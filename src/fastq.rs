//! High-throughput FASTQ / FASTA stream parser supporting gzip compression,
//! Phred quality score auto-detection (Phred-33 vs Phred-64), and linked-read barcode parsing.

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Quality score encoding format.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum PhredEncoding {
    /// Phred+33 (Sanger, Illumina 1.8+, standard).
    Phred33,
    /// Phred+64 (Solexa, Illumina 1.3-1.7).
    Phred64,
}

impl PhredEncoding {
    #[inline(always)]
    pub fn offset(&self) -> u8 {
        match self {
            Self::Phred33 => 33,
            Self::Phred64 => 64,
        }
    }
}

/// Auto-detects Phred quality score offset from a slice of ASCII quality characters.
pub fn detect_phred_offset(qualities: &[u8]) -> PhredEncoding {
    let mut min_q = 255u8;
    let mut max_q = 0u8;

    for &q in qualities {
        if q >= 33 {
            min_q = min_q.min(q);
            max_q = max_q.max(q);
        }
    }

    // Characters between 33 and 58 ('!' through ':') can only exist in Phred+33.
    if min_q < 59 {
        PhredEncoding::Phred33
    } else if max_q > 74 {
        // High ASCII qualities with min >= 64 indicate Phred+64.
        PhredEncoding::Phred64
    } else {
        // Default to the modern genomics standard (Phred+33)
        PhredEncoding::Phred33
    }
}

/// Configuration for read quality filtering and adapter/tail trimming.
#[derive(Clone, Debug)]
pub struct QualityFilterConfig {
    /// Minimum Phred quality score (e.g. Q10 or Q15). Bases below this at the 3' end are trimmed.
    pub min_quality: u8,
    /// Minimum read length after trimming to keep.
    pub min_read_len: usize,
    /// Maximum number of uncalled 'N' bases allowed in a read.
    pub max_ns: usize,
    /// Manual override for Phred encoding, or None to auto-detect.
    pub phred_override: Option<PhredEncoding>,
}

impl Default for QualityFilterConfig {
    fn default() -> Self {
        Self {
            min_quality: 10,
            min_read_len: 30,
            max_ns: 3,
            phred_override: None,
        }
    }
}

/// A sequencing read accompanied by an optional linked-read barcode (e.g. 10x Genomics or TELL-Seq).
#[derive(Clone, Debug)]
pub struct LinkedReadRecord {
    pub barcode: Option<u64>,
    pub seq: Vec<u8>,
}

/// Extracts a 10x / TELL-Seq barcode from a FASTQ header line (e.g. `BX:Z:ACGT...` or `BC:Z:ACGT...`).
pub fn extract_barcode_from_header(header: &[u8]) -> Option<u64> {
    let header_str = std::str::from_utf8(header).ok()?;
    for token in header_str.split_whitespace() {
        if let Some(bc_str) = token.strip_prefix("BX:Z:") {
            return encode_barcode_str(bc_str);
        }
        if let Some(bc_str) = token.strip_prefix("BC:Z:") {
            return encode_barcode_str(bc_str);
        }
    }
    // Check delimiter format `@read_id#ACGT...`
    if let Some(pos) = header_str.find('#') {
        let tail = &header_str[pos + 1..];
        let bc_str = tail.split('/').next().unwrap_or(tail);
        return encode_barcode_str(bc_str);
    }
    None
}

/// Encodes an ASCII barcode (up to 32 bp) into a 64-bit integer, or hashes it if longer.
fn encode_barcode_str(bc_str: &str) -> Option<u64> {
    let clean_str = bc_str.trim_end_matches("-1");
    if clean_str.is_empty() {
        return None;
    }
    if clean_str.len() <= 32
        && clean_str
            .bytes()
            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
    {
        let mut val = 0u64;
        for b in clean_str.bytes() {
            let two_bit = match b {
                b'A' | b'a' => 0u64,
                b'C' | b'c' => 1u64,
                b'G' | b'g' => 2u64,
                b'T' | b't' => 3u64,
                _ => return None,
            };
            val = (val << 2) | two_bit;
        }
        Some(val)
    } else {
        // 64-bit FNV-1a hash fallback for non-canonical barcodes
        let mut hash = 0xcbf29ce484222325u64;
        for b in clean_str.bytes() {
            hash ^= b as u64;
            hash = hash.wrapping_mul(0x100000001b3);
        }
        Some(hash)
    }
}

/// Reads all valid DNA sequences from a FASTQ or FASTA file (plain text or gzip).
pub fn parse_reads_from_file<P: AsRef<Path>>(path: P) -> Result<Vec<Vec<u8>>> {
    parse_reads_with_quality(path, None)
}

/// Reads all valid DNA sequences with optional quality-based trimming and filtering.
pub fn parse_reads_with_quality<P: AsRef<Path>>(
    path: P,
    filter: Option<&QualityFilterConfig>,
) -> Result<Vec<Vec<u8>>> {
    let path = path.as_ref();
    let file = File::open(path).with_context(|| format!("Failed to open file: {:?}", path))?;

    let is_gz = path.extension().is_some_and(|ext| ext == "gz");
    let reader: Box<dyn Read + Send> = if is_gz {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    let mut buf_reader = BufReader::with_capacity(512 * 1024, reader);
    let mut sequences = Vec::with_capacity(100_000);

    // Peak at first byte to determine format: '@' for FASTQ, '>' for FASTA
    let mut first_byte = [0u8; 1];
    let n = buf_reader.read(&mut first_byte)?;
    if n == 0 {
        return Ok(sequences); // empty file
    }

    let is_fastq = first_byte[0] == b'@';

    if is_fastq {
        let mut header = Vec::with_capacity(512);
        let mut seq_buf = Vec::with_capacity(512);
        let mut qual_buf = Vec::with_capacity(512);
        let mut sep = Vec::with_capacity(128);

        // Discard rest of first header line
        buf_reader.read_until(b'\n', &mut header)?;

        // Phred offset: use explicit override if set, otherwise auto-detect from the first
        // quality line encountered (distinguishes Phred+33 from legacy Phred+64 data).
        let phred_override = filter.and_then(|cfg| cfg.phred_override);
        let mut phred_offset = phred_override.unwrap_or(PhredEncoding::Phred33).offset();
        let mut phred_detected = phred_override.is_some();

        loop {
            // Line 2: Sequence
            seq_buf.clear();
            if buf_reader.read_until(b'\n', &mut seq_buf)? == 0 {
                break;
            }
            while seq_buf.ends_with(b"\n") || seq_buf.ends_with(b"\r") {
                seq_buf.pop();
            }

            // Line 3: '+' header
            sep.clear();
            if buf_reader.read_until(b'\n', &mut sep)? == 0 {
                break;
            }

            // Line 4: Quality string
            qual_buf.clear();
            if buf_reader.read_until(b'\n', &mut qual_buf)? == 0 {
                break;
            }
            while qual_buf.ends_with(b"\n") || qual_buf.ends_with(b"\r") {
                qual_buf.pop();
            }

            // Auto-detect Phred encoding from first non-empty quality line when no override given
            if !phred_detected && !qual_buf.is_empty() {
                phred_offset = detect_phred_offset(&qual_buf).offset();
                phred_detected = true;
            }

            if !seq_buf.is_empty() {
                if let Some(cfg) = filter {
                    // Check N base count
                    let n_count = seq_buf.iter().filter(|&&b| b == b'N' || b == b'n').count();
                    if n_count <= cfg.max_ns {
                        // 3' quality score trimming
                        let mut valid_len = seq_buf.len().min(qual_buf.len());
                        while valid_len > 0 {
                            let q_char = qual_buf[valid_len - 1];
                            let q_score = q_char.saturating_sub(phred_offset);
                            if q_score >= cfg.min_quality {
                                break;
                            }
                            valid_len -= 1;
                        }

                        if valid_len >= cfg.min_read_len {
                            sequences.push(seq_buf[..valid_len].to_vec());
                        }
                    }
                } else {
                    sequences.push(seq_buf.clone());
                }
            }

            // Read Line 1 (Header) of the next record
            header.clear();
            if buf_reader.read_until(b'\n', &mut header)? == 0 {
                break;
            }
        }
    } else {
        // FASTA format
        let mut line = String::with_capacity(1024);
        let mut current_seq = Vec::with_capacity(1024);

        // Finish reading remainder of initial header line
        buf_reader.read_line(&mut line)?;

        loop {
            line.clear();
            if buf_reader.read_line(&mut line)? == 0 {
                break;
            }
            let trimmed = line.trim_end();
            if trimmed.starts_with('>') {
                if !current_seq.is_empty() {
                    sequences.push(std::mem::take(&mut current_seq));
                }
            } else {
                current_seq.extend_from_slice(trimmed.as_bytes());
            }
        }
        if !current_seq.is_empty() {
            sequences.push(current_seq);
        }
    }

    Ok(sequences)
}

/// Parses linked reads with barcode extraction for SPlitteR repeat resolution.
pub fn parse_linked_reads_from_file<P: AsRef<Path>>(path: P) -> Result<Vec<LinkedReadRecord>> {
    let path = path.as_ref();
    let file = File::open(path).with_context(|| format!("Failed to open file: {:?}", path))?;

    let is_gz = path.extension().is_some_and(|ext| ext == "gz");
    let reader: Box<dyn Read + Send> = if is_gz {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    let mut buf_reader = BufReader::with_capacity(512 * 1024, reader);
    let mut records = Vec::with_capacity(100_000);

    let mut header = Vec::with_capacity(512);
    let mut seq = Vec::with_capacity(512);
    let mut sep = Vec::with_capacity(128);
    let mut qual = Vec::with_capacity(512);

    loop {
        header.clear();
        if buf_reader.read_until(b'\n', &mut header)? == 0 {
            break;
        }
        seq.clear();
        if buf_reader.read_until(b'\n', &mut seq)? == 0 {
            break;
        }
        sep.clear();
        if buf_reader.read_until(b'\n', &mut sep)? == 0 {
            break;
        }
        qual.clear();
        if buf_reader.read_until(b'\n', &mut qual)? == 0 {
            break;
        }

        while seq.ends_with(b"\n") || seq.ends_with(b"\r") {
            seq.pop();
        }

        if !seq.is_empty() {
            let barcode = extract_barcode_from_header(&header);
            records.push(LinkedReadRecord {
                barcode,
                seq: seq.clone(),
            });
        }
    }

    Ok(records)
}
