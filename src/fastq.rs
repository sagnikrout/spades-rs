//! High-throughput FASTQ / FASTA stream parser supporting gzip compression.

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// A simple DNA sequencing read.
#[derive(Clone, Debug)]
pub struct ReadRecord {
    pub id: String,
    pub seq: Vec<u8>,
}

/// Reads all valid DNA sequences from a FASTQ or FASTA file (plain text or gzip).
pub fn parse_reads_from_file<P: AsRef<Path>>(path: P) -> Result<Vec<Vec<u8>>> {
    let path = path.as_ref();
    let file = File::open(path).with_context(|| format!("Failed to open file: {:?}", path))?;

    let is_gz = path.extension().map_or(false, |ext| ext == "gz");
    let reader: Box<dyn Read + Send> = if is_gz {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    let mut buf_reader = BufReader::with_capacity(256 * 1024, reader);
    let mut sequences = Vec::with_capacity(100_000);

    let mut line = String::with_capacity(1024);

    // Peak at first byte to determine format: '@' for FASTQ, '>' for FASTA
    let mut first_byte = [0u8; 1];
    let n = buf_reader.read(&mut first_byte)?;
    if n == 0 {
        return Ok(sequences); // empty file
    }

    let is_fastq = first_byte[0] == b'@';
    line.push(first_byte[0] as char);

    if is_fastq {
        // FASTQ format: 4 lines per record:
        // @ID
        // SEQ
        // +
        // QUAL
        loop {
            // Finish reading remainder of the header line
            line.clear();
            if buf_reader.read_line(&mut line)? == 0 {
                break;
            }

            // Line 2: Sequence
            let mut seq_line = String::with_capacity(256);
            if buf_reader.read_line(&mut seq_line)? == 0 {
                break;
            }
            let trimmed_seq = seq_line.trim_end().as_bytes().to_vec();

            // Line 3: '+' header
            line.clear();
            if buf_reader.read_line(&mut line)? == 0 {
                break;
            }

            // Line 4: Quality string
            line.clear();
            if buf_reader.read_line(&mut line)? == 0 {
                break;
            }

            if !trimmed_seq.is_empty() {
                sequences.push(trimmed_seq);
            }
        }
    } else {
        // FASTA format
        let mut current_seq = Vec::with_capacity(1024);
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
