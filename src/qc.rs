//! Stage 2: Read QC & Preprocessing Module.
//!
//! High-throughput paired-end and single-end FASTQ preprocessor performing:
//! - Adapter trimming (Illumina Universal, TruSeq, Nextera, or custom adapters)
//! - Poly-G tail trimming (two-color Illumina NextSeq/NovaSeq artifact)
//! - Poly-N and 3'/5' low-quality end trimming
//! - Sliding-window quality filtering (Phred-33 / Phred-64 auto-detected)
//! - Length and ambiguous base filtering
//! - Paired-end read synchronization (synchronized output or singleton routing)
//! - QC summary metrics: Q20%, Q30%, GC%, total bases, read survival rate

use crate::fastq::{detect_phred_offset, PhredEncoding};
use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

/// Standard known adapter sequences.
pub const ADAPTER_ILLUMINA_TRUSEQ_R1: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
pub const ADAPTER_ILLUMINA_TRUSEQ_R2: &[u8] = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
pub const ADAPTER_NEXTERA: &[u8] = b"CTGTCTCTTATACACATCT";

/// A single FASTQ record.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FastqRecord {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

impl FastqRecord {
    #[inline]
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

/// Configuration options for read preprocessing and quality control.
#[derive(Clone, Debug)]
pub struct QcConfig {
    /// Minimum Phred quality score for window trimming (default 15)
    pub min_quality: u8,
    /// Sliding window size for quality trimming (default 4)
    pub window_size: usize,
    /// Minimum read length after trimming to keep (default 30 bp)
    pub min_length: usize,
    /// Maximum number of uncalled 'N' bases allowed (default 3)
    pub max_ns: usize,
    /// Minimum length of trailing poly-G run to trim (e.g. 5). Set to 0 to disable.
    pub trim_poly_g: usize,
    /// Minimum overlap length to detect and trim an adapter at read 3' end (default 10)
    pub min_adapter_overlap: usize,
    /// Custom adapter sequences to trim from forward/R1 reads
    pub adapters_r1: Vec<Vec<u8>>,
    /// Custom adapter sequences to trim from reverse/R2 reads
    pub adapters_r2: Vec<Vec<u8>>,
    /// Manual override for Phred encoding format (None = auto-detect)
    pub phred_override: Option<PhredEncoding>,
    /// Number of worker threads (None = Rayon default)
    pub threads: Option<usize>,
}

impl Default for QcConfig {
    fn default() -> Self {
        Self {
            min_quality: 15,
            window_size: 4,
            min_length: 30,
            max_ns: 3,
            trim_poly_g: 5,
            min_adapter_overlap: 10,
            adapters_r1: vec![
                ADAPTER_ILLUMINA_TRUSEQ_R1.to_vec(),
                ADAPTER_NEXTERA.to_vec(),
            ],
            adapters_r2: vec![
                ADAPTER_ILLUMINA_TRUSEQ_R2.to_vec(),
                ADAPTER_NEXTERA.to_vec(),
            ],
            phred_override: None,
            threads: None,
        }
    }
}

/// Statistics collected before and after quality control filtering.
#[derive(Default, Clone, Debug)]
pub struct QcReport {
    pub raw_reads: usize,
    pub raw_bases: usize,
    pub raw_q20_bases: usize,
    pub raw_q30_bases: usize,
    pub raw_gc_bases: usize,

    pub clean_reads: usize,
    pub clean_bases: usize,
    pub clean_q20_bases: usize,
    pub clean_q30_bases: usize,
    pub clean_gc_bases: usize,

    pub adapter_trimmed_reads: usize,
    pub poly_g_trimmed_reads: usize,
    pub quality_trimmed_reads: usize,
    pub dropped_too_short: usize,
    pub dropped_too_many_ns: usize,
    pub dropped_low_quality: usize,

    pub detected_phred: String,
}

impl QcReport {
    pub fn raw_q20_pct(&self) -> f64 {
        if self.raw_bases == 0 {
            0.0
        } else {
            (self.raw_q20_bases as f64 / self.raw_bases as f64) * 100.0
        }
    }

    pub fn raw_q30_pct(&self) -> f64 {
        if self.raw_bases == 0 {
            0.0
        } else {
            (self.raw_q30_bases as f64 / self.raw_bases as f64) * 100.0
        }
    }

    pub fn raw_gc_pct(&self) -> f64 {
        if self.raw_bases == 0 {
            0.0
        } else {
            (self.raw_gc_bases as f64 / self.raw_bases as f64) * 100.0
        }
    }

    pub fn clean_q20_pct(&self) -> f64 {
        if self.clean_bases == 0 {
            0.0
        } else {
            (self.clean_q20_bases as f64 / self.clean_bases as f64) * 100.0
        }
    }

    pub fn clean_q30_pct(&self) -> f64 {
        if self.clean_bases == 0 {
            0.0
        } else {
            (self.clean_q30_bases as f64 / self.clean_bases as f64) * 100.0
        }
    }

    pub fn clean_gc_pct(&self) -> f64 {
        if self.clean_bases == 0 {
            0.0
        } else {
            (self.clean_gc_bases as f64 / self.clean_bases as f64) * 100.0
        }
    }

    pub fn survival_rate_pct(&self) -> f64 {
        if self.raw_reads == 0 {
            0.0
        } else {
            (self.clean_reads as f64 / self.raw_reads as f64) * 100.0
        }
    }

    pub fn to_json(&self) -> String {
        format!(
            "{{\n  \
              \"detected_phred\": \"{}\",\n  \
              \"raw\": {{\n    \
                \"reads\": {},\n    \
                \"bases\": {},\n    \
                \"q20_pct\": {:.2},\n    \
                \"q30_pct\": {:.2},\n    \
                \"gc_pct\": {:.2}\n  \
              }},\n  \
              \"clean\": {{\n    \
                \"reads\": {},\n    \
                \"bases\": {},\n    \
                \"q20_pct\": {:.2},\n    \
                \"q30_pct\": {:.2},\n    \
                \"gc_pct\": {:.2},\n    \
                \"survival_rate_pct\": {:.2}\n  \
              }},\n  \
              \"trimming\": {{\n    \
                \"adapter_trimmed_reads\": {},\n    \
                \"poly_g_trimmed_reads\": {},\n    \
                \"quality_trimmed_reads\": {},\n    \
                \"dropped_too_short\": {},\n    \
                \"dropped_too_many_ns\": {},\n    \
                \"dropped_low_quality\": {}\n  \
              }}\n\
            }}",
            self.detected_phred,
            self.raw_reads,
            self.raw_bases,
            self.raw_q20_pct(),
            self.raw_q30_pct(),
            self.raw_gc_pct(),
            self.clean_reads,
            self.clean_bases,
            self.clean_q20_pct(),
            self.clean_q30_pct(),
            self.clean_gc_pct(),
            self.survival_rate_pct(),
            self.adapter_trimmed_reads,
            self.poly_g_trimmed_reads,
            self.quality_trimmed_reads,
            self.dropped_too_short,
            self.dropped_too_many_ns,
            self.dropped_low_quality
        )
    }

    pub fn print_summary(&self) {
        println!("============================================================");
        println!("           STAGE 2: READ PREPROCESSING & QC SUMMARY         ");
        println!("============================================================");
        println!("  Phred Quality Encoding:     {}", self.detected_phred);
        println!("  Raw Reads:                  {:>12}", self.raw_reads);
        println!("  Raw Bases:                  {:>12} bp", self.raw_bases);
        println!(
            "  Raw Q20 / Q30:              {:>10.2}% / {:>.2}%",
            self.raw_q20_pct(),
            self.raw_q30_pct()
        );
        println!("  Raw GC Content:             {:>10.2}%", self.raw_gc_pct());
        println!("------------------------------------------------------------");
        println!(
            "  Clean Reads Passed:         {:>12} ({:.2}%)",
            self.clean_reads,
            self.survival_rate_pct()
        );
        println!("  Clean Bases Passed:         {:>12} bp", self.clean_bases);
        println!(
            "  Clean Q20 / Q30:            {:>10.2}% / {:>.2}%",
            self.clean_q20_pct(),
            self.clean_q30_pct()
        );
        println!(
            "  Clean GC Content:           {:>10.2}%",
            self.clean_gc_pct()
        );
        println!("------------------------------------------------------------");
        println!(
            "  Reads with Adapters Trimmed:{:>12}",
            self.adapter_trimmed_reads
        );
        println!(
            "  Reads with Poly-G Trimmed:  {:>12}",
            self.poly_g_trimmed_reads
        );
        println!(
            "  Reads Quality Trimmed:      {:>12}",
            self.quality_trimmed_reads
        );
        println!(
            "  Dropped (Length < min):     {:>12}",
            self.dropped_too_short
        );
        println!(
            "  Dropped (Excessive Ns):     {:>12}",
            self.dropped_too_many_ns
        );
        println!(
            "  Dropped (Mean Q < threshold):{:>11}",
            self.dropped_low_quality
        );
        println!("============================================================");
    }
}

/// Trims leading and trailing ambiguous 'N' / 'n' bases and their qualities.
pub fn trim_poly_n(record: &mut FastqRecord) {
    let mut start = 0;
    while start < record.seq.len() && (record.seq[start] == b'N' || record.seq[start] == b'n') {
        start += 1;
    }
    let mut end = record.seq.len();
    while end > start && (record.seq[end - 1] == b'N' || record.seq[end - 1] == b'n') {
        end -= 1;
    }
    if start > 0 || end < record.seq.len() {
        record.seq = record.seq[start..end].to_vec();
        if record.qual.len() >= end {
            record.qual = record.qual[start..end].to_vec();
        }
    }
}

/// Trims trailing poly-G runs (NextSeq/NovaSeq two-color artifact) if run length >= min_g_len.
pub fn trim_poly_g(record: &mut FastqRecord, min_g_len: usize) -> bool {
    if min_g_len == 0 || record.seq.len() < min_g_len {
        return false;
    }
    let mut g_count = 0;
    for &b in record.seq.iter().rev() {
        if b == b'G' || b == b'g' {
            g_count += 1;
        } else {
            break;
        }
    }
    if g_count >= min_g_len {
        let new_len = record.seq.len() - g_count;
        record.seq.truncate(new_len);
        if record.qual.len() >= new_len {
            record.qual.truncate(new_len);
        }
        true
    } else {
        false
    }
}

/// Trims adapter sequences from the 3' end of the read.
/// Returns true if an adapter was detected and trimmed.
pub fn trim_adapters(record: &mut FastqRecord, adapters: &[Vec<u8>], min_overlap: usize) -> bool {
    if record.seq.len() < min_overlap || adapters.is_empty() {
        return false;
    }

    let mut best_trim_pos: Option<usize> = None;
    let r_len = record.seq.len();

    for adapter in adapters {
        let ad_len = adapter.len();
        if ad_len < min_overlap {
            continue;
        }

        // 1. Partial adapter at 3' end (overlap of read suffix with adapter prefix)
        let min_start = r_len.saturating_sub(ad_len);
        let max_start = r_len.saturating_sub(min_overlap);

        for start in min_start..=max_start {
            let overlap = r_len - start;
            let max_allowed = (overlap / 10).max(if overlap >= 12 { 1 } else { 0 });
            let mut mismatches = 0;

            for (i, &ab_byte) in adapter.iter().take(overlap).enumerate() {
                let rb = record.seq[start + i].to_ascii_uppercase();
                let ab = ab_byte.to_ascii_uppercase();
                if rb != ab {
                    mismatches += 1;
                    if mismatches > max_allowed {
                        break;
                    }
                }
            }

            if mismatches <= max_allowed {
                best_trim_pos = Some(match best_trim_pos {
                    Some(prev) => prev.min(start),
                    None => start,
                });
                break;
            }
        }

        // 2. Full adapter embedded in read (read length > insert + ad_len)
        if r_len > ad_len {
            let max_full_start = r_len - ad_len;
            let max_allowed = (ad_len / 10).max(1);

            for start in (0..max_full_start).rev() {
                let mut mismatches = 0;
                for (i, &ab_byte) in adapter.iter().take(ad_len).enumerate() {
                    let rb = record.seq[start + i].to_ascii_uppercase();
                    let ab = ab_byte.to_ascii_uppercase();
                    if rb != ab {
                        mismatches += 1;
                        if mismatches > max_allowed {
                            break;
                        }
                    }
                }

                if mismatches <= max_allowed {
                    best_trim_pos = Some(match best_trim_pos {
                        Some(prev) => prev.min(start),
                        None => start,
                    });
                    break;
                }
            }
        }
    }

    if let Some(pos) = best_trim_pos {
        record.seq.truncate(pos);
        if record.qual.len() >= pos {
            record.qual.truncate(pos);
        }
        true
    } else {
        false
    }
}

/// Performs sliding-window quality trimming from the 3' end toward the 5' end.
/// If the average quality in a window drops below min_qual, clips the read.
pub fn sliding_window_trim(
    record: &mut FastqRecord,
    window_size: usize,
    min_qual: u8,
    phred_offset: u8,
) -> bool {
    let len = record.seq.len().min(record.qual.len());
    if len < window_size || window_size == 0 {
        return false;
    }

    let scores: Vec<u8> = record.qual[..len]
        .iter()
        .map(|&q| q.saturating_sub(phred_offset))
        .collect();

    let target_sum = (min_qual as usize) * window_size;
    let mut end = len;

    // 1. Trim trailing individual bases below min_qual
    while end > 0 && scores[end - 1] < min_qual {
        end -= 1;
    }

    // 2. Sliding window check from right to left
    while end >= window_size {
        let win_sum: usize = scores[end - window_size..end]
            .iter()
            .map(|&s| s as usize)
            .sum();
        if win_sum < target_sum {
            end -= 1;
            while end > 0 && scores[end - 1] < min_qual {
                end -= 1;
            }
        } else {
            break;
        }
    }

    if end < len {
        record.seq.truncate(end);
        record.qual.truncate(end);
        true
    } else {
        false
    }
}

/// Evaluates if a processed record passes length and quality thresholds.
pub fn passes_filter(
    record: &FastqRecord,
    cfg: &QcConfig,
    phred_offset: u8,
) -> Result<(), &'static str> {
    if record.seq.len() < cfg.min_length {
        return Err("too_short");
    }

    let n_count = record
        .seq
        .iter()
        .filter(|&&b| b == b'N' || b == b'n')
        .count();
    if n_count > cfg.max_ns {
        return Err("too_many_ns");
    }

    if !record.qual.is_empty() {
        let total_q: usize = record
            .qual
            .iter()
            .map(|&q| q.saturating_sub(phred_offset) as usize)
            .sum();
        let mean_q = total_q as f64 / record.qual.len() as f64;
        if mean_q < (cfg.min_quality as f64 * 0.75) {
            return Err("low_quality");
        }
    }

    Ok(())
}

/// Reads all FastqRecord entries from a file (FASTQ / FASTQ.GZ).
pub fn read_fastq_records<P: AsRef<Path>>(path: P) -> Result<(Vec<FastqRecord>, PhredEncoding)> {
    let path = path.as_ref();
    let file = File::open(path).with_context(|| format!("Failed to open file: {:?}", path))?;
    let is_gz = path.extension().is_some_and(|ext| ext == "gz");

    let reader: Box<dyn Read + Send> = if is_gz {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    let mut buf = BufReader::with_capacity(1024 * 1024, reader);
    let mut records = Vec::with_capacity(65_536);

    let mut header = Vec::with_capacity(256);
    let mut seq = Vec::with_capacity(256);
    let mut sep = Vec::with_capacity(32);
    let mut qual = Vec::with_capacity(256);

    let mut sample_qualities = Vec::new();

    loop {
        header.clear();
        if buf.read_until(b'\n', &mut header)? == 0 {
            break;
        }
        while header.ends_with(b"\n") || header.ends_with(b"\r") {
            header.pop();
        }
        if header.is_empty() {
            continue;
        }

        seq.clear();
        if buf.read_until(b'\n', &mut seq)? == 0 {
            break;
        }
        while seq.ends_with(b"\n") || seq.ends_with(b"\r") {
            seq.pop();
        }

        sep.clear();
        if buf.read_until(b'\n', &mut sep)? == 0 {
            break;
        }

        qual.clear();
        if buf.read_until(b'\n', &mut qual)? == 0 {
            break;
        }
        while qual.ends_with(b"\n") || qual.ends_with(b"\r") {
            qual.pop();
        }

        if sample_qualities.len() < 10_000 && !qual.is_empty() {
            sample_qualities.extend_from_slice(&qual);
        }

        records.push(FastqRecord {
            id: header.clone(),
            seq: seq.clone(),
            qual: qual.clone(),
        });
    }

    let encoding = detect_phred_offset(&sample_qualities);
    Ok((records, encoding))
}

/// Opens a FASTQ writer that automatically gzips if path ends with `.gz`.
pub fn open_fastq_writer<P: AsRef<Path>>(path: P) -> Result<Box<dyn Write>> {
    let path = path.as_ref();
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let file = File::create(path).with_context(|| format!("Failed to create: {:?}", path))?;
    let is_gz = path.extension().is_some_and(|ext| ext == "gz");

    if is_gz {
        let enc = GzEncoder::new(file, Compression::fast());
        Ok(Box::new(BufWriter::with_capacity(512 * 1024, enc)))
    } else {
        Ok(Box::new(BufWriter::with_capacity(512 * 1024, file)))
    }
}

/// Writes a single FASTQ record to the output writer.
#[inline]
pub fn write_fastq_record<W: Write>(writer: &mut W, rec: &FastqRecord) -> Result<()> {
    writer.write_all(&rec.id)?;
    writer.write_all(b"\n")?;
    writer.write_all(&rec.seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(&rec.qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}

/// Runs Stage 2 preprocessing on paired-end reads.
pub fn process_paired_reads(
    r1_path: &Path,
    r2_path: &Path,
    out1_path: &Path,
    out2_path: &Path,
    unpaired_out_path: Option<&Path>,
    cfg: &QcConfig,
) -> Result<QcReport> {
    let (mut reads1, enc1) = read_fastq_records(r1_path)?;
    let (mut reads2, _enc2) = read_fastq_records(r2_path)?;

    if reads1.len() != reads2.len() {
        anyhow::bail!(
            "Paired FASTQ record counts do not match: R1 has {} reads, R2 has {} reads",
            reads1.len(),
            reads2.len()
        );
    }

    let encoding = cfg.phred_override.unwrap_or(enc1);
    let offset = encoding.offset();

    let mut report = QcReport {
        detected_phred: format!("{:?} (offset {})", encoding, offset),
        ..Default::default()
    };

    // Calculate raw statistics
    for rec in &reads1 {
        accumulate_stats(
            &mut report.raw_reads,
            &mut report.raw_bases,
            &mut report.raw_q20_bases,
            &mut report.raw_q30_bases,
            &mut report.raw_gc_bases,
            rec,
            offset,
        );
    }
    for rec in &reads2 {
        accumulate_stats(
            &mut report.raw_reads,
            &mut report.raw_bases,
            &mut report.raw_q20_bases,
            &mut report.raw_q30_bases,
            &mut report.raw_gc_bases,
            rec,
            offset,
        );
    }

    let mut out1 = open_fastq_writer(out1_path)?;
    let mut out2 = open_fastq_writer(out2_path)?;
    let mut out_unpaired = match unpaired_out_path {
        Some(p) => Some(open_fastq_writer(p)?),
        None => None,
    };

    // Process records
    for i in 0..reads1.len() {
        let rec1 = &mut reads1[i];
        let rec2 = &mut reads2[i];

        // 1. Poly-N trimming
        trim_poly_n(rec1);
        trim_poly_n(rec2);

        // 2. Poly-G trimming
        if trim_poly_g(rec1, cfg.trim_poly_g) {
            report.poly_g_trimmed_reads += 1;
        }
        if trim_poly_g(rec2, cfg.trim_poly_g) {
            report.poly_g_trimmed_reads += 1;
        }

        // 3. Adapter trimming
        if trim_adapters(rec1, &cfg.adapters_r1, cfg.min_adapter_overlap) {
            report.adapter_trimmed_reads += 1;
        }
        if trim_adapters(rec2, &cfg.adapters_r2, cfg.min_adapter_overlap) {
            report.adapter_trimmed_reads += 1;
        }

        // 4. Sliding-window quality trimming
        if sliding_window_trim(rec1, cfg.window_size, cfg.min_quality, offset) {
            report.quality_trimmed_reads += 1;
        }
        if sliding_window_trim(rec2, cfg.window_size, cfg.min_quality, offset) {
            report.quality_trimmed_reads += 1;
        }

        // 5. Check filters
        let res1 = passes_filter(rec1, cfg, offset);
        let res2 = passes_filter(rec2, cfg, offset);

        match (res1, res2) {
            (Ok(()), Ok(())) => {
                write_fastq_record(&mut out1, rec1)?;
                write_fastq_record(&mut out2, rec2)?;
                accumulate_stats(
                    &mut report.clean_reads,
                    &mut report.clean_bases,
                    &mut report.clean_q20_bases,
                    &mut report.clean_q30_bases,
                    &mut report.clean_gc_bases,
                    rec1,
                    offset,
                );
                accumulate_stats(
                    &mut report.clean_reads,
                    &mut report.clean_bases,
                    &mut report.clean_q20_bases,
                    &mut report.clean_q30_bases,
                    &mut report.clean_gc_bases,
                    rec2,
                    offset,
                );
            }
            (Ok(()), Err(reason2)) => {
                record_drop(&mut report, reason2);
                if let Some(ref mut u_writer) = out_unpaired {
                    write_fastq_record(u_writer, rec1)?;
                    accumulate_stats(
                        &mut report.clean_reads,
                        &mut report.clean_bases,
                        &mut report.clean_q20_bases,
                        &mut report.clean_q30_bases,
                        &mut report.clean_gc_bases,
                        rec1,
                        offset,
                    );
                } else {
                    report.dropped_too_short += 1; // R1 dropped because mate failed
                }
            }
            (Err(reason1), Ok(())) => {
                record_drop(&mut report, reason1);
                if let Some(ref mut u_writer) = out_unpaired {
                    write_fastq_record(u_writer, rec2)?;
                    accumulate_stats(
                        &mut report.clean_reads,
                        &mut report.clean_bases,
                        &mut report.clean_q20_bases,
                        &mut report.clean_q30_bases,
                        &mut report.clean_gc_bases,
                        rec2,
                        offset,
                    );
                } else {
                    report.dropped_too_short += 1; // R2 dropped because mate failed
                }
            }
            (Err(r1), Err(r2)) => {
                record_drop(&mut report, r1);
                record_drop(&mut report, r2);
            }
        }
    }

    Ok(report)
}

/// Runs Stage 2 preprocessing on single-end reads.
pub fn process_single_reads(
    input_path: &Path,
    output_path: &Path,
    cfg: &QcConfig,
) -> Result<QcReport> {
    let (mut reads, enc) = read_fastq_records(input_path)?;
    let encoding = cfg.phred_override.unwrap_or(enc);
    let offset = encoding.offset();

    let mut report = QcReport {
        detected_phred: format!("{:?} (offset {})", encoding, offset),
        ..Default::default()
    };

    for rec in &reads {
        accumulate_stats(
            &mut report.raw_reads,
            &mut report.raw_bases,
            &mut report.raw_q20_bases,
            &mut report.raw_q30_bases,
            &mut report.raw_gc_bases,
            rec,
            offset,
        );
    }

    let mut writer = open_fastq_writer(output_path)?;

    for rec in &mut reads {
        trim_poly_n(rec);
        if trim_poly_g(rec, cfg.trim_poly_g) {
            report.poly_g_trimmed_reads += 1;
        }
        if trim_adapters(rec, &cfg.adapters_r1, cfg.min_adapter_overlap) {
            report.adapter_trimmed_reads += 1;
        }
        if sliding_window_trim(rec, cfg.window_size, cfg.min_quality, offset) {
            report.quality_trimmed_reads += 1;
        }

        match passes_filter(rec, cfg, offset) {
            Ok(()) => {
                write_fastq_record(&mut writer, rec)?;
                accumulate_stats(
                    &mut report.clean_reads,
                    &mut report.clean_bases,
                    &mut report.clean_q20_bases,
                    &mut report.clean_q30_bases,
                    &mut report.clean_gc_bases,
                    rec,
                    offset,
                );
            }
            Err(reason) => {
                record_drop(&mut report, reason);
            }
        }
    }

    Ok(report)
}

#[inline]
fn accumulate_stats(
    n_reads: &mut usize,
    n_bases: &mut usize,
    n_q20: &mut usize,
    n_q30: &mut usize,
    n_gc: &mut usize,
    rec: &FastqRecord,
    phred_offset: u8,
) {
    *n_reads += 1;
    *n_bases += rec.seq.len();
    for &b in &rec.seq {
        if b == b'G' || b == b'C' || b == b'g' || b == b'c' {
            *n_gc += 1;
        }
    }
    for &q in &rec.qual {
        let score = q.saturating_sub(phred_offset);
        if score >= 20 {
            *n_q20 += 1;
        }
        if score >= 30 {
            *n_q30 += 1;
        }
    }
}

#[inline]
fn record_drop(report: &mut QcReport, reason: &str) {
    match reason {
        "too_short" => report.dropped_too_short += 1,
        "too_many_ns" => report.dropped_too_many_ns += 1,
        "low_quality" => report.dropped_low_quality += 1,
        _ => report.dropped_too_short += 1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim_poly_n() {
        let mut rec = FastqRecord {
            id: b"@read1".to_vec(),
            seq: b"NNNACGTNN".to_vec(),
            qual: b"IIIIIIIII".to_vec(),
        };
        trim_poly_n(&mut rec);
        assert_eq!(rec.seq, b"ACGT");
        assert_eq!(rec.qual, b"IIII");
    }

    #[test]
    fn test_trim_poly_g() {
        let mut rec = FastqRecord {
            id: b"@read1".to_vec(),
            seq: b"ACGTGGGGGG".to_vec(),
            qual: b"IIIIIIIIII".to_vec(),
        };
        assert!(trim_poly_g(&mut rec, 5));
        assert_eq!(rec.seq, b"ACGT");
        assert_eq!(rec.qual, b"IIII");
    }

    #[test]
    fn test_trim_adapters() {
        let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec();
        let mut rec = FastqRecord {
            id: b"@read1".to_vec(),
            seq: b"GATTACAGATCGGAAGAGCACACGTCTGAA".to_vec(),
            qual: vec![b'I'; 30],
        };
        let trimmed = trim_adapters(&mut rec, &[adapter], 10);
        assert!(trimmed);
        assert_eq!(rec.seq, b"GATTAC");
    }

    #[test]
    fn test_sliding_window_trim() {
        let mut rec = FastqRecord {
            id: b"@read1".to_vec(),
            seq: b"ACGTACGTACGT".to_vec(),
            // First 8 high qual (score 40, 'I'=73), last 4 low qual (score 2, '$'=35)
            qual: b"IIIIIIII$$$$".to_vec(),
        };
        // phred offset 33. Window 4, min qual 15
        assert!(sliding_window_trim(&mut rec, 4, 15, 33));
        assert_eq!(rec.seq, b"ACGTACGT");
    }

    #[test]
    fn test_passes_filter() {
        let cfg = QcConfig {
            min_length: 20,
            max_ns: 2,
            min_quality: 10,
            ..Default::default()
        };
        let good_rec = FastqRecord {
            id: b"@r".to_vec(),
            seq: b"ACGTACGTACGTACGTACGT".to_vec(),
            qual: b"IIIIIIIIIIIIIIIIIIII".to_vec(),
        };
        assert!(passes_filter(&good_rec, &cfg, 33).is_ok());

        let short_rec = FastqRecord {
            id: b"@r".to_vec(),
            seq: b"ACGT".to_vec(),
            qual: b"IIII".to_vec(),
        };
        assert_eq!(passes_filter(&short_rec, &cfg, 33), Err("too_short"));
    }
}
