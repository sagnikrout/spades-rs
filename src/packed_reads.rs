//! Contiguous 2-bit packed read store for high-scale memory optimization.
//!
//! Packs 4 DNA bases per byte (A=00, C=01, G=10, T=11).
//! Stores millions of reads in a single contiguous buffer to eliminate
//! multi-gigabyte heap fragmentation and glibc arena bloat.
//!
//! # Three-Block Layout & N-Base Ambiguity Model
//!
//! Reads are ingested through a three-gate pipeline and organized into a deterministic
//! three-block contiguous layout:
//!
//! ```text
//! ┌───────────────────────────────┬───────────────────────────────┬─────────────────────────┐
//! │ BLOCK 1: Clean R1 Reads       │ BLOCK 2: Clean R2 Reads       │ BLOCK 3:                │
//! │ Indices: 0 .. pe_count - 1    │ Indices: pe_count .. 2*pe-1   │ Unpaired & Sub-reads    │
//! │ (Intact, zero internal Ns)    │ (Intact, zero internal Ns)    │ (from split reads ≥21bp)│
//! └───────────────────────────────┴───────────────────────────────┴─────────────────────────┘
//! ```
//!
//! 1. **Gate 1 (Terminal Trim)**: Leading and trailing `N`/`n` bases are stripped first (zero allocation).
//!    If terminal trimming leaves a clean sequence ≥ `MIN_K` (21 bp) without internal `N` breaks
//!    in both mates, the pair is retained in Blocks 1 and 2 with $O(1)$ mate lookup: `mate(p) = p + pe_count`.
//! 2. **Gate 2 (Internal N-Split)**: Reads containing internal `N` breaks are split into contiguous
//!    ACGT sub-reads. Sub-reads ≥ `MIN_K` are routed to Block 3 so all valid genomic k-mers enter
//!    the Bloom filter and de Bruijn graph without desynchronizing paired-end indices.
//! 3. **Gate 3 (Ambiguous Sidecar)**: Reads whose sub-reads are all shorter than `MIN_K` are
//!    preserved in `ambiguous_sidecar` as raw ASCII for Stage 8 consensus base voting.

use crate::dna::{base_to_2bit, bit2_to_base};
use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Minimum k-mer size: sub-reads shorter than this cannot produce any valid k-mer.
const MIN_K: usize = 21;

/// Contiguous 2-bit in-memory read store with synchronized paired-end blocks and an ambiguous-base sidecar.
#[derive(Clone, Default, Debug)]
pub struct PackedReads {
    /// Contiguous packed 2-bit bases (4 bases per byte).
    pub data: Vec<u8>,
    /// Byte start offset of each read in `data`.
    pub offsets: Vec<usize>,
    /// Exact base length of each read.
    pub lengths: Vec<u32>,
    /// Number of synchronized paired-end read pairs (N_pairs).
    /// - Block 1 (R1): `0 .. pe_count`
    /// - Block 2 (R2): `pe_count .. 2 * pe_count`
    /// - Block 3 (Unpaired & split sub-reads): `2 * pe_count .. len()`
    pub pe_count: usize,
    /// Backwards-compatible alias for `pe_count`.
    pub pe_boundary: usize,
    /// Reads with internal Ns whose all ACGT sub-reads were shorter than MIN_K.
    /// Stored as trimmed raw ASCII. Used only for polishing coverage votes.
    pub ambiguous_sidecar: Vec<Vec<u8>>,
}

/// Streaming FASTQ/FASTA reader supporting plain text, multi-line FASTA, and `.gz` archives.
struct FastxReader {
    reader: BufReader<Box<dyn Read + Send>>,
    is_fastq: bool,
    eof: bool,
    line_buf: Vec<u8>,
}

impl FastxReader {
    fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path).with_context(|| format!("Failed to open file: {:?}", path))?;

        let is_gz = path.extension().is_some_and(|ext| ext == "gz");
        let reader: Box<dyn Read + Send> = if is_gz {
            Box::new(MultiGzDecoder::new(file))
        } else {
            Box::new(file)
        };

        let mut buf_reader = BufReader::with_capacity(1024 * 1024, reader);
        let mut line_buf = Vec::with_capacity(1024);

        // Skip any leading blank lines before inspecting the first record header
        let is_fastq;
        loop {
            let buf = buf_reader.fill_buf()?;
            if buf.is_empty() {
                return Ok(Self {
                    reader: buf_reader,
                    is_fastq: true,
                    eof: true,
                    line_buf,
                });
            }
            if buf[0] == b'\n' || buf[0] == b'\r' {
                buf_reader.read_until(b'\n', &mut line_buf)?;
                line_buf.clear();
                continue;
            }
            is_fastq = buf[0] == b'@';
            break;
        }

        // For FASTA, consume the initial '>header' line so next_seq reads sequence lines until the next '>'.
        // For FASTQ, leave the '@header' line in the stream so next_seq validates all 4 lines per record.
        if !is_fastq {
            buf_reader.read_until(b'\n', &mut line_buf)?;
            line_buf.clear();
        }

        Ok(Self {
            reader: buf_reader,
            is_fastq,
            eof: false,
            line_buf,
        })
    }

    fn next_seq(&mut self, seq_out: &mut Vec<u8>) -> Result<bool> {
        seq_out.clear();
        if self.eof {
            return Ok(false);
        }

        if self.is_fastq {
            loop {
                // Line 1: '@header' (skip blank lines between records)
                self.line_buf.clear();
                if self.reader.read_until(b'\n', &mut self.line_buf)? == 0 {
                    self.eof = true;
                    return Ok(false);
                }
                while self.line_buf.ends_with(b"\n") || self.line_buf.ends_with(b"\r") {
                    self.line_buf.pop();
                }
                if self.line_buf.is_empty() {
                    continue;
                }
                if !self.line_buf.starts_with(b"@") {
                    anyhow::bail!(
                        "Malformed FASTQ record: expected '@' header, found '{}'",
                        String::from_utf8_lossy(&self.line_buf[..self.line_buf.len().min(32)])
                    );
                }

                // Line 2: Sequence
                self.line_buf.clear();
                if self.reader.read_until(b'\n', &mut self.line_buf)? == 0 {
                    self.eof = true;
                    anyhow::bail!(
                        "Malformed FASTQ record: unexpected EOF waiting for sequence line"
                    );
                }
                while self.line_buf.ends_with(b"\n") || self.line_buf.ends_with(b"\r") {
                    self.line_buf.pop();
                }
                seq_out.extend_from_slice(&self.line_buf);

                // Line 3: '+' separator
                self.line_buf.clear();
                if self.reader.read_until(b'\n', &mut self.line_buf)? == 0 {
                    self.eof = true;
                    anyhow::bail!(
                        "Malformed FASTQ record: unexpected EOF waiting for '+' separator"
                    );
                }
                if !self.line_buf.starts_with(b"+") {
                    anyhow::bail!(
                        "Malformed FASTQ record: expected '+' separator line, found '{}'",
                        String::from_utf8_lossy(&self.line_buf[..self.line_buf.len().min(32)])
                    );
                }

                // Line 4: Quality scores
                self.line_buf.clear();
                if self.reader.read_until(b'\n', &mut self.line_buf)? == 0 {
                    self.eof = true;
                    anyhow::bail!(
                        "Malformed FASTQ record: unexpected EOF waiting for quality scores"
                    );
                }
                while self.line_buf.ends_with(b"\n") || self.line_buf.ends_with(b"\r") {
                    self.line_buf.pop();
                }
                if self.line_buf.len() != seq_out.len() {
                    anyhow::bail!(
                        "Malformed FASTQ record: sequence length ({}) != quality length ({})",
                        seq_out.len(),
                        self.line_buf.len()
                    );
                }

                return Ok(true);
            }
        } else {
            // FASTA (single-line or multi-line wrapped)
            loop {
                self.line_buf.clear();
                if self.reader.read_until(b'\n', &mut self.line_buf)? == 0 {
                    self.eof = true;
                    return Ok(!seq_out.is_empty());
                }
                while self.line_buf.ends_with(b"\n") || self.line_buf.ends_with(b"\r") {
                    self.line_buf.pop();
                }
                if self.line_buf.starts_with(b">") {
                    if !seq_out.is_empty() {
                        return Ok(true);
                    }
                } else {
                    seq_out.extend_from_slice(&self.line_buf);
                }
            }
        }
    }
}

impl PackedReads {
    /// Creates an empty PackedReads store with preallocated capacity.
    pub fn with_capacity(num_reads: usize, total_bases: usize) -> Self {
        Self {
            data: Vec::with_capacity(total_bases.div_ceil(4)),
            offsets: Vec::with_capacity(num_reads),
            lengths: Vec::with_capacity(num_reads),
            pe_count: 0,
            pe_boundary: 0,
            ambiguous_sidecar: Vec::new(),
        }
    }

    /// Total number of packed reads stored across all 3 blocks.
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.lengths.len()
    }

    /// True if store contains zero reads.
    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.lengths.is_empty()
    }

    /// Number of synchronized clean paired-end read pairs ($N_{\text{pairs}}$).
    #[inline(always)]
    pub fn pe_count(&self) -> usize {
        self.pe_count.max(self.pe_boundary)
    }

    /// True if the store holds synchronized paired-end reads.
    #[inline(always)]
    pub fn is_paired(&self) -> bool {
        self.pe_count() > 0
    }

    /// Retrieves both mates of clean pair `p` ($0 \le p < \text{pe\_count}$) in $O(1)$ time:
    /// `r1_idx = p`, `r2_idx = p + pe_count`.
    #[inline(always)]
    pub fn get_pair(&self, p: usize, buf1: &mut Vec<u8>, buf2: &mut Vec<u8>) {
        let count = self.pe_count();
        debug_assert!(
            p < count,
            "Pair index {} out of range (pe_count: {})",
            p,
            count
        );
        self.get_read(p, buf1);
        self.get_read(p + count, buf2);
    }

    /// Appends all packed reads from `other` into `self` with zero re-encoding.
    fn append_packed(&mut self, mut other: Self) {
        let base_offset = self.data.len();
        self.data.append(&mut other.data);
        self.offsets
            .extend(other.offsets.into_iter().map(|o| o + base_offset));
        self.lengths.append(&mut other.lengths);
        self.ambiguous_sidecar.append(&mut other.ambiguous_sidecar);
    }

    /// Ingests and packs FASTQ/FASTA files into the Three-Block layout.
    /// - If `files.len() >= 2`, `files[0]` and `files[1]` are ingested synchronously as paired-end
    ///   into Blocks 1 and 2 (with internal-N sub-reads routed to Block 3), and any additional
    ///   files (`files[2..]`) are appended to Block 3.
    /// - If `files.len() == 1`, the file is ingested into Block 3 (`pe_count = 0`).
    pub fn from_files<P: AsRef<Path>>(files: &[P]) -> Result<Self> {
        if files.is_empty() {
            return Ok(Self::default());
        }
        if files.len() >= 2 {
            let mut store = Self::from_paired_files(&files[0], &files[1])?;
            for extra in &files[2..] {
                store.append_from_file(extra)?;
            }
            Ok(store)
        } else {
            let mut store = Self::default();
            store.append_from_file(&files[0])?;
            Ok(store)
        }
    }

    /// Synchronously ingests forward (`r1_path`) and reverse (`r2_path`) paired-end files
    /// directly into 2-bit Three-Block buffers in a single streaming pass.
    pub fn from_paired_files<P1: AsRef<Path>, P2: AsRef<Path>>(
        r1_path: P1,
        r2_path: P2,
    ) -> Result<Self> {
        let mut r1_reader = FastxReader::open(&r1_path)?;
        let mut r2_reader = FastxReader::open(&r2_path)?;

        let mut b1_store = Self::with_capacity(65_536, 65_536 * 150);
        let mut b2_store = Self::with_capacity(65_536, 65_536 * 150);
        let mut b3_store = Self::default();

        let mut seq1 = Vec::with_capacity(512);
        let mut seq2 = Vec::with_capacity(512);
        let mut r1_records = 0usize;
        let mut r2_records = 0usize;

        loop {
            let has1 = r1_reader.next_seq(&mut seq1)?;
            let has2 = r2_reader.next_seq(&mut seq2)?;
            if has1 {
                r1_records += 1;
            }
            if has2 {
                r2_records += 1;
            }

            match (has1, has2) {
                (true, true) => {
                    partition_pair_streaming(
                        &seq1,
                        &seq2,
                        &mut b1_store,
                        &mut b2_store,
                        &mut b3_store,
                    );
                }
                (true, false) => {
                    b3_store.ingest_read(&seq1);
                }
                (false, true) => {
                    b3_store.ingest_read(&seq2);
                }
                (false, false) => break,
            }
        }

        if r1_records != r2_records {
            eprintln!(
                "  [PackedReads Warning] Paired files {:?} ({} reads) and {:?} ({} reads) have unequal record counts; trailing reads routed to Block 3.",
                r1_path.as_ref(),
                r1_records,
                r2_path.as_ref(),
                r2_records
            );
        }

        Ok(Self::finalize_three_blocks(b1_store, b2_store, b3_store))
    }

    /// Synchronously ingests an interleaved paired-end file (`--12`) into the Three-Block layout.
    pub fn from_interleaved_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut reader = FastxReader::open(path)?;

        let mut b1_store = Self::with_capacity(65_536, 65_536 * 150);
        let mut b2_store = Self::with_capacity(65_536, 65_536 * 150);
        let mut b3_store = Self::default();

        let mut seq1 = Vec::with_capacity(512);
        let mut seq2 = Vec::with_capacity(512);

        while reader.next_seq(&mut seq1)? {
            if reader.next_seq(&mut seq2)? {
                partition_pair_streaming(&seq1, &seq2, &mut b1_store, &mut b2_store, &mut b3_store);
            } else {
                // Trailing orphan read at end of interleaved stream
                b3_store.ingest_read(&seq1);
                break;
            }
        }

        Ok(Self::finalize_three_blocks(b1_store, b2_store, b3_store))
    }

    /// Synchronously ingests multiple interleaved paired-end files into the Three-Block layout.
    pub fn from_interleaved_files<P: AsRef<Path>>(paths: &[P]) -> Result<Self> {
        let mut b1_store = Self::with_capacity(65_536, 65_536 * 150);
        let mut b2_store = Self::with_capacity(65_536, 65_536 * 150);
        let mut b3_store = Self::default();

        let mut seq1 = Vec::with_capacity(512);
        let mut seq2 = Vec::with_capacity(512);

        for path in paths {
            let mut reader = FastxReader::open(path)?;
            while reader.next_seq(&mut seq1)? {
                if reader.next_seq(&mut seq2)? {
                    partition_pair_streaming(
                        &seq1,
                        &seq2,
                        &mut b1_store,
                        &mut b2_store,
                        &mut b3_store,
                    );
                } else {
                    b3_store.ingest_read(&seq1);
                    break;
                }
            }
        }

        Ok(Self::finalize_three_blocks(b1_store, b2_store, b3_store))
    }

    /// Constructs a Three-Block `PackedReads` store from in-memory paired read slices.
    pub fn from_paired_slices(r1_reads: &[Vec<u8>], r2_reads: &[Vec<u8>]) -> Self {
        let mut b1_store = Self::with_capacity(r1_reads.len(), r1_reads.len() * 150);
        let mut b2_store = Self::with_capacity(r2_reads.len(), r2_reads.len() * 150);
        let mut b3_store = Self::default();

        let common = r1_reads.len().min(r2_reads.len());
        for i in 0..common {
            partition_pair_streaming(
                &r1_reads[i],
                &r2_reads[i],
                &mut b1_store,
                &mut b2_store,
                &mut b3_store,
            );
        }
        for r1 in &r1_reads[common..] {
            b3_store.ingest_read(r1);
        }
        for r2 in &r2_reads[common..] {
            b3_store.ingest_read(r2);
        }

        Self::finalize_three_blocks(b1_store, b2_store, b3_store)
    }

    /// Constructs a Three-Block `PackedReads` store from an in-memory interleaved slice.
    pub fn from_interleaved_slices(reads: &[Vec<u8>]) -> Self {
        let mut b1_store = Self::with_capacity(reads.len() / 2, (reads.len() / 2) * 150);
        let mut b2_store = Self::with_capacity(reads.len() / 2, (reads.len() / 2) * 150);
        let mut b3_store = Self::default();

        let mut iter = reads.iter();
        while let Some(r1) = iter.next() {
            if let Some(r2) = iter.next() {
                partition_pair_streaming(r1, r2, &mut b1_store, &mut b2_store, &mut b3_store);
            } else {
                b3_store.ingest_read(r1);
            }
        }

        Self::finalize_three_blocks(b1_store, b2_store, b3_store)
    }

    fn finalize_three_blocks(mut b1_store: Self, b2_store: Self, b3_store: Self) -> Self {
        debug_assert_eq!(b1_store.len(), b2_store.len());
        let pe_count = b1_store.len();
        b1_store.pe_count = pe_count;
        b1_store.pe_boundary = pe_count;
        b1_store.append_packed(b2_store);
        b1_store.append_packed(b3_store);
        b1_store
    }

    /// Streams unpaired reads directly from a FASTQ/FASTA file into Block 3 of the 2-bit store.
    pub fn append_from_file<P: AsRef<Path>>(&mut self, path: P) -> Result<usize> {
        let mut reader = FastxReader::open(path)?;
        let mut seq = Vec::with_capacity(512);
        let mut count = 0;

        while reader.next_seq(&mut seq)? {
            self.ingest_read(&seq);
            count += 1;
        }

        Ok(count)
    }

    // -------------------------------------------------------------------------
    // Three-Gate Single/Unpaired Ingestion
    // -------------------------------------------------------------------------

    /// Ingests a raw unpaired read through the three-gate pipeline:
    /// 1. End-trim leading/trailing Ns (zero allocation — slice bounds only).
    /// 2. Split on internal Ns; pack each ACGT sub-read ≥ MIN_K into the primary store.
    /// 3. If no sub-read was long enough, push the trimmed read into the ambiguous sidecar.
    pub fn ingest_read(&mut self, read: &[u8]) {
        // Gate 1: End trim.
        let trimmed = trim_n_ends(read);
        self.ingest_trimmed_unpaired(trimmed);
    }

    /// Routes an already terminal-trimmed read through Gate 2 (internal N-split) and Gate 3 (sidecar).
    #[inline]
    fn ingest_trimmed_unpaired(&mut self, trimmed: &[u8]) {
        if trimmed.is_empty() {
            return;
        }

        // Gate 2: Split on internal Ns and pack qualifying sub-reads.
        let mut any_packed = false;
        for segment in split_on_n(trimmed) {
            if segment.len() >= MIN_K {
                self.pack_clean(segment);
                any_packed = true;
            }
            // Sub-reads < MIN_K cannot produce any k-mer — silently discard.
        }

        // Gate 3: Entire read had no long-enough sub-read → sidecar.
        if !any_packed {
            self.ambiguous_sidecar.push(trimmed.to_vec());
        }
    }

    /// Packs and appends a clean ACGT read into the 2-bit buffer in a single pass.
    ///
    /// Precondition: every byte in `read` is A, C, G, or T (upper or lowercase).
    /// In debug builds a non-ACGT byte triggers a panic assertion.
    /// In release builds a non-ACGT byte is skipped (no `unwrap_or(0)` mutation).
    pub fn pack_clean(&mut self, read: &[u8]) {
        let offset = self.data.len();
        self.offsets.push(offset);

        let mut clean_len = 0u32;
        let mut byte = 0u8;
        let mut shift = 6i32;

        for &b in read {
            let b2 = match base_to_2bit(b) {
                Some(v) => v,
                None => {
                    debug_assert!(false, "pack_clean received non-ACGT byte: {}", b as char);
                    continue; // release: skip instead of corrupting
                }
            };
            clean_len += 1;
            byte |= b2 << shift;
            if shift == 0 {
                self.data.push(byte);
                byte = 0;
                shift = 6;
            } else {
                shift -= 2;
            }
        }

        if shift != 6 {
            self.data.push(byte);
        }
        self.lengths.push(clean_len);
    }

    /// Backwards-compatible entry point for callers that already hold clean ACGT sequences.
    /// Routes directly to `pack_clean`.
    pub fn add_read(&mut self, read: &[u8]) {
        self.pack_clean(read);
    }

    // -------------------------------------------------------------------------
    // Read Decoding
    // -------------------------------------------------------------------------

    /// Decodes read `idx` into the caller's reusable scratch buffer without heap reallocation.
    #[inline(always)]
    pub fn get_read(&self, idx: usize, buf: &mut Vec<u8>) {
        buf.clear();
        let len = self.lengths[idx] as usize;
        if len == 0 {
            return;
        }
        buf.reserve(len);

        let start_byte = self.offsets[idx];
        let num_bytes = len.div_ceil(4);
        let slice = &self.data[start_byte..start_byte + num_bytes];
        let mut rem = len;

        for &byte in slice {
            let take = rem.min(4);
            for s in (4 - take..4).rev() {
                let code = (byte >> (s * 2)) & 3;
                buf.push(bit2_to_base(code));
            }
            rem -= take;
        }
    }

    /// Total allocated memory footprint of packed store in bytes.
    pub fn memory_usage_bytes(&self) -> usize {
        self.data.capacity()
            + self.offsets.capacity() * std::mem::size_of::<usize>()
            + self.lengths.capacity() * std::mem::size_of::<u32>()
            + self.ambiguous_sidecar.capacity() * std::mem::size_of::<Vec<u8>>()
            + self
                .ambiguous_sidecar
                .iter()
                .map(|v| v.capacity())
                .sum::<usize>()
    }
}

// -------------------------------------------------------------------------
// N-Handling & Partitioning Helpers (module-private)
// -------------------------------------------------------------------------

/// Strips leading and trailing `N`/`n` (or non-ACGT ambiguous) bases from a slice — zero allocation.
#[inline]
fn trim_n_ends(read: &[u8]) -> &[u8] {
    let is_ambig = |b: &u8| base_to_2bit(*b).is_none();
    let start = read.iter().position(|b| !is_ambig(b)).unwrap_or(read.len());
    let end = read
        .iter()
        .rposition(|b| !is_ambig(b))
        .map(|i| i + 1)
        .unwrap_or(0);
    if start >= end {
        &[]
    } else {
        &read[start..end]
    }
}

/// Checks if a terminal-trimmed sequence contains any internal `N`/`n` (or non-ACGT) bases.
#[inline]
fn has_internal_n(trimmed: &[u8]) -> bool {
    trimmed.iter().any(|&b| base_to_2bit(b).is_none())
}

/// Iterates over contiguous ACGT sub-slices by splitting at every `N`/`n` (or non-ACGT base).
/// Zero allocation — yields non-empty sub-slices of the original `read`.
#[inline]
fn split_on_n(read: &[u8]) -> impl Iterator<Item = &[u8]> {
    read.split(|&b| base_to_2bit(b).is_none())
        .filter(|s| !s.is_empty())
}

/// Routes a paired-end read pair directly into the streaming 2-bit stores:
/// - Clean pairs (both ≥ MIN_K, zero internal Ns after terminal trim) -> `b1_store` & `b2_store`.
/// - Broken/internal-N pairs -> `b3_store` (sub-reads ≥ MIN_K and sidecar < MIN_K).
#[inline]
fn partition_pair_streaming(
    raw1: &[u8],
    raw2: &[u8],
    b1_store: &mut PackedReads,
    b2_store: &mut PackedReads,
    b3_store: &mut PackedReads,
) {
    // Gate 1: Terminal trimming of poly-N ends BEFORE checking for internal breaks.
    let t1 = trim_n_ends(raw1);
    let t2 = trim_n_ends(raw2);

    // Gate 2: Check whether remainder contains internal Ns.
    let int_n1 = has_internal_n(t1);
    let int_n2 = has_internal_n(t2);

    if !int_n1 && !int_n2 && t1.len() >= MIN_K && t2.len() >= MIN_K {
        b1_store.pack_clean(t1);
        b2_store.pack_clean(t2);
    } else {
        b3_store.ingest_trimmed_unpaired(t1);
        b3_store.ingest_trimmed_unpaired(t2);
    }
}

// -------------------------------------------------------------------------
// Tests
// -------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_packed_reads_roundtrip() {
        let original_reads = vec![
            b"ACGTACGT".to_vec(),
            b"TTTTAAAACCCCGGGG".to_vec(),
            b"A".to_vec(),
            b"CGA".to_vec(),
            b"GATTACA".to_vec(),
        ];

        let mut store = PackedReads::default();
        for r in &original_reads {
            store.add_read(r);
        }

        assert_eq!(store.len(), 5);

        let mut buf = Vec::new();
        for (i, orig) in original_reads.iter().enumerate() {
            store.get_read(i, &mut buf);
            assert_eq!(&buf, orig, "Mismatch at read index {}", i);
        }
    }

    #[test]
    fn test_trim_n_ends_leading() {
        assert_eq!(trim_n_ends(b"NNNACGT"), b"ACGT");
    }

    #[test]
    fn test_trim_n_ends_trailing() {
        assert_eq!(trim_n_ends(b"ACGTNNN"), b"ACGT");
    }

    #[test]
    fn test_trim_n_ends_both() {
        assert_eq!(trim_n_ends(b"NNACGTNN"), b"ACGT");
    }

    #[test]
    fn test_trim_n_ends_all_n() {
        assert_eq!(trim_n_ends(b"NNNNN"), b"");
    }

    #[test]
    fn test_trim_n_ends_no_n() {
        assert_eq!(trim_n_ends(b"ACGTACGT"), b"ACGTACGT");
    }

    #[test]
    fn test_split_on_n_basic() {
        let segs: Vec<&[u8]> = split_on_n(b"ACGTNNNACGT").collect();
        assert_eq!(segs, vec![b"ACGT" as &[u8], b"ACGT"]);
    }

    #[test]
    fn test_split_on_n_no_n() {
        let segs: Vec<&[u8]> = split_on_n(b"ACGTACGT").collect();
        assert_eq!(segs, vec![b"ACGTACGT" as &[u8]]);
    }

    #[test]
    fn test_split_on_n_multiple() {
        let segs: Vec<&[u8]> = split_on_n(b"ACGTNACGTNACGT").collect();
        assert_eq!(segs, vec![b"ACGT" as &[u8], b"ACGT", b"ACGT"]);
    }

    #[test]
    fn test_ingest_n_tail_packs_trimmed_clean() {
        // 21 clean ACGT bases + 3 trailing Ns → trimmed to 21 bp, packed to primary store.
        let read = b"ACGTACGTACGTACGTACGTANNN";
        let mut store = PackedReads::default();
        store.ingest_read(read);

        assert_eq!(store.len(), 1, "Trimmed read must enter primary store");
        assert!(store.ambiguous_sidecar.is_empty());

        let mut buf = Vec::new();
        store.get_read(0, &mut buf);
        assert_eq!(buf, b"ACGTACGTACGTACGTACGTA");
        assert!(!buf.contains(&b'N'));
    }

    #[test]
    fn test_ingest_internal_n_splits_into_two_primary_reads() {
        // 21 bp + internal N + 21 bp → two reads in primary store, sidecar empty.
        let read = b"ACGTACGTACGTACGTACGTANNNACGTACGTACGTACGTACGTA";
        let mut store = PackedReads::default();
        store.ingest_read(read);

        assert_eq!(store.len(), 2);
        assert!(store.ambiguous_sidecar.is_empty());
    }

    #[test]
    fn test_ingest_all_n_discarded() {
        let mut store = PackedReads::default();
        store.ingest_read(b"NNNNNN");
        assert_eq!(store.len(), 0);
        assert!(store.ambiguous_sidecar.is_empty());
    }

    #[test]
    fn test_ingest_sub_reads_all_short_goes_to_sidecar() {
        // "ACGT" (4 bp) + N + "CGTA" (4 bp) — neither half ≥ 21 bp → sidecar.
        let mut store = PackedReads::default();
        store.ingest_read(b"ACGTNNCGTA");

        assert_eq!(
            store.len(),
            0,
            "No sub-read ≥ 21 bp; primary store must be empty"
        );
        assert_eq!(
            store.ambiguous_sidecar.len(),
            1,
            "Short-sub-read read must land in sidecar"
        );
    }

    #[test]
    fn test_no_artificial_a_from_n() {
        // The original bug: N → unwrap_or(0) → 'A'. Confirm this no longer happens.
        let mut store = PackedReads::default();
        store.ingest_read(b"ACGTNACGT");

        assert_eq!(store.len(), 0);
        assert!(!store.ambiguous_sidecar.is_empty());
        assert!(
            store.ambiguous_sidecar[0].contains(&b'N'),
            "Sidecar must preserve N, not mutate it to A"
        );
    }

    #[test]
    fn test_sidecar_read_not_in_primary_kmer_stream() {
        let mut store = PackedReads::default();
        store.ingest_read(b"ACGTNACGT");

        assert_eq!(
            store.len(),
            0,
            "Sidecar reads must not enter primary 2-bit store"
        );
    }

    #[test]
    fn test_three_block_paired_synchronization() {
        // Pair 0: Clean (24 bp)
        let r1_0 = b"AAAAACCCCCGGGGGTTTTTAAAA".to_vec();
        let r2_0 = b"TTTTTAAAAAGGGGGCCCCCTTTT".to_vec();

        // Pair 1: R1 has internal N splitting into TWO 22-bp sub-reads; R2 is clean (24 bp)
        // If not handled by Three-Block layout, this split would shift Pair 2!
        let r1_1 = b"ACGTACGTACGTACGTACGTACNNACGTACGTACGTACGTACGTAC".to_vec();
        let r2_1 = b"GGGGGGGGGGGGCCCCCCCCCCCC".to_vec();

        // Pair 2: Terminal Ns on R1 and R2, but clean 24-bp interior -> MUST stay in Blocks 1 & 2!
        let r1_2 = b"NNNCCCCCGGGGTTTTTAAAAACCCCCNN".to_vec();
        let r2_2 = b"NNNGGGGGCCCCTTTTTAAAAAGGGGGNN".to_vec();

        let store = PackedReads::from_paired_slices(
            &[r1_0.clone(), r1_1, r1_2],
            &[r2_0.clone(), r2_1.clone(), r2_2],
        );

        // Clean pairs: Pair 0 and Pair 2 -> pe_count = 2
        assert_eq!(store.pe_count(), 2);
        // Block 1 (2) + Block 2 (2) + Block 3 (2 sub-reads from r1_1 + 1 from r2_1 = 3) = 7 total reads
        assert_eq!(store.len(), 7);

        let mut b1 = Vec::new();
        let mut b2 = Vec::new();

        // Pair p=0 must be r1_0 and r2_0
        store.get_pair(0, &mut b1, &mut b2);
        assert_eq!(b1, r1_0);
        assert_eq!(b2, r2_0);

        // Pair p=1 must be trimmed r1_2 and r2_2 (zero index drift from Pair 1's split!)
        store.get_pair(1, &mut b1, &mut b2);
        assert_eq!(b1, b"CCCCCGGGGTTTTTAAAAACCCCC");
        assert_eq!(b2, b"GGGGGCCCCTTTTTAAAAAGGGGG");

        // Block 3 (indices 4, 5, 6) holds the sub-reads from broken Pair 1
        store.get_read(4, &mut b1);
        assert_eq!(b1, b"ACGTACGTACGTACGTACGTAC");
        store.get_read(5, &mut b1);
        assert_eq!(b1, b"ACGTACGTACGTACGTACGTAC");
        store.get_read(6, &mut b1);
        assert_eq!(b1, r2_1);
    }
}
