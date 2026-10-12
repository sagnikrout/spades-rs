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

        let buf = buf_reader.fill_buf()?;
        if buf.is_empty() {
            return Ok(Self {
                reader: buf_reader,
                is_fastq: true,
                eof: true,
                line_buf,
            });
        }

        let is_fastq = buf[0] == b'@';
        // Consume initial header line ('@...' or '>...')
        buf_reader.read_until(b'\n', &mut line_buf)?;
        line_buf.clear();

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
                // Line 2: Sequence
                self.line_buf.clear();
                if self.reader.read_until(b'\n', &mut self.line_buf)? == 0 {
                    self.eof = true;
                    return Ok(false);
                }
                while self.line_buf.ends_with(b"\n") || self.line_buf.ends_with(b"\r") {
                    self.line_buf.pop();
                }
                seq_out.extend_from_slice(&self.line_buf);

                // Line 3: '+' separator
                self.line_buf.clear();
                if self.reader.read_until(b'\n', &mut self.line_buf)? == 0 {
                    self.eof = true;
                    return Ok(!seq_out.is_empty());
                }

                // Line 4: Quality scores
                self.line_buf.clear();
                if self.reader.read_until(b'\n', &mut self.line_buf)? == 0 {
                    self.eof = true;
                    return Ok(!seq_out.is_empty());
                }

                // Line 1 of next record: '@header'
                self.line_buf.clear();
                if self.reader.read_until(b'\n', &mut self.line_buf)? == 0 {
                    self.eof = true;
                }

                if !seq_out.is_empty() {
                    return Ok(true);
                }
                if self.eof {
                    return Ok(false);
                }
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
        if self.pe_count > 0 {
            self.pe_count
        } else {
            self.pe_boundary
        }
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
    /// into the Three-Block layout.
    pub fn from_paired_files<P1: AsRef<Path>, P2: AsRef<Path>>(
        r1_path: P1,
        r2_path: P2,
    ) -> Result<Self> {
        let mut r1_reader = FastxReader::open(r1_path)?;
        let mut r2_reader = FastxReader::open(r2_path)?;

        let mut clean_r1: Vec<Vec<u8>> = Vec::with_capacity(65_536);
        let mut clean_r2: Vec<Vec<u8>> = Vec::with_capacity(65_536);
        let mut unpaired_subreads: Vec<Vec<u8>> = Vec::new();
        let mut sidecar: Vec<Vec<u8>> = Vec::new();

        let mut seq1 = Vec::with_capacity(512);
        let mut seq2 = Vec::with_capacity(512);

        loop {
            let has1 = r1_reader.next_seq(&mut seq1)?;
            let has2 = r2_reader.next_seq(&mut seq2)?;

            match (has1, has2) {
                (true, true) => {
                    partition_pair(
                        &seq1,
                        &seq2,
                        &mut clean_r1,
                        &mut clean_r2,
                        &mut unpaired_subreads,
                        &mut sidecar,
                    );
                }
                (true, false) => {
                    route_unpaired_or_split(
                        trim_n_ends(&seq1),
                        &mut unpaired_subreads,
                        &mut sidecar,
                    );
                }
                (false, true) => {
                    route_unpaired_or_split(
                        trim_n_ends(&seq2),
                        &mut unpaired_subreads,
                        &mut sidecar,
                    );
                }
                (false, false) => break,
            }
        }

        Ok(Self::build_three_block_store(
            clean_r1,
            clean_r2,
            unpaired_subreads,
            sidecar,
        ))
    }

    /// Synchronously ingests an interleaved paired-end file (`--12`) into the Three-Block layout.
    pub fn from_interleaved_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut reader = FastxReader::open(path)?;

        let mut clean_r1: Vec<Vec<u8>> = Vec::with_capacity(65_536);
        let mut clean_r2: Vec<Vec<u8>> = Vec::with_capacity(65_536);
        let mut unpaired_subreads: Vec<Vec<u8>> = Vec::new();
        let mut sidecar: Vec<Vec<u8>> = Vec::new();

        let mut seq1 = Vec::with_capacity(512);
        let mut seq2 = Vec::with_capacity(512);

        while reader.next_seq(&mut seq1)? {
            if reader.next_seq(&mut seq2)? {
                partition_pair(
                    &seq1,
                    &seq2,
                    &mut clean_r1,
                    &mut clean_r2,
                    &mut unpaired_subreads,
                    &mut sidecar,
                );
            } else {
                // Trailing orphan read at end of interleaved stream
                route_unpaired_or_split(trim_n_ends(&seq1), &mut unpaired_subreads, &mut sidecar);
                break;
            }
        }

        Ok(Self::build_three_block_store(
            clean_r1,
            clean_r2,
            unpaired_subreads,
            sidecar,
        ))
    }

    /// Constructs a Three-Block `PackedReads` store from in-memory paired read slices.
    pub fn from_paired_slices(r1_reads: &[Vec<u8>], r2_reads: &[Vec<u8>]) -> Self {
        let mut clean_r1 = Vec::with_capacity(r1_reads.len());
        let mut clean_r2 = Vec::with_capacity(r2_reads.len());
        let mut unpaired_subreads = Vec::new();
        let mut sidecar = Vec::new();

        let common = r1_reads.len().min(r2_reads.len());
        for i in 0..common {
            partition_pair(
                &r1_reads[i],
                &r2_reads[i],
                &mut clean_r1,
                &mut clean_r2,
                &mut unpaired_subreads,
                &mut sidecar,
            );
        }
        for r1 in &r1_reads[common..] {
            route_unpaired_or_split(trim_n_ends(r1), &mut unpaired_subreads, &mut sidecar);
        }
        for r2 in &r2_reads[common..] {
            route_unpaired_or_split(trim_n_ends(r2), &mut unpaired_subreads, &mut sidecar);
        }

        Self::build_three_block_store(clean_r1, clean_r2, unpaired_subreads, sidecar)
    }

    /// Constructs a Three-Block `PackedReads` store from an in-memory interleaved slice.
    pub fn from_interleaved_slices(reads: &[Vec<u8>]) -> Self {
        let mut clean_r1 = Vec::with_capacity(reads.len() / 2);
        let mut clean_r2 = Vec::with_capacity(reads.len() / 2);
        let mut unpaired_subreads = Vec::new();
        let mut sidecar = Vec::new();

        let mut iter = reads.iter();
        while let Some(r1) = iter.next() {
            if let Some(r2) = iter.next() {
                partition_pair(
                    r1,
                    r2,
                    &mut clean_r1,
                    &mut clean_r2,
                    &mut unpaired_subreads,
                    &mut sidecar,
                );
            } else {
                route_unpaired_or_split(trim_n_ends(r1), &mut unpaired_subreads, &mut sidecar);
            }
        }

        Self::build_three_block_store(clean_r1, clean_r2, unpaired_subreads, sidecar)
    }

    fn build_three_block_store(
        clean_r1: Vec<Vec<u8>>,
        clean_r2: Vec<Vec<u8>>,
        unpaired_subreads: Vec<Vec<u8>>,
        sidecar: Vec<Vec<u8>>,
    ) -> Self {
        debug_assert_eq!(clean_r1.len(), clean_r2.len());
        let pe_count = clean_r1.len();
        let total_reads = pe_count * 2 + unpaired_subreads.len();
        let total_bases: usize = clean_r1.iter().map(|r| r.len()).sum::<usize>()
            + clean_r2.iter().map(|r| r.len()).sum::<usize>()
            + unpaired_subreads.iter().map(|r| r.len()).sum::<usize>();

        let mut store = Self::with_capacity(total_reads, total_bases);
        store.ambiguous_sidecar = sidecar;

        // Block 1: Clean R1 [0 .. pe_count]
        for r1 in &clean_r1 {
            store.pack_clean(r1);
        }
        store.pe_count = pe_count;
        store.pe_boundary = pe_count;

        // Block 2: Clean R2 [pe_count .. 2 * pe_count]
        for r2 in &clean_r2 {
            store.pack_clean(r2);
        }

        // Block 3: Unpaired & split sub-reads [2 * pe_count .. total_reads]
        for sub in &unpaired_subreads {
            store.pack_clean(sub);
        }

        store
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
        let read = trim_n_ends(read);
        if read.is_empty() {
            return;
        }

        // Gate 2: Split on internal Ns and pack qualifying sub-reads.
        let mut any_packed = false;
        for segment in split_on_n(read) {
            if segment.len() >= MIN_K {
                self.pack_clean(segment);
                any_packed = true;
            }
            // Sub-reads < MIN_K cannot produce any k-mer — silently discard.
        }

        // Gate 3: Entire read had no long-enough sub-read → sidecar.
        if !any_packed {
            self.ambiguous_sidecar.push(read.to_vec());
        }
    }

    /// Packs and appends a clean ACGT read into the 2-bit buffer.
    ///
    /// Precondition: every byte in `read` is A, C, G, or T (upper or lowercase).
    /// In debug builds a non-ACGT byte triggers a panic assertion.
    /// In release builds a non-ACGT byte is skipped (no `unwrap_or(0)` mutation).
    pub fn pack_clean(&mut self, read: &[u8]) {
        let clean_len = read.iter().filter(|&&b| base_to_2bit(b).is_some()).count();
        let offset = self.data.len();
        self.offsets.push(offset);
        self.lengths.push(clean_len as u32);
        if clean_len == 0 {
            return;
        }

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

    /// Total memory footprint of packed store in bytes.
    pub fn memory_usage_bytes(&self) -> usize {
        self.data.capacity()
            + self.offsets.capacity() * std::mem::size_of::<usize>()
            + self.lengths.capacity() * std::mem::size_of::<u32>()
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

/// Strips leading and trailing `N`/`n` bases from a slice — zero allocation.
#[inline]
fn trim_n_ends(read: &[u8]) -> &[u8] {
    let is_n = |b: &u8| *b == b'N' || *b == b'n';
    let start = read.iter().position(|b| !is_n(b)).unwrap_or(read.len());
    let end = read
        .iter()
        .rposition(|b| !is_n(b))
        .map(|i| i + 1)
        .unwrap_or(0);
    if start >= end {
        &[]
    } else {
        &read[start..end]
    }
}

/// Checks if a terminal-trimmed sequence contains any internal `N`/`n` bases.
#[inline]
fn has_internal_n(trimmed: &[u8]) -> bool {
    trimmed.iter().any(|&b| b == b'N' || b == b'n')
}

/// Iterates over contiguous ACGT sub-slices by splitting at every `N`/`n`.
/// Zero allocation — yields sub-slices of the original `read`.
#[inline]
fn split_on_n(read: &[u8]) -> impl Iterator<Item = &[u8]> {
    read.split(|&b| b == b'N' || b == b'n')
        .filter(|s| !s.is_empty())
}

/// Routes a paired-end read pair into either the synchronized paired blocks (Blocks 1 & 2)
/// or demotes broken/internal-N reads into Block 3 (sub-reads ≥ MIN_K) and `sidecar` (< MIN_K).
#[inline]
fn partition_pair(
    raw1: &[u8],
    raw2: &[u8],
    clean_r1: &mut Vec<Vec<u8>>,
    clean_r2: &mut Vec<Vec<u8>>,
    unpaired_subreads: &mut Vec<Vec<u8>>,
    sidecar: &mut Vec<Vec<u8>>,
) {
    // Gate 1: Terminal trimming of poly-N ends BEFORE checking for internal breaks.
    let t1 = trim_n_ends(raw1);
    let t2 = trim_n_ends(raw2);

    // Gate 2: Check whether remainder contains internal Ns.
    let int_n1 = has_internal_n(t1);
    let int_n2 = has_internal_n(t2);

    if !int_n1 && !int_n2 && t1.len() >= MIN_K && t2.len() >= MIN_K {
        clean_r1.push(t1.to_vec());
        clean_r2.push(t2.to_vec());
    } else {
        route_unpaired_or_split(t1, unpaired_subreads, sidecar);
        route_unpaired_or_split(t2, unpaired_subreads, sidecar);
    }
}

/// Splits a trimmed read on internal Ns: segments ≥ MIN_K go to Block 3 (`unpaired_subreads`),
/// and if no segment is ≥ MIN_K, the read is preserved in `sidecar` for consensus polishing.
#[inline]
fn route_unpaired_or_split(
    trimmed: &[u8],
    unpaired_subreads: &mut Vec<Vec<u8>>,
    sidecar: &mut Vec<Vec<u8>>,
) {
    if trimmed.is_empty() {
        return;
    }
    let mut any_packed = false;
    for segment in split_on_n(trimmed) {
        if segment.len() >= MIN_K {
            unpaired_subreads.push(segment.to_vec());
            any_packed = true;
        }
    }
    if !any_packed {
        sidecar.push(trimmed.to_vec());
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
