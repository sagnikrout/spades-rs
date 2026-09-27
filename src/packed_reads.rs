//! Contiguous 2-bit packed read store for high-scale memory optimization.
//!
//! Packs 4 DNA bases per byte (A=00, C=01, G=10, T=11).
//! Stores millions of reads in a single contiguous buffer to eliminate
//! multi-gigabyte heap fragmentation and glibc arena bloat.
//!
//! # N-Base Handling
//!
//! Reads are ingested through a three-gate pipeline:
//! 1. **End Trim**: Leading and trailing `N`/`n` bases are stripped in-place (zero allocation).
//! 2. **N-Split**: Internal `N`/`n` runs split the read into contiguous ACGT sub-reads.
//!    Sub-reads ≥ `MIN_K` (21 bp) are packed into the primary 2-bit store.
//! 3. **Sidecar**: Reads whose every sub-read is shorter than `MIN_K` (extremely rare, ~0.5%)
//!    are kept in `ambiguous_sidecar` as raw ASCII for polishing coverage votes only.
//!    They never enter the Bloom filter or the de Bruijn graph.

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
    /// Index separating Read 1 and Read 2 (if paired).
    pub pe_boundary: usize,
    /// Reads with internal Ns whose all ACGT sub-reads were shorter than MIN_K.
    /// Stored as trimmed raw ASCII. Used only for polishing coverage votes.
    pub ambiguous_sidecar: Vec<Vec<u8>>,
}

impl PackedReads {
    /// Creates an empty PackedReads store with preallocated capacity.
    pub fn with_capacity(num_reads: usize, total_bases: usize) -> Self {
        Self {
            data: Vec::with_capacity(total_bases.div_ceil(4)),
            offsets: Vec::with_capacity(num_reads),
            lengths: Vec::with_capacity(num_reads),
            pe_boundary: 0,
            ambiguous_sidecar: Vec::new(),
        }
    }

    /// Number of packed reads stored.
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.lengths.len()
    }

    /// True if store contains zero reads.
    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.lengths.is_empty()
    }

    /// Ingests and packs FASTQ/FASTA files into contiguous memory with zero intermediate vector allocations.
    pub fn from_files<P: AsRef<Path>>(files: &[P]) -> Result<Self> {
        let mut store = Self::default();
        let mut r1_count = 0;

        for (idx, f) in files.iter().enumerate() {
            let count = store.append_from_file(f)?;
            if idx == 0 {
                r1_count = count;
            }
        }
        store.pe_boundary = r1_count;
        Ok(store)
    }

    /// Streams reads directly from a FASTQ/FASTA file (plain text or gzip) into the 2-bit packed store.
    pub fn append_from_file<P: AsRef<Path>>(&mut self, path: P) -> Result<usize> {
        let path = path.as_ref();
        let file = File::open(path).with_context(|| format!("Failed to open file: {:?}", path))?;

        let is_gz = path.extension().is_some_and(|ext| ext == "gz");
        let reader: Box<dyn Read + Send> = if is_gz {
            Box::new(MultiGzDecoder::new(file))
        } else {
            Box::new(file)
        };

        let mut buf_reader = BufReader::with_capacity(1024 * 1024, reader);
        let mut line = Vec::with_capacity(1024);
        let mut count = 0;

        // Peek at first byte to determine format: '@' for FASTQ, '>' for FASTA
        let mut first_byte = [0u8; 1];
        if buf_reader.read(&mut first_byte)? == 0 {
            return Ok(0); // empty file
        }
        let is_fastq = first_byte[0] == b'@';

        if is_fastq {
            // Read remainder of line 1 (header)
            buf_reader.read_until(b'\n', &mut line)?;
            line.clear();

            loop {
                // Line 2: Sequence
                if buf_reader.read_until(b'\n', &mut line)? == 0 {
                    break;
                }
                while line.ends_with(b"\n") || line.ends_with(b"\r") {
                    line.pop();
                }
                if !line.is_empty() {
                    self.ingest_read(&line);
                    count += 1;
                }
                line.clear();

                // Line 3: +
                if buf_reader.read_until(b'\n', &mut line)? == 0 {
                    break;
                }
                line.clear();

                // Line 4: Qual
                if buf_reader.read_until(b'\n', &mut line)? == 0 {
                    break;
                }
                line.clear();

                // Line 1 of next record (@ID)
                if buf_reader.read_until(b'\n', &mut line)? == 0 {
                    break;
                }
                line.clear();
            }
        } else {
            // FASTA format
            let mut seq_buf = Vec::with_capacity(4096);
            loop {
                line.clear();
                if buf_reader.read_until(b'\n', &mut line)? == 0 {
                    break;
                }
                while line.ends_with(b"\n") || line.ends_with(b"\r") {
                    line.pop();
                }
                if line.starts_with(b">") {
                    if !seq_buf.is_empty() {
                        self.ingest_read(&seq_buf);
                        count += 1;
                        seq_buf.clear();
                    }
                } else {
                    seq_buf.extend_from_slice(&line);
                }
            }
            if !seq_buf.is_empty() {
                self.ingest_read(&seq_buf);
                count += 1;
            }
        }

        Ok(count)
    }

    // -------------------------------------------------------------------------
    // Three-Gate Ingestion
    // -------------------------------------------------------------------------

    /// Ingests a raw read through the three-gate pipeline:
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
        let offset = self.data.len();
        self.offsets.push(offset);
        self.lengths.push(read.len() as u32);

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
// N-Handling Helpers (module-private, zero allocation)
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

/// Iterates over contiguous ACGT sub-slices by splitting at every `N`/`n`.
/// Zero allocation — yields sub-slices of the original `read`.
#[inline]
fn split_on_n(read: &[u8]) -> impl Iterator<Item = &[u8]> {
    read.split(|&b| b == b'N' || b == b'n')
        .filter(|s| !s.is_empty())
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
        // Both halves of "ACGTNACGT" are 4 bp each (<21), so the read goes to the sidecar.
        let mut store = PackedReads::default();
        store.ingest_read(b"ACGTNACGT");

        // Primary 2-bit store is empty — no corruption possible.
        assert_eq!(store.len(), 0);

        // Sidecar must preserve the raw bytes including N — not substitute A.
        assert!(!store.ambiguous_sidecar.is_empty());
        assert!(
            store.ambiguous_sidecar[0].contains(&b'N'),
            "Sidecar must preserve N, not mutate it to A"
        );
    }

    #[test]
    fn test_sidecar_read_not_in_primary_kmer_stream() {
        // Reads in the sidecar must not appear in the primary index.
        let mut store = PackedReads::default();
        store.ingest_read(b"ACGTNACGT"); // both halves < 21 bp → sidecar

        assert_eq!(
            store.len(),
            0,
            "Sidecar reads must not enter primary 2-bit store"
        );
    }
}
