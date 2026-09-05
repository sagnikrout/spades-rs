//! Contiguous 2-bit packed read store for high-scale memory optimization.
//!
//! Packs 4 DNA bases per byte (A=00, C=01, G=10, T=11). 
//! Stores millions of reads in a single contiguous buffer to eliminate 
//! multi-gigabyte heap fragmentation and glibc arena bloat.

use crate::dna::{base_to_2bit, bit2_to_base};
use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

#[derive(Clone, Default, Debug)]
pub struct PackedReads {
    /// Contiguous packed 2-bit bases (4 bases per byte)
    pub data: Vec<u8>,
    /// Byte start offset of each read in `data`
    pub offsets: Vec<usize>,
    /// Exact base length of each read
    pub lengths: Vec<u32>,
    /// Index separating Read 1 and Read 2 (if paired)
    pub pe_boundary: usize,
}

impl PackedReads {
    /// Creates an empty PackedReads store with preallocated capacity.
    pub fn with_capacity(num_reads: usize, total_bases: usize) -> Self {
        Self {
            data: Vec::with_capacity((total_bases + 3) / 4),
            offsets: Vec::with_capacity(num_reads),
            lengths: Vec::with_capacity(num_reads),
            pe_boundary: 0,
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
                    self.add_read(&line);
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
                        self.add_read(&seq_buf);
                        count += 1;
                        seq_buf.clear();
                    }
                } else {
                    seq_buf.extend_from_slice(&line);
                }
            }
            if !seq_buf.is_empty() {
                self.add_read(&seq_buf);
                count += 1;
            }
        }

        Ok(count)
    }

    /// Packs and appends a single read into the contiguous byte buffer.
    pub fn add_read(&mut self, read: &[u8]) {
        let offset = self.data.len();
        self.offsets.push(offset);
        self.lengths.push(read.len() as u32);

        let mut byte = 0u8;
        let mut shift = 6i32;

        for &b in read {
            let b2 = base_to_2bit(b).unwrap_or(0); // non-ACGT default to 00
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
        let num_bytes = (len + 3) / 4;
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
    }
}

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
}
