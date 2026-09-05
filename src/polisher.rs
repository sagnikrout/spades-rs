//! Mismatch Corrector & Base Polisher (Module 4.2).
//!
//! Re-aligns reads against assembled contigs to calculate position-wise consensus,
//! correcting residual substitution errors and homopolymer indels.

use crate::dna::{canonical_kmer_u64, string_to_kmer};
use crate::graph::Unitig;
use hashbrown::HashMap;
use rayon::prelude::*;
use std::sync::atomic::{AtomicU32, Ordering};

pub struct Polisher {
    pub k: usize,
    pub min_coverage_support: u32,
}

impl Default for Polisher {
    fn default() -> Self {
        Self {
            k: 21,
            min_coverage_support: 5,
        }
    }
}

impl Polisher {
    /// Polishes contigs by consensus base-calling.
    pub fn polish_contigs(
        &self,
        mut contigs: Vec<Unitig>,
        reads: &[Vec<u8>],
    ) -> (Vec<Unitig>, usize) {
        let k = self.k;
        if contigs.is_empty() || reads.is_empty() {
            return (contigs, 0);
        }

        // 1. Build a global k-mer lookup table for all contigs:
        // kmer -> (contig_idx, offset, is_rc)
        let total_kmers: usize = contigs.iter().map(|c| c.sequence.len().saturating_sub(k - 1)).sum();
        let mut kmer_pos: HashMap<u64, (u32, u32, bool)> = HashMap::with_capacity(total_kmers);

        for (c_idx, contig) in contigs.iter().enumerate() {
            let len = contig.sequence.len();
            if len < k {
                continue;
            }
            for i in 0..=(len - k) {
                if let Some(km) = string_to_kmer(&contig.sequence[i..i + k], k) {
                    let (can, is_rc) = canonical_kmer_u64(km, k);
                    kmer_pos.entry(can).or_insert((c_idx as u32, i as u32, is_rc));
                }
            }
        }

        // 2. Prepare atomic base tally structures for all contigs
        let tallies: Vec<Vec<[AtomicU32; 4]>> = contigs
            .iter()
            .map(|c| {
                (0..c.sequence.len())
                    .map(|_| [
                        AtomicU32::new(0),
                        AtomicU32::new(0),
                        AtomicU32::new(0),
                        AtomicU32::new(0),
                    ])
                    .collect()
            })
            .collect();

        // 3. Process reads in parallel chunks
        reads.par_chunks(4096).for_each(|chunk| {
            for read in chunk {
                if read.len() < k {
                    continue;
                }
                // Try finding an anchor in the contigs
                let step = (k / 2 + 1).max(1);
                for r_idx in (0..=(read.len() - k)).step_by(step) {
                    if let Some(km) = string_to_kmer(&read[r_idx..r_idx + k], k) {
                        let (can, r_is_rc) = canonical_kmer_u64(km, k);
                        if let Some(&(c_idx, c_pos, c_is_rc)) = kmer_pos.get(&can) {
                            let c_idx = c_idx as usize;
                            let c_pos = c_pos as usize;
                            let c_len = tallies[c_idx].len();

                            if r_is_rc == c_is_rc {
                                // Align read to contig starting from anchor offset (forward)
                                let start_contig = c_pos as isize - r_idx as isize;
                                for (read_offset, &base) in read.iter().enumerate() {
                                    let curr_c_pos = start_contig + read_offset as isize;
                                    if curr_c_pos >= 0 && (curr_c_pos as usize) < c_len {
                                        if let Some(code) = crate::dna::base_to_2bit(base) {
                                            tallies[c_idx][curr_c_pos as usize][code as usize]
                                                .fetch_add(1, Ordering::Relaxed);
                                        }
                                    }
                                }
                            } else {
                                // Align read in reverse complement
                                let start_contig = c_pos as isize + k as isize - 1 + r_idx as isize;
                                for (read_offset, &base) in read.iter().enumerate() {
                                    let curr_c_pos = start_contig - read_offset as isize;
                                    if curr_c_pos >= 0 && (curr_c_pos as usize) < c_len {
                                        if let Some(code) = crate::dna::base_to_2bit(base) {
                                            let rc_code = (!code) & 3;
                                            tallies[c_idx][curr_c_pos as usize][rc_code as usize]
                                                .fetch_add(1, Ordering::Relaxed);
                                        }
                                    }
                                }
                            }
                            break;
                        }
                    }
                }
            }
        });

        // 4. Apply consensus where support >= threshold and > 80% dominant
        let mut total_corrections = 0;
        for (c_idx, contig) in contigs.iter_mut().enumerate() {
            let len = contig.sequence.len();
            for pos in 0..len {
                let counts = &tallies[c_idx][pos];
                let c0 = counts[0].load(Ordering::Relaxed);
                let c1 = counts[1].load(Ordering::Relaxed);
                let c2 = counts[2].load(Ordering::Relaxed);
                let c3 = counts[3].load(Ordering::Relaxed);
                let total = c0 + c1 + c2 + c3;

                if total >= self.min_coverage_support {
                    let ary = [c0, c1, c2, c3];
                    let mut max_code = 0;
                    let mut max_c = 0;
                    for (code, &c) in ary.iter().enumerate() {
                        if c > max_c {
                            max_c = c;
                            max_code = code;
                        }
                    }

                    if (max_c as f64 / total as f64) >= 0.8 {
                        let consensus_base = crate::dna::bit2_to_base(max_code as u8);
                        if contig.sequence[pos] != consensus_base {
                            contig.sequence[pos] = consensus_base;
                            total_corrections += 1;
                        }
                    }
                }
            }
        }

        (contigs, total_corrections)
    }

    /// Polishes contigs by consensus base-calling directly from PackedReads without intermediate vectors.
    pub fn polish_contigs_packed(
        &self,
        mut contigs: Vec<Unitig>,
        packed: &crate::packed_reads::PackedReads,
    ) -> (Vec<Unitig>, usize) {
        let k = self.k;
        if contigs.is_empty() || packed.is_empty() {
            return (contigs, 0);
        }

        let total_kmers: usize = contigs.iter().map(|c| c.sequence.len().saturating_sub(k - 1)).sum();
        let mut kmer_pos: HashMap<u64, (u32, u32, bool)> = HashMap::with_capacity(total_kmers);

        for (c_idx, contig) in contigs.iter().enumerate() {
            let len = contig.sequence.len();
            if len < k {
                continue;
            }
            for i in 0..=(len - k) {
                if let Some(km) = string_to_kmer(&contig.sequence[i..i + k], k) {
                    let (can, is_rc) = canonical_kmer_u64(km, k);
                    kmer_pos.entry(can).or_insert((c_idx as u32, i as u32, is_rc));
                }
            }
        }

        let tallies: Vec<Vec<[AtomicU32; 4]>> = contigs
            .iter()
            .map(|c| {
                (0..c.sequence.len())
                    .map(|_| [
                        AtomicU32::new(0),
                        AtomicU32::new(0),
                        AtomicU32::new(0),
                        AtomicU32::new(0),
                    ])
                    .collect()
            })
            .collect();

        let num_reads = packed.len();
        let indices: Vec<usize> = (0..num_reads).collect();

        indices.par_chunks(4096).for_each(|chunk| {
            let mut read = Vec::with_capacity(512);
            for &idx in chunk {
                packed.get_read(idx, &mut read);
                if read.len() < k {
                    continue;
                }
                let step = (k / 2 + 1).max(1);
                for r_idx in (0..=(read.len() - k)).step_by(step) {
                    if let Some(km) = string_to_kmer(&read[r_idx..r_idx + k], k) {
                        let (can, r_is_rc) = canonical_kmer_u64(km, k);
                        if let Some(&(c_idx, c_pos, c_is_rc)) = kmer_pos.get(&can) {
                            let c_idx = c_idx as usize;
                            let c_pos = c_pos as usize;
                            let c_len = tallies[c_idx].len();

                            if r_is_rc == c_is_rc {
                                let start_contig = c_pos as isize - r_idx as isize;
                                for (read_offset, &base) in read.iter().enumerate() {
                                    let curr_c_pos = start_contig + read_offset as isize;
                                    if curr_c_pos >= 0 && (curr_c_pos as usize) < c_len {
                                        if let Some(code) = crate::dna::base_to_2bit(base) {
                                            tallies[c_idx][curr_c_pos as usize][code as usize]
                                                .fetch_add(1, Ordering::Relaxed);
                                        }
                                    }
                                }
                            } else {
                                let start_contig = c_pos as isize + k as isize - 1 + r_idx as isize;
                                for (read_offset, &base) in read.iter().enumerate() {
                                    let curr_c_pos = start_contig - read_offset as isize;
                                    if curr_c_pos >= 0 && (curr_c_pos as usize) < c_len {
                                        if let Some(code) = crate::dna::base_to_2bit(base) {
                                            let rc_code = (!code) & 3;
                                            tallies[c_idx][curr_c_pos as usize][rc_code as usize]
                                                .fetch_add(1, Ordering::Relaxed);
                                        }
                                    }
                                }
                            }
                            break;
                        }
                    }
                }
            }
        });

        let mut total_corrections = 0;
        for (c_idx, contig) in contigs.iter_mut().enumerate() {
            let len = contig.sequence.len();
            for pos in 0..len {
                let counts = &tallies[c_idx][pos];
                let c0 = counts[0].load(Ordering::Relaxed);
                let c1 = counts[1].load(Ordering::Relaxed);
                let c2 = counts[2].load(Ordering::Relaxed);
                let c3 = counts[3].load(Ordering::Relaxed);
                let total = c0 + c1 + c2 + c3;

                if total >= self.min_coverage_support {
                    let ary = [c0, c1, c2, c3];
                    let mut max_code = 0;
                    let mut max_c = 0;
                    for (code, &c) in ary.iter().enumerate() {
                        if c > max_c {
                            max_c = c;
                            max_code = code;
                        }
                    }

                    if (max_c as f64 / total as f64) >= 0.8 {
                        let consensus_base = crate::dna::bit2_to_base(max_code as u8);
                        if contig.sequence[pos] != consensus_base {
                            contig.sequence[pos] = consensus_base;
                            total_corrections += 1;
                        }
                    }
                }
            }
        }

        (contigs, total_corrections)
    }
}
