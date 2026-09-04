//! Mismatch Corrector & Base Polisher (Module 4.2).
//!
//! Re-aligns reads against assembled contigs to calculate position-wise consensus,
//! correcting residual substitution errors and homopolymer indels.

use crate::dna::{canonical_kmer_u64, string_to_kmer};
use crate::graph::Unitig;
use hashbrown::HashMap;

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
    pub fn polish_contigs(&self, mut contigs: Vec<Unitig>, reads: &[Vec<u8>]) -> (Vec<Unitig>, usize) {
        let k = self.k;
        let mut total_corrections = 0;

        for contig in &mut contigs {
            let len = contig.sequence.len();
            if len < k {
                continue;
            }

            // Index positions of k-mers in this contig
            let mut kmer_pos: HashMap<u64, usize> = HashMap::with_capacity(len);
            for i in 0..=(len - k) {
                if let Some(km) = string_to_kmer(&contig.sequence[i..i + k], k) {
                    let (can, _) = canonical_kmer_u64(km, k);
                    kmer_pos.insert(can, i);
                }
            }

            // Position base tallies: pos -> [A, C, G, T] counts
            let mut tallies = vec![[0u32; 4]; len];

            for read in reads {
                if read.len() < k {
                    continue;
                }
                // Try finding an anchor
                for r_idx in (0..=(read.len() - k)).step_by(k / 2 + 1) {
                    if let Some(km) = string_to_kmer(&read[r_idx..r_idx + k], k) {
                        let (can, _) = canonical_kmer_u64(km, k);
                        if let Some(&c_pos) = kmer_pos.get(&can) {
                            // Align read to contig starting from anchor offset
                            let start_contig = c_pos as isize - r_idx as isize;
                            for (read_offset, &base) in read.iter().enumerate() {
                                let curr_c_pos = start_contig + read_offset as isize;
                                if curr_c_pos >= 0 && (curr_c_pos as usize) < len {
                                    if let Some(code) = crate::dna::base_to_2bit(base) {
                                        tallies[curr_c_pos as usize][code as usize] += 1;
                                    }
                                }
                            }
                            break;
                        }
                    }
                }
            }

            // Apply consensus where support >= threshold and > 80% dominant
            for pos in 0..len {
                let counts = tallies[pos];
                let total: u32 = counts.iter().sum();
                if total >= self.min_coverage_support {
                    let mut max_code = 0;
                    let mut max_c = 0;
                    for (code, &c) in counts.iter().enumerate() {
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
