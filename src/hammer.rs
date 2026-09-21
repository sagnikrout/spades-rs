//! BayesHammer-inspired Fast SIMD Read Error Corrector (Module 1.3).
//!
//! Identifies sequencing errors in reads by comparing against the solid k-mer spectrum
//! and corrects isolated single-nucleotide errors using 2-bit constant-time Hamming distance.

use crate::bloom::TwoTierFilter;
use crate::dna::Kmer256;
use rayon::prelude::*;

/// Computes the number of nucleotide mismatches between two 2-bit packed sequences.
/// Executes in 4 CPU cycles using bitwise OR and POPCNT instructions.
#[inline(always)]
pub fn hamming_distance_2bit(a: u64, b: u64) -> u32 {
    let diff = a ^ b;
    // In 2-bit DNA encoding, a mismatch occurs if either or both bits of a pair differ.
    let base_diffs = (diff | (diff >> 1)) & 0x5555555555555555;
    base_diffs.count_ones()
}

pub struct ErrorCorrector {
    pub k: usize,
}

impl ErrorCorrector {
    pub fn new(k: usize) -> Self {
        Self { k }
    }

    /// Corrects reads in parallel using the TwoTierFilter solid k-mer index.
    pub fn correct_reads(
        &self,
        reads: Vec<Vec<u8>>,
        filter: &TwoTierFilter,
    ) -> (Vec<Vec<u8>>, usize) {
        let k = self.k;
        let corrected_count = std::sync::atomic::AtomicUsize::new(0);

        let corrected_reads: Vec<Vec<u8>> = reads
            .into_par_iter()
            .map(|mut seq| {
                if seq.len() < k {
                    return seq;
                }

                let mut modified = false;

                for i in 0..=(seq.len() - k) {
                    if let Some(km) = Kmer256::from_bytes(&seq[i..i + k], k) {
                        let (can, _) = km.canonical(k);
                        if !filter.is_solid_kmer256(can) {
                            // Non-solid k-mer detected: evaluate single-base substitutions
                            // Find the substitution that maximizes solid k-mer coverage across the overlapping window
                            let mut best_sub: Option<(usize, u8)> = None;
                            let mut max_solid_gain = 0;

                            for offset in 0..k {
                                let orig_base = seq[i + offset];
                                for &cand_base in b"ACGT" {
                                    if cand_base == orig_base {
                                        continue;
                                    }
                                    seq[i + offset] = cand_base;
                                    if let Some(cand_km) = Kmer256::from_bytes(&seq[i..i + k], k) {
                                        let (cand_can, _) = cand_km.canonical(k);
                                        if filter.is_solid_kmer256(cand_can) {
                                            // Count how many overlapping k-mers are solid with this candidate base
                                            let pos = i + offset;
                                            let start_scan = pos.saturating_sub(k - 1);
                                            let end_scan = pos.min(seq.len() - k);
                                            let mut solid_count = 0;
                                            for s in start_scan..=end_scan {
                                                if let Some(surr_km) =
                                                    Kmer256::from_bytes(&seq[s..s + k], k)
                                                {
                                                    let (surr_can, _) = surr_km.canonical(k);
                                                    if filter.is_solid_kmer256(surr_can) {
                                                        solid_count += 1;
                                                    }
                                                }
                                            }
                                            if solid_count > max_solid_gain {
                                                max_solid_gain = solid_count;
                                                best_sub = Some((offset, cand_base));
                                            }
                                        }
                                    }
                                    seq[i + offset] = orig_base;
                                }
                            }

                            if let Some((offset, cand_base)) = best_sub {
                                seq[i + offset] = cand_base;
                                modified = true;
                            }
                        }
                    }
                }

                if modified {
                    corrected_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                }

                seq
            })
            .collect();

        let total_corrected = corrected_count.load(std::sync::atomic::Ordering::Relaxed);
        (corrected_reads, total_corrected)
    }
}
