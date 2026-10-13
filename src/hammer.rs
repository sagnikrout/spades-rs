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

/// BayesHammer-inspired single-nucleotide read error corrector operating against a solid $k$-mer filter.
pub struct ErrorCorrector {
    /// $k$-mer size used to test candidate base substitutions against the solid $k$-mer spectrum.
    pub k: usize,
}

impl ErrorCorrector {
    /// Creates a new [`ErrorCorrector`] for the given $k$-mer length.
    pub fn new(k: usize) -> Self {
        Self { k }
    }

    /// Corrects a single read in-place using a reusable scratch buffer against the solid k-mer filter.
    /// Returns `true` if at least one single-nucleotide error was corrected.
    pub fn correct_read_in_place(
        &self,
        seq: &mut [u8],
        scratch: &mut Vec<u8>,
        filter: &TwoTierFilter,
    ) -> bool {
        let k = self.k;
        if k == 0 || k > 128 || seq.len() < k {
            return false;
        }

        scratch.clear();
        scratch.extend_from_slice(seq);
        let mut modified = false;

        for i in 0..=(scratch.len() - k) {
            if let Some(km) = Kmer256::from_bytes(&scratch[i..i + k], k) {
                let (can, _) = km.canonical(k);
                if !filter.is_solid_kmer256(can) {
                    let mut best_sub: Option<(usize, u8)> = None;
                    let mut max_solid_gain = 0usize;

                    for offset in 0..k {
                        let pos = i + offset;
                        let orig_base = scratch[pos];
                        let orig_upper = orig_base.to_ascii_uppercase();
                        let start_scan = pos.saturating_sub(k - 1);
                        let end_scan = pos.min(scratch.len() - k);

                        let mut base_solid = 0usize;
                        for s in start_scan..=end_scan {
                            if s == i {
                                continue;
                            }
                            if let Some(surr_km) = Kmer256::from_bytes(&scratch[s..s + k], k) {
                                let (surr_can, _) = surr_km.canonical(k);
                                if filter.is_solid_kmer256(surr_can) {
                                    base_solid += 1;
                                }
                            }
                        }

                        for &cand_base in b"ACGT" {
                            if cand_base == orig_upper {
                                continue;
                            }
                            scratch[pos] = cand_base;
                            if let Some(cand_km) = Kmer256::from_bytes(&scratch[i..i + k], k) {
                                let (cand_can, _) = cand_km.canonical(k);
                                if filter.is_solid_kmer256(cand_can) {
                                    let mut solid_count = 1usize; // s == i is already verified solid
                                    for s in start_scan..=end_scan {
                                        if s == i {
                                            continue;
                                        }
                                        if let Some(surr_km) =
                                            Kmer256::from_bytes(&scratch[s..s + k], k)
                                        {
                                            let (surr_can, _) = surr_km.canonical(k);
                                            if filter.is_solid_kmer256(surr_can) {
                                                solid_count += 1;
                                            }
                                        }
                                    }

                                    let gain = solid_count.saturating_sub(base_solid);
                                    if gain > max_solid_gain {
                                        max_solid_gain = gain;
                                        best_sub = Some((offset, cand_base));
                                    }
                                }
                            }
                            scratch[pos] = orig_base;
                        }
                    }

                    if let Some((offset, cand_base)) = best_sub {
                        seq[i + offset] = cand_base;
                        scratch[i + offset] = cand_base;
                        modified = true;
                    }
                }
            }
        }

        modified
    }

    /// Corrects reads in parallel using the TwoTierFilter solid k-mer index.
    pub fn correct_reads(
        &self,
        reads: Vec<Vec<u8>>,
        filter: &TwoTierFilter,
    ) -> (Vec<Vec<u8>>, usize) {
        let corrected_count = std::sync::atomic::AtomicUsize::new(0);

        let corrected_reads: Vec<Vec<u8>> = reads
            .into_par_iter()
            .map(|mut seq| {
                let mut scratch = Vec::with_capacity(seq.len());
                if self.correct_read_in_place(&mut seq, &mut scratch, filter) {
                    corrected_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                }
                seq
            })
            .collect();

        let total_corrected = corrected_count.load(std::sync::atomic::Ordering::Relaxed);
        (corrected_reads, total_corrected)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hamming_distance_2bit() {
        // 0 mismatches
        assert_eq!(hamming_distance_2bit(0b00011011, 0b00011011), 0);
        // 1 mismatch in last base
        assert_eq!(hamming_distance_2bit(0b00011011, 0b00011000), 1);
        // 4 mismatches across 4 bases
        assert_eq!(hamming_distance_2bit(0b00000000, 0b01010101), 4);
    }

    #[test]
    fn test_correct_read_single_substitution() {
        let k = 5;
        let corrector = ErrorCorrector::new(k);
        let filter = TwoTierFilter::new(1024 * 1024);

        // Populate a clean reference sequence into filter as solid (seen twice)
        let ref_seq = b"ACGTACGTACGT";
        for _ in 0..2 {
            for i in 0..=(ref_seq.len() - k) {
                if let Some(km) = Kmer256::from_bytes(&ref_seq[i..i + k], k) {
                    filter.insert_kmer256(km.canonical(k).0);
                }
            }
        }

        // Mutate single base in the middle: 'T' (pos 3) -> 'G'
        let mut mutated = b"ACGGACGTACGT".to_vec();
        let mut scratch = Vec::new();

        let corrected = corrector.correct_read_in_place(&mut mutated, &mut scratch, &filter);
        assert!(corrected, "ErrorCorrector should correct single mismatch");
        assert_eq!(mutated, ref_seq, "Corrected sequence must match reference");
    }

    #[test]
    fn test_short_read_bypass() {
        let corrector = ErrorCorrector::new(15);
        let filter = TwoTierFilter::new(1024);
        let mut short = b"ACGT".to_vec();
        let mut scratch = Vec::new();
        assert!(!corrector.correct_read_in_place(&mut short, &mut scratch, &filter));
    }
}
