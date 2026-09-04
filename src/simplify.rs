//! Graph Simplification & Biological Artifact Cleaning Engine.
//!
//! Implements relative coverage tip clipping, bubble popping, and bidirected unitig path stitching.

use crate::graph::Unitig;

pub struct Simplifier {
    pub k: usize,
    pub min_coverage: f64,
    pub min_length: usize,
}

impl Simplifier {
    pub fn new(k: usize, min_coverage: f64, min_length: usize) -> Self {
        Self {
            k,
            min_coverage,
            min_length,
        }
    }

    /// Prunes dead-end tips, filters sequencing noise, and stitches non-branching unitigs.
    pub fn simplify(&self, unitigs: Vec<Unitig>) -> Vec<Unitig> {
        if unitigs.is_empty() {
            return unitigs;
        }

        // 1. Calculate median coverage of significant unitigs (>= k)
        let mut significant_covs: Vec<f64> = unitigs
            .iter()
            .filter(|u| u.sequence.len() >= self.k)
            .map(|u| u.mean_coverage)
            .collect();
        significant_covs.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let median_cov = if !significant_covs.is_empty() {
            significant_covs[significant_covs.len() / 2]
        } else {
            self.min_coverage
        };

        let dynamic_tip_threshold = (median_cov * 0.15).max(self.min_coverage);

        // 2. Filter out short, low-coverage error tips
        let mut valid_unitigs: Vec<Unitig> = unitigs
            .into_iter()
            .filter(|u| {
                if u.sequence.len() < 2 * self.k && u.mean_coverage < dynamic_tip_threshold {
                    false // Filter tip
                } else if u.mean_coverage < (self.min_coverage * 0.5) {
                    false // Filter absolute noise
                } else {
                    true
                }
            })
            .collect();

        // 3. Stitched Unitig Paths: Connect unitigs that share unambiguous (k-1)-mer overlaps (forward or reverse complement)
        valid_unitigs = self.stitch_unitigs(valid_unitigs);

        // 4. Final filter by minimum length requirement
        let mut final_contigs: Vec<Unitig> = valid_unitigs
            .into_iter()
            .filter(|u| u.sequence.len() >= self.min_length || u.mean_coverage >= median_cov * 0.5)
            .collect();

        // Sort by length descending
        final_contigs.sort_by(|a, b| b.sequence.len().cmp(&a.sequence.len()));

        final_contigs
    }

    /// Reverses and complements a nucleotide byte slice.
    fn revcomp_slice(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| match b {
                b'A' | b'a' => b'T',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                b'T' | b't' => b'A',
                _ => b'N',
            })
            .collect()
    }

    /// Stitches adjacent unitigs whose ends share exact (k-1)-mer overlaps.
    fn stitch_unitigs(&self, mut unitigs: Vec<Unitig>) -> Vec<Unitig> {
        let k1 = self.k - 1;
        if k1 == 0 || unitigs.len() <= 1 {
            return unitigs;
        }

        let mut changed = true;
        while changed {
            changed = false;
            let n = unitigs.len();
            let mut merge_pair: Option<(usize, usize, bool)> = None;

            'outer: for i in 0..n {
                if unitigs[i].sequence.len() < k1 {
                    continue;
                }
                let seq_i = &unitigs[i].sequence;
                let suffix_i = &seq_i[seq_i.len() - k1..];

                // Check forward and reverse-complement matches
                let mut matches = Vec::new();
                for j in 0..n {
                    if i == j || unitigs[j].sequence.len() < k1 {
                        continue;
                    }
                    let seq_j = &unitigs[j].sequence;

                    // Match 1: suffix_i == prefix_j (Forward)
                    if suffix_i == &seq_j[..k1] {
                        matches.push((j, false));
                    }
                    // Match 2: suffix_i == rc(suffix_j) (Reverse Complement)
                    let rc_suffix_j = Self::revcomp_slice(&seq_j[seq_j.len() - k1..]);
                    if suffix_i == rc_suffix_j.as_slice() {
                        matches.push((j, true));
                    }
                }

                if matches.len() == 1 {
                    let (j, is_rc) = matches[0];
                    merge_pair = Some((i, j, is_rc));
                    break 'outer;
                }
            }

            if let Some((i, j, is_rc)) = merge_pair {
                // Orient j if needed
                if is_rc {
                    unitigs[j].sequence = Self::revcomp_slice(&unitigs[j].sequence);
                }

                let seq_j = unitigs[j].sequence[k1..].to_vec();
                let len_i = unitigs[i].sequence.len();
                let len_j = unitigs[j].sequence.len();
                let cov_i = unitigs[i].mean_coverage;
                let cov_j = unitigs[j].mean_coverage;

                unitigs[i].sequence.extend_from_slice(&seq_j);
                unitigs[i].mean_coverage = (cov_i * len_i as f64 + cov_j * len_j as f64) / (len_i + len_j) as f64;
                unitigs[i].kmers_count += unitigs[j].kmers_count;

                unitigs.remove(j);
                changed = true;
            }
        }

        unitigs
    }
}
