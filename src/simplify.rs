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

    /// Prunes dead-end tips, pops error bubbles, and stitches non-branching unitigs iteratively.
    pub fn simplify(&self, mut unitigs: Vec<Unitig>) -> Vec<Unitig> {
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

        // Iterative simplification loop: alternate tip clipping, bubble popping, and path stitching
        let mut iteration = 0;
        let max_iterations = 8;

        while iteration < max_iterations {
            iteration += 1;
            let count_before = unitigs.len();

            // Step A: Clip short low-coverage tips and absolute noise
            unitigs = self.clip_tips(unitigs, dynamic_tip_threshold);

            // Step B: Pop parallel sequencing error bubbles (Bulge Remover)
            unitigs = self.pop_bubbles(unitigs);

            // Step C: Stitch unambiguous unitig paths (in-degree == 1 and out-degree == 1)
            unitigs = self.stitch_unitigs(unitigs);

            if unitigs.len() == count_before {
                break; // Converged to fixed point
            }
        }

        // Final filter by minimum length requirement
        let mut final_contigs: Vec<Unitig> = unitigs
            .into_iter()
            .filter(|u| u.sequence.len() >= self.min_length)
            .collect();

        // Sort by length descending
        final_contigs.sort_by_key(|a| std::cmp::Reverse(a.sequence.len()));

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

    /// Filters dead-end tips and sequencing artifacts with length-awareness.
    pub fn clip_tips(&self, unitigs: Vec<Unitig>, dynamic_tip_threshold: f64) -> Vec<Unitig> {
        let tip_max_len = 2 * self.k;
        unitigs
            .into_iter()
            .filter(|u| {
                // If a dead-end branch is longer than or equal to 2*k base pairs,
                // do not prune it based on low coverage alone.
                // True sequencing errors virtually never generate an unbroken, unbranched error chain of 2*k bases.
                // Keep the low-coverage and noise cutoff strictly for short tips (< 2*k).
                if u.sequence.len() >= tip_max_len {
                    true
                } else {
                    let is_tip = u.mean_coverage < dynamic_tip_threshold;
                    let is_noise = u.mean_coverage < (self.min_coverage * 0.5);
                    !is_tip && !is_noise
                }
            })
            .collect()
    }

    /// Identifies and removes parallel bubble paths (bulges) between common junctions.
    pub fn pop_bubbles(&self, unitigs: Vec<Unitig>) -> Vec<Unitig> {
        let k1 = self.k - 1;
        if k1 == 0 || unitigs.len() <= 1 {
            return unitigs;
        }

        // Map canonical endpoint pair (start_(k-1), end_(k-1)) -> list of unitig indices
        let mut bubble_groups: hashbrown::HashMap<(Vec<u8>, Vec<u8>), Vec<usize>> =
            hashbrown::HashMap::new();

        for (idx, u) in unitigs.iter().enumerate() {
            if u.sequence.len() < k1 {
                continue;
            }
            let p = u.sequence[..k1].to_vec();
            let q = u.sequence[u.sequence.len() - k1..].to_vec();
            let rc_p = Self::revcomp_slice(&p);
            let rc_q = Self::revcomp_slice(&q);

            let key_fwd = (p, q);
            let key_rev = (rc_q, rc_p);

            let canonical_key = if key_fwd <= key_rev { key_fwd } else { key_rev };

            bubble_groups.entry(canonical_key).or_default().push(idx);
        }

        let mut to_remove = hashbrown::HashSet::new();

        for (_endpoints, members) in bubble_groups {
            if members.len() <= 1 {
                continue;
            }

            // Find member with highest coverage (tie-break with length)
            let mut best_idx = members[0];
            let mut best_cov = unitigs[best_idx].mean_coverage;
            let mut best_len = unitigs[best_idx].sequence.len();

            for &idx in &members[1..] {
                let cov = unitigs[idx].mean_coverage;
                let len = unitigs[idx].sequence.len();
                if cov > best_cov || (cov == best_cov && len > best_len) {
                    best_idx = idx;
                    best_cov = cov;
                    best_len = len;
                }
            }

            // Pop inferior members if they meet bubble criteria
            for &idx in &members {
                if idx == best_idx {
                    continue;
                }
                let len = unitigs[idx].sequence.len();
                let len_diff = (len as isize - best_len as isize).unsigned_abs();

                // Bubble criteria: length difference within 3k or coverage ratio >= 2.0
                if len_diff <= 3 * self.k || (best_cov / unitigs[idx].mean_coverage.max(0.1)) >= 2.0
                {
                    to_remove.insert(idx);
                }
            }
        }

        if to_remove.is_empty() {
            return unitigs;
        }

        unitigs
            .into_iter()
            .enumerate()
            .filter(|(idx, _)| !to_remove.contains(idx))
            .map(|(_, u)| u)
            .collect()
    }

    /// Stitches adjacent unitigs whose ends share exact (k-1)-mer overlaps with in-degree == 1 and out-degree == 1.
    pub fn stitch_unitigs(&self, mut unitigs: Vec<Unitig>) -> Vec<Unitig> {
        let k1 = self.k - 1;
        if k1 == 0 || unitigs.len() <= 1 {
            return unitigs;
        }

        let mut changed = true;
        while changed {
            changed = false;
            let n = unitigs.len();

            // Build index of prefixes and suffixes for degree checking
            let mut prefix_map: hashbrown::HashMap<Vec<u8>, Vec<(usize, bool)>> =
                hashbrown::HashMap::new();
            let mut suffix_map: hashbrown::HashMap<Vec<u8>, Vec<(usize, bool)>> =
                hashbrown::HashMap::new();

            for (i, u_i) in unitigs.iter().enumerate().take(n) {
                if u_i.sequence.len() < k1 {
                    continue;
                }
                let seq = &u_i.sequence;
                let prefix = seq[..k1].to_vec();
                let suffix = seq[seq.len() - k1..].to_vec();
                let rc_prefix = Self::revcomp_slice(&prefix);
                let rc_suffix = Self::revcomp_slice(&suffix);

                // Forward orientations
                prefix_map.entry(prefix).or_default().push((i, false));
                suffix_map.entry(suffix).or_default().push((i, false));

                // RC orientations
                prefix_map.entry(rc_suffix).or_default().push((i, true));
                suffix_map.entry(rc_prefix).or_default().push((i, true));
            }

            let mut merge_pair: Option<(usize, usize, bool)> = None;

            for (i, u_i) in unitigs.iter().enumerate().take(n) {
                if u_i.sequence.len() < k1 {
                    continue;
                }
                let suffix_i = u_i.sequence[u_i.sequence.len() - k1..].to_vec();

                if let Some(matches) = prefix_map.get(&suffix_i) {
                    let valid_matches: Vec<_> = matches.iter().filter(|&&(j, _)| i != j).collect();
                    let mut unique_targets: Vec<(usize, bool)> = Vec::new();
                    for &&(j, is_rc) in &valid_matches {
                        if !unique_targets.iter().any(|&(uj, _)| uj == j) {
                            unique_targets.push((j, is_rc));
                        }
                    }

                    if unique_targets.len() == 1 {
                        let (j, is_rc) = unique_targets[0];

                        // Degree check: ensure j's incoming endpoint also has exactly 1 predecessor (i)
                        let j_prefix = if !is_rc {
                            unitigs[j].sequence[..k1].to_vec()
                        } else {
                            Self::revcomp_slice(
                                &unitigs[j].sequence[unitigs[j].sequence.len() - k1..],
                            )
                        };

                        if let Some(in_matches) = suffix_map.get(&j_prefix) {
                            let valid_in: Vec<_> =
                                in_matches.iter().filter(|&&(src, _)| src != j).collect();
                            let mut unique_in: Vec<usize> = Vec::new();
                            for &&(src, _) in &valid_in {
                                if !unique_in.contains(&src) {
                                    unique_in.push(src);
                                }
                            }

                            if unique_in.len() == 1 && unique_in[0] == i {
                                merge_pair = Some((i, j, is_rc));
                                break;
                            }
                        }
                    }
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
                unitigs[i].mean_coverage =
                    (cov_i * len_i as f64 + cov_j * len_j as f64) / (len_i + len_j) as f64;
                unitigs[i].kmers_count += unitigs[j].kmers_count;

                unitigs.remove(j);
                changed = true;
            }
        }

        unitigs
    }
}
