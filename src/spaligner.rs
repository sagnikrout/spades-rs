//! Hybrid Long-Read Resolver (Module 3.3: Spaligner).
//!
//! Maps PacBio HiFi and Oxford Nanopore (ONT) long reads across assembly graph unitigs
//! using seed chaining and copy-number repeat unrolling to resolve multi-kilobase
//! repeats (such as bacterial ribosomal RNA operons).

use crate::dna::{canonical_kmer_u64, revcomp_bytes, string_to_kmer};
use crate::graph::Unitig;
use hashbrown::{HashMap, HashSet};
use rayon::prelude::*;

pub struct LongReadResolver {
    pub k: usize,
    pub min_seed_matches: usize,
}

impl Default for LongReadResolver {
    fn default() -> Self {
        Self {
            k: 31,
            min_seed_matches: 3,
        }
    }
}

type BridgePair = (usize, usize);
type BridgeConfig = (bool, Vec<(usize, bool)>, bool);

#[derive(Clone, Debug)]
struct ReadAnchor {
    u_idx: usize,
    is_rev: bool,
    start: usize,
    end: usize,
    hits: usize,
}

impl LongReadResolver {
    /// Bridges short-read unitigs using long reads by unrolling copy-number repeats
    /// and resolving multi-copy unitig chains.
    pub fn bridge_with_long_reads(
        &self,
        unitigs: Vec<Unitig>,
        long_reads: &[Vec<u8>],
    ) -> Vec<Unitig> {
        if long_reads.is_empty() || unitigs.len() <= 1 {
            return unitigs;
        }

        let k = self.k;
        let seed_k = 15.min(k).max(4);
        let step = 5.min((seed_k / 3).max(1));
        let k1 = k.saturating_sub(1);
        let min_ov = if unitigs.len() <= 3 { (k1 / 2).max(1) } else { 15 };
        let min_read_len = if unitigs.len() <= 3 { seed_k * 2 } else { 1500 };

        // 1. Calculate expected coverage depth to distinguish repeats from unique contigs
        let mut covs: Vec<f64> = unitigs
            .iter()
            .filter(|u| u.sequence.len() >= 500)
            .map(|u| u.mean_coverage)
            .collect();
        if covs.is_empty() {
            covs = unitigs.iter().map(|u| u.mean_coverage).collect();
        }
        covs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median_cov = if !covs.is_empty() {
            covs[covs.len() / 2]
        } else {
            10.0
        };

        let repeat_threshold = median_cov * 1.8;
        let is_unique_flank: Vec<bool> = if unitigs.len() <= 3 {
            vec![true; unitigs.len()]
        } else {
            unitigs
                .iter()
                .map(|u| u.sequence.len() >= 1500 && u.mean_coverage < median_cov * 1.35)
                .collect()
        };
        let is_repeat: Vec<bool> = unitigs
            .iter()
            .enumerate()
            .map(|(i, _)| !is_unique_flank[i])
            .collect();
        let is_major_repeat: Vec<bool> = unitigs
            .iter()
            .map(|u| u.mean_coverage >= median_cov * 2.0 && u.sequence.len() >= 1000)
            .collect();

        println!(
            "  [Spaligner] Median coverage: {:.1}x, Repeat cutoff (>=1.8x): {:.1}x (Detected {} repeat unitigs, {} unique flanks)",
            median_cov,
            repeat_threshold,
            is_repeat.iter().filter(|&&r| r).count(),
            is_unique_flank.iter().filter(|&&u| u).count()
        );

        // 2. Build seed k-mer index mapping canonical seed -> Vec<(unitig_idx, pos, is_rc)>
        let mut seed_to_unitig: HashMap<u64, Vec<(u32, u32, bool)>> =
            HashMap::with_capacity(unitigs.len() * 50);

        for (idx, u) in unitigs.iter().enumerate() {
            if u.sequence.len() < seed_k {
                continue;
            }
            for i in (0..=(u.sequence.len() - seed_k)).step_by(step) {
                if let Some(km) = string_to_kmer(&u.sequence[i..i + seed_k], seed_k) {
                    let (can, is_rc) = canonical_kmer_u64(km, seed_k);
                    seed_to_unitig.entry(can).or_default().push((idx as u32, i as u32, is_rc));
                }
            }
        }

        // 3. Trace ordered unitig anchors along each long read in parallel
        let min_support = self.min_seed_matches;
        let bridge_lists: Vec<HashMap<BridgePair, Vec<BridgeConfig>>> = long_reads
            .par_chunks(256)
            .map(|chunk| {
                let mut local_bridges: HashMap<BridgePair, Vec<BridgeConfig>> = HashMap::new();

                for read in chunk {
                    if read.len() < min_read_len {
                        continue;
                    }

                    // Collect matching seeds on this read
                    let mut hits: Vec<(usize, bool, usize)> = Vec::new();
                    for i in (0..=(read.len() - seed_k)).step_by(step) {
                        if let Some(km) = string_to_kmer(&read[i..i + seed_k], seed_k) {
                            let (can, r_is_rc) = canonical_kmer_u64(km, seed_k);
                            if let Some(matches) = seed_to_unitig.get(&can) {
                                for &(u_idx, _, u_is_rc) in matches {
                                    hits.push((u_idx as usize, r_is_rc != u_is_rc, i));
                                }
                            }
                        }
                    }

                    if hits.is_empty() {
                        continue;
                    }

                    // Cluster hits into continuous anchors
                    hits.sort_by_key(|h| (h.0, h.1, h.2));
                    let mut raw_anchors: Vec<ReadAnchor> = Vec::new();

                    let mut cur_u = hits[0].0;
                    let mut cur_rev = hits[0].1;
                    let mut min_pos = hits[0].2;
                    let mut max_pos = hits[0].2;
                    let mut h_count = 1;

                    for &(u, rev, pos) in &hits[1..] {
                        if u == cur_u && rev == cur_rev && pos.saturating_sub(max_pos) <= 1500 {
                            max_pos = pos;
                            h_count += 1;
                        } else {
                            let span = max_pos - min_pos + seed_k;
                            let min_span = if unitigs.len() <= 3 { seed_k } else if is_unique_flank[cur_u] { 500 } else { 200 };
                            if h_count >= min_support && span >= min_span {
                                raw_anchors.push(ReadAnchor {
                                    u_idx: cur_u,
                                    is_rev: cur_rev,
                                    start: min_pos,
                                    end: max_pos,
                                    hits: h_count,
                                });
                            }

                            cur_u = u;
                            cur_rev = rev;
                            min_pos = pos;
                            max_pos = pos;
                            h_count = 1;
                        }
                    }

                    let span = max_pos - min_pos + seed_k;
                    let min_span = if unitigs.len() <= 3 { seed_k } else if is_unique_flank[cur_u] { 500 } else { 200 };
                    if h_count >= min_support && span >= min_span {
                        raw_anchors.push(ReadAnchor {
                            u_idx: cur_u,
                            is_rev: cur_rev,
                            start: min_pos,
                            end: max_pos,
                            hits: h_count,
                        });
                    }

                    if raw_anchors.len() < 2 {
                        continue;
                    }

                    // Sort anchors along the read
                    raw_anchors.sort_by_key(|a| a.start);

                    // Collapse consecutive anchors of same unitig
                    let mut ordered_anchors: Vec<ReadAnchor> = Vec::with_capacity(raw_anchors.len());
                    for a in raw_anchors {
                        if let Some(last) = ordered_anchors.last_mut() {
                            if last.u_idx == a.u_idx && a.start.saturating_sub(last.end) <= 1500 {
                                if a.hits > last.hits {
                                    last.is_rev = a.is_rev;
                                }
                                last.end = a.end.max(last.end);
                                last.hits += a.hits;
                                continue;
                            }
                        }
                        ordered_anchors.push(a);
                    }

                    // Extract bridges from unique flank -> repeats -> unique flank
                    let n_anchors = ordered_anchors.len();
                    for i in 0..n_anchors {
                        let mut left_anchor = &ordered_anchors[i];
                        if !is_unique_flank[left_anchor.u_idx] {
                            continue;
                        }

                        let mut reps_list = Vec::new();
                        for j in (i + 1)..n_anchors {
                            let cand = &ordered_anchors[j];
                            if cand.u_idx == left_anchor.u_idx {
                                left_anchor = cand;
                                reps_list.clear();
                                continue;
                            }
                            if is_repeat[cand.u_idx] {
                                if reps_list.last().map(|&(u, _)| u) != Some(cand.u_idx) {
                                    reps_list.push((cand.u_idx, cand.is_rev));
                                }
                            } else {
                                // Found right unique flank
                                let u_l = left_anchor.u_idx;
                                let u_r = cand.u_idx;
                                let rc_l = left_anchor.is_rev;
                                let rc_r = cand.is_rev;

                                if u_l <= u_r {
                                    let pair_key = (u_l, u_r);
                                    let config_val = (rc_l, reps_list, rc_r);
                                    local_bridges.entry(pair_key).or_default().push(config_val);
                                } else {
                                    let pair_key = (u_r, u_l);
                                    let reps_rc: Vec<(usize, bool)> = reps_list
                                        .into_iter()
                                        .rev()
                                        .map(|(u, rev)| (u, !rev))
                                        .collect();
                                    let config_val = (!rc_r, reps_rc, !rc_l);
                                    local_bridges.entry(pair_key).or_default().push(config_val);
                                }
                                break;
                            }
                        }
                    }
                }

                local_bridges
            })
            .collect();

        let mut all_bridges: HashMap<BridgePair, Vec<BridgeConfig>> = HashMap::new();
        for l_bridges in bridge_lists {
            for (pair, configs) in l_bridges {
                all_bridges.entry(pair).or_default().extend(configs);
            }
        }

        println!(
            "  [Spaligner] Evaluated {} long reads, identified {} unique flank bridge pairs",
            long_reads.len(),
            all_bridges.len()
        );

        // 4. Resolve bridges by unrolling repeats with rRNA prioritization
        let mut sorted_bridges: Vec<((usize, usize), Vec<BridgeConfig>)> = all_bridges
            .into_iter()
            .filter(|(_, obs)| obs.len() >= min_support)
            .collect();

        sorted_bridges.sort_by(|a, b| {
            let a_spans_major = a.1.iter().any(|c| c.1.iter().any(|&(u, _)| is_major_repeat[u]));
            let b_spans_major = b.1.iter().any(|c| c.1.iter().any(|&(u, _)| is_major_repeat[u]));
            b_spans_major
                .cmp(&a_spans_major)
                .then_with(|| b.1.len().cmp(&a.1.len()))
        });

        let mut consumed_unique: HashSet<usize> = HashSet::new();
        let mut unrolled_repeats: HashSet<usize> = HashSet::new();
        let mut stitched_contigs: Vec<Unitig> = Vec::new();

        for ((u_l, u_r), obs) in sorted_bridges {
            if consumed_unique.contains(&u_l) || consumed_unique.contains(&u_r) {
                continue;
            }

            // Select best configuration (dominant orientation and most complete repeat chain)
            let mut counts: HashMap<BridgeConfig, usize> = HashMap::new();
            for cfg in &obs {
                *counts.entry(cfg.clone()).or_insert(0) += 1;
            }
            let (best_config, _) = counts
                .into_iter()
                .max_by(|(c1, sup1), (c2, sup2)| {
                    c1.1.len().cmp(&c2.1.len()).then_with(|| sup1.cmp(sup2))
                })
                .unwrap();

            let (rc_l, reps, rc_r) = best_config;

            let left_seq = if rc_l {
                revcomp_bytes(&unitigs[u_l].sequence)
            } else {
                unitigs[u_l].sequence.clone()
            };

            let mut assembled = left_seq;
            let mut total_kmers = unitigs[u_l].kmers_count;
            let mut cov_sum = unitigs[u_l].mean_coverage * unitigs[u_l].sequence.len() as f64;
            let mut total_bp = unitigs[u_l].sequence.len();

            // Intermediate repeats
            for &(rep_idx, rep_rev) in &reps {
                let rep_seq = if rep_rev {
                    revcomp_bytes(&unitigs[rep_idx].sequence)
                } else {
                    unitigs[rep_idx].sequence.clone()
                };

                let mut best_overlap = 0;
                let max_ov = k1.min(assembled.len()).min(rep_seq.len());
                if max_ov >= min_ov {
                    for ov in (min_ov..=max_ov).rev() {
                        if assembled[assembled.len() - ov..] == rep_seq[..ov] {
                            best_overlap = ov;
                            break;
                        }
                    }
                }

                if best_overlap >= min_ov {
                    assembled.extend_from_slice(&rep_seq[best_overlap..]);
                    total_bp += rep_seq.len() - best_overlap;
                    cov_sum += unitigs[rep_idx].mean_coverage * (rep_seq.len() - best_overlap) as f64;
                } else {
                    assembled.extend_from_slice(b"NNNNNNNNNNNNNNNNNNNN");
                    assembled.extend_from_slice(&rep_seq);
                    total_bp += rep_seq.len() + 20;
                    cov_sum += unitigs[rep_idx].mean_coverage * rep_seq.len() as f64;
                }
                total_kmers += unitigs[rep_idx].kmers_count;
                unrolled_repeats.insert(rep_idx);
            }

            // Right flank
            let right_seq = if rc_r {
                revcomp_bytes(&unitigs[u_r].sequence)
            } else {
                unitigs[u_r].sequence.clone()
            };

            let mut best_overlap = 0;
            let max_ov = k1.min(assembled.len()).min(right_seq.len());
            if max_ov >= min_ov {
                for ov in (min_ov..=max_ov).rev() {
                    if assembled[assembled.len() - ov..] == right_seq[..ov] {
                        best_overlap = ov;
                        break;
                    }
                }
            }

            if best_overlap >= min_ov {
                assembled.extend_from_slice(&right_seq[best_overlap..]);
                total_bp += right_seq.len() - best_overlap;
                cov_sum += unitigs[u_r].mean_coverage * (right_seq.len() - best_overlap) as f64;
            } else {
                assembled.extend_from_slice(b"NNNNNNNNNNNNNNNNNNNN");
                assembled.extend_from_slice(&right_seq);
                total_bp += right_seq.len() + 20;
                cov_sum += unitigs[u_r].mean_coverage * right_seq.len() as f64;
            }
            total_kmers += unitigs[u_r].kmers_count;

            consumed_unique.insert(u_l);
            consumed_unique.insert(u_r);

            stitched_contigs.push(Unitig {
                id: stitched_contigs.len(),
                sequence: assembled,
                mean_coverage: if total_bp > 0 {
                    cov_sum / total_bp as f64
                } else {
                    median_cov
                },
                kmers_count: total_kmers,
            });
        }

        println!(
            "  [Spaligner] Successfully unrolled {} repeat bridges ({} unique flanks stitched, {} repeat nodes unrolled)",
            stitched_contigs.len(),
            consumed_unique.len(),
            unrolled_repeats.len()
        );

        // Collect final unitigs
        let mut final_unitigs = stitched_contigs;
        for (idx, u) in unitigs.into_iter().enumerate() {
            if !consumed_unique.contains(&idx) {
                final_unitigs.push(u);
            }
        }

        final_unitigs.sort_by_key(|u| std::cmp::Reverse(u.sequence.len()));
        final_unitigs
    }
}

