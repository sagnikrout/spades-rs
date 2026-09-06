//! Hybrid Long-Read Resolver (Module 3.3: Spaligner).
//!
//! Maps PacBio HiFi and Oxford Nanopore (ONT) long reads across assembly graph unitigs
//! using seed chaining, repeat seed filtering, and bidirected port involution chaining
//! to resolve multi-kilobase repeats and bridge unitigs into chromosome-scale scaffolds.

use crate::dna::{canonical_kmer_u64, revcomp_bytes, string_to_kmer};
use crate::graph::Unitig;
use hashbrown::HashMap;
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

struct DisjointSet {
    parent: Vec<usize>,
}

impl DisjointSet {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
        }
    }

    fn find(&mut self, i: usize) -> usize {
        let mut root = i;
        while root != self.parent[root] {
            root = self.parent[root];
        }
        let mut curr = i;
        while curr != root {
            let nxt = self.parent[curr];
            self.parent[curr] = root;
            curr = nxt;
        }
        root
    }

    fn union(&mut self, i: usize, j: usize) -> bool {
        let root_i = self.find(i);
        let root_j = self.find(j);
        if root_i == root_j {
            false
        } else {
            self.parent[root_i] = root_j;
            true
        }
    }
}

impl LongReadResolver {
    /// Bridges short-read unitigs using long reads by unrolling copy-number repeats
    /// and resolving multi-copy unitig chains into chromosome scaffolds.
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
        let min_ov = if unitigs.len() <= 3 {
            (k1 / 2).max(1)
        } else {
            10
        };
        let min_read_len = if unitigs.len() <= 3 {
            seed_k * 2
        } else {
            500.max(seed_k * 2)
        };
        let min_span = if unitigs.len() <= 3 {
            seed_k
        } else {
            150.max(seed_k * 2)
        };
        let n_unitigs = unitigs.len();

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
                .map(|u| u.sequence.len() >= 500 && u.mean_coverage < median_cov * 1.5)
                .collect()
        };
        let is_repeat: Vec<bool> = is_unique_flank.iter().map(|&u| !u).collect();

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
                    seed_to_unitig
                        .entry(can)
                        .or_default()
                        .push((idx as u32, i as u32, is_rc));
                }
            }
        }

        // Filter high-frequency repeat seeds for large eukaryotic graphs
        if unitigs.len() > 10 {
            let initial_seeds = seed_to_unitig.len();
            seed_to_unitig.retain(|_, matches| matches.len() <= 5);
            let kept = seed_to_unitig.len();
            if initial_seeds > 0 {
                println!(
                    "  [Spaligner] Filtered repeat seeds: retained {} / {} ({:.1}%)",
                    kept,
                    initial_seeds,
                    (kept as f64 / initial_seeds as f64) * 100.0
                );
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
                            let required_span = if is_unique_flank[cur_u] {
                                min_span
                            } else {
                                (min_span / 2).max(seed_k)
                            };
                            if h_count >= min_support.min(3) && span >= required_span {
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
                    let required_span = if is_unique_flank[cur_u] {
                        min_span
                    } else {
                        (min_span / 2).max(seed_k)
                    };
                    if h_count >= min_support.min(3) && span >= required_span {
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

                    raw_anchors.sort_by_key(|a| a.start);

                    let mut ordered_anchors: Vec<ReadAnchor> =
                        Vec::with_capacity(raw_anchors.len());
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

                    let n_anchors = ordered_anchors.len();
                    for i in 0..n_anchors {
                        let mut left_anchor = &ordered_anchors[i];
                        if !is_unique_flank[left_anchor.u_idx] {
                            continue;
                        }

                        let mut reps_list = Vec::new();
                        #[allow(clippy::needless_range_loop)]
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

        // 4. Resolve bridges using Bidirected Port Involution Chaining
        let mut sorted_bridges: Vec<((usize, usize), Vec<BridgeConfig>)> = all_bridges
            .into_iter()
            .filter(|(_, obs)| obs.len() >= min_support)
            .collect();

        sorted_bridges.sort_by_key(|a| std::cmp::Reverse(a.1.len()));

        let total_ports = 2 * n_unitigs;
        #[allow(clippy::type_complexity)]
        let mut partner: Vec<Option<(usize, Vec<(usize, bool)>, usize)>> = vec![None; total_ports];
        let mut dsu = DisjointSet::new(n_unitigs);
        let mut chained_bridges = 0;

        for ((u_l, u_r), obs) in &sorted_bridges {
            let mut counts: HashMap<&BridgeConfig, usize> = HashMap::new();
            for cfg in obs {
                *counts.entry(cfg).or_insert(0) += 1;
            }
            let (best_cfg, sup) = counts.into_iter().max_by_key(|&(_, c)| c).unwrap();

            // For larger graphs, enforce dominant orientation consistency >= 70%
            if unitigs.len() > 10 && (sup as f64) / (obs.len() as f64) < 0.70 {
                continue;
            }

            let &(rc_l, ref reps, rc_r) = best_cfg;

            let p_l = if rc_l { 2 * u_l } else { 2 * u_l + 1 };
            let p_r = if rc_r { 2 * u_r + 1 } else { 2 * u_r };

            if partner[p_l].is_some() || partner[p_r].is_some() {
                continue;
            }

            if !dsu.union(*u_l, *u_r) {
                continue;
            }

            partner[p_l] = Some((p_r, reps.clone(), sup));
            let reps_rc: Vec<(usize, bool)> =
                reps.iter().rev().map(|&(u, rev)| (u, !rev)).collect();
            partner[p_r] = Some((p_l, reps_rc, sup));
            chained_bridges += 1;
        }

        println!(
            "  [Spaligner] Successfully connected {} bridges into scaffolds",
            chained_bridges
        );

        // 5. Reconstruct linear scaffold paths
        let mut visited_unitigs = vec![false; n_unitigs];
        let mut stitched_contigs: Vec<Unitig> = Vec::new();

        for port in 0..total_ports {
            let start_u = port / 2;
            if visited_unitigs[start_u] {
                continue;
            }

            if partner[port].is_some() {
                continue;
            }

            let opposite_port = port ^ 1;
            if partner[opposite_port].is_none() {
                continue;
            }

            let mut cur_u = start_u;
            let mut cur_is_rev = port % 2 == 1;

            let mut path_unitigs: Vec<(usize, bool)> = Vec::new();
            let mut path_reps: Vec<Vec<(usize, bool)>> = Vec::new();

            loop {
                visited_unitigs[cur_u] = true;
                path_unitigs.push((cur_u, cur_is_rev));

                let out_port = if cur_is_rev { 2 * cur_u } else { 2 * cur_u + 1 };

                if let Some((next_port, ref reps, _)) = partner[out_port] {
                    let next_u = next_port / 2;
                    if visited_unitigs[next_u] {
                        break;
                    }
                    path_reps.push(reps.clone());
                    let next_is_rev = next_port % 2 == 1;
                    cur_u = next_u;
                    cur_is_rev = next_is_rev;
                } else {
                    break;
                }
            }

            if path_unitigs.len() >= 2 {
                let first_u = path_unitigs[0].0;
                let first_rev = path_unitigs[0].1;
                let mut seq = if first_rev {
                    revcomp_bytes(&unitigs[first_u].sequence)
                } else {
                    unitigs[first_u].sequence.clone()
                };

                let mut cov_sum =
                    unitigs[first_u].mean_coverage * unitigs[first_u].sequence.len() as f64;
                let mut total_bp = unitigs[first_u].sequence.len();
                let mut total_kmers = unitigs[first_u].kmers_count;

                for step_idx in 0..path_reps.len() {
                    for &(rep_u, rep_rev) in &path_reps[step_idx] {
                        let r_seq = if rep_rev {
                            revcomp_bytes(&unitigs[rep_u].sequence)
                        } else {
                            unitigs[rep_u].sequence.clone()
                        };
                        seq.extend_from_slice(b"NNNNNNNNNNNNNNNNNNNN");
                        seq.extend_from_slice(&r_seq);
                        total_bp += r_seq.len() + 20;
                        cov_sum += unitigs[rep_u].mean_coverage * r_seq.len() as f64;
                        total_kmers += unitigs[rep_u].kmers_count;
                    }

                    let (nxt_u, nxt_rev) = path_unitigs[step_idx + 1];
                    let n_seq = if nxt_rev {
                        revcomp_bytes(&unitigs[nxt_u].sequence)
                    } else {
                        unitigs[nxt_u].sequence.clone()
                    };

                    let max_ov = k1.min(seq.len()).min(n_seq.len());
                    let mut ov_found = 0;
                    if max_ov >= min_ov {
                        for ov in (min_ov..=max_ov).rev() {
                            if seq[seq.len() - ov..] == n_seq[..ov] {
                                ov_found = ov;
                                break;
                            }
                        }
                    }

                    if ov_found >= min_ov {
                        seq.extend_from_slice(&n_seq[ov_found..]);
                        total_bp += n_seq.len() - ov_found;
                    } else {
                        seq.extend_from_slice(b"NNNNNNNNNNNNNNNNNNNN");
                        seq.extend_from_slice(&n_seq);
                        total_bp += n_seq.len() + 20;
                    }
                    cov_sum += unitigs[nxt_u].mean_coverage * n_seq.len() as f64;
                    total_kmers += unitigs[nxt_u].kmers_count;
                }

                stitched_contigs.push(Unitig {
                    id: stitched_contigs.len(),
                    sequence: seq,
                    mean_coverage: if total_bp > 0 {
                        cov_sum / total_bp as f64
                    } else {
                        median_cov
                    },
                    kmers_count: total_kmers,
                });
            }
        }

        let mut final_unitigs = stitched_contigs;
        for (i, u) in unitigs.into_iter().enumerate() {
            if !visited_unitigs[i] {
                final_unitigs.push(u);
            }
        }

        final_unitigs.sort_by_key(|u| std::cmp::Reverse(u.sequence.len()));
        final_unitigs
    }
}
