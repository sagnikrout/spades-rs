//! ExSPAnder (Exact SPath Assembler) Repeat Resolution Engine (Module 3.2).
//!
//! Resolves repetitive graph knots using paired-end read distance constraints.

use crate::graph::Unitig;
use crate::paired_info::PairedInfoIndex;

pub struct ExSPAnder {
    pub min_support: u32,
    pub confidence_ratio: f64,
}

impl Default for ExSPAnder {
    fn default() -> Self {
        Self {
            min_support: 3,
            confidence_ratio: 2.5,
        }
    }
}

impl ExSPAnder {
    /// Resolves repeat tangles in the unitig set using paired-end read linkages.
    pub fn resolve_repeats(
        &self,
        k: usize,
        unitigs: Vec<Unitig>,
        paired_info: &PairedInfoIndex,
    ) -> Vec<Unitig> {
        let k1 = k - 1;
        if unitigs.len() <= 1 {
            return unitigs;
        }

        let n = unitigs.len();

        // 1. Build fast prefix index for candidate successors in both orientations:
        // A successor j can match suffix_i either in forward orientation (prefix_j)
        // or reverse-complement orientation (rc(suffix_j))
        let mut prefix_map: hashbrown::HashMap<Vec<u8>, Vec<(usize, bool)>> =
            hashbrown::HashMap::with_capacity(n * 2);

        for (j, u_j) in unitigs.iter().enumerate() {
            if u_j.sequence.len() < k1 {
                continue;
            }
            // Forward orientation: prefix_j
            let p_fwd = u_j.sequence[..k1].to_vec();
            prefix_map.entry(p_fwd).or_default().push((j, false));

            // Reverse complement orientation: rc(suffix_j)
            let p_rc = crate::dna::revcomp_bytes(&u_j.sequence[u_j.sequence.len() - k1..]);
            prefix_map.entry(p_rc).or_default().push((j, true));
        }

        // 2. Maps unitig_i -> list of candidate successors (unitig_j, is_rc)
        let mut adj: Vec<Vec<(usize, bool)>> = vec![Vec::new(); n];

        for i in 0..n {
            if unitigs[i].sequence.len() < k1 {
                continue;
            }
            let suffix_i = unitigs[i].sequence[unitigs[i].sequence.len() - k1..].to_vec();

            if let Some(candidates) = prefix_map.get(&suffix_i) {
                for &(j, is_rc) in candidates {
                    if i != j {
                        adj[i].push((j, is_rc));
                    }
                }
            }
        }

        // If no branches exist, no repeat resolution is needed
        let has_branches = adj.iter().any(|succs| succs.len() > 1);
        if !has_branches {
            return unitigs;
        }

        println!("  [ExSPAnder] Resolving bidirected repeat junctions with paired-end support...");

        // Resolve branches using paired-end linkages
        let mut resolved_adj: Vec<Option<(usize, bool)>> = vec![None; n];

        for i in 0..n {
            if adj[i].len() == 1 {
                resolved_adj[i] = Some(adj[i][0]);
            } else if adj[i].len() > 1 {
                // Repeat knot! Evaluate paired-end support for each candidate
                let mut best_target = None;
                let mut best_is_rc = false;
                let mut best_support = 0u32;
                let mut second_best_support = 0u32;

                for &(target, is_rc) in &adj[i] {
                    let support = paired_info.get_support(i, target);
                    if support > best_support {
                        second_best_support = best_support;
                        best_support = support;
                        best_target = Some(target);
                        best_is_rc = is_rc;
                    } else if support > second_best_support {
                        second_best_support = support;
                    }
                }

                // Check ExSPAnder confidence criterion
                if best_support >= self.min_support {
                    let ratio = if second_best_support > 0 {
                        best_support as f64 / second_best_support as f64
                    } else {
                        best_support as f64
                    };

                    if ratio >= self.confidence_ratio {
                        println!(
                            "    Resolved repeat branch {} -> {:?} (RC: {}, Support: {}, Ratio: {:.1}x)",
                            i, best_target, best_is_rc, best_support, ratio
                        );
                        resolved_adj[i] = best_target.map(|t| (t, best_is_rc));
                    }
                }
            }
        }

        // Disambiguate incoming repeat collisions:
        // If multiple unitigs point to the same target, evaluate paired-end support
        let mut in_candidates: Vec<Vec<(usize, bool)>> = vec![Vec::new(); n];
        for (i, &adj) in resolved_adj.iter().enumerate() {
            if let Some((next_node, is_rc)) = adj {
                in_candidates[next_node].push((i, is_rc));
            }
        }

        for (target, candidates) in in_candidates.iter().enumerate() {
            if candidates.len() > 1 {
                let mut best_src = None;
                let mut best_sup = 0u32;
                let mut second_sup = 0u32;

                for &(src, _) in candidates {
                    let sup = paired_info.get_support(src, target);
                    if sup > best_sup {
                        second_sup = best_sup;
                        best_sup = sup;
                        best_src = Some(src);
                    } else if sup > second_sup {
                        second_sup = sup;
                    }
                }

                let winner = if best_sup >= self.min_support
                    && (second_sup == 0
                        || best_sup as f64 / second_sup as f64 >= self.confidence_ratio)
                {
                    best_src
                } else {
                    None
                };

                for &(src, _) in candidates {
                    if Some(src) != winner {
                        resolved_adj[src] = None;
                    }
                }
            }
        }

        // Compute final in-degrees to identify true path sources
        let mut final_in_degree = vec![0usize; n];
        for &adj in &resolved_adj {
            if let Some((next_node, _)) = adj {
                final_in_degree[next_node] += 1;
            }
        }

        // Traverse resolved paths from true source nodes (in-degree == 0) first
        let mut start_order: Vec<usize> = (0..n).collect();
        start_order.sort_unstable_by_key(|&idx| final_in_degree[idx]);

        let mut visited = vec![false; n];
        let mut extended_contigs = Vec::new();
        let mut new_id = 0;

        for start_node in start_order {
            if visited[start_node] {
                continue;
            }

            let mut path: Vec<(usize, bool)> = vec![(start_node, false)];
            visited[start_node] = true;
            let mut curr = start_node;

            while let Some((next_node, is_rc)) = resolved_adj[curr] {
                if visited[next_node] {
                    break; // Prevent cycles
                }
                visited[next_node] = true;
                path.push((next_node, is_rc));
                curr = next_node;
            }

            // Assemble path sequence
            let mut seq = unitigs[path[0].0].sequence.clone();
            let mut total_cov =
                unitigs[path[0].0].mean_coverage * unitigs[path[0].0].sequence.len() as f64;
            let mut total_len = unitigs[path[0].0].sequence.len();

            for &(node, is_rc) in &path[1..] {
                let node_seq = if !is_rc {
                    unitigs[node].sequence.clone()
                } else {
                    crate::dna::revcomp_bytes(&unitigs[node].sequence)
                };

                if node_seq.len() > k1 {
                    seq.extend_from_slice(&node_seq[k1..]);
                    total_cov += unitigs[node].mean_coverage * (node_seq.len() - k1) as f64;
                    total_len += node_seq.len() - k1;
                }
            }

            extended_contigs.push(Unitig {
                id: new_id,
                sequence: seq,
                mean_coverage: total_cov / total_len as f64,
                kmers_count: total_len.saturating_sub(k1),
            });
            new_id += 1;
        }

        // Sort descending by length
        extended_contigs.sort_by_key(|a| std::cmp::Reverse(a.sequence.len()));

        extended_contigs
    }
}
