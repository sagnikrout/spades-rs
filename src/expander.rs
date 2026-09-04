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

        // Build adjacency from sequence overlaps
        let n = unitigs.len();
        // Maps unitig_i -> list of candidate successors (unitig_j, is_rc)
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];

        for i in 0..n {
            if unitigs[i].sequence.len() < k1 {
                continue;
            }
            let suffix_i = &unitigs[i].sequence[unitigs[i].sequence.len() - k1..];

            for j in 0..n {
                if i == j || unitigs[j].sequence.len() < k1 {
                    continue;
                }
                if suffix_i == &unitigs[j].sequence[..k1] {
                    adj[i].push(j);
                }
            }
        }

        // If no branches exist, no repeat resolution is needed
        let has_branches = adj.iter().any(|succs| succs.len() > 1);
        if !has_branches {
            return unitigs;
        }

        println!("  [ExSPAnder] Resolving repeat junctions with paired-end support...");

        // Resolve branches using paired-end linkages
        let mut resolved_adj: Vec<Option<usize>> = vec![None; n];

        for i in 0..n {
            if adj[i].len() == 1 {
                resolved_adj[i] = Some(adj[i][0]);
            } else if adj[i].len() > 1 {
                // Repeat knot! Evaluate paired-end support for each candidate
                let mut best_target = None;
                let mut best_support = 0u32;
                let mut second_best_support = 0u32;

                for &target in &adj[i] {
                    let support = paired_info.get_support(i, target);
                    if support > best_support {
                        second_best_support = best_support;
                        best_support = support;
                        best_target = Some(target);
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
                            "    Resolved repeat branch {} -> {:?} (Support: {}, Ratio: {:.1}x)",
                            i, best_target, best_support, ratio
                        );
                        resolved_adj[i] = best_target;
                    }
                }
            }
        }

        // Traverse resolved paths and construct new extended contigs
        let mut visited = vec![false; n];
        let mut extended_contigs = Vec::new();
        let mut new_id = 0;

        for start_node in 0..n {
            if visited[start_node] {
                continue;
            }

            let mut path = vec![start_node];
            visited[start_node] = true;
            let mut curr = start_node;

            while let Some(next_node) = resolved_adj[curr] {
                if visited[next_node] {
                    break; // Prevent cycles
                }
                visited[next_node] = true;
                path.push(next_node);
                curr = next_node;
            }

            // Assemble path sequence
            let mut seq = unitigs[path[0]].sequence.clone();
            let mut total_cov = unitigs[path[0]].mean_coverage * unitigs[path[0]].sequence.len() as f64;
            let mut total_len = unitigs[path[0]].sequence.len();

            for &node in &path[1..] {
                seq.extend_from_slice(&unitigs[node].sequence[k1..]);
                total_cov += unitigs[node].mean_coverage * (unitigs[node].sequence.len() - k1) as f64;
                total_len += unitigs[node].sequence.len() - k1;
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
        extended_contigs.sort_by(|a, b| b.sequence.len().cmp(&a.sequence.len()));

        extended_contigs
    }
}
