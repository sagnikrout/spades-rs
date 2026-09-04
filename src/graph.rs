//! Bidirected Compacted de Bruijn Graph (cDBG) Engine.
//!
//! Models double-stranded DNA where each vertex is a canonical k-mer with two orientations:
//! - Forward (false): sequence = kmer
//! - Reverse (true): sequence = revcomp(kmer)
//!
//! Compacts non-branching paths into maximal Unitigs.

use crate::dna::{canonical_kmer_u64, revcomp_kmer_u64};
use hashbrown::{HashMap, HashSet};

/// An oriented k-mer in the bidirected de Bruijn graph.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct OrientedKmer {
    pub kmer: u64,     // Canonical k-mer
    pub is_rc: bool,    // false = Forward, true = Reverse Complement
}

impl OrientedKmer {
    #[inline(always)]
    pub fn explicit_seq(&self, k: usize) -> u64 {
        if !self.is_rc {
            self.kmer
        } else {
            revcomp_kmer_u64(self.kmer, k)
        }
    }

    #[inline(always)]
    pub fn twin(&self) -> Self {
        Self {
            kmer: self.kmer,
            is_rc: !self.is_rc,
        }
    }
}

/// A maximal non-branching path of oriented k-mers.
#[derive(Clone, Debug)]
pub struct Unitig {
    pub id: usize,
    pub sequence: Vec<u8>,
    pub mean_coverage: f64,
    pub kmers_count: usize,
}

pub struct CompactedGraph {
    pub k: usize,
    pub unitigs: Vec<Unitig>,
}

impl CompactedGraph {
    /// Builds maximal unitigs from the solid canonical k-mers.
    pub fn build(
        k: usize,
        solid_kmers: &HashSet<u64>,
        kmer_cov: &HashMap<u64, u32>,
    ) -> Self {
        let kmer_mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };

        // Helper closures for getting successors and predecessors
        let get_successors = |ok: OrientedKmer| -> Vec<(OrientedKmer, u8)> {
            let s = ok.explicit_seq(k);
            let mut succs = Vec::with_capacity(4);
            for b in 0..4u8 {
                let next_val = ((s << 2) | (b as u64)) & kmer_mask;
                let (can, is_rc) = canonical_kmer_u64(next_val, k);
                if solid_kmers.contains(&can) {
                    succs.push((OrientedKmer { kmer: can, is_rc }, b));
                }
            }
            succs
        };

        let get_predecessors = |ok: OrientedKmer| -> Vec<(OrientedKmer, u8)> {
            let s = ok.explicit_seq(k);
            let mut preds = Vec::with_capacity(4);
            for a in 0..4u8 {
                let prev_val = (s >> 2) | ((a as u64) << (2 * (k - 1)));
                let (can, is_rc) = canonical_kmer_u64(prev_val, k);
                if solid_kmers.contains(&can) {
                    preds.push((OrientedKmer { kmer: can, is_rc }, a));
                }
            }
            preds
        };

        let mut visited = HashSet::<u64>::with_capacity(solid_kmers.len());
        let mut unitigs = Vec::new();
        let mut unitig_id = 0;

        for &kmer in solid_kmers {
            if visited.contains(&kmer) {
                continue;
            }

            // Test if either orientation is a unitig start (i.e. in-degree != 1 or cycle)
            let mut start_ok = OrientedKmer { kmer, is_rc: false };
            let preds_fwd = get_predecessors(start_ok);
            let preds_rev = get_predecessors(start_ok.twin());

            // If start_ok has pred == 1, but its predecessor has out-degree == 1, we shouldn't start here
            // unless we can't find a boundary
            if preds_fwd.len() == 1 {
                let (pred_node, _) = preds_fwd[0];
                let succs_of_pred = get_successors(pred_node);
                if succs_of_pred.len() == 1 && !visited.contains(&pred_node.kmer) {
                    // There is an unvisited predecessor that leads uniquely to this node
                    // Check if reverse orientation is better
                    if preds_rev.len() == 1 {
                        let (pred_rev_node, _) = preds_rev[0];
                        let succs_of_rev_pred = get_successors(pred_rev_node);
                        if succs_of_rev_pred.len() == 1 && !visited.contains(&pred_rev_node.kmer) {
                            continue; // In the middle of an unvisited path in both directions
                        }
                    } else {
                        start_ok = start_ok.twin();
                    }
                }
            }

            // Traverse forward from start_ok
            let mut path: Vec<OrientedKmer> = Vec::new();
            let mut curr = start_ok;

            while !visited.contains(&curr.kmer) {
                visited.insert(curr.kmer);
                path.push(curr);

                let succs = get_successors(curr);
                if succs.len() != 1 {
                    break; // Branch point or dead end
                }

                let (next_ok, _) = succs[0];
                if visited.contains(&next_ok.kmer) {
                    break; // Cycle or visited
                }

                let preds_of_next = get_predecessors(next_ok);
                if preds_of_next.len() != 1 {
                    break; // Merge point
                }

                curr = next_ok;
            }

            if path.is_empty() {
                continue;
            }

            // Also check if we can extend backwards from start_ok if in-degree == 1
            // (e.g. if we started in the middle of a path)
            let mut backward_path: Vec<OrientedKmer> = Vec::new();
            let mut back_curr = path[0];
            loop {
                let preds = get_predecessors(back_curr);
                if preds.len() != 1 {
                    break;
                }
                let (prev_ok, _) = preds[0];
                if visited.contains(&prev_ok.kmer) {
                    break;
                }
                let succs_of_prev = get_successors(prev_ok);
                if succs_of_prev.len() != 1 {
                    break;
                }
                visited.insert(prev_ok.kmer);
                backward_path.push(prev_ok);
                back_curr = prev_ok;
            }

            // Stitch backward_path (reversed) + path
            backward_path.reverse();
            backward_path.extend(path);
            let full_path = backward_path;

            // Reconstruct full unitig sequence
            let first_seq = full_path[0].explicit_seq(k);
            let mut unitig_bytes = crate::dna::kmer_to_string(first_seq, k).into_bytes();
            let mut total_cov = *kmer_cov.get(&full_path[0].kmer).unwrap_or(&1) as f64;

            for ok in &full_path[1..] {
                let explicit = ok.explicit_seq(k);
                let last_base = crate::dna::bit2_to_base((explicit & 3) as u8);
                unitig_bytes.push(last_base);
                total_cov += *kmer_cov.get(&ok.kmer).unwrap_or(&1) as f64;
            }

            let mean_cov = total_cov / full_path.len() as f64;

            unitigs.push(Unitig {
                id: unitig_id,
                sequence: unitig_bytes,
                mean_coverage: mean_cov,
                kmers_count: full_path.len(),
            });
            unitig_id += 1;
        }

        Self { k, unitigs }
    }
}
