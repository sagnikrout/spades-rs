//! Bidirected Compacted de Bruijn Graph (cDBG) Engine.
//!
//! Models double-stranded DNA where each vertex is a canonical k-mer with two orientations:
//! - Forward (false): sequence = kmer
//! - Reverse (true): sequence = revcomp(kmer)
//!
//! Compacts non-branching paths into maximal Unitigs.

use crate::dna::Kmer256;
use hashbrown::{HashMap, HashSet};

/// An oriented k-mer in the bidirected de Bruijn graph.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct OrientedKmer {
    /// Canonical (lexicographically smaller of forward and reverse-complement) packed $k$-mer.
    pub kmer: Kmer256,
    /// Orientation flag: `false` for forward (`kmer`), `true` for reverse-complement (`kmer.revcomp(k)`).
    pub is_rc: bool,
}

impl OrientedKmer {
    /// Returns the explicit 2-bit $k$-mer sequence in its current orientation.
    #[inline(always)]
    pub fn explicit_seq(&self, k: usize) -> Kmer256 {
        if !self.is_rc {
            self.kmer
        } else {
            self.kmer.revcomp(k)
        }
    }

    /// Returns the opposite-strand twin of this oriented $k$-mer (same canonical $k$-mer with flipped `is_rc`).
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
    /// Unique numeric identifier of the unitig or contig within the current graph.
    pub id: usize,
    /// Spelled ASCII nucleotide sequence (`A`, `C`, `G`, `T`, or `N` in scaffold gaps).
    pub sequence: Vec<u8>,
    /// Mean $k$-mer multiplicity (coverage depth) across the unitig.
    pub mean_coverage: f64,
    /// Total number of $k$-mers compacted into this unitig.
    pub kmers_count: usize,
}

/// Bidirected compacted de Bruijn graph represented as a collection of maximal [`Unitig`]s.
pub struct CompactedGraph {
    /// $k$-mer size used to construct the graph.
    pub k: usize,
    /// Maximal non-branching unitigs extracted from the graph.
    pub unitigs: Vec<Unitig>,
}

impl CompactedGraph {
    /// Builds maximal unitigs from the solid canonical k-mers.
    pub fn build(
        k: usize,
        solid_kmers: &HashSet<Kmer256>,
        kmer_cov: &HashMap<Kmer256, u32>,
    ) -> Self {
        // Helper closures for getting successors and predecessors
        let get_successors = |ok: OrientedKmer| -> Vec<(OrientedKmer, u8)> {
            let s = ok.explicit_seq(k);
            let mut succs = Vec::with_capacity(4);
            for b in 0..4u8 {
                let next_val = s.extend_right(b, k);
                let (can, is_rc) = next_val.canonical(k);
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
                let prev_val = s.prepend_left(a, k);
                let (can, is_rc) = prev_val.canonical(k);
                if solid_kmers.contains(&can) {
                    preds.push((OrientedKmer { kmer: can, is_rc }, a));
                }
            }
            preds
        };

        let extract_maximal_path =
            |start_ok: OrientedKmer, visited: &mut HashSet<Kmer256>| -> Vec<OrientedKmer> {
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
                    return path;
                }

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

                backward_path.reverse();
                backward_path.extend(path);
                backward_path
            };

        let build_unitig = |id: usize, path: &[OrientedKmer]| -> Unitig {
            let first_seq = path[0].explicit_seq(k);
            let mut unitig_bytes = first_seq.to_string(k).into_bytes();
            let mut total_cov = *kmer_cov.get(&path[0].kmer).unwrap_or(&2) as f64;

            for ok in &path[1..] {
                let explicit = ok.explicit_seq(k);
                let last_base = crate::dna::bit2_to_base(explicit.last_base(k));
                unitig_bytes.push(last_base);
                total_cov += *kmer_cov.get(&ok.kmer).unwrap_or(&2) as f64;
            }

            Unitig {
                id,
                sequence: unitig_bytes,
                mean_coverage: total_cov / path.len() as f64,
                kmers_count: path.len(),
            }
        };

        let mut visited = HashSet::<Kmer256>::with_capacity(solid_kmers.len());
        let mut unitigs = Vec::new();
        let mut unitig_id = 0;

        // Pass 1: Traverse paths starting preferentially from boundary nodes
        for &kmer in solid_kmers {
            if visited.contains(&kmer) {
                continue;
            }

            let mut start_ok = OrientedKmer { kmer, is_rc: false };
            let preds_fwd = get_predecessors(start_ok);
            let preds_rev = get_predecessors(start_ok.twin());

            if preds_fwd.len() == 1 {
                let (pred_node, _) = preds_fwd[0];
                let succs_of_pred = get_successors(pred_node);
                if succs_of_pred.len() == 1 && !visited.contains(&pred_node.kmer) {
                    if preds_rev.len() == 1 {
                        let (pred_rev_node, _) = preds_rev[0];
                        let succs_of_rev_pred = get_successors(pred_rev_node);
                        if succs_of_rev_pred.len() == 1 && !visited.contains(&pred_rev_node.kmer) {
                            continue; // Interior node of an unvisited path; defer to boundary or Pass 2
                        }
                    } else {
                        start_ok = start_ok.twin();
                    }
                }
            }

            let full_path = extract_maximal_path(start_ok, &mut visited);
            if !full_path.is_empty() {
                unitigs.push(build_unitig(unitig_id, &full_path));
                unitig_id += 1;
            }
        }

        // Pass 2: Recover isolated cycles and any remaining unvisited paths
        for &kmer in solid_kmers {
            if visited.contains(&kmer) {
                continue;
            }

            let start_ok = OrientedKmer { kmer, is_rc: false };
            let full_path = extract_maximal_path(start_ok, &mut visited);
            if !full_path.is_empty() {
                unitigs.push(build_unitig(unitig_id, &full_path));
                unitig_id += 1;
            }
        }

        Self { k, unitigs }
    }
}
