//! Hybrid Long-Read Resolver (Module 3.3: Spaligner).
//!
//! Maps PacBio HiFi and Oxford Nanopore (ONT) long reads across assembly graph unitigs
//! using seed chaining to bridge multi-kilobase genomic repeats.

use crate::dna::{canonical_kmer_u64, string_to_kmer};
use crate::graph::Unitig;
use hashbrown::HashMap;

pub struct LongReadResolver {
    pub k: usize,
    pub min_seed_matches: usize,
}

impl Default for LongReadResolver {
    fn default() -> Self {
        Self {
            k: 21,
            min_seed_matches: 3,
        }
    }
}

impl LongReadResolver {
    /// Bridges unitigs using long reads.
    pub fn bridge_with_long_reads(
        &self,
        mut unitigs: Vec<Unitig>,
        long_reads: &[Vec<u8>],
    ) -> Vec<Unitig> {
        if long_reads.is_empty() || unitigs.len() <= 1 {
            return unitigs;
        }

        let k = self.k;
        let k1 = k - 1;

        // Build k-mer to unitig index
        let mut kmer_to_unitig: HashMap<u64, usize> = HashMap::with_capacity(unitigs.len() * 50);
        for (idx, u) in unitigs.iter().enumerate() {
            if u.sequence.len() < k {
                continue;
            }
            for i in (0..=(u.sequence.len() - k)).step_by(k / 2 + 1) {
                if let Some(km) = string_to_kmer(&u.sequence[i..i + k], k) {
                    let (can, _) = canonical_kmer_u64(km, k);
                    kmer_to_unitig.insert(can, idx);
                }
            }
        }

        // Trace ordered unitigs along each long read
        let mut bridge_counts: HashMap<(usize, usize), usize> = HashMap::new();

        for read in long_reads {
            if read.len() < k * 2 {
                continue;
            }

            let mut ordered_hits: Vec<usize> = Vec::new();
            let step = (k / 2).max(1);

            for i in (0..=(read.len() - k)).step_by(step) {
                if let Some(km) = string_to_kmer(&read[i..i + k], k) {
                    let (can, _) = canonical_kmer_u64(km, k);
                    if let Some(&u_id) = kmer_to_unitig.get(&can) {
                        if ordered_hits.last() != Some(&u_id) {
                            ordered_hits.push(u_id);
                        }
                    }
                }
            }

            // Record transitions
            for window in ordered_hits.windows(2) {
                let u1 = window[0];
                let u2 = window[1];
                if u1 != u2 {
                    *bridge_counts.entry((u1, u2)).or_insert(0) += 1;
                }
            }
        }

        println!("  [Spaligner] Evaluated {} long reads, identified {} unitig bridge candidates", long_reads.len(), bridge_counts.len());

        // Merge unitigs supported by >= min_seed_matches long reads
        let mut merged = true;
        while merged {
            merged = false;
            let mut best_pair = None;
            let mut max_support = 0;

            for (&(u1, u2), &support) in &bridge_counts {
                if support >= self.min_seed_matches && support > max_support {
                    if u1 < unitigs.len() && u2 < unitigs.len() && u1 != u2 {
                        best_pair = Some((u1, u2));
                        max_support = support;
                    }
                }
            }

            if let Some((u1, u2)) = best_pair {
                // If u1 suffix overlaps u2 prefix by k1, stitch directly
                let seq2 = &unitigs[u2].sequence;
                let overlap = if unitigs[u1].sequence.len() >= k1 && seq2.len() >= k1 {
                    let s1 = &unitigs[u1].sequence[unitigs[u1].sequence.len() - k1..];
                    if s1 == &seq2[..k1] {
                        k1
                    } else {
                        0
                    }
                } else {
                    0
                };

                let tail = unitigs[u2].sequence[overlap..].to_vec();
                unitigs[u1].sequence.extend_from_slice(&tail);
                unitigs.remove(u2);
                bridge_counts.remove(&(u1, u2));
                merged = true;
            }
        }

        unitigs
    }
}
