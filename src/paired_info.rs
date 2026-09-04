//! Paired-End Read Distance Estimator & Graph Mapping Engine (Module 3.1).
//!
//! Maps paired-end reads to assembly graph unitigs to estimate library insert size
//! and construct the paired-end transition matrix for repeat resolution.

use crate::dna::{canonical_kmer_u64, string_to_kmer};
use crate::graph::Unitig;
use hashbrown::HashMap;

/// A paired-end connection between two unitigs supported by read pairs.
#[derive(Clone, Debug, Default)]
pub struct PairedLink {
    pub target_unitig: usize,
    pub support_count: u32,
    pub total_distance: f64,
}

/// Index of paired-end connections between unitigs.
#[derive(Default, Debug)]
pub struct PairedInfoIndex {
    pub mean_insert_size: f64,
    pub insert_size_stdev: f64,
    // Maps source_unitig_id -> (target_unitig_id -> PairedLink)
    pub links: HashMap<usize, HashMap<usize, PairedLink>>,
}

impl PairedInfoIndex {
    /// Builds paired-end connections from paired reads and unitigs.
    pub fn build(
        k: usize,
        unitigs: &[Unitig],
        reads1: &[Vec<u8>],
        reads2: &[Vec<u8>],
    ) -> Self {
        // 1. Build a k-mer to unitig lookup table
        // kmer -> (unitig_id, offset)
        let mut kmer_to_unitig: HashMap<u64, (usize, usize)> = HashMap::with_capacity(unitigs.len() * 100);
        for (u_id, u) in unitigs.iter().enumerate() {
            if u.sequence.len() < k {
                continue;
            }
            for i in 0..=(u.sequence.len() - k) {
                if let Some(km) = string_to_kmer(&u.sequence[i..i + k], k) {
                    let (can, _) = canonical_kmer_u64(km, k);
                    kmer_to_unitig.entry(can).or_insert((u_id, i));
                }
            }
        }

        // Helper to map a read to a unitig
        let map_read = |seq: &[u8]| -> Option<(usize, usize)> {
            if seq.len() < k {
                return None;
            }
            // Sample a few k-mers along the read
            let step = (seq.len() - k).max(1) / 3;
            for i in (0..=(seq.len() - k)).step_by(step.max(1)) {
                if let Some(km) = string_to_kmer(&seq[i..i + k], k) {
                    let (can, _) = canonical_kmer_u64(km, k);
                    if let Some(&loc) = kmer_to_unitig.get(&can) {
                        return Some(loc);
                    }
                }
            }
            None
        };

        let mut internal_distances: Vec<f64> = Vec::new();
        let mut links: HashMap<usize, HashMap<usize, PairedLink>> = HashMap::new();

        let num_pairs = reads1.len().min(reads2.len());
        for i in 0..num_pairs {
            let loc1 = map_read(&reads1[i]);
            let loc2 = map_read(&reads2[i]);

            match (loc1, loc2) {
                (Some((u1, pos1)), Some((u2, pos2))) => {
                    if u1 == u2 {
                        let dist = (pos1 as f64 - pos2 as f64).abs();
                        if dist > 0.0 && dist < 2000.0 {
                            internal_distances.push(dist);
                        }
                    } else {
                        // Read 1 on u1, Read 2 on u2 -> paired support between u1 and u2
                        let entry1 = links.entry(u1).or_default();
                        let link1 = entry1.entry(u2).or_insert(PairedLink {
                            target_unitig: u2,
                            support_count: 0,
                            total_distance: 0.0,
                        });
                        link1.support_count += 1;

                        // Symmetric link
                        let entry2 = links.entry(u2).or_default();
                        let link2 = entry2.entry(u1).or_insert(PairedLink {
                            target_unitig: u1,
                            support_count: 0,
                            total_distance: 0.0,
                        });
                        link2.support_count += 1;
                    }
                }
                _ => {}
            }
        }

        // Calculate insert size distribution
        let (mean_insert, stdev) = if !internal_distances.is_empty() {
            let sum: f64 = internal_distances.iter().sum();
            let mean = sum / internal_distances.len() as f64;
            let variance: f64 = internal_distances.iter().map(|&d| (d - mean).powi(2)).sum::<f64>()
                / internal_distances.len() as f64;
            (mean, variance.sqrt())
        } else {
            (300.0, 50.0) // Default typical Illumina insert size
        };

        Self {
            mean_insert_size: mean_insert,
            insert_size_stdev: stdev,
            links,
        }
    }

    /// Queries paired-end support count between unitig u1 and unitig u2.
    pub fn get_support(&self, u1: usize, u2: usize) -> u32 {
        self.links
            .get(&u1)
            .and_then(|targets| targets.get(&u2))
            .map(|l| l.support_count)
            .unwrap_or(0)
    }
}
