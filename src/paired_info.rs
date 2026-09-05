//! Paired-End Read Distance Estimator & Graph Mapping Engine (Module 3.1).
//!
//! Maps paired-end reads to assembly graph unitigs to estimate library insert size
//! and construct the paired-end transition matrix for repeat resolution.

use crate::dna::Kmer256;
use crate::graph::Unitig;
use hashbrown::HashMap;
use rayon::prelude::*;

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
    pub fn build(k: usize, unitigs: &[Unitig], reads1: &[Vec<u8>], reads2: &[Vec<u8>]) -> Self {
        // 1. Build a k-mer to unitig lookup table
        // kmer -> (unitig_id, offset)
        let mut kmer_to_unitig: HashMap<Kmer256, (usize, usize)> =
            HashMap::with_capacity(unitigs.len() * 100);
        for (u_id, u) in unitigs.iter().enumerate() {
            if u.sequence.len() < k {
                continue;
            }
            for i in 0..=(u.sequence.len() - k) {
                if let Some(km) = Kmer256::from_bytes(&u.sequence[i..i + k], k) {
                    let (can, _) = km.canonical(k);
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
                if let Some(km) = Kmer256::from_bytes(&seq[i..i + k], k) {
                    let (can, _) = km.canonical(k);
                    if let Some(&loc) = kmer_to_unitig.get(&can) {
                        return Some(loc);
                    }
                }
            }
            None
        };

        let num_pairs = reads1.len().min(reads2.len());
        let pair_indices: Vec<usize> = (0..num_pairs).collect();

        // 2. Map paired reads in parallel chunks
        let thread_results: Vec<(Vec<f64>, HashMap<(usize, usize), u32>)> = pair_indices
            .par_chunks(4096)
            .map(|chunk| {
                let mut local_dists = Vec::new();
                let mut local_links: HashMap<(usize, usize), u32> = HashMap::new();
                for &i in chunk {
                    let loc1 = map_read(&reads1[i]);
                    let loc2 = map_read(&reads2[i]);

                    if let (Some((u1, pos1)), Some((u2, pos2))) = (loc1, loc2) {
                        if u1 == u2 {
                            let dist = (pos1 as f64 - pos2 as f64).abs();
                            if dist > 0.0 && dist < 2000.0 {
                                local_dists.push(dist);
                            }
                        } else {
                            *local_links.entry((u1, u2)).or_insert(0) += 1;
                            *local_links.entry((u2, u1)).or_insert(0) += 1;
                        }
                    }
                }
                (local_dists, local_links)
            })
            .collect();

        // 3. Aggregate thread results
        let mut internal_distances: Vec<f64> = Vec::new();
        let mut links: HashMap<usize, HashMap<usize, PairedLink>> = HashMap::new();

        for (dists, l_map) in thread_results {
            internal_distances.extend(dists);
            for ((u1, u2), cnt) in l_map {
                let entry = links.entry(u1).or_default();
                let link = entry.entry(u2).or_insert(PairedLink {
                    target_unitig: u2,
                    support_count: 0,
                    total_distance: 0.0,
                });
                link.support_count += cnt;
            }
        }

        // Calculate insert size distribution
        let (mean_insert, stdev) = if !internal_distances.is_empty() {
            let sum: f64 = internal_distances.iter().sum();
            let mean = sum / internal_distances.len() as f64;
            let variance: f64 = internal_distances
                .iter()
                .map(|&d| (d - mean).powi(2))
                .sum::<f64>()
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

    /// Builds paired-end connections directly from PackedReads without unpacking vectors.
    pub fn build_from_packed(
        k: usize,
        unitigs: &[Unitig],
        packed: &crate::packed_reads::PackedReads,
    ) -> Self {
        let mut kmer_to_unitig: HashMap<Kmer256, (usize, usize)> =
            HashMap::with_capacity(unitigs.len() * 100);
        for (u_id, u) in unitigs.iter().enumerate() {
            if u.sequence.len() < k {
                continue;
            }
            for i in 0..=(u.sequence.len() - k) {
                if let Some(km) = Kmer256::from_bytes(&u.sequence[i..i + k], k) {
                    let (can, _) = km.canonical(k);
                    kmer_to_unitig.entry(can).or_insert((u_id, i));
                }
            }
        }

        let map_read = |seq: &[u8]| -> Option<(usize, usize)> {
            if seq.len() < k {
                return None;
            }
            let step = (seq.len() - k).max(1) / 3;
            for i in (0..=(seq.len() - k)).step_by(step.max(1)) {
                if let Some(km) = Kmer256::from_bytes(&seq[i..i + k], k) {
                    let (can, _) = km.canonical(k);
                    if let Some(&loc) = kmer_to_unitig.get(&can) {
                        return Some(loc);
                    }
                }
            }
            None
        };

        let n1 = packed.pe_boundary;
        let n2 = packed.len().saturating_sub(n1);
        let num_pairs = n1.min(n2);
        let pair_indices: Vec<usize> = (0..num_pairs).collect();

        let thread_results: Vec<(Vec<f64>, HashMap<(usize, usize), u32>)> = pair_indices
            .par_chunks(4096)
            .map(|chunk| {
                let mut buf1 = Vec::with_capacity(512);
                let mut buf2 = Vec::with_capacity(512);
                let mut local_dists = Vec::new();
                let mut local_links: HashMap<(usize, usize), u32> = HashMap::new();

                for &i in chunk {
                    packed.get_read(i, &mut buf1);
                    packed.get_read(n1 + i, &mut buf2);
                    let loc1 = map_read(&buf1);
                    let loc2 = map_read(&buf2);

                    if let (Some((u1, pos1)), Some((u2, pos2))) = (loc1, loc2) {
                        if u1 == u2 {
                            let dist = (pos1 as f64 - pos2 as f64).abs();
                            if dist > 0.0 && dist < 2000.0 {
                                local_dists.push(dist);
                            }
                        } else {
                            *local_links.entry((u1, u2)).or_insert(0) += 1;
                            *local_links.entry((u2, u1)).or_insert(0) += 1;
                        }
                    }
                }
                (local_dists, local_links)
            })
            .collect();

        let mut internal_distances: Vec<f64> = Vec::new();
        let mut links: HashMap<usize, HashMap<usize, PairedLink>> = HashMap::new();

        for (dists, l_map) in thread_results {
            internal_distances.extend(dists);
            for ((u1, u2), cnt) in l_map {
                let entry = links.entry(u1).or_default();
                let link = entry.entry(u2).or_insert(PairedLink {
                    target_unitig: u2,
                    support_count: 0,
                    total_distance: 0.0,
                });
                link.support_count += cnt;
            }
        }

        let (mean_insert, stdev) = if !internal_distances.is_empty() {
            let sum: f64 = internal_distances.iter().sum();
            let mean = sum / internal_distances.len() as f64;
            let variance: f64 = internal_distances
                .iter()
                .map(|&d| (d - mean).powi(2))
                .sum::<f64>()
                / internal_distances.len() as f64;
            (mean, variance.sqrt())
        } else {
            (300.0, 50.0)
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
