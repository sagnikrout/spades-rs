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
    fn build_kmer_index(k: usize, unitigs: &[Unitig]) -> HashMap<Kmer256, (usize, usize)> {
        let total_kmers: usize = unitigs
            .iter()
            .map(|u| u.sequence.len().saturating_sub(k - 1))
            .sum();
        let mut kmer_to_unitig: HashMap<Kmer256, (usize, usize)> =
            HashMap::with_capacity(total_kmers);
        let mut seen_repeats: hashbrown::HashSet<Kmer256> = hashbrown::HashSet::new();

        for (u_id, u) in unitigs.iter().enumerate() {
            if u.sequence.len() < k {
                continue;
            }
            for i in 0..=(u.sequence.len() - k) {
                if let Some(km) = Kmer256::from_bytes(&u.sequence[i..i + k], k) {
                    let (can, _) = km.canonical(k);
                    if seen_repeats.contains(&can) {
                        continue;
                    }
                    if kmer_to_unitig.contains_key(&can) {
                        kmer_to_unitig.remove(&can);
                        seen_repeats.insert(can);
                    } else {
                        kmer_to_unitig.insert(can, (u_id, i));
                    }
                }
            }
        }
        kmer_to_unitig
    }

    #[inline]
    fn map_read(
        seq: &[u8],
        k: usize,
        kmer_to_unitig: &HashMap<Kmer256, (usize, usize)>,
    ) -> Option<(usize, usize)> {
        if seq.len() < k {
            return None;
        }
        let last_start = seq.len() - k;
        let step = (last_start / 8).clamp(1, 8);
        for i in (0..=last_start).step_by(step) {
            if let Some(km) = Kmer256::from_bytes(&seq[i..i + k], k) {
                let (can, _) = km.canonical(k);
                if let Some(&loc) = kmer_to_unitig.get(&can) {
                    return Some(loc);
                }
            }
        }
        if !last_start.is_multiple_of(step) {
            if let Some(km) = Kmer256::from_bytes(&seq[last_start..last_start + k], k) {
                let (can, _) = km.canonical(k);
                if let Some(&loc) = kmer_to_unitig.get(&can) {
                    return Some(loc);
                }
            }
        }
        None
    }

    #[inline]
    fn record_pair_mapping(
        loc1: Option<(usize, usize)>,
        loc2: Option<(usize, usize)>,
        read2_len: usize,
        local_dists: &mut Vec<f64>,
        local_links: &mut HashMap<(usize, usize), u32>,
    ) {
        if let (Some((u1, pos1)), Some((u2, pos2))) = (loc1, loc2) {
            if u1 == u2 {
                let offset_diff = (pos1 as f64 - pos2 as f64).abs();
                if offset_diff > 0.0 {
                    let insert_span = offset_diff + read2_len as f64;
                    if insert_span <= 2_000.0 {
                        local_dists.push(insert_span);
                    }
                }
            } else {
                *local_links.entry((u1, u2)).or_insert(0) += 1;
                *local_links.entry((u2, u1)).or_insert(0) += 1;
            }
        }
    }

    #[allow(clippy::type_complexity)]
    fn finalize_from_threads(
        thread_results: Vec<(Vec<f64>, HashMap<(usize, usize), u32>)>,
    ) -> Self {
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
            (300.0, 50.0) // Default typical Illumina insert size
        };

        Self {
            mean_insert_size: mean_insert,
            insert_size_stdev: stdev,
            links,
        }
    }

    /// Builds paired-end connections from paired reads and unitigs.
    pub fn build(k: usize, unitigs: &[Unitig], reads1: &[Vec<u8>], reads2: &[Vec<u8>]) -> Self {
        let kmer_to_unitig = Self::build_kmer_index(k, unitigs);
        let num_pairs = reads1.len().min(reads2.len());
        let pair_indices: Vec<usize> = (0..num_pairs).collect();

        #[allow(clippy::type_complexity)]
        let thread_results: Vec<(Vec<f64>, HashMap<(usize, usize), u32>)> = pair_indices
            .par_chunks(4096)
            .map(|chunk| {
                let mut local_dists = Vec::new();
                let mut local_links: HashMap<(usize, usize), u32> = HashMap::new();
                for &i in chunk {
                    let loc1 = Self::map_read(&reads1[i], k, &kmer_to_unitig);
                    let loc2 = Self::map_read(&reads2[i], k, &kmer_to_unitig);
                    Self::record_pair_mapping(
                        loc1,
                        loc2,
                        reads2[i].len(),
                        &mut local_dists,
                        &mut local_links,
                    );
                }
                (local_dists, local_links)
            })
            .collect();

        Self::finalize_from_threads(thread_results)
    }

    /// Builds paired-end connections directly from PackedReads without unpacking vectors.
    pub fn build_from_packed(
        k: usize,
        unitigs: &[Unitig],
        packed: &crate::packed_reads::PackedReads,
    ) -> Self {
        let kmer_to_unitig = Self::build_kmer_index(k, unitigs);
        let num_pairs = packed.pe_count();
        let pair_indices: Vec<usize> = (0..num_pairs).collect();

        #[allow(clippy::type_complexity)]
        let thread_results: Vec<(Vec<f64>, HashMap<(usize, usize), u32>)> = pair_indices
            .par_chunks(4096)
            .map(|chunk| {
                let mut buf1 = Vec::with_capacity(512);
                let mut buf2 = Vec::with_capacity(512);
                let mut local_dists = Vec::new();
                let mut local_links: HashMap<(usize, usize), u32> = HashMap::new();

                for &p in chunk {
                    packed.get_pair(p, &mut buf1, &mut buf2);
                    let loc1 = Self::map_read(&buf1, k, &kmer_to_unitig);
                    let loc2 = Self::map_read(&buf2, k, &kmer_to_unitig);
                    Self::record_pair_mapping(
                        loc1,
                        loc2,
                        buf2.len(),
                        &mut local_dists,
                        &mut local_links,
                    );
                }
                (local_dists, local_links)
            })
            .collect();

        Self::finalize_from_threads(thread_results)
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
