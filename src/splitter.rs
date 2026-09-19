//! SpLitteR: Linked-Read / Synthetic Long-Read Scaffolding Engine (Module 6.6).
//!
//! Reconstructs long-range genomic scaffolds from barcoded linked reads (e.g., 10x Genomics, TELL-Seq)
//! by tracing shared barcode co-occurrence across assembly graph unitigs to resolve repeat bifurcations.
//! References: Tolstoganov et al. (2024) PeerJ 10.7717/peerj.18050.

use crate::dna::Kmer256;
use crate::fastq::LinkedReadRecord;
use crate::graph::Unitig;
use hashbrown::{HashMap, HashSet};
use rayon::prelude::*;

/// Configuration and engine for linked-read repeat resolution.
#[derive(Clone, Debug)]
pub struct LinkedReadResolver {
    /// K-mer size used for graph anchoring.
    pub k: usize,
    /// Minimum number of distinct barcoded molecules supporting a bridge between unitigs.
    pub min_shared_barcodes: usize,
    /// Minimum read count per barcode cluster to be considered informative.
    pub min_reads_per_barcode: usize,
}

impl Default for LinkedReadResolver {
    fn default() -> Self {
        Self {
            k: 31,
            min_shared_barcodes: 3,
            min_reads_per_barcode: 2,
        }
    }
}

impl LinkedReadResolver {
    pub fn new(k: usize) -> Self {
        Self {
            k,
            ..Default::default()
        }
    }

    /// Resolves repeat bifurcations and scaffolds unitigs using linked-read barcode co-occurrence.
    pub fn bridge_with_linked_reads(
        &self,
        contigs: Vec<Unitig>,
        linked_reads: &[LinkedReadRecord],
    ) -> Vec<Unitig> {
        if contigs.len() <= 1 || linked_reads.is_empty() {
            return contigs;
        }

        let k = self.k;

        // Step 1: Index unique k-mer anchors across contigs
        // Build map: canonical k-mer -> (contig_idx, is_rc, pos)
        let mut kmer_to_contig: HashMap<Kmer256, (usize, bool, usize)> = HashMap::new();
        let mut ambiguous_kmers: HashSet<Kmer256> = HashSet::new();

        for (c_idx, c) in contigs.iter().enumerate() {
            if c.sequence.len() < k {
                continue;
            }
            for i in 0..=(c.sequence.len() - k) {
                if let Some(kmer) = Kmer256::from_bytes(&c.sequence[i..i + k], k) {
                    let (can, is_rc) = kmer.canonical(k);
                    if ambiguous_kmers.contains(&can) {
                        continue;
                    }
                    if let Some(&(existing_c, _, _)) = kmer_to_contig.get(&can) {
                        if existing_c != c_idx {
                            kmer_to_contig.remove(&can);
                            ambiguous_kmers.insert(can);
                        }
                    } else {
                        kmer_to_contig.insert(can, (c_idx, is_rc, i));
                    }
                }
            }
        }

        if kmer_to_contig.is_empty() {
            return contigs;
        }

        // Step 2: Group reads by barcode and collect visited contigs
        // barcode -> (set of visited contigs)
        let mut barcode_groups: HashMap<u64, Vec<&[u8]>> = HashMap::new();
        for r in linked_reads {
            if let Some(bc) = r.barcode {
                barcode_groups.entry(bc).or_default().push(&r.seq);
            }
        }

        // Filter barcodes by minimum read count and map to contigs
        let barcode_contigs: Vec<HashSet<usize>> = barcode_groups
            .into_par_iter()
            .filter(|(_, seqs)| seqs.len() >= self.min_reads_per_barcode)
            .map(|(_, seqs)| {
                let mut visited = HashSet::new();
                for seq in seqs {
                    if seq.len() < k {
                        continue;
                    }
                    // Sample k-mers along read
                    let step = (k / 2).max(1);
                    for i in (0..=(seq.len() - k)).step_by(step) {
                        if let Some(kmer) = Kmer256::from_bytes(&seq[i..i + k], k) {
                            let (can, _) = kmer.canonical(k);
                            if let Some(&(c_idx, _, _)) = kmer_to_contig.get(&can) {
                                visited.insert(c_idx);
                            }
                        }
                    }
                }
                visited
            })
            .collect();

        // Step 3: Count pairwise barcode co-occurrences between contigs
        let mut pair_counts: HashMap<(usize, usize), usize> = HashMap::new();
        for visited in barcode_contigs {
            if visited.len() < 2 || visited.len() > 20 {
                // Ignore empty or overly promiscuous barcode droplets
                continue;
            }
            let list: Vec<usize> = visited.into_iter().collect();
            for i in 0..list.len() {
                for j in (i + 1)..list.len() {
                    let u = list[i].min(list[j]);
                    let v = list[i].max(list[j]);
                    *pair_counts.entry((u, v)).or_insert(0) += 1;
                }
            }
        }

        // Step 4: Filter connections with sufficient shared barcodes
        let mut candidates: Vec<((usize, usize), usize)> = pair_counts
            .into_iter()
            .filter(|&(_, count)| count >= self.min_shared_barcodes)
            .collect();

        // Sort descending by barcode support
        candidates.sort_unstable_by_key(|a| std::cmp::Reverse(a.1));

        if candidates.is_empty() {
            return contigs;
        }

        // Step 5: Greedy non-conflicting scaffold chaining
        let n = contigs.len();
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(i: usize, parent: &mut [usize]) -> usize {
            if parent[i] == i {
                i
            } else {
                let root = find(parent[i], parent);
                parent[i] = root;
                root
            }
        }

        let mut right_neighbor: Vec<Option<usize>> = vec![None; n];
        let mut has_left_neighbor = vec![false; n];

        let mut connections = 0;
        for &((u, v), _) in &candidates {
            let root_u = find(u, &mut parent);
            let root_v = find(v, &mut parent);
            if root_u == root_v {
                continue; // Avoid cycles
            }

            if right_neighbor[u].is_none() && !has_left_neighbor[v] {
                right_neighbor[u] = Some(v);
                has_left_neighbor[v] = true;
                parent[root_u] = root_v;
                connections += 1;
            } else if right_neighbor[v].is_none() && !has_left_neighbor[u] {
                right_neighbor[v] = Some(u);
                has_left_neighbor[u] = true;
                parent[root_v] = root_u;
                connections += 1;
            }
        }

        if connections == 0 {
            return contigs;
        }

        // Reconstruct scaffolds
        let mut chained = Vec::with_capacity(n);
        let mut visited = vec![false; n];

        for i in 0..n {
            if has_left_neighbor[i] || visited[i] {
                continue;
            }

            // Walk along chain
            let mut curr = i;
            let mut scaffold_seq = contigs[curr].sequence.clone();
            let mut total_coverage =
                contigs[curr].mean_coverage * (contigs[curr].sequence.len() as f64);
            let mut total_len = contigs[curr].sequence.len();
            let mut total_kmers = contigs[curr].kmers_count;
            visited[curr] = true;

            while let Some(next) = right_neighbor[curr] {
                visited[next] = true;
                // Insert standard 100-bp N scaffold gap
                scaffold_seq.extend_from_slice(&[b'N'; 100]);
                scaffold_seq.extend_from_slice(&contigs[next].sequence);
                total_coverage +=
                    contigs[next].mean_coverage * (contigs[next].sequence.len() as f64);
                total_len += contigs[next].sequence.len();
                total_kmers += contigs[next].kmers_count;
                curr = next;
            }

            let avg_cov = if total_len > 0 {
                total_coverage / (total_len as f64)
            } else {
                contigs[i].mean_coverage
            };

            chained.push(Unitig {
                id: contigs[i].id,
                sequence: scaffold_seq,
                mean_coverage: avg_cov,
                kmers_count: total_kmers,
            });
        }

        // Add any remaining unchained contigs
        for (i, c) in contigs.into_iter().enumerate() {
            if !visited[i] {
                chained.push(c);
            }
        }

        chained
    }
}
