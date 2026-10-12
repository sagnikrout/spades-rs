//! Mismatch Corrector & Base Polisher (Module 4.2).
//!
//! Re-aligns reads against assembled contigs to calculate position-wise consensus,
//! correcting residual substitution errors and homopolymer indels.

use crate::dna::Kmer256;
use crate::graph::Unitig;
use hashbrown::HashMap;
use rayon::prelude::*;
use std::sync::atomic::{AtomicU32, Ordering};

pub struct Polisher {
    pub k: usize,
    pub min_coverage_support: u32,
    pub careful: bool,
}

impl Default for Polisher {
    fn default() -> Self {
        Self {
            k: 21,
            min_coverage_support: 5,
            careful: false,
        }
    }
}

impl Polisher {
    fn build_kmer_index(k: usize, contigs: &[Unitig]) -> HashMap<Kmer256, (u32, u32, bool)> {
        let total_kmers: usize = contigs
            .iter()
            .map(|c| c.sequence.len().saturating_sub(k - 1))
            .sum();
        let mut kmer_pos: HashMap<Kmer256, (u32, u32, bool)> = HashMap::with_capacity(total_kmers);
        let mut seen_repeats: hashbrown::HashSet<Kmer256> = hashbrown::HashSet::new();

        for (c_idx, contig) in contigs.iter().enumerate() {
            let len = contig.sequence.len();
            if len < k {
                continue;
            }
            for i in 0..=(len - k) {
                if let Some(km) = Kmer256::from_bytes(&contig.sequence[i..i + k], k) {
                    let (can, is_rc) = km.canonical(k);
                    if seen_repeats.contains(&can) {
                        continue;
                    }
                    if let Some(&(existing_c, _, _)) = kmer_pos.get(&can) {
                        // Multi-mapped k-mer across different contigs
                        if k >= 11 || (existing_c as usize) != c_idx {
                            kmer_pos.remove(&can);
                            seen_repeats.insert(can);
                        }
                    } else {
                        kmer_pos.insert(can, (c_idx as u32, i as u32, is_rc));
                    }
                }
            }
        }
        kmer_pos
    }

    fn init_tallies(contigs: &[Unitig]) -> Vec<Vec<[AtomicU32; 4]>> {
        contigs
            .iter()
            .map(|c| {
                (0..c.sequence.len())
                    .map(|_| {
                        [
                            AtomicU32::new(0),
                            AtomicU32::new(0),
                            AtomicU32::new(0),
                            AtomicU32::new(0),
                        ]
                    })
                    .collect()
            })
            .collect()
    }

    #[inline]
    fn align_read_to_tallies(
        read: &[u8],
        k: usize,
        kmer_pos: &HashMap<Kmer256, (u32, u32, bool)>,
        tallies: &[Vec<[AtomicU32; 4]>],
    ) {
        if read.len() < k {
            return;
        }
        let step = (k / 2 + 1).max(1);
        for r_idx in (0..=(read.len() - k)).step_by(step) {
            if let Some(km) = Kmer256::from_bytes(&read[r_idx..r_idx + k], k) {
                let (can, r_is_rc) = km.canonical(k);
                if let Some(&(c_idx, c_pos, c_is_rc)) = kmer_pos.get(&can) {
                    let c_idx = c_idx as usize;
                    let c_pos = c_pos as usize;
                    let c_len = tallies[c_idx].len();

                    if r_is_rc == c_is_rc {
                        // Forward alignment: read[r_idx] aligns to contig[c_pos],
                        // so read[read_offset] aligns to contig[c_pos - r_idx + read_offset].
                        let start_contig = c_pos as isize - r_idx as isize;
                        for (read_offset, &base) in read.iter().enumerate() {
                            let curr_c_pos = start_contig + read_offset as isize;
                            if curr_c_pos >= 0 && (curr_c_pos as usize) < c_len {
                                if let Some(code) = crate::dna::base_to_2bit(base) {
                                    tallies[c_idx][curr_c_pos as usize][code as usize]
                                        .fetch_add(1, Ordering::Relaxed);
                                }
                            }
                        }
                    } else {
                        // Reverse-complement alignment: read[r_idx..r_idx+k] is RC of contig[c_pos..c_pos+k].
                        // Therefore read[r_idx] aligns to contig[c_pos + k - 1], and read[r_idx + k - 1] aligns to contig[c_pos].
                        // In general, read[read_offset] aligns to contig[(c_pos + k - 1 + r_idx) - read_offset].
                        let start_contig = c_pos as isize + k as isize - 1 + r_idx as isize;
                        for (read_offset, &base) in read.iter().enumerate() {
                            let curr_c_pos = start_contig - read_offset as isize;
                            if curr_c_pos >= 0 && (curr_c_pos as usize) < c_len {
                                if let Some(code) = crate::dna::base_to_2bit(base) {
                                    let rc_code = (!code) & 3;
                                    tallies[c_idx][curr_c_pos as usize][rc_code as usize]
                                        .fetch_add(1, Ordering::Relaxed);
                                }
                            }
                        }
                    }
                    break;
                }
            }
        }
    }

    fn apply_consensus(
        &self,
        mut contigs: Vec<Unitig>,
        tallies: &[Vec<[AtomicU32; 4]>],
    ) -> (Vec<Unitig>, usize) {
        let (min_support, dominance_ratio) = if self.careful {
            (self.min_coverage_support.max(8), 0.85)
        } else {
            (self.min_coverage_support, 0.80)
        };

        let mut total_corrections = 0;
        for (c_idx, contig) in contigs.iter_mut().enumerate() {
            let len = contig.sequence.len();
            #[allow(clippy::needless_range_loop)]
            for pos in 0..len {
                let counts = &tallies[c_idx][pos];
                let c0 = counts[0].load(Ordering::Relaxed);
                let c1 = counts[1].load(Ordering::Relaxed);
                let c2 = counts[2].load(Ordering::Relaxed);
                let c3 = counts[3].load(Ordering::Relaxed);
                let total = c0 + c1 + c2 + c3;

                if total >= min_support {
                    let ary = [c0, c1, c2, c3];
                    let mut max_code = 0;
                    let mut max_c = 0;
                    for (code, &c) in ary.iter().enumerate() {
                        if c > max_c {
                            max_c = c;
                            max_code = code;
                        }
                    }

                    if (max_c as f64 / total as f64) >= dominance_ratio {
                        let consensus_base = crate::dna::bit2_to_base(max_code as u8);
                        if contig.sequence[pos] != consensus_base {
                            contig.sequence[pos] = consensus_base;
                            total_corrections += 1;
                        }
                    }
                }
            }
        }

        (contigs, total_corrections)
    }

    /// Polishes contigs by consensus base-calling.
    pub fn polish_contigs(&self, contigs: Vec<Unitig>, reads: &[Vec<u8>]) -> (Vec<Unitig>, usize) {
        let k = self.k;
        if contigs.is_empty() || reads.is_empty() {
            return (contigs, 0);
        }

        let kmer_pos = Self::build_kmer_index(k, &contigs);
        let tallies = Self::init_tallies(&contigs);

        reads.par_chunks(4096).for_each(|chunk| {
            for read in chunk {
                Self::align_read_to_tallies(read, k, &kmer_pos, &tallies);
            }
        });

        self.apply_consensus(contigs, &tallies)
    }

    /// Polishes contigs by consensus base-calling directly from PackedReads without intermediate vectors.
    pub fn polish_contigs_packed(
        &self,
        contigs: Vec<Unitig>,
        packed: &crate::packed_reads::PackedReads,
    ) -> (Vec<Unitig>, usize) {
        let k = self.k;
        if contigs.is_empty() || (packed.is_empty() && packed.ambiguous_sidecar.is_empty()) {
            return (contigs, 0);
        }

        let kmer_pos = Self::build_kmer_index(k, &contigs);
        let tallies = Self::init_tallies(&contigs);

        let num_reads = packed.len();
        let chunk_size = 4096;
        let num_chunks = num_reads.div_ceil(chunk_size);

        (0..num_chunks).into_par_iter().for_each(|chunk_idx| {
            let start = chunk_idx * chunk_size;
            let end = (start + chunk_size).min(num_reads);
            let mut read = Vec::with_capacity(512);
            for idx in start..end {
                packed.get_read(idx, &mut read);
                Self::align_read_to_tallies(&read, k, &kmer_pos, &tallies);
            }
        });

        // Sidecar pass: ambiguous reads with N positions vote with their clean flanking bases.
        // N bases return None from base_to_2bit and contribute zero to the tally.
        for raw_read in &packed.ambiguous_sidecar {
            Self::align_read_to_tallies(raw_read, k, &kmer_pos, &tallies);
        }

        self.apply_consensus(contigs, &tallies)
    }
}
