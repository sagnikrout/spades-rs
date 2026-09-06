//! Late Gap Closer & Scaffolder (Module 4.1).
//!
//! Orders and orients contigs into scaffolds across unresolved gaps using
//! paired-end insert distance estimates, inserting estimated 'N' runs.

use crate::graph::Unitig;
use crate::paired_info::PairedInfoIndex;
use anyhow::Result;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::dna::{revcomp_bytes, string_to_kmer};
use hashbrown::HashMap;

/// Bounded local de Bruijn graph gap closer for bridging contig ends.
pub struct LocalGapCloser {
    pub k: usize,
    pub max_gap_len: usize,
}

impl Default for LocalGapCloser {
    fn default() -> Self {
        Self {
            k: 21,
            max_gap_len: 2000,
        }
    }
}

impl LocalGapCloser {
    /// Attempts to close the gap between left_seq and right_seq.
    /// Returns Some(joined_sequence) if closed, or None if unclosed.
    pub fn close_gap(
        &self,
        left_seq: &[u8],
        right_seq: &[u8],
        candidate_reads: &[&[u8]],
    ) -> Option<Vec<u8>> {
        // 1. Check direct prefix/suffix overlap
        let min_ov = self.k.min(4);
        let max_ov = 31.min(left_seq.len()).min(right_seq.len());
        if max_ov >= min_ov {
            for ov in (min_ov..=max_ov).rev() {
                if left_seq[left_seq.len() - ov..] == right_seq[..ov] {
                    let mut merged = left_seq.to_vec();
                    merged.extend_from_slice(&right_seq[ov..]);
                    return Some(merged);
                }
            }
        }

        if candidate_reads.is_empty() {
            return None;
        }

        let k = self.k;
        let start_kmer = string_to_kmer(&left_seq[left_seq.len() - k..], k)?;
        let end_kmer = string_to_kmer(&right_seq[..k], k)?;
        let rc_start = crate::dna::revcomp_kmer_u64(start_kmer, k);
        let rc_end = crate::dna::revcomp_kmer_u64(end_kmer, k);

        // 2. Build local directed k-mer count table from oriented candidate reads
        let mut kmer_counts: HashMap<u64, u32> = HashMap::with_capacity(candidate_reads.len() * 50);
        for read in candidate_reads {
            if read.len() < k {
                continue;
            }
            let mut is_rc = false;
            let mut has_fwd = false;
            for i in 0..=(read.len() - k) {
                if let Some(km) = string_to_kmer(&read[i..i + k], k) {
                    if km == rc_start || km == rc_end {
                        is_rc = true;
                    }
                    if km == start_kmer || km == end_kmer {
                        has_fwd = true;
                    }
                }
            }

            let oriented_seq = if is_rc && !has_fwd {
                revcomp_bytes(read)
            } else {
                read.to_vec()
            };

            for i in 0..=(oriented_seq.len() - k) {
                if let Some(km) = string_to_kmer(&oriented_seq[i..i + k], k) {
                    *kmer_counts.entry(km).or_insert(0) += 1;
                }
            }
        }

        let mut cur = start_kmer;
        let mut walked = Vec::new();
        let mut visited = hashbrown::HashSet::new();
        visited.insert(cur);

        let mask = if k >= 32 {
            !0u64
        } else {
            (1u64 << (2 * (k - 1))) - 1
        };

        for _ in 0..self.max_gap_len {
            if cur == end_kmer {
                if walked.len() >= k {
                    let insert_len = walked.len() - k;
                    let mut merged = left_seq.to_vec();
                    merged.extend_from_slice(&walked[..insert_len]);
                    merged.extend_from_slice(right_seq);
                    return Some(merged);
                } else {
                    let mut merged = left_seq.to_vec();
                    merged.extend_from_slice(right_seq);
                    return Some(merged);
                }
            }

            let prefix = (cur & mask) << 2;
            let mut candidates: Vec<(u64, u32, u8)> = Vec::with_capacity(4);

            for &(code, byte) in [(0u64, b'A'), (1u64, b'C'), (2u64, b'G'), (3u64, b'T')].iter() {
                let nxt = prefix | code;
                if let Some(&cnt) = kmer_counts.get(&nxt) {
                    if cnt >= 2 {
                        candidates.push((nxt, cnt, byte));
                    }
                }
            }

            if candidates.is_empty() {
                break;
            }

            candidates.sort_by_key(|c| std::cmp::Reverse(c.1));
            let (best_nxt, best_cnt, best_byte) = candidates[0];

            if candidates.len() > 1 && candidates[1].1 * 2 > best_cnt {
                break; // Ambiguous branch
            }

            if visited.contains(&best_nxt) {
                break; // Cycle
            }

            visited.insert(best_nxt);
            walked.push(best_byte);
            cur = best_nxt;
        }

        None
    }
}

pub struct Scaffolder {
    pub min_links: u32,
    pub default_gap_len: usize,
}

impl Default for Scaffolder {
    fn default() -> Self {
        Self {
            min_links: 3,
            default_gap_len: 50,
        }
    }
}

impl Scaffolder {
    /// Connects contigs into scaffolds across unresolved gaps, attempting local gap closure.
    pub fn build_scaffolds(
        &self,
        contigs: &[Unitig],
        paired_info: &PairedInfoIndex,
    ) -> Vec<Unitig> {
        if contigs.len() <= 1 {
            return contigs.to_vec();
        }

        let n = contigs.len();
        let mut visited = vec![false; n];
        let mut scaffolds = Vec::new();
        let mut scaffold_id = 0;
        let gap_closer = LocalGapCloser::default();

        for i in 0..n {
            if visited[i] {
                continue;
            }

            let mut current_seq = contigs[i].sequence.clone();
            let mut current_cov = contigs[i].mean_coverage * contigs[i].sequence.len() as f64;
            let mut current_len = contigs[i].sequence.len();
            visited[i] = true;
            let mut curr = i;

            // Greedily find best forward link
            while let Some((next_idx, support)) =
                self.find_best_forward_link(curr, n, &visited, paired_info)
            {
                if support < self.min_links {
                    break;
                }

                let next_seq = &contigs[next_idx].sequence;
                if let Some(merged) = gap_closer.close_gap(&current_seq, next_seq, &[]) {
                    let added_len = merged.len().saturating_sub(current_seq.len());
                    current_cov +=
                        contigs[next_idx].mean_coverage * contigs[next_idx].sequence.len() as f64;
                    current_len += added_len;
                    current_seq = merged;
                } else {
                    // Insert gap
                    current_seq.extend(std::iter::repeat_n(b'N', self.default_gap_len));
                    current_seq.extend_from_slice(next_seq);
                    current_cov +=
                        contigs[next_idx].mean_coverage * contigs[next_idx].sequence.len() as f64;
                    current_len += contigs[next_idx].sequence.len() + self.default_gap_len;
                }

                visited[next_idx] = true;
                curr = next_idx;
            }

            scaffolds.push(Unitig {
                id: scaffold_id,
                sequence: current_seq,
                mean_coverage: if current_len > 0 {
                    current_cov / current_len as f64
                } else {
                    10.0
                },
                kmers_count: current_len,
            });
            scaffold_id += 1;
        }

        scaffolds.sort_by_key(|a| std::cmp::Reverse(a.sequence.len()));
        scaffolds
    }

    fn find_best_forward_link(
        &self,
        curr: usize,
        _n: usize,
        visited: &[bool],
        paired_info: &PairedInfoIndex,
    ) -> Option<(usize, u32)> {
        let mut best = None;
        let mut max_support = 0;

        if let Some(targets) = paired_info.links.get(&curr) {
            for (&target, link) in targets {
                if target < visited.len() && !visited[target] && link.support_count > max_support {
                    max_support = link.support_count;
                    best = Some((target, link.support_count));
                }
            }
        }
        best
    }
}

/// Writes scaffolds to a FASTA file.
pub fn write_scaffolds_fasta<P: AsRef<Path>>(scaffolds: &[Unitig], out_path: P) -> Result<()> {
    let file = File::create(out_path)?;
    let mut writer = std::io::BufWriter::with_capacity(1024 * 1024, file);
    for (i, scaf) in scaffolds.iter().enumerate() {
        writeln!(
            writer,
            ">scaffold_{}_len_{}_cov_{:.1}",
            i + 1,
            scaf.sequence.len(),
            scaf.mean_coverage
        )?;
        for chunk in scaf.sequence.chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }
    writer.flush()?;
    Ok(())
}
