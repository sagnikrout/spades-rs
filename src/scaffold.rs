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
    /// Connects contigs into scaffolds separated by 'N' runs.
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
            while let Some((next_idx, support)) = self.find_best_forward_link(curr, n, &visited, paired_info) {
                if support < self.min_links {
                    break;
                }

                // Insert gap
                current_seq.extend(std::iter::repeat(b'N').take(self.default_gap_len));
                current_seq.extend_from_slice(&contigs[next_idx].sequence);
                current_cov += contigs[next_idx].mean_coverage * contigs[next_idx].sequence.len() as f64;
                current_len += contigs[next_idx].sequence.len() + self.default_gap_len;

                visited[next_idx] = true;
                curr = next_idx;
            }

            scaffolds.push(Unitig {
                id: scaffold_id,
                sequence: current_seq,
                mean_coverage: current_cov / current_len as f64,
                kmers_count: current_len,
            });
            scaffold_id += 1;
        }

        scaffolds.sort_by(|a, b| b.sequence.len().cmp(&a.sequence.len()));
        scaffolds
    }

    fn find_best_forward_link(
        &self,
        curr: usize,
        n: usize,
        visited: &[bool],
        paired_info: &PairedInfoIndex,
    ) -> Option<(usize, u32)> {
        let mut best = None;
        let mut max_support = 0;

        for target in 0..n {
            if target == curr || visited[target] {
                continue;
            }
            let support = paired_info.get_support(curr, target);
            if support > max_support {
                max_support = support;
                best = Some((target, support));
            }
        }
        best
    }
}

/// Writes scaffolds to a FASTA file.
pub fn write_scaffolds_fasta<P: AsRef<Path>>(scaffolds: &[Unitig], out_path: P) -> Result<()> {
    let mut file = File::create(out_path)?;
    for (i, scaf) in scaffolds.iter().enumerate() {
        writeln!(
            file,
            ">scaffold_{}_len_{}_cov_{:.1}",
            i + 1,
            scaf.sequence.len(),
            scaf.mean_coverage
        )?;
        for chunk in scaf.sequence.chunks(80) {
            file.write_all(chunk)?;
            file.write_all(b"\n")?;
        }
    }
    Ok(())
}
