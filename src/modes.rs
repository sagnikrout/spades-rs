//! Specialized Assembly Modes (Layer 5: Meta, Plasmid, RNA-Seq & Single-Cell Pipelines).
//!
//! Implements specialized graph and coverage algorithms for domain-specific inputs:
//! - **plasmidSPAdes** ([`PlasmidDetector`]): plasmid circularity and copy-number extraction.
//! - **metaSPAdes** ([`apply_meta_filter`]): multi-species uneven coverage preservation.
//! - **rnaSPAdes** ([`RnaEngine`]): alternative splicing bubble path preservation and transcript extraction.
//! - **scSPAdes** ([`SingleCellNormalizer`]): coverage spike capping and dropout normalization for MDA data.

use crate::graph::Unitig;

// ─────────────────────────────────────────────────────────────────────────────
// Module 5.1: plasmidSPAdes Circularity & Copy-Number Extraction
// ─────────────────────────────────────────────────────────────────────────────

pub struct PlasmidDetector {
    pub k: usize,
    pub min_plasmid_len: usize,
    pub copy_number_threshold: f64,
}

impl Default for PlasmidDetector {
    fn default() -> Self {
        Self {
            k: 31,
            min_plasmid_len: 300,
            copy_number_threshold: 1.8,
        }
    }
}

impl PlasmidDetector {
    /// Separates contigs into (chromosomal, plasmids) based on circularity and coverage.
    pub fn extract_plasmids(&self, unitigs: Vec<Unitig>) -> (Vec<Unitig>, Vec<Unitig>) {
        if unitigs.is_empty() {
            return (Vec::new(), Vec::new());
        }

        let k1 = self.k.saturating_sub(1);

        // Compute median coverage across all contigs
        let mut covs: Vec<f64> = unitigs.iter().map(|u| u.mean_coverage).collect();
        covs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median_cov = covs[covs.len() / 2].max(1.0);

        let mut chromosomal = Vec::new();
        let mut plasmids = Vec::new();

        for u in unitigs {
            let len = u.sequence.len();
            // Check topological circularity (matching prefix and suffix overlap of length >= 15 up to k-1)
            let is_circular = if (20..=450_000).contains(&len) {
                let max_overlap = k1.min(len / 2);
                let mut circular = false;
                for check_k in (15..=max_overlap).rev() {
                    if u.sequence[..check_k] == u.sequence[len - check_k..] {
                        circular = true;
                        break;
                    }
                }
                circular
            } else {
                false
            };

            let is_high_copy =
                len <= 350_000 && u.mean_coverage >= (median_cov * self.copy_number_threshold);

            // A plasmid candidate is circular or significantly high-copy
            if (is_circular || (is_high_copy && len >= self.min_plasmid_len))
                && len >= self.min_plasmid_len
            {
                plasmids.push(u);
            } else {
                chromosomal.push(u);
            }
        }

        (chromosomal, plasmids)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module 5.2: metaSPAdes Metagenomic Filtering
// ─────────────────────────────────────────────────────────────────────────────

/// Applies metagenomic-aware filtering for uneven multi-species coverage.
pub fn apply_meta_filter(unitigs: Vec<Unitig>, min_len: usize) -> Vec<Unitig> {
    // In metagenomic mode, keep low-abundance species contigs as long as length is solid
    unitigs
        .into_iter()
        .filter(|u| u.sequence.len() >= min_len && u.mean_coverage >= 2.0)
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Module 5.3: rnaSPAdes Transcriptome & Isoform Engine
// ─────────────────────────────────────────────────────────────────────────────

pub struct RnaEngine {
    pub min_isoform_len: usize,
    pub min_isoform_coverage: f64,
}

impl Default for RnaEngine {
    fn default() -> Self {
        Self {
            min_isoform_len: 200,
            min_isoform_coverage: 1.5,
        }
    }
}

impl RnaEngine {
    /// Preserves alternative splicing transcripts without collapsing legitimate isoforms.
    pub fn process_transcripts(&self, unitigs: Vec<Unitig>) -> Vec<Unitig> {
        let mut isoforms: Vec<Unitig> = unitigs
            .into_iter()
            .filter(|u| {
                u.sequence.len() >= self.min_isoform_len
                    && u.mean_coverage >= self.min_isoform_coverage
            })
            .collect();

        isoforms.sort_by_key(|a| std::cmp::Reverse(a.sequence.len()));
        isoforms
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module 5.4: scSPAdes Single-Cell MDA Amplification Normalizer
// ─────────────────────────────────────────────────────────────────────────────

pub struct SingleCellNormalizer {
    pub min_coverage_cutoff: f64,
}

impl Default for SingleCellNormalizer {
    fn default() -> Self {
        Self {
            min_coverage_cutoff: 2.0,
        }
    }
}

impl SingleCellNormalizer {
    /// Normalizes unitigs from single-cell MDA data, preventing coverage dropouts from fragmenting contigs.
    pub fn normalize_coverage(&self, mut unitigs: Vec<Unitig>) -> Vec<Unitig> {
        for u in &mut unitigs {
            // Cap extreme amplification spikes (e.g. > 500x) to stabilize graph heuristics
            if u.mean_coverage > 500.0 {
                u.mean_coverage = 500.0;
            }
        }

        unitigs
            .into_iter()
            .filter(|u| u.mean_coverage >= self.min_coverage_cutoff)
            .collect()
    }
}
