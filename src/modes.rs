//! Specialized Assembly Modes (Layer 5: Meta & Plasmid Pipelines).
//!
//! Implements plasmid circularity and copy-number detection (plasmidSPAdes)
//! and multi-species uneven coverage preservation (metaSPAdes).

use crate::graph::Unitig;

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
            copy_number_threshold: 2.0,
        }
    }
}

impl PlasmidDetector {
    /// Separates contigs into (chromosomal, plasmids) based on circularity and coverage.
    pub fn extract_plasmids(&self, unitigs: Vec<Unitig>) -> (Vec<Unitig>, Vec<Unitig>) {
        if unitigs.is_empty() {
            return (Vec::new(), Vec::new());
        }

        let k1 = self.k - 1;

        // Compute median coverage across all contigs
        let mut covs: Vec<f64> = unitigs.iter().map(|u| u.mean_coverage).collect();
        covs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median_cov = covs[covs.len() / 2].max(1.0);

        let mut chromosomal = Vec::new();
        let mut plasmids = Vec::new();

        for u in unitigs {
            let len = u.sequence.len();
            let is_circular = if len >= k1 && k1 > 0 {
                let prefix = &u.sequence[..k1];
                let suffix = &u.sequence[len - k1..];
                prefix == suffix
            } else {
                false
            };

            let is_high_copy = u.mean_coverage >= (median_cov * self.copy_number_threshold);

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

/// Applies metagenomic-aware filtering for uneven multi-species coverage.
pub fn apply_meta_filter(unitigs: Vec<Unitig>, min_len: usize) -> Vec<Unitig> {
    // In metagenomic mode, keep low-abundance species contigs as long as length is solid
    unitigs
        .into_iter()
        .filter(|u| u.sequence.len() >= min_len && u.mean_coverage >= 2.0)
        .collect()
}
