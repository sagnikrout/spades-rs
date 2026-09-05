//! scSPAdes Single-Cell MDA Amplification Normalizer (Module 5.4).
//!
//! Handles extreme 10x-1000x coverage fluctuations typical of Multiple Displacement Amplification (MDA).

use crate::graph::Unitig;

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
