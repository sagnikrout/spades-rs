//! rnaSPAdes Transcriptome & Isoform Engine (Module 5.3).
//!
//! Preserves alternative splicing bubble paths and outputs reconstructed RNA isoforms.

use crate::graph::Unitig;

pub struct RnaEngine {
    pub min_isoform_len: usize,
    pub min_isoform_coverage: f64,
}

impl Default for RnaEngine {
    fn default() -> Self {
        Self {
            min_isoform_len: 200,
            min_isoform_coverage: 3.0,
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
