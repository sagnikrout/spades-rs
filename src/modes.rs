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

/// Detector that isolates circular and high-copy-number plasmid contigs (`plasmidSPAdes` mode).
pub struct PlasmidDetector {
    /// $k$-mer size used when checking prefix/suffix circularity overlaps.
    pub k: usize,
    /// Minimum sequence length (in base pairs) required for a contig to be classified as a plasmid.
    pub min_plasmid_len: usize,
    /// Minimum coverage ratio relative to the median chromosomal coverage for high-copy plasmid detection.
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
        covs.sort_by(|a, b| a.total_cmp(b));
        let median_cov = covs[covs.len() / 2].max(1.0);

        let mut chromosomal = Vec::new();
        let mut plasmids = Vec::new();

        for u in unitigs {
            let len = u.sequence.len();
            // Check topological circularity (matching prefix and suffix overlap up to k-1)
            let is_circular = if k1 > 0 && (20..=450_000).contains(&len) {
                let min_overlap = 15.min(k1);
                let max_overlap = k1.min(len / 2);
                let mut circular = false;
                if max_overlap >= min_overlap {
                    for check_k in (min_overlap..=max_overlap).rev() {
                        if u.sequence[..check_k] == u.sequence[len - check_k..] {
                            circular = true;
                            break;
                        }
                    }
                }
                circular
            } else {
                false
            };

            let is_high_copy =
                len <= 350_000 && u.mean_coverage >= (median_cov * self.copy_number_threshold);

            // A plasmid candidate is circular or significantly high-copy
            if (is_circular || is_high_copy) && len >= self.min_plasmid_len {
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

/// Transcriptome isoform filter that preserves alternative splicing paths (`rnaSPAdes` mode).
pub struct RnaEngine {
    /// Minimum transcript length (in base pairs) required to retain an assembled isoform.
    pub min_isoform_len: usize,
    /// Minimum mean $k$-mer coverage required to retain a low-expression transcript isoform.
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

/// Coverage normalizer for single-cell Multiple Displacement Amplification (MDA) assemblies (`scSPAdes` mode).
pub struct SingleCellNormalizer {
    /// Minimum mean coverage threshold required to retain low-coverage dropout regions.
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

// ─────────────────────────────────────────────────────────────────────────────
// Module 5.5: dipSPAdes Diploid Haplotype Separation & Haplotig Extraction
// ─────────────────────────────────────────────────────────────────────────────

/// Diploid haplotype phaser and bubble separator (`dipSPAdes` mode).
#[derive(Clone, Debug)]
pub struct DiploidPhaser {
    /// Maximum length difference ratio between homologous bubble alleles.
    pub max_allele_len_diff: f64,
    /// Maximum coverage asymmetry ratio between homologous alleles (e.g., 2.5).
    pub max_cov_ratio: f64,
}

impl Default for DiploidPhaser {
    fn default() -> Self {
        Self {
            max_allele_len_diff: 0.30,
            max_cov_ratio: 2.5,
        }
    }
}

impl DiploidPhaser {
    /// Separates parallel polymorphic bubble branches into primary reference contigs and secondary haplotigs.
    pub fn phase_haplotigs(&self, k: usize, unitigs: Vec<Unitig>) -> (Vec<Unitig>, Vec<Unitig>) {
        let k1 = k.saturating_sub(1);
        if k1 == 0 || unitigs.len() <= 1 {
            return (unitigs, Vec::new());
        }

        let mut bubble_groups: hashbrown::HashMap<(Vec<u8>, Vec<u8>), Vec<usize>> =
            hashbrown::HashMap::new();

        for (idx, u) in unitigs.iter().enumerate() {
            if u.sequence.len() < k1 {
                continue;
            }
            let p = u.sequence[..k1].to_vec();
            let q = u.sequence[u.sequence.len() - k1..].to_vec();
            let rc_p = crate::dna::revcomp_bytes(&p);
            let rc_q = crate::dna::revcomp_bytes(&q);

            let key_fwd = (p, q);
            let key_rev = (rc_q, rc_p);
            let canonical_key = if key_fwd <= key_rev { key_fwd } else { key_rev };

            bubble_groups.entry(canonical_key).or_default().push(idx);
        }

        let mut haplotig_indices = hashbrown::HashSet::new();

        for members in bubble_groups.into_values() {
            if members.len() == 2 {
                let u0 = &unitigs[members[0]];
                let u1 = &unitigs[members[1]];

                let len0 = u0.sequence.len();
                let len1 = u1.sequence.len();
                let max_len = len0.max(len1) as f64;
                let min_len = len0.min(len1) as f64;

                let cov0 = u0.mean_coverage.max(0.1);
                let cov1 = u1.mean_coverage.max(0.1);
                let max_cov = cov0.max(cov1);
                let min_cov = cov0.min(cov1);

                let len_ratio = (max_len - min_len) / max_len.max(1.0);
                let cov_ratio = max_cov / min_cov;

                if len_ratio <= self.max_allele_len_diff && cov_ratio <= self.max_cov_ratio {
                    // Legitimate heterozygous diploid locus:
                    // Retain higher-coverage or longer branch as primary, and route the other to haplotigs.
                    let haplotig_idx = if cov0 >= cov1 { members[1] } else { members[0] };
                    haplotig_indices.insert(haplotig_idx);
                }
            }
        }

        let mut primary = Vec::new();
        let mut haplotigs = Vec::new();

        for (idx, u) in unitigs.into_iter().enumerate() {
            if haplotig_indices.contains(&idx) {
                haplotigs.push(u);
            } else {
                primary.push(u);
            }
        }

        (primary, haplotigs)
    }
}
