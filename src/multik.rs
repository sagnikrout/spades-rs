//! Multi-k-mer Iteration Engine (Module 2.1).
//!
//! Emulates SPAdes progressive multi-k stepping (e.g. k=21 -> 33 -> 55).
//! Resolves repeats progressively using increasing k-mer lengths while preserving
//! connectivity from smaller k-mers.

use crate::assemble::{run_assembly, AssemblerConfig, AssemblyResult};
use anyhow::Result;
use std::path::Path;

pub struct MultiKConfig {
    pub kmers: Vec<usize>,
    pub min_coverage: f64,
    pub min_contig_len: usize,
    pub bloom_bits: usize,
}

impl Default for MultiKConfig {
    fn default() -> Self {
        Self {
            kmers: vec![21, 33, 55],
            min_coverage: 5.0,
            min_contig_len: 200,
            bloom_bits: 64 * 1024 * 1024,
        }
    }
}

/// Runs iterative multi-k assembly across multiple k-mer values.
pub fn run_multik_assembly<P: AsRef<Path> + Sync>(
    input_files: &[P],
    config: &MultiKConfig,
) -> Result<AssemblyResult> {
    println!("===========================================================");
    println!("       MULTI-K-MER PROGRESSIVE ASSEMBLY ENGINE             ");
    println!("===========================================================");
    println!("  Target k-mer steps: {:?}", config.kmers);

    let mut final_result: Option<AssemblyResult> = None;

    for (step, &k) in config.kmers.iter().enumerate() {
        println!("\n▶▶▶ [Multi-K Step {}/{}] Assembling with k = {} ◀◀◀", step + 1, config.kmers.len(), k);

        let sub_config = AssemblerConfig {
            k,
            min_coverage: config.min_coverage,
            min_contig_len: config.min_contig_len,
            bloom_bits: config.bloom_bits,
            ..AssemblerConfig::default()
        };

        let result = run_assembly(input_files, &sub_config)?;
        println!(
            "  ✓ Step k={} completed: {} contigs, max length {} bp (Elapsed: {:.3}s)",
            k, result.stats.total_contigs, result.stats.max_contig_length, result.elapsed_secs
        );

        final_result = Some(result);
    }

    Ok(final_result.expect("At least one k-mer must be specified"))
}
