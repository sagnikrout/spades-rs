//! Multi-k-mer Iteration Engine (Module 2.1).
//!
//! Emulates SPAdes progressive multi-k stepping (e.g. k=21 -> 33 -> 55).
//! Resolves repeats progressively using increasing k-mer lengths while preserving
//! connectivity from smaller k-mers.

use crate::assemble::{run_assembly_with_packed_reads, AssemblerConfig, AssemblyResult};
use crate::packed_reads::PackedReads;
use anyhow::Result;
use std::path::{Path, PathBuf};

pub struct MultiKConfig {
    pub kmers: Vec<usize>,
    pub min_coverage: f64,
    pub min_contig_len: usize,
    pub bloom_bits: usize,
    pub error_correct: bool,
    pub is_meta: bool,
    pub is_plasmid: bool,
    pub is_rna: bool,
    pub is_sc: bool,
    pub polish: bool,
    pub long_reads: Option<Vec<PathBuf>>,
}

impl Default for MultiKConfig {
    fn default() -> Self {
        Self {
            kmers: vec![21, 33, 55],
            min_coverage: 5.0,
            min_contig_len: 200,
            bloom_bits: 512 * 1024 * 1024,
            error_correct: false,
            is_meta: false,
            is_plasmid: false,
            is_rna: false,
            is_sc: false,
            polish: true,
            long_reads: None,
        }
    }
}

/// Runs iterative multi-k assembly across multiple k-mer values using 2-bit packed reads.
pub fn run_multik_assembly<P: AsRef<Path> + Sync>(
    input_files: &[P],
    config: &MultiKConfig,
) -> Result<AssemblyResult> {
    println!("===========================================================");
    println!("       MULTI-K-MER PROGRESSIVE ASSEMBLY ENGINE             ");
    println!("===========================================================");
    println!("  Target k-mer steps: {:?}", config.kmers);

    println!(
        "─── [Multi-K Preload] Ingesting and packing reads from {} input file(s) ───",
        input_files.len()
    );
    let packed = PackedReads::from_files(input_files)?;
    println!(
        "  Total reads packed for Multi-K: {} ({:.2} MB in RAM)",
        packed.len(),
        packed.memory_usage_bytes() as f64 / 1_048_576.0
    );

    let mut final_result: Option<AssemblyResult> = None;
    let mut prior_contigs: Option<Vec<Vec<u8>>> = None;
    let mut intermediate_unitigs: Vec<crate::graph::Unitig> = Vec::new();

    for (step, &k) in config.kmers.iter().enumerate() {
        let is_last_step = step + 1 == config.kmers.len();

        println!(
            "\n▶▶▶ [Multi-K Step {}/{}] Assembling with k = {} (Backbone priors: {}) ◀◀◀",
            step + 1,
            config.kmers.len(),
            k,
            prior_contigs.as_ref().map(|p| p.len()).unwrap_or(0)
        );

        let sub_config = AssemblerConfig {
            k,
            min_coverage: config.min_coverage,
            min_contig_len: config.min_contig_len,
            bloom_bits: config.bloom_bits,
            error_correct: config.error_correct && step == 0, // Error correct once at the beginning
            is_meta: config.is_meta,
            is_plasmid: config.is_plasmid,
            is_rna: config.is_rna,
            is_sc: config.is_sc,
            polish: config.polish && is_last_step, // Polish on final assembly
            long_reads: if is_last_step { config.long_reads.clone() } else { None },
            prior_contigs: prior_contigs.clone(),
            skip_repeat_resolution: !is_last_step,
        };

        let result = run_assembly_with_packed_reads(&packed, &sub_config)?;
        println!(
            "  ✓ Step k={} completed: {} contigs, max length {} bp, N50 {} bp (Elapsed: {:.3}s)",
            k, result.stats.total_contigs, result.stats.max_contig_length, result.stats.n50, result.elapsed_secs
        );

        // Save high-confidence intermediate unitigs from earlier steps for potential rescue
        if !is_last_step {
            for u in &result.contigs {
                if u.sequence.len() >= 500 && u.mean_coverage >= config.min_coverage {
                    intermediate_unitigs.push(u.clone());
                }
            }
        }

        // Forward assembled unitigs into next step
        prior_contigs = Some(result.contigs.iter().map(|u| u.sequence.clone()).collect());
        final_result = Some(result);

        #[cfg(target_os = "linux")]
        unsafe {
            extern "C" {
                fn malloc_trim(pad: usize) -> i32;
            }
            malloc_trim(0);
        }
    }

    let mut result = final_result.expect("At least one k-mer must be specified");

    // Rescue valid unitigs from lower-k steps that were dropped due to higher-k coverage thinning
    if !intermediate_unitigs.is_empty() {
        let final_kmers: hashbrown::HashSet<u64> = result
            .contigs
            .iter()
            .flat_map(|u| {
                let k = 31;
                if u.sequence.len() < k {
                    vec![]
                } else {
                    let mut kmers = Vec::with_capacity(u.sequence.len() - k + 1);
                    for i in 0..=(u.sequence.len() - k) {
                        if let Some(km) = crate::dna::string_to_kmer(&u.sequence[i..i + k], k) {
                            let (can, _) = crate::dna::canonical_kmer_u64(km, k);
                            kmers.push(can);
                        }
                    }
                    kmers
                }
            })
            .collect();

        let mut rescued_count = 0;
        let k = 31;
        for prior_u in intermediate_unitigs {
            if prior_u.sequence.len() < 500 {
                continue;
            }
            let mut missing_kmers = 0;
            let mut total_kmers = 0;
            for i in 0..=(prior_u.sequence.len() - k) {
                if let Some(km) = crate::dna::string_to_kmer(&prior_u.sequence[i..i + k], k) {
                    let (can, _) = crate::dna::canonical_kmer_u64(km, k);
                    total_kmers += 1;
                    if !final_kmers.contains(&can) {
                        missing_kmers += 1;
                    }
                }
            }
            if total_kmers > 0 && (missing_kmers as f64 / total_kmers as f64) >= 0.5 {
                result.contigs.push(prior_u);
                rescued_count += 1;
            }
        }
        if rescued_count > 0 {
            println!(
                "  [Multi-K Rescue] Rescued {} valid unitigs from lower-k steps dropped at final k",
                rescued_count
            );
            result.contigs.sort_by_key(|u| std::cmp::Reverse(u.sequence.len()));
            result.stats = crate::assemble::calculate_stats(&result.contigs);
        }
    }

    Ok(result)
}
