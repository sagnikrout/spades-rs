//! Multi-k-mer Iteration Engine (Module 2.1).
//!
//! Emulates SPAdes progressive multi-k stepping (e.g. k=21 -> 33 -> 55).
//! Resolves repeats progressively using increasing k-mer lengths while preserving
//! connectivity from smaller k-mers.

use crate::assemble::{run_assembly_with_packed_reads, AssemblerConfig, AssemblyResult};
use crate::packed_reads::PackedReads;
use anyhow::Result;
use std::path::{Path, PathBuf};

/// Configuration parameters for progressive multi-$k$-mer de novo assembly.
pub struct MultiKConfig {
    /// Ordered sequence of odd $k$-mer sizes (e.g. `[21, 33, 55]`) stepped through progressively.
    pub kmers: Vec<usize>,
    /// Minimum mean $k$-mer coverage required to retain unitigs during graph simplification.
    pub min_coverage: f64,
    /// Minimum sequence length (in base pairs) required for output contigs and scaffolds.
    pub min_contig_len: usize,
    /// Bit capacity of each tier of the two-tier Bloom filter used to screen singleton $k$-mers.
    pub bloom_bits: usize,
    /// Whether to run BayesHammer-style single-base read error correction before graph construction.
    pub error_correct: bool,
    /// Enables metagenomic assembly mode (`metaSPAdes`).
    pub is_meta: bool,
    /// Enables plasmid extraction mode (`plasmidSPAdes`).
    pub is_plasmid: bool,
    /// Enables transcriptome assembly mode (`rnaSPAdes`).
    pub is_rna: bool,
    /// Enables single-cell MDA assembly mode (`scSPAdes`).
    pub is_sc: bool,
    /// Whether to run post-assembly read-alignment consensus polishing on the final $k$-mer step.
    pub polish: bool,
    /// Enables mismatch-intolerant diploid-aware consensus polishing (`--careful`).
    pub careful: bool,
    /// Indicates that each input FASTQ file contains interleaved forward (`R1`) and reverse (`R2`) mates.
    pub interleaved: bool,
    /// Optional paths to long-read FASTQ/FASTA files (PacBio or Oxford Nanopore) for hybrid repeat resolution.
    pub long_reads: Option<Vec<PathBuf>>,
    /// Optional paths to barcoded linked-read FASTQ files (10x Genomics or TELL-Seq) for long-range scaffolding.
    pub linked_reads: Option<Vec<PathBuf>>,
    /// Optional paths to trusted reference or prior assembly FASTA contigs injected into the graph.
    pub trusted_contigs: Option<Vec<PathBuf>>,
    /// Optional system memory budget and hardware limits.
    pub memory_limits: Option<crate::memory::MemoryLimits>,
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
            careful: false,
            interleaved: false,
            long_reads: None,
            linked_reads: None,
            trusted_contigs: None,
            memory_limits: None,
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
    if config.kmers.is_empty() {
        anyhow::bail!("At least one k-mer must be specified for multi-k assembly");
    }
    println!("  Target k-mer steps: {:?}", config.kmers);

    println!(
        "─── [Multi-K Preload] Ingesting and packing reads from {} input file(s) ───",
        input_files.len()
    );
    let packed = if config.interleaved {
        PackedReads::from_interleaved_files(input_files)?
    } else {
        PackedReads::from_files(input_files)?
    };
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
            error_correct: config.error_correct,
            is_meta: config.is_meta,
            is_plasmid: config.is_plasmid && is_last_step,
            is_rna: config.is_rna,
            is_sc: config.is_sc,
            polish: config.polish && is_last_step, // Polish on final assembly
            careful: config.careful && is_last_step,
            interleaved: config.interleaved,
            long_reads: if is_last_step {
                config.long_reads.clone()
            } else {
                None
            },
            linked_reads: if is_last_step {
                config.linked_reads.clone()
            } else {
                None
            },
            trusted_contigs: config.trusted_contigs.clone(),
            prior_contigs: prior_contigs.clone(),
            skip_repeat_resolution: !is_last_step,
            memory_limits: config.memory_limits,
        };

        let result = run_assembly_with_packed_reads(&packed, &sub_config)?;
        println!(
            "  ✓ Step k={} completed: {} contigs, max length {} bp, N50 {} bp (Elapsed: {:.3}s)",
            k,
            result.stats.total_contigs,
            result.stats.max_contig_length,
            result.stats.n50,
            result.elapsed_secs
        );

        // Save high-confidence intermediate unitigs from earlier steps for potential rescue
        if !is_last_step {
            for u in &result.contigs {
                if u.sequence.len() >= 300 && u.mean_coverage >= config.min_coverage {
                    intermediate_unitigs.push(u.clone());
                }
            }
        }

        // Forward assembled unitigs into next step
        prior_contigs = Some(result.contigs.iter().map(|u| u.sequence.clone()).collect());
        final_result = Some(result);
        crate::memory::trim_memory();
    }

    let mut result =
        final_result.ok_or_else(|| anyhow::anyhow!("At least one k-mer must be specified"))?;

    // Rescue valid unitigs from lower-k steps that were dropped due to higher-k coverage thinning
    if !intermediate_unitigs.is_empty() {
        // Use the actual final k for fingerprinting, capped at 32 because string_to_kmer /
        // canonical_kmer_u64 pack into u64 and only support k ≤ 32.
        let rescue_k = config.kmers.last().copied().unwrap_or(31).min(32);

        let mut final_kmers: hashbrown::HashSet<u64> = hashbrown::HashSet::new();
        for u in result.contigs.iter().chain(result.plasmids.iter()) {
            if u.sequence.len() >= rescue_k {
                for i in 0..=(u.sequence.len() - rescue_k) {
                    if let Some(km) =
                        crate::dna::string_to_kmer(&u.sequence[i..i + rescue_k], rescue_k)
                    {
                        let (can, _) = crate::dna::canonical_kmer_u64(km, rescue_k);
                        final_kmers.insert(can);
                    }
                }
            }
        }

        let mut rescued_count = 0;
        let k = rescue_k;
        for prior_u in intermediate_unitigs.into_iter().rev() {
            if prior_u.sequence.len() < 300 {
                continue;
            }
            let mut missing_kmers = 0;
            let mut total_kmers = 0;
            let mut u_kmers = Vec::new();
            for i in 0..=(prior_u.sequence.len() - k) {
                if let Some(km) = crate::dna::string_to_kmer(&prior_u.sequence[i..i + k], k) {
                    let (can, _) = crate::dna::canonical_kmer_u64(km, k);
                    total_kmers += 1;
                    u_kmers.push(can);
                    if !final_kmers.contains(&can) {
                        missing_kmers += 1;
                    }
                }
            }
            if total_kmers > 0 && (missing_kmers as f64 / total_kmers as f64) >= 0.5 {
                final_kmers.extend(u_kmers);
                let mut rescued_u = prior_u;
                rescued_u.id = result.contigs.len();
                result.scaffolds.push(rescued_u.clone());
                result.contigs.push(rescued_u);
                rescued_count += 1;
            }
        }
        if rescued_count > 0 {
            println!(
                "  [Multi-K Rescue] Rescued {} valid unitigs from lower-k steps dropped at final k",
                rescued_count
            );
            result
                .contigs
                .sort_by_key(|u| std::cmp::Reverse(u.sequence.len()));
            for (idx, u) in result.contigs.iter_mut().enumerate() {
                u.id = idx;
            }
            result
                .scaffolds
                .sort_by_key(|u| std::cmp::Reverse(u.sequence.len()));
            for (idx, u) in result.scaffolds.iter_mut().enumerate() {
                u.id = idx;
            }
            result.stats = crate::assemble::calculate_stats(&result.contigs);
        }
    }

    Ok(result)
}
