//! End-to-End De Novo Assembly Pipeline & Assembly QC.

use crate::bloom::TwoTierFilter;
use crate::dna::Kmer256;
use crate::fastq::parse_reads_from_file;
use crate::graph::CompactedGraph;
use crate::simplify::Simplifier;
use anyhow::Result;
use hashbrown::HashMap;
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::time::Instant;

/// Assembly quality control statistics.
#[derive(Debug, Clone, Default)]
pub struct AssemblyStats {
    pub total_contigs: usize,
    pub total_length: usize,
    pub max_contig_length: usize,
    pub n50: usize,
    pub l50: usize,
    pub gc_content: f64,
}

pub struct AssemblerConfig {
    pub k: usize,
    pub min_coverage: f64,
    pub min_contig_len: usize,
    pub bloom_bits: usize,
    pub error_correct: bool,
    pub is_meta: bool,
    pub is_plasmid: bool,
    pub is_rna: bool,
    pub is_sc: bool,
    pub polish: bool,
    pub long_reads: Option<Vec<std::path::PathBuf>>,
    pub prior_contigs: Option<Vec<Vec<u8>>>,
    pub skip_repeat_resolution: bool,
    pub memory_limits: Option<crate::memory::MemoryLimits>,
}

impl Default for AssemblerConfig {
    fn default() -> Self {
        Self {
            k: 31,
            min_coverage: 5.0,
            min_contig_len: 200,
            bloom_bits: 512 * 1024 * 1024, // 512M bits = 64 MB RAM
            error_correct: false,
            is_meta: false,
            is_plasmid: false,
            is_rna: false,
            is_sc: false,
            polish: true,
            long_reads: None,
            prior_contigs: None,
            skip_repeat_resolution: false,
            memory_limits: None,
        }
    }
}

pub struct AssemblyResult {
    pub contigs: Vec<crate::graph::Unitig>,
    pub scaffolds: Vec<crate::graph::Unitig>,
    pub plasmids: Vec<crate::graph::Unitig>,
    pub stats: AssemblyStats,
    pub elapsed_secs: f64,
}

/// Assembles pre-loaded reads into contigs, enabling multi-k iteration without repeated disk I/O.
pub fn run_assembly_with_loaded_reads(
    all_reads: &[Vec<u8>],
    pe_reads1: Option<&[Vec<u8>]>,
    pe_reads2: Option<&[Vec<u8>]>,
    config: &AssemblerConfig,
) -> Result<AssemblyResult> {
    let start_time = Instant::now();
    let k = config.k;

    println!("─── [Stage 2] Streaming reads into Two-Tier Bloom Filter (Memory Shield) ───");
    let filter = TwoTierFilter::new(config.bloom_bits);
    println!(
        "  Bloom filter memory: {:.2} MB",
        filter.memory_usage_bytes() as f64 / 1_048_576.0
    );

    all_reads.par_iter().for_each(|seq| {
        if seq.len() < k {
            return;
        }
        for i in 0..=(seq.len() - k) {
            if let Some(kmer) = Kmer256::from_bytes(&seq[i..i + k], k) {
                let (can, _) = kmer.canonical(k);
                filter.insert_kmer256(can);
            }
        }
    });
    println!(
        "  Bloom filter populated. (Elapsed: {:.3}s)",
        start_time.elapsed().as_secs_f64()
    );

    println!("─── [Stage 3] Building Solid K-mer Index ───");
    let num_threads = rayon::current_num_threads().max(1);
    let chunk_size = (all_reads.len() / num_threads).max(1000);
    let mut solid_maps: Vec<HashMap<Kmer256, u32>> = all_reads
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut local_counts: HashMap<Kmer256, u32> = HashMap::with_capacity(chunk.len() * 4);

            for seq in chunk {
                if seq.len() < k {
                    continue;
                }
                for i in 0..=(seq.len() - k) {
                    if let Some(kmer) = Kmer256::from_bytes(&seq[i..i + k], k) {
                        let (can, _) = kmer.canonical(k);
                        if filter.is_solid_kmer256(can) {
                            *local_counts.entry(can).or_insert(0) += 1;
                        }
                    }
                }
            }
            local_counts
        })
        .collect();

    // Drop filter immediately to free memory shield
    drop(filter);
    #[cfg(target_os = "linux")]
    unsafe {
        extern "C" {
            fn malloc_trim(pad: usize) -> i32;
        }
        malloc_trim(0);
    }
    if let Some(ref limits) = config.memory_limits {
        if let Err(e) = limits.check_headroom() {
            eprintln!("  [Memory Governor Warning] {}", e);
        }
    }

    // Merge thread-local maps by draining and popping to immediately free memory
    let mut global_counts: HashMap<Kmer256, u32> = HashMap::new();
    while let Some(mut l_counts) = solid_maps.pop() {
        for (kmer, cnt) in l_counts.drain() {
            *global_counts.entry(kmer).or_insert(0) += cnt;
        }
    }

    // Inject prior contigs for progressive multi-k continuity
    if let Some(ref priors) = config.prior_contigs {
        println!(
            "  [Multi-K] Injecting {} prior unitigs as high-confidence backbone paths...",
            priors.len()
        );
        for seq in priors {
            if seq.len() < k {
                continue;
            }
            for i in 0..=(seq.len() - k) {
                if let Some(kmer) = Kmer256::from_bytes(&seq[i..i + k], k) {
                    let (can, _) = kmer.canonical(k);
                    *global_counts.entry(can).or_insert(0) += 50;
                }
            }
        }
    }

    // Determine solid k-mer cutoff dynamically (robust noise valley detection)
    let min_kmer_cov = if config.min_coverage > 5.0 {
        (config.min_coverage * 0.2).clamp(2.0, 10.0) as u32
    } else {
        let mut sample_covs: Vec<u32> = global_counts.values().copied().collect();
        if sample_covs.len() > 100 {
            sample_covs.sort_unstable();
            let top_cov = sample_covs[sample_covs.len() * 95 / 100];
            if top_cov >= 30 {
                let max_valley_search = (top_cov / 4).clamp(10, 35) as usize;
                let mut hist = vec![0usize; max_valley_search + 1];
                for &c in &sample_covs {
                    if (c as usize) <= max_valley_search {
                        hist[c as usize] += 1;
                    }
                }
                let mut valley = 2u32;
                let mut min_val = usize::MAX;
                let mut found_valley = false;
                for (c, &h_val) in hist.iter().enumerate().take(max_valley_search + 1).skip(2) {
                    if h_val <= min_val {
                        min_val = h_val;
                        valley = c as u32;
                    } else if h_val > min_val * 2 && c > valley as usize + 2 {
                        found_valley = true;
                        break;
                    }
                }
                let cutoff = if found_valley {
                    valley
                } else {
                    2u32
                };
                if cutoff > 2 {
                    println!(
                        "  [Auto-Cutoff] Robust noise valley detected at {}x (Top: {}x)",
                        cutoff, top_cov
                    );
                }
                cutoff
            } else {
                2u32
            }
        } else {
            2u32
        }
    };

    let solid_kmers: hashbrown::HashSet<Kmer256> = global_counts
        .iter()
        .filter(|(_, &cov)| cov >= min_kmer_cov)
        .map(|(&kmer, _)| kmer)
        .collect();

    println!(
        "  Total solid k-mers retained: {} (Elapsed: {:.3}s)",
        solid_kmers.len(),
        start_time.elapsed().as_secs_f64()
    );

    println!("─── [Stage 4] Compacting de Bruijn Graph into Unitigs ───");
    let cdbg = CompactedGraph::build(k, &solid_kmers, &global_counts);
    drop(solid_kmers);
    drop(global_counts);
    #[cfg(target_os = "linux")]
    unsafe {
        extern "C" {
            fn malloc_trim(pad: usize) -> i32;
        }
        malloc_trim(0);
    }
    println!(
        "  Raw unitigs constructed: {} (Elapsed: {:.3}s)",
        cdbg.unitigs.len(),
        start_time.elapsed().as_secs_f64()
    );
    if let Some(ref limits) = config.memory_limits {
        if let Err(e) = limits.check_headroom() {
            eprintln!("  [Memory Governor Warning] {}", e);
        }
    }

    let raw_unitigs = if config.is_sc {
        println!("─── [Single-Cell Mode] Normalizing MDA Coverage Discrepancies ───");
        crate::single_cell::SingleCellNormalizer::default().normalize_coverage(cdbg.unitigs)
    } else {
        cdbg.unitigs
    };

    println!("─── [Stage 5] Simplification (Tip Clipping & Artifact Cleaning) ───");
    let mut simplifier = Simplifier::new(k, config.min_coverage, config.min_contig_len);
    simplifier.is_rna = config.is_rna;
    let simplified_contigs = simplifier.simplify(raw_unitigs);
    println!(
        "  Assembled contigs after simplification: {} (Elapsed: {:.3}s)",
        simplified_contigs.len(),
        start_time.elapsed().as_secs_f64()
    );

    println!("─── [Stage 6] ExSPAnder Repeat Resolution (Paired-End Linkages) ───");
    let contigs = if !config.skip_repeat_resolution && pe_reads1.is_some() && pe_reads2.is_some() {
        let (r1, r2) = (pe_reads1.unwrap(), pe_reads2.unwrap());
        let paired_info =
            crate::paired_info::PairedInfoIndex::build(k, &simplified_contigs, r1, r2);
        println!(
            "  Paired library estimated insert size: {:.1} ± {:.1} bp",
            paired_info.mean_insert_size, paired_info.insert_size_stdev
        );
        let expander = crate::expander::ExSPAnder::default();
        let resolved = expander.resolve_repeats(k, simplified_contigs, &paired_info);
        println!(
            "  Contigs after repeat resolution: {} (Elapsed: {:.3}s)",
            resolved.len(),
            start_time.elapsed().as_secs_f64()
        );
        resolved
    } else {
        simplified_contigs
    };

    let contigs = if !config.skip_repeat_resolution {
        if let Some(ref lr_paths) = config.long_reads {
            if !lr_paths.is_empty() {
                println!("─── [Stage 6.5] Spaligner Hybrid Long-Read Bridging ───");
                let mut lr_seqs = Vec::new();
                for p in lr_paths {
                    if let Ok(reads) = parse_reads_from_file(p) {
                        lr_seqs.extend(reads);
                    }
                }
                let bridged = crate::spaligner::LongReadResolver {
                    k,
                    min_seed_matches: 3,
                }
                .bridge_with_long_reads(contigs, &lr_seqs);
                println!(
                    "  Contigs after long-read bridging: {} (Elapsed: {:.3}s)",
                    bridged.len(),
                    start_time.elapsed().as_secs_f64()
                );
                bridged
            } else {
                contigs
            }
        } else {
            contigs
        }
    } else {
        contigs
    };

    let contigs = if config.is_rna {
        println!("─── [RNA Mode] Preserving Alternative Splicing Isoforms ───");
        crate::rna::RnaEngine::default().process_transcripts(contigs)
    } else {
        contigs
    };

    println!("─── [Stage 7] Scaffolding Across Unresolved Gaps ───");
    let scaffolds = if !config.skip_repeat_resolution && pe_reads1.is_some() && pe_reads2.is_some() {
        let (r1, r2) = (pe_reads1.unwrap(), pe_reads2.unwrap());
        let paired_info = crate::paired_info::PairedInfoIndex::build(k, &contigs, r1, r2);
        crate::scaffold::Scaffolder::default().build_scaffolds(&contigs, &paired_info)
    } else {
        contigs.clone()
    };
    println!(
        "  Scaffolds generated: {} (Elapsed: {:.3}s)",
        scaffolds.len(),
        start_time.elapsed().as_secs_f64()
    );

    let contigs = if config.polish && !config.skip_repeat_resolution {
        println!("─── [Stage 8] Consensus Base Polishing ───");
        let (polished, fixes) =
            crate::polisher::Polisher::default().polish_contigs(contigs, all_reads);
        println!(
            "  Polished {} base discrepancies (Elapsed: {:.3}s)",
            fixes,
            start_time.elapsed().as_secs_f64()
        );
        polished
    } else {
        contigs
    };

    let (contigs, plasmids) = if config.is_plasmid {
        let (chrom, plas) = crate::modes::PlasmidDetector {
            k: config.k,
            ..Default::default()
        }
        .extract_plasmids(contigs);
        println!(
            "  [Plasmid Mode] Separated {} circular/high-copy plasmids and {} chromosomal contigs",
            plas.len(),
            chrom.len()
        );
        (chrom, plas)
    } else {
        (contigs, Vec::new())
    };

    let contigs = if config.is_meta {
        let meta_contigs = crate::modes::apply_meta_filter(contigs, config.min_contig_len);
        println!(
            "  [Meta Mode] Retained {} contigs across uneven coverage depths",
            meta_contigs.len()
        );
        meta_contigs
    } else {
        contigs
    };

    // Calculate QC statistics
    let stats = calculate_stats(&contigs);
    let elapsed = start_time.elapsed().as_secs_f64();

    Ok(AssemblyResult {
        contigs,
        scaffolds,
        plasmids,
        stats,
        elapsed_secs: elapsed,
    })
}

/// Assembles contigs directly from 2-bit packed reads without intermediate read vector allocations.
pub fn run_assembly_with_packed_reads(
    packed: &crate::packed_reads::PackedReads,
    config: &AssemblerConfig,
) -> Result<AssemblyResult> {
    let start_time = Instant::now();
    let k = config.k;

    println!("─── [Stage 2] Streaming reads into Two-Tier Bloom Filter (Memory Shield) ───");
    let filter = TwoTierFilter::new(config.bloom_bits);
    println!(
        "  Bloom filter memory: {:.2} MB",
        filter.memory_usage_bytes() as f64 / 1_048_576.0
    );

    let num_reads = packed.len();
    let num_threads = rayon::current_num_threads().max(1);
    let chunk_size = (num_reads / num_threads).max(1000);
    let read_indices: Vec<usize> = (0..num_reads).collect();

    read_indices.par_chunks(chunk_size).for_each(|chunk| {
        let mut buf = Vec::with_capacity(512);
        for &idx in chunk {
            packed.get_read(idx, &mut buf);
            if buf.len() < k {
                continue;
            }
            for i in 0..=(buf.len() - k) {
                if let Some(kmer) = Kmer256::from_bytes(&buf[i..i + k], k) {
                    let (can, _) = kmer.canonical(k);
                    filter.insert_kmer256(can);
                }
            }
        }
    });
    println!(
        "  Bloom filter populated. (Elapsed: {:.3}s)",
        start_time.elapsed().as_secs_f64()
    );

    println!("─── [Stage 3] Building Solid K-mer Index ───");
    const NUM_SHARDS: usize = 64;
    let shards: Vec<std::sync::Mutex<HashMap<Kmer256, u32>>> = (0..NUM_SHARDS)
        .map(|_| std::sync::Mutex::new(HashMap::with_capacity(75_000)))
        .collect();

    read_indices.par_chunks(chunk_size).for_each(|chunk| {
        let mut thread_buffers: Vec<Vec<Kmer256>> = (0..NUM_SHARDS)
            .map(|_| Vec::with_capacity(256))
            .collect();
        let mut buf = Vec::with_capacity(512);

        for &idx in chunk {
            packed.get_read(idx, &mut buf);
            if buf.len() < k {
                continue;
            }
            for i in 0..=(buf.len() - k) {
                if let Some(kmer) = Kmer256::from_bytes(&buf[i..i + k], k) {
                    let (can, _) = kmer.canonical(k);
                    if filter.is_solid_kmer256(can) {
                        let shard_idx = (can.0 as usize) % NUM_SHARDS;
                        thread_buffers[shard_idx].push(can);
                        if thread_buffers[shard_idx].len() >= 256 {
                            let mut guard = shards[shard_idx].lock().unwrap();
                            for km in thread_buffers[shard_idx].drain(..) {
                                *guard.entry(km).or_insert(0) += 1;
                            }
                        }
                    }
                }
            }
        }

        // Flush remaining thread buffers
        for (shard_idx, tb) in thread_buffers.iter_mut().enumerate() {
            if !tb.is_empty() {
                let mut guard = shards[shard_idx].lock().unwrap();
                for km in tb.drain(..) {
                    *guard.entry(km).or_insert(0) += 1;
                }
            }
        }
    });

    // Drop filter immediately
    drop(filter);
    #[cfg(target_os = "linux")]
    unsafe {
        extern "C" {
            fn malloc_trim(pad: usize) -> i32;
        }
        malloc_trim(0);
    }
    if let Some(ref limits) = config.memory_limits {
        if let Err(e) = limits.check_headroom() {
            eprintln!("  [Memory Governor Warning] {}", e);
        }
    }

    // Merge non-overlapping shards into global_counts (zero key collisions!)
    let mut global_counts: HashMap<Kmer256, u32> = HashMap::with_capacity(5_000_000);
    for shard in shards {
        let mut map = shard.into_inner().unwrap();
        global_counts.extend(map.drain());
    }

    // Inject prior contigs for progressive multi-k continuity
    if let Some(ref priors) = config.prior_contigs {
        println!(
            "  [Multi-K] Injecting {} prior unitigs as high-confidence backbone paths...",
            priors.len()
        );
        for seq in priors {
            if seq.len() < k {
                continue;
            }
            for i in 0..=(seq.len() - k) {
                if let Some(kmer) = Kmer256::from_bytes(&seq[i..i + k], k) {
                    let (can, _) = kmer.canonical(k);
                    *global_counts.entry(can).or_insert(0) += 50;
                }
            }
        }
    }

    // Determine solid k-mer cutoff dynamically (robust noise valley detection)
    let min_kmer_cov = if config.min_coverage > 5.0 {
        (config.min_coverage * 0.2).clamp(2.0, 10.0) as u32
    } else {
        let mut sample_covs: Vec<u32> = global_counts.values().copied().collect();
        if sample_covs.len() > 100 {
            sample_covs.sort_unstable();
            let top_cov = sample_covs[sample_covs.len() * 95 / 100];
            if top_cov >= 30 {
                let max_valley_search = (top_cov / 4).clamp(10, 35) as usize;
                let mut hist = vec![0usize; max_valley_search + 1];
                for &c in &sample_covs {
                    if (c as usize) <= max_valley_search {
                        hist[c as usize] += 1;
                    }
                }
                let mut valley = 2u32;
                let mut min_val = usize::MAX;
                let mut found_valley = false;
                for (c, &h_val) in hist.iter().enumerate().take(max_valley_search + 1).skip(2) {
                    if h_val <= min_val {
                        min_val = h_val;
                        valley = c as u32;
                    } else if h_val > min_val * 2 && c > valley as usize + 2 {
                        found_valley = true;
                        break;
                    }
                }
                let cutoff = if found_valley {
                    valley
                } else {
                    2u32
                };
                if cutoff > 2 {
                    println!(
                        "  [Auto-Cutoff] Robust noise valley detected at {}x (Top: {}x)",
                        cutoff, top_cov
                    );
                }
                cutoff
            } else {
                2u32
            }
        } else {
            2u32
        }
    };

    let solid_kmers: hashbrown::HashSet<Kmer256> = global_counts
        .iter()
        .filter(|(_, &cov)| cov >= min_kmer_cov)
        .map(|(&kmer, _)| kmer)
        .collect();

    println!(
        "  Total solid k-mers retained: {} (Elapsed: {:.3}s)",
        solid_kmers.len(),
        start_time.elapsed().as_secs_f64()
    );

    println!("─── [Stage 4] Compacting de Bruijn Graph into Unitigs ───");
    let cdbg = CompactedGraph::build(k, &solid_kmers, &global_counts);
    drop(solid_kmers);
    drop(global_counts);
    #[cfg(target_os = "linux")]
    unsafe {
        extern "C" {
            fn malloc_trim(pad: usize) -> i32;
        }
        malloc_trim(0);
    }
    println!(
        "  Raw unitigs constructed: {} (Elapsed: {:.3}s)",
        cdbg.unitigs.len(),
        start_time.elapsed().as_secs_f64()
    );
    if let Some(ref limits) = config.memory_limits {
        if let Err(e) = limits.check_headroom() {
            eprintln!("  [Memory Governor Warning] {}", e);
        }
    }

    let raw_unitigs = if config.is_sc {
        println!("─── [Single-Cell Mode] Normalizing MDA Coverage Discrepancies ───");
        crate::single_cell::SingleCellNormalizer::default().normalize_coverage(cdbg.unitigs)
    } else {
        cdbg.unitigs
    };

    println!("─── [Stage 5] Simplification (Tip Clipping & Artifact Cleaning) ───");
    let mut simplifier = Simplifier::new(k, config.min_coverage, config.min_contig_len);
    simplifier.is_rna = config.is_rna;
    let simplified_contigs = simplifier.simplify(raw_unitigs);
    println!(
        "  Assembled contigs after simplification: {} (Elapsed: {:.3}s)",
        simplified_contigs.len(),
        start_time.elapsed().as_secs_f64()
    );

    println!("─── [Stage 6] ExSPAnder Repeat Resolution (Paired-End Linkages) ───");
    let contigs = if !config.skip_repeat_resolution && packed.pe_boundary > 0 {
        let paired_info =
            crate::paired_info::PairedInfoIndex::build_from_packed(k, &simplified_contigs, packed);
        println!(
            "  Paired library estimated insert size: {:.1} ± {:.1} bp",
            paired_info.mean_insert_size, paired_info.insert_size_stdev
        );
        let expander = crate::expander::ExSPAnder::default();
        let resolved = expander.resolve_repeats(k, simplified_contigs, &paired_info);
        println!(
            "  Contigs after repeat resolution: {} (Elapsed: {:.3}s)",
            resolved.len(),
            start_time.elapsed().as_secs_f64()
        );
        resolved
    } else {
        simplified_contigs
    };

    let contigs = if !config.skip_repeat_resolution {
        if let Some(ref lr_paths) = config.long_reads {
            if !lr_paths.is_empty() {
                println!("─── [Stage 6.5] Spaligner Hybrid Long-Read Bridging ───");
                let mut lr_seqs = Vec::new();
                for p in lr_paths {
                    if let Ok(reads) = parse_reads_from_file(p) {
                        lr_seqs.extend(reads);
                    }
                }
                let bridged = crate::spaligner::LongReadResolver {
                    k,
                    min_seed_matches: 3,
                }
                .bridge_with_long_reads(contigs, &lr_seqs);
                println!(
                    "  Contigs after long-read bridging: {} (Elapsed: {:.3}s)",
                    bridged.len(),
                    start_time.elapsed().as_secs_f64()
                );
                bridged
            } else {
                contigs
            }
        } else {
            contigs
        }
    } else {
        contigs
    };

    let contigs = if config.is_rna {
        println!("─── [RNA Mode] Preserving Alternative Splicing Isoforms ───");
        crate::rna::RnaEngine::default().process_transcripts(contigs)
    } else {
        contigs
    };

    println!("─── [Stage 7] Scaffolding Across Unresolved Gaps ───");
    let scaffolds = if !config.skip_repeat_resolution && packed.pe_boundary > 0 {
        let paired_info = crate::paired_info::PairedInfoIndex::build_from_packed(k, &contigs, packed);
        crate::scaffold::Scaffolder::default().build_scaffolds(&contigs, &paired_info)
    } else {
        contigs.clone()
    };
    println!(
        "  Scaffolds generated: {} (Elapsed: {:.3}s)",
        scaffolds.len(),
        start_time.elapsed().as_secs_f64()
    );

    let contigs = if config.polish && !config.skip_repeat_resolution {
        println!("─── [Stage 8] Consensus Base Polishing ───");
        let (polished, fixes) =
            crate::polisher::Polisher::default().polish_contigs_packed(contigs, packed);
        println!(
            "  Polished {} base discrepancies (Elapsed: {:.3}s)",
            fixes,
            start_time.elapsed().as_secs_f64()
        );
        polished
    } else {
        contigs
    };

    let (contigs, plasmids) = if config.is_plasmid {
        let (chrom, plas) = crate::modes::PlasmidDetector {
            k: config.k,
            ..Default::default()
        }
        .extract_plasmids(contigs);
        println!(
            "  [Plasmid Mode] Separated {} circular/high-copy plasmids and {} chromosomal contigs",
            plas.len(),
            chrom.len()
        );
        (chrom, plas)
    } else {
        (contigs, Vec::new())
    };

    let contigs = if config.is_meta {
        let meta_contigs = crate::modes::apply_meta_filter(contigs, config.min_contig_len);
        println!(
            "  [Meta Mode] Retained {} contigs across uneven coverage depths",
            meta_contigs.len()
        );
        meta_contigs
    } else {
        contigs
    };

    let stats = calculate_stats(&contigs);
    let elapsed = start_time.elapsed().as_secs_f64();

    Ok(AssemblyResult {
        contigs,
        scaffolds,
        plasmids,
        stats,
        elapsed_secs: elapsed,
    })
}

/// Assembles input FASTQ/FASTA files into contigs.
pub fn run_assembly<P: AsRef<Path> + Sync>(
    input_files: &[P],
    config: &AssemblerConfig,
) -> Result<AssemblyResult> {
    println!(
        "─── [Stage 1] Ingesting and packing reads from {} input file(s) ───",
        input_files.len()
    );
    let packed = crate::packed_reads::PackedReads::from_files(input_files)?;
    println!(
        "  Total reads packed: {} ({:.2} MB in RAM)",
        packed.len(),
        packed.memory_usage_bytes() as f64 / 1_048_576.0
    );

    run_assembly_with_packed_reads(&packed, config)
}

/// Computes assembly QC metrics: N50, L50, Max length, total bases, GC%.
pub fn calculate_stats(contigs: &[crate::graph::Unitig]) -> AssemblyStats {
    if contigs.is_empty() {
        return AssemblyStats::default();
    }

    let mut lengths: Vec<usize> = contigs.iter().map(|c| c.sequence.len()).collect();
    lengths.sort_unstable_by(|a, b| b.cmp(a));

    let total_length: usize = lengths.iter().sum();
    let max_contig_length = lengths[0];

    let half_total = total_length / 2;
    let mut cumulative = 0;
    let mut n50 = 0;
    let mut l50 = 0;

    for (i, &len) in lengths.iter().enumerate() {
        cumulative += len;
        if cumulative >= half_total {
            n50 = len;
            l50 = i + 1;
            break;
        }
    }

    let mut gc_count = 0usize;
    for c in contigs {
        for &b in &c.sequence {
            if b == b'G' || b == b'g' || b == b'C' || b == b'c' {
                gc_count += 1;
            }
        }
    }
    let gc_content = if total_length > 0 {
        (gc_count as f64 / total_length as f64) * 100.0
    } else {
        0.0
    };

    AssemblyStats {
        total_contigs: contigs.len(),
        total_length,
        max_contig_length,
        n50,
        l50,
        gc_content,
    }
}

/// Writes contigs to a FASTA file.
pub fn write_contigs_fasta<P: AsRef<Path>>(
    contigs: &[crate::graph::Unitig],
    out_path: P,
) -> Result<()> {
    let file = File::create(out_path)?;
    let mut writer = std::io::BufWriter::with_capacity(1024 * 1024, file);
    for (i, c) in contigs.iter().enumerate() {
        writeln!(
            writer,
            ">contig_{}_len_{}_cov_{:.1}",
            i + 1,
            c.sequence.len(),
            c.mean_coverage
        )?;
        // Write 80 characters per line
        for chunk in c.sequence.chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }
    writer.flush()?;
    Ok(())
}
