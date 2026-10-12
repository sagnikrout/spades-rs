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
    /// Total number of assembled contigs passing the minimum length filter.
    pub total_contigs: usize,
    /// Cumulative length (in base pairs) of all assembled contigs.
    pub total_length: usize,
    /// Length (in base pairs) of the longest assembled contig.
    pub max_contig_length: usize,
    /// Sequence length at which 50% of the total assembly length is contained in contigs of at least this size.
    pub n50: usize,
    /// Minimum number of contigs whose summed length is at least 50% of the total assembly length.
    pub l50: usize,
    /// Percentage of unambiguous bases (`A`, `C`, `G`, `T`) that are `G` or `C` (0.0 to 100.0).
    pub gc_content: f64,
}

/// Configuration parameters for a single-$k$ de novo assembly run.
#[derive(Clone, Debug)]
pub struct AssemblerConfig {
    /// Odd $k$-mer size ($15 \le k \le 127$) used for de Bruijn graph construction.
    pub k: usize,
    /// Minimum mean $k$-mer coverage required to retain unitigs during graph simplification.
    pub min_coverage: f64,
    /// Minimum sequence length (in base pairs) required for output contigs and scaffolds.
    pub min_contig_len: usize,
    /// Bit capacity of each tier of the two-tier Bloom filter used to screen singleton $k$-mers.
    pub bloom_bits: usize,
    /// Whether to run BayesHammer-style single-base read error correction before graph construction.
    pub error_correct: bool,
    /// Enables metagenomic assembly mode (`metaSPAdes`), preserving low-abundance species contigs.
    pub is_meta: bool,
    /// Enables plasmid extraction mode (`plasmidSPAdes`), isolating circular and high-copy extrachromosomal contigs.
    pub is_plasmid: bool,
    /// Enables transcriptome assembly mode (`rnaSPAdes`), preserving alternative splicing isoforms.
    pub is_rna: bool,
    /// Enables single-cell MDA assembly mode (`scSPAdes`), normalizing extreme amplification spikes and dropouts.
    pub is_sc: bool,
    /// Whether to run post-assembly read-alignment consensus polishing.
    pub polish: bool,
    /// Enables mismatch-intolerant diploid-aware consensus polishing (`--careful`).
    pub careful: bool,
    /// Indicates that each input FASTQ file contains interleaved forward (`R1`) and reverse (`R2`) mates.
    pub interleaved: bool,
    /// Optional paths to long-read FASTQ/FASTA files (PacBio or Oxford Nanopore) for hybrid repeat resolution.
    pub long_reads: Option<Vec<std::path::PathBuf>>,
    /// Optional paths to barcoded linked-read FASTQ files (10x Genomics or TELL-Seq) for long-range scaffolding.
    pub linked_reads: Option<Vec<std::path::PathBuf>>,
    /// Optional paths to trusted reference or prior assembly FASTA contigs injected into the graph.
    pub trusted_contigs: Option<Vec<std::path::PathBuf>>,
    /// Contigs carried forward from smaller $k$-mer iterations during multi-$k$ assembly.
    pub prior_contigs: Option<Vec<Vec<u8>>>,
    /// Skips paired-end repeat resolution and scaffolding (used on intermediate multi-$k$ steps).
    pub skip_repeat_resolution: bool,
    /// Optional system memory budget and hardware limits.
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
            careful: false,
            interleaved: false,
            long_reads: None,
            linked_reads: None,
            trusted_contigs: None,
            prior_contigs: None,
            skip_repeat_resolution: false,
            memory_limits: None,
        }
    }
}

/// Output contigs, scaffolds, extracted plasmids, and summary statistics from an assembly run.
pub struct AssemblyResult {
    /// Assembled contigs sorted by descending sequence length.
    pub contigs: Vec<crate::graph::Unitig>,
    /// Assembled scaffolds bridged across paired-end gaps, sorted by descending sequence length.
    pub scaffolds: Vec<crate::graph::Unitig>,
    /// Extrachromosomal plasmid contigs isolated when `is_plasmid` is enabled.
    pub plasmids: Vec<crate::graph::Unitig>,
    /// Contiguity and composition metrics computed over `contigs`.
    pub stats: AssemblyStats,
    /// Total wall-clock runtime of the assembly pipeline in seconds.
    pub elapsed_secs: f64,
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
    let num_chunks = num_reads.div_ceil(chunk_size);

    (0..num_chunks).into_par_iter().for_each(|c_idx| {
        let start = c_idx * chunk_size;
        let end = (start + chunk_size).min(num_reads);
        let mut buf = Vec::with_capacity(512);
        for idx in start..end {
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

    if config.error_correct {
        println!("─── [Stage 2.5] BayesHammer Single-Base Solid K-mer Error Correction ───");
    }
    println!("─── [Stage 3] Building Solid K-mer Index ───");
    let num_shards = (num_threads * 4).next_power_of_two().clamp(8, 64);
    let shard_mask = num_shards - 1;
    let shards: Vec<std::sync::Mutex<HashMap<Kmer256, u32>>> = (0..num_shards)
        .map(|_| std::sync::Mutex::new(HashMap::with_capacity(75_000)))
        .collect();
    let corrected_reads_count = std::sync::atomic::AtomicUsize::new(0);
    let corrector = crate::hammer::ErrorCorrector::new(k);

    (0..num_chunks).into_par_iter().for_each(|c_idx| {
        let start = c_idx * chunk_size;
        let end = (start + chunk_size).min(num_reads);
        let mut thread_buffers: Vec<Vec<Kmer256>> =
            (0..num_shards).map(|_| Vec::with_capacity(256)).collect();
        let mut buf = Vec::with_capacity(512);
        let mut scratch = Vec::with_capacity(512);

        for idx in start..end {
            packed.get_read(idx, &mut buf);
            if buf.len() < k {
                continue;
            }
            if config.error_correct
                && corrector.correct_read_in_place(&mut buf, &mut scratch, &filter)
            {
                corrected_reads_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            }
            for i in 0..=(buf.len() - k) {
                if let Some(kmer) = Kmer256::from_bytes(&buf[i..i + k], k) {
                    let (can, _) = kmer.canonical(k);
                    if filter.is_solid_kmer256(can) {
                        let shard_idx = (can.0 as usize) & shard_mask;
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

    if config.error_correct {
        println!(
            "  [BayesHammer] Corrected single-nucleotide errors in {} reads",
            corrected_reads_count.load(std::sync::atomic::Ordering::Relaxed)
        );
    }

    // Drop filter immediately
    drop(filter);
    crate::memory::trim_memory();
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

    // Inject trusted contigs and prior contigs for progressive multi-k continuity
    let mut all_priors: Vec<Vec<u8>> = config.prior_contigs.clone().unwrap_or_default();
    if let Some(ref trusted) = config.trusted_contigs {
        for tp in trusted {
            match parse_reads_from_file(tp) {
                Ok(seqs) => {
                    println!(
                        "  [Trusted Contigs] Loaded {} high-confidence backbone priors from {:?}",
                        seqs.len(),
                        tp
                    );
                    all_priors.extend(seqs);
                }
                Err(e) => {
                    eprintln!(
                        "  [Warning] Failed to parse trusted contigs from {:?}: {}",
                        tp, e
                    );
                }
            }
        }
    }

    if !all_priors.is_empty() {
        println!(
            "  [Multi-K] Injecting {} prior unitigs as high-confidence backbone paths...",
            all_priors.len()
        );
        let prior_boost = (config.min_coverage as u32 * 4).max(50);
        for seq in &all_priors {
            if seq.len() < k {
                continue;
            }
            for i in 0..=(seq.len() - k) {
                if let Some(kmer) = Kmer256::from_bytes(&seq[i..i + k], k) {
                    let (can, _) = kmer.canonical(k);
                    *global_counts.entry(can).or_insert(0) += prior_boost;
                }
            }
        }
    }

    let min_kmer_cov = determine_solid_kmer_cutoff(&global_counts, config.min_coverage);

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
    crate::memory::trim_memory();
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
        crate::modes::SingleCellNormalizer::default().normalize_coverage(cdbg.unitigs)
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
    let contigs = if !config.skip_repeat_resolution && packed.is_paired() {
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
                    match parse_reads_from_file(p) {
                        Ok(reads) => lr_seqs.extend(reads),
                        Err(e) => eprintln!("  Warning: Failed to read long reads {:?}: {}", p, e),
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

    let contigs = if !config.skip_repeat_resolution {
        if let Some(ref lr_paths) = config.linked_reads {
            if !lr_paths.is_empty() {
                println!("─── [Stage 6.6] SpLitteR Linked-Read Barcode Repeat Resolution ───");
                let mut all_linked = Vec::new();
                for p in lr_paths {
                    match crate::fastq::parse_linked_reads_from_file(p) {
                        Ok(records) => all_linked.extend(records),
                        Err(e) => {
                            eprintln!("  Warning: Failed to read linked reads {:?}: {}", p, e)
                        }
                    }
                }
                let bridged = crate::splitter::LinkedReadResolver::new(k)
                    .bridge_with_linked_reads(contigs, &all_linked);
                println!(
                    "  Contigs after linked-read repeat resolution: {} (Elapsed: {:.3}s)",
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
        crate::modes::RnaEngine::default().process_transcripts(contigs)
    } else {
        contigs
    };

    println!("─── [Stage 7] Scaffolding Across Unresolved Gaps ───");
    let scaffolds = if !config.skip_repeat_resolution && packed.is_paired() {
        let paired_info =
            crate::paired_info::PairedInfoIndex::build_from_packed(k, &contigs, packed);
        crate::scaffold::Scaffolder::default().build_scaffolds_with_reads(
            &contigs,
            &paired_info,
            Some(packed),
        )
    } else {
        contigs.clone()
    };
    println!(
        "  Scaffolds generated: {} (Elapsed: {:.3}s)",
        scaffolds.len(),
        start_time.elapsed().as_secs_f64()
    );

    let (contigs, scaffolds) = if config.polish && !config.skip_repeat_resolution {
        println!(
            "─── [Stage 8] Consensus Base Polishing{} ───",
            if config.careful {
                " (Careful Mode: Active)"
            } else {
                ""
            }
        );
        let polisher = crate::polisher::Polisher {
            k: 21,
            min_coverage_support: if config.careful { 8 } else { 5 },
            careful: config.careful,
        };
        let (polished_c, fixes_c) = polisher.polish_contigs_packed(contigs, packed);
        let (polished_s, fixes_s) = if packed.is_paired() {
            polisher.polish_contigs_packed(scaffolds, packed)
        } else {
            (polished_c.clone(), 0)
        };
        println!(
            "  Polished {} base discrepancies across contigs and scaffolds (Elapsed: {:.3}s)",
            fixes_c + fixes_s,
            start_time.elapsed().as_secs_f64()
        );
        (polished_c, polished_s)
    } else {
        (contigs, scaffolds)
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

    let mut contigs = if config.is_meta {
        let meta_contigs = crate::modes::apply_meta_filter(contigs, config.min_contig_len);
        println!(
            "  [Meta Mode] Retained {} contigs across uneven coverage depths",
            meta_contigs.len()
        );
        meta_contigs
    } else {
        contigs
    };

    let mut scaffolds = if config.is_meta {
        crate::modes::apply_meta_filter(scaffolds, config.min_contig_len)
    } else {
        scaffolds
    };

    let mut plasmids = plasmids;
    for (idx, u) in contigs.iter_mut().enumerate() {
        u.id = idx;
    }
    for (idx, u) in scaffolds.iter_mut().enumerate() {
        u.id = idx;
    }
    for (idx, u) in plasmids.iter_mut().enumerate() {
        u.id = idx;
    }

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
    let packed = if config.interleaved {
        crate::packed_reads::PackedReads::from_interleaved_files(input_files)?
    } else {
        crate::packed_reads::PackedReads::from_files(input_files)?
    };
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
    if total_length == 0 {
        return AssemblyStats::default();
    }
    let max_contig_length = lengths[0];

    let half_total = total_length.div_ceil(2);
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
    let mut atgc_count = 0usize;
    for c in contigs {
        for &b in &c.sequence {
            match b {
                b'G' | b'g' | b'C' | b'c' => {
                    gc_count += 1;
                    atgc_count += 1;
                }
                b'A' | b'a' | b'T' | b't' => {
                    atgc_count += 1;
                }
                _ => {}
            }
        }
    }
    let gc_content = if atgc_count > 0 {
        (gc_count as f64 / atgc_count as f64) * 100.0
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
        if c.sequence.is_empty() {
            writer.write_all(b"\n")?;
        } else {
            for chunk in c.sequence.chunks(80) {
                writer.write_all(chunk)?;
                writer.write_all(b"\n")?;
            }
        }
    }
    writer.flush()?;
    Ok(())
}

/// Determines the minimum k-mer coverage threshold dynamically using robust noise valley detection.
fn determine_solid_kmer_cutoff(global_counts: &HashMap<Kmer256, u32>, min_coverage: f64) -> u32 {
    if min_coverage > 5.0 {
        (min_coverage * 0.35).max(2.0) as u32
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
                let cutoff = if found_valley { valley } else { 2u32 };
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
    }
}
