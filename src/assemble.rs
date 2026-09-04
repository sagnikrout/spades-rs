//! End-to-End De Novo Assembly Pipeline & Assembly QC.

use crate::bloom::TwoTierFilter;
use crate::dna::{canonical_kmer_u64, string_to_kmer};
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
}

impl Default for AssemblerConfig {
    fn default() -> Self {
        Self {
            k: 31,
            min_coverage: 5.0,
            min_contig_len: 200,
            bloom_bits: 64 * 1024 * 1024, // 64M bits = 8 MB RAM
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

pub struct AssemblyResult {
    pub contigs: Vec<crate::graph::Unitig>,
    pub scaffolds: Vec<crate::graph::Unitig>,
    pub plasmids: Vec<crate::graph::Unitig>,
    pub stats: AssemblyStats,
    pub elapsed_secs: f64,
}

/// Assembles input FASTQ/FASTA files into contigs.
pub fn run_assembly<P: AsRef<Path> + Sync>(
    input_files: &[P],
    config: &AssemblerConfig,
) -> Result<AssemblyResult> {
    let start_time = Instant::now();
    let k = config.k;

    println!("─── [Stage 1] Ingesting reads from {} input file(s) ───", input_files.len());
    let mut all_reads: Vec<Vec<u8>> = Vec::new();
    for f in input_files {
        let reads = parse_reads_from_file(f)?;
        println!("  • Loaded {} reads from {:?}", reads.len(), f.as_ref().file_name().unwrap_or_default());
        all_reads.extend(reads);
    }
    println!("  Total reads loaded: {} (Elapsed: {:.3}s)", all_reads.len(), start_time.elapsed().as_secs_f64());

    println!("─── [Stage 2] Streaming reads into Two-Tier Bloom Filter (Memory Shield) ───");
    let filter = TwoTierFilter::new(config.bloom_bits);
    println!("  Bloom filter memory: {:.2} MB", filter.memory_usage_bytes() as f64 / 1_048_576.0);

    all_reads.par_iter().for_each(|seq| {
        if seq.len() < k {
            return;
        }
        for i in 0..=(seq.len() - k) {
            if let Some(kmer) = string_to_kmer(&seq[i..i + k], k) {
                let (can, _) = canonical_kmer_u64(kmer, k);
                filter.insert(can);
            }
        }
    });
    println!("  Bloom filter populated. (Elapsed: {:.3}s)", start_time.elapsed().as_secs_f64());

    let all_reads = if config.error_correct {
        println!("─── [Stage 2.5] BayesHammer Read Error Correction ───");
        let corrector = crate::hammer::ErrorCorrector::new(k);
        let (corrected, count) = corrector.correct_reads(all_reads, &filter);
        println!("  Corrected {} reads via solid consensus (Elapsed: {:.3}s)", count, start_time.elapsed().as_secs_f64());
        corrected
    } else {
        all_reads
    };

    println!("─── [Stage 3] Building Solid K-mer Index ───");
    // Thread-local accumulation of solid k-mers to avoid mutex locks
    let solid_maps: Vec<HashMap<u64, u32>> = all_reads
        .par_chunks(2000)
        .map(|chunk| {
            let mut local_counts: HashMap<u64, u32> = HashMap::with_capacity(4096);

            for seq in chunk {
                if seq.len() < k {
                    continue;
                }
                for i in 0..=(seq.len() - k) {
                    if let Some(kmer) = string_to_kmer(&seq[i..i + k], k) {
                        let (can, _) = canonical_kmer_u64(kmer, k);
                        if filter.is_solid(can) {
                            *local_counts.entry(can).or_insert(0) += 1;
                        }
                    }
                }
            }
            local_counts
        })
        .collect();

    // Merge thread-local maps
    let mut global_counts: HashMap<u64, u32> = HashMap::new();
    for l_counts in solid_maps {
        for (kmer, cnt) in l_counts {
            *global_counts.entry(kmer).or_insert(0) += cnt;
        }
    }

    // Filter out k-mers with coverage below minimum threshold (e.g. noise filter)
    let min_kmer_cov = (config.min_coverage * 0.2).max(2.0) as u32;
    let solid_kmers: hashbrown::HashSet<u64> = global_counts
        .iter()
        .filter(|(_, &cov)| cov >= min_kmer_cov)
        .map(|(&kmer, _)| kmer)
        .collect();

    println!("  Total solid k-mers retained: {} (Elapsed: {:.3}s)", solid_kmers.len(), start_time.elapsed().as_secs_f64());

    println!("─── [Stage 4] Compacting de Bruijn Graph into Unitigs ───");
    let cdbg = CompactedGraph::build(k, &solid_kmers, &global_counts);
    println!("  Raw unitigs constructed: {} (Elapsed: {:.3}s)", cdbg.unitigs.len(), start_time.elapsed().as_secs_f64());

    let raw_unitigs = if config.is_sc {
        println!("─── [Single-Cell Mode] Normalizing MDA Coverage Discrepancies ───");
        crate::single_cell::SingleCellNormalizer::default().normalize_coverage(cdbg.unitigs)
    } else {
        cdbg.unitigs
    };

    println!("─── [Stage 5] Simplification (Tip Clipping & Artifact Cleaning) ───");
    let simplifier = Simplifier::new(k, config.min_coverage, config.min_contig_len);
    let simplified_contigs = simplifier.simplify(raw_unitigs);
    println!("  Assembled contigs after simplification: {} (Elapsed: {:.3}s)", simplified_contigs.len(), start_time.elapsed().as_secs_f64());

    println!("─── [Stage 6] ExSPAnder Repeat Resolution (Paired-End Linkages) ───");
    let contigs = if input_files.len() >= 2 {
        let reads1 = parse_reads_from_file(&input_files[0])?;
        let reads2 = parse_reads_from_file(&input_files[1])?;
        let paired_info = crate::paired_info::PairedInfoIndex::build(k, &simplified_contigs, &reads1, &reads2);
        println!("  Paired library estimated insert size: {:.1} ± {:.1} bp", paired_info.mean_insert_size, paired_info.insert_size_stdev);
        let expander = crate::expander::ExSPAnder::default();
        let resolved = expander.resolve_repeats(k, simplified_contigs, &paired_info);
        println!("  Contigs after repeat resolution: {} (Elapsed: {:.3}s)", resolved.len(), start_time.elapsed().as_secs_f64());
        resolved
    } else {
        simplified_contigs
    };

    let contigs = if let Some(ref lr_paths) = config.long_reads {
        if !lr_paths.is_empty() {
            println!("─── [Stage 6.5] Spaligner Hybrid Long-Read Bridging ───");
            let mut lr_seqs = Vec::new();
            for p in lr_paths {
                if let Ok(reads) = parse_reads_from_file(p) {
                    lr_seqs.extend(reads);
                }
            }
            let bridged = crate::spaligner::LongReadResolver::default().bridge_with_long_reads(contigs, &lr_seqs);
            println!("  Contigs after long-read bridging: {} (Elapsed: {:.3}s)", bridged.len(), start_time.elapsed().as_secs_f64());
            bridged
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
    let scaffolds = if input_files.len() >= 2 {
        let reads1 = parse_reads_from_file(&input_files[0])?;
        let reads2 = parse_reads_from_file(&input_files[1])?;
        let paired_info = crate::paired_info::PairedInfoIndex::build(k, &contigs, &reads1, &reads2);
        crate::scaffold::Scaffolder::default().build_scaffolds(&contigs, &paired_info)
    } else {
        contigs.clone()
    };
    println!("  Scaffolds generated: {} (Elapsed: {:.3}s)", scaffolds.len(), start_time.elapsed().as_secs_f64());

    let contigs = if config.polish {
        println!("─── [Stage 8] Consensus Base Polishing ───");
        let (polished, fixes) = crate::polisher::Polisher::default().polish_contigs(contigs, &all_reads);
        println!("  Polished {} base discrepancies (Elapsed: {:.3}s)", fixes, start_time.elapsed().as_secs_f64());
        polished
    } else {
        contigs
    };

    let (contigs, plasmids) = if config.is_plasmid {
        let (chrom, plas) = crate::modes::PlasmidDetector::default().extract_plasmids(contigs);
        println!("  [Plasmid Mode] Separated {} circular/high-copy plasmids and {} chromosomal contigs", plas.len(), chrom.len());
        (chrom, plas)
    } else {
        (contigs, Vec::new())
    };

    let contigs = if config.is_meta {
        let meta_contigs = crate::modes::apply_meta_filter(contigs, config.min_contig_len);
        println!("  [Meta Mode] Retained {} contigs across uneven coverage depths", meta_contigs.len());
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
pub fn write_contigs_fasta<P: AsRef<Path>>(contigs: &[crate::graph::Unitig], out_path: P) -> Result<()> {
    let mut file = File::create(out_path)?;
    for (i, c) in contigs.iter().enumerate() {
        writeln!(
            file,
            ">contig_{}_len_{}_cov_{:.1}",
            i + 1,
            c.sequence.len(),
            c.mean_coverage
        )?;
        // Write 80 characters per line
        for chunk in c.sequence.chunks(80) {
            file.write_all(chunk)?;
            file.write_all(b"\n")?;
        }
    }
    Ok(())
}
