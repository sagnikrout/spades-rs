use clap::{Parser, Subcommand};
use spades_rs::assemble::{run_assembly, write_contigs_fasta, AssemblerConfig};
use spades_rs::multik::{run_multik_assembly, MultiKConfig};
use spades_rs::scaffold::write_scaffolds_fasta;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(name = "spades-rs")]
#[command(about = "Ultra-Fast, Low-Memory De Novo Genome Assembler in Pure Rust")]
#[command(after_help = "Citations:\n  SPAdes: Bankevich et al. (2012) J Comput Biol 19(5):455-477\n  Protocol: Prjibelski et al. (2020) Curr Protoc Bioinformatics 70(1):e102\n  See README.md for full citations & BibTeX entries.")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Assemble reads from FASTQ / FASTA files into contigs and scaffolds
    Assemble {
        /// Input FASTQ / FASTA files (can specify multiple or gzipped .fq.gz)
        #[arg(short, long, required = true, num_args = 1..)]
        inputs: Vec<PathBuf>,

        /// Output contigs FASTA file
        #[arg(short, long, default_value = "contigs.fasta")]
        output: PathBuf,

        /// K-mer size (must be odd and <= 127)
        #[arg(short, long, default_value_t = 31)]
        k: usize,

        /// Minimum average coverage threshold for contigs
        #[arg(short, long, default_value_t = 5.0)]
        coverage: f64,

        /// Minimum contig length to output (bp)
        #[arg(short, long, default_value_t = 200)]
        min_len: usize,

        /// Number of worker threads (defaults to system logical threads)
        #[arg(short, long)]
        threads: Option<usize>,

        /// Run BayesHammer read error correction before assembling
        #[arg(long)]
        error_correct: bool,

        /// Enable metaSPAdes metagenomics mode for uneven coverage datasets
        #[arg(long)]
        meta: bool,

        /// Enable plasmidSPAdes mode to extract plasmids into a separate FASTA
        #[arg(long)]
        plasmid: bool,

        /// Enable rnaSPAdes mode for transcriptome and alternative splicing isoforms
        #[arg(long)]
        rna: bool,

        /// Enable scSPAdes mode for single-cell MDA coverage normalization
        #[arg(long)]
        sc: bool,

        /// Disable consensus base polishing
        #[arg(long)]
        no_polish: bool,

        /// Oxford Nanopore reads for hybrid long-read repeat bridging
        #[arg(long)]
        nanopore: Option<PathBuf>,

        /// PacBio HiFi reads for hybrid long-read repeat bridging
        #[arg(long)]
        pacbio: Option<PathBuf>,

        /// Optional list of k-mers for multi-k iterative assembly (e.g. 21,33,55)
        #[arg(long, value_delimiter = ',', num_args = 1..)]
        multik: Option<Vec<usize>>,

        /// Maximum physical memory budget in GB (defaults to Available System RAM - 20%)
        #[arg(long)]
        max_memory: Option<f64>,
    },

    /// Run automatic benchmark on SPAdes reference test dataset
    Benchmark {
        /// Path to SPAdes test dataset directory (defaults to data/)
        #[arg(long)]
        test_dir: Option<PathBuf>,

        /// K-mer size
        #[arg(short, long, default_value_t = 31)]
        k: usize,
    },
}

fn main() -> anyhow::Result<()> {
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    unsafe {
        extern "C" {
            fn mallopt(param: i32, value: i32) -> i32;
        }
        const M_ARENA_MAX: i32 = -8;
        mallopt(M_ARENA_MAX, 2);
    }

    let cli = Cli::parse();

    match cli.command {
        Commands::Assemble {
            inputs,
            output,
            k,
            coverage,
            min_len,
            threads,
            error_correct,
            meta,
            plasmid,
            rna,
            sc,
            no_polish,
            nanopore,
            pacbio,
            multik,
            max_memory,
        } => {
            if let Some(t) = threads {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(t)
                    .build_global()?;
            }

            let memory_limits = spades_rs::memory::MemoryLimits::determine(max_memory);

            let mut long_reads = Vec::new();
            if let Some(np) = nanopore {
                long_reads.push(np);
            }
            if let Some(pb) = pacbio {
                long_reads.push(pb);
            }
            let long_reads_opt = if !long_reads.is_empty() {
                Some(long_reads)
            } else {
                None
            };

            println!("===========================================================");
            println!("      SPADES-RS: ULTRA-FAST DE NOVO GENOME ASSEMBLER       ");
            println!("===========================================================");
            println!(
                "  Hardware Concurrency: {} threads active",
                rayon::current_num_threads()
            );
            println!(
                "  Memory Governor:      {:.2} GB budget ({}) [System: {:.2} GB avail / {:.2} GB total]",
                memory_limits.budget_gb(),
                if memory_limits.is_user_specified {
                    "User Specified: --max-memory"
                } else {
                    "Auto: 80% of available RAM, 20% reserved for OS"
                },
                memory_limits.available_gb(),
                memory_limits.total_gb(),
            );
            println!("  K-mer size: {}", k);
            println!("  Min Coverage: {:.1}x", coverage);
            println!("  Min Contig Length: {} bp", min_len);
            println!(
                "  Error Correction: {}",
                if error_correct {
                    "ENABLED (BayesHammer)"
                } else {
                    "DISABLED"
                }
            );
            println!("  Meta Mode: {}", if meta { "ENABLED" } else { "DISABLED" });
            println!(
                "  Plasmid Mode: {}",
                if plasmid { "ENABLED" } else { "DISABLED" }
            );
            println!(
                "  RNA Mode: {}",
                if rna {
                    "ENABLED (Isoform Preserver)"
                } else {
                    "DISABLED"
                }
            );
            println!(
                "  Single-Cell MDA: {}",
                if sc {
                    "ENABLED (MDA Normalizer)"
                } else {
                    "DISABLED"
                }
            );
            println!(
                "  Consensus Polishing: {}",
                if !no_polish { "ENABLED" } else { "DISABLED" }
            );
            if let Some(ref kms) = multik {
                println!("  Multi-K Iteration: {:?}", kms);
            }
            println!("  Output path: {:?}", output);
            println!("-----------------------------------------------------------");

            let result = if let Some(kms) = multik {
                let mk_config = MultiKConfig {
                    kmers: kms,
                    min_coverage: coverage,
                    min_contig_len: min_len,
                    bloom_bits: memory_limits.optimal_bloom_bits(),
                    error_correct,
                    is_meta: meta,
                    is_plasmid: plasmid,
                    is_rna: rna,
                    is_sc: sc,
                    polish: !no_polish,
                    long_reads: long_reads_opt,
                    memory_limits: Some(memory_limits),
                };
                run_multik_assembly(&inputs, &mk_config)?
            } else {
                let config = AssemblerConfig {
                    k,
                    min_coverage: coverage,
                    min_contig_len: min_len,
                    bloom_bits: memory_limits.optimal_bloom_bits(),
                    error_correct,
                    is_meta: meta,
                    is_plasmid: plasmid,
                    is_rna: rna,
                    is_sc: sc,
                    polish: !no_polish,
                    long_reads: long_reads_opt,
                    prior_contigs: None,
                    skip_repeat_resolution: false,
                    memory_limits: Some(memory_limits),
                };
                run_assembly(&inputs, &config)?
            };

            println!("-----------------------------------------------------------");
            println!("               ASSEMBLY QUALITY SUMMARY                    ");
            println!("-----------------------------------------------------------");
            println!("  Total Contigs:        {}", result.stats.total_contigs);
            println!("  Total Scaffolds:      {}", result.scaffolds.len());
            println!("  Total Assembled bp:   {} bp", result.stats.total_length);
            println!(
                "  Max Contig Length:    {} bp",
                result.stats.max_contig_length
            );
            println!("  N50:                  {} bp", result.stats.n50);
            println!("  L50:                  {}", result.stats.l50);
            println!("  GC Content:           {:.2}%", result.stats.gc_content);
            println!("  Total Wall-Clock:     {:.4} seconds", result.elapsed_secs);
            println!("-----------------------------------------------------------");

            let (contig_path, scaffold_path, gfa_path, plasmid_path) = if output.is_dir()
                || output.extension().is_none()
            {
                std::fs::create_dir_all(&output)?;
                (
                    output.join("contigs.fasta"),
                    output.join("scaffolds.fasta"),
                    output.join("assembly_graph.gfa"),
                    output.join("plasmids.fasta"),
                )
            } else {
                if let Some(parent) = output.parent() {
                    let _ = std::fs::create_dir_all(parent);
                }
                (
                    output.clone(),
                    output.with_file_name("scaffolds.fasta"),
                    output.with_extension("gfa"),
                    output.with_file_name("plasmids.fasta"),
                )
            };

            write_contigs_fasta(&result.contigs, &contig_path)?;
            println!("  Contigs successfully exported to: {:?}", contig_path);

            write_scaffolds_fasta(&result.scaffolds, &scaffold_path)?;
            println!("  Scaffolds successfully exported to: {:?}", scaffold_path);

            spades_rs::gfa::write_graph_gfa(k, &result.contigs, &gfa_path)?;
            println!("  Assembly graph (GFA v1.1) exported to: {:?}", gfa_path);

            if plasmid && !result.plasmids.is_empty() {
                write_contigs_fasta(&result.plasmids, &plasmid_path)?;
                println!("  Plasmids successfully exported to: {:?}", plasmid_path);
            }
        }

        Commands::Benchmark { test_dir, k } => {
            let default_path = PathBuf::from("data");
            let dir = test_dir.unwrap_or(default_path);

            let r1 = dir.join("ecoli_1K_1.fq.gz");
            let r2 = dir.join("ecoli_1K_2.fq.gz");

            if !r1.exists() || !r2.exists() {
                anyhow::bail!("Test files not found in {:?}.", dir);
            }

            println!("===========================================================");
            println!("   RUNNING BENCHMARK ON E. COLI 1K TEST DATASET            ");
            println!("===========================================================");
            println!("  Read 1: {:?}", r1);
            println!("  Read 2: {:?}", r2);
            println!("  Active Threads: {}", rayon::current_num_threads());

            let config = AssemblerConfig {
                k,
                min_coverage: 5.0,
                min_contig_len: 200,
                bloom_bits: 16 * 1024 * 1024,
                error_correct: true,
                is_meta: false,
                is_plasmid: false,
                is_rna: false,
                is_sc: false,
                polish: true,
                long_reads: None,
                prior_contigs: None,
                skip_repeat_resolution: false,
                memory_limits: None,
            };

            let t0 = Instant::now();
            let result = run_assembly(&[r1, r2], &config)?;
            let total_wall = t0.elapsed().as_secs_f64();

            println!("-----------------------------------------------------------");
            println!("               BENCHMARK RESULTS                           ");
            println!("-----------------------------------------------------------");
            println!("  Total Assembled Contigs: {}", result.stats.total_contigs);
            println!(
                "  Max Contig Length:       {} bp",
                result.stats.max_contig_length
            );
            println!("  Total Wall-Clock Time:   {:.4} seconds", total_wall);
            println!("===========================================================");
        }
    }

    Ok(())
}
