use clap::{Parser, Subcommand};
use spades_rs::assemble::{run_assembly, write_contigs_fasta, AssemblerConfig};
use spades_rs::multik::{run_multik_assembly, MultiKConfig};
use spades_rs::scaffold::write_scaffolds_fasta;
use std::path::PathBuf;
use std::time::Instant;

/// Validates the k-mer size: must be odd, between 11 and 127 (inclusive).
fn validate_k(s: &str) -> Result<usize, String> {
    let k: usize = s
        .parse()
        .map_err(|_| format!("'{}' is not a valid integer", s))?;
    if k < 11 {
        return Err(format!("k={} is too small; minimum is 11", k));
    }
    if k > 127 {
        return Err(format!("k={} exceeds maximum of 127 (Kmer256 limit)", k));
    }
    if k.is_multiple_of(2) {
        return Err(format!(
            "k={} is even; SPAdes requires an odd k-mer size",
            k
        ));
    }
    Ok(k)
}

#[derive(Parser, Debug)]
#[command(name = "spades-rs")]
#[command(about = "A Rust-based de novo genome assembler designed for low-memory environments")]
#[command(
    after_help = "Citations:\n  SPAdes: Bankevich et al. (2012) J Comput Biol 19(5):455-477\n  Protocol: Prjibelski et al. (2020) Curr Protoc Bioinformatics 70(1):e102\n  SpLitteR: Tolstoganov et al. (2024) PeerJ 12:e18050\n  See README.md for full citations & BibTeX entries."
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
#[allow(clippy::large_enum_variant)]
enum Commands {
    /// Assemble reads from FASTQ / FASTA files into contigs and scaffolds
    Assemble {
        /// Input FASTQ / FASTA files (can specify multiple or gzipped .fq.gz)
        #[arg(short, long, num_args = 0..)]
        inputs: Vec<PathBuf>,

        /// File with forward paired-end reads (SPAdes -1 option)
        #[arg(short = '1', long = "pe1-1", value_name = "FILE")]
        pe1_1: Option<PathBuf>,

        /// File with reverse paired-end reads (SPAdes -2 option)
        #[arg(short = '2', long = "pe1-2", value_name = "FILE")]
        pe1_2: Option<PathBuf>,

        /// File with unpaired / single reads (SPAdes -s option)
        #[arg(short = 's', long = "pe1-s", value_name = "FILE")]
        pe1_s: Option<PathBuf>,

        /// File with interleaved paired-end reads (SPAdes --12 option)
        #[arg(long = "12", value_name = "FILE")]
        pe1_12: Option<PathBuf>,

        /// Output contigs FASTA file
        #[arg(short, long, default_value = "contigs.fasta")]
        output: PathBuf,

        /// K-mer size (must be odd, between 11 and 127)
        #[arg(short, long, default_value_t = 31, value_parser = validate_k)]
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

        /// Run careful mode (tries to reduce number of mismatches and short indels)
        #[arg(long)]
        careful: bool,

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

        /// Path to trusted contigs / prior high-confidence backbones
        #[arg(long, value_delimiter = ',', num_args = 1..)]
        trusted_contigs: Option<Vec<PathBuf>>,

        /// File with barcoded linked reads for SpLitteR repeat resolution (TELL-Seq or 10x)
        #[arg(long = "splitter", aliases = ["linked-reads"], value_delimiter = ',', num_args = 1..)]
        splitter: Option<Vec<PathBuf>>,

        /// Optional list of k-mers for multi-k iterative assembly (e.g. 21,33,55)
        #[arg(long, value_delimiter = ',', num_args = 1..)]
        multik: Option<Vec<usize>>,

        /// Maximum physical memory budget in GB (defaults to Available System RAM - 20%)
        #[arg(long)]
        max_memory: Option<f64>,
    },

    /// Stage 2: Quality control, adapter clipping, and read preprocessing
    #[command(name = "qc", aliases = ["preprocess"])]
    Qc {
        /// Forward paired-end reads file (FASTQ / FASTQ.GZ)
        #[arg(short = '1', long = "pe1-1", value_name = "FILE")]
        pe1_1: Option<PathBuf>,

        /// Reverse paired-end reads file (FASTQ / FASTQ.GZ)
        #[arg(short = '2', long = "pe1-2", value_name = "FILE")]
        pe1_2: Option<PathBuf>,

        /// Single-end reads file (FASTQ / FASTQ.GZ)
        #[arg(short = 's', long = "single", value_name = "FILE")]
        single: Option<PathBuf>,

        /// Clean output file for forward reads (or single reads)
        #[arg(short = 'o', long = "out1", value_name = "FILE")]
        out1: Option<PathBuf>,

        /// Clean output file for reverse reads
        #[arg(short = 'O', long = "out2", value_name = "FILE")]
        out2: Option<PathBuf>,

        /// Output file for surviving unpaired / singleton reads from filtered pairs
        #[arg(long = "unpaired", value_name = "FILE")]
        unpaired: Option<PathBuf>,

        /// Minimum Phred quality threshold for sliding-window trimming
        #[arg(long, default_value_t = 15)]
        min_quality: u8,

        /// Sliding window size for quality trimming
        #[arg(long, default_value_t = 4)]
        window_size: usize,

        /// Minimum read length to keep after trimming
        #[arg(short = 'l', long = "min-len", default_value_t = 30)]
        min_len: usize,

        /// Maximum allowed uncalled 'N' bases per read
        #[arg(long, default_value_t = 3)]
        max_ns: usize,

        /// Minimum trailing poly-G length to trim (0 to disable)
        #[arg(long, default_value_t = 5)]
        trim_poly_g: usize,

        /// Minimum adapter overlap length for 3' adapter clipping
        #[arg(long, default_value_t = 10)]
        min_adapter_overlap: usize,

        /// Additional custom adapter sequence(s) to trim
        #[arg(long = "adapter", num_args = 1..)]
        adapters: Option<Vec<String>>,

        /// Force Phred+64 encoding format
        #[arg(long)]
        phred64: bool,

        /// Force Phred+33 encoding format
        #[arg(long)]
        phred33: bool,

        /// Number of worker threads
        #[arg(short = 't', long = "threads")]
        threads: Option<usize>,

        /// Save machine-readable QC report to JSON file
        #[arg(long = "json", value_name = "FILE")]
        json: Option<PathBuf>,
    },

    /// Stage 5: QUAST-equivalent assembly evaluation
    #[command(name = "eval")]
    Eval {
        /// Assembly contigs or scaffolds FASTA file (plain text or .gz)
        #[arg(value_name = "CONTIGS")]
        contigs: PathBuf,

        /// Reference genome FASTA file for genome fraction and mismatch scoring
        #[arg(short = 'r', long = "reference", value_name = "FILE")]
        reference: Option<PathBuf>,

        /// Minimum contig length cutoff (bp)
        #[arg(short = 'l', long = "min-len", default_value_t = 200)]
        min_len: usize,

        /// Save machine-readable evaluation report to JSON file
        #[arg(long = "json", value_name = "FILE")]
        json: Option<PathBuf>,
    },

    /// Stage 5: Standalone consensus base polishing for draft contigs
    #[command(name = "polish")]
    Polish {
        /// Assembly contigs or scaffolds FASTA file to polish
        #[arg(value_name = "CONTIGS")]
        contigs: PathBuf,

        /// Read files for polishing (comma-separated or multiple flags)
        #[arg(short = 'i', long = "reads", value_delimiter = ',', num_args = 1..)]
        reads: Vec<PathBuf>,

        /// Output path for polished contigs FASTA
        #[arg(short = 'o', long = "output", default_value = "polished_contigs.fasta")]
        output: PathBuf,

        /// K-mer size for polishing alignment anchors
        #[arg(short = 'k', default_value_t = 21)]
        k: usize,

        /// Minimum coverage threshold required to alter a consensus base
        #[arg(long = "min-coverage", default_value_t = 5)]
        min_coverage: u32,

        /// Enable careful polishing mode
        #[arg(long)]
        careful: bool,
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

    let mut args: Vec<String> = std::env::args().collect();
    if args.len() > 1 {
        let first = &args[1];
        if first != "assemble"
            && first != "benchmark"
            && first != "qc"
            && first != "preprocess"
            && first != "eval"
            && first != "polish"
            && first != "help"
            && first != "-h"
            && first != "--help"
            && first != "-V"
            && first != "--version"
        {
            args.insert(1, "assemble".to_string());
        }
    }
    let cli = Cli::parse_from(args);

    match cli.command {
        Commands::Assemble {
            inputs,
            pe1_1,
            pe1_2,
            pe1_s,
            pe1_12,
            output,
            k,
            coverage,
            min_len,
            threads,
            error_correct,
            careful,
            meta,
            plasmid,
            rna,
            sc,
            no_polish,
            nanopore,
            pacbio,
            trusted_contigs,
            splitter,
            multik,
            max_memory,
        } => {
            if let Some(t) = threads {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(t)
                    .build_global()?;
            }

            let mut all_inputs = inputs;
            if let Some(r1) = pe1_1 {
                all_inputs.insert(0, r1);
                if let Some(r2) = pe1_2 {
                    all_inputs.insert(1, r2);
                }
            } else if let Some(r2) = pe1_2 {
                all_inputs.push(r2);
            }
            if let Some(s) = pe1_s {
                all_inputs.push(s);
            }
            if let Some(interleaved) = pe1_12 {
                all_inputs.push(interleaved);
            }

            if all_inputs.is_empty() {
                anyhow::bail!(
                    "No input reads provided. Specify reads via -1 / -2, -s, --12, or -i / --inputs."
                );
            }

            let effective_error_correct = error_correct || careful;

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
            println!("            SPADES-RS: DE NOVO GENOME ASSEMBLER            ");
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
                if effective_error_correct {
                    "ENABLED (BayesHammer)"
                } else {
                    "DISABLED"
                }
            );
            println!(
                "  Careful Mode:     {}",
                if careful {
                    "ENABLED (Mismatch Corrections Active)"
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
            if let Some(ref trusted) = trusted_contigs {
                println!("  Trusted Contigs:      {} file(s)", trusted.len());
            }
            if let Some(ref sp) = splitter {
                println!("  SpLitteR Barcoded:    {} file(s)", sp.len());
            }
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
                    error_correct: effective_error_correct,
                    is_meta: meta,
                    is_plasmid: plasmid,
                    is_rna: rna,
                    is_sc: sc,
                    polish: !no_polish,
                    careful,
                    long_reads: long_reads_opt,
                    linked_reads: splitter,
                    trusted_contigs,
                    memory_limits: Some(memory_limits),
                };
                run_multik_assembly(&all_inputs, &mk_config)?
            } else {
                let config = AssemblerConfig {
                    k,
                    min_coverage: coverage,
                    min_contig_len: min_len,
                    bloom_bits: memory_limits.optimal_bloom_bits(),
                    error_correct: effective_error_correct,
                    is_meta: meta,
                    is_plasmid: plasmid,
                    is_rna: rna,
                    is_sc: sc,
                    polish: !no_polish,
                    careful,
                    long_reads: long_reads_opt,
                    linked_reads: splitter,
                    trusted_contigs,
                    prior_contigs: None,
                    skip_repeat_resolution: false,
                    memory_limits: Some(memory_limits),
                };
                run_assembly(&all_inputs, &config)?
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

            let (contig_path, scaffold_path, gfa_path, plasmid_path) =
                if output.is_dir() || output.extension().is_none() {
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

            // Stage 5: Automatic assembly evaluation summary
            if let Ok(eval_metrics) =
                spades_rs::eval::evaluate_assembly(&contig_path, None, min_len)
            {
                println!();
                eval_metrics.print_summary();
            }
        }

        Commands::Qc {
            pe1_1,
            pe1_2,
            single,
            out1,
            out2,
            unpaired,
            min_quality,
            window_size,
            min_len,
            max_ns,
            trim_poly_g,
            min_adapter_overlap,
            adapters,
            phred64,
            phred33,
            threads,
            json,
        } => {
            if let Some(t) = threads {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(t)
                    .build_global()?;
            }

            let mut cfg = spades_rs::qc::QcConfig {
                min_quality,
                window_size,
                min_length: min_len,
                max_ns,
                trim_poly_g,
                min_adapter_overlap,
                threads,
                ..Default::default()
            };

            if phred64 {
                cfg.phred_override = Some(spades_rs::fastq::PhredEncoding::Phred64);
            } else if phred33 {
                cfg.phred_override = Some(spades_rs::fastq::PhredEncoding::Phred33);
            }

            if let Some(user_adapters) = adapters {
                for ad in user_adapters {
                    let bytes = ad.into_bytes();
                    cfg.adapters_r1.push(bytes.clone());
                    cfg.adapters_r2.push(bytes);
                }
            }

            let report = match (pe1_1, pe1_2, single) {
                (Some(r1), Some(r2), None) => {
                    let o1 = out1.unwrap_or_else(|| PathBuf::from("clean_1.fastq.gz"));
                    let o2 = out2.unwrap_or_else(|| PathBuf::from("clean_2.fastq.gz"));
                    spades_rs::qc::process_paired_reads(
                        &r1,
                        &r2,
                        &o1,
                        &o2,
                        unpaired.as_deref(),
                        &cfg,
                    )?
                }
                (None, None, Some(s)) => {
                    let o = out1.unwrap_or_else(|| PathBuf::from("clean_single.fastq.gz"));
                    spades_rs::qc::process_single_reads(&s, &o, &cfg)?
                }
                (Some(_), None, _) | (None, Some(_), _) => {
                    anyhow::bail!("Paired-end processing requires both -1/--pe1-1 and -2/--pe1-2.");
                }
                (None, None, None) => {
                    anyhow::bail!("No input reads provided. Specify -1 and -2 for paired reads, or -s for single reads.");
                }
                (Some(_), Some(_), Some(_)) => {
                    anyhow::bail!(
                        "Cannot mix paired reads (-1/-2) and single reads (-s) in one run."
                    );
                }
            };

            report.print_summary();

            if let Some(json_path) = json {
                std::fs::write(&json_path, report.to_json())?;
                println!("  QC JSON report saved to: {:?}", json_path);
            }
        }

        Commands::Eval {
            contigs,
            reference,
            min_len,
            json,
        } => {
            let metrics =
                spades_rs::eval::evaluate_assembly(&contigs, reference.as_deref(), min_len)?;
            metrics.print_summary();
            if let Some(json_path) = json {
                std::fs::write(&json_path, metrics.to_json())?;
                println!("  Evaluation JSON report saved to: {:?}", json_path);
            }
        }

        Commands::Polish {
            contigs,
            reads,
            output,
            k,
            min_coverage,
            careful,
        } => {
            println!("Running Stage 5 standalone consensus base polishing...");
            let corrected = spades_rs::eval::polish_assembly_file(
                &contigs,
                &reads,
                &output,
                k,
                min_coverage,
                careful,
            )?;
            println!(
                "Consensus polishing completed: {} bases modified. Output saved to: {:?}",
                corrected, output
            );
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
                careful: false,
                long_reads: None,
                linked_reads: None,
                trusted_contigs: None,
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
