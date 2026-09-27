use clap::Parser;
use spades_rs::fastq::PhredEncoding;
use spades_rs::qc::{process_paired_reads, process_single_reads, QcConfig};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(name = "spades-qc")]
#[command(about = "Stage 2: High-throughput FASTQ Read QC, Adapter Trimming & Preprocessing")]
#[command(
    after_help = "Examples:\n  spades-qc -1 r1.fq.gz -2 r2.fq.gz -o clean_1.fq.gz -O clean_2.fq.gz\n  spades-qc -s single.fq -o clean_single.fq --min-quality 20 --json qc.json"
)]
struct Cli {
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

    /// Minimum trailing poly-G length to trim (two-color NextSeq/NovaSeq artifact; 0 to disable)
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

    /// Number of worker threads (defaults to system logical threads)
    #[arg(short = 't', long = "threads")]
    threads: Option<usize>,

    /// Save machine-readable QC report to JSON file
    #[arg(long = "json", value_name = "FILE")]
    json: Option<PathBuf>,
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    if let Some(t) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()?;
    }

    let mut cfg = QcConfig {
        min_quality: cli.min_quality,
        window_size: cli.window_size,
        min_length: cli.min_len,
        max_ns: cli.max_ns,
        trim_poly_g: cli.trim_poly_g,
        min_adapter_overlap: cli.min_adapter_overlap,
        threads: cli.threads,
        ..Default::default()
    };

    if cli.phred64 {
        cfg.phred_override = Some(PhredEncoding::Phred64);
    } else if cli.phred33 {
        cfg.phred_override = Some(PhredEncoding::Phred33);
    }

    if let Some(user_adapters) = cli.adapters {
        for ad in user_adapters {
            let bytes = ad.into_bytes();
            cfg.adapters_r1.push(bytes.clone());
            cfg.adapters_r2.push(bytes);
        }
    }

    let report = match (cli.pe1_1, cli.pe1_2, cli.single) {
        (Some(r1), Some(r2), None) => {
            let o1 = cli
                .out1
                .unwrap_or_else(|| PathBuf::from("clean_1.fastq.gz"));
            let o2 = cli
                .out2
                .unwrap_or_else(|| PathBuf::from("clean_2.fastq.gz"));
            process_paired_reads(&r1, &r2, &o1, &o2, cli.unpaired.as_deref(), &cfg)?
        }
        (None, None, Some(s)) => {
            let o = cli
                .out1
                .unwrap_or_else(|| PathBuf::from("clean_single.fastq.gz"));
            process_single_reads(&s, &o, &cfg)?
        }
        (Some(_), None, _) | (None, Some(_), _) => {
            anyhow::bail!("Paired-end processing requires both -1/--pe1-1 and -2/--pe1-2.");
        }
        (None, None, None) => {
            anyhow::bail!("No input reads provided. Specify -1 and -2 for paired reads, or -s for single reads.");
        }
        (Some(_), Some(_), Some(_)) => {
            anyhow::bail!("Cannot mix paired reads (-1/-2) and single reads (-s) in one run.");
        }
    };

    report.print_summary();

    if let Some(json_path) = cli.json {
        std::fs::write(&json_path, report.to_json())?;
        println!("  QC JSON report saved to: {:?}", json_path);
    }

    Ok(())
}
