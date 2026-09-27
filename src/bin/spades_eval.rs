use clap::Parser;
use spades_rs::eval::{evaluate_assembly, polish_assembly_file};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(name = "spades-eval")]
#[command(about = "Stage 5: QUAST-equivalent Assembly Evaluation & Standalone Polishing")]
#[command(
    after_help = "Examples:\n  spades-eval contigs.fasta\n  spades-eval contigs.fasta -r reference.fasta --json report.json\n  spades-eval contigs.fasta --polish --reads r1.fq.gz,r2.fq.gz -o polished.fasta"
)]
struct Cli {
    /// Assembly contigs or scaffolds FASTA file (plain text or .gz)
    #[arg(value_name = "CONTIGS")]
    contigs: PathBuf,

    /// Reference genome FASTA file for genome fraction and mismatch scoring
    #[arg(short = 'r', long = "reference", value_name = "FILE")]
    reference: Option<PathBuf>,

    /// Minimum contig length cutoff (bp)
    #[arg(short = 'l', long = "min-len", default_value_t = 200)]
    min_len: usize,

    /// Enable read-based consensus base polishing before evaluation
    #[arg(long)]
    polish: bool,

    /// Read files for polishing (comma-separated or multiple flags)
    #[arg(long = "reads", value_delimiter = ',', num_args = 1..)]
    reads: Option<Vec<PathBuf>>,

    /// Output path for polished contigs FASTA (used when --polish is set)
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    output: Option<PathBuf>,

    /// K-mer size for polishing alignment anchors
    #[arg(short = 'k', default_value_t = 21)]
    k: usize,

    /// Minimum coverage threshold required to alter a consensus base
    #[arg(long = "min-coverage", default_value_t = 5)]
    min_coverage: u32,

    /// Enable careful polishing mode
    #[arg(long)]
    careful: bool,

    /// Save machine-readable evaluation report to JSON file
    #[arg(long = "json", value_name = "FILE")]
    json: Option<PathBuf>,
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    let target_eval_file = if cli.polish {
        let reads = cli.reads.ok_or_else(|| {
            anyhow::anyhow!("Polishing requires read files via --reads <FILE1,FILE2>.")
        })?;
        let out_contigs = cli
            .output
            .unwrap_or_else(|| PathBuf::from("polished_contigs.fasta"));

        println!("Running Stage 5 standalone consensus base polishing...");
        let corrected = polish_assembly_file(
            &cli.contigs,
            &reads,
            &out_contigs,
            cli.k,
            cli.min_coverage,
            cli.careful,
        )?;
        println!(
            "Consensus polishing completed: {} bases modified. Output saved to: {:?}",
            corrected, out_contigs
        );
        out_contigs
    } else {
        cli.contigs
    };

    let metrics = evaluate_assembly(&target_eval_file, cli.reference.as_deref(), cli.min_len)?;

    metrics.print_summary();

    if let Some(json_path) = cli.json {
        std::fs::write(&json_path, metrics.to_json())?;
        println!("  Evaluation JSON report saved to: {:?}", json_path);
    }

    Ok(())
}
