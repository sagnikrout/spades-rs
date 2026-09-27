# spades-rs

A Rust implementation of the SPAdes genome assembly pipeline designed for lower memory consumption.

[![Rust CI](https://github.com/sagnikrout/spades-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/sagnikrout/spades-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Release: v1.0.0](https://img.shields.io/badge/Release-v1.0.0-teal.svg)](https://github.com/sagnikrout/spades-rs/releases)

`spades-rs` is a standalone *de novo* genome assembler written in Rust. It implements the primary algorithms of the SPAdes pipeline (Bankevich et al., 2012) as a single static executable, with a focus on avoiding disk churn and keeping peak memory low enough to run on standard computers.

## Practical expectations

Assembly performance depends heavily on library quality, coverage depth, repeat content, and error profiles. In everyday practice, expect the following:

* **Memory usage:** For standard bacterial genomes (2 to 6 Mb) at 30x to 100x coverage, peak RAM typically remains between 1.5 GB and 3.5 GB. Classical SPAdes often allocates 8 GB to 16 GB for the same data. On larger datasets, memory scales with the number of distinct k-mers.
* **Execution speed:** On modern multi-core processors, runs are generally 1.5x to 3x faster than legacy SPAdes on clean bacterial isolates, mainly because multi-k iterations proceed in memory without saving intermediate graph states to disk.
* **Assembly quality:** On clean short-read bacterial data with adequate coverage, contiguity (N50) and genome fraction are generally comparable to classical SPAdes (typically 95% to 98% reference coverage). Repeat regions longer than the read length or insert size will fragment contigs unless bridged by long reads.
* **Specialized modes:** Basic implementations of `--meta` (metagenomics), `--plasmid` (plasmid extraction), `--rna` (transcriptomes), and `--sc` (single-cell MDA) are included. For large, complex environmental metagenomes or specialized clinical pipelines, results should be evaluated alongside the mature C++ SPAdes release.

## Installation

### Installer script (Linux, WSL, macOS, Google Colab)

Install the compiled binary to `/usr/local/bin` (or `~/.local/bin`):

```bash
curl -fsSL https://raw.githubusercontent.com/sagnikrout/spades-rs/master/install.sh | bash
```

**In Google Colab**, run in any notebook cell:
```python
!curl -fsSL https://raw.githubusercontent.com/sagnikrout/spades-rs/master/install.sh | bash
!spades-rs -1 reads_1.fq.gz -2 reads_2.fq.gz -o out_dir
```

### Precompiled binaries

Download directly from [GitHub Releases](https://github.com/sagnikrout/spades-rs/releases):

| Platform | Target triple | Binary |
| :--- | :--- | :--- |
| Linux x86_64 (musl static) | `x86_64-unknown-linux-musl` | [`spades-rs-linux-x86_64-musl`](https://github.com/sagnikrout/spades-rs/releases/latest/download/spades-rs-linux-x86_64-musl) |
| Linux ARM64 | `aarch64-unknown-linux-gnu` | [`spades-rs-linux-arm64`](https://github.com/sagnikrout/spades-rs/releases/latest/download/spades-rs-linux-arm64) |
| macOS Apple Silicon | `aarch64-apple-darwin` | [`spades-rs-macos-arm64`](https://github.com/sagnikrout/spades-rs/releases/latest/download/spades-rs-macos-arm64) |
| Windows x86_64 | `x86_64-pc-windows-msvc` | [`spades-rs-windows-x86_64.exe`](https://github.com/sagnikrout/spades-rs/releases/latest/download/spades-rs-windows-x86_64.exe) |

```bash
# Manual installation on 64-bit Linux
curl -fsSL https://github.com/sagnikrout/spades-rs/releases/latest/download/spades-rs-linux-x86_64-musl -o /usr/local/bin/spades-rs
chmod +x /usr/local/bin/spades-rs
```

### Build from source

Building requires Rust 1.75 or newer on a 64-bit platform:

```bash
git clone https://github.com/sagnikrout/spades-rs.git
cd spades-rs
cargo build --release
```

The compiled binary will be placed at `target/release/spades-rs`.

## Pipeline architecture

`spades-rs` structures the assembly process into sequential in-memory stages:

```
[ Raw FASTQ / FASTA (.gz or plain) ]
                  │
                  ▼
┌───────────────────────────────────────────────────────────────┐
│ Stage 1: Ingestion & 2-Bit Packing                            │
│ • End trimming of poly-N runs and low-quality bases           │
│ • Internal N splitting into clean sub-reads (>= 21 bp)        │
│ • Ambiguous sidecar preservation for consensus base voting    │
│ • Contiguous 2-bit buffer (4 bases per byte)                  │
└───────────────────────────────┬───────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────┐
│ Stage 2 & 3: Counting & Solid K-mers                          │
│ • Lock-free Two-Tier Atomic Bloom filter (64 MB to 4 GB)      │
│ • Sharded solid k-mer frequency index                         │
│ • Automatic noise cutoff valley detection                     │
└───────────────────────────────┬───────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────┐
│ Stage 4 & 5: Graph Construction & Simplification              │
│ • Compacted de Bruijn graph (cDBG) with 256-bit SIMD k-mers   │
│ • Tip clipping and bubble popping                             │
│ • Disjoint path stitching into unitigs                        │
└───────────────────────────────┬───────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────┐
│ Stage 6 & 7: Repeat Resolution & Scaffolding                  │
│ • ExSPAnder paired-end insert size estimation and pathing     │
│ • Spaligner hybrid repeat unrolling (Nanopore / PacBio HiFi)  │
│ • SpLitteR linked-read barcode repeat resolution (10x)        │
│ • Local de Bruijn gap closure with N-bridge insertion         │
└───────────────────────────────┬───────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────┐
│ Stage 8: Polishing & Output Export                            │
│ • Consensus base voting across primary reads and sidecar      │
│ • FASTA export (contigs.fasta, scaffolds.fasta)               │
│ • Assembly graph export in standard GFA v1.1 format           │
└───────────────────────────────────────────────────────────────┘
```

## Command-line options

`spades-rs` supports direct invocation where `assemble` is the default subcommand:

```bash
# Direct execution
spades-rs -1 R1.fq.gz -2 R2.fq.gz -o output_dir

# Explicit subcommand syntax
spades-rs assemble -1 R1.fq.gz -2 R2.fq.gz -o output_dir
```

| Option | Argument | Default | Description |
| :--- | :--- | :--- | :--- |
| `-1, --pe1-1` | Path | None | Forward paired-end reads |
| `-2, --pe1-2` | Path | None | Reverse paired-end reads |
| `-s, --pe1-s` | Path | None | Unpaired / single-end reads |
| `--12` | Path | None | Interleaved paired-end reads |
| `-i, --inputs` | Paths | None | General FASTQ/FASTA files (plain or `.gz`) |
| `-o, --output` | Path | `contigs.fasta` | Output directory or contig FASTA path |
| `-k, --k` | Integer | `31` | K-mer size (must be odd, 11 to 127) |
| `-c, --coverage` | Float | `5.0` | Minimum unitig k-mer coverage threshold |
| `-m, --min-len` | Integer | `200` | Minimum contig length in output (bp) |
| `-t, --threads` | Integer | System threads | Worker threads (defaults to logical core count) |
| `--error-correct` | Flag | `false` | Run BayesHammer read error correction before assembly |
| `--careful` | Flag | `false` | Run careful mode with stricter mismatch correction |
| `--multik` | List | None | Comma-separated k-mer sizes (e.g., `21,33,55`) |
| `--nanopore` | Path | None | Oxford Nanopore reads for hybrid bridging |
| `--pacbio` | Path | None | PacBio HiFi reads for hybrid bridging |
| `--splitter` | Paths | None | Barcoded linked reads for SpLitteR repeat resolution |
| `--trusted-contigs` | Paths | None | High-confidence contigs to guide assembly paths |
| `--max-memory` | Float | 80% RAM | Maximum physical RAM budget in GB |
| `--meta` | Flag | `false` | Enable metagenomic uneven coverage filtering |
| `--plasmid` | Flag | `false` | Segregate circular and high-copy elements into `plasmids.fasta` |
| `--rna` | Flag | `false` | Enable transcriptome mode preserving alternative splicing |
| `--sc` | Flag | `false` | Enable single-cell MDA normalization |
| `--no-polish` | Flag | `false` | Skip consensus base polishing pass |

## Usage examples

### 1. Standard bacterial isolate assembly

```bash
spades-rs -1 reads_1.fq.gz -2 reads_2.fq.gz --careful -o output_dir/
```

Files produced in `output_dir/`:
* `contigs.fasta`: Primary assembled genomic contigs.
* `scaffolds.fasta`: Scaffolds linked across repeat gaps.
* `assembly_graph.gfa`: Standard Graphical Fragment Assembly v1.1 file, viewable in Bandage.

### 2. Multi-k iteration

Stepping through multiple k-mer sizes helps resolve both short repeats and low-coverage regions:

```bash
spades-rs -1 reads_1.fq.gz -2 reads_2.fq.gz --multik 21,33,55 -o output_dir/
```

### 3. Setting a memory ceiling

By default, `spades-rs` inspects host RAM and caps itself at 80% of available memory. A hard ceiling can be set explicitly:

```bash
spades-rs -1 reads_1.fq.gz -2 reads_2.fq.gz --max-memory 2.0 -o output_dir/
```

### 4. Hybrid assembly with long reads

```bash
spades-rs -1 short_1.fq.gz -2 short_2.fq.gz --nanopore ont_reads.fq.gz -o hybrid_out/
```

### 5. Metagenomic and plasmid assembly

```bash
# Metagenomic community
spades-rs -1 meta_1.fq.gz -2 meta_2.fq.gz --meta -o meta_out/

# Plasmid extraction
spades-rs -1 isolate_1.fq.gz -2 isolate_2.fq.gz --plasmid -o plasmid_out/
```

## Testing and validation

The repository includes both native Rust tests and a real-world scenario validation harness:

### Native Rust test suite

```bash
cargo test --release
```

Includes 53 tests covering DNA primitives, 256-bit SIMD k-mers, Two-Tier Bloom filters, graph simplification, 2-bit read packing, N-base handling, and module workflows.

### Real-world scenario validation

```bash
# Fast validation across all library types and edge cases (~20 seconds)
python3 tools/run_realworld_validation.py --quick

# Full matrix including multi-million read benchmarks
python3 tools/run_realworld_validation.py
```

Validates 19 distinct operational scenarios:
1. Viral control ground truth (PhiX174, circular ssDNA)
2. Single-end reads (`-s`)
3. Interleaved paired-end reads (`--12`)
4. Legacy Phred+64 auto-detection
5. Multi-line wrapped FASTA input
6. Deep multi-K stepping (k=21 to k=99, testing 256-bit SIMD math)
7. Hybrid PacBio HiFi repeat bridging (`--pacbio`)
8. 10x Genomics barcoded linked reads (`--splitter`)
9. Multi-plasmid segregation (`--plasmid`)
10. Single-cell MDA depth normalization (`--sc`)
11. RNA-Seq transcriptome assembly (`--rna`)
12. Metagenome multi-coverage filtering (`--meta`)
13. Memory budget enforcement (`--max-memory 0.5`)
14. Thread concurrency scaling (`-t 1`)
15. Reads with degraded poly-N tails
16. Reads with 100% N bases
17. Even k-mer rejection fault injection
18. Out-of-bounds k-mer rejection fault injection
19. Missing input file error handling

## Codebase layout

The Rust source code is contained in `src/` (4,600 lines, zero dynamic runtime dependencies):

| File | Primary role |
| :--- | :--- |
| [`src/main.rs`](src/main.rs) | CLI parsing, argument validation, thread pool setup |
| [`src/lib.rs`](src/lib.rs) | Crate root, 64-bit architecture assertion, public exports |
| [`src/dna.rs`](src/dna.rs) | 2-bit base encoding, 64-bit and 256-bit canonical k-mer arithmetic |
| [`src/bloom.rs`](src/bloom.rs) | Two-Tier Atomic Bloom filter for memory-bounded k-mer filtering |
| [`src/fastq.rs`](src/fastq.rs) | FASTQ/FASTA reader, Phred-33/64 auto-detection, quality trimming |
| [`src/packed_reads.rs`](src/packed_reads.rs) | Contiguous 2-bit read buffer with end-trimming and N-handling |
| [`src/hammer.rs`](src/hammer.rs) | BayesHammer Hamming error correction |
| [`src/graph.rs`](src/graph.rs) | Compacted de Bruijn graph with bidirected port involution |
| [`src/simplify.rs`](src/simplify.rs) | Tip clipping, bubble popping, and unitig stitching |
| [`src/paired_info.rs`](src/paired_info.rs) | Paired-end library insert size estimation |
| [`src/expander.rs`](src/expander.rs) | ExSPAnder paired-end repeat navigation |
| [`src/spaligner.rs`](src/spaligner.rs) | Spaligner long-read graph alignment and repeat bridging |
| [`src/splitter.rs`](src/splitter.rs) | SpLitteR linked-read barcode repeat resolution |
| [`src/multik.rs`](src/multik.rs) | Progressive in-memory multi-k unitig seeding |
| [`src/scaffold.rs`](src/scaffold.rs) | Scaffolding path walker and gap closing |
| [`src/polisher.rs`](src/polisher.rs) | Consensus base voting across primary reads and sidecar |
| [`src/modes.rs`](src/modes.rs) | Specialized pipelines: meta, plasmid, RNA, and single-cell |
| [`src/memory.rs`](src/memory.rs) | Memory budget calculations and dynamic Bloom sizing |
| [`src/gfa.rs`](src/gfa.rs) | Graphical Fragment Assembly (GFA v1.1) exporter |
| [`src/assemble.rs`](src/assemble.rs) | End-to-end 8-stage assembly pipeline orchestrator |

## Attribution and citations

`spades-rs` is an independent Rust reimplementation based on the algorithms developed by the SPAdes research team (Algorithmic Biology Lab, St. Petersburg Academic University / Center for Algorithmic Biotechnology).

When using `spades-rs`, please cite both this repository and the original publications:

* **SPAdes:** Bankevich, A. et al. (2012). *Journal of Computational Biology*, 19(5), 455–477. [doi:10.1089/cmb.2012.0021](https://doi.org/10.1089/cmb.2012.0021).
* **Using SPAdes:** Prjibelski, A. et al. (2020). *Current Protocols in Bioinformatics*, 70(1), e102. [doi:10.1002/cpbi.102](https://doi.org/10.1002/cpbi.102).
* **BayesHammer:** Nikolenko, S. I. et al. (2013). *BMC Genomics*, 14(Suppl 1), S7. [doi:10.1186/1471-2164-14-S1-S7](https://doi.org/10.1186/1471-2164-14-S1-S7).
* **ExSPAnder:** Prjibelski, A. D. et al. (2014). *Bioinformatics*, 30(12), i293–i301. [doi:10.1093/bioinformatics/btu266](https://doi.org/10.1093/bioinformatics/btu266).
* **metaSPAdes:** Nurk, S. et al. (2017). *Genome Research*, 27(5), 824–834. [doi:10.1101/gr.213959.116](https://doi.org/10.1101/gr.213959.116).
* **plasmidSPAdes:** Antipov, D. et al. (2016). *Bioinformatics*, 32(22), 3380–3387. [doi:10.1093/bioinformatics/btw493](https://doi.org/10.1093/bioinformatics/btw493).
* **rnaSPAdes:** Bushmanova, E. et al. (2019). *GigaScience*, 8(9), giz100. [doi:10.1093/gigascience/giz100](https://doi.org/10.1093/gigascience/giz100).
* **hybridSPAdes / Spaligner:** Antipov, D. et al. (2016). *Bioinformatics*, 32(7), 1009–1015; Dvorkina, T. et al. (2020). *Bioinformatics*, 36(Suppl 1), i188–i195.
* **SpLitteR:** Tolstoganov, I. et al. (2024). *PeerJ*, 12, e18050. [doi:10.7717/peerj.18050](https://doi.org/10.7717/peerj.18050).

## License

MIT License. Free for academic, non-commercial, and commercial research.
