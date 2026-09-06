# spades-rs

**Ultra-Fast, Memory-Shielded De Novo Genome Assembler in Pure Rust**

[![Rust CI](https://github.com/sagnikrout/spades-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/sagnikrout/spades-rs)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Release: v1.0.0](https://img.shields.io/badge/Release-v1.0.0-teal.svg)](https://github.com/sagnikrout/spades-rs/releases)

`spades-rs` is a high-performance, single-binary rewrite of the core algorithms behind classical *de novo* genome assemblers (like SPAdes). Built from the ground up in 100% safe, modern Rust for multi-core SIMD architectures, it eliminates the notorious memory bloat, process crashes, and CPU thrashing of legacy bioinformatics tools.

* **10× – 32× Less Memory**: Runs on standard laptops and workstations via lock-free Two-Tier Bloom filters and 2-bit packed streaming.
* **Automatic Hardware Memory Governor**: Dynamically locks budget to `Available System RAM - 20%`, preventing OOM crashes.
* **3× – 4× Faster Wall-Clock**: Assembles full bacterial isolates in ~4 minutes and viral genomes in under 1 second.
* **Complete Multi-Omics Suite**: Built-in support for WGS isolates, metagenomes (`--meta`), plasmids (`--plasmid`), RNA transcriptomes (`--rna`), single-cell MDA (`--sc`), hybrid long-read repeat bridging (`--nanopore`, `--pacbio`), and progressive multi-K iteration (`--multik`).

---

## 🚀 Key Performance Benchmarks

### 1. Full Bacterial Isolate Benchmark (*E. coli* MG1655, 100x depth, 1.28M PE reads, Multi-K 33,55,77,99,111)
Measured empirically on an **Intel Core Ultra 9 185H (22 logical threads, AVX2, 32 GB RAM)** against NCBI Reference `NC_000913.3` (4.64 Mb):

| Metric | Legacy SPAdes v4.3.0 | spades-rs (Rust) | Advantage / Delta |
| :--- | :--- | :--- | :--- |
| **Elapsed Wall-Clock** | 834.93 s (13m 55s) | **259.86 s (4m 20s)** | **3.21x Faster** |
| **User CPU Time** | 5,567.95 s (92m 48s) | **2,275.41 s (37m 55s)** | **2.45x Less CPU Compute** |
| **Peak RAM (Max RSS)** | 5,560.73 MB (**5.56 GB**) | **2,157.60 MB (2.16 GB)** | **2.58x Less Memory** |
| **Contigs ($\ge$ 200 bp)** | 700 contigs | **244 contigs** | **2.87x Fewer Fragments** |
| **Genome Fraction** | 99.19% | **98.82%** | **High-fidelity de novo representation** |

### 2. Hybrid Assembly with Long Reads (*E. coli* MG1655 + Oxford Nanopore `DRR242214`)
| Metric | Legacy SPAdes v4.3.0 | spades-rs (Rust) | Advantage / Delta |
| :--- | :--- | :--- | :--- |
| **Elapsed Wall-Clock** | 961.68 s (16m 01s) | **256.63 s (4m 17s)** | **3.75x Faster** |
| **Peak RAM (Max RSS)** | 5,278.59 MB (**5.28 GB**) | **2,115.82 MB (2.12 GB)** | **2.50x Less Memory** |
| **Scaffold N50 / L50** | N/A | **212,344 bp / 7** | **Extensive structural contiguity** |
| **Longest Scaffold** | 469,088 bp | **811,852 bp** | **+342.7 kb Longer Scaffold** |
| **Unaligned Contigs** | 468 contigs (374.5 kb) | **8 contigs (131.2 kb)** | **58.5x Fewer Junk Contigs** |
| **Base Accuracy (Mismatches/100kb)** | 9.15 (Q40) | **6.15** (3.78 in scaffolds) | **Up to 2.4x Higher Base Accuracy** |
| **7 Ribosomal RNA Operons** | Collapsed / Fragmented | **7 / 7 (100% BRIDGED)** | **Zero collapsed multi-copy repeats** |
| **Executable Size** | ~640 MB across 15+ binaries | **1.1 MB standalone native binary** | **580x Smaller, Zero dependencies** |

---

## 🏛️ Modular Architecture

The engine is decomposed into 5 distinct architectural layers across 16 self-contained modules:

```
[ Raw Gzipped Reads (.fq.gz) ]
              │
              ▼
┌───────────────────────────────────────────────────────────────┐
│ Layer 1: Ingestion & Error Correction                         │
│ • [Module 1.1] SIMD 2-Bit Streaming Fastq Reader              │
│ • [Module 1.2] Lock-Free Two-Tier Bloom Filter (Memory Shield)│
│ • [Module 1.3] BayesHammer 2-Bit POPCNT Error Corrector       │
└───────────────────────────────┬───────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────┐
│ Layer 2: Graph Topology & Simplification                      │
│ • [Module 2.1] Bidirected Compacted de Bruijn Graph (cDBG)    │
│ • [Module 2.2] Relative Coverage Tip-Clipper & Path Stitcher  │
└───────────────────────────────┬───────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────┐
│ Layer 3: Repeat Resolution & Hybrid Bridging                  │
│ • [Module 3.1] Library Insert Size Distance Estimator         │
│ • [Module 3.2] ExSPAnder Paired-End Repeat Resolver           │
│ • [Module 3.3] Spaligner Nanopore / PacBio Long-Read Bridging │
└───────────────────────────────┬───────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────┐
│ Layer 4: Scaffolding, Polishing & Export                      │
│ • [Module 4.1] Late Gap Closer & Scaffolder ('N' filling)     │
│ • [Module 4.2] Consensus Base Polisher (WFA consensus)        │
│ • [Module 4.3] Standard GFA v1.1 & FASTA Exporters            │
└───────────────────────────────┬───────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────┐
│ Layer 5: Specialized Biological Pipelines                     │
│ • [Module 5.1] metaSPAdes: Uneven multi-species preservation  │
│ • [Module 5.2] plasmidSPAdes: Circular & copy-number detector │
│ • [Module 5.3] rnaSPAdes: Transcriptome isoform preserver     │
│ • [Module 5.4] scSPAdes: Single-cell MDA normalizer           │
└───────────────────────────────────────────────────────────────┘
```

---

## 📦 Building & Installation

### Prerequisites
* Rust 1.75+ (`rustc` and `cargo`)

### Build from Source
```bash
git clone https://github.com/sagnikrout/spades-rs.git
cd spades-rs

# Optimized release build with LTO and native SIMD vectorization:
cargo build --release
```
The compiled standalone executable will be located at `target/release/spades-rs` (~1.1 MB).

---

## 💻 CLI Usage Guide

### 1. Standard High-Speed De Novo Assembly
```bash
./target/release/spades-rs assemble \
    -i reads_1.fq.gz reads_2.fq.gz \
    -o contigs.fasta \
    -k 31 \
    --error-correct
```
Outputs:
* `contigs.fasta`: Primary assembled genomic contigs.
* `scaffolds.fasta`: Scaffolds linked across unresolved repeat gaps.
* `assembly_graph.gfa`: Standard Graphical Fragment Assembly v1.1 (viewable in [Bandage](https://rrwick.github.io/Bandage/)).

### 2. Hardware Memory Governor (Available RAM - 20%)
By default, `spades-rs` auto-detects real-time system memory and bounds its working envelope to **80% of available RAM** (leaving 20% headroom for OS, disk cache, and glibc arenas). To specify a custom limit:
```bash
./target/release/spades-rs assemble \
    -i reads_1.fq.gz reads_2.fq.gz \
    --max-memory 4.5 \
    -o contigs.fasta
```

### 3. Hybrid Assembly (Illumina + Oxford Nanopore / PacBio)
```bash
./target/release/spades-rs assemble \
    -i short_reads_1.fq.gz short_reads_2.fq.gz \
    --nanopore ont_long_reads.fq.gz \
    -o hybrid_contigs.fasta
```

### 4. Metagenomic & Plasmid Assembly
```bash
./target/release/spades-rs assemble \
    -i metagenome_1.fq.gz metagenome_2.fq.gz \
    --meta \
    --plasmid \
    -o meta_assembly.fasta
```
Extracts circular and high-copy elements into `plasmids.fasta` while retaining low-abundance species.

### 5. Transcriptome Mode (RNA-Seq)
```bash
./target/release/spades-rs assemble \
    -i rna_1.fq.gz rna_2.fq.gz \
    --rna \
    -o transcripts.fasta
```

### 6. Multi-K Progressive Iterative Stepping
```bash
./target/release/spades-rs assemble \
    -i reads_1.fq.gz reads_2.fq.gz \
    --multik 21,33,55 \
    -o contigs.fasta
```

---

## 🔬 Benchmark Verification & Technical Report

All 11 sequential ground-truth benchmarks against published biological references are documented in detail with QUAST scorecards in:
* **[`TECHNICAL_REPORT_AND_ROADMAP.md`](TECHNICAL_REPORT_AND_ROADMAP.md)**

---

## 🧪 Automated Testing

Run the full test suite (39 unit tests, stress tests, and end-to-end reference genome verifications):
```bash
cargo test
```
All 39 tests pass with 0 failures:
```text
test result: ok. 39 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 18.80s
```

---

## 📄 License
MIT License. Free for academic, non-commercial, and commercial genomics workflows.
