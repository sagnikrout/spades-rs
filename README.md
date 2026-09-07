# spades-rs

A Rust-based de novo genome assembler designed for low-memory environments.

[![Rust CI](https://github.com/sagnikrout/spades-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/sagnikrout/spades-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Release: v1.0.0](https://img.shields.io/badge/Release-v1.0.0-teal.svg)](https://github.com/sagnikrout/spades-rs/releases)

`spades-rs` is an experimental *de novo* genome assembler written in Rust. It implements core algorithmic concepts from the SPAdes assembly pipeline (Bankevich et al., 2012), with a focus on reduced memory consumption and single-binary deployment on commodity 64-bit hardware.

* **Bounded memory footprint**: Uses lock-free Two-Tier Bloom filters and 2-bit packed reads to operate within limited RAM environments.
* **Automatic memory budgeting**: Inspects available physical memory at launch and defaults to an 80% budget ceiling, leaving 20% headroom for operating system processes and file cache.
* **Single binary**: Compiles into an independent executable with no external dynamic library dependencies.
* **Pipeline modes**: Supports short-read isolate assembly, progressive multi-K stepping, long-read hybrid bridging (`--nanopore`, `--pacbio`), and preliminary support for metagenomic (`--meta`), plasmid (`--plasmid`), RNA-Seq (`--rna`), and single-cell (`--sc`) datasets.

> Note: `spades-rs` is actively evolving software. While it reproduces key SPAdes heuristics, legacy SPAdes remains the mature reference standard for production genomics pipelines.

## Quickstart

### Precompiled binaries

Precompiled standalone binaries for 64-bit platforms are available from GitHub releases:

| Platform | Target triple | Download link |
| :--- | :--- | :--- |
| Linux x86_64 | `x86_64-unknown-linux-musl` | [`spades-rs-linux-musl`](https://github.com/sagnikrout/spades-rs/releases/latest/download/spades-rs-linux-musl) |
| Linux ARM64 | `aarch64-unknown-linux-gnu` | [`spades-rs-linux-arm64`](https://github.com/sagnikrout/spades-rs/releases/latest/download/spades-rs-linux-arm64) |
| macOS Apple Silicon | `aarch64-apple-darwin` | [`spades-rs-macos-arm64`](https://github.com/sagnikrout/spades-rs/releases/latest/download/spades-rs-macos-arm64) |
| Windows x86_64 | `x86_64-pc-windows-msvc` | [`spades-rs-windows-x86_64.exe`](https://github.com/sagnikrout/spades-rs/releases/latest/download/spades-rs-windows-x86_64.exe) |

```bash
# Download and execute on 64-bit Linux
curl -L -o spades-rs https://github.com/sagnikrout/spades-rs/releases/latest/download/spades-rs-linux-musl
chmod +x spades-rs
./spades-rs assemble -i reads_1.fq.gz reads_2.fq.gz -o contigs.fasta
```

### Build from source

Building from source requires Rust 1.75 or later on a 64-bit operating system:

```bash
git clone https://github.com/sagnikrout/spades-rs.git
cd spades-rs
cargo build --release
```
The resulting executable is located at `target/release/spades-rs`.

## Codebase and file structure

The codebase is partitioned into 19 modules in `src/`, automated tests in `tests/`, and verification scripts in `tools/`:

| Path | Primary responsibility | Key algorithms and mechanisms |
| :--- | :--- | :--- |
| [`src/main.rs`](src/main.rs) | CLI entry point and thread pool configuration | Argument parsing (`clap`), Rayon thread pool initialization |
| [`src/lib.rs`](src/lib.rs) | Crate root and architecture checks | Compile-time 64-bit pointer assertion, module re-exports |
| [`src/dna.rs`](src/dna.rs) | Nucleotide encoding and k-mer hashing | 2-bit representation (A=0, C=1, G=2, T=3), 64-bit and 256-bit canonical k-mers |
| [`src/bloom.rs`](src/bloom.rs) | Memory shield for k-mer counting | Lock-free Two-Tier Atomic Bloom filter (512 MB fixed bitset) |
| [`src/fastq.rs`](src/fastq.rs) | FASTQ and FASTA sequence ingestion | Multi-threaded streaming gzip and plain text sequence parser |
| [`src/packed_reads.rs`](src/packed_reads.rs) | Compact read storage | Cache-aligned 2-bit packed array (8.0M reads in 352 MB RAM) |
| [`src/hammer.rs`](src/hammer.rs) | Read error correction | BayesHammer algorithm: bit-parallel Hamming clustering and quality voting |
| [`src/graph.rs`](src/graph.rs) | Compacted de Bruijn graph (cDBG) | Unitig topology, bidirected port involution (2N ports), adjacency indices |
| [`src/simplify.rs`](src/simplify.rs) | Graph cleaning and topology simplification | Length-aware tip clipping (>2k protected), bubble popping, linear unitig stitching |
| [`src/paired_info.rs`](src/paired_info.rs) | Paired-end insert size estimation | Gaussian distance distribution estimation (mean and standard deviation) |
| [`src/expander.rs`](src/expander.rs) | Paired-end repeat resolution | ExSPAnder algorithm: Dijkstra path extension with paired-read voting |
| [`src/spaligner.rs`](src/spaligner.rs) | Long-read hybrid repeat bridging | Spaligner algorithm: ONT and PacBio graph alignment and repeat unrolling |
| [`src/multik.rs`](src/multik.rs) | Multi-K progressive iteration | Progressive unitig-to-reads seeding loop across increasing k-mer sizes |
| [`src/scaffold.rs`](src/scaffold.rs) | Scaffolding and gap closing | Local de Bruijn path walker with cycle guards, insertion of 'N' bridges |
| [`src/polisher.rs`](src/polisher.rs) | Consensus base polishing | Wavefront-style consensus alignment against raw reads |
| [`src/modes.rs`](src/modes.rs) | Metagenomics and plasmid pipelines | metaSPAdes coverage filtering and plasmidSPAdes circularity extraction |
| [`src/rna.rs`](src/rna.rs) | Transcriptome assembly pipeline | rnaSPAdes alternative isoform preservation and transcript extraction |
| [`src/single_cell.rs`](src/single_cell.rs) | Single-cell MDA normalization | scSPAdes local coverage normalization for severe amplification bias |
| [`src/memory.rs`](src/memory.rs) | Hardware memory governor | Real-time RAM detection, 20% OS headroom protection, dynamic Bloom sizing |
| [`src/gfa.rs`](src/gfa.rs) | Graph visualization export | Graphical Fragment Assembly (GFA v1.1) exporter for Bandage |
| [`src/assemble.rs`](src/assemble.rs) | Top-level assembly orchestration | Coordinates ingestion, correction, graph construction, resolution, and output |
| [`tests/`](tests/) | Integration test suite | 39 automated tests covering unitigs, graph simplification, and full genomes |
| [`tools/`](tools/) | Biological audit scripts | Reference-based QUAST evaluation, AMR gene checks, and rRNA synteny audits |

## Assembly pipeline

```
[ Raw Gzipped Reads (.fq.gz) ]
              │
              ▼
┌───────────────────────────────────────────────────────────────┐
│ Layer 1: Ingestion and Error Correction                       │
│ • [Module 1.1] SIMD 2-Bit Streaming Fastq Reader              │
│ • [Module 1.2] Lock-Free Two-Tier Bloom Filter (Memory Shield)│
│ • [Module 1.3] BayesHammer 2-Bit POPCNT Error Corrector       │
└───────────────────────────────┬───────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────┐
│ Layer 2: Graph Topology and Simplification                    │
│ • [Module 2.1] Bidirected Compacted de Bruijn Graph (cDBG)    │
│ • [Module 2.2] Relative Coverage Tip-Clipper & Path Stitcher  │
└───────────────────────────────┬───────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────┐
│ Layer 3: Repeat Resolution and Hybrid Bridging                │
│ • [Module 3.1] Library Insert Size Distance Estimator         │
│ • [Module 3.2] ExSPAnder Paired-End Repeat Resolver           │
│ • [Module 3.3] Spaligner Nanopore / PacBio Long-Read Bridging │
└───────────────────────────────┬───────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────┐
│ Layer 4: Scaffolding, Polishing and Export                    │
│ • [Module 4.1] Late Gap Closer & Scaffolder ('N' filling)     │
│ • [Module 4.2] Consensus Base Polisher (WFA consensus)        │
│ • [Module 4.3] Standard GFA v1.1 and FASTA Exporters          │
└───────────────────────────────┬───────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────┐
│ Layer 5: Specialized Biological Pipelines                     │
│ • [Module 5.1] metaSPAdes: Uneven multi-species preservation  │
│ • [Module 5.2] plasmidSPAdes: Circular and copy-number finder │
│ • [Module 5.3] rnaSPAdes: Transcriptome isoform preserver     │
│ • [Module 5.4] scSPAdes: Single-cell MDA normalizer           │
└───────────────────────────────────────────────────────────────┘
```

## Performance comparisons

These tests were performed on an Intel Core Ultra 9 185H (22 logical threads, AVX2, 32 GB RAM) running Ubuntu under WSL2.

### 1. Bacterial isolate (Escherichia coli MG1655, 100x depth, 1.28M paired reads, Multi-K 33,55,77,99,111)
Evaluated against NCBI Reference `NC_000913.3` (4.64 Mb):

| Metric | SPAdes v4.3.0 | spades-rs (v1.0.0) | Observation |
| :--- | :--- | :--- | :--- |
| Wall-Clock Time | 834.93 s (13m 55s) | 259.86 s (4m 20s) | Shorter runtime |
| User CPU Time | 5,567.95 s (92m 48s) | 2,275.41 s (37m 55s) | Reduced CPU time |
| Peak RAM (Max RSS) | 5,560.73 MB (5.56 GB) | 2,157.60 MB (2.16 GB) | Lower peak memory |
| Contigs (≥ 200 bp) | 700 contigs | 244 contigs | Fewer small fragments |
| Genome Fraction | 99.19% | 98.82% | Comparable recovery |

### 2. Hybrid assembly (Escherichia coli MG1655 paired-end + Oxford Nanopore DRR242214)

| Metric | SPAdes v4.3.0 | spades-rs (v1.0.0) | Observation |
| :--- | :--- | :--- | :--- |
| Wall-Clock Time | 961.68 s (16m 01s) | 256.63 s (4m 17s) | Shorter runtime |
| Peak RAM (Max RSS) | 5,278.59 MB (5.28 GB) | 2,115.82 MB (2.12 GB) | Lower peak memory |
| Scaffold N50 / L50 | Not produced | 212,344 bp / 7 | Produced scaffold output |
| Longest Scaffold | 469,088 bp | 811,852 bp | Longer primary scaffold |
| Unaligned Contigs | 468 contigs (374.5 kb) | 8 contigs (131.2 kb) | Fewer unaligned contigs |
| Base Error Rate | 9.15 mismatches / 100 kb | 6.15 mismatches / 100 kb | Lower mismatch rate |
| 7 rRNA Operons | Fragmented | 7 / 7 bridged | Bridged repeat copies |
| Installation Size | Multiple binaries (~640 MB) | Single binary (1.1 MB) | Self-contained binary |

*Note: Assembly results, runtime, and memory consumption vary with dataset characteristics, sequencing depth, error profiles, and hardware. Users should evaluate results against their own quality criteria.*

## Command-line options

```
spades-rs assemble [OPTIONS] -i <INPUTS>...
```

| Option | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `-i, --inputs <PATHS>` | Paths | Required | Input paired or single-end FASTQ/FASTA files (plain or `.gz`) |
| `-o, --output <PATH>` | Path | `contigs.fasta` | Output contigs FASTA file |
| `-k, --k <INT>` | Integer | `31` | Primary k-mer size (must be odd, <= 127) |
| `-c, --coverage <FLOAT>` | Float | `5.0` | Minimum unitig k-mer coverage threshold |
| `-m, --min-len <INT>` | Integer | `200` | Minimum contig length in output (bp) |
| `-t, --threads <INT>` | Integer | Auto | Worker thread count (defaults to logical CPU cores) |
| `--error-correct` | Flag | `false` | Run BayesHammer read error correction prior to assembly |
| `--max-memory <GB>` | Float | Auto (80% RAM) | Hard RAM budget ceiling in GB |
| `--multik <LIST>` | Comma-separated | None | Progressive multi-K sizes (e.g., `21,33,55,77`) |
| `--nanopore <PATH>` | Path | None | Oxford Nanopore reads for hybrid bridging |
| `--pacbio <PATH>` | Path | None | PacBio reads for hybrid bridging |
| `--meta` | Flag | `false` | Enable metaSPAdes mode for uneven metagenomic communities |
| `--plasmid` | Flag | `false` | Enable plasmidSPAdes mode to extract plasmids to `plasmids.fasta` |
| `--rna` | Flag | `false` | Enable rnaSPAdes mode for transcriptomes and alternative isoforms |
| `--sc` | Flag | `false` | Enable scSPAdes mode for single-cell MDA normalization |
| `--no-polish` | Flag | `false` | Skip final consensus base polishing pass |

## Usage examples

### 1. Standard isolate assembly with error correction
```bash
spades-rs assemble \
    -i reads_1.fq.gz reads_2.fq.gz \
    -o contigs.fasta \
    -k 31 \
    --error-correct
```
Outputs produced:
* `contigs.fasta`: Primary assembled genomic contigs.
* `scaffolds.fasta`: Scaffolds linked across repeat gaps.
* `assembly_graph.gfa`: Graphical Fragment Assembly v1.1 format (compatible with Bandage).

### 2. Specifying a memory ceiling
By default, `spades-rs` uses up to 80% of detected available RAM. A specific limit can be set manually:
```bash
spades-rs assemble \
    -i reads_1.fq.gz reads_2.fq.gz \
    --max-memory 4.0 \
    -o contigs.fasta
```

### 3. Hybrid assembly with Oxford Nanopore or PacBio
```bash
spades-rs assemble \
    -i short_reads_1.fq.gz short_reads_2.fq.gz \
    --nanopore ont_reads.fq.gz \
    -o hybrid_contigs.fasta
```

### 4. Metagenomic and plasmid extraction
```bash
spades-rs assemble \
    -i metagenome_1.fq.gz metagenome_2.fq.gz \
    --meta \
    --plasmid \
    -o meta_assembly.fasta
```

### 5. Transcriptome assembly (RNA-Seq)
```bash
spades-rs assemble \
    -i rna_1.fq.gz rna_2.fq.gz \
    --rna \
    -o transcripts.fasta
```

### 6. Progressive multi-K iterations
```bash
spades-rs assemble \
    -i reads_1.fq.gz reads_2.fq.gz \
    --multik 21,33,55 \
    -o contigs.fasta
```

## Technical report and verification

Additional evaluation details, including per-species metrics and biological feature recovery across 11 test sets, are documented in:
* [`TECHNICAL_REPORT_AND_ROADMAP.md`](TECHNICAL_REPORT_AND_ROADMAP.md)

## Automated testing

The test suite includes 39 unit and integration tests:
```bash
cargo test
```
All tests pass:
```text
test result: ok. 39 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 18.80s
```

## Academic citations and attribution

`spades-rs` is a Rust reimplementation and optimization of algorithmic concepts developed by the SPAdes research group (Algorithmic Biology Lab, St. Petersburg Academic University / Center for Algorithmic Biotechnology).

When using `spades-rs`, please cite both this repository and the original literature describing the foundational algorithms:

### Primary SPAdes algorithms
* **SPAdes Core Assembler:**  
  Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., Lesin, V. M., Nikolenko, S. I., Pham, S., Prjibelski, A. D., Pyshkin, A. V., Sirotkin, A. V., Vyahhi, N., Tesler, G., Alekseyev, M. A., & Pevzner, P. A. (2012). **SPAdes: A new genome assembly algorithm and its applications to single-cell sequencing.** *Journal of Computational Biology*, 19(5), 455–477. [doi:10.1089/cmb.2012.0021](https://doi.org/10.1089/cmb.2012.0021).

* **Comprehensive SPAdes Ecosystem:**  
  Prjibelski, A., Antipov, D., Meleshko, D., Lapidus, A., & Korobeynikov, A. (2020). **Using SPAdes De Novo Assembler.** *Current Protocols in Bioinformatics*, 70(1), e102. [doi:10.1002/cpbi.102](https://doi.org/10.1002/cpbi.102).

### Specialized algorithmic pipelines
* **BayesHammer Error Correction (`--error-correct`):**  
  Nikolenko, S. I., Korobeynikov, A. I., & Alekseyev, M. A. (2013). **BayesHammer: Bayesian clustering for error correction in single-cell sequencing.** *BMC Genomics*, 14(Suppl 1), S7. [doi:10.1186/1471-2164-14-S1-S7](https://doi.org/10.1186/1471-2164-14-S1-S7).

* **ExSPAnder Repeat Resolution (`src/expander.rs`):**  
  Prjibelski, A. D., Vasilinetc, I., Bankevich, A., Gurevich, A., Krivosheev, T., Nurk, S., Pham, S., & Pevzner, P. A. (2014). **ExSPAnder: a universal repeat resolver for DNA fragment assembly.** *Bioinformatics*, 30(12), i293–i301. [doi:10.1093/bioinformatics/btu266](https://doi.org/10.1093/bioinformatics/btu266).

* **metaSPAdes Metagenomics Pipeline (`--meta`):**  
  Nurk, S., Meleshko, D., Korobeynikov, A., & Pevzner, P. A. (2017). **metaSPAdes: a new versatile metagenomic assembler.** *Genome Research*, 27(5), 824–834. [doi:10.1101/gr.213959.116](https://doi.org/10.1101/gr.213959.116).

* **plasmidSPAdes Plasmid Extraction (`--plasmid`):**  
  Antipov, D., Hartwick, N., Shen, M., & Pevzner, P. A. (2016). **plasmidSPAdes: assembling plasmids from whole genome sequencing data.** *Bioinformatics*, 32(22), 3380–3387. [doi:10.1093/bioinformatics/btw493](https://doi.org/10.1093/bioinformatics/btw493).

* **rnaSPAdes Transcriptome Assembly (`--rna`):**  
  Bushmanova, E., Antipov, D., Lapidus, A., & Prjibelski, A. D. (2019). **rnaSPAdes: a de novo transcriptome assembler and its application to RNA-Seq data.** *GigaScience*, 8(9), giz100. [doi:10.1093/gigascience/giz100](https://doi.org/10.1093/gigascience/giz100).

* **hybridSPAdes and Spaligner Long-Read Graph Alignment (`--nanopore`, `--pacbio`):**  
  Antipov, D., Korobeynikov, A., McLean, J. S., & Pevzner, P. A. (2016). **hybridSPAdes: an algorithm for hybrid assembly of short and long reads.** *Bioinformatics*, 32(7), 1009–1015. [doi:10.1093/bioinformatics/btv688](https://doi.org/10.1093/bioinformatics/btv688).  
  Dvorkina, T., Antipov, D., & Korobeynikov, A. (2020). **Spaligner: alignment of long reads to assembly graphs.** *Bioinformatics*, 36(Suppl 1), i188–i195. [doi:10.1093/bioinformatics/btaa444](https://doi.org/10.1093/bioinformatics/btaa444).

### Computational and model assistance
* **Gemini 3.8 Flash (Google DeepMind):**  
  Interactive LLM assistance was used for code translation, test authoring, and documentation auditing.

### BibTeX entries
```bibtex
@article{bankevich2012spades,
  author    = {Bankevich, Anton and Nurk, Sergey and Antipov, Dmitry and Gurevich, Alexey A. and Dvorkin, Mikhail and Kulikov, Alexander S. and Lesin, Valery M. and Nikolenko, Sergey I. and Pham, Son and Prjibelski, Andrey D. and Pyshkin, Alexey V. and Sirotkin, Alexander V. and Vyahhi, Nikolay and Tesler, Glenn and Alekseyev, Max A. and Pevzner, Pavel A.},
  title     = {{SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing}},
  journal   = {Journal of Computational Biology},
  volume    = {19},
  number    = {5},
  pages     = {455--477},
  year      = {2012},
  doi       = {10.1089/cmb.2012.0021}
}

@article{prjibelski2020using,
  author    = {Prjibelski, Andrey and Antipov, Dmitry and Meleshko, Dmitry and Lapidus, Alla and Korobeynikov, Anton},
  title     = {{Using SPAdes De Novo Assembler}},
  journal   = {Current Protocols in Bioinformatics},
  volume    = {70},
  number    = {1},
  pages     = {e102},
  year      = {2020},
  doi       = {10.1002/cpbi.102}
}

@software{spades_rs2026,
  author    = {Rout, Sagnik and {Gemini 3.8 Flash (Google DeepMind)}},
  title     = {{spades-rs: A Rust-based de novo genome assembler designed for low-memory environments}},
  url       = {https://github.com/sagnikrout/spades-rs},
  version   = {1.0.0},
  year      = {2026}
}

@misc{gemini38flash,
  author    = {{Google DeepMind}},
  title     = {{Gemini 3.8 Flash}},
  year      = {2026},
  url       = {https://deepmind.google/technologies/gemini/},
  note      = {Interactive LLM assistance for code translation, test authoring, and documentation auditing}
}
```

## Authors and contributors

* **Sagnik Rout** (Lead developer, architecture and algorithms)
* **Gemini 3.8 Flash** (Google DeepMind; LLM co-author for code translation, test authoring, and documentation auditing)

## License
MIT License. Free for academic, non-commercial, and commercial genomics workflows.
