//! `spades-rs`: De novo genome, metagenome, plasmid, single-cell, and transcriptome assembler.
//!
//! # Architecture Overview
//! - [`dna`]: 2-bit DNA base encoding, hardware bit-reversal reverse-complement, and [`Kmer256`] up to $k = 127$.
//! - [`packed_reads`]: 2-bit packed read storage in memory with paired-end synchronization and ambiguous (`N`) base handling.
//! - [`bloom`]: Lock-free two-tier atomic Bloom filter ([`TwoTierFilter`]) for filtering singleton sequencing errors.
//! - [`hammer`]: BayesHammer-style single-base error correction against the solid $k$-mer spectrum.
//! - [`graph`]: Bidirected compacted de Bruijn graph (cDBG) construction and maximal unitig extraction.
//! - [`simplify`]: Topological graph simplification (tip clipping, bubble popping, and unitig stitching).
//! - [`paired_info`] & [`expander`]: Paired-end insert size estimation and repeat resolution (ExSPAnder).
//! - [`spaligner`] & [`splitter`]: Hybrid long-read (PacBio/ONT) and barcoded linked-read (10x/TELL-Seq) repeat bridging.
//! - [`scaffold`] & [`polisher`]: Paired-end scaffolding, local micro-assembly gap closure, and read-alignment consensus polishing.
//! - [`modes`]: Specialized filtering and normalization for metagenomic, plasmid, RNA-Seq, and single-cell assemblies.
//! - [`qc`], [`eval`], & [`gfa`]: Read quality control/trimming, QUAST-style assembly evaluation, and GFA v1.1 graph export.

#![warn(missing_docs)]

#[cfg(not(target_pointer_width = "64"))]
compile_error!("spades-rs requires a 64-bit architecture (x86_64, aarch64, riscv64). 32-bit systems are not supported due to the 4 GB address space barrier and 64-bit atomic requirements for de novo genome assembly.");

pub mod assemble;
pub mod bloom;
pub mod dna;
pub mod eval;
pub mod expander;
pub mod fastq;
pub mod gfa;
pub mod graph;
pub mod hammer;
pub mod memory;
pub mod modes;
pub mod multik;
pub mod packed_reads;
pub mod paired_info;
pub mod polisher;
pub mod qc;
pub mod scaffold;
pub mod simplify;
pub mod spaligner;
pub mod splitter;

// Re-exports for convenience and backward compatibility
pub use assemble::{run_assembly, AssemblerConfig, AssemblyResult, AssemblyStats};
pub use bloom::TwoTierFilter;
pub use dna::Kmer256;
pub use graph::Unitig;
pub use modes::{apply_meta_filter, PlasmidDetector, RnaEngine, SingleCellNormalizer};
pub use multik::{run_multik_assembly, MultiKConfig};
pub use packed_reads::PackedReads;
