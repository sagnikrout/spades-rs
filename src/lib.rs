#[cfg(not(target_pointer_width = "64"))]
compile_error!("spades-rs requires a 64-bit architecture (x86_64, aarch64, riscv64). 32-bit systems are not supported due to the 4 GB address space barrier and 64-bit atomic requirements for de novo genome assembly.");

pub mod assemble;
pub mod bloom;
pub mod dna;
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
pub mod scaffold;
pub mod simplify;
pub mod spaligner;
pub mod splitter;

// Re-exports for convenience and backward compatibility
pub use modes::{apply_meta_filter, PlasmidDetector, RnaEngine, SingleCellNormalizer};
