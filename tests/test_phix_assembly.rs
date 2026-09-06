//! End-to-end Integration Test for Real Public Sequencing Data: Bacteriophage PhiX174 (NC_001422.1).
//!
//! Validates:
//! 1. Ingestion and assembly of genuine Illumina HiSeq 2000 paired-end reads (SRR2057028).
//! 2. Circular genome traversal in the compacted de Bruijn graph.
//! 3. Parallel error bubble popping and automatic noise valley thresholding.
//! 4. Generation of a >= 5.0 kb single contig matching the NCBI reference with >= 99.8% identity.

use spades_rs::assemble::{run_assembly, AssemblerConfig};
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;

#[test]
fn test_phix174_wgs_public_assembly() {
    let r1 = PathBuf::from("data/phix/phix_std_1.fq.gz");
    let r2 = PathBuf::from("data/phix/phix_std_2.fq.gz");
    let ref_path = PathBuf::from("data/phix/phix174_ref.fa");

    if !r1.exists() || !r2.exists() || !ref_path.exists() {
        eprintln!("PhiX174 dataset not present, skipping public benchmark test.");
        return;
    }

    let config = AssemblerConfig {
        k: 31,
        min_coverage: 5.0,
        min_contig_len: 200,
        bloom_bits: 16 * 1024 * 1024,
        error_correct: false,
        ..AssemblerConfig::default()
    };

    let result = run_assembly(&[r1, r2], &config).expect("PhiX assembly failed");

    println!("Assembled contigs: {}", result.stats.total_contigs);
    println!("Max contig length: {} bp", result.stats.max_contig_length);
    println!("N50: {} bp", result.stats.n50);

    // Verify contiguity: Primary contig must represent the entire circular genome (>= 5,300 bp)
    assert!(
        result.stats.max_contig_length >= 5300,
        "Expected single contig covering the full ~5.4 kb PhiX genome, got: {} bp",
        result.stats.max_contig_length
    );

    // Verify sequence accuracy against NCBI NC_001422.1
    let mut ref_file = File::open(ref_path).expect("Failed to open PhiX reference");
    let mut ref_str = String::new();
    ref_file.read_to_string(&mut ref_str).unwrap();

    let ref_seq: Vec<u8> = ref_str
        .lines()
        .filter(|l| !l.starts_with('>'))
        .flat_map(|l| l.bytes())
        .collect();

    let primary_contig = &result.contigs[0].sequence;
    assert!(
        primary_contig.len() >= ref_seq.len(),
        "Primary contig length {} is shorter than reference genome {}",
        primary_contig.len(),
        ref_seq.len()
    );

    // Compute circular alignment identity
    let mut ref_doubled = ref_seq.clone();
    ref_doubled.extend_from_slice(&ref_seq);

    let rc_ref: Vec<u8> = ref_seq
        .iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect();
    let mut rc_ref_doubled = rc_ref.clone();
    rc_ref_doubled.extend_from_slice(&rc_ref);

    let probe = &primary_contig[..50];
    let find_sub = |haystack: &[u8], needle: &[u8]| -> Option<usize> {
        haystack
            .windows(needle.len())
            .position(|window| window == needle)
    };

    let (offset, is_rc) = if let Some(pos) = find_sub(&ref_doubled, probe) {
        (pos, false)
    } else if let Some(pos) = find_sub(&rc_ref_doubled, probe) {
        (pos, true)
    } else {
        panic!("Probe could not be aligned to reference genome!");
    };

    let aligned_ref = if is_rc {
        &rc_ref_doubled[offset..offset + ref_seq.len()]
    } else {
        &ref_doubled[offset..offset + ref_seq.len()]
    };

    let evaluated_contig = &primary_contig[..ref_seq.len()];
    let matches = evaluated_contig
        .iter()
        .zip(aligned_ref.iter())
        .filter(|(&a, &b)| a == b)
        .count();

    let identity = (matches as f64 / ref_seq.len() as f64) * 100.0;
    println!("Matching bases: {} / {}", matches, ref_seq.len());
    println!("Sequence identity: {:.4}%", identity);

    assert!(
        identity >= 99.8,
        "Assembly identity must be >= 99.8%, got: {:.4}%",
        identity
    );
}
