use spades_rs::assemble::{run_assembly, AssemblerConfig};
use spades_rs::dna::{canonical_kmer_u64, revcomp_kmer_u64, string_to_kmer};
use spades_rs::hammer::hamming_distance_2bit;
use std::path::PathBuf;

#[test]
fn test_dna_primitives() {
    let kmer = string_to_kmer(b"ACGTACGT", 8).unwrap();
    let rc = revcomp_kmer_u64(kmer, 8);
    let (can, _) = canonical_kmer_u64(kmer, 8);
    assert_eq!(can, kmer.min(rc));
}

#[test]
fn test_hamming_distance() {
    // Exact match
    let k1 = string_to_kmer(b"ACGT", 4).unwrap();
    let k2 = string_to_kmer(b"ACGT", 4).unwrap();
    assert_eq!(hamming_distance_2bit(k1, k2), 0);

    // 1 base mutation (A -> T)
    let k3 = string_to_kmer(b"TCGT", 4).unwrap();
    assert_eq!(hamming_distance_2bit(k1, k3), 1);

    // 2 base mutations
    let k4 = string_to_kmer(b"TGGT", 4).unwrap();
    assert_eq!(hamming_distance_2bit(k1, k4), 2);
}

#[test]
fn test_assembly_accuracy_100_percent() {
    let r1 = PathBuf::from("data/ecoli_1K_1.fq.gz");
    let r2 = PathBuf::from("data/ecoli_1K_2.fq.gz");
    let ref_path = PathBuf::from("data/reference_1K.fa.gz");

    if !r1.exists() || !r2.exists() || !ref_path.exists() {
        eprintln!("Test dataset not present, skipping integration test");
        return;
    }

    let config = AssemblerConfig {
        k: 31,
        min_coverage: 5.0,
        min_contig_len: 200,
        bloom_bits: 16 * 1024 * 1024,
        error_correct: true,
        ..AssemblerConfig::default()
    };

    let result = run_assembly(&[r1, r2], &config).expect("Assembly failed");
    assert_eq!(
        result.stats.total_contigs, 1,
        "Expected exactly 1 assembled contig"
    );
    assert_eq!(
        result.stats.max_contig_length, 1000,
        "Expected exactly 1,000 bp assembled length"
    );

    // Verify 100% identity to reference
    use flate2::read::MultiGzDecoder;
    use std::fs::File;
    use std::io::Read;

    let mut gz = MultiGzDecoder::new(File::open(ref_path).unwrap());
    let mut ref_str = String::new();
    gz.read_to_string(&mut ref_str).unwrap();
    let ref_seq: Vec<u8> = ref_str
        .lines()
        .filter(|l| !l.starts_with('>'))
        .flat_map(|l| l.bytes())
        .collect();

    let seq = &result.contigs[0].sequence;
    let rc_seq: Vec<u8> = seq
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

    let is_match = (seq == &ref_seq) || (rc_seq == ref_seq);
    assert!(is_match, "Assembled sequence must match reference with 100% identity in forward or reverse orientation!");
}
