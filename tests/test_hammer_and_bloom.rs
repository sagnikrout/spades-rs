use spades_rs::bloom::TwoTierFilter;
use spades_rs::dna::{base_to_2bit, canonical_kmer_u64, string_to_kmer};
use spades_rs::hammer::{hamming_distance_2bit, ErrorCorrector};

#[test]
fn test_hamming_distance_exhaustive_single_base() {
    let bases = [b'A', b'C', b'G', b'T'];
    for (i, &b1) in bases.iter().enumerate() {
        for (j, &b2) in bases.iter().enumerate() {
            let val1 = base_to_2bit(b1).unwrap() as u64;
            let val2 = base_to_2bit(b2).unwrap() as u64;
            let dist = hamming_distance_2bit(val1, val2);
            if i == j {
                assert_eq!(dist, 0, "Same base {} must have distance 0", b1 as char);
            } else {
                assert_eq!(dist, 1, "Different bases {} and {} must have distance 1", b1 as char, b2 as char);
            }
        }
    }
}

#[test]
fn test_hamming_distance_multi_base() {
    let s1 = b"ACGTACGTACGTACGT";
    let s2 = b"ACCTACATACGAACGT"; // mismatches at index 2 (G->C), 6 (G->A), 11 (T->A) = 3 mismatches

    let k1 = string_to_kmer(s1, 16).unwrap();
    let k2 = string_to_kmer(s2, 16).unwrap();

    assert_eq!(hamming_distance_2bit(k1, k2), 3);
}

#[test]
fn test_two_tier_bloom_transitions() {
    let filter = TwoTierFilter::new(1024 * 1024);

    let kmer_a = 0x123456789abcdef0u64;
    let kmer_b = 0xfedcba9876543210u64;
    let kmer_c = 0xa5a5a5a5a5a5a5a5u64;

    // Initially none are solid
    assert!(!filter.is_solid(kmer_a));
    assert!(!filter.is_solid(kmer_b));
    assert!(!filter.is_solid(kmer_c));

    // First insertion of A: should return false (not solid yet)
    let became_solid_1 = filter.insert(kmer_a);
    assert!(!became_solid_1);
    assert!(!filter.is_solid(kmer_a));

    // Second insertion of A: should return true (promoted to solid)
    let became_solid_2 = filter.insert(kmer_a);
    assert!(became_solid_2);
    assert!(filter.is_solid(kmer_a));

    // Third insertion of A: still solid
    let became_solid_3 = filter.insert(kmer_a);
    assert!(became_solid_3);
    assert!(filter.is_solid(kmer_a));

    // B and C should still not be solid
    assert!(!filter.is_solid(kmer_b));
    assert!(!filter.is_solid(kmer_c));
}

#[test]
fn test_error_corrector_single_base_correction() {
    let k = 15;
    let filter = TwoTierFilter::new(1024 * 1024);

    // Non-repetitive reference sequence: 41 bp
    let ref_seq = b"GATTACATAGGATACAGGATACCAGATTAGACATAGACAG";

    // Populate filter with solid k-mers from reference (insert twice)
    for i in 0..=(ref_seq.len() - k) {
        let km = string_to_kmer(&ref_seq[i..i + k], k).unwrap();
        let (can, _) = canonical_kmer_u64(km, k);
        filter.insert(can);
        filter.insert(can);
    }

    // Corrupt one base at position 20: 'A' -> 'T'
    let mut corrupted = ref_seq.to_vec();
    assert_eq!(corrupted[20], b'A');
    corrupted[20] = b'T';

    let corrector = ErrorCorrector::new(k);
    let (corrected_reads, count) = corrector.correct_reads(vec![corrupted], &filter);

    assert_eq!(count, 1, "Expected exactly 1 read corrected");
    assert_eq!(corrected_reads[0], ref_seq.to_vec(), "Corrupted read should be restored to reference");
}
