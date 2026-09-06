use spades_rs::dna::{
    base_to_2bit, bit2_to_base, canonical_kmer_u64, kmer_to_string, revcomp_bytes,
    revcomp_kmer_u64, string_to_kmer, Kmer256,
};

#[test]
fn test_base_encoding_and_invalid_chars() {
    // Standard nucleotides
    assert_eq!(base_to_2bit(b'A'), Some(0));
    assert_eq!(base_to_2bit(b'a'), Some(0));
    assert_eq!(base_to_2bit(b'C'), Some(1));
    assert_eq!(base_to_2bit(b'c'), Some(1));
    assert_eq!(base_to_2bit(b'G'), Some(2));
    assert_eq!(base_to_2bit(b'g'), Some(2));
    assert_eq!(base_to_2bit(b'T'), Some(3));
    assert_eq!(base_to_2bit(b't'), Some(3));

    // Decode check
    assert_eq!(bit2_to_base(0), b'A');
    assert_eq!(bit2_to_base(1), b'C');
    assert_eq!(bit2_to_base(2), b'G');
    assert_eq!(bit2_to_base(3), b'T');

    // Invalid nucleotides must return None in string_to_kmer
    assert!(string_to_kmer(b"ACGTN", 5).is_none());
    assert!(string_to_kmer(b"ACGT1", 5).is_none());
    assert!(string_to_kmer(b"AC GT", 5).is_none());
    assert!(string_to_kmer(b"AC-GT", 5).is_none());
    assert!(string_to_kmer(b"ACGT", 5).is_none()); // string shorter than k
}

#[test]
fn test_revcomp_bytes_properties() {
    let original = b"ACGTACGTNNACGT";
    let rc = revcomp_bytes(original);
    let rc_rc = revcomp_bytes(&rc);
    assert_eq!(original.to_vec(), rc_rc);

    assert_eq!(revcomp_bytes(b"AAAA"), b"TTTT");
    assert_eq!(revcomp_bytes(b"CCCC"), b"GGGG");
    assert_eq!(revcomp_bytes(b"ACGT"), b"ACGT"); // Palindrome
    assert_eq!(revcomp_bytes(b"CGCG"), b"CGCG"); // Palindrome
}

#[test]
fn test_revcomp_u64_invariance() {
    // For all k from 1 to 32, revcomp(revcomp(kmer)) == kmer
    for k in 1..=32 {
        let seq: Vec<u8> = (0..k)
            .map(|i| match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();
        let kmer = string_to_kmer(&seq, k).expect("valid kmer");
        let rc = revcomp_kmer_u64(kmer, k);
        let rc_rc = revcomp_kmer_u64(rc, k);
        assert_eq!(kmer, rc_rc, "Failed double revcomp invariance for k={}", k);

        // String roundtrip
        let decoded = kmer_to_string(kmer, k);
        assert_eq!(decoded.as_bytes(), &seq[..], "Failed string decode for k={}", k);
    }
}

#[test]
fn test_canonical_kmer_u64_idempotence() {
    for k in [1, 7, 15, 21, 31, 32] {
        let seq: Vec<u8> = (0..k)
            .map(|i| match (i * 7) % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();
        let kmer = string_to_kmer(&seq, k).unwrap();
        let (can1, is_rc1) = canonical_kmer_u64(kmer, k);
        let (can2, _) = canonical_kmer_u64(can1, k);
        assert_eq!(can1, can2, "Canonical representation must be idempotent for k={}", k);

        let rc = revcomp_kmer_u64(kmer, k);
        let (can_from_rc, is_rc2) = canonical_kmer_u64(rc, k);
        assert_eq!(can1, can_from_rc, "Canonical must be identical from forward or RC for k={}", k);
        if kmer != rc {
            assert_ne!(is_rc1, is_rc2, "One must be RC and one forward");
        }
    }
}

#[test]
fn test_kmer256_large_k_invariance() {
    // Test large k-mers spanning across 64-bit boundaries (k = 33, 55, 64, 77, 99, 111, 127, 128)
    for &k in &[33, 55, 64, 77, 99, 111, 127, 128] {
        let seq: Vec<u8> = (0..k)
            .map(|i| match (i * 13 + 5) % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();

        let kmer = Kmer256::from_bytes(&seq, k).expect("valid Kmer256");
        let rc = kmer.revcomp(k);
        let rc_rc = rc.revcomp(k);
        assert_eq!(kmer, rc_rc, "Double revcomp failed for Kmer256 with k={}", k);

        let (can, _) = kmer.canonical(k);
        let (can_rc, _) = rc.canonical(k);
        assert_eq!(can, can_rc, "Canonical mismatch between fwd and rc in Kmer256 for k={}", k);

        let decoded = kmer.to_string(k);
        assert_eq!(decoded.as_bytes(), seq.as_slice(), "String round-trip failed in Kmer256 for k={}", k);
    }
}

#[test]
fn test_kmer256_extend_and_prepend_exactness() {
    let k = 63;
    let base_seq: Vec<u8> = (0..k)
        .map(|i| match i % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        })
        .collect();

    let kmer = Kmer256::from_bytes(&base_seq, k).unwrap();

    // Extend right with 'T' (base 3)
    let next_base = b'T';
    let extended = kmer.extend_right(base_to_2bit(next_base).unwrap(), k);
    let mut expected_seq = base_seq[1..].to_vec();
    expected_seq.push(next_base);
    assert_eq!(extended.to_string(k).as_bytes(), expected_seq.as_slice());

    // Prepend left with 'G' (base 2)
    let prev_base = b'G';
    let prepended = kmer.prepend_left(base_to_2bit(prev_base).unwrap(), k);
    let mut expected_prep = vec![prev_base];
    expected_prep.extend_from_slice(&base_seq[..k - 1]);
    assert_eq!(prepended.to_string(k).as_bytes(), expected_prep.as_slice());
}
