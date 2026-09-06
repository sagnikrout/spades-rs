use spades_rs::graph::Unitig;
use spades_rs::simplify::Simplifier;

#[test]
fn test_length_aware_tip_clipping() {
    let k = 21;
    let simplifier = Simplifier::new(k, 5.0, 50);
    let dynamic_tip_threshold = 10.0;

    let tip_short_low_cov = Unitig {
        id: 0,
        sequence: vec![b'A'; 30], // len 30 < 2*k (42)
        mean_coverage: 4.0,       // < dynamic_tip_threshold (10.0)
        kmers_count: 10,
    };

    let tip_long_low_cov = Unitig {
        id: 1,
        sequence: vec![b'C'; 100], // len 100 >= 2*k (42)
        mean_coverage: 4.0,        // < dynamic_tip_threshold (10.0)
        kmers_count: 80,
    };

    let unitig_short_high_cov = Unitig {
        id: 2,
        sequence: vec![b'G'; 30], // len 30 < 2*k
        mean_coverage: 25.0,      // > dynamic_tip_threshold
        kmers_count: 10,
    };

    let unitig_noise = Unitig {
        id: 3,
        sequence: vec![b'T'; 30],
        mean_coverage: 1.0, // < 2.5 (min_coverage * 0.5)
        kmers_count: 10,
    };

    let input = vec![
        tip_short_low_cov,
        tip_long_low_cov,
        unitig_short_high_cov,
        unitig_noise,
    ];

    let clipped = simplifier.clip_tips(input, dynamic_tip_threshold);
    let kept_ids: Vec<usize> = clipped.iter().map(|u| u.id).collect();

    // id 0 (short low cov) clipped
    // id 3 (noise) clipped
    // id 1 (long low cov) PRESERVED by length-aware guard
    // id 2 (short high cov) PRESERVED
    assert_eq!(kept_ids, vec![1, 2]);
}

#[test]
fn test_bubble_popping_snp_and_indels() {
    let k = 5;
    let simplifier = Simplifier::new(k, 5.0, 10);

    // Endpoints: prefix "ACGT", suffix "TGCA"
    // Allele 1: high coverage true path
    let major_allele = Unitig {
        id: 0,
        sequence: b"ACGTAAAATGCA".to_vec(),
        mean_coverage: 50.0,
        kmers_count: 8,
    };

    // Allele 2: low coverage sequencing error bubble
    let minor_allele = Unitig {
        id: 1,
        sequence: b"ACGTTTTATGCA".to_vec(),
        mean_coverage: 2.0,
        kmers_count: 8,
    };

    // Independent unitig with different endpoints
    let other_unitig = Unitig {
        id: 2,
        sequence: b"GGGGCCCCAAAA".to_vec(),
        mean_coverage: 20.0,
        kmers_count: 8,
    };

    let input = vec![major_allele, minor_allele, other_unitig];
    let popped = simplifier.pop_bubbles(input);

    assert_eq!(popped.len(), 2);
    let kept_ids: Vec<usize> = popped.iter().map(|u| u.id).collect();
    assert!(kept_ids.contains(&0));
    assert!(kept_ids.contains(&2));
    assert!(!kept_ids.contains(&1), "Minor allele bubble must be popped");
}

#[test]
fn test_stitch_unitigs_linear_chain() {
    let k = 5;
    let simplifier = Simplifier::new(k, 5.0, 10);

    // Three unitigs in series sharing 4-bp overlaps:
    // u0: TTTAACGT -> suffix "ACGT"
    // u1: ACGTCCAT -> prefix "ACGT", suffix "CCAT"
    // u2: CCATGAAA -> prefix "CCAT"
    let u0 = Unitig {
        id: 0,
        sequence: b"TTTAACGT".to_vec(),
        mean_coverage: 20.0,
        kmers_count: 4,
    };
    let u1 = Unitig {
        id: 1,
        sequence: b"ACGTCCAT".to_vec(),
        mean_coverage: 20.0,
        kmers_count: 4,
    };
    let u2 = Unitig {
        id: 2,
        sequence: b"CCATGAAA".to_vec(),
        mean_coverage: 20.0,
        kmers_count: 4,
    };

    let stitched = simplifier.stitch_unitigs(vec![u0, u1, u2]);
    assert_eq!(stitched.len(), 1, "Linear chain should stitch into exactly 1 unitig");
    assert_eq!(stitched[0].sequence, b"TTTAACGTCCATGAAA");
}

#[test]
fn test_simplify_end_to_end() {
    let k = 5;
    let simplifier = Simplifier::new(k, 5.0, 10);

    // Linear chain with a bubble in the middle and a short dead-end tip
    let u_start = Unitig {
        id: 0,
        sequence: b"CCCCACGT".to_vec(),
        mean_coverage: 30.0,
        kmers_count: 4,
    };
    let u_bubble_major = Unitig {
        id: 1,
        sequence: b"ACGTAAAATGCA".to_vec(),
        mean_coverage: 30.0,
        kmers_count: 8,
    };
    let u_bubble_minor = Unitig {
        id: 2,
        sequence: b"ACGTTTTATGCA".to_vec(),
        mean_coverage: 1.5,
        kmers_count: 8,
    };
    let u_end = Unitig {
        id: 3,
        sequence: b"TGCAGGGG".to_vec(),
        mean_coverage: 30.0,
        kmers_count: 4,
    };
    let u_tip = Unitig {
        id: 4,
        sequence: b"ACGTTA".to_vec(), // len 6 < 2*k (10)
        mean_coverage: 1.0,
        kmers_count: 2,
    };

    let result = simplifier.simplify(vec![u_start, u_bubble_major, u_bubble_minor, u_end, u_tip]);

    // Should pop bubble minor, clip tip, and stitch u_start -> u_bubble_major -> u_end
    assert_eq!(result.len(), 1);
    assert_eq!(result[0].sequence, b"CCCCACGTAAAATGCAGGGG");
}
