//! Integration tests for specialized assembly modes (plasmid, metagenomic, RNA, single-cell), scaffolding, and polishing.

use hashbrown::HashMap;
use spades_rs::graph::Unitig;
use spades_rs::modes::{apply_meta_filter, PlasmidDetector, RnaEngine, SingleCellNormalizer};
use spades_rs::paired_info::{PairedInfoIndex, PairedLink};
use spades_rs::polisher::Polisher;
use spades_rs::scaffold::Scaffolder;

#[test]
fn test_plasmid_detector() {
    let detector = PlasmidDetector {
        k: 21,
        min_plasmid_len: 20,
        copy_number_threshold: 2.0,
    };

    // Circular contig (prefix == suffix)
    let seq_circular = b"ACGTACGTACGTACGTACGTACGTA".to_vec(); // length 25 > 20
    let mut circ = seq_circular.clone();
    circ.extend_from_slice(&seq_circular[..20]); // ensure prefix == suffix of length 20

    // High copy but short (len=19 < 20)
    let seq_short = b"ACGTAAAAACGTAAAAACG".to_vec();
    // High copy and long -> Plasmid (make it high copy by giving it 100 cov)
    let mut seq_plasmid = vec![b'A'; 200];
    seq_plasmid[0] = b'T'; // break circularity

    // Chromosomal (low copy)
    let mut seq_chrom = vec![b'C'; 200];
    seq_chrom[0] = b'G'; // break circularity

    let unitigs = vec![
        Unitig {
            id: 0,
            sequence: circ,
            mean_coverage: 10.0,
            kmers_count: 50,
        }, // Circular -> Plasmid
        Unitig {
            id: 1,
            sequence: seq_short,
            mean_coverage: 100.0,
            kmers_count: 50,
        }, // High copy but len 19 < min_plasmid_len (20) -> Chromosomal
        Unitig {
            id: 2,
            sequence: seq_plasmid,
            mean_coverage: 300.0,
            kmers_count: 180,
        }, // High copy and long -> Plasmid
        Unitig {
            id: 3,
            sequence: seq_chrom,
            mean_coverage: 1.0,
            kmers_count: 180,
        }, // Chromosomal
    ];

    let (chromosomal, plasmids) = detector.extract_plasmids(unitigs);

    assert_eq!(plasmids.len(), 2);
    assert_eq!(chromosomal.len(), 2);

    let plasmid_ids: Vec<_> = plasmids.iter().map(|u| u.id).collect();
    assert!(plasmid_ids.contains(&0));
    assert!(plasmid_ids.contains(&2));
}

#[test]
fn test_meta_filter() {
    let unitigs = vec![
        Unitig {
            id: 0,
            sequence: vec![b'A'; 1000],
            mean_coverage: 1.5,
            kmers_count: 980,
        }, // Filtered (cov < 2.0)
        Unitig {
            id: 1,
            sequence: vec![b'C'; 1000],
            mean_coverage: 3.0,
            kmers_count: 980,
        }, // Kept
        Unitig {
            id: 2,
            sequence: vec![b'G'; 100],
            mean_coverage: 10.0,
            kmers_count: 80,
        }, // Filtered (len < 500)
    ];
    let filtered = apply_meta_filter(unitigs, 500);
    assert_eq!(filtered.len(), 1);
    assert_eq!(filtered[0].id, 1);
}

#[test]
fn test_scaffolder_basic() {
    let scaffolder = Scaffolder {
        min_links: 2,
        default_gap_len: 5,
    };

    let contigs = vec![
        Unitig {
            id: 0,
            sequence: b"ACGTACGT".to_vec(),
            mean_coverage: 10.0,
            kmers_count: 8,
        },
        Unitig {
            id: 1,
            sequence: b"TGCATGCA".to_vec(),
            mean_coverage: 10.0,
            kmers_count: 8,
        },
    ];

    let mut links = HashMap::new();
    let mut link_map = HashMap::new();
    link_map.insert(
        1,
        PairedLink {
            target_unitig: 1,
            support_count: 3,
            total_distance: 5.0,
        },
    );
    links.insert(0, link_map);

    let paired_info = PairedInfoIndex {
        mean_insert_size: 100.0,
        insert_size_stdev: 10.0,
        links,
    };

    let scaffolds = scaffolder.build_scaffolds(&contigs, &paired_info);

    assert_eq!(scaffolds.len(), 1);
    let combined_len = 8 + 5 + 8; // len(contig0) + gap + len(contig1)
    assert_eq!(scaffolds[0].sequence.len(), combined_len);
    assert_eq!(&scaffolds[0].sequence[0..8], b"ACGTACGT");
    assert_eq!(&scaffolds[0].sequence[8..13], b"NNNNN");
    assert_eq!(&scaffolds[0].sequence[13..21], b"TGCATGCA");
}

#[test]
fn test_rna_engine() {
    let engine = RnaEngine {
        min_isoform_len: 50,
        min_isoform_coverage: 5.0,
    };
    let unitigs = vec![
        Unitig {
            id: 0,
            sequence: vec![b'A'; 100],
            mean_coverage: 50.0,
            kmers_count: 80,
        }, // Kept
        Unitig {
            id: 1,
            sequence: vec![b'C'; 40],
            mean_coverage: 50.0,
            kmers_count: 20,
        }, // Filtered (len)
        Unitig {
            id: 2,
            sequence: vec![b'G'; 100],
            mean_coverage: 2.0,
            kmers_count: 80,
        }, // Filtered (cov)
    ];
    let transcripts = engine.process_transcripts(unitigs);
    assert_eq!(transcripts.len(), 1);
    assert_eq!(transcripts[0].id, 0);
}

#[test]
fn test_single_cell_normalizer() {
    let normalizer = SingleCellNormalizer {
        min_coverage_cutoff: 5.0,
    };
    let unitigs = vec![
        Unitig {
            id: 0,
            sequence: b"AAAA".to_vec(),
            mean_coverage: 1000.0,
            kmers_count: 4,
        }, // Capped
        Unitig {
            id: 1,
            sequence: b"CCCC".to_vec(),
            mean_coverage: 2.0,
            kmers_count: 4,
        }, // Filtered
        Unitig {
            id: 2,
            sequence: b"GGGG".to_vec(),
            mean_coverage: 10.0,
            kmers_count: 4,
        }, // Unchanged
    ];
    let normalized = normalizer.normalize_coverage(unitigs);
    assert_eq!(normalized.len(), 2);

    let u0 = normalized.iter().find(|u| u.id == 0).unwrap();
    assert_eq!(u0.mean_coverage, 500.0); // SC mode caps at 500.0

    let u2 = normalized.iter().find(|u| u.id == 2).unwrap();
    assert_eq!(u2.mean_coverage, 10.0);
}

#[test]
fn test_polisher() {
    let polisher = Polisher {
        k: 3,
        min_coverage_support: 1,
        careful: false,
    };

    // Contig with an error 'T' instead of 'C' at index 3
    let contigs = vec![Unitig {
        id: 0,
        sequence: b"ACGTTA".to_vec(),
        mean_coverage: 10.0,
        kmers_count: 6,
    }];

    let reads = vec![b"ACGTCA".to_vec(), b"ACGTCA".to_vec(), b"ACGTCA".to_vec()];

    let polished = polisher.polish_contigs(contigs, &reads);

    assert_eq!(polished.0.len(), 1);
    assert_eq!(polished.0[0].sequence, b"ACGTCA");
}

#[test]
fn test_spaligner_hybrid() {
    use spades_rs::spaligner::LongReadResolver;

    let resolver = LongReadResolver {
        k: 5,
        min_seed_matches: 1,
    };

    // Two unitigs that share a 4-bp overlap (k-1 = 4)
    // u0: GGCTACGT (8 bp) -> ends with ACGT
    // u1: ACGTCCTA (8 bp) -> starts with ACGT
    let u0 = Unitig {
        id: 0,
        sequence: b"GGCTACGT".to_vec(),
        mean_coverage: 20.0,
        kmers_count: 4,
    };
    let u1 = Unitig {
        id: 1,
        sequence: b"ACGTCCTA".to_vec(),
        mean_coverage: 20.0,
        kmers_count: 4,
    };

    // Long read spanning both: "GGCTACGT" + "CCTA" = "GGCTACGTCCTA"
    let long_reads = vec![b"GGCTACGTCCTA".to_vec()];

    let bridged = resolver.bridge_with_long_reads(vec![u0, u1], &long_reads);
    assert_eq!(bridged.len(), 1);
    assert_eq!(bridged[0].sequence, b"GGCTACGTCCTA");
}

#[test]
fn test_local_gap_closer() {
    use spades_rs::scaffold::LocalGapCloser;

    let closer = LocalGapCloser {
        k: 5,
        max_gap_len: 100,
    };

    // Test direct overlap
    let left = b"ACGTACGTACGT";
    let right = b"ACGTGGGGGG";
    let closed = closer.close_gap(left, right, &[]);
    assert!(closed.is_some());
    assert_eq!(closed.unwrap(), b"ACGTACGTACGTGGGGGG");

    // Test k-mer walk across gap
    let left2 = b"ACGTACGTAG";
    let right2 = b"TGCATGCATG";
    let read = b"ACGTACGTAGAATTGCATGCATG".to_vec();
    let read_slice: &[u8] = &read;
    let closed2 = closer.close_gap(left2, right2, &[read_slice, read_slice]);
    assert!(closed2.is_some());
    assert_eq!(closed2.unwrap(), b"ACGTACGTAGAATTGCATGCATG");
}

#[test]
fn test_diploid_phaser() {
    use spades_rs::modes::DiploidPhaser;

    let phaser = DiploidPhaser::default();
    let k = 5;

    // Heterozygous bubble alleles with shared (k-1) prefix ACGT and suffix TGCA
    let allele_a = Unitig {
        id: 0,
        sequence: b"ACGTAAAATGCA".to_vec(),
        mean_coverage: 32.0,
        kmers_count: 8,
    };
    let allele_b = Unitig {
        id: 1,
        sequence: b"ACGTAATATGCA".to_vec(),
        mean_coverage: 28.0,
        kmers_count: 8,
    };
    let homozygous_flank = Unitig {
        id: 2,
        sequence: b"CCCCCCCCCCCC".to_vec(),
        mean_coverage: 60.0,
        kmers_count: 8,
    };

    let (primary, haplotigs) =
        phaser.phase_haplotigs(k, vec![allele_a, allele_b, homozygous_flank]);

    assert_eq!(primary.len(), 2);
    assert_eq!(haplotigs.len(), 1);
    // allele_a has higher coverage, so allele_b is phased into haplotigs
    assert_eq!(primary[0].sequence, b"ACGTAAAATGCA");
    assert_eq!(primary[1].sequence, b"CCCCCCCCCCCC");
    assert_eq!(haplotigs[0].sequence, b"ACGTAATATGCA");
}

#[test]
fn test_diploid_phaser_identical_and_revcomp_not_phased() {
    use spades_rs::modes::DiploidPhaser;

    let phaser = DiploidPhaser::default();
    let k = 5;

    // Identical sequences (homozygous locus): must NOT be split
    let u0 = Unitig {
        id: 0,
        sequence: b"ACGTAAAATGCA".to_vec(),
        mean_coverage: 30.0,
        kmers_count: 8,
    };
    let u1 = Unitig {
        id: 1,
        sequence: b"ACGTAAAATGCA".to_vec(),
        mean_coverage: 30.0,
        kmers_count: 8,
    };

    let (primary, haplotigs) = phaser.phase_haplotigs(k, vec![u0, u1]);
    assert_eq!(primary.len(), 2);
    assert!(
        haplotigs.is_empty(),
        "Identical homozygous sequences should not be phased into haplotigs"
    );

    // Reverse-complement identical sequences (opposite strand homozygous duplicate): must NOT be split
    let rc_seq = spades_rs::dna::revcomp_bytes(b"ACGTAAAATGCA");
    let u2 = Unitig {
        id: 2,
        sequence: b"ACGTAAAATGCA".to_vec(),
        mean_coverage: 30.0,
        kmers_count: 8,
    };
    let u3 = Unitig {
        id: 3,
        sequence: rc_seq,
        mean_coverage: 30.0,
        kmers_count: 8,
    };

    let (primary2, haplotigs2) = phaser.phase_haplotigs(k, vec![u2, u3]);
    assert_eq!(primary2.len(), 2);
    assert!(
        haplotigs2.is_empty(),
        "Reverse-complement duplicate should not be phased into haplotigs"
    );
}

#[test]
fn test_diploid_phaser_coverage_and_length_ratio_thresholds() {
    use spades_rs::modes::DiploidPhaser;

    let phaser = DiploidPhaser::default();
    let k = 5;

    // Asymmetric coverage: 100x vs 10x (ratio = 10.0 > 2.5 cutoff) -> sequencing artifact, not diploid allele
    let u0 = Unitig {
        id: 0,
        sequence: b"ACGTAAAATGCA".to_vec(),
        mean_coverage: 100.0,
        kmers_count: 8,
    };
    let u1 = Unitig {
        id: 1,
        sequence: b"ACGTAATATGCA".to_vec(),
        mean_coverage: 10.0,
        kmers_count: 8,
    };

    let (primary, haplotigs) = phaser.phase_haplotigs(k, vec![u0, u1]);
    assert_eq!(primary.len(), 2);
    assert!(
        haplotigs.is_empty(),
        "Coverage asymmetry > 2.5x must not be phased as diploid"
    );

    // Length divergence: 100 bp vs 200 bp (ratio = 0.50 > 0.30 cutoff) -> structural variant / repeat, not homologous bubble
    let mut seq_short = b"ACGT".to_vec();
    seq_short.extend_from_slice(&[b'A'; 92]);
    seq_short.extend_from_slice(b"TGCA");

    let mut seq_long = b"ACGT".to_vec();
    seq_long.extend_from_slice(&[b'C'; 192]);
    seq_long.extend_from_slice(b"TGCA");

    let u2 = Unitig {
        id: 2,
        sequence: seq_short,
        mean_coverage: 30.0,
        kmers_count: 8,
    };
    let u3 = Unitig {
        id: 3,
        sequence: seq_long,
        mean_coverage: 30.0,
        kmers_count: 8,
    };

    let (primary2, haplotigs2) = phaser.phase_haplotigs(k, vec![u2, u3]);
    assert_eq!(primary2.len(), 2);
    assert!(
        haplotigs2.is_empty(),
        "Length difference > 30% must not be phased as diploid bubble"
    );
}

#[test]
fn test_diploid_phaser_multibranch_repeat_and_boundaries() {
    use spades_rs::modes::DiploidPhaser;

    let phaser = DiploidPhaser::default();
    let k = 5;

    // Triploid / repeat multi-branch (3 branches sharing same endpoints): must not arbitrarily pair
    let u0 = Unitig {
        id: 0,
        sequence: b"ACGTAAAATGCA".to_vec(),
        mean_coverage: 30.0,
        kmers_count: 8,
    };
    let u1 = Unitig {
        id: 1,
        sequence: b"ACGTAATATGCA".to_vec(),
        mean_coverage: 30.0,
        kmers_count: 8,
    };
    let u2 = Unitig {
        id: 2,
        sequence: b"ACGTAACATGCA".to_vec(),
        mean_coverage: 30.0,
        kmers_count: 8,
    };

    let (primary, haplotigs) = phaser.phase_haplotigs(k, vec![u0, u1, u2]);
    assert_eq!(primary.len(), 3);
    assert!(
        haplotigs.is_empty(),
        "3-way repeat/multibranch locus should not be phased as diploid"
    );

    // Boundary cases: empty unitigs, unitigs shorter than k-1, k=0, k=1
    let (p_empty, h_empty) = phaser.phase_haplotigs(k, Vec::new());
    assert!(p_empty.is_empty());
    assert!(h_empty.is_empty());

    let (p_short, h_short) = phaser.phase_haplotigs(
        31,
        vec![Unitig {
            id: 0,
            sequence: b"ACGT".to_vec(),
            mean_coverage: 10.0,
            kmers_count: 0,
        }],
    );
    assert_eq!(p_short.len(), 1);
    assert!(h_short.is_empty());

    let (p_k0, h_k0) = phaser.phase_haplotigs(
        0,
        vec![Unitig {
            id: 0,
            sequence: b"ACGTACGT".to_vec(),
            mean_coverage: 10.0,
            kmers_count: 0,
        }],
    );
    assert_eq!(p_k0.len(), 1);
    assert!(h_k0.is_empty());
}
