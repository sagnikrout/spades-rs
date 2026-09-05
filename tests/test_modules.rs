use hashbrown::HashMap;
use intelligent_pascal::graph::Unitig;
use intelligent_pascal::modes::{apply_meta_filter, PlasmidDetector};
use intelligent_pascal::paired_info::{PairedInfoIndex, PairedLink};
use intelligent_pascal::polisher::Polisher;
use intelligent_pascal::rna::RnaEngine;
use intelligent_pascal::scaffold::Scaffolder;
use intelligent_pascal::single_cell::SingleCellNormalizer;

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
        }, // High copy but len 24 < min_plasmid_len (20? wait, I changed min to 20, so 24 >= 20. I should make it 19)
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
    use intelligent_pascal::spaligner::LongReadResolver;

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
    use intelligent_pascal::scaffold::LocalGapCloser;

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
