use hashbrown::HashMap;
use spades_rs::expander::ExSPAnder;
use spades_rs::graph::Unitig;
use spades_rs::paired_info::{PairedInfoIndex, PairedLink};
use spades_rs::spaligner::LongReadResolver;

#[test]
fn test_expander_bifurcation_resolution() {
    let k = 5;
    let expander = ExSPAnder {
        min_support: 3,
        confidence_ratio: 2.0,
    };

    // Unitig 0 ends with "ACGT"
    // Unitig 1 starts with "ACGT", ends with "CGCG"
    // Unitig 2 starts with "ACGT", ends with "TTTA"
    let u0 = Unitig {
        id: 0,
        sequence: b"TTTTACGT".to_vec(),
        mean_coverage: 20.0,
        kmers_count: 4,
    };
    let u1 = Unitig {
        id: 1,
        sequence: b"ACGTAAAACGCG".to_vec(),
        mean_coverage: 20.0,
        kmers_count: 8,
    };
    let u2 = Unitig {
        id: 2,
        sequence: b"ACGTTTTATTTA".to_vec(),
        mean_coverage: 20.0,
        kmers_count: 8,
    };

    // Paired-end support heavily favors 0 -> 1 over 0 -> 2
    let mut links = HashMap::new();
    let mut links_0 = HashMap::new();
    links_0.insert(
        1,
        PairedLink {
            target_unitig: 1,
            support_count: 10,
            total_distance: 100.0,
        },
    );
    links_0.insert(
        2,
        PairedLink {
            target_unitig: 2,
            support_count: 1,
            total_distance: 100.0,
        },
    );
    links.insert(0, links_0);

    let paired_info = PairedInfoIndex {
        mean_insert_size: 200.0,
        insert_size_stdev: 20.0,
        links,
    };

    let resolved = expander.resolve_repeats(k, vec![u0, u1, u2], &paired_info);

    // Should resolve and stitch 0 -> 1 into "TTTTACGTAAAACGCG", leaving u2 as a separate unitig
    assert_eq!(resolved.len(), 2);
    let has_stitched = resolved.iter().any(|u| u.sequence == b"TTTTACGTAAAACGCG");
    assert!(has_stitched, "Unitig 0 and 1 should be stitched by ExSPAnder");
}

#[test]
fn test_spaligner_repeat_unrolling_two_copies() {
    let k = 21;
    let resolver = LongReadResolver {
        k,
        min_seed_matches: 1,
    };

    // Generate unique pseudo-random sequences for flanks and repeat
    let make_seq = |seed: u64, len: usize| -> Vec<u8> {
        let mut s = seed;
        let mut res = Vec::with_capacity(len);
        for _ in 0..len {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
            let b = match (s >> 32) % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
            res.push(b);
        }
        res
    };

    let flank_len = 1600;
    let repeat_len = 2000;

    let mut flank_a_seq = make_seq(111, flank_len);
    let mut flank_b_seq = make_seq(222, flank_len);
    let mut flank_c_seq = make_seq(333, flank_len);
    let mut flank_d_seq = make_seq(444, flank_len);
    let mut repeat_seq = make_seq(999, repeat_len);

    let junction_start = b"GATTACAGATTACAGATTAC"; // 20 bp (k-1)
    let junction_end = b"TGCAATGCAATGCAATGCAA";   // 20 bp (k-1)

    // Ensure exact (k-1)-mer overlaps at junctions
    flank_a_seq[flank_len - 20..].copy_from_slice(junction_start);
    flank_c_seq[flank_len - 20..].copy_from_slice(junction_start);
    repeat_seq[..20].copy_from_slice(junction_start);

    repeat_seq[repeat_len - 20..].copy_from_slice(junction_end);
    flank_b_seq[..20].copy_from_slice(junction_end);
    flank_d_seq[..20].copy_from_slice(junction_end);

    let u_flank_a = Unitig {
        id: 0,
        sequence: flank_a_seq.clone(),
        mean_coverage: 20.0,
        kmers_count: flank_len - 20,
    };
    let u_repeat = Unitig {
        id: 1,
        sequence: repeat_seq.clone(),
        mean_coverage: 40.0, // >= 1.8x median
        kmers_count: repeat_len - 20,
    };
    let u_flank_b = Unitig {
        id: 2,
        sequence: flank_b_seq.clone(),
        mean_coverage: 20.0,
        kmers_count: flank_len - 20,
    };
    let u_flank_c = Unitig {
        id: 3,
        sequence: flank_c_seq.clone(),
        mean_coverage: 20.0,
        kmers_count: flank_len - 20,
    };
    let u_flank_d = Unitig {
        id: 4,
        sequence: flank_d_seq.clone(),
        mean_coverage: 20.0,
        kmers_count: flank_len - 20,
    };

    // Long read 1: spans Flank A (last 800 bp) + Repeat R (2000 bp) + Flank B (first 800 bp)
    let mut lr1 = flank_a_seq[flank_len - 800..].to_vec();
    lr1.extend_from_slice(&repeat_seq[20..]);
    lr1.extend_from_slice(&flank_b_seq[20..800]);

    // Long read 2: spans Flank C (last 800 bp) + Repeat R (2000 bp) + Flank D (first 800 bp)
    let mut lr2 = flank_c_seq[flank_len - 800..].to_vec();
    lr2.extend_from_slice(&repeat_seq[20..]);
    lr2.extend_from_slice(&flank_d_seq[20..800]);

    let bridged = resolver.bridge_with_long_reads(
        vec![u_flank_a, u_repeat, u_flank_b, u_flank_c, u_flank_d],
        &[lr1, lr2],
    );

    // Verify that repeat node is successfully resolved into both contexts
    let probe_a = &flank_a_seq[flank_len - 100..flank_len - 20];
    let probe_b = &flank_b_seq[20..100];
    let probe_c = &flank_c_seq[flank_len - 100..flank_len - 20];
    let probe_d = &flank_d_seq[20..100];

    let has_arb = bridged.iter().any(|u| {
        u.sequence.windows(probe_a.len()).any(|w| w == probe_a)
            && u.sequence.windows(probe_b.len()).any(|w| w == probe_b)
    });
    let has_crd = bridged.iter().any(|u| {
        u.sequence.windows(probe_c.len()).any(|w| w == probe_c)
            && u.sequence.windows(probe_d.len()).any(|w| w == probe_d)
    });

    assert!(has_arb, "Bridge A-R-B must be resolved");
    assert!(has_crd, "Bridge C-R-D must be resolved");
}
