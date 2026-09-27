use spades_rs::packed_reads::PackedReads;

#[test]
fn test_packed_reads_all_modulo4_lengths() {
    let mut store = PackedReads::default();
    let mut expected_reads = Vec::new();

    // Test every length from 0 to 64
    for len in 0..=64 {
        let seq: Vec<u8> = (0..len)
            .map(|i| match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();
        store.add_read(&seq);
        expected_reads.push(seq);
    }

    assert_eq!(store.len(), 65);
    assert!(!store.is_empty());

    let mut buf = Vec::new();
    for (i, expected) in expected_reads.iter().enumerate() {
        store.get_read(i, &mut buf);
        assert_eq!(
            &buf,
            expected,
            "Mismatch at read index {} of length {}",
            i,
            expected.len()
        );
    }
}

#[test]
fn test_packed_reads_synthetic_workload() {
    let mut store = PackedReads::with_capacity(2000, 2000 * 150);
    let mut expected_reads = Vec::with_capacity(2000);

    // Generate 2,000 reads with varied lengths (50 to 250 bp)
    let mut rng_state = 0x12345678u64;
    for _i in 0..2000 {
        rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
        let len = 50 + (rng_state % 200) as usize;

        let mut seq = Vec::with_capacity(len);
        for j in 0..len {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let base = match (rng_state ^ (j as u64)) % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
            seq.push(base);
        }

        store.add_read(&seq);
        expected_reads.push(seq);
    }

    assert_eq!(store.len(), 2000);
    assert!(store.memory_usage_bytes() > 0);

    let mut scratch = Vec::new();
    for (idx, expected) in expected_reads.iter().enumerate() {
        store.get_read(idx, &mut scratch);
        assert_eq!(
            &scratch, expected,
            "Data corruption detected in read #{}",
            idx
        );
    }
}

#[test]
fn test_packed_reads_boundary_cases() {
    let mut store = PackedReads::default();

    // 0-length read
    store.add_read(&[]);
    // 1-bp read
    store.add_read(b"G");
    // 2-bp read
    store.add_read(b"TA");
    // 3-bp read
    store.add_read(b"CGC");
    // 4-bp read (exactly 1 byte)
    store.add_read(b"ACGT");
    // 5-bp read (1 byte + 1 base)
    store.add_read(b"ACGTA");

    let mut buf = Vec::new();
    store.get_read(0, &mut buf);
    assert_eq!(buf, b"");

    store.get_read(1, &mut buf);
    assert_eq!(buf, b"G");

    store.get_read(2, &mut buf);
    assert_eq!(buf, b"TA");

    store.get_read(3, &mut buf);
    assert_eq!(buf, b"CGC");

    store.get_read(4, &mut buf);
    assert_eq!(buf, b"ACGT");

    store.get_read(5, &mut buf);
    assert_eq!(buf, b"ACGTA");
}

// ─────────────────────────────────────────────────────────────────────────────
// N-handling integration tests (tests/test_packed_reads.rs)
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_n_end_trim_leading_enters_primary_store() {
    // "NNN" + 21 clean bases → trimmed to 21 bp, packed to primary store.
    let read = b"NNNACGTACGTACGTACGTACGTA"; // NNN + 21 bp clean
    let mut store = PackedReads::default();
    store.ingest_read(read);

    assert_eq!(store.len(), 1, "Trimmed read must be in primary store");
    assert!(store.ambiguous_sidecar.is_empty());

    let mut buf = Vec::new();
    store.get_read(0, &mut buf);
    assert_eq!(buf, b"ACGTACGTACGTACGTACGTA");
    assert!(!buf.contains(&b'N'));
}

#[test]
fn test_n_end_trim_trailing_enters_primary_store() {
    let read = b"ACGTACGTACGTACGTACGTANNN";
    let mut store = PackedReads::default();
    store.ingest_read(read);

    assert_eq!(store.len(), 1);
    assert!(store.ambiguous_sidecar.is_empty());

    let mut buf = Vec::new();
    store.get_read(0, &mut buf);
    // 21 bytes, last of which is 'A'
    assert_eq!(buf.len(), 21);
    assert!(!buf.contains(&b'N'));
}

#[test]
fn test_n_internal_split_both_halves_long_enough() {
    // 21 bp + NNN + 21 bp → two sub-reads in primary store, sidecar empty.
    let read = b"ACGTACGTACGTACGTACGTANNNACGTACGTACGTACGTACGTA";
    let mut store = PackedReads::default();
    store.ingest_read(read);

    assert_eq!(store.len(), 2, "Both 21+ bp sub-reads must be packed");
    assert!(store.ambiguous_sidecar.is_empty());
}

#[test]
fn test_n_internal_split_one_half_short_one_long() {
    // 4 bp short sub-read + N + 21 bp long sub-read.
    // Short part (<21) is discarded. Long part goes to primary store. No sidecar.
    let mut read = b"ACGTN".to_vec();
    read.extend_from_slice(b"ACGTACGTACGTACGTACGTA"); // 21 bp
    let mut store = PackedReads::default();
    store.ingest_read(&read);

    assert_eq!(store.len(), 1, "Only the 21+ bp sub-read should be packed");
    assert!(store.ambiguous_sidecar.is_empty());
}

#[test]
fn test_n_all_sub_reads_short_goes_to_sidecar() {
    // "ACGT"(4) + N + "CGTA"(4) — neither ≥ 21 → sidecar, primary store empty.
    let mut store = PackedReads::default();
    store.ingest_read(b"ACGTNNCGTA");

    assert_eq!(store.len(), 0);
    assert_eq!(store.ambiguous_sidecar.len(), 1);
}

#[test]
fn test_n_all_n_read_discarded_completely() {
    let mut store = PackedReads::default();
    store.ingest_read(b"NNNNNNNNNNN");

    assert_eq!(store.len(), 0);
    assert!(
        store.ambiguous_sidecar.is_empty(),
        "All-N reads must be fully discarded"
    );
}

#[test]
fn test_no_artificial_a_from_n_bug() {
    // Regression test: original `unwrap_or(0)` mapped N → 'A' (00 in 2-bit).
    // Both halves of "ACGTNACGT" are 4 bp (<21), so the read goes to sidecar.
    let mut store = PackedReads::default();
    store.ingest_read(b"ACGTNACGT");

    // Primary store must be empty → no A-corruption possible at the bit level.
    assert_eq!(store.len(), 0);

    // Sidecar entry must retain the raw N byte — not substitute A.
    assert!(!store.ambiguous_sidecar.is_empty());
    assert!(
        store.ambiguous_sidecar[0].contains(&b'N'),
        "Sidecar must preserve the N character, not corrupt it to A"
    );
}

#[test]
fn test_sidecar_n_position_contributes_no_primary_kmer() {
    // Reads in the sidecar are absent from the primary 2-bit index, so they
    // cannot contaminate the Bloom filter or de Bruijn graph with spurious k-mers.
    let mut store = PackedReads::default();
    store.ingest_read(b"ACGTNACGT"); // both halves < 21 bp → sidecar

    let num_primary = store.len();
    assert_eq!(
        num_primary, 0,
        "Sidecar reads must not be iterable via the primary packed store"
    );
}

#[test]
fn test_ingest_pure_acgt_read_bypasses_sidecar() {
    // A completely clean read (≥21 bp) should go straight to primary, sidecar untouched.
    let read = b"ACGTACGTACGTACGTACGTACGT"; // 24 bp, all clean
    let mut store = PackedReads::default();
    store.ingest_read(read);

    assert_eq!(store.len(), 1);
    assert!(store.ambiguous_sidecar.is_empty());

    let mut buf = Vec::new();
    store.get_read(0, &mut buf);
    assert_eq!(buf, read);
}

#[test]
fn test_sidecar_memory_included_in_usage_bytes() {
    let mut store = PackedReads::default();
    store.ingest_read(b"ACGTNACGT"); // → sidecar
    let usage = store.memory_usage_bytes();
    assert!(
        usage > 0,
        "Sidecar memory must be counted in memory_usage_bytes()"
    );
}
