use intelligent_pascal::packed_reads::PackedReads;

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
        assert_eq!(&buf, expected, "Mismatch at read index {} of length {}", i, expected.len());
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
        assert_eq!(&scratch, expected, "Data corruption detected in read #{}", idx);
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
