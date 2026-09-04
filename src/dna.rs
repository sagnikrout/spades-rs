//! Fast 2-bit DNA encoding and SIMD/bit-twiddling primitives.
//!
//! A = 00 (0)
//! C = 01 (1)
//! G = 10 (2)
//! T = 11 (3)
//!
//! Under this representation, complement is simply: base ^ 3 (or !base & 3).

/// Lookup table mapping ASCII byte to 2-bit value (0..3), or 0xFF if invalid.
pub static ASCII_TO_2BIT: [u8; 256] = {
    let mut table = [0xFF; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

pub static BIT2_TO_ASCII: [u8; 4] = [b'A', b'C', b'G', b'T'];

#[inline(always)]
pub fn base_to_2bit(b: u8) -> Option<u8> {
    let val = ASCII_TO_2BIT[b as usize];
    if val != 0xFF {
        Some(val)
    } else {
        None
    }
}

#[inline(always)]
pub fn bit2_to_base(val: u8) -> u8 {
    BIT2_TO_ASCII[(val & 3) as usize]
}

/// Reverse-complement of a packed 2-bit k-mer stored in the lower 2*k bits of a u64.
/// Executes in constant time without branches using bit-permutation swaps.
#[inline(always)]
pub fn revcomp_kmer_u64(mut kmer: u64, k: usize) -> u64 {
    debug_assert!(k > 0 && k <= 32);
    // Complement all 2-bit pairs
    kmer = !kmer;
    // Swap adjacent 2-bit pairs
    kmer = ((kmer >> 2) & 0x3333333333333333) | ((kmer & 0x3333333333333333) << 2);
    // Swap 4-bit nibbles
    kmer = ((kmer >> 4) & 0x0F0F0F0F0F0F0F0F) | ((kmer & 0x0F0F0F0F0F0F0F0F) << 4);
    // Swap bytes
    kmer = ((kmer >> 8) & 0x00FF00FF00FF00FF) | ((kmer & 0x00FF00FF00FF00FF) << 8);
    // Swap 16-bit words
    kmer = ((kmer >> 16) & 0x0000FFFF0000FFFF) | ((kmer & 0x0000FFFF0000FFFF) << 16);
    // Swap 32-bit halves
    kmer = (kmer >> 32) | (kmer << 32);

    // Shift down so the k-mer occupies the lowest 2*k bits
    kmer >> (2 * (32 - k))
}

/// Returns the canonical k-mer (lexicographical minimum of forward and reverse-complement).
#[inline(always)]
pub fn canonical_kmer_u64(kmer: u64, k: usize) -> (u64, bool) {
    let rc = revcomp_kmer_u64(kmer, k);
    if kmer <= rc {
        (kmer, false) // forward
    } else {
        (rc, true) // reverse complement
    }
}

/// Decodes a packed 2-bit k-mer into an ASCII String.
pub fn kmer_to_string(mut kmer: u64, k: usize) -> String {
    let mut bytes = vec![0u8; k];
    for i in (0..k).rev() {
        bytes[i] = bit2_to_base((kmer & 3) as u8);
        kmer >>= 2;
    }
    String::from_utf8(bytes).unwrap()
}

/// Encodes an ASCII DNA string slice of length k into a packed u64.
/// Returns None if any non-ACGT character is encountered.
pub fn string_to_kmer(s: &[u8], k: usize) -> Option<u64> {
    if s.len() < k {
        return None;
    }
    let mut val = 0u64;
    for &b in &s[..k] {
        let code = base_to_2bit(b)?;
        val = (val << 2) | (code as u64);
    }
    Some(val)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_revcomp_kmer() {
        // Forward: "ACGT" -> A=00, C=01, G=10, T=11 => 00 01 10 11 (0x1B)
        // Reverse complement of ACGT is ACGT (it's a palindrome)
        let kmer = string_to_kmer(b"ACGT", 4).unwrap();
        let rc = revcomp_kmer_u64(kmer, 4);
        assert_eq!(kmer, rc);

        // Forward: "AAAA" -> 00 00 00 00
        // RC: "TTTT" -> 11 11 11 11
        let a4 = string_to_kmer(b"AAAA", 4).unwrap();
        let t4 = string_to_kmer(b"TTTT", 4).unwrap();
        assert_eq!(revcomp_kmer_u64(a4, 4), t4);
        assert_eq!(revcomp_kmer_u64(t4, 4), a4);
    }

    #[test]
    fn test_encode_decode_roundtrip() {
        let original = "ATGCGATCGATCGAT";
        let k = original.len();
        let encoded = string_to_kmer(original.as_bytes(), k).unwrap();
        let decoded = kmer_to_string(encoded, k);
        assert_eq!(original, decoded);
    }
}
