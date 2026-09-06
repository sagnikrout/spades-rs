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

pub static BIT2_TO_ASCII: [u8; 4] = *b"ACGT";

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

/// Computes reverse complement of an ASCII DNA byte slice.
pub fn revcomp_bytes(seq: &[u8]) -> Vec<u8> {
    let mut rc = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        let comp = match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            b'N' | b'n' => b'N',
            other => other,
        };
        rc.push(comp);
    }
    rc
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
    kmer = kmer.rotate_left(32);

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

/// Reverse-complement of a packed 2-bit k-mer in u128 for k <= 64.
#[inline(always)]
pub fn revcomp_kmer_u128(mut kmer: u128, k: usize) -> u128 {
    debug_assert!(k > 0 && k <= 64);
    kmer = !kmer;
    kmer = ((kmer >> 2) & 0x33333333333333333333333333333333)
        | ((kmer & 0x33333333333333333333333333333333) << 2);
    kmer = ((kmer >> 4) & 0x0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F)
        | ((kmer & 0x0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F) << 4);
    kmer = ((kmer >> 8) & 0x00FF00FF00FF00FF00FF00FF00FF00FF)
        | ((kmer & 0x00FF00FF00FF00FF00FF00FF00FF00FF) << 8);
    kmer = ((kmer >> 16) & 0x0000FFFF0000FFFF0000FFFF0000FFFF)
        | ((kmer & 0x0000FFFF0000FFFF0000FFFF0000FFFF) << 16);
    kmer = ((kmer >> 32) & 0x00000000FFFFFFFF00000000FFFFFFFF)
        | ((kmer & 0x00000000FFFFFFFF00000000FFFFFFFF) << 32);
    kmer = (kmer >> 64) | (kmer << 64);

    if k == 64 {
        kmer
    } else {
        kmer >> (2 * (64 - k))
    }
}

/// 256-bit packed 2-bit k-mer supporting any k up to 127.
/// (lo, hi): for k <= 64, hi is 0 and lo holds the k-mer in lower 2*k bits.
/// For k > 64, lo holds bases 0..64 (128 bits), and hi holds bases 64..k in lower 2*(k-64) bits.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default, Debug)]
pub struct Kmer256(pub u128, pub u128);

impl Kmer256 {
    pub const ZERO: Self = Kmer256(0, 0);

    #[inline(always)]
    pub fn is_zero(&self) -> bool {
        self.0 == 0 && self.1 == 0
    }

    #[inline(always)]
    pub fn revcomp(&self, k: usize) -> Self {
        debug_assert!(k > 0 && k <= 128);
        if k <= 64 {
            Kmer256(revcomp_kmer_u128(self.0, k), 0)
        } else {
            let k2 = k - 64;
            let lo_rc = revcomp_kmer_u128(self.0, 64);
            let hi_rc = revcomp_kmer_u128(self.1, k2);

            let shift_lo = 2 * (64 - k2);
            let rc_lo = if shift_lo >= 128 { 0 } else { hi_rc << shift_lo }
                | if k2 >= 64 { 0 } else { lo_rc >> (2 * k2) };
            let mask_hi = if k2 >= 64 { u128::MAX } else { (1u128 << (2 * k2)) - 1 };
            let rc_hi = lo_rc & mask_hi;
            Kmer256(rc_lo, rc_hi)
        }
    }

    #[inline(always)]
    pub fn canonical(&self, k: usize) -> (Self, bool) {
        let rc = self.revcomp(k);
        if *self <= rc {
            (*self, false)
        } else {
            (rc, true)
        }
    }

    #[inline(always)]
    pub fn extend_right(&self, b: u8, k: usize) -> Self {
        debug_assert!(k > 0 && k <= 128);
        if k <= 64 {
            let mask = if k == 64 { u128::MAX } else { (1u128 << (2 * k)) - 1 };
            Kmer256(((self.0 << 2) | (b as u128)) & mask, 0)
        } else {
            let k2 = k - 64;
            let mask_hi = if k2 == 64 { u128::MAX } else { (1u128 << (2 * k2)) - 1 };
            let s64 = (self.1 >> (2 * (k2 - 1))) & 3;
            let new_lo = (self.0 << 2) | s64;
            let new_hi = ((self.1 << 2) | (b as u128)) & mask_hi;
            Kmer256(new_lo, new_hi)
        }
    }

    #[inline(always)]
    pub fn prepend_left(&self, a: u8, k: usize) -> Self {
        debug_assert!(k > 0 && k <= 128);
        if k <= 64 {
            Kmer256((self.0 >> 2) | ((a as u128) << (2 * (k - 1))), 0)
        } else {
            let k2 = k - 64;
            let mask_hi = if k2 == 64 { u128::MAX } else { (1u128 << (2 * k2)) - 1 };
            let s63 = self.0 & 3;
            let new_lo = (self.0 >> 2) | ((a as u128) << 126);
            let new_hi = ((self.1 >> 2) | (s63 << (2 * (k2 - 1)))) & mask_hi;
            Kmer256(new_lo, new_hi)
        }
    }

    #[inline(always)]
    pub fn last_base(&self, k: usize) -> u8 {
        if k <= 64 {
            (self.0 & 3) as u8
        } else {
            (self.1 & 3) as u8
        }
    }

    #[inline(always)]
    pub fn first_base(&self, k: usize) -> u8 {
        if k <= 64 {
            ((self.0 >> (2 * (k - 1))) & 3) as u8
        } else {
            ((self.0 >> 126) & 3) as u8
        }
    }

    pub fn to_string(&self, k: usize) -> String {
        let mut bytes = vec![0u8; k];
        if k <= 64 {
            let mut val = self.0;
            for i in (0..k).rev() {
                bytes[i] = bit2_to_base((val & 3) as u8);
                val >>= 2;
            }
        } else {
            let mut val_lo = self.0;
            for i in (0..64).rev() {
                bytes[i] = bit2_to_base((val_lo & 3) as u8);
                val_lo >>= 2;
            }
            let mut val_hi = self.1;
            for i in (64..k).rev() {
                bytes[i] = bit2_to_base((val_hi & 3) as u8);
                val_hi >>= 2;
            }
        }
        String::from_utf8(bytes).unwrap()
    }

    pub fn from_bytes(s: &[u8], k: usize) -> Option<Self> {
        if s.len() < k || k > 128 {
            return None;
        }
        let mut lo = 0u128;
        let mut hi = 0u128;
        if k <= 64 {
            for &b in &s[..k] {
                let code = base_to_2bit(b)?;
                lo = (lo << 2) | (code as u128);
            }
        } else {
            for &b in &s[..64] {
                let code = base_to_2bit(b)?;
                lo = (lo << 2) | (code as u128);
            }
            for &b in &s[64..k] {
                let code = base_to_2bit(b)?;
                hi = (hi << 2) | (code as u128);
            }
        }
        Some(Kmer256(lo, hi))
    }
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

    #[test]
    fn test_kmer256_roundtrip_and_revcomp() {
        let template = "ATGCGATCGATCGATAGCTAGCTAGCTAGCTAAGCTAGCTAGCTAGCTAATGCGATCGATCGATAGCTAGCTAGCTAGCTAAGCTAGCTAGCTAGCTAATGCGATCGATCGATAGCTAGCTAGCTAGCTAAGCTAGCTAGCTAGCTAATGCGATCGATCGATAGCTAGCTAGCTAGCTAAGCTAGCTAGCTAGCTA";
        for &k in &[21, 33, 55, 64, 77, 99, 127] {
            let seq = template[..k].as_bytes();
            let km = Kmer256::from_bytes(seq, k).expect("encoding failed");
            let s_dec = km.to_string(k);
            assert_eq!(std::str::from_utf8(seq).unwrap(), s_dec);

            // Double reverse complement is identity
            let rc = km.revcomp(k);
            let rc_rc = rc.revcomp(k);
            assert_eq!(km, rc_rc);
        }
    }

    #[test]
    fn test_kmer256_extend_and_prepend() {
        let template = "ATGCGATCGATCGATAGCTAGCTAGCTAGCTAAGCTAGCTAGCTAGCTAATGCGATCGATCGATAGCTAGCTAGCTAGCTAAGCTAGCTAGCTAGCTAATGCGATCGATCGATAGCTAGCTAGCTAGCTAAGCTAGCTAGCTAGCTAATGCGATCGATCGATAGCTAGCTAGCTAGCTAAGCTAGCTAGCTAGCTA";
        for &k in &[21, 33, 55, 64, 77, 99, 127] {
            let seq = &template[..k];
            let km = Kmer256::from_bytes(seq.as_bytes(), k).unwrap();

            // Extend right with 'C' (1)
            let ext = km.extend_right(1, k);
            let expected_ext = format!("{}C", &seq[1..]);
            assert_eq!(ext.to_string(k), expected_ext);

            // Prepend left with 'G' (2)
            let prep = km.prepend_left(2, k);
            let expected_prep = format!("G{}", &seq[..k - 1]);
            assert_eq!(prep.to_string(k), expected_prep);
        }
    }
}
