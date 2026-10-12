//! Lock-free concurrent two-tier Bloom filter for $k$-mer abundance screening.
//!
//! Filters out singleton sequencing errors (k-mers seen only once) before allocating
//! the hash-based de Bruijn graph:
//! - Tier 1 (`seen_bits`): records whether a $k$-mer has been observed at least once.
//! - Tier 2 (`solid_bits`): promotes a $k$-mer when all 3 independent hash positions were
//!   already set in `seen_bits` (indicating coverage $\ge 2$).

use crate::dna::Kmer256;
use std::sync::atomic::{AtomicU64, Ordering};

/// Lock-free two-tier bit-vector Bloom filter backed by 64-bit atomic words.
pub struct TwoTierFilter {
    /// Total number of addressable bits in each tier (always a multiple of 64).
    size_bits: usize,
    /// Number of 64-bit atomic words in each tier.
    num_words: usize,
    /// Tier 1 bit-vector: bits set on first observation of a k-mer.
    seen_bits: Vec<AtomicU64>,
    /// Tier 2 bit-vector: bits set on second and subsequent observations of a k-mer.
    solid_bits: Vec<AtomicU64>,
}

impl TwoTierFilter {
    /// Creates a two-tier filter with `capacity_bits` per tier (rounded up to 64-bit words).
    pub fn new(capacity_bits: usize) -> Self {
        let num_words = capacity_bits.div_ceil(64);
        let size_bits = num_words * 64;

        let mut seen = Vec::with_capacity(num_words);
        let mut solid = Vec::with_capacity(num_words);
        for _ in 0..num_words {
            seen.push(AtomicU64::new(0));
            solid.push(AtomicU64::new(0));
        }

        Self {
            size_bits,
            num_words,
            seen_bits: seen,
            solid_bits: solid,
        }
    }

    /// Observes a 64-bit k-mer. If already observed, marks it as 'solid' (coverage >= 2).
    /// Returns true if it was already observed at least once before this call.
    #[inline(always)]
    pub fn insert(&self, kmer: u64) -> bool {
        self.insert_kmer256(Kmer256(kmer as u128, 0))
    }

    /// Checks if a 64-bit k-mer is considered solid (seen at least twice).
    #[inline(always)]
    pub fn is_solid(&self, kmer: u64) -> bool {
        self.is_solid_kmer256(Kmer256(kmer as u128, 0))
    }

    /// Folds a 256-bit k-mer into a 64-bit hash using MurmurHash3 finalizer mixing.
    /// For `k <= 64`, `kmer.1` is 0 and only `kmer.0` contributes; for `k > 64`, both halves are mixed.
    #[inline(always)]
    fn hash_kmer256(kmer: Kmer256, seed: u64) -> u64 {
        let x0 = (kmer.0 as u64) ^ ((kmer.0 >> 64) as u64);
        let x1 = (kmer.1 as u64) ^ ((kmer.1 >> 64) as u64);
        let mut x = x0 ^ x1 ^ seed;
        x = x.wrapping_mul(0xff51afd7ed558ccd);
        x ^= x >> 33;
        x = x.wrapping_mul(0xc4ceb9fe1a85ec53);
        x ^= x >> 33;
        x
    }

    /// Probabilistically observes a `Kmer256`: sets bits in `seen_bits` on first sight and promotes
    /// to `solid_bits` when all 3 hash positions were already set in `seen_bits` (coverage >= 2).
    #[inline(always)]
    pub fn insert_kmer256(&self, kmer: Kmer256) -> bool {
        let h1 = Self::hash_kmer256(kmer, 0x517cc1b727220a95) as usize % self.size_bits;
        let h2 = Self::hash_kmer256(kmer, 0x9e3779b97f4a7c15) as usize % self.size_bits;
        let h3 = Self::hash_kmer256(kmer, 0x6c62272e07bb0142) as usize % self.size_bits;

        let w1 = h1 / 64;
        let b1 = 1u64 << (h1 % 64);
        let w2 = h2 / 64;
        let b2 = 1u64 << (h2 % 64);
        let w3 = h3 / 64;
        let b3 = 1u64 << (h3 % 64);

        let prev1 = self.seen_bits[w1].fetch_or(b1, Ordering::Relaxed);
        let prev2 = self.seen_bits[w2].fetch_or(b2, Ordering::Relaxed);
        let prev3 = self.seen_bits[w3].fetch_or(b3, Ordering::Relaxed);

        if (prev1 & b1 != 0) && (prev2 & b2 != 0) && (prev3 & b3 != 0) {
            self.solid_bits[w1].fetch_or(b1, Ordering::Relaxed);
            self.solid_bits[w2].fetch_or(b2, Ordering::Relaxed);
            self.solid_bits[w3].fetch_or(b3, Ordering::Relaxed);
            true
        } else {
            false
        }
    }

    /// Returns `true` if all three hash positions for `kmer` are set in `solid_bits` (coverage >= 2).
    #[inline(always)]
    pub fn is_solid_kmer256(&self, kmer: Kmer256) -> bool {
        let h1 = Self::hash_kmer256(kmer, 0x517cc1b727220a95) as usize % self.size_bits;
        let h2 = Self::hash_kmer256(kmer, 0x9e3779b97f4a7c15) as usize % self.size_bits;
        let h3 = Self::hash_kmer256(kmer, 0x6c62272e07bb0142) as usize % self.size_bits;

        let w1 = h1 / 64;
        let b1 = 1u64 << (h1 % 64);
        let w2 = h2 / 64;
        let b2 = 1u64 << (h2 % 64);
        let w3 = h3 / 64;
        let b3 = 1u64 << (h3 % 64);

        (self.solid_bits[w1].load(Ordering::Relaxed) & b1 != 0)
            && (self.solid_bits[w2].load(Ordering::Relaxed) & b2 != 0)
            && (self.solid_bits[w3].load(Ordering::Relaxed) & b3 != 0)
    }

    /// Total memory occupied by the filter in bytes.
    pub fn memory_usage_bytes(&self) -> usize {
        self.num_words * 8 * 2
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_two_tier_filter() {
        let filter = TwoTierFilter::new(1024 * 1024);
        let kmer1 = 123456789u64;
        let kmer2 = 987654321u64;

        // First insertion: not yet solid
        assert!(!filter.is_solid(kmer1));
        filter.insert(kmer1);
        // Still not solid after 1 observation
        assert!(!filter.is_solid(kmer1));

        // Second insertion: becomes solid
        filter.insert(kmer1);
        assert!(filter.is_solid(kmer1));

        // kmer2 was never inserted
        assert!(!filter.is_solid(kmer2));
    }
}
