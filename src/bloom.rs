//! Lock-free concurrent Two-Tier Bloom Filter (The Memory Shield).
//!
//! Eliminates singleton sequencing error k-mers before graph allocation.
//! Uses atomic bit-vectors for concurrent access across all 22 threads without lock contention.

use std::sync::atomic::{AtomicU64, Ordering};

pub struct TwoTierFilter {
    size_bits: usize,
    num_words: usize,
    seen_bits: Vec<AtomicU64>,
    solid_bits: Vec<AtomicU64>,
}

impl TwoTierFilter {
    /// Creates a filter with `capacity` bits (rounded up to 64-bit words).
    pub fn new(capacity_bits: usize) -> Self {
        let num_words = (capacity_bits + 63) / 64;
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

    /// Fast 64-bit integer mix hash (Murmur / SplitMix).
    #[inline(always)]
    fn hash_kmer(mut x: u64, seed: u64) -> u64 {
        x ^= seed;
        x = x.wrapping_mul(0xff51afd7ed558ccd);
        x ^= x >> 33;
        x = x.wrapping_mul(0xc4ceb9fe1a85ec53);
        x ^= x >> 33;
        x
    }

    /// Observes a k-mer. If already observed, marks it as 'solid' (coverage >= 2).
    /// Returns true if it was already marked solid.
    #[inline(always)]
    pub fn insert(&self, kmer: u64) -> bool {
        let h1 = Self::hash_kmer(kmer, 0x517cc1b727220a95) as usize % self.size_bits;
        let h2 = Self::hash_kmer(kmer, 0x9e3779b97f4a7c15) as usize % self.size_bits;

        let w1 = h1 / 64;
        let b1 = 1u64 << (h1 % 64);

        let w2 = h2 / 64;
        let b2 = 1u64 << (h2 % 64);

        // Check if both bits are set in seen
        let prev1 = self.seen_bits[w1].fetch_or(b1, Ordering::Relaxed);
        let prev2 = self.seen_bits[w2].fetch_or(b2, Ordering::Relaxed);

        if (prev1 & b1 != 0) && (prev2 & b2 != 0) {
            // Already seen at least once! Now mark solid
            self.solid_bits[w1].fetch_or(b1, Ordering::Relaxed);
            self.solid_bits[w2].fetch_or(b2, Ordering::Relaxed);
            true
        } else {
            false
        }
    }

    /// Checks if a k-mer is considered solid (seen at least twice).
    #[inline(always)]
    pub fn is_solid(&self, kmer: u64) -> bool {
        let h1 = Self::hash_kmer(kmer, 0x517cc1b727220a95) as usize % self.size_bits;
        let h2 = Self::hash_kmer(kmer, 0x9e3779b97f4a7c15) as usize % self.size_bits;

        let w1 = h1 / 64;
        let b1 = 1u64 << (h1 % 64);

        let w2 = h2 / 64;
        let b2 = 1u64 << (h2 % 64);

        (self.solid_bits[w1].load(Ordering::Relaxed) & b1 != 0)
            && (self.solid_bits[w2].load(Ordering::Relaxed) & b2 != 0)
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
