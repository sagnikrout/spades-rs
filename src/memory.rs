//! Memory Management & Governor Engine.
//!
//! Automatically detects available system RAM and sets the maximum memory
//! budget to (Available RAM - 20%), leaving a 20% safety margin for the OS
//! and disk cache. Supports user override via `--max-memory <GB>`.

use std::fs;

/// System memory hardware profile and configured assembly RAM budget.
#[derive(Debug, Clone, Copy)]
pub struct MemoryLimits {
    /// Total physical RAM installed on the host system, in bytes.
    pub total_system_bytes: usize,
    /// Currently available (free + reclaimable cache) physical RAM on the host system, in bytes.
    pub available_system_bytes: usize,
    /// Maximum memory budget allocated to the assembly run, in bytes.
    pub max_budget_bytes: usize,
    /// Whether `max_budget_bytes` was explicitly overridden by the user (`--max-memory`).
    pub is_user_specified: bool,
}

impl MemoryLimits {
    /// Detects system memory and computes default budget (Available RAM - 20%)
    /// or applies user-specified override.
    pub fn determine(user_override_gb: Option<f64>) -> Self {
        let (total, available) = detect_system_memory();

        let (budget, is_user) = match user_override_gb {
            Some(gb) if gb.is_finite() && gb > 0.0 => {
                let bytes = (gb * 1024.0 * 1024.0 * 1024.0) as usize;
                (bytes, true)
            }
            _ => {
                // Default: Available RAM - 20% (i.e. 80% of available memory)
                let budget = (available as f64 * 0.80) as usize;
                // Floor at 512 MB (capped by available RAM) to ensure minimal working buffer
                (budget.max(512 * 1024 * 1024).min(available), false)
            }
        };

        Self {
            total_system_bytes: total,
            available_system_bytes: available,
            max_budget_bytes: budget,
            is_user_specified: is_user,
        }
    }

    /// Returns the configured maximum memory budget in gigabytes (GiB).
    pub fn budget_gb(&self) -> f64 {
        self.max_budget_bytes as f64 / (1024.0 * 1024.0 * 1024.0)
    }

    /// Returns the detected available system memory in gigabytes (GiB).
    pub fn available_gb(&self) -> f64 {
        self.available_system_bytes as f64 / (1024.0 * 1024.0 * 1024.0)
    }

    /// Returns the total installed system memory in gigabytes (GiB).
    pub fn total_gb(&self) -> f64 {
        self.total_system_bytes as f64 / (1024.0 * 1024.0 * 1024.0)
    }

    /// Computes optimal Bloom filter bits based on memory budget.
    /// Allocates an optimal ~3-6% of memory to the transient Bloom filter (which is
    /// dropped immediately after k-mer screening), preventing false-positive noise
    /// while leaving >90% of budget for graph compaction.
    pub fn optimal_bloom_bits(&self) -> usize {
        let budget_mb = self.max_budget_bytes / (1024 * 1024);
        if budget_mb < 2048 {
            // Under 2 GB budget: 32 MB per tier (256M bits = 64 MB total)
            256 * 1024 * 1024
        } else if budget_mb < 8192 {
            // 2 - 8 GB budget: 64 MB per tier (512M bits = 128 MB total)
            512 * 1024 * 1024
        } else if budget_mb < 16384 {
            // 8 - 16 GB budget: 128 MB per tier (1024M bits = 256 MB total)
            1024 * 1024 * 1024
        } else if budget_mb < 32768 {
            // 16 - 32 GB budget: 256 MB per tier (2048M bits = 512 MB total)
            2048 * 1024 * 1024
        } else if budget_mb < 65536 {
            // 32 - 64 GB workstation: 512 MB per tier (4096M bits = 1 GB total)
            4096 * 1024 * 1024
        } else if budget_mb < 131072 {
            // 64 - 128 GB server: 1 GB per tier (8192M bits = 2 GB total)
            8192 * 1024 * 1024
        } else {
            // > 128 GB high-memory server: 2 GB per tier (16384M bits = 4 GB total)
            16384 * 1024 * 1024
        }
    }

    /// Gets current physical resident set size (RSS) in bytes.
    pub fn current_rss_bytes() -> usize {
        #[cfg(target_os = "linux")]
        {
            if let Ok(content) = fs::read_to_string("/proc/self/status") {
                for line in content.lines() {
                    if line.starts_with("VmRSS:") {
                        let kb = line
                            .split_whitespace()
                            .nth(1)
                            .and_then(|v| v.parse::<usize>().ok())
                            .unwrap_or(0);
                        return kb * 1024;
                    }
                }
            }
        }
        #[cfg(target_os = "windows")]
        {
            #[repr(C)]
            struct ProcessMemoryCounters {
                cb: u32,
                page_fault_count: u32,
                peak_working_set_size: usize,
                working_set_size: usize,
                quota_peak_paged_pool_usage: usize,
                quota_paged_pool_usage: usize,
                quota_peak_non_paged_pool_usage: usize,
                quota_non_paged_pool_usage: usize,
                pagefile_usage: usize,
                peak_pagefile_usage: usize,
            }
            extern "system" {
                fn GetCurrentProcess() -> isize;
                fn K32GetProcessMemoryInfo(
                    process: isize,
                    counters: *mut ProcessMemoryCounters,
                    cb: u32,
                ) -> i32;
            }
            unsafe {
                let mut counters = std::mem::zeroed::<ProcessMemoryCounters>();
                counters.cb = std::mem::size_of::<ProcessMemoryCounters>() as u32;
                if K32GetProcessMemoryInfo(GetCurrentProcess(), &mut counters, counters.cb) != 0 {
                    return counters.working_set_size;
                }
            }
        }
        0
    }

    /// Checks if current memory exceeds warning or critical thresholds.
    pub fn check_headroom(&self) -> Result<(), String> {
        let rss = Self::current_rss_bytes();
        if rss > 0 && rss >= self.max_budget_bytes {
            return Err(format!(
                "Memory budget exceeded: current RSS {:.2} GB >= budget {:.2} GB",
                rss as f64 / (1024.0 * 1024.0 * 1024.0),
                self.budget_gb()
            ));
        }
        Ok(())
    }
}

/// Trims the process heap back to the OS if supported (glibc Linux).
#[inline]
pub fn trim_memory() {
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    unsafe {
        extern "C" {
            fn malloc_trim(pad: usize) -> i32;
        }
        malloc_trim(0);
    }
}

#[cfg(target_os = "windows")]
fn detect_windows_memory() -> Option<(usize, usize)> {
    #[repr(C)]
    struct MemoryStatusEx {
        length: u32,
        memory_load: u32,
        total_phys: u64,
        avail_phys: u64,
        total_page_file: u64,
        avail_page_file: u64,
        total_virtual: u64,
        avail_virtual: u64,
        avail_extended_virtual: u64,
    }
    extern "system" {
        fn GlobalMemoryStatusEx(stat: *mut MemoryStatusEx) -> i32;
    }
    unsafe {
        let mut stat = std::mem::zeroed::<MemoryStatusEx>();
        stat.length = std::mem::size_of::<MemoryStatusEx>() as u32;
        if GlobalMemoryStatusEx(&mut stat) != 0 {
            return Some((stat.total_phys as usize, stat.avail_phys as usize));
        }
    }
    None
}

fn detect_system_memory() -> (usize, usize) {
    // 1. Try Windows API on Windows
    #[cfg(target_os = "windows")]
    {
        if let Some(mem) = detect_windows_memory() {
            return mem;
        }
    }

    // 2. Try /proc/meminfo (Linux / WSL)
    if let Ok(content) = fs::read_to_string("/proc/meminfo") {
        let mut total_kb = 0;
        let mut avail_kb = 0;
        for line in content.lines() {
            if line.starts_with("MemTotal:") {
                total_kb = parse_meminfo_kb(line);
            } else if line.starts_with("MemAvailable:") {
                avail_kb = parse_meminfo_kb(line);
            }
        }
        if total_kb > 0 {
            if avail_kb == 0 {
                // If MemAvailable not present on older kernels, estimate ~70% of total
                avail_kb = (total_kb as f64 * 0.70) as usize;
            }
            return (total_kb * 1024, avail_kb * 1024);
        }
    }

    // 3. Safe fallback if /proc/meminfo not available (e.g. macOS)
    let fallback_total = 16 * 1024 * 1024 * 1024; // 16 GB
    let fallback_avail = 12 * 1024 * 1024 * 1024; // 12 GB
    (fallback_total, fallback_avail)
}

fn parse_meminfo_kb(line: &str) -> usize {
    line.split_whitespace()
        .nth(1)
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_memory_limits_calculation() {
        let limits = MemoryLimits::determine(None);
        assert!(limits.total_system_bytes > 0);
        assert!(limits.available_system_bytes > 0);
        assert!(limits.max_budget_bytes > 0);
        assert!(limits.max_budget_bytes <= limits.available_system_bytes);
        assert!(!limits.is_user_specified);
    }

    #[test]
    fn test_memory_limits_override() {
        let limits = MemoryLimits::determine(Some(4.5));
        assert_eq!(limits.budget_gb(), 4.5);
        assert!(limits.is_user_specified);
    }
}
