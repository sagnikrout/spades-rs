//! Memory Management & Governor Engine.
//!
//! Automatically detects available system RAM and sets the maximum memory
//! budget to (Available RAM - 20%), leaving a 20% safety margin for the OS
//! and disk cache. Supports user override via `--max-memory <GB>`.

use std::fs;

#[derive(Debug, Clone, Copy)]
pub struct MemoryLimits {
    pub total_system_bytes: usize,
    pub available_system_bytes: usize,
    pub max_budget_bytes: usize,
    pub is_user_specified: bool,
}

impl MemoryLimits {
    /// Detects system memory and computes default budget (Available RAM - 20%)
    /// or applies user-specified override.
    pub fn determine(user_override_gb: Option<f64>) -> Self {
        let (total, available) = detect_system_memory();

        let (budget, is_user) = match user_override_gb {
            Some(gb) => {
                let bytes = (gb * 1024.0 * 1024.0 * 1024.0) as usize;
                (bytes, true)
            }
            None => {
                // Default: Available RAM - 20% (i.e. 80% of available memory)
                let budget = (available as f64 * 0.80) as usize;
                // Floor at 512 MB to ensure minimal working buffer
                (budget.max(512 * 1024 * 1024), false)
            }
        };

        Self {
            total_system_bytes: total,
            available_system_bytes: available,
            max_budget_bytes: budget,
            is_user_specified: is_user,
        }
    }

    pub fn budget_gb(&self) -> f64 {
        self.max_budget_bytes as f64 / (1024.0 * 1024.0 * 1024.0)
    }

    pub fn available_gb(&self) -> f64 {
        self.available_system_bytes as f64 / (1024.0 * 1024.0 * 1024.0)
    }

    pub fn total_gb(&self) -> f64 {
        self.total_system_bytes as f64 / (1024.0 * 1024.0 * 1024.0)
    }

    /// Computes optimal Bloom filter bits based on memory budget.
    pub fn optimal_bloom_bits(&self) -> usize {
        let budget_mb = self.max_budget_bytes / (1024 * 1024);
        if budget_mb < 2048 {
            // Under 2 GB budget: 64 MB per tier (256M bits)
            256 * 1024 * 1024
        } else if budget_mb < 8192 {
            // 2 - 8 GB budget: 128 MB per tier (512M bits) - standard default
            512 * 1024 * 1024
        } else if budget_mb < 32768 {
            // 8 - 32 GB budget: 256 MB per tier (1024M bits)
            1024 * 1024 * 1024
        } else {
            // > 32 GB server: 512 MB per tier (2048M bits)
            2048 * 1024 * 1024
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

fn detect_system_memory() -> (usize, usize) {
    // 1. Try /proc/meminfo (Linux / WSL)
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

    // 2. Safe fallback if /proc/meminfo not available (e.g. native Windows / macOS)
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
