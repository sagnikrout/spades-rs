# Contributing to spades-rs

Contributions are welcome. Please follow these guidelines to keep the codebase fast, safe, and maintainable.

## Development environment

* Rust toolchain: 1.75 or later.
* Target architecture: 64-bit only (`x86_64`, `aarch64`). 32-bit targets are explicitly disallowed via compile-time assertions in `src/lib.rs`.

## Development workflow

### 1. Building the project
```bash
cargo build
```

For release builds with SIMD vectorization and Link-Time Optimization (LTO):
```bash
cargo build --release
```

### 2. Running tests
All 39 unit and integration tests must pass:
```bash
cargo test
```

### 3. Formatting and linting
Format code before opening a pull request:
```bash
cargo fmt --check
cargo clippy --all-targets -- -D warnings
```

## Architectural guidelines

1. **Safety and Dependencies**: The codebase is 100% safe Rust (`#![forbid(unsafe_code)]` where possible). Avoid adding external dependencies if the standard library or existing crates already provide the functionality.
2. **Memory Efficiency**: Never hold raw ASCII string sequences in memory for large datasets. Use the 2-bit packed storage (`src/packed_reads.rs`) and Two-Tier Bloom filters (`src/bloom.rs`).
3. **Memory Governor**: Respect the dynamic memory governor in `src/memory.rs`. Long-lived allocations should remain within the calculated hardware budget.
4. **Graph Algorithms**: Follow the bidirected port-involution conventions in `src/graph.rs` and `src/spaligner.rs`.
