# Contributing to spades-rs

Contributions are welcome. Please follow these guidelines to keep the codebase reliable, safe, and maintainable.

## Development environment

* Rust toolchain: 1.75 or later.
* Target architecture: 64-bit only (`x86_64`, `aarch64`). 32-bit targets are explicitly disallowed by compile-time assertions in `src/lib.rs`.

## Development workflow

### 1. Building the project

```bash
cargo build
```

For release builds with vectorization and optimizations:
```bash
cargo build --release
```

### 2. Running the test suite

All 53 unit and integration tests must pass:
```bash
cargo test --release
```

To run the real-world scenario validation suite across sequencing modalities and edge cases:
```bash
python3 tools/run_realworld_validation.py --quick
```

### 3. Formatting and linting

Code must be formatted and free of compiler and linter warnings:
```bash
cargo fmt --check
cargo clippy --all-targets -- -D warnings
```

## Architectural conventions

1. **Safety and dependencies**: Avoid adding external crate dependencies if standard library primitives or existing dependencies provide the required capability.
2. **Compact read representation**: Do not store millions of raw ASCII sequence strings in memory. Use the contiguous 2-bit packed read store (`src/packed_reads.rs`) with three-gate N-base handling.
3. **Memory governor**: Respect the dynamic memory governor in `src/memory.rs`. Allocations must remain bounded within the calculated hardware budget.
4. **Graph conventions**: Follow the bidirected port-involution topology in `src/graph.rs` and `src/spaligner.rs`.
