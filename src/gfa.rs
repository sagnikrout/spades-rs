//! Graphical Fragment Assembly (GFA v1.1) Exporter (Module 4.3).
//!
//! Exports assembly graph topology and unitig sequences into standard GFA format
//! compatible with Bandage, IGV, and other graph genome viewers.

use crate::graph::Unitig;
use anyhow::Result;
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// Exports assembly graph unitigs and links into GFA v1.1 format.
pub fn write_graph_gfa<P: AsRef<Path>>(
    k: usize,
    unitigs: &[Unitig],
    out_path: P,
) -> Result<()> {
    let mut file = File::create(out_path)?;
    let k1 = k - 1;

    // 1. Header
    writeln!(file, "H\tVN:Z:1.0")?;

    // 2. Segment lines: S <id> <seq> LN:i:<len> RC:f:<coverage>
    for u in unitigs {
        let seq_str = String::from_utf8_lossy(&u.sequence);
        writeln!(
            file,
            "S\t{}\t{}\tLN:i:{}\tRC:f:{:.1}",
            u.id,
            seq_str,
            u.sequence.len(),
            u.mean_coverage
        )?;
    }

    // 3. Link lines: L <from> + <to> + <k-1>M
    for (i, u_i) in unitigs.iter().enumerate() {
        if u_i.sequence.len() < k1 {
            continue;
        }
        let suffix_i = &u_i.sequence[u_i.sequence.len() - k1..];

        for (j, u_j) in unitigs.iter().enumerate() {
            if i == j || u_j.sequence.len() < k1 {
                continue;
            }
            if suffix_i == &u_j.sequence[..k1] {
                writeln!(
                    file,
                    "L\t{}\t+\t{}\t+\t{}M",
                    u_i.id,
                    u_j.id,
                    k1
                )?;
            }
        }
    }

    Ok(())
}
