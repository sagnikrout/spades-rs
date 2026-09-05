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
pub fn write_graph_gfa<P: AsRef<Path>>(k: usize, unitigs: &[Unitig], out_path: P) -> Result<()> {
    let file = File::create(out_path)?;
    let mut writer = std::io::BufWriter::with_capacity(1024 * 1024, file);
    let k1 = k - 1;

    // 1. Header
    writeln!(writer, "H\tVN:Z:1.0")?;

    // 2. Segment lines: S <id> <seq> LN:i:<len> RC:f:<coverage>
    for u in unitigs {
        let seq_str = String::from_utf8_lossy(&u.sequence);
        writeln!(
            writer,
            "S\t{}\t{}\tLN:i:{}\tRC:f:{:.1}",
            u.id,
            seq_str,
            u.sequence.len(),
            u.mean_coverage
        )?;
    }

    // 3. Link lines: L <from> + <to> + <k-1>M
    // Index prefixes of length k1
    let mut prefix_map: hashbrown::HashMap<&[u8], Vec<usize>> = hashbrown::HashMap::new();
    for u in unitigs {
        if u.sequence.len() >= k1 {
            prefix_map.entry(&u.sequence[..k1]).or_default().push(u.id);
        }
    }

    for u_i in unitigs {
        if u_i.sequence.len() < k1 {
            continue;
        }
        let suffix_i = &u_i.sequence[u_i.sequence.len() - k1..];
        if let Some(target_ids) = prefix_map.get(suffix_i) {
            for &target_id in target_ids {
                if target_id != u_i.id {
                    writeln!(writer, "L\t{}\t+\t{}\t+\t{}M", u_i.id, target_id, k1)?;
                }
            }
        }
    }

    writer.flush()?;
    Ok(())
}
