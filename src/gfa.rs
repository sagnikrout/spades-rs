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
    if k < 2 {
        anyhow::bail!("k must be at least 2 for GFA link generation, got {}", k);
    }
    let k1 = k - 1;

    // 1. Header
    writeln!(writer, "H\tVN:Z:1.1")?;

    // 2. Segment lines: S <id> <seq> LN:i:<len> cv:f:<coverage> RC:f:<coverage>
    for u in unitigs {
        let seq_str = String::from_utf8_lossy(&u.sequence);
        writeln!(
            writer,
            "S\t{}\t{}\tLN:i:{}\tcv:f:{:.1}\tRC:f:{:.1}",
            u.id,
            seq_str,
            u.sequence.len(),
            u.mean_coverage,
            u.mean_coverage
        )?;
    }

    // 3. Link lines: L <from> <strand> <to> <strand> <k-1>M
    // Index forward and reverse-complement prefixes of length k1
    // (prefix_seq, target_id, is_target_rc)
    let mut prefix_map: hashbrown::HashMap<Vec<u8>, Vec<(usize, char)>> = hashbrown::HashMap::new();
    for u in unitigs {
        if u.sequence.len() >= k1 {
            // Forward prefix: entering target with '+'
            prefix_map
                .entry(u.sequence[..k1].to_vec())
                .or_default()
                .push((u.id, '+'));
            // Reverse-complement prefix (i.e. RC of suffix): entering target with '-'
            let rc_prefix = crate::dna::revcomp_bytes(&u.sequence[u.sequence.len() - k1..]);
            prefix_map.entry(rc_prefix).or_default().push((u.id, '-'));
        }
    }

    for u_i in unitigs {
        if u_i.sequence.len() < k1 {
            continue;
        }
        // Outgoing from forward 3' end ('+')
        let suffix_fwd = &u_i.sequence[u_i.sequence.len() - k1..];
        if let Some(targets) = prefix_map.get(suffix_fwd) {
            for &(target_id, target_strand) in targets {
                writeln!(
                    writer,
                    "L\t{}\t+\t{}\t{}\t{}M",
                    u_i.id, target_id, target_strand, k1
                )?;
            }
        }

        // Outgoing from reverse-complement 3' end ('-') = RC of forward prefix
        let suffix_rev = crate::dna::revcomp_bytes(&u_i.sequence[..k1]);
        if let Some(targets) = prefix_map.get(&suffix_rev) {
            for &(target_id, target_strand) in targets {
                // To avoid duplicate mirror links, only emit if target_id != u_i.id or target_strand == '+'
                if target_id != u_i.id || target_strand == '+' {
                    writeln!(
                        writer,
                        "L\t{}\t-\t{}\t{}\t{}M",
                        u_i.id, target_id, target_strand, k1
                    )?;
                }
            }
        }
    }

    writer.flush()?;
    Ok(())
}
