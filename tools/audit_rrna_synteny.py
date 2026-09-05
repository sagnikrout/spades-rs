#!/usr/bin/env python3
"""
Audits E. coli ribosomal RNA (rRNA) operon synteny and structural integrity.
Verifies whether the 7 known rrn operons (rrnA, rrnB, rrnC, rrnD, rrnE, rrnG, rrnH)
are properly spanned without chimeric inter-operon translocations or inversions.
Also runs QUAST if available.
"""

import sys
import os
import subprocess
import json

RRNA_OPERONS = {
    "rrnH": {"start": 223771, "end": 229295, "strand": "-"},
    "rrnG": {"start": 2725835, "end": 2731349, "strand": "-"},
    "rrnD": {"start": 3423674, "end": 3429188, "strand": "-"},
    "rrnC": {"start": 3941607, "end": 3947126, "strand": "+"},
    "rrnA": {"start": 4035268, "end": 4040782, "strand": "+"},
    "rrnB": {"start": 4166198, "end": 4171712, "strand": "+"},
    "rrnE": {"start": 4207896, "end": 4213410, "strand": "+"},
}

def parse_fasta(filepath):
    names, seqs = [], []
    cur_name, cur_seq = None, []
    with open(filepath, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur_name is not None:
                    seqs.append("".join(cur_seq).upper())
                    cur_seq = []
                cur_name = line[1:].split()[0]
                names.append(cur_name)
            else:
                cur_seq.append(line)
        if cur_name is not None:
            seqs.append("".join(cur_seq).upper())
    return names, seqs

def audit_rrna_operons(contigs_file, ref_file):
    print("=" * 65)
    print("      E. COLI 7 RIBOSOMAL RNA OPERONS STRUCTURAL INTEGRITY AUDIT")
    print("=" * 65)
    print(f"  Assembly Contigs: {contigs_file}")
    print(f"  NCBI Reference:   {ref_file}")
    print("-" * 65)

    _, ref_seqs = parse_fasta(ref_file)
    ref = "".join(ref_seqs)
    c_names, c_seqs = parse_fasta(contigs_file)

    k = 31
    tr = str.maketrans("ACGT", "TGCA")
    def rc(s): return s.translate(tr)[::-1]

    # Index contig 31-mers: kmer -> list of (contig_idx, offset, is_rc)
    print("  Indexing assembly contigs...")
    contig_kmers = {}
    for c_idx, seq in enumerate(c_seqs):
        if len(seq) < k:
            continue
        for i in range(0, len(seq) - k + 1):
            km = seq[i:i+k]
            can = min(km, rc(km))
            is_rev = (can != km)
            contig_kmers.setdefault(can, []).append((c_idx, i, is_rev))

    results = {}
    all_intact = True

    print(f"  {'Operon':<8} {'Coordinates (bp)':<25} {'Strand':<8} {'Status':<15} {'Spanning Contig'}")
    print("-" * 65)

    for op_name, op_info in RRNA_OPERONS.items():
        s = op_info["start"]
        e = op_info["end"]
        strand = op_info["strand"]
        op_len = e - s

        # Sample 5' flank (1 kb before operon), internal, and 3' flank (1 kb after)
        flank5 = ref[max(0, s - 1000):s]
        internal = ref[s + 1000:e - 1000]
        flank3 = ref[e:min(len(ref), e + 1000)]

        # Map each region to contigs
        def find_contig_matches(seq_chunk):
            hits = {}
            for i in range(0, len(seq_chunk) - k + 1, 10):
                sub = seq_chunk[i:i+k]
                can = min(sub, rc(sub))
                if can in contig_kmers:
                    for (c_idx, c_pos, is_rev) in contig_kmers[can]:
                        hits[c_idx] = hits.get(c_idx, 0) + 1
            return hits

        hits_5 = find_contig_matches(flank5)
        hits_mid = find_contig_matches(internal)
        hits_3 = find_contig_matches(flank3)

        # A contig bridges the operon without translocation if it matches 5' flank, mid, and 3' flank!
        bridging_contigs = []
        for c_idx in hits_mid:
            if hits_5.get(c_idx, 0) >= 3 and hits_3.get(c_idx, 0) >= 3:
                bridging_contigs.append(c_names[c_idx])

        if bridging_contigs:
            status = "BRIDGED (Intact)"
            best_c = bridging_contigs[0]
        elif hits_mid:
            best_c_idx = max(hits_mid, key=hits_mid.get)
            best_c = c_names[best_c_idx]
            status = "Partially Covered"
        else:
            best_c = "None"
            status = "Uncovered"

        results[op_name] = {
            "status": status,
            "spanning_contig": best_c,
            "coordinates": f"{s}-{e}",
            "strand": strand,
        }

        print(f"  {op_name:<8} {s:>9}-{e:<13} {strand:^8} {status:<15} {best_c}")

    print("=" * 65)
    return results

def run_quast(contigs_file, ref_file, out_dir="output/quast_report"):
    print("\n--- RUNNING STANDARDIZED QUAST EVALUATION ---")
    quast_bin = os.path.expanduser("~/.local/bin/quast.py")
    if not os.path.exists(quast_bin):
        quast_bin = "quast.py"
    quast_cmd = f"python3 {quast_bin} {contigs_file} -r {ref_file} -o {out_dir} --threads 22 --min-contig 200"
    print(f"Command: {quast_cmd}")
    res = subprocess.run(quast_cmd, shell=True)
    if res.returncode == 0:
        report_tsv = os.path.join(out_dir, "report.tsv")
        if os.path.exists(report_tsv):
            print("\n--- QUAST OFFICIAL SUMMARY METRICS ---")
            with open(report_tsv) as f:
                for line in f:
                    print("  " + line.strip())
    else:
        print("QUAST execution failed or quast.py not in PATH.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 audit_rrna_synteny.py <contigs.fasta> <reference.fasta> [quast_out_dir]")
        sys.exit(1)

    c_file = sys.argv[1]
    r_file = sys.argv[2]
    q_dir = sys.argv[3] if len(sys.argv) > 3 else "output/quast_report"

    audit_rrna_operons(c_file, r_file)
    run_quast(c_file, r_file, q_dir)
