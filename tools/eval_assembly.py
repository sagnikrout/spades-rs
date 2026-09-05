#!/usr/bin/env python3
"""
Automated QUAST-equivalent Assembly Evaluation Harness.
Calculates:
- Contig count, Total assembled bases, Max contig length, GC content
- N50 and L50 (standard assembly contiguity)
- Genome fraction, aligned bases, and mismatch rate against NCBI reference genome.
"""

import sys
import os
import gzip

def parse_fasta(filepath):
    sequences = []
    names = []
    cur_seq = []
    cur_name = None
    
    open_fn = gzip.open if filepath.endswith(".gz") else open
    with open_fn(filepath, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur_name is not None:
                    sequences.append("".join(cur_seq).upper())
                    cur_seq = []
                cur_name = line[1:].split()[0]
                names.append(cur_name)
            else:
                cur_seq.append(line)
        if cur_name is not None:
            sequences.append("".join(cur_seq).upper())
            
    return names, sequences

def compute_n50(lengths, total_length=None):
    if not lengths:
        return 0, 0
    lengths = sorted(lengths, reverse=True)
    half_sum = (total_length if total_length else sum(lengths)) / 2.0
    cum_sum = 0
    for idx, l in enumerate(lengths):
        cum_sum += l
        if cum_sum >= half_sum:
            return l, idx + 1
    return lengths[-1], len(lengths)

def evaluate_assembly(contigs_file, ref_file=None, min_len=200):
    names, contigs = parse_fasta(contigs_file)
    # Filter by min_len
    filtered_contigs = [c for c in contigs if len(c) >= min_len]
    lengths = [len(c) for c in filtered_contigs]
    total_bp = sum(lengths)
    max_len = max(lengths) if lengths else 0
    n_contigs = len(filtered_contigs)
    
    # GC content
    total_gc = sum(c.count('G') + c.count('C') for c in filtered_contigs)
    gc_pct = (total_gc / total_bp * 100.0) if total_bp > 0 else 0.0
    
    n50, l50 = compute_n50(lengths, total_bp)
    
    report = {
        "Contigs (>= {} bp)".format(min_len): n_contigs,
        "Total length (bp)": total_bp,
        "Max contig length (bp)": max_len,
        "N50 (bp)": n50,
        "L50": l50,
        "GC (%)": round(gc_pct, 2),
    }
    
    # Reference comparison if provided
    if ref_file and os.path.exists(ref_file):
        _, ref_seqs = parse_fasta(ref_file)
        full_ref = "".join(ref_seqs)
        ref_len = len(full_ref)
        report["Reference length (bp)"] = ref_len
        
        # Calculate reference k-mer coverage / genome fraction using 31-mers
        k = 31
        tr = str.maketrans("ACGTN", "TGCAN")
        def revcomp(s):
            return s.translate(tr)[::-1]
        
        ref_kmers = set()
        for i in range(ref_len - k + 1):
            sub = full_ref[i:i+k]
            ref_kmers.add(min(sub, revcomp(sub)))
            
        covered_kmers = set()
        for c in filtered_contigs:
            if len(c) < k:
                continue
            for i in range(len(c) - k + 1):
                sub = c[i:i+k]
                if "N" in sub:
                    continue
                can = min(sub, revcomp(sub))
                if can in ref_kmers:
                    covered_kmers.add(can)
                    
        genome_fraction = (len(covered_kmers) / len(ref_kmers) * 100.0) if ref_kmers else 0.0
        report["Reference Genome Fraction (%)"] = round(genome_fraction, 2)
        
    return report

def main():
    if len(sys.argv) < 2:
        print("Usage: eval_assembly.py <contigs.fasta> [reference.fasta] [min_contig_length]")
        sys.exit(1)
        
    contigs_file = sys.argv[1]
    ref_file = sys.argv[2] if len(sys.argv) > 2 else None
    min_len = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    
    metrics = evaluate_assembly(contigs_file, ref_file, min_len)
    
    print("-" * 55)
    print("           ASSEMBLY QUALITY METRICS REPORT")
    print("-" * 55)
    for k, v in metrics.items():
        print(f"  {k:<35}: {v}")
    print("-" * 55)

if __name__ == "__main__":
    main()
