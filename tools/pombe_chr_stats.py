import csv
from audit_pombe_loci import load_fasta

ref = load_fasta('data/pombe/pombe_ref.fa')
chr_lens = {k: len(v) for k, v in ref.items()}

covered_intervals = {k: [] for k in chr_lens}

tsv_file = 'output/bench_pombe/quast_report/contigs_reports/all_alignments_contigs.tsv'
with open(tsv_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 6 and parts[0].isdigit() and parts[1].isdigit():
            s = int(parts[0])
            e = int(parts[1])
            ref_name = parts[4]
            if s > e:
                s, e = e, s
            if ref_name in covered_intervals:
                covered_intervals[ref_name].append((s, e))

print("=" * 75)
print(f"{'Chromosome':<25} | {'Ref Length':>14} | {'Aligned Bp':>14} | {'Coverage %':>10}")
print("-" * 75)

total_ref = sum(chr_lens.values())
total_aligned = 0

for chrom, length in chr_lens.items():
    intervals = sorted(covered_intervals[chrom])
    merged = []
    for s, e in intervals:
        if not merged:
            merged.append([s, e])
        else:
            if s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
    aligned = sum(e - s + 1 for s, e in merged)
    total_aligned += aligned
    pct = (aligned / length) * 100.0
    label = chrom
    if 'NC_003424' in chrom: label = "Chromosome I (5.58 Mb)"
    elif 'NC_003423' in chrom: label = "Chromosome II (4.54 Mb)"
    elif 'NC_003421' in chrom: label = "Chromosome III (2.45 Mb)"
    elif 'NC_001326' in chrom: label = "Mitochondrion (19.4 kb)"
    print(f"{label:<25} | {length:>12,d} bp | {aligned:>12,d} bp | {pct:>9.2f}%")

print("-" * 75)
total_pct = (total_aligned / total_ref) * 100.0
print(f"{'TOTAL GENOME':<25} | {total_ref:>12,d} bp | {total_aligned:>12,d} bp | {total_pct:>9.2f}%")
