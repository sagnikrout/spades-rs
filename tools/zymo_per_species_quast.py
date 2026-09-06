import csv
import glob
import os

ref_dir = 'data/zymo/references/ZymoBIOMICS.STD.refseq.v2/Genomes/'
files = sorted(glob.glob(os.path.join(ref_dir, '*.fasta')))

species_lens = {}
for f in files:
    sp_name = os.path.basename(f).replace('_complete_genome.fasta', '').replace('_draft_genome.fasta', '')
    bp = 0
    with open(f, 'r') as fh:
        for line in fh:
            if not line.startswith('>'):
                bp += len(line.strip())
    species_lens[sp_name] = bp

# Extract covered intervals per species
tsv_file = 'output/bench_zymo/quast_report/contigs_reports/all_alignments_contigs.tsv'
covered = {sp: {} for sp in species_lens}

with open(tsv_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 6 and parts[0].isdigit() and parts[1].isdigit():
            s = int(parts[0])
            e = int(parts[1])
            ref_name = parts[4] # e.g. Bacillus_subtilis_contig_1
            # extract species prefix
            sp = None
            for s_name in species_lens:
                if ref_name.startswith(s_name):
                    sp = s_name
                    break
            if sp:
                if ref_name not in covered[sp]:
                    covered[sp][ref_name] = []
                if s > e: s, e = e, s
                covered[sp][ref_name].append((s, e))

print("=" * 85)
print(f"{'Species':<28} | {'Ref Length':>12} | {'Assembled Bp':>14} | {'Recovery %':>10} | {'Status'}")
print("-" * 85)

total_ref = sum(species_lens.values())
total_cov = 0

bacterial_ref = 0
bacterial_cov = 0

for sp, ref_len in sorted(species_lens.items(), key=lambda x: x[1]):
    sp_cov = 0
    for ref_chr, intervals in covered[sp].items():
        merged = []
        for s, e in sorted(intervals):
            if not merged:
                merged.append([s, e])
            else:
                if s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
        sp_cov += sum(e - s + 1 for s, e in merged)
    
    total_cov += sp_cov
    pct = (sp_cov / ref_len) * 100.0
    
    is_yeast = "cerevisiae" in sp or "neoformans" in sp
    if not is_yeast:
        bacterial_ref += ref_len
        bacterial_cov += sp_cov
        
    status = "HIGH RECOVERY" if pct >= 70.0 else ("MODERATE" if pct >= 30.0 else "LOW ABUNDANCE (2%)")
    print(f"{sp:<28} | {ref_len:>10,d} bp | {sp_cov:>12,d} bp | {pct:>9.2f}% | {status}")

print("-" * 85)
bact_pct = (bacterial_cov / bacterial_ref) * 100.0
print(f"{'8 BACTERIAL SPECIES':<28} | {bacterial_ref:>10,d} bp | {bacterial_cov:>12,d} bp | {bact_pct:>9.2f}% | DOMINANT COMMUNITY")
total_pct = (total_cov / total_ref) * 100.0
print(f"{'COMBINED 10 SPECIES':<28} | {total_ref:>10,d} bp | {total_cov:>12,d} bp | {total_pct:>9.2f}% | TOTAL METAGENOME")
