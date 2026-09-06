import gzip
import glob
import os
import sys

ref_dir = 'data/zymo/references/ZymoBIOMICS.STD.refseq.v2/Genomes/'
files = sorted(glob.glob(os.path.join(ref_dir, '*.fasta')))

print("Building reference 31-mer tables for all 10 species...")
species_kmers = {}
for f in files:
    sp_name = os.path.basename(f).replace('_complete_genome.fasta', '').replace('_draft_genome.fasta', '')
    km_set = set()
    with open(f, 'r') as fh:
        seq = []
        for line in fh:
            if line.startswith('>'):
                s = "".join(seq).upper()
                for i in range(len(s) - 31 + 1):
                    km_set.add(s[i:i+31])
                seq = []
            else:
                seq.append(line.strip())
        s = "".join(seq).upper()
        for i in range(len(s) - 31 + 1):
            km_set.add(s[i:i+31])
    species_kmers[sp_name] = km_set
    print(f"  {sp_name:25s}: {len(km_set):10,d} 31-mers")

def revcomp(s):
    tr = str.maketrans('ACGTN', 'TGCAN')
    return s.translate(tr)[::-1]

print("\nSampling 10,000 raw reads from data/zymo/zymo_R1.fq.gz...")
counts = {sp: 0 for sp in species_kmers}
unassigned = 0
total_tested = 0

with gzip.open('data/zymo/zymo_R1.fq.gz', 'rt') as f:
    for i, line in enumerate(f):
        if i % 4 == 1:
            read = line.strip().upper()
            total_tested += 1
            assigned = False
            for j in range(0, len(read) - 31 + 1, 10):
                km = read[j:j+31]
                rc = revcomp(km)
                for sp, km_set in species_kmers.items():
                    if km in km_set or rc in km_set:
                        counts[sp] += 1
                        assigned = True
                        break
                if assigned:
                    break
            if not assigned:
                unassigned += 1
            if total_tested >= 10000:
                break

print("=" * 65)
print(f"{'Species':<28} | {'Read Count':>12} | {'Percentage':>10}")
print("-" * 65)
for sp, cnt in sorted(counts.items(), key=lambda x: x[1], reverse=True):
    pct = (cnt / total_tested) * 100.0
    print(f"{sp:<28} | {cnt:>12,d} | {pct:>9.2f}%")
print("-" * 65)
unass_pct = (unassigned / total_tested) * 100.0
print(f"{'Unassigned / Ambiguous':<28} | {unassigned:>12,d} | {unass_pct:>9.2f}%")
print(f"{'TOTAL TESTED':<28} | {total_tested:>12,d} | 100.00%")
print("=" * 65)

# Verify all 10 species are present
present = sum(1 for cnt in counts.values() if cnt > 10)
print(f"Species detected (>10 reads): {present}/10 species!")
if present == 10:
    print("COMMUNITY IDENTITY CONFIRMED: Authentic ZymoBIOMICS 10-species mock community!")
