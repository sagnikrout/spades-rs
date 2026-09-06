import glob
import os

ref_dir = 'data/zymo/references/ZymoBIOMICS.STD.refseq.v2/Genomes/'
files = glob.glob(os.path.join(ref_dir, '*.fasta'))

species = {}
for f in files:
    name = os.path.basename(f).replace('_complete_genome.fasta', '').replace('_draft_genome.fasta', '')
    seqs = []
    with open(f, 'r') as fh:
        for line in fh:
            if not line.startswith('>'):
                seqs.append(line.strip().upper())
    species[name] = "".join(seqs)
    print(f"Loaded {name}: {len(species[name]):,} bp")

def revcomp(s):
    tr = str.maketrans('ACGTN', 'TGCAN')
    return s.translate(tr)[::-1]

# Test reads
r1 = "AGCCTGCAAAAAAACTTGGGAAAAGTGCCTCTGCACGAACTTTTCCAATCGCCAAAACTTCTAATCCGAAAACCGTTCCCGCTAATGGTGTGCCAAAAACAGAACTAAAACCTGCACTAATTCCGCTAATAATAATTACTTGTCGTTCCAA"
r2 = "TGGTTGGTCCGCTGATCGACGTTATGGGCTCTGCGGGTGAAGACCTGAAACGCCAGCAGGCGCAGGTTGAGCAGGTGCTGAAGACTGAAGAAGAGCAGTTTGCTCGTACTCTGGAGCGCGGTCTGGCGTTGCTGGATGAAGAGCTGGCAAA"

for r_name, r in [("Read 1", r1), ("Read 2", r2)]:
    print(f"\nMatching {r_name}:")
    for name, s in species.items():
        if r in s:
            print(f"  -> EXACT FORWARD in {name}!")
        elif revcomp(r) in s:
            print(f"  -> EXACT REVERSE in {name}!")
        elif r[:50] in s or revcomp(r[:50]) in s:
            print(f"  -> Prefix 50bp match in {name}!")
