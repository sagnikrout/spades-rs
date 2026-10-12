import gzip
import os
import sys

def load_fasta_kmers(fa_path, k=31):
    kmers = set()
    with open(fa_path, 'r') as f:
        seq = []
        for line in f:
            if line.startswith('>'):
                s = "".join(seq).upper()
                for i in range(len(s) - k + 1):
                    kmers.add(s[i:i+k])
                seq = []
            else:
                seq.append(line.strip())
        s = "".join(seq).upper()
        for i in range(len(s) - k + 1):
            kmers.add(s[i:i+k])
    return kmers

def revcomp(s):
    tr = str.maketrans('ACGTN', 'TGCAN')
    return s.translate(tr)[::-1]

ref_path = sys.argv[1] if len(sys.argv) > 1 else 'data/pombe/pombe_ref.fa'
fq_path = sys.argv[2] if len(sys.argv) > 2 else 'data/pombe/pombe_R1.fq.gz'

if not os.path.exists(ref_path):
    print(f"Error: Reference file '{ref_path}' not found.")
    sys.exit(1)
if not os.path.exists(fq_path):
    print(f"Error: FASTQ file '{fq_path}' not found.")
    sys.exit(1)

print("Loading reference kmers...")
ref_kmers = load_fasta_kmers(ref_path, k=31)
print(f"Loaded {len(ref_kmers):,} 31-mers from S. pombe reference.")

tested = 0
matched = 0
open_fn = gzip.open if fq_path.endswith('.gz') else open
with open_fn(fq_path, 'rt') as f:
    for i, line in enumerate(f):
        if i % 4 == 1:
            read = line.strip().upper()
            tested += 1
            # check any 31-mer in read
            found = False
            for j in range(0, len(read) - 31 + 1, 10):
                km = read[j:j+31]
                rc = revcomp(km)
                if km in ref_kmers or rc in ref_kmers:
                    found = True
                    break
            if found:
                matched += 1
            if tested >= 2000:
                break

pct = (matched / max(1, tested)) * 100
print(f"Sample verification: {matched}/{tested} reads ({pct:.2f}%) match S. pombe reference!")
assert pct > 90.0, f"Low identity: {pct}%"
print("IDENTITY VERIFIED: Authentic Schizosaccharomyces pombe reads!")
