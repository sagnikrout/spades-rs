from verify_pombe_identity import load_fasta_kmers, revcomp
import gzip
import urllib.request
import io

print("Loading E. coli reference kmers...")
ref_kmers = load_fasta_kmers("data/ecoli_mg1655/ref_NC_000913.3.fa", k=31)
print(f"Loaded {len(ref_kmers):,} 31-mers from E. coli reference.")

print("Fetching initial chunk from SRR31677630...")
req = urllib.request.Request(
    "http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR316/030/SRR31677630/SRR31677630_1.fastq.gz",
    headers={"Range": "bytes=0-500000"}
)
with urllib.request.urlopen(req) as resp:
    data = resp.read()

f = gzip.GzipFile(fileobj=io.BytesIO(data))
matched = 0
total = 0
for i, line in enumerate(f):
    if i % 4 == 1:
        read = line.decode().strip().upper()
        total += 1
        found = any(read[j:j+31] in ref_kmers or revcomp(read[j:j+31]) in ref_kmers for j in range(0, len(read)-31+1, 10))
        if found:
            matched += 1
        if total >= 500:
            break

pct = (matched / total) * 100.0
print(f"Single-cell identity verification: {matched}/{total} ({pct:.2f}%) match E. coli reference!")
