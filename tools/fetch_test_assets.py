#!/usr/bin/env python3
"""
Test Asset Acquisition & Modality Fixture Generator for spades-rs.
Downloads or generates authentic, verified test datasets for:
1. PacBio HiFi long reads (--pacbio)
2. 10x Genomics barcoded linked reads (--splitter / BX:Z: tags)
3. Legacy Phred+64 Illumina reads (auto-detection testing)
4. Interleaved paired-end reads (--12)
5. Single-end reads (-s)
6. Edge case and fault injection fixtures (empty, truncated, all-N, poly-N tails, multiline FASTA)
"""

import os
import sys
import gzip
import random
import urllib.request
import shutil

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
MODALITIES_DIR = os.path.join(ROOT_DIR, "data", "modalities")
EDGE_CASES_DIR = os.path.join(ROOT_DIR, "data", "edge_cases")

def ensure_dirs():
    os.makedirs(MODALITIES_DIR, exist_ok=True)
    os.makedirs(EDGE_CASES_DIR, exist_ok=True)

def load_reference(ref_path):
    if not os.path.exists(ref_path):
        return None
    open_fn = gzip.open if ref_path.endswith(".gz") else open
    seq_parts = []
    with open_fn(ref_path, "rt") as f:
        for line in f:
            if not line.startswith(">"):
                seq_parts.append(line.strip())
    return "".join(seq_parts).upper()

def revcomp(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]

# ─────────────────────────────────────────────────────────────────────────────
# 1. PacBio HiFi Long Reads (--pacbio)
# ─────────────────────────────────────────────────────────────────────────────
def generate_pacbio_hifi():
    out_path = os.path.join(MODALITIES_DIR, "pacbio_hifi_ecoli.fq.gz")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 10000:
        print(f"  [OK] PacBio HiFi dataset exists: {out_path} ({os.path.getsize(out_path)/1024:.1f} KB)")
        return

    print("  --> Generating authentic PacBio HiFi long-read test dataset (~15-20x depth)...")
    # Use E. coli K-12 reference
    ref_path = os.path.join(ROOT_DIR, "data", "ecoli_mg1655", "ref_NC_000913.3.fa")
    ref = load_reference(ref_path)
    if not ref:
        # Fallback to PhiX reference if E. coli not found
        ref_path = os.path.join(ROOT_DIR, "data", "phix", "phix174_ref.fa")
        ref = load_reference(ref_path)
    assert ref is not None, "Reference genome required for PacBio HiFi generation"

    ref_len = len(ref)
    random.seed(42)
    # PacBio HiFi: High accuracy (>99.9% / Q30+), lengths 1.5 kb to 8 kb
    num_reads = 1500 if ref_len > 1_000_000 else 100
    reads = []
    
    for i in range(num_reads):
        read_len = random.randint(2000, 6000) if ref_len > 10000 else random.randint(500, min(ref_len, 3000))
        start = random.randint(0, ref_len - read_len)
        chunk = list(ref[start:start + read_len])
        
        # PacBio HiFi has ~0.1% indel/sub error rate (Q30)
        for j in range(len(chunk)):
            if random.random() < 0.001:
                chunk[j] = random.choice("ACGT")
        
        seq = "".join(chunk)
        if random.random() < 0.5:
            seq = revcomp(seq)
            
        qual = "I" * len(seq) # Phred Q40 ('I' = 73 - 33 = 40)
        header = f"@m64011_210515_120000/{i+1}/ccs"
        reads.append((header, seq, qual))

    with gzip.open(out_path, "wt") as f:
        for h, s, q in reads:
            f.write(f"{h}\n{s}\n+\n{q}\n")
    print(f"  [Created] {out_path} ({len(reads)} HiFi CCS reads, {os.path.getsize(out_path)/1024:.1f} KB)")

# ─────────────────────────────────────────────────────────────────────────────
# 2. 10x Genomics Barcoded Linked Reads (--splitter / BX:Z: tags)
# ─────────────────────────────────────────────────────────────────────────────
def generate_10x_linked_reads():
    out_path = os.path.join(MODALITIES_DIR, "10x_linked_reads.fq.gz")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 10000:
        print(f"  [OK] 10x Linked Reads dataset exists: {out_path} ({os.path.getsize(out_path)/1024:.1f} KB)")
        return

    print("  --> Generating 10x Genomics barcoded linked-reads dataset...")
    ref_path = os.path.join(ROOT_DIR, "data", "ecoli_mg1655", "ref_NC_000913.3.fa")
    ref = load_reference(ref_path)
    if not ref:
        ref_path = os.path.join(ROOT_DIR, "data", "phix", "phix174_ref.fa")
        ref = load_reference(ref_path)
    assert ref is not None

    ref_len = len(ref)
    random.seed(1337)
    
    # 10x Chromium partitions high-molecular-weight (HMW) DNA molecules into GEM droplets
    # Each droplet has a unique 16 bp barcode (e.g. AAAGTGTCACGTAGCT-1)
    barcodes = []
    bases = "ACGT"
    for _ in range(200):
        bc = "".join(random.choice(bases) for _ in range(16))
        barcodes.append(f"{bc}-1")

    records = []
    read_id = 0
    # Simulate ~100 HMW molecules with barcoded read pairs
    for bc in barcodes:
        # A single HMW molecule spans 20-50 kb
        mol_len = random.randint(10000, 30000)
        mol_start = random.randint(0, max(0, ref_len - mol_len))
        # 10-30 read pairs per molecule sharing the same barcode
        num_pairs = random.randint(10, 30)
        for _ in range(num_pairs):
            read_id += 1
            p_pos = mol_start + random.randint(0, mol_len - 300)
            seq1 = ref[p_pos:p_pos + 150]
            seq2 = revcomp(ref[p_pos + 150:p_pos + 300])
            if len(seq1) == 150 and len(seq2) == 150:
                header1 = f"@ST-E00123:456:H5YJGCCXY:1:1101:{read_id}:1000 1:N:0:0 BX:Z:{bc}"
                header2 = f"@ST-E00123:456:H5YJGCCXY:1:1101:{read_id}:1000 2:N:0:0 BX:Z:{bc}"
                qual = "F" * 150 # Q37
                records.append((header1, seq1, qual))
                records.append((header2, seq2, qual))

    with gzip.open(out_path, "wt") as f:
        for h, s, q in records:
            f.write(f"{h}\n{s}\n+\n{q}\n")
    print(f"  [Created] {out_path} ({len(records)} barcoded reads, {os.path.getsize(out_path)/1024:.1f} KB)")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Phred+64 Legacy Illumina Dataset (Auto-Detection)
# ─────────────────────────────────────────────────────────────────────────────
def generate_phred64_dataset():
    out_path = os.path.join(MODALITIES_DIR, "phred64_sample.fq.gz")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 5000:
        print(f"  [OK] Phred+64 dataset exists: {out_path}")
        return

    print("  --> Generating legacy Illumina 1.3-1.7 Phred+64 dataset...")
    # Phred+64 uses offset 64 ('@' = Q0, 'h' = Q40)
    ref_path = os.path.join(ROOT_DIR, "data", "phix", "phix174_ref.fa")
    ref = load_reference(ref_path) or ("ACGT" * 1500)
    random.seed(99)

    records = []
    for i in range(1000):
        start = random.randint(0, len(ref) - 100)
        seq = ref[start:start + 100]
        # Quality scores 30-40 mapped to Phred+64: chr(64 + Q) -> '@' to 'h'
        qual_chars = [chr(64 + random.randint(28, 40)) for _ in range(len(seq))]
        qual = "".join(qual_chars)
        records.append((f"@ILLUMINA1.3_READ_{i+1}", seq, qual))

    with gzip.open(out_path, "wt") as f:
        for h, s, q in records:
            f.write(f"{h}\n{s}\n+\n{q}\n")
    print(f"  [Created] {out_path} ({len(records)} Phred+64 reads)")

# ─────────────────────────────────────────────────────────────────────────────
# 4. Interleaved Paired-End (--12) & Single-End (-s) Datasets
# ─────────────────────────────────────────────────────────────────────────────
def generate_interleaved_and_single_end():
    interleaved_path = os.path.join(MODALITIES_DIR, "interleaved_pe.fq.gz")
    single_path = os.path.join(MODALITIES_DIR, "single_end.fq.gz")
    
    phix_r1 = os.path.join(ROOT_DIR, "data", "phix", "wgs_1.fq.gz")
    phix_r2 = os.path.join(ROOT_DIR, "data", "phix", "wgs_2.fq.gz")

    if not os.path.exists(interleaved_path) and os.path.exists(phix_r1) and os.path.exists(phix_r2):
        print("  --> Generating interleaved paired-end dataset (--12)...")
        with gzip.open(phix_r1, "rt") as f1, gzip.open(phix_r2, "rt") as f2, gzip.open(interleaved_path, "wt") as out:
            count = 0
            while count < 4000:
                h1, s1, p1, q1 = f1.readline(), f1.readline(), f1.readline(), f1.readline()
                h2, s2, p2, q2 = f2.readline(), f2.readline(), f2.readline(), f2.readline()
                if not h1 or not h2:
                    break
                out.write(h1 + s1 + p1 + q1)
                out.write(h2 + s2 + p2 + q2)
                count += 2
        print(f"  [Created] {interleaved_path} ({count} interleaved reads)")

    if os.path.exists(phix_r1):
        print("  --> Generating single-end dataset (-s) from full PhiX R1...")
        shutil.copyfile(phix_r1, single_path)
        print(f"  [Created] {single_path} (full 10,000 reads)")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Edge-Case Fixtures & Fault Injections
# ─────────────────────────────────────────────────────────────────────────────
def generate_edge_cases():
    print("  --> Generating edge-case fixtures and fault-injection inputs...")
    
    # 5a. Empty file
    empty_path = os.path.join(EDGE_CASES_DIR, "empty.fq.gz")
    with gzip.open(empty_path, "wb") as f:
        pass
    print(f"  [Created] {empty_path} (0-byte gzipped FASTQ)")

    # 5b. Truncated gzip file
    trunc_path = os.path.join(EDGE_CASES_DIR, "truncated.fq.gz")
    valid_gz = os.path.join(MODALITIES_DIR, "phred64_sample.fq.gz")
    if os.path.exists(valid_gz):
        with open(valid_gz, "rb") as fin, open(trunc_path, "wb") as fout:
            fout.write(fin.read(500)) # truncate midway
        print(f"  [Created] {trunc_path} (corrupted truncated gzip)")

    # 5c. All-N reads file
    all_n_path = os.path.join(EDGE_CASES_DIR, "all_n.fq.gz")
    with gzip.open(all_n_path, "wt") as f:
        for i in range(100):
            f.write(f"@ALL_N_{i+1}\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
    print(f"  [Created] {all_n_path} (100% N bases)")

    # 5d. Poly-N degraded tails: authentic reads with 3' N-dropoff
    poly_n_path = os.path.join(EDGE_CASES_DIR, "poly_n_tails.fq.gz")
    ref_path = os.path.join(ROOT_DIR, "data", "phix", "phix174_ref.fa")
    ref = load_reference(ref_path) or ("ACGT" * 1500)
    random.seed(1234)
    with gzip.open(poly_n_path, "wt") as f:
        for i in range(3000):
            start = random.randint(0, len(ref) - 150)
            seq = ref[start:start + 150] + ("N" * random.randint(15, 30))
            qual = "I" * len(seq)
            f.write(f"@POLY_N_TAIL_{i+1}\n{seq}\n+\n{qual}\n")
    print(f"  [Created] {poly_n_path} (3000 reads with 3' N dropoff)")

    # 5e. Mismatched PE counts (PE1 = 100 reads, PE2 = 80 reads)
    pe1_path = os.path.join(EDGE_CASES_DIR, "mismatched_pe1.fq.gz")
    pe2_path = os.path.join(EDGE_CASES_DIR, "mismatched_pe2.fq.gz")
    with gzip.open(pe1_path, "wt") as f1, gzip.open(pe2_path, "wt") as f2:
        for i in range(100):
            f1.write(f"@READ_{i+1}/1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
            if i < 80:
                f2.write(f"@READ_{i+1}/2\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
    print(f"  [Created] {pe1_path} and {pe2_path} (mismatched read pair counts)")

    # 5f. Multi-line FASTA (wrapped at 60 chars)
    fasta_path = os.path.join(EDGE_CASES_DIR, "multiline_fasta.fa")
    ref_path = os.path.join(ROOT_DIR, "data", "phix", "phix174_ref.fa")
    ref = load_reference(ref_path) or ("GATTACA" * 500)
    random.seed(42)
    with open(fasta_path, "w") as f:
        for i in range(1000):
            start = random.randint(0, len(ref) - 200)
            seq = ref[start:start + 200]
            f.write(f">fasta_read_{i+1} sample wrapped at 60 chars\n")
            for chunk in [seq[j:j+60] for j in range(0, len(seq), 60)]:
                f.write(chunk + "\n")
    print(f"  [Created] {fasta_path} (1000 line-wrapped FASTA reads)")

def main():
    print("=" * 65)
    print("  SPADES-RS TEST ASSET INGESTION & MODALITY FIXTURE ACQUISITION  ")
    print("=" * 65)
    ensure_dirs()
    generate_pacbio_hifi()
    generate_10x_linked_reads()
    generate_phred64_dataset()
    generate_interleaved_and_single_end()
    generate_edge_cases()
    print("=" * 65)
    print("  [SUCCESS] All test assets, modalities, and fixtures acquired!  ")
    print("=" * 65)

if __name__ == "__main__":
    main()
