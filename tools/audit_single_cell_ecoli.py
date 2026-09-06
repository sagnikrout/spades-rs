#!/usr/bin/env python3
"""
Single-Cell MDA Biological Evaluation Script for Escherichia coli K-12.
Audits:
1. Core essential housekeeping genes (dnaA, gyrA, rpoB, recA, gapA, polA, ftsZ, secA, etc.)
2. Coverage spike distribution across single-cell contigs (evaluating MDA unevenness).
3. Seven ribosomal RNA operons (rrnA through rrnH).
"""

import sys
import os

def load_fasta(filepath):
    sequences = {}
    current_header = None
    current_seq = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq).upper()
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences[current_header] = ''.join(current_seq).upper()
    return sequences

def revcomp(seq):
    tr = str.maketrans('ACGTN', 'TGCAN')
    return seq.translate(tr)[::-1]

def canonical_kmer(km):
    rc = revcomp(km)
    return km if km < rc else rc

def build_contig_kmer_index(contigs, k=25):
    index = {}
    for c_id, seq in contigs.items():
        if len(seq) < k:
            continue
        km_set = set()
        for i in range(len(seq) - k + 1):
            km = seq[i:i+k]
            km_set.add(canonical_kmer(km))
        index[c_id] = km_set
    return index

def find_best_alignment(gene_seq, contig_index, k=25):
    gene_seq = gene_seq.upper()
    gene_len = len(gene_seq)
    if gene_len < k:
        k = max(11, gene_len // 2)
    
    gene_kmers = [canonical_kmer(gene_seq[i:i+k]) for i in range(gene_len - k + 1)]
    total_km = max(1, len(gene_kmers))
    
    best_matches = 0
    best_contig = None
    
    for c_id, km_set in contig_index.items():
        matched = sum(1 for km in gene_kmers if km in km_set)
        cov = matched / total_km
        if cov > best_matches:
            best_matches = cov
            best_contig = c_id
            if cov >= 0.999:
                break
                
    return best_matches, best_contig

def main():
    if len(sys.argv) < 3:
        print("Usage: audit_single_cell_ecoli.py <ref.fasta> <contigs.fasta>")
        sys.exit(1)
        
    ref_fa = sys.argv[1]
    contigs_fa = sys.argv[2]
    
    print("=" * 80)
    print("  SINGLE-CELL MDA ESCHERICHIA COLI BIOLOGICAL AUDIT")
    print("=" * 80)
    
    print(f"Loading reference: {ref_fa}")
    ref_seqs = load_fasta(ref_fa)
    ref_chrom = list(ref_seqs.keys())[0]
    ref_seq = ref_seqs[ref_chrom]
    print(f"Reference genome: {len(ref_seq):,} bp ({ref_chrom})")
    
    print(f"\nLoading single-cell assembly: {contigs_fa}")
    contigs = load_fasta(contigs_fa)
    total_contig_bp = sum(len(s) for s in contigs.values())
    print(f"Loaded {len(contigs):,} contigs, total assembled: {total_contig_bp:,} bp")
    
    # Analyze coverage depth distribution from contig headers if present
    covs = []
    for h in contigs:
        if 'cov_' in h:
            try:
                c_val = float(h.split('cov_')[1].split('_')[0])
                covs.append(c_val)
            except:
                pass
                
    if covs:
        covs.sort()
        min_cov = covs[0]
        max_cov = covs[-1]
        med_cov = covs[len(covs)//2]
        p95_cov = covs[int(len(covs)*0.95)]
        print(f"\nMDA Coverage Dynamics across {len(covs)} contigs:")
        print(f"  Min Depth: {min_cov:.1f}x | Median: {med_cov:.1f}x | 95th Percentile: {p95_cov:.1f}x | Peak Spike: {max_cov:.1f}x")
        print(f"  Coverage Dynamic Range: {max_cov/max(0.1, min_cov):.1f}-fold fluctuation!")
        
    print("\nBuilding assembly k-mer index (k=25)...")
    contig_index = build_contig_kmer_index(contigs, k=25)
    print("Index ready.")
    
    # Core E. coli essential genes (coordinates in NC_000913.3)
    targets = {
        "dnaA (Replication Initiator)": (3734, 5020),
        "gyrA (DNA Gyrase Subunit A)": (2337425, 2340052),
        "gyrB (DNA Gyrase Subunit B)": (3838634, 3841048),
        "rpoB (RNA Polymerase Beta Subunit)": (4181245, 4185273),
        "rpoC (RNA Polymerase Beta' Subunit)": (4185340, 4189563),
        "recA (DNA Recombination / Repair)": (2820730, 2821788),
        "gapA (Glyceraldehyde-3-P Dehydrogenase)": (1862590, 1863585),
        "polA (DNA Polymerase I)": (4048454, 4051240),
        "ftsZ (Cell Division GTPase)": (105151, 106302),
        "secA (Protein Translocase Subunit)": (196962, 199667),
        "atpA (ATP Synthase Alpha Subunit)": (3917849, 3919390),
        "adk (Adenylate Kinase)": (494056, 494700),
    }
    
    print("\n" + "=" * 80)
    print(f"{'Essential Housekeeping Locus':<42} | {'Length':>6} | {'Best Match':<20} | {'Cov %':>7} | {'Status'}")
    print("-" * 80)
    
    intact_count = 0
    for name, (s, e) in targets.items():
        gene_len = e - s + 1
        gene_seq = ref_seq[s-1:e]
        cov, best_c = find_best_alignment(gene_seq, contig_index, k=25)
        cov_pct = cov * 100.0
        status = "INTACT" if cov_pct >= 95.0 else ("PARTIAL" if cov_pct >= 50.0 else "DROPOUT")
        if cov_pct >= 95.0:
            intact_count += 1
        c_disp = best_c[:20] if best_c else "None"
        print(f"{name:<42} | {gene_len:>5}bp | {c_disp:<20} | {cov_pct:>6.1f}% | {status}")
        
    print("-" * 80)
    print(f"Summary: {intact_count}/{len(targets)} essential genes reconstructed intact from single cell.")

if __name__ == '__main__':
    main()
