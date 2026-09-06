#!/usr/bin/env python3
"""
Biological ground-truth evaluation script for Schizosaccharomyces pombe 972h-.
Audits:
1. Master cell-cycle regulators (cdc2, cdc13, cdc25, wee1, rad3, chk1, pol1, tor1, act1, ade6, ura4, leu1)
2. Mating-type cassette (mat1, mat3) on Chromosome II
3. Regional centromeres (cen1, cen2, cen3) on Chromosomes I, II, III
against de novo assembled contigs/scaffolds.
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
        print("Usage: audit_pombe_loci.py <ref.fasta> <contigs.fasta>")
        sys.exit(1)
        
    ref_fa = sys.argv[1]
    contigs_fa = sys.argv[2]
    
    print("=" * 80)
    print("  SCHIZOSACCHAROMYCES POMBE 972h- BIOLOGICAL & GENOMIC AUDIT")
    print("=" * 80)
    
    print(f"Loading reference: {ref_fa}")
    ref_seqs = load_fasta(ref_fa)
    for h, s in ref_seqs.items():
        print(f"  {h:15s}: {len(s):,} bp")
    total_ref_bp = sum(len(s) for s in ref_seqs.values())
    print(f"Total reference genome: {total_ref_bp:,} bp across {len(ref_seqs)} chromosomes")
    
    print(f"\nLoading assembly: {contigs_fa}")
    contigs = load_fasta(contigs_fa)
    total_contig_bp = sum(len(s) for s in contigs.values())
    print(f"Loaded {len(contigs)} sequences, total: {total_contig_bp:,} bp")
    
    print("Building k-mer search index (k=25)...")
    contig_index = build_contig_kmer_index(contigs, k=25)
    print("Index ready.")
    
    # Map reference accession IDs
    chr1_id = [k for k in ref_seqs if 'NC_003424' in k][0]
    chr2_id = [k for k in ref_seqs if 'NC_003423' in k][0]
    chr3_id = [k for k in ref_seqs if 'NC_003421' in k][0]
    
    targets = {
        # Master Cell Cycle Regulators
        "cdc2 (CDK1 Master Kinase)": (chr2_id, 1500208, 1502095, "SPBC11B10.09"),
        "cdc13 (B-type Cyclin)": (chr2_id, 420111, 423276, "SPBC582.03"),
        "cdc25 (Mitotic Phosphatase)": (chr1_id, 479228, 482373, "SPAC24H6.05"),
        "wee1 (Mitotic Inhibitor Kinase)": (chr3_id, 721804, 725835, "SPCC18B5.03"),
        "rad3 (ATR Checkpoint Master)": (chr2_id, 904184, 913629, "SPBC216.05"),
        "chk1 (DNA Damage Effector Kinase)": (chr3_id, 1059796, 1064443, "SPCC1259.13"),
        "pol1 (DNA Pol Alpha Subunit)": (chr1_id, 3430955, 3436241, "SPAC3H5.06C"),
        "tor1 (TOR Complex Kinase)": (chr2_id, 3075748, 3083248, "SPBC30D10.10C"),
        # Classical Genetic Markers
        "act1 (Actin)": (chr2_id, 1475485, 1477333, "SPBC32H8.12C"),
        "ura4 (OMP Decarboxylase)": (chr3_id, 115589, 116625, "SPCC330.05C"),
        "ade6 (Purine Biosynthesis)": (chr3_id, 1316281, 1318035, "SPCC1322.13"),
        "leu1 (Leucine Biosynthesis)": (chr2_id, 1974390, 1975674, "SPBC1A4.02C"),
        # Mating Type Locus & Cassettes
        "mat1-Mi/Mc (Expressed Mating Cassette)": (chr2_id, 2114219, 2115027, "SPBC23G7.17C/09"),
        "mat3-Mi/Mc (Silenced Mating Cassette)": (chr2_id, 2132586, 2133394, "SPBC1711.01C/02"),
    }
    
    print("\n" + "=" * 80)
    print(f"{'Target Locus':<38} | {'Ref Chrom':<12} | {'Length':>6} | {'Contig Match':<32} | {'Cov %':>7} | {'Status'}")
    print("-" * 80)
    
    intact_count = 0
    total_targets = len(targets)
    
    for name, (chrom, s, e, locus_id) in targets.items():
        gene_len = e - s + 1
        gene_seq = ref_seqs[chrom][s-1:e]
        cov, best_c = find_best_alignment(gene_seq, contig_index, k=25)
        cov_pct = cov * 100.0
        
        if cov_pct >= 95.0:
            status = "INTACT"
            intact_count += 1
        elif cov_pct >= 50.0:
            status = "FRAGMENTED"
        else:
            status = "MISSING"
            
        c_disp = best_c[:32] if best_c else "None"
        print(f"{name:<38} | {chrom:<12} | {gene_len:>5}bp | {c_disp:<32} | {cov_pct:>6.1f}% | {status}")
        
    print("-" * 80)
    print(f"Summary: {intact_count}/{total_targets} loci reconstructed intact (>=95% identity across full open reading frame).")
    
    print("\n" + "=" * 80)
    print("  REGIONAL CENTROMERIC HETEROCHROMATIN AUDIT")
    print("=" * 80)
    
    centromeres = {
        "cen1 (Chr I Regional Centromere ~35 kb)": (chr1_id, 3752000, 3788000),
        "cen2 (Chr II Regional Centromere ~65 kb)": (chr2_id, 1602000, 1667000),
        "cen3 (Chr III Regional Centromere ~110 kb)": (chr3_id, 1070000, 1142000),
    }
    
    for cen_name, (chrom, s, e) in centromeres.items():
        cen_len = e - s + 1
        cen_seq = ref_seqs[chrom][s-1:e]
        cov, best_c = find_best_alignment(cen_seq, contig_index, k=25)
        cov_pct = cov * 100.0
        c_disp = best_c[:32] if best_c else "None"
        print(f"  {cen_name}: {cen_len:,} bp -> Best contig {c_disp}, span coverage: {cov_pct:.1f}%")

if __name__ == '__main__':
    main()
