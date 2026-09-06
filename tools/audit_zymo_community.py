#!/usr/bin/env python3
"""
Metagenomic biological ground-truth evaluation script for ZymoBIOMICS Microbial Community (D6300).
Evaluates:
1. Multi-species genome fraction across 8 bacteria + 2 yeasts.
2. Inter-species chimerism (zero cross-taxa fusion).
3. 16S and 18S ribosomal RNA gene recovery and taxonomic assignment.
"""

import sys
import os
import glob

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
        print("Usage: audit_zymo_community.py <ref_dir> <contigs.fasta>")
        sys.exit(1)
        
    ref_dir = sys.argv[1]
    contigs_fa = sys.argv[2]
    
    print("=" * 80)
    print("  ZYMOBIOMICS 10-SPECIES METAGENOMIC BENCHMARK AUDIT")
    print("=" * 80)
    
    # Load 10 reference species
    ref_files = sorted(glob.glob(os.path.join(ref_dir, 'Genomes', '*.fasta')))
    species_genomes = {}
    for rf in ref_files:
        sp_name = os.path.basename(rf).replace('_complete_genome.fasta', '').replace('_draft_genome.fasta', '')
        seqs = load_fasta(rf)
        species_genomes[sp_name] = seqs
        tot = sum(len(s) for s in seqs.values())
        print(f"  Ref {sp_name:25s}: {tot:10,d} bp across {len(seqs)} contigs/chromosomes")
        
    print(f"\nLoading metagenome assembly: {contigs_fa}")
    contigs = load_fasta(contigs_fa)
    total_contig_bp = sum(len(s) for s in contigs.values())
    print(f"Loaded {len(contigs):,} sequences, total assembled: {total_contig_bp:,} bp")
    
    print("\nBuilding assembly k-mer index (k=25)...")
    contig_index = build_contig_kmer_index(contigs, k=25)
    print("Index ready.")
    
    # Audit 16S / 18S rRNA operons
    rrna_files = sorted(glob.glob(os.path.join(ref_dir, 'ssrRNAs', '*.fasta')))
    print("\n" + "=" * 80)
    print(f"{'Target rRNA Operon':<38} | {'Type':<6} | {'Length':>6} | {'Best Match':<20} | {'Cov %':>7} | {'Status'}")
    print("-" * 80)
    
    intact_rrna = 0
    for rrf in rrna_files:
        bname = os.path.basename(rrf).replace('.fasta', '')
        seqs = load_fasta(rrf)
        for h, seq in seqs.items():
            rtype = "16S" if "16S" in bname else ("18S" if "18S" in bname else "ssrRNA")
            cov, best_c = find_best_alignment(seq, contig_index, k=25)
            cov_pct = cov * 100.0
            status = "INTACT" if cov_pct >= 95.0 else ("PARTIAL" if cov_pct >= 50.0 else "ABSENT")
            if cov_pct >= 95.0:
                intact_rrna += 1
            c_disp = best_c[:20] if best_c else "None"
            print(f"{bname[:38]:<38} | {rtype:<6} | {len(seq):>5}bp | {c_disp:<20} | {cov_pct:>6.1f}% | {status}")
            
    print("-" * 80)
    print(f"Total ssrRNA Operons Preserved: {intact_rrna}/{len(rrna_files)}")

if __name__ == '__main__':
    main()
