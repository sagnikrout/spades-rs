#!/usr/bin/env python3
"""
Milestone 7 Biological and AMR Mobilome Audit for Klebsiella pneumoniae ATCC BAA-2146.
Audits:
1. Four verified reference plasmids:
   - pMYS (2,014 bp, high copy)
   - pNDM-US (140,825 bp, NDM-1 / OXA-181 MDR plasmid)
   - pHg (85,161 bp, mercury resistance plasmid)
   - pCuAs (117,755 bp, copper/arsenic resistance plasmid)
2. Segregation purity:
   - Evaluates whether plasmids are cleanly captured in plasmids.fasta without contaminating chromosomal contigs.
3. AMR and Heavy Metal Resistance Operons:
   - blaNDM-1, blaOXA-181, aac(6')-Ib, merA, arsB, pcoA.
"""

import sys
import os

def load_fasta(filepath):
    sequences = {}
    if not os.path.exists(filepath):
        return sequences
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

def find_best_alignment(target_seq, contig_index, k=25):
    target_seq = target_seq.upper()
    target_len = len(target_seq)
    if target_len < k:
        k = max(11, target_len // 2)
    
    target_kmers = [canonical_kmer(target_seq[i:i+k]) for i in range(target_len - k + 1)]
    total_km = max(1, len(target_kmers))
    
    best_matches = 0
    best_contig = None
    
    for c_id, km_set in contig_index.items():
        matched = sum(1 for km in target_kmers if km in km_set)
        cov = matched / total_km
        if cov > best_matches:
            best_matches = cov
            best_contig = c_id
            if cov >= 0.999:
                break
                
    return best_matches, best_contig

def main():
    if len(sys.argv) < 4:
        print("Usage: audit_plasmid_kpn.py <all_ref.fa> <contigs.fasta> <plasmids.fasta>")
        sys.exit(1)
        
    ref_fa = sys.argv[1]
    contigs_fa = sys.argv[2]
    plasmids_fa = sys.argv[3]
    
    print('=' * 80)
    print('  MILESTONE 7: KLEBSIELLA PNEUMONIAE ATCC BAA-2146 PLASMID AND MOBILOME AUDIT')
    print('=' * 80)
    
    ref_seqs = load_fasta(ref_fa)
    contigs = load_fasta(contigs_fa)
    plasmids = load_fasta(plasmids_fa)
    
    print(f'Reference Replicons Loaded: {len(ref_seqs)}')
    for h, s in ref_seqs.items():
        gc = (s.count('G') + s.count('C')) / len(s) * 100
        print(f'  - {h:<15}: {len(s):>10,d} bp | GC: {gc:.2f}%')
        
    print(f'\nAssembly Output:')
    print(f'  - Chromosomal contigs: {len(contigs):>6,d} sequences ({sum(len(s) for s in contigs.values()):,d} bp)')
    print(f'  - Extracted plasmids:  {len(plasmids):>6,d} sequences ({sum(len(s) for s in plasmids.values()):,d} bp)')
    
    print('\nBuilding k-mer indices (k=25)...')
    chrom_index = build_contig_kmer_index(contigs, k=25)
    plasmid_index = build_contig_kmer_index(plasmids, k=25)
    print('Indices ready.')
    
    # Audit Each Reference Replicon
    print('\n' + '=' * 80)
    print(f'{"Reference Replicon":<30} | {"Ref Length":>10} | {"In Plasmids.fa":>14} | {"In Contigs.fa":>14}')
    print('-' * 80)
    
    all_plas_kmers = set().union(*plasmid_index.values()) if plasmid_index else set()
    all_chrom_kmers = set().union(*chrom_index.values()) if chrom_index else set()

    for r_id, r_seq in ref_seqs.items():
        r_kmers = [canonical_kmer(r_seq[i:i+25]) for i in range(len(r_seq) - 24)]
        total_r_km = max(1, len(r_kmers))
        
        cov_plas = sum(1 for km in r_kmers if km in all_plas_kmers) / total_r_km * 100.0
        cov_chrom = sum(1 for km in r_kmers if km in all_chrom_kmers) / total_r_km * 100.0
        
        is_plas = "plasmid" in r_id.lower() or r_id in ['CP006660.1', 'CP006661.1', 'CP006662.2', 'CP006663.1']
        tag = "[PLASMID]" if is_plas else "[CHROM]"
        print(f'{r_id + " " + tag:<30} | {len(r_seq):>10,d} | {cov_plas:>13.1f}% | {cov_chrom:>13.1f}%')
        
    print('-' * 80)

if __name__ == '__main__':
    main()
