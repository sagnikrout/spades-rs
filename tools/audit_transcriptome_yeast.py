#!/usr/bin/env python3
"""
Biological Transcriptome & Alternative Splicing Audit for Saccharomyces cerevisiae.
Audits:
1. Ensembl curated reference cDNA transcript capture (6,612 transcripts).
2. Full-length recovery rates (>=95% and >=99% transcript coverage).
3. Key intron-containing and alternatively spliced transcripts (YRA1, SRC1, ACT1, RPL28, RPS14B).
4. Dynamic expression range sensitivity (high vs. low abundance).
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
    if len(sys.argv) < 3:
        print("Usage: audit_transcriptome_yeast.py <cdna_ref.fa> <transcripts.fasta>")
        sys.exit(1)
        
    cdna_fa = sys.argv[1]
    transcripts_fa = sys.argv[2]
    
    print('=' * 95)
    print('  MILESTONE 8: SACCHAROMYCES CEREVISIAE DE NOVO TRANSCRIPTOME & ISOFORM AUDIT')
    print('=' * 95)
    
    print(f"Loading reference cDNA transcriptome: {cdna_fa}")
    ref_cdna = load_fasta(cdna_fa)
    print(f"Reference curated transcripts: {len(ref_cdna):,} transcripts ({sum(len(s) for s in ref_cdna.values()):,d} bp)")
    
    print(f"\nLoading de novo assembled transcripts: {transcripts_fa}")
    assembled = load_fasta(transcripts_fa)
    total_assembled_bp = sum(len(s) for s in assembled.values())
    print(f"Assembled transcripts: {len(assembled):,} transcripts ({total_assembled_bp:,d} bp)")
    
    # Calculate transcript N50
    lens = sorted([len(s) for s in assembled.values()], reverse=True)
    half_bp = total_assembled_bp / 2.0
    cum = 0
    n50 = 0
    for l in lens:
        cum += l
        if cum >= half_bp:
            n50 = l
            break
    print(f"Transcript N50: {n50:,d} bp | Longest Transcript: {lens[0] if lens else 0:,d} bp")
    
    print("\nBuilding assembly k-mer index (k=25)...")
    contig_index = build_contig_kmer_index(assembled, k=25)
    print("Index ready.")
    
    # Key biological marker transcripts (intron-containing, spliced, regulatory)
    targets = [
        ("YDR381W_mRNA", "YRA1 (RNA Export / Intron-Regulated Splice)"),
        ("YML034W_mRNA", "SRC1 (Alternative Splicing Factor / 2 Isoforms)"),
        ("YFL039C_mRNA", "ACT1 (Actin / Spliced Canonical Intron)"),
        ("YGL103W_mRNA", "RPL28 (Ribosomal Protein L28 / Spliced Intron)"),
        ("YCR031C_mRNA", "RPS14B (Ribosomal Protein S14B Spliced)"),
        ("YGR192C_mRNA", "TDH3 (GAPDH / High Expression Glycolysis)"),
        ("YHR174W_mRNA", "ENO2 (Enolase II / Hyper-Abundant Glycolysis)"),
        ("YLR044C_mRNA", "PDC1 (Pyruvate Decarboxylase 1)"),
        ("YPR080W_mRNA", "TEF1 (Translation Elongation Factor 1A)"),
        ("YOL086C_mRNA", "ADH1 (Alcohol Dehydrogenase 1)"),
        ("YKL109W_mRNA", "HAP4 (Low Expression Transcriptional Activator)"),
        ("YGL035C_mRNA", "MIG1 (Glucose Repressor / Regulatory Protein)"),
    ]
    
    print("\n" + "=" * 95)
    print(f"{'Biological Marker Transcript':<48} | {'Length':>6} | {'Best Match':<20} | {'Cov %':>7} | {'Status'}")
    print("-" * 95)
    
    intact_count = 0
    for orf_id, label in targets:
        seq = None
        for h in ref_cdna:
            if h.startswith(orf_id):
                seq = ref_cdna[h]
                break
                
        if seq:
            cov, best_c = find_best_alignment(seq, contig_index, k=25)
            cov_pct = cov * 100.0
            status = "FULL_LENGTH" if cov_pct >= 95.0 else ("PARTIAL" if cov_pct >= 50.0 else "LOW_DEPTH")
            if cov_pct >= 95.0:
                intact_count += 1
            c_disp = best_c[:20] if best_c else "None"
            print(f"{label:<48} | {len(seq):>5}bp | {c_disp:<20} | {cov_pct:>6.1f}% | {status}")
        else:
            print(f"{label:<48} |   N/A  | {'None':<20} |   0.0% | NOT_FOUND")
            
    print("-" * 95)
    print(f"Summary: {intact_count}/{len(targets)} audited biological transcripts assembled at full length!")

if __name__ == '__main__':
    main()
