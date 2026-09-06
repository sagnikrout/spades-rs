#!/usr/bin/env python3
"""
Biological ground-truth evaluation script for Plasmodium falciparum 3D7.
Audits antimalarial drug resistance genes (pfcrt, pfmdr1, k13, dhfr, dhps),
essential invasion antigens (csp, msp1, ama1, eba175), and housekeeping loci
across 14 nuclear chromosomes against assembled contigs/scaffolds.
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
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences[current_header] = ''.join(current_seq)
    return sequences

def revcomp(seq):
    tr = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(tr)[::-1]

def find_best_alignment(gene_seq, contigs, k=25):
    gene_len = len(gene_seq)
    if gene_len < k:
        k = max(11, gene_len // 2)
    
    gene_kmers = {}
    for i in range(gene_len - k + 1):
        km = gene_seq[i:i+k]
        if km not in gene_kmers:
            gene_kmers[km] = []
        gene_kmers[km].append(i)
        
    best_matches = 0
    best_contig = None
    
    for c_id, c_seq in contigs.items():
        if len(c_seq) < k:
            continue
        rc_seq = revcomp(c_seq)
        matched_positions = set()
        for strand_seq in [c_seq, rc_seq]:
            for i in range(len(strand_seq) - k + 1):
                km = strand_seq[i:i+k]
                if km in gene_kmers:
                    for pos in gene_kmers[km]:
                        matched_positions.add(pos)
        
        cov = len(matched_positions) / max(1, (gene_len - k + 1))
        if cov > best_matches:
            best_matches = cov
            best_contig = c_id
            
    return best_matches, best_contig

def parse_gff_target_genes(gff_path, targets):
    genes = {}
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            feat_type = parts[2]
            if feat_type != 'gene':
                continue
            chrom, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
            attrs = parts[8]
            attr_dict = {}
            for item in attrs.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    attr_dict[k] = v
                    
            locus = attr_dict.get('locus_tag', attr_dict.get('Name', ''))
            product = attr_dict.get('product', '')
            
            for t_name, t_locus in targets.items():
                if t_locus == locus or t_name.lower() in attrs.lower():
                    genes[t_name] = {
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'locus': locus,
                        'product': product
                    }
    return genes

def main():
    if len(sys.argv) < 4:
        print("Usage: audit_pf_genes.py <ref.fasta> <annotations.gff> <contigs.fasta>")
        sys.exit(1)
        
    ref_fa = sys.argv[1]
    gff_file = sys.argv[2]
    contigs_fa = sys.argv[3]
    
    print("=" * 75)
    print("  PLASMODIUM FALCIPARUM 3D7 CLINICAL & ANTIMALARIAL AUDIT")
    print("=" * 75)
    
    print(f"Loading reference: {ref_fa}")
    ref_seqs = load_fasta(ref_fa)
    total_ref_bp = sum(len(s) for s in ref_seqs.values())
    print(f"Reference genome length: {total_ref_bp:,} bp across {len(ref_seqs)} chromosomes")
    
    print(f"Loading assembly: {contigs_fa}")
    contigs = load_fasta(contigs_fa)
    total_contig_bp = sum(len(s) for s in contigs.values())
    print(f"Loaded {len(contigs)} contigs, total: {total_contig_bp:,} bp")
    
    targets = {
        # Clinical Drug Resistance Targets
        "pfcrt (Chloroquine)": "PF3D7_0709000",
        "pfmdr1 (Multidrug)": "PF3D7_0523000",
        "k13 (Artemisinin)": "PF3D7_1343700",
        "pfdhfr (Pyrimethamine)": "PF3D7_0417200",
        "pfdhps (Sulfadoxine)": "PF3D7_0810800",
        # Key Vaccine & Surface Antigens
        "csp (Circumsporozoite)": "PF3D7_0304600",
        "msp1 (Merozoite Surface 1)": "PF3D7_0930300",
        "ama1 (Apical Membrane 1)": "PF3D7_1133400",
        "eba175 (Erythrocyte Binding)": "PF3D7_0731500",
        # Housekeeping / Replication Markers
        "gapdh": "PF3D7_1462800",
        "act1 (Actin-1)": "PF3D7_1246200",
        "rpoB (RNA Pol Beta)": "PF3D7_1017000",
        "dnaA (Replication)": "PF3D7_1342600"
    }
    
    annot_genes = parse_gff_target_genes(gff_file, targets)
    print(f"Found {len(annot_genes)} target loci in GFF annotations.")
    
    print("\n%-28s %-16s %-12s %-8s %-10s %s" % ("Gene / Clinical Role", "Locus Tag", "Chromosome", "Length", "Recovery", "Status"))
    print("-" * 85)
    
    intact_count = 0
    total_audited = 0
    
    for gname, locus_tag in targets.items():
        if gname not in annot_genes:
            continue
        info = annot_genes[gname]
        chrom = info['chrom']
        if chrom not in ref_seqs:
            continue
        g_seq = ref_seqs[chrom][info['start']-1:info['end']]
        if info['strand'] == '-':
            g_seq = revcomp(g_seq)
            
        recov, c_id = find_best_alignment(g_seq, contigs, k=25)
        total_audited += 1
        
        status = "INTACT" if recov >= 0.99 else ("FRAGMENTED" if recov >= 0.50 else "ABSENT")
        if recov >= 0.99:
            intact_count += 1
            
        print("%-28s %-16s %-12s %-8d %6.1f%%    %s" % (
            gname, info['locus'], chrom, len(g_seq), recov * 100, status
        ))
        
    print("-" * 85)
    print(f"OVERALL SUMMARY: {intact_count} / {total_audited} ({intact_count / max(1, total_audited) * 100:.1f}%) loci 100% INTACT in assembly.")
    print("=" * 85)

if __name__ == '__main__':
    main()
