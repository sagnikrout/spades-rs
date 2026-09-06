#!/usr/bin/env python3
"""
Biological ground-truth evaluation script for Pseudomonas aeruginosa PAO1.
Audits key virulence clusters, efflux pumps, alginate operon, quorum sensing,
and ribosomal RNA operons against assembled contigs.
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

def find_best_alignment(gene_seq, contigs, k=31):
    gene_len = len(gene_seq)
    if gene_len < k:
        k = max(11, gene_len // 2)
    
    # Build k-mer index of gene
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
        c_matches = 0
        rc_seq = revcomp(c_seq)
        
        # Count shared k-mers forward and reverse
        matched_positions = set()
        for strand_seq in [c_seq, rc_seq]:
            for i in range(0, len(strand_seq) - k + 1): # step 3 for speed
                km = strand_seq[i:i+k]
                if km in gene_kmers:
                    for pos in gene_kmers[km]:
                        matched_positions.add(pos)
        
        cov = len(matched_positions) / max(1, (gene_len - k + 1))
        if cov > best_matches:
            best_matches = cov
            best_contig = c_id
            
    return best_matches, best_contig

def parse_gff_target_genes(gff_path, target_names):
    genes = {}
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            feat_type = parts[2]
            if feat_type not in ('gene', 'CDS'):
                continue
            start, end, strand = int(parts[3]), int(parts[4]), parts[6]
            attrs = parts[8]
            attr_dict = {}
            for item in attrs.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    attr_dict[k] = v
            gene_name = attr_dict.get('Name', attr_dict.get('gene', ''))
            locus = attr_dict.get('locus_tag', '')
            product = attr_dict.get('product', '')
            
            for target in target_names:
                if (gene_name and target.lower() == gene_name.lower()) or (locus and target.lower() == locus.lower()):
                    if target not in genes or feat_type == 'CDS':
                        genes[target] = {
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'name': gene_name or locus,
                            'locus': locus,
                            'product': product
                        }
    return genes

def main():
    if len(sys.argv) < 4:
        print("Usage: audit_pa_genes.py <ref.fasta> <annotations.gff> <contigs.fasta>")
        sys.exit(1)
        
    ref_fa = sys.argv[1]
    gff_file = sys.argv[2]
    contigs_fa = sys.argv[3]
    
    print("=" * 65)
    print("  PAO1 HIGH-GC CLINICAL & VIRULENCE LOCI AUDIT")
    print("=" * 65)
    
    print(f"Loading reference: {ref_fa}")
    ref_seqs = load_fasta(ref_fa)
    ref_chrom = list(ref_seqs.values())[0]
    ref_len = len(ref_chrom)
    print(f"Reference genome length: {ref_len:,} bp")
    
    print(f"Loading assembly contigs: {contigs_fa}")
    contigs = load_fasta(contigs_fa)
    total_contig_bp = sum(len(s) for s in contigs.values())
    print(f"Loaded {len(contigs)} contigs, total: {total_contig_bp:,} bp")
    
    target_genes = [
        # Essential Housekeeping
        "dnaA", "gyrA", "gyrB", "parC", "parE", "rpoB", "recA",
        # Multidrug Efflux Pumps
        "mexA", "mexB", "oprM",
        "mexC", "mexD", "oprJ",
        "mexE", "mexF", "oprN",
        "mexX", "mexY",
        # Alginate Operon (Mucoidy / Exopolysaccharide)
        "algD", "alg8", "alg44", "algK", "algE", "algG", "algX", "algL", "algI", "algJ", "algF", "algA",
        # Pyoverdine Siderophore Cluster (Virulence)
        "pvdA", "pvdE", "pvdD", "pvdJ", "pvdI", "pvdL", "fpvA",
        # Quorum Sensing Systems
        "lasR", "lasI", "rhlR", "rhlI", "pqsA", "pqsR"
    ]
    
    annot_genes = parse_gff_target_genes(gff_file, target_genes)
    print(f"Found {len(annot_genes)} target loci in GFF annotations.")
    
    print("\n%-10s %-10s %-8s %-25s %-12s %s" % ("Gene", "Locus", "Length", "Category", "Recovery", "Status"))
    print("-" * 75)
    
    categories = {
        "Housekeeping": ["dnaA", "gyrA", "gyrB", "parC", "parE", "rpoB", "recA"],
        "Efflux Pumps": ["mexA", "mexB", "oprM", "mexC", "mexD", "oprJ", "mexE", "mexF", "oprN", "mexX", "mexY"],
        "Alginate Operon": ["algD", "alg8", "alg44", "algK", "algE", "algG", "algX", "algL", "algI", "algJ", "algF", "algA"],
        "Pyoverdine Siderophore": ["pvdA", "pvdE", "pvdD", "pvdJ", "pvdI", "pvdL", "fpvA"],
        "Quorum Sensing": ["lasR", "lasI", "rhlR", "rhlI", "pqsA", "pqsR"]
    }
    
    gene_cat = {}
    for cat, glist in categories.items():
        for g in glist:
            gene_cat[g] = cat
            
    intact_count = 0
    total_audited = 0
    
    for gname in target_genes:
        if gname not in annot_genes:
            continue
        info = annot_genes[gname]
        g_seq = ref_chrom[info['start']-1:info['end']]
        if info['strand'] == '-':
            g_seq = revcomp(g_seq)
            
        recov, c_id = find_best_alignment(g_seq, contigs, k=31)
        total_audited += 1
        
        status = "INTACT" if recov >= 0.99 else ("FRAGMENTED" if recov >= 0.50 else "ABSENT")
        if recov >= 0.99:
            intact_count += 1
            
        cat = gene_cat.get(gname, "Other")
        print("%-10s %-10s %-8d %-25s %6.1f%%     %s" % (
            gname, info['locus'], len(g_seq), cat, recov * 100, status
        ))
        
    print("-" * 75)
    print(f"OVERALL SUMMARY: {intact_count} / {total_audited} ({intact_count / max(1, total_audited) * 100:.1f}%) loci 100% INTACT in assembly.")
    print("=" * 75)

if __name__ == '__main__':
    main()
