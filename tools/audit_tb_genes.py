#!/usr/bin/env python3
"""
Audit Mycobacterium tuberculosis H37Rv Assembly for:
1. Core AMR diagnostic loci (rpoB, katG, inhA, gyrA, gyrB, pncA, embB).
2. PE/PPE multigene repeat family integrity.
"""

import sys
import os
import re

AMR_GENES = {
    "Rv0667": ("rpoB", "Rifampicin resistance (RRDR)"),
    "Rv1908c": ("katG", "Isoniazid high-level resistance"),
    "Rv1484": ("inhA", "Isoniazid / Ethionamide resistance"),
    "Rv0006": ("gyrA", "Fluoroquinolone resistance (QRDR)"),
    "Rv0005": ("gyrB", "Fluoroquinolone resistance"),
    "Rv2043c": ("pncA", "Pyrazinamide resistance"),
    "Rv3795": ("embB", "Ethambutol resistance"),
    "Rv2447c": ("folC", "p-Aminosalicylic acid resistance"),
    "Rv3919c": ("gidB", "Streptomycin resistance"),
    "Rv0682": ("rpsL", "Streptomycin high-level resistance"),
}

def parse_gff(gff_path):
    genes = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            if fields[2] == "gene":
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                info = fields[8]
                locus_match = re.search(r"locus_tag=([^;]+)", info)
                gene_match = re.search(r"gene=([^;]+)", info)
                
                if locus_match:
                    locus = locus_match.group(1)
                    name = gene_match.group(1) if gene_match else locus
                    is_pe_ppe = "PE" in name or "PPE" in name or "PE_PGRS" in info or "PPE family" in info
                    genes[locus] = {
                        "name": name,
                        "start": start,
                        "end": end,
                        "length": end - start + 1,
                        "strand": strand,
                        "is_pe_ppe": is_pe_ppe,
                    }
    return genes

def parse_alignments(align_tsv):
    alignments = []
    if not os.path.exists(align_tsv):
        return alignments
    with open(align_tsv) as f:
        for line in f:
            if line.startswith("CONTIG") or line.startswith("S1"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 7 and parts[0].isdigit():
                r1 = int(parts[0])
                r2 = int(parts[1])
                alignments.append({
                    "r_start": min(r1, r2),
                    "r_end": max(r1, r2),
                    "q_start": int(parts[2]),
                    "q_end": int(parts[3]),
                    "contig": parts[5],
                    "idy": float(parts[6]),
                })
    return alignments

def main():
    gff_path = "data/tb/h37rv_annotations.gff"
    align_tsv = sys.argv[1] if len(sys.argv) > 1 else "output/quast_tb_h37rv/contigs_reports/all_alignments_scaffolds.tsv"
    
    if not os.path.exists(gff_path):
        print(f"Error: Missing {gff_path}")
        sys.exit(1)
        
    genes = parse_gff(gff_path)
    alignments = parse_alignments(align_tsv)
    
    print("================================================================================")
    print("      MYCOBACTERIUM TUBERCULOSIS H37RV: AMR & PE/PPE GENE INTEGRITY AUDIT       ")
    print("================================================================================")
    
    # 1. Audit AMR Diagnostic Genes
    print("\n[SECTION 1: CLINICAL ANTIMICROBIAL RESISTANCE (AMR) DIAGNOSTIC LOCI]")
    print(f"{'Gene':<8} {'Locus':<10} {'Position (bp)':<22} {'Clinical Target':<38} {'Status':<10}")
    print("-" * 88)
    
    amr_intact = 0
    for locus, (gene_name, target_desc) in AMR_GENES.items():
        if locus in genes:
            g = genes[locus]
            g_start, g_end = g["start"], g["end"]
            
            # Check if any alignment completely spans the gene with >= 50 bp flanks
            spanning = [
                a for a in alignments
                if a["r_start"] <= g_start - 30 and a["r_end"] >= g_end + 30
            ]
            
            if spanning:
                status = "INTACT"
                amr_intact += 1
            else:
                # Check coverage
                cov = max((min(a["r_end"], g_end) - max(a["r_start"], g_start) + 1 for a in alignments if a["r_end"] >= g_start and a["r_start"] <= g_end), default=0)
                frac = cov / g["length"] * 100
                status = f"{frac:.1f}%"
                
            pos_str = f"{g_start:,} - {g_end:,}"
            print(f"{gene_name:<8} {locus:<10} {pos_str:<22} {target_desc:<38} {status:<10}")
            
    print("-" * 88)
    print(f"AMR Diagnostic Loci Summary: {amr_intact} / {len(AMR_GENES)} genes 100% intact with flanking context")
    
    # 2. Audit PE/PPE Multigene Family
    pe_ppe_genes = {k: v for k, v in genes.items() if v["is_pe_ppe"]}
    print(f"\n[SECTION 2: REPETITIVE PE / PPE MULTIGENE FAMILY RECOVERY (Total: {len(pe_ppe_genes)} genes)]")
    
    pe_intact_99 = 0
    pe_intact_95 = 0
    pe_partial = 0
    
    for locus, g in pe_ppe_genes.items():
        g_start, g_end = g["start"], g["end"]
        
        # Calculate maximum continuous alignment coverage of the gene
        max_span = 0
        for a in alignments:
            if a["r_end"] >= g_start and a["r_start"] <= g_end:
                overlap = min(a["r_end"], g_end) - max(a["r_start"], g_start) + 1
                if overlap > max_span:
                    max_span = overlap
                    
        frac = max_span / g["length"] * 100
        if frac >= 99.0:
            pe_intact_99 += 1
            pe_intact_95 += 1
        elif frac >= 95.0:
            pe_intact_95 += 1
        elif frac >= 50.0:
            pe_partial += 1
            
    print(f"  >= 99% Complete (Zero-break full length): {pe_intact_99} / {len(pe_ppe_genes)} ({pe_intact_99/len(pe_ppe_genes)*100:.1f}%)")
    print(f"  >= 95% Complete (Coding domain intact):   {pe_intact_95} / {len(pe_ppe_genes)} ({pe_intact_95/len(pe_ppe_genes)*100:.1f}%)")
    print(f"  Partial (>50%):                           {pe_partial} / {len(pe_ppe_genes)} ({pe_partial/len(pe_ppe_genes)*100:.1f}%)")
    print("================================================================================\n")

if __name__ == "__main__":
    main()
