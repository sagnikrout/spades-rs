import csv
import glob
import os

species_list = [
    "Bacillus_subtilis", "Cryptococcus_neoformans", "Enterococcus_faecalis",
    "Escherichia_coli", "Lactobacillus_fermentum", "Listeria_monocytogenes",
    "Pseudomonas_aeruginosa", "Saccharomyces_cerevisiae", "Salmonella_enterica",
    "Staphylococcus_aureus"
]

info_file = 'output/bench_zymo/quast_report/contigs_reports/contigs_report_contigs.mis_contigs.info'

inter_species_chimeras = []
intra_species_mis = []

if os.path.exists(info_file):
    current_contig = None
    with open(info_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('contig_'):
                current_contig = line.split()[0]
            elif 'Extensive misassembly' in line or 'translocation' in line or 'relocation' in line or 'inversion' in line:
                # check alignments for current_contig in all_alignments
                pass

# Let's inspect all_alignments_contigs.tsv to find contigs matching >1 species
tsv_file = 'output/bench_zymo/quast_report/contigs_reports/all_alignments_contigs.tsv'
contig_species = {}

with open(tsv_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 6 and parts[0].isdigit():
            c_name = parts[5]
            ref_name = parts[4]
            sp = None
            for s in species_list:
                if ref_name.startswith(s):
                    sp = s
                    break
            if sp:
                if c_name not in contig_species:
                    contig_species[c_name] = set()
                contig_species[c_name].add(sp)

multi_species_contigs = {c: sps for c, sps in contig_species.items() if len(sps) > 1}
total_aligned_contigs = len(contig_species)

print("=" * 70)
print("  ZYMOBIOMICS INTER-SPECIES CHIMERISM AUDIT")
print("=" * 70)
print(f"Total aligned contigs evaluated: {total_aligned_contigs:,}")
print(f"Contigs mapping to >1 species (Chimeras): {len(multi_species_contigs)}")
if multi_species_contigs:
    print("Examples of multi-species alignments:")
    for c, sps in list(multi_species_contigs.items())[:10]:
        print(f"  {c}: spans {', '.join(sps)}")
else:
    print("ZERO INTER-SPECIES CHIMERAS: All contigs cleanly map to a single species!")

purity = (1.0 - (len(multi_species_contigs) / max(1, total_aligned_contigs))) * 100.0
print(f"Species Taxonomic Separation Purity: {purity:.2f}%")
