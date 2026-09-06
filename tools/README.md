# Biological verification and benchmark audit tools

This directory contains standalone evaluation scripts used to audit `spades-rs` assemblies against published reference genomes and NCBI ground truths.

## Benchmark runners and metrics

* `run_benchmark.py`: Runs wall-clock, CPU usage, and resident set size (RSS) memory profiling comparisons between `spades-rs` and legacy SPAdes.
* `eval_assembly.py`: Runs QUAST and Minimap2 to calculate genome fraction, N50, mismatch rate, indel rate, and structural misassemblies.

## Biological audit scripts

### Bacterial and fungal isolates
* `audit_rrna_synteny.py`: Checks complete bridging of all 7 identical ribosomal RNA operons (*rrnA* through *rrnH*) in *E. coli* K-12 MG1655 without collapse.
* `audit_yeast_chromosomes.py`: Aligns contigs against the 16 chromosomes of *Saccharomyces cerevisiae* S288C, checking centromeric synteny (*CEN1*–*CEN16*) and retrotransposon boundaries.
* `audit_tb_genes.py`: Evaluates clinical drug-resistance loci (*rpoB*, *inhA*, *gyrA*, *gyrB*, *pncA*, *embB*, *folC*, *gidB*, *rpsL*) and 155 repetitive *PE/PPE* multigene families in *Mycobacterium tuberculosis* H37Rv (65.6% GC).
* `audit_pa_genes.py`: Verifies multidrug efflux operons (*mexAB-oprM*, *mexCD-oprJ*, *mexEF-oprN*), alginate regulatory clusters, and pyoverdine siderophore loci in *Pseudomonas aeruginosa* PAO1 (66.6% GC).
* `audit_pf_genes.py`: Evaluates assembly across the 14 AT-rich chromosomes (19.4% GC) and repetitive *var* / *rifin* / *stevor* multigene families in *Plasmodium falciparum* 3D7.
* `audit_pombe_loci.py`, `find_pombe_genes.py`, `pombe_chr_stats.py`, `verify_pombe_identity.py`: Measures recovery across the three chromosomes (3.5 Mb to 5.7 Mb) of *Schizosaccharomyces pombe* 972h-, verifying core cell-cycle regulators (*cdc2*, *cdc25*, *wee1*).

### Metagenomics
* `audit_zymo_community.py`: Computes per-species genome fraction across the 10-species ZymoBIOMICS mock community (8 bacteria, 2 yeasts).
* `check_zymo_chimeras.py`: Detects inter-species chimeric contigs across divergent bacterial genomes.
* `zymo_per_species_quast.py`: Runs automated per-species QUAST metrics against individual reference assemblies.
* `prepare_zymo_ref.py`: Downloads and validates the reference genomes for the ZymoBIOMICS standard.

### Single-cell and plasmid pipelines
* `audit_single_cell_ecoli.py`: Analyzes coverage depth fluctuations, dropouts, and structural misassemblies in single-cell multiple displacement amplification (MDA) datasets.
* `audit_plasmid_kpn.py`: Assesses extraction of circular plasmids (*pKPN-182*, *pKPN-116*, *pKPN-66*, *pKPN-4*) in *Klebsiella pneumoniae* BAA-2146.
* `audit_transcriptome_yeast.py`: Compares assembled transcripts against SGD curated isoforms in *S. cerevisiae* RNA-Seq, checking transcript completeness, chimera rates, and expression correlation.
