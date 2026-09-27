# spades-rs evaluation and testing tools

This directory contains validation harnesses, test asset generators, and biological evaluation scripts for `spades-rs`.

## Test harnesses and benchmark scripts

* `run_realworld_validation.py`: Automated test runner executing 19 operational scenarios across sequencing modalities, library formats, biological modes, and fault injections. Supports `--quick` for rapid regression testing (~20 seconds) and full execution for deep benchmarks.
* `fetch_test_assets.py`: Generates and acquires verified test datasets under `data/modalities/` and `data/edge_cases/` (PacBio HiFi, 10x Genomics barcoded reads, legacy Phred+64, interleaved paired-end, single-end, poly-N tails, and line-wrapped FASTA).
* `eval_assembly.py`: Standalone assembly evaluation script calculating contig counts, N50, L50, GC content, genome fraction, and reference coverage against ground-truth genomes without requiring external dependencies.
* `run_benchmark.py`: Measures execution time, CPU utilization, and peak resident memory (RSS) using `/usr/bin/time -v`.

## Biological evaluation scripts

### Bacterial and fungal isolates
* `audit_rrna_synteny.py`: Verifies bridging of all 7 identical ribosomal RNA operons (*rrnA* through *rrnH*) in *E. coli* K-12 MG1655.
* `audit_yeast_chromosomes.py`: Aligns contigs against the 16 chromosomes of *Saccharomyces cerevisiae* S288C, checking centromeric synteny (*CEN1*–*CEN16*) and repeat boundaries.
* `audit_tb_genes.py`: Evaluates clinical drug-resistance loci (*rpoB*, *inhA*, *gyrA*, *gyrB*, *pncA*, *embB*, *folC*, *gidB*, *rpsL*) and 155 *PE/PPE* repeat families in *Mycobacterium tuberculosis* H37Rv (65.6% GC).
* `audit_pa_genes.py`: Checks multidrug efflux operons (*mexAB-oprM*, *mexCD-oprJ*, *mexEF-oprN*), alginate regulatory genes, and pyoverdine loci in *Pseudomonas aeruginosa* PAO1 (66.6% GC).
* `audit_pf_genes.py`: Measures assembly across the 14 AT-rich chromosomes (19.4% GC) in *Plasmodium falciparum* 3D7.
* `audit_pombe_loci.py`, `find_pombe_genes.py`, `pombe_chr_stats.py`, `verify_pombe_identity.py`: Evaluates sequence recovery across the three chromosomes of *Schizosaccharomyces pombe* 972h-, verifying core cell-cycle regulators (*cdc2*, *cdc25*, *wee1*).

### Metagenomics
* `audit_zymo_community.py`: Computes per-species genome fraction across the 10-species ZymoBIOMICS mock community (8 bacteria, 2 yeasts).
* `check_zymo_chimeras.py`: Checks for inter-species chimeric contigs across co-assembled genomes.
* `zymo_per_species_quast.py`: Computes per-species assembly metrics against individual reference assemblies.
* `prepare_zymo_ref.py`: Validates reference genomes for the ZymoBIOMICS community standard.

### Single-cell and plasmid pipelines
* `audit_single_cell_ecoli.py`: Measures coverage fluctuations and dropout patterns in single-cell multiple displacement amplification (MDA) datasets.
* `audit_plasmid_kpn.py`: Assesses recovery and circularity of drug-resistance plasmids in *Klebsiella pneumoniae* BAA-2146.
* `audit_transcriptome_yeast.py`: Compares assembled transcripts against curated cDNA isoforms in *S. cerevisiae* RNA-Seq datasets.
