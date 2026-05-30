# ChickEmbryo_ATAC_RNA

## Project description
ATAC-seq and RNA-seq integration pipeline for chicken embryo
(Gallus gallus, GalGal6/GRCg6a) antagomir study.
Comparisons: AM133, AM1/206, AMAll vs AMScr scramble control.
Genome annotation: Ensembl v106.

## Authors
Shahrzad Moradi Fard
shahrzadmf94@gmail.com
University of East Anglia

## Reproducing this analysis

### Requirements
- R version 4.5.3
- Quarto version 1.9.38
- All R package dependencies managed via renv

### Setup (run once before anything else)
1. Clone this repository
2. Open ChickEmbryo_ATAC_RNA.Rproj in RStudio
3. Run in R console: renv::restore()
   This installs all required packages at exact versions used
4. Place input data files in data/ as described below
5. Run scripts in order

## Data
Raw ATAC-seq data: deposited at NCBI SRA (accession pending)
RNA-seq data: processed by Novomagic, filtered at
|log2FC| >= 0.322 AND padj <= 0.05

Raw data files are not included in this repository.
Place files as follows before running scripts:
- data/raw/atac/deseq2/      DiffBind DESeq2 output TSV files
- data/raw/atac/narrowpeaks/ ATAC-seq narrow peak files
- data/raw/atac/bigwig/      BigWig coverage files
- data/processed/rna/        RNA-seq DEG files from Novomagic

## Scripts

### Script 0: ATAC-seq Direction Correction
File: scripts/Script0_ATAC_direction_correction.R
Documentation: quarto/Script0_ATAC_direction_correction.qmd
Purpose: Corrects fold change direction in raw DiffBind output.
DiffBind calculates Fold as log2(Group2/Group1) instead of
log2(Group1/Group2). Five files are flipped; one is renamed only.
Outputs: data/processed/atac/deseq2/ and deseq2_filtered/

### Script 1: ATAC-seq Peak Annotation
File: scripts/Script1_ATAC_annotation.R
Documentation: quarto/Script1_ATAC_annotation.qmd
Purpose: Annotates all ATAC-seq peaks with genomic features
using ChIPseeker and the Ensembl v106 EnsDb for Gallus gallus
GRCg6a (AnnotationHub record AH100636). All peaks annotated
regardless of significance.
Outputs: results/Script1_annotation/annotated/

### Script 2: ATAC-seq and RNA-seq Integration (coming soon)
### Script 3: ATAC-seq Visualisation (coming soon)
### Script 4: Integration Visualisation (coming soon)

## Reproducibility
Session information and MD5 checksums for each script run
are saved in reproducibility/
Package versions are locked in renv.lock
