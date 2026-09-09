# ============================================================
# Script: Script1_ATAC_annotation.R
# Project: Chromatin accessibility and gene expression
#          integration in chicken embryo antagomir study
# Author: Shahrzad Moradi Fard
# Email: shahrzadmf94@gmail.com
#        s.moradi-fard@uea.ac.uk
# GitHub: Shahrzad-94
# Date: 2026-05-30
# R version: 4.5.3
# Bioconductor version: 3.22
#
# Purpose:
#   Annotates all ATAC-seq peaks from the 6 corrected DiffBind
#   output files with genomic features using ChIPseeker and
#   the Ensembl v106 EnsDb for Gallus gallus GRCg6a (AH100636).
#   All peaks annotated regardless of significance — downstream
#   scripts apply their own filters as needed.
#
# Inputs:
#   data/processed/atac/deseq2/AM133_vs_AMScr_atac.tsv
#   data/processed/atac/deseq2/AMAll_vs_AMScr_atac.tsv
#   data/processed/atac/deseq2/AM1_206_vs_AMScr_atac.tsv
#   data/processed/atac/deseq2/AM133_vs_AM1_206_atac.tsv
#   data/processed/atac/deseq2/AM133_vs_AMAll_atac.tsv
#   data/processed/atac/deseq2/AM1_206_vs_AMAll_atac.tsv
#
# Outputs:
#   results/Script1_annotation/annotated/[name]_annotated.tsv
#   reproducibility/Script1_annotation/sessionInfo/
#   reproducibility/Script1_annotation/md5/
# ============================================================

# ------------------------------------------------------------
# REPRODUCIBILITY SETUP
# set.seed(): ensures identical results for any random process
# here(): verifies correct project root before anything runs
# Package version checks: warns if versions differ from those
#   used when this script was written
# Output folders: created here so script is fully self-contained
# ------------------------------------------------------------
set.seed(42)
script_start <- Sys.time()
cat("Script started:", format(script_start), "\n")
cat("R version:", R.version.string, "\n")
cat("Platform:", .Platform$OS.type, "\n")

# 1. DEPENDENCIES & CONFIGURATION ----
library(here)
library(data.table)
library(AnnotationHub)
library(ensembldb)
library(ChIPseeker)
library(GenomicRanges)
library(GenomicFeatures)

cat("Working directory:", here(), "\n")

# Verify package versions match those used when script
# was written — warns but does not stop if mismatch found
expected_versions <- list(
  "here"            = "1.0.2",
  "data.table"      = "1.18.4",
  "AnnotationHub"   = "4.0.0",
  "ensembldb"       = "2.34.0",
  "ChIPseeker"      = "1.46.1",
  "GenomicRanges"   = "1.62.1",
  "GenomicFeatures" = "1.62.0"
)
for (pkg in names(expected_versions)) {
  installed <- as.character(packageVersion(pkg))
  expected  <- expected_versions[[pkg]]
  if (installed != expected) {
    cat("WARNING: Package", pkg,
        "version", installed,
        "does not match expected", expected, "\n")
    cat("Run renv::restore() to restore correct versions\n")
  } else {
    cat("OK:", pkg, installed, "\n")
  }
}

# Create all output folders needed by this script
dir.create(here("results", "Script1_annotation", "annotated"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(here("reproducibility", "Script1_annotation",
                "sessionInfo"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(here("reproducibility", "Script1_annotation",
                "md5"),
           recursive = TRUE, showWarnings = FALSE)
cat("Output folders created\n")

# 2. LOAD ENSEMBL V106 ENSDB ----
# AH100636: Ensembl v106 EnsDb for Gallus gallus GRCg6a
# Downloaded from AnnotationHub and cached locally on first run
# Subsequent runs use the local cache — no internet needed
# This replaces both TxDb.Ggallus.UCSC.galGal6.refGene and
# org.Gg.eg.db used in previous pipeline versions, ensuring
# complete Ensembl v106 consistency throughout the pipeline
#
# NOTE on chromosome naming:
# seqlevelsStyle() is not available for Gallus gallus in
# GenomeInfoDb — this is a known limitation of the package.
# We handle the UCSC (chr1) vs Ensembl (1) naming mismatch
# manually: strip chr prefix before annotation, restore after.
cat("\nLoading Ensembl v106 EnsDb (AH100636)...\n")
ah       <- AnnotationHub()
edb_v106 <- ah[["AH100636"]]
cat("EnsDb loaded successfully\n")
cat("Organism:", metadata(edb_v106)[
  metadata(edb_v106)$name == "Organism", "value"], "\n")
cat("Genome build:", metadata(edb_v106)[
  metadata(edb_v106)$name == "genome_build", "value"], "\n")
cat("Ensembl version:", metadata(edb_v106)[
  metadata(edb_v106)$name == "ensembl_version", "value"], "\n")

# 3. BUILD GENE ANNOTATION LOOKUP TABLE ----
# Query all genes from EnsDb v106 once before the annotation
# loop — retrieves gene_id, gene_name, gene_biotype, and
# chromosome for all 24,356 genes in Ensembl v106.
# gene_chromosome is retrieved directly from EnsDb to ensure
# consistent Ensembl v106 sourcing throughout the pipeline.
# This lookup is merged with ChIPseeker output in the loop
# to add gene symbols, biotypes, and chromosome to each peak.
# gene_name is empty for unannotated genes — in these cases
# the Ensembl gene ID serves as the gene identifier,
# consistent with the RNA-seq annotation from Novomagic
# which also uses Ensembl v106 identifiers.
cat("\nBuilding gene name lookup table from EnsDb v106...\n")

all_genes_edb <- genes(
  edb_v106,
  columns = c("gene_id", "gene_name", "gene_biotype", "seq_name")
)

gene_lookup <- as.data.table(as.data.frame(all_genes_edb))
gene_lookup <- gene_lookup[, .(
  ensembl_gene_id    = gene_id,
  external_gene_name = gene_name,
  gene_biotype       = gene_biotype,
  gene_chromosome    = paste0("chr", as.character(seqnames))
)]

# Where gene_name is empty use ensembl_gene_id as placeholder
# This ensures every peak has a non-empty gene identifier
# for matching with RNA-seq data in Script 2
gene_lookup[
  is.na(external_gene_name) | external_gene_name == "",
  external_gene_name := ensembl_gene_id
]

cat("Total genes in lookup table:", nrow(gene_lookup), "\n")
cat("With gene symbol:",
    sum(gene_lookup$external_gene_name !=
          gene_lookup$ensembl_gene_id), "\n")
cat("Using Ensembl ID as identifier:",
    sum(gene_lookup$external_gene_name ==
          gene_lookup$ensembl_gene_id), "\n")

# 4. DEFINE FILE LIST ----
# One entry per comparison:
# input:  full corrected file from Script 0
# output: annotated file saved by this script
# name:   comparison label used in messages and file names
atac_files <- list(
  
  list(input  = here("data", "processed", "atac", "deseq2",
                     "AM133_vs_AMScr_atac.tsv"),
       output = here("results", "Script1_annotation",
                     "annotated",
                     "AM133_vs_AMScr_annotated.tsv"),
       name   = "AM133_vs_AMScr"),
  
  list(input  = here("data", "processed", "atac", "deseq2",
                     "AMAll_vs_AMScr_atac.tsv"),
       output = here("results", "Script1_annotation",
                     "annotated",
                     "AMAll_vs_AMScr_annotated.tsv"),
       name   = "AMAll_vs_AMScr"),
  
  list(input  = here("data", "processed", "atac", "deseq2",
                     "AM1_206_vs_AMScr_atac.tsv"),
       output = here("results", "Script1_annotation",
                     "annotated",
                     "AM1_206_vs_AMScr_annotated.tsv"),
       name   = "AM1_206_vs_AMScr"),
  
  list(input  = here("data", "processed", "atac", "deseq2",
                     "AM133_vs_AM1_206_atac.tsv"),
       output = here("results", "Script1_annotation",
                     "annotated",
                     "AM133_vs_AM1_206_annotated.tsv"),
       name   = "AM133_vs_AM1_206"),
  
  list(input  = here("data", "processed", "atac", "deseq2",
                     "AM133_vs_AMAll_atac.tsv"),
       output = here("results", "Script1_annotation",
                     "annotated",
                     "AM133_vs_AMAll_annotated.tsv"),
       name   = "AM133_vs_AMAll"),
  
  list(input  = here("data", "processed", "atac", "deseq2",
                     "AM1_206_vs_AMAll_atac.tsv"),
       output = here("results", "Script1_annotation",
                     "annotated",
                     "AM1_206_vs_AMAll_annotated.tsv"),
       name   = "AM1_206_vs_AMAll")
)

# 5. ANNOTATE EACH FILE ----
# For each file:
# Step 5a: verify input exists
# Step 5b: read corrected file
# Step 5c: strip chr prefix — converts UCSC to Ensembl style
# Step 5d: create GRanges object
# Step 5e: annotate with ChIPseeker using EnsDb v106
# Step 5f: convert back to data table
# Step 5g: restore chr prefix in output
# Step 5h: add gene annotations from lookup table
#          (gene name, biotype, and chromosome)
# Step 5i: rename and reorder columns
# Step 5j: save annotated TSV

cat("\n============================================================\n")
cat("Annotating ATAC-seq peaks...\n")
cat("============================================================\n")

# Storage for summary table
summary_records <- list()

for (af in atac_files) {
  
  cat("\n========================================\n")
  cat("Comparison:", af$name, "\n")
  cat("========================================\n")
  
  # ----------------------------------------------------------
  # Step 5a: Verify input file exists
  # stopifnot() stops immediately with clear error if missing
  # ----------------------------------------------------------
  stopifnot(
    "Input file not found — check data/processed/atac/deseq2/" =
      file.exists(af$input)
  )
  cat("Input verified:", basename(af$input), "\n")
  
  # ----------------------------------------------------------
  # Step 5b: Read corrected file
  # Contains all peaks with standardised column names
  # from Script 0: log2FoldChange, padj, pvalue
  # ----------------------------------------------------------
  dt <- fread(af$input, sep = "\t")
  cat("Rows loaded:", nrow(dt), "\n")
  cat("Columns:", paste(colnames(dt), collapse = ", "), "\n")
  
  # ----------------------------------------------------------
  # Step 5c: Strip chr prefix from chromosome names
  # EnsDb v106 uses Ensembl-style names (1, 2... Z, W)
  # DiffBind files use UCSC-style names (chr1, chr2... chrZ)
  # seqlevelsStyle() not available for Gallus gallus in
  # GenomeInfoDb — handled manually here
  # Unplaced scaffolds (chrUn_NW_...) become Un_NW_...
  # These will not match EnsDb and receive NA annotation
  # which is scientifically correct for unplaced scaffolds
  # ----------------------------------------------------------
  dt[, Chr_ensembl := sub("^chr", "", Chr)]
  cat("Chr prefix stripped for EnsDb matching\n")
  cat("Example:", dt$Chr[1], "->", dt$Chr_ensembl[1], "\n")
  
  # ----------------------------------------------------------
  # Step 5d: Create GRanges object
  # GRanges is the standard Bioconductor genomic data structure
  # required by ChIPseeker's annotatePeak() function
  # All original columns carried as metadata so they are
  # preserved through annotation and appear in output
  # ----------------------------------------------------------
  peaks_gr <- makeGRangesFromDataFrame(
    dt,
    keep.extra.columns  = TRUE,
    seqnames.field      = "Chr_ensembl",
    start.field         = "Start",
    end.field           = "End"
  )
  cat("GRanges created:", length(peaks_gr), "peaks\n")
  # After creating peaks_gr, explicitly set seqlevels
  # to match EnsDb seqlevels order exactly
  # This prevents factor level confusion in ChIPseeker
  peaks_gr <- keepSeqlevels(
    peaks_gr,
    intersect(seqlevels(peaks_gr), seqlevels(edb_v106)),
    pruning.mode = "coarse"
  )
  seqlevels(peaks_gr) <- seqlevels(edb_v106)[
    seqlevels(edb_v106) %in% seqlevels(peaks_gr)]
  # ----------------------------------------------------------
  # Step 5e: Annotate peaks with ChIPseeker
  # annotatePeak() assigns each peak to nearest genomic feature
  # TxDb = edb_v106: uses Ensembl v106 transcript models
  #   ensuring annotation matches RNA-seq reference exactly
  # tssRegion = c(-3000, 3000): peaks within 3kb of a TSS
  #   are classified as promoter peaks — standard definition
  #   in published ATAC-seq studies
  # verbose = FALSE: suppresses internal progress messages
  # ----------------------------------------------------------
  cat("Annotating with ChIPseeker + Ensembl v106...\n")
  peak_anno <- annotatePeak(
    peaks_gr,
    TxDb      = edb_v106,
    tssRegion = c(-3000, 3000),
    verbose   = FALSE
  )
  
  # ----------------------------------------------------------
  # Step 5f: Convert annotation back to data table
  # as.data.frame() extracts the full annotated data frame
  # New columns added by ChIPseeker:
  #   annotation:        genomic feature category
  #   geneChr:           chromosome of nearest gene
  #   geneStart/End:     coordinates of nearest gene
  #   geneLength:        length of nearest gene
  #   geneStrand:        strand of nearest gene
  #   geneId:            Ensembl gene ID of nearest gene
  #   transcriptId:      Ensembl transcript ID
  #   transcriptBiotype: transcript biotype
  #   distanceToTSS:     distance in bp to nearest TSS
  # ----------------------------------------------------------
  anno_dt <- as.data.table(as.data.frame(peak_anno))
  cat("Annotation complete:", nrow(anno_dt), "peaks\n")
  
  # ----------------------------------------------------------
  # Step 5g: Restore chr prefix to peak chromosome column
  # seqnames contains Ensembl-style chromosome names (1, 2...)
  # restored to UCSC-style (chr1, chr2...) to match the
  # original DiffBind input files throughout the pipeline
  # gene_chromosome is added from the lookup table in Step 5h
  # ----------------------------------------------------------
  anno_dt[, seqnames := paste0("chr", seqnames)]
  cat("Chr prefix restored in output\n")
  # ----------------------------------------------------------
  # Step 5h: Add gene annotations from lookup table
  # Merge by geneId (Ensembl gene ID from ChIPseeker)
  # Adds external_gene_name, gene_biotype, and
  # gene_chromosome columns from Ensembl v106 EnsDb
  # all.x = TRUE keeps all peaks even without a gene match
  # Peaks on unplaced scaffolds and truly unannotated peaks
  # receive NA — this is scientifically correct and expected
  # ----------------------------------------------------------
  anno_dt <- merge(
    anno_dt,
    gene_lookup,
    by.x  = "geneId",
    by.y  = "ensembl_gene_id",
    all.x = TRUE
  )
  
  # ----------------------------------------------------------
  # Step 5i: Rename and reorder columns
  # Rename seqnames -> Chr, start -> Start, end -> End
  # to match original DiffBind column naming convention
  # Rename geneId -> ensembl_gene_id for consistency with
  # RNA-seq column naming used in Script 2 integration
  # Reorder: original DiffBind columns first, then annotation
  # ----------------------------------------------------------
  setnames(anno_dt,
           old = c("seqnames", "start",  "end",   "geneId"),
           new = c("Chr",      "Start",  "End",   "ensembl_gene_id"))
  
  # Define column order — original columns first
  original_cols <- c("Chr", "Start", "End",
                     "Conc", "log2FoldChange", "padj", "pvalue")
  
  # Identify which original columns exist in this file
  # (Conc columns vary by comparison)
  conc_cols <- grep("^Conc", colnames(anno_dt), value = TRUE)
  replicate_cols <- grep(
    "^(AM133|AMAll|AM1_206|AMScr)_[0-9]",
    colnames(anno_dt), value = TRUE)
  
  # Annotation columns from ChIPseeker and EnsDb v106 lookup
  anno_cols <- c("annotation", "gene_chromosome",
                 "geneStart", "geneEnd",
                 "geneLength", "geneStrand",
                 "ensembl_gene_id", "external_gene_name",
                 "gene_biotype", "transcriptId",
                 "transcriptBiotype", "distanceToTSS")
  
  # Build final column order
  final_cols <- c("Chr", "Start", "End",
                  conc_cols,
                  "log2FoldChange", "padj", "pvalue",
                  replicate_cols,
                  anno_cols)
  
  # Keep only columns that exist
  final_cols <- final_cols[final_cols %in% colnames(anno_dt)]
  anno_dt    <- anno_dt[, ..final_cols]
  
  # ----------------------------------------------------------
  # Step 5j: Save annotated TSV
  # ALL peaks saved — no filtering applied
  # Used by Scripts 2, 3, and 4 which apply their own filters
  # ----------------------------------------------------------
  fwrite(anno_dt, af$output, sep = "\t")
  cat("Saved annotated TSV:", basename(af$output), "\n")
  cat("Total peaks saved:", nrow(anno_dt), "\n")
  
  # Summary statistics for this comparison
  n_total     <- nrow(anno_dt)
  n_standard  <- sum(!grepl("^chrUn", anno_dt$Chr))
  n_unplaced  <- sum(grepl("^chrUn", anno_dt$Chr))
  n_with_name <- sum(!is.na(anno_dt$external_gene_name))
  n_no_name   <- sum(is.na(anno_dt$external_gene_name))
  
  cat("Standard chromosomes:", n_standard, "\n")
  cat("Unplaced scaffolds:", n_unplaced, "\n")
  cat("Peaks with gene identifier:", n_with_name,
      "(", round(n_with_name / n_total * 100, 1), "%)\n")
  cat("Peaks without gene identifier:", n_no_name,
      "(", round(n_no_name / n_total * 100, 1), "%)\n")
  
  summary_records[[af$name]] <- list(
    comparison   = af$name,
    total_peaks  = n_total,
    standard_chr = n_standard,
    unplaced     = n_unplaced,
    with_gene    = n_with_name,
    without_gene = n_no_name,
    pct_annotated = round(n_with_name / n_total * 100, 1)
  )
  
  rm(dt, peaks_gr, peak_anno, anno_dt)
  gc()
}

# 6. SUMMARY TABLE ----
cat("\n============================================================\n")
cat("SCRIPT 1 COMPLETE — ANNOTATION SUMMARY\n")
cat("============================================================\n")

summary_dt <- rbindlist(
  lapply(summary_records, as.data.table),
  fill = TRUE
)
print(summary_dt)

# 7. REPRODUCIBILITY RECORDS ----
cat("\nSaving reproducibility records...\n")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Session information
session_info <- capture.output(sessionInfo())
session_file <- here(
  "reproducibility", "Script1_annotation",
  "sessionInfo",
  paste0("sessionInfo_Script1_", timestamp, ".txt")
)
writeLines(session_info, session_file)
cat("Session info saved:", basename(session_file), "\n")

# MD5 checksums for all input and output files
all_files <- c(
  sapply(atac_files, function(x) x$input),
  sapply(atac_files, function(x) x$output)
)
existing_files <- all_files[file.exists(all_files)]

md5_records <- data.table(
  file   = basename(existing_files),
  path   = existing_files,
  md5sum = sapply(existing_files, tools::md5sum,
                  USE.NAMES = FALSE),
  type   = c(rep("input",  length(atac_files)),
             rep("output", length(atac_files))
  )[seq_along(existing_files)],
  script = "Script1_ATAC_annotation",
  run_at = format(Sys.time())
)

md5_file <- here(
  "reproducibility", "Script1_annotation",
  "md5",
  paste0("md5_Script1_", timestamp, ".csv")
)
fwrite(md5_records, md5_file)
cat("MD5 checksums saved:", basename(md5_file), "\n")

# Runtime
script_end <- Sys.time()
cat("\nScript completed:", format(script_end), "\n")
cat("Total runtime:",
    format(round(script_end - script_start, 2)), "\n")

# Save workspace
save.image(here("Script1_ATAC_annotation.RData"))
cat("Workspace saved\n")
cat("Script 1 complete!\n")