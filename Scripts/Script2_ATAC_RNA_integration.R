# ============================================================
# Script: Script2_ATAC_RNA_integration.R
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
#   Integrates annotated ATAC-seq differential peaks with
#   RNA-seq differential expression data by matching Ensembl
#   gene IDs. Identifies genes with concordant or discordant
#   changes in chromatin accessibility and gene expression.
#   Calculates integration scores for downstream visualisation.
#
# Inputs:
#   results/Script1_annotation/annotated/[name]_annotated.tsv
#   data/processed/rna/[name]_deg_all.tsv
#
# Outputs:
#   results/Script2_integration/ATAC_padj0.05_RNA_FC1.25_padj0.05/
#     per_comparison/integration_[name].tsv
#     per_comparison/integration_[name].xlsx
#     per_comparison/venn_genes_[name].xlsx
#     summary/integration_summary.tsv
#     summary/integration_summary.xlsx
#   reproducibility/Script2_integration/sessionInfo/
#   reproducibility/Script2_integration/md5/
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
library(openxlsx)

cat("Working directory:", here(), "\n")

# Verify package versions match those used when script
# was written — warns but does not stop if mismatch found
expected_versions <- list(
  "here"       = "1.0.2",
  "data.table" = "1.18.4",
  "openxlsx"   = "4.2.8.1"
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

# 2. DEFINE PARAMETERS ----
# Significance thresholds defined as named variables
# Never hardcoded silently inside filtering code
# ATAC threshold: padj <= 0.05, no fold change filter
# RNA threshold: pre-applied by Novomagic at
#   |log2FC| >= 0.322 (FC >= 1.25) AND padj <= 0.05
atac_padj_threshold <- 0.05
rna_threshold_label <- "FC1.25_padj0.05"
results_folder      <- here(
  "results", "Script2_integration",
  paste0("ATAC_padj", atac_padj_threshold,
         "_RNA_", rna_threshold_label)
)

cat("\nATAC threshold: padj <=", atac_padj_threshold,
    "(no FC filter)\n")
cat("RNA threshold (pre-applied by Novomagic):",
    rna_threshold_label, "\n")
cat("Results folder:", results_folder, "\n")

# 3. CREATE OUTPUT FOLDER STRUCTURE ----
# per_comparison: integration TSV, Excel, Venn Excel
# summary: summary table across all comparisons
dir.create(here(results_folder, "per_comparison"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(here(results_folder, "summary"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(here("reproducibility", "Script2_integration",
                "sessionInfo"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(here("reproducibility", "Script2_integration",
                "md5"),
           recursive = TRUE, showWarnings = FALSE)
cat("Output folders created\n")

# 4. DEFINE COMPARISONS ----
# atac:    annotated ATAC file from Script 1
# rna:     RNA DEG file from Novomagic (pre-filtered)
# name:    comparison label for file names
# display: human-readable label for Excel sheets
# group:   biological grouping
comparisons <- list(
  
  list(atac    = here("results", "Script1_annotation",
                      "annotated",
                      "AM133_vs_AMScr_annotated.tsv"),
       rna     = here("data", "processed", "rna",
                      "AM133vsAMScr_deg_all.tsv"),
       name    = "AM133_vs_AMScr",
       display = "AM133 vs AMScr",
       group   = "Group1_vs_Scramble"),
  
  list(atac    = here("results", "Script1_annotation",
                      "annotated",
                      "AMAll_vs_AMScr_annotated.tsv"),
       rna     = here("data", "processed", "rna",
                      "AMAllvsAMScr_deg_all.tsv"),
       name    = "AMAll_vs_AMScr",
       display = "AMAll vs AMScr",
       group   = "Group1_vs_Scramble"),
  
  list(atac    = here("results", "Script1_annotation",
                      "annotated",
                      "AM1_206_vs_AMScr_annotated.tsv"),
       rna     = here("data", "processed", "rna",
                      "AM1_206vsAMScr_deg_all.tsv"),
       name    = "AM1_206_vs_AMScr",
       display = "AM1/206 vs AMScr",
       group   = "Group1_vs_Scramble"),
  
  list(atac    = here("results", "Script1_annotation",
                      "annotated",
                      "AM133_vs_AM1_206_annotated.tsv"),
       rna     = here("data", "processed", "rna",
                      "AM133vsAM1_206_deg_all.tsv"),
       name    = "AM133_vs_AM1_206",
       display = "AM133 vs AM1/206",
       group   = "Group2_vs_Each_Other"),
  
  list(atac    = here("results", "Script1_annotation",
                      "annotated",
                      "AM133_vs_AMAll_annotated.tsv"),
       rna     = here("data", "processed", "rna",
                      "AM133vsAMAll_deg_all.tsv"),
       name    = "AM133_vs_AMAll",
       display = "AM133 vs AMAll",
       group   = "Group2_vs_Each_Other"),
  
  list(atac    = here("results", "Script1_annotation",
                      "annotated",
                      "AM1_206_vs_AMAll_annotated.tsv"),
       rna     = here("data", "processed", "rna",
                      "AM1_206vsAMAll_deg_all.tsv"),
       name    = "AM1_206_vs_AMAll",
       display = "AM1/206 vs AMAll",
       group   = "Group2_vs_Each_Other")
)

# 5. MAIN INTEGRATION LOOP ----
# For each comparison:
# Step 5a: verify input files exist
# Step 5b: read annotated ATAC file
# Step 5c: apply ATAC significance filter
# Step 5d: read RNA file
# Step 5e: find overlapping Ensembl gene IDs
# Step 5f: print overlapping genes
# Step 5g: prepare slim data tables for merging
# Step 5h: merge ATAC and RNA by Ensembl gene ID
# Step 5i: assign concordance category
# Step 5j: calculate integration scores
# Step 5k: Pearson correlation
# Step 5l: save integration TSV
# Step 5m: save integration Excel
# Step 5n: save Venn gene list Excel
# Step 5o: store summary record

cat("\n============================================================\n")
cat("Starting integration...\n")
cat("============================================================\n")

# Storage for summary table
summary_records <- list()

# Storage for overlap genes — used in 3-way Venn in Script 4
overlap_genes_by_comp <- list()

for (comp in comparisons) {
  
  cat("\n========================================\n")
  cat("Comparison:", comp$name, "\n")
  cat("========================================\n")
  
  # ----------------------------------------------------------
  # Step 5a: Verify input files exist
  # Both ATAC and RNA files must exist before proceeding
  # stopifnot() stops immediately with clear error if missing
  # ----------------------------------------------------------
  stopifnot(
    "ATAC file not found — check results/Script1_annotation/" =
      file.exists(comp$atac),
    "RNA file not found — check data/processed/rna/" =
      file.exists(comp$rna)
  )
  cat("Input files verified\n")
  
  # ----------------------------------------------------------
  # Step 5b: Read annotated ATAC file
  # Full annotated file from Script 1 — all peaks, all columns
  # Contains: Chr, Start, End, Conc columns, log2FoldChange,
  # padj, pvalue, replicate counts, annotation, gene_chromosome,
  # geneStart, geneEnd, geneLength, geneStrand,
  # ensembl_gene_id, external_gene_name, gene_biotype,
  # transcriptId, transcriptBiotype, distanceToTSS
  # ----------------------------------------------------------
  atac <- fread(comp$atac, sep = "\t")
  cat("ATAC peaks loaded:", nrow(atac), "\n")
  
  # ----------------------------------------------------------
  # Step 5c: Apply ATAC significance filter
  # Filter to padj <= atac_padj_threshold (0.05)
  # No fold change filter — chromatin accessibility changes
  # are often modest in absolute fold change but biologically
  # meaningful, so no FC threshold is applied
  # ----------------------------------------------------------
  sig_atac <- atac[padj <= atac_padj_threshold]
  cat("Significant ATAC peaks (padj <=",
      atac_padj_threshold, "):", nrow(sig_atac), "\n")
  
  atac_genes <- unique(sig_atac$ensembl_gene_id)
  atac_genes <- atac_genes[!is.na(atac_genes) &
                             atac_genes != ""]
  cat("Unique ATAC Ensembl gene IDs:", length(atac_genes), "\n")
  
  # ----------------------------------------------------------
  # Step 5d: Read RNA file
  # Use all rows directly — pre-filtered by Novomagic at
  # |log2FC| >= 0.322 (FC >= 1.25) AND padj <= 0.05
  # No additional filtering applied here
  # ----------------------------------------------------------
  rna <- fread(comp$rna, sep = "\t")
  cat("RNA genes loaded:", nrow(rna), "\n")
  
  rna_genes <- unique(rna$gene_id)
  rna_genes <- rna_genes[!is.na(rna_genes) &
                           rna_genes != ""]
  cat("Unique RNA Ensembl gene IDs:", length(rna_genes), "\n")
  
  # ----------------------------------------------------------
  # Step 5e: Find overlapping Ensembl gene IDs
  # Match ensembl_gene_id (ATAC) vs gene_id (RNA)
  # Ensembl IDs used as matching key — stable identifiers
  # that do not change between annotation versions unlike
  # gene symbols which can be renamed between releases
  # ----------------------------------------------------------
  overlap_genes <- intersect(atac_genes, rna_genes)
  cat("Overlapping Ensembl gene IDs:", length(overlap_genes), "\n")
  
  # Store for 3-way Venn diagrams in Script 4
  overlap_genes_by_comp[[comp$name]] <- overlap_genes
  
  # Store summary record — updated with correlation later
  summary_records[[comp$name]] <- list(
    comparison = comp$name,
    display    = comp$display,
    group      = comp$group,
    atac_genes = length(atac_genes),
    rna_genes  = length(rna_genes),
    overlap    = length(overlap_genes),
    pearson_r  = NA_real_,
    pearson_p  = NA_real_
  )
  
  # Handle no overlap case
  if (length(overlap_genes) == 0) {
    cat("No overlapping genes — skipping integration\n")
    rm(atac, sig_atac, rna)
    gc()
    next
  }
  
  # ----------------------------------------------------------
  # Step 5f: Print overlapping genes
  # Shows Ensembl ID and gene names from both datasets
  # for immediate biological assessment
  # ----------------------------------------------------------
  cat("\nOverlapping genes:\n")
  for (gid in overlap_genes) {
    atac_name <- sig_atac[ensembl_gene_id == gid,
                          external_gene_name][1]
    rna_name  <- rna[gene_id == gid, gene_name][1]
    cat(" ", gid, "| ATAC:", atac_name,
        "| RNA:", rna_name, "\n")
  }
  
  # ----------------------------------------------------------
  # Step 5g: Prepare slim data tables for merging
  # Keep only columns needed for integration output
  # Rename with atac_/rna_ prefixes to avoid conflicts
  # gene_name_atac: gene symbol from ATAC annotation (EnsDb v106)
  # gene_name_rna:  gene symbol from RNA-seq (Novomagic/EnsDb v106)
  # Both kept for transparency — discrepancies reflect
  # ongoing gene nomenclature revision between sources
  # ----------------------------------------------------------
  
  # ATAC slim table
  atac_slim <- sig_atac[, .(
    Chr,
    Start,
    End,
    annotation,
    gene_chromosome,
    distanceToTSS,
    ensembl_gene_id,
    gene_name_atac  = external_gene_name,
    atac_log2FC     = log2FoldChange,
    atac_padj       = padj
  )]
  
  # RNA slim table
  rna_slim <- rna[, .(
    ensembl_gene_id = gene_id,
    gene_name_rna   = gene_name,
    gene_biotype    = gene_biotype,
    rna_log2FC      = log2FoldChange,
    rna_padj        = padj,
    rna_pvalue      = pvalue
  )]
  
  # ----------------------------------------------------------
  # Step 5h: Merge ATAC and RNA by Ensembl gene ID
  # allow.cartesian = TRUE is necessary because one gene
  # can have multiple ATAC peaks — the RNA row for a gene
  # is duplicated for each ATAC peak of that gene
  # This preserves ALL peak-gene associations — scientifically
  # correct since multiple peaks near one gene all contribute
  # to the chromatin accessibility landscape of that gene
  # ----------------------------------------------------------
  integrated <- merge(
    atac_slim,
    rna_slim,
    by              = "ensembl_gene_id",
    allow.cartesian = TRUE
  )
  
  cat("Integrated rows (all peaks x genes):",
      nrow(integrated), "\n")
  
  # Add comparison metadata
  integrated[, comparison := comp$name]
  integrated[, group      := comp$group]
  
  # ----------------------------------------------------------
  # Step 5i: Assign concordance category
  # Based on sign of atac_log2FC and rna_log2FC:
  # Concordant up:        both positive — gene activated
  # Concordant down:      both negative — gene repressed
  # Discordant (open/down): ATAC+ RNA- — complex regulation
  # Discordant (closed/up): ATAC- RNA+ — indirect regulation
  # ----------------------------------------------------------
  integrated[, concordance := fifelse(
    atac_log2FC > 0 & rna_log2FC > 0,
    "Concordant up",
    fifelse(
      atac_log2FC < 0 & rna_log2FC < 0,
      "Concordant down",
      fifelse(
        atac_log2FC > 0 & rna_log2FC < 0,
        "Discordant (open/down)",
        fifelse(
          atac_log2FC < 0 & rna_log2FC > 0,
          "Discordant (closed/up)",
          "Other"
        )
      )
    )
  )]
  
  cat("\nConcordance summary:\n")
  print(integrated[, .N, by = concordance])
  
  # ----------------------------------------------------------
  # Step 5j: Calculate integration scores
  # All scores are ranking and reference metrics only —
  # NOT statistical tests
  #
  # combined_score:          atac_log2FC x rna_log2FC
  #                          positive = concordant
  #                          negative = discordant
  # combined_pvalue:         atac_padj x rna_padj
  #                          product of both padj values
  # neg_log10_pval:          -log10(combined_pvalue)
  #                          Y axis for volcano plot in Script 4
  # weighted_score:          combined_score x neg_log10_pval
  #                          primary ranking metric for plots
  #                          rewards strong + significant genes
  # z_combined_score:        scale(atac) x scale(rna)
  #                          equal-weight normalised score
  # sig_score:               -log10(atac_padj) + -log10(rna_padj)
  #                          additive significance score
  # distance_weight:         exp(-min(|distToTSS|,100kb)/10kb)
  #                          ~0.95 at 500bp, ~0.007 at 50kb
  # distance_weighted_score: combined_score x distance_weight
  # ----------------------------------------------------------
  integrated[, `:=`(
    combined_score          = atac_log2FC * rna_log2FC,
    combined_pvalue         = atac_padj   * rna_padj,
    neg_log10_pval          = -log10(atac_padj * rna_padj),
    weighted_score          = (atac_log2FC * rna_log2FC) *
      (-log10(atac_padj * rna_padj)),
    z_combined_score        = scale(atac_log2FC)[, 1] *
      scale(rna_log2FC)[, 1],
    sig_score               = -log10(atac_padj) +
      -log10(rna_padj),
    distance_weight         = exp(
      -pmin(abs(distanceToTSS),
            100000) / 10000),
    distance_weighted_score = (atac_log2FC * rna_log2FC) *
      exp(-pmin(abs(distanceToTSS),
                100000) / 10000)
  )]
  
  # ----------------------------------------------------------
  # Step 5k: Pearson correlation
  # Measures linear relationship between atac_log2FC and
  # rna_log2FC across all overlapping gene-peak pairs
  # r value: +1 = perfect positive, 0 = none, -1 = negative
  # p value: probability the correlation is due to chance
  # Minimum 3 integrated rows required to compute correlation
  # If fewer than 3 rows exist store NA and skip
  # ----------------------------------------------------------
  if (nrow(integrated) >= 3) {
    correlation <- cor.test(integrated$atac_log2FC,
                            integrated$rna_log2FC)
    cat("\nPearson r:", round(correlation$estimate, 3),
        "| p =", format(correlation$p.value, digits = 3),
        "\n")
    summary_records[[comp$name]]$pearson_r <-
      round(correlation$estimate, 3)
    summary_records[[comp$name]]$pearson_p <-
      round(correlation$p.value, 4)
  } else {
    cat("\nPearson correlation skipped —",
        "fewer than 3 observations\n")
  }
  
  # ----------------------------------------------------------
  # Step 5l: Save integration TSV
  # Full merged table — all peaks, all columns
  # Primary input for Scripts 3 and 4
  # ----------------------------------------------------------
  tsv_out <- here(results_folder, "per_comparison",
                  paste0("integration_", comp$name, ".tsv"))
  fwrite(integrated, tsv_out, sep = "\t")
  cat("Saved TSV:", basename(tsv_out), "\n")
  
  # ----------------------------------------------------------
  # Step 5m: Save integration Excel
  # Same content as TSV — for browsing and sharing
  # ----------------------------------------------------------
  xlsx_out <- here(results_folder, "per_comparison",
                   paste0("integration_", comp$name, ".xlsx"))
  write.xlsx(integrated, xlsx_out)
  cat("Saved Excel:", basename(xlsx_out), "\n")
  
  # ----------------------------------------------------------
  # Step 5n: Save Venn gene list Excel
  # Multi-sheet Excel showing gene lists per category
  # Enables biological interpretation of overlap results
  # Sheet 1: Overlap_genes — full details of all overlaps
  # Sheet 2-5: stratified by concordance category
  # Sheet 6: ATAC_only — significant in ATAC not in RNA
  # Sheet 7: RNA_only  — significant in RNA not in ATAC
  # ----------------------------------------------------------
  
  # ATAC only genes
  atac_only_genes <- setdiff(atac_genes, rna_genes)
  atac_only_dt    <- sig_atac[
    ensembl_gene_id %in% atac_only_genes,
    .(ensembl_gene_id, gene_name_atac = external_gene_name,
      Chr, Start, End, annotation,
      gene_chromosome, distanceToTSS,
      atac_log2FC = log2FoldChange,
      atac_padj   = padj)]
  
  # RNA only genes
  rna_only_genes <- setdiff(rna_genes, atac_genes)
  rna_only_dt    <- rna[
    gene_id %in% rna_only_genes,
    .(ensembl_gene_id = gene_id,
      gene_name_rna  = gene_name,
      gene_biotype,
      rna_log2FC     = log2FoldChange,
      rna_padj       = padj)]
  
  wb <- createWorkbook()
  
  addWorksheet(wb, "Overlap_genes")
  writeData(wb, "Overlap_genes", integrated)
  
  addWorksheet(wb, "Concordant_up")
  writeData(wb, "Concordant_up",
            integrated[concordance == "Concordant up"])
  
  addWorksheet(wb, "Concordant_down")
  writeData(wb, "Concordant_down",
            integrated[concordance == "Concordant down"])
  
  addWorksheet(wb, "Discordant_open_down")
  writeData(wb, "Discordant_open_down",
            integrated[concordance == "Discordant (open/down)"])
  
  addWorksheet(wb, "Discordant_closed_up")
  writeData(wb, "Discordant_closed_up",
            integrated[concordance == "Discordant (closed/up)"])
  
  addWorksheet(wb, "ATAC_only")
  writeData(wb, "ATAC_only", atac_only_dt)
  
  addWorksheet(wb, "RNA_only")
  writeData(wb, "RNA_only", rna_only_dt)
  
  venn_out <- here(results_folder, "per_comparison",
                   paste0("venn_genes_", comp$name, ".xlsx"))
  saveWorkbook(wb, venn_out, overwrite = TRUE)
  cat("Saved Venn Excel:", basename(venn_out), "\n")
  
  # Free memory before next comparison
  rm(atac, sig_atac, rna, atac_slim, rna_slim,
     integrated, atac_only_dt, rna_only_dt, wb)
  gc()
}

# 6. SUMMARY TABLE ----
cat("\n============================================================\n")
cat("INTEGRATION SUMMARY\n")
cat("============================================================\n")

summary_dt <- rbindlist(
  lapply(summary_records, as.data.table),
  fill = TRUE
)
print(summary_dt)

# Save summary TSV
summary_tsv <- here(results_folder, "summary",
                    "integration_summary.tsv")
fwrite(summary_dt, summary_tsv, sep = "\t")
cat("Saved summary TSV:", basename(summary_tsv), "\n")

# Save summary Excel
summary_xlsx <- here(results_folder, "summary",
                     "integration_summary.xlsx")
write.xlsx(summary_dt, summary_xlsx)
cat("Saved summary Excel:", basename(summary_xlsx), "\n")

# 7. REPRODUCIBILITY RECORDS ----
cat("\nSaving reproducibility records...\n")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Session information
session_info <- capture.output(sessionInfo())
session_file <- here(
  "reproducibility", "Script2_integration",
  "sessionInfo",
  paste0("sessionInfo_Script2_", timestamp, ".txt")
)
writeLines(session_info, session_file)
cat("Session info saved:", basename(session_file), "\n")

# MD5 checksums for all input and output files
input_files <- c(
  sapply(comparisons, function(x) x$atac),
  sapply(comparisons, function(x) x$rna)
)
output_files <- c(
  list.files(here(results_folder, "per_comparison"),
             full.names = TRUE),
  list.files(here(results_folder, "summary"),
             full.names = TRUE)
)
all_files      <- c(input_files, output_files)
existing_files <- all_files[file.exists(all_files)]

md5_records <- data.table(
  file   = basename(existing_files),
  path   = existing_files,
  md5sum = sapply(existing_files, tools::md5sum,
                  USE.NAMES = FALSE),
  type   = ifelse(existing_files %in% input_files,
                  "input", "output"),
  script = "Script2_ATAC_RNA_integration",
  run_at = format(Sys.time())
)

md5_file <- here(
  "reproducibility", "Script2_integration",
  "md5",
  paste0("md5_Script2_", timestamp, ".csv")
)
fwrite(md5_records, md5_file)
cat("MD5 checksums saved:", basename(md5_file), "\n")

# Runtime
script_end <- Sys.time()
cat("\nScript completed:", format(script_end), "\n")
cat("Total runtime:",
    format(round(script_end - script_start, 2)), "\n")

# Save workspace
save.image(here("Script2_ATAC_RNA_integration.RData"))
cat("Workspace saved\n")
cat("Script 2 complete!\n")