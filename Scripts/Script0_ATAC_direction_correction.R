# ============================================================
# Script: Script0_ATAC_direction_correction.R
# Project: Chromatin accessibility and gene expression
#          integration in chicken embryo antagomir study
# Author: Shahrzad Moradi Fard
# Email: shahrzadmf94@gmail.com
#        s.moradi-fard@uea.ac.uk
# GitHub: Shahrzad-94
# Date: 2026-05-29
# R version: 4.5.3
# Bioconductor version: not required for this script
#
# Purpose:
#   Corrects fold change direction in raw DiffBind ATAC-seq
#   output files. DiffBind calculates Fold as
#   log2(Group2/Group1) instead of log2(Group1/Group2).
#   This script flips the direction for 5 of the 6 files
#   and renames the 6th to reflect its true direction.
#   All columns are standardised to DESeq2 naming convention.
#   Filtered versions (padj <= 0.05) are also saved.
#
# Inputs:
#   data/raw/atac/deseq2/DEseq2_AM1_206_v_AM133.tsv
#   data/raw/atac/deseq2/DEseq2_AM1_206_v_AMAll.tsv
#   data/raw/atac/deseq2/DEseq2_AM1_206_v_AMScr.tsv
#   data/raw/atac/deseq2/DEseq2_AM133_v_AMAll.tsv
#   data/raw/atac/deseq2/DEseq2_AM133_v_AMScr.tsv
#   data/raw/atac/deseq2/DEseq2_AMAll_v_AMScr.tsv
#
# Outputs:
#   data/processed/atac/deseq2/[name]_atac.tsv (6 files)
#   data/processed/atac/deseq2_filtered/[name]_atac_padj0.05.tsv
#   reproducibility/Script0_corrections/sessionInfo/
#   reproducibility/Script0_corrections/md5/
# ============================================================

# ------------------------------------------------------------
# REPRODUCIBILITY SETUP
# set.seed(): ensures identical results for any random process
# here(): verifies correct project root before anything runs
# Package version checks: warns if versions differ from those
#   used when this script was written — does not stop script
#   because correction is deterministic and version-independent
# Output folders: created here so script is fully self-contained
#   recursive = TRUE  — creates all parent folders if needed
#   showWarnings = FALSE — no error if folder already exists,
#   allowing safe reruns without manual cleanup
# ------------------------------------------------------------
set.seed(42)
script_start <- Sys.time()
cat("Script started:", format(script_start), "\n")
cat("R version:", R.version.string, "\n")
cat("Platform:", .Platform$OS.type, "\n")

# 1. DEPENDENCIES & CONFIGURATION ----
library(here)
library(data.table)

cat("Working directory:", here(), "\n")

# Verify package versions match those used when script
# was written — warns but does not stop if mismatch found
expected_versions <- list(
  "here"       = "1.0.2",
  "data.table" = "1.18.4"
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
# All other folders created by their respective scripts
dir.create(here("data", "processed", "atac", "deseq2"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(here("data", "processed", "atac", "deseq2_filtered"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(here("reproducibility", "Script0_corrections",
                "sessionInfo"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(here("reproducibility", "Script0_corrections",
                "md5"),
           recursive = TRUE, showWarnings = FALSE)
cat("Output folders created\n")

# 2. PARAMETERS ----
# Significance threshold for filtered output files
# Defined as a named variable — never hardcoded silently
# No fold change filter applied — padj only
# This matches the ATAC threshold agreed for this study
atac_padj_threshold <- 0.05

cat("\nATAC significance threshold: padj <=",
    atac_padj_threshold, "\n")
cat("Fold change filter: none\n")

# 3. DEFINE ATAC FILE LIST ----
# One entry per comparison with all processing information:
# input:        raw DiffBind TSV in data/raw/atac/deseq2/
# output_full:  corrected full file in data/processed/atac/deseq2/
# output_filt:  filtered file in data/processed/atac/deseq2_filtered/
# name:         comparison label used in messages and file names
# flip:         TRUE = multiply log2FoldChange by -1
#               FALSE = rename only (direction already correct)
# conc_group1:  Conc column for first-named group
#               used in verification — higher value = more
#               accessible in group1 = should be positive FC
# conc_group2:  Conc column for second-named group
#
# NOTE on AM133_vs_AM1_206:
#   DiffBind file is named AM1_206_v_AM133 but calculates
#   Fold as log2(AM133/AM1_206) — so positive Fold already
#   means higher in AM133. RNA file AM133vsAM1_206_deg_all
#   also has positive log2FC = higher in AM133. They already
#   agree — no flip needed. File is renamed only to reflect
#   its true biological direction.
atac_files <- list(
  
  list(input       = here("data", "raw", "atac", "deseq2",
                          "DEseq2_AM133_v_AMScr.tsv"),
       output_full = here("data", "processed", "atac", "deseq2",
                          "AM133_vs_AMScr_atac.tsv"),
       output_filt = here("data", "processed", "atac",
                          "deseq2_filtered",
                          "AM133_vs_AMScr_atac_padj0.05.tsv"),
       name        = "AM133_vs_AMScr",
       flip        = TRUE,
       conc_group1 = "Conc_AM133",
       conc_group2 = "Conc_AMScr"),
  
  list(input       = here("data", "raw", "atac", "deseq2",
                          "DEseq2_AMAll_v_AMScr.tsv"),
       output_full = here("data", "processed", "atac", "deseq2",
                          "AMAll_vs_AMScr_atac.tsv"),
       output_filt = here("data", "processed", "atac",
                          "deseq2_filtered",
                          "AMAll_vs_AMScr_atac_padj0.05.tsv"),
       name        = "AMAll_vs_AMScr",
       flip        = TRUE,
       conc_group1 = "Conc_AMAll",
       conc_group2 = "Conc_AMScr"),
  
  list(input       = here("data", "raw", "atac", "deseq2",
                          "DEseq2_AM1_206_v_AMScr.tsv"),
       output_full = here("data", "processed", "atac", "deseq2",
                          "AM1_206_vs_AMScr_atac.tsv"),
       output_filt = here("data", "processed", "atac",
                          "deseq2_filtered",
                          "AM1_206_vs_AMScr_atac_padj0.05.tsv"),
       name        = "AM1_206_vs_AMScr",
       flip        = TRUE,
       conc_group1 = "Conc_AM1_206",
       conc_group2 = "Conc_AMScr"),
  
  list(input       = here("data", "raw", "atac", "deseq2",
                          "DEseq2_AM1_206_v_AM133.tsv"),
       output_full = here("data", "processed", "atac", "deseq2",
                          "AM133_vs_AM1_206_atac.tsv"),
       output_filt = here("data", "processed", "atac",
                          "deseq2_filtered",
                          "AM133_vs_AM1_206_atac_padj0.05.tsv"),
       name        = "AM133_vs_AM1_206",
       flip        = FALSE,
       conc_group1 = "Conc_AM133",
       conc_group2 = "Conc_AM1_206"),
  
  list(input       = here("data", "raw", "atac", "deseq2",
                          "DEseq2_AM133_v_AMAll.tsv"),
       output_full = here("data", "processed", "atac", "deseq2",
                          "AM133_vs_AMAll_atac.tsv"),
       output_filt = here("data", "processed", "atac",
                          "deseq2_filtered",
                          "AM133_vs_AMAll_atac_padj0.05.tsv"),
       name        = "AM133_vs_AMAll",
       flip        = TRUE,
       conc_group1 = "Conc_AM133",
       conc_group2 = "Conc_AMAll"),
  
  list(input       = here("data", "raw", "atac", "deseq2",
                          "DEseq2_AM1_206_v_AMAll.tsv"),
       output_full = here("data", "processed", "atac", "deseq2",
                          "AM1_206_vs_AMAll_atac.tsv"),
       output_filt = here("data", "processed", "atac",
                          "deseq2_filtered",
                          "AM1_206_vs_AMAll_atac_padj0.05.tsv"),
       name        = "AM1_206_vs_AMAll",
       flip        = TRUE,
       conc_group1 = "Conc_AM1_206",
       conc_group2 = "Conc_AMAll")
)

# 4. PROCESS EACH ATAC FILE ----
# For each file:
# Step 4a: verify input file exists before reading
# Step 4b: read raw file and print dimensions
# Step 4c: rename columns to DESeq2 standard naming FIRST
#          before any other operation — cleaner and consistent
# Step 4d: verify fold change direction BEFORE correction
#          using Conc columns as ground truth
# Step 4e: apply fold change flip if needed
# Step 4f: verify fold change direction AFTER correction
# Step 4g: save full corrected TSV — all peaks, no filter
# Step 4h: save filtered TSV — padj <= threshold only

cat("\n============================================================\n")
cat("Processing ATAC files...\n")
cat("============================================================\n")

# Storage for summary table built after all files processed
summary_records <- list()

for (af in atac_files) {
  
  cat("\n========================================\n")
  cat("Comparison:", af$name, "\n")
  cat("========================================\n")
  
  # ----------------------------------------------------------
  # Step 4a: Verify input file exists
  # stopifnot() stops immediately with clear error message
  # if the raw file is missing — catches path problems
  # before any processing begins rather than crashing mid-loop
  # ----------------------------------------------------------
  stopifnot(
    "Input file not found — check data/raw/atac/deseq2/" =
      file.exists(af$input)
  )
  cat("Input file verified:", basename(af$input), "\n")
  
  # ----------------------------------------------------------
  # Step 4b: Read raw DiffBind TSV file
  # fread() is used for speed — files have ~70,000-95,000 rows
  # sep = "\t" explicitly specifies tab separation
  # Print dimensions as confirmation of successful load
  # ----------------------------------------------------------
  dt <- fread(af$input, sep = "\t")
  cat("Rows loaded:", nrow(dt), "\n")
  cat("Columns:", paste(colnames(dt), collapse = ", "), "\n")
  
  # ----------------------------------------------------------
  # Step 4c: Rename columns to DESeq2 standard naming
  # Done FIRST before any other operation for consistency
  # Fold    -> log2FoldChange : standard DESeq2 column name
  # FDR     -> padj           : standard adjusted p-value name
  # p.value -> pvalue         : removes dot to avoid backtick
  #                             issues when referencing in R
  # All subsequent steps use the new standardised names
  # ----------------------------------------------------------
  setnames(dt,
           old = c("Fold",            "FDR",  "p-value"),
           new = c("log2FoldChange",  "padj", "pvalue"))
  cat("Columns renamed:",
      "Fold -> log2FoldChange,",
      "FDR -> padj,",
      "p-value -> pvalue\n")
  
  # ----------------------------------------------------------
  # Step 4d: Verify fold change direction BEFORE correction
  # Logic: find a peak where Group1 concentration is clearly
  # higher than Group2 concentration (2x threshold ensures
  # an unambiguous test case). For this peak:
  # - If flip = TRUE: log2FoldChange should be NEGATIVE
  #   (DiffBind calculated log2(Group2/Group1) so higher
  #   Group1 gives negative value — confirming flip is needed)
  # - If flip = FALSE: log2FoldChange should be POSITIVE
  #   (direction already correct — confirming no flip needed)
  # Conc columns contain mean log2 normalised read
  # concentrations per group — reliable ground truth
  # ----------------------------------------------------------
  cat("\nVerification BEFORE correction:\n")
  
  test_peak <- dt[get(af$conc_group1) >
                    get(af$conc_group2) * 2][1]
  
  if (nrow(test_peak) == 0) {
    cat("WARNING: No clear test peak found for verification\n")
  } else {
    cat("Test peak:", test_peak$Chr, ":",
        test_peak$Start, "-", test_peak$End, "\n")
    cat(af$conc_group1, "=",
        round(test_peak[[af$conc_group1]], 3), "\n")
    cat(af$conc_group2, "=",
        round(test_peak[[af$conc_group2]], 3), "\n")
    cat("log2FoldChange =",
        round(test_peak$log2FoldChange, 3), "\n")
    
    if (af$flip) {
      cat("Expected BEFORE flip: NEGATIVE\n")
      cat("Result:",
          ifelse(test_peak$log2FoldChange < 0,
                 "NEGATIVE - flip confirmed needed",
                 "WARNING: POSITIVE - check file direction"),
          "\n")
    } else {
      cat("Expected (no flip): POSITIVE\n")
      cat("Result:",
          ifelse(test_peak$log2FoldChange > 0,
                 "POSITIVE - direction already correct",
                 "WARNING: NEGATIVE - check file direction"),
          "\n")
    }
  }
  
  # Store test peak fold change for after-correction check
  fc_before <- test_peak$log2FoldChange
  
  # ----------------------------------------------------------
  # Step 4e: Apply fold change flip if needed
  # Only applied when flip = TRUE
  # Multiplies log2FoldChange by -1 to reverse direction
  # so positive log2FoldChange = higher in first-named group
  # For flip = FALSE (AM133_vs_AM1_206): direction already
  # correct — column renamed only, no calculation applied
  # ----------------------------------------------------------
  if (af$flip) {
    dt[, log2FoldChange := log2FoldChange * -1]
    cat("\nFold change direction corrected (x -1)\n")
  } else {
    cat("\nNo flip applied — direction already correct\n")
    cat("File renamed to reflect true biological direction\n")
  }
  
  # ----------------------------------------------------------
  # Step 4f: Verify fold change direction AFTER correction
  # Find same test peak by chromosome coordinates
  # For flipped files: log2FoldChange should now be POSITIVE
  # For non-flipped file: log2FoldChange should still be
  # POSITIVE — unchanged from before
  # Signs should be opposite for flipped files
  # ----------------------------------------------------------
  cat("\nVerification AFTER correction:\n")
  
  fc_after <- dt[Chr   == test_peak$Chr  &
                   Start == test_peak$Start &
                   End   == test_peak$End,
                 log2FoldChange][1]
  
  cat("log2FoldChange after =", round(fc_after, 3), "\n")
  
  if (af$flip) {
    cat("Expected: POSITIVE (opposite sign to before)\n")
    cat("Result:",
        ifelse(fc_after > 0 &
                 sign(fc_before) != sign(fc_after),
               "POSITIVE and OPPOSITE SIGN - correction confirmed",
               "WARNING: check correction"),
        "\n")
  } else {
    cat("Expected: POSITIVE (unchanged)\n")
    cat("Result:",
        ifelse(fc_after > 0,
               "POSITIVE - direction confirmed correct",
               "WARNING: check file"),
        "\n")
  }
  
  # ----------------------------------------------------------
  # Step 4g: Save full corrected TSV
  # ALL peaks saved — no filtering applied
  # Preserves complete dataset for downstream annotation
  # in Script 1 which needs all peaks regardless of significance
  # fwrite() used for speed — consistent with fread()
  # ----------------------------------------------------------
  fwrite(dt, af$output_full, sep = "\t")
  cat("\nSaved full corrected TSV:", basename(af$output_full), "\n")
  cat("Total peaks saved:", nrow(dt), "\n")
  
  # ----------------------------------------------------------
  # Step 4h: Save filtered TSV
  # Filter to padj <= atac_padj_threshold (0.05)
  # No fold change filter — padj only as agreed
  # These filtered files are used directly in Script 2
  # integration and Script 3 visualisation
  # Filtering applied AFTER saving full file to ensure
  # complete dataset is always preserved
  # ----------------------------------------------------------
  dt_filtered <- dt[padj <= atac_padj_threshold]
  fwrite(dt_filtered, af$output_filt, sep = "\t")
  cat("Saved filtered TSV:", basename(af$output_filt), "\n")
  cat("Peaks passing padj <=", atac_padj_threshold, ":",
      nrow(dt_filtered), "of", nrow(dt),
      paste0("(", round(nrow(dt_filtered) / nrow(dt) * 100, 1),
             "%)\n"))
  
  # Store summary record for this file
  summary_records[[af$name]] <- list(
    comparison       = af$name,
    action           = ifelse(af$flip,
                              "Flipped + renamed",
                              "Renamed only"),
    total_peaks      = nrow(dt),
    significant_peaks = nrow(dt_filtered),
    pct_significant  = round(
      nrow(dt_filtered) / nrow(dt) * 100, 1),
    fc_before        = round(fc_before, 3),
    fc_after         = round(fc_after, 3)
  )
  
  rm(dt, dt_filtered, test_peak)
  gc()
}

# 5. SUMMARY TABLE ----
# Complete record of all corrections made
# Printed to console and serves as permanent run record
# alongside the sessionInfo and MD5 files saved below
cat("\n============================================================\n")
cat("SCRIPT 0 COMPLETE — SUMMARY OF ALL CORRECTIONS\n")
cat("============================================================\n")

summary_dt <- rbindlist(
  lapply(summary_records, as.data.table),
  fill = TRUE
)
print(summary_dt)

cat("\nPositive log2FoldChange now means:\n")
positive_means <- data.table(
  comparison = c("AM133_vs_AMScr",
                 "AMAll_vs_AMScr",
                 "AM1_206_vs_AMScr",
                 "AM133_vs_AM1_206",
                 "AM133_vs_AMAll",
                 "AM1_206_vs_AMAll"),
  positive_log2FC_means = c(
    "More open in AM133",
    "More open in AMAll",
    "More open in AM1_206",
    "More open in AM133",
    "More open in AM133",
    "More open in AM1_206"
  )
)
print(positive_means)

cat("\nAll corrected files saved in:\n")
cat("  Full:     data/processed/atac/deseq2/\n")
cat("  Filtered: data/processed/atac/deseq2_filtered/\n")
cat("  Filter applied: padj <=", atac_padj_threshold,
    "(no FC filter)\n")

# 6. REPRODUCIBILITY RECORDS ----
# Timestamped so each run creates a new file
# Never overwrites previous runs — permanent record per run
# MD5 checksums verify file integrity for all inputs
# and outputs — anyone can confirm they have identical files
cat("\nSaving reproducibility records...\n")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Session information — complete R environment record
session_info <- capture.output(sessionInfo())
session_file <- here(
  "reproducibility", "Script0_corrections",
  "sessionInfo",
  paste0("sessionInfo_Script0_", timestamp, ".txt")
)
writeLines(session_info, session_file)
cat("Session info saved:", basename(session_file), "\n")

# MD5 checksums for all input and output files
all_files <- c(
  # Input files
  sapply(atac_files, function(x) x$input),
  # Full corrected output files
  sapply(atac_files, function(x) x$output_full),
  # Filtered output files
  sapply(atac_files, function(x) x$output_filt)
)

existing_files <- all_files[file.exists(all_files)]

md5_records <- data.table(
  file    = basename(existing_files),
  path    = existing_files,
  md5sum  = sapply(existing_files, tools::md5sum,
                   USE.NAMES = FALSE),
  type    = c(
    rep("input",    length(atac_files)),
    rep("output_full",     length(atac_files)),
    rep("output_filtered", length(atac_files))
  )[seq_along(existing_files)],
  script  = "Script0_ATAC_direction_correction",
  run_at  = format(Sys.time())
)

md5_file <- here(
  "reproducibility", "Script0_corrections",
  "md5",
  paste0("md5_Script0_", timestamp, ".csv")
)
fwrite(md5_records, md5_file)
cat("MD5 checksums saved:", basename(md5_file), "\n")

# Runtime
script_end <- Sys.time()
cat("\nScript completed:", format(script_end), "\n")
cat("Total runtime:",
    format(round(script_end - script_start, 2)), "\n")

# Save workspace
save.image(here("Script0_ATAC_direction_correction.RData"))
cat("Workspace saved\n")
cat("Script 0 complete!\n")



