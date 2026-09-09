# ============================================================
# Script: Script4_integration_visualisation.R
# Project: Chromatin accessibility and gene expression
#          integration in chicken embryo antagomir study
# Author: Shahrzad Moradi Fard
# Email: shahrzadmf94@gmail.com
#        s.moradi-fard@uea.ac.uk
# GitHub: Shahrzad-94
# Date: 2026-05-31
# R version: 4.5.3
# Bioconductor version: 3.22
#
# Purpose:
#   Generates all integration visualisations from the
#   ATAC-seq and RNA-seq integration files produced by
#   Script 2. Produces heatmaps, scatter plots, volcano
#   plots, Venn diagrams, distance histograms, GO/KEGG
#   enrichment bar charts, and protein interaction networks.
#
# Input:
#   results/Script2_integration/
#     ATAC_padj0.05_RNA_FC1.25_padj0.05/per_comparison/
#     integration_[name].tsv
#   results/Script2_integration/
#     ATAC_padj0.05_RNA_FC1.25_padj0.05/summary/
#     integration_summary.tsv
#   results/Script1_annotation/annotated/
#     [name]_annotated.tsv (for GO/KEGG background)
#
# Output:
#   results/Script4_integration_visualisation/
#     heatmap/, scatter/, volcano/, venn/,
#     distance_histogram/, go_kegg/, network/
#   cache/STRINGdb_cache/
#   reproducibility/Script4_integration_visualisation/
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

# 1. LOAD LIBRARIES ----
library(here)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(ggpubr)
library(openxlsx)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(enrichplot)
library(org.Gg.eg.db)
library(STRINGdb)
library(igraph)
library(ggraph)
library(tidygraph)

cat("Working directory:", here(), "\n")

# Verify package versions
expected_versions <- list(
  "here"           = "1.0.2",
  "data.table"     = "1.18.4",
  "ggplot2"        = "4.0.3",
  "ggrepel"        = "0.9.8",
  "ggforce"        = "0.5.0",
  "ggpubr"         = "0.6.3",
  "openxlsx"       = "4.2.8.1",
  "dplyr"          = "1.2.1",
  "ComplexHeatmap" = "2.26.1",
  "clusterProfiler"= "4.18.4",
  "org.Gg.eg.db"   = "3.22.0",
  "STRINGdb"       = "2.22.0",
  "igraph"         = "2.3.1",
  "ggraph"         = "2.2.2",
  "tidygraph"      = "1.3.1"
)
for (pkg in names(expected_versions)) {
  installed <- as.character(packageVersion(pkg))
  expected  <- expected_versions[[pkg]]
  if (installed != expected) {
    cat("WARNING: Package", pkg,
        "version", installed,
        "does not match expected", expected, "\n")
  } else {
    cat("OK:", pkg, installed, "\n")
  }
}

# 2. DEFINE COLOUR PALETTES ----
# All colours from the Okabe-Ito colourblind-safe palette.
# Defined once and referenced throughout — guarantees
# consistency across all plots of the same type.
# Concordance colours chosen to be maximally distinct
# from heatmap ATAC colours (Vermillion, Blue) and
# heatmap RNA colours (Orange, Bluish Green).

concordance_colours <- c(
  "Concordant up"          = "#FF8C00",  # Orange
  "Concordant down"        = "#8A2BE2",  # Purple
  "Discordant (open/down)" = "#8B4513",  # Brown
  "Discordant (closed/up)" = "#808080",  # Gray
  "Other"                  = "#000000"   # Black
)

# GO/KEGG ontology colours
ontology_colours <- c(
  "BP"   = "#009E73",  # Bluish Green
  "MF"   = "#D55E00",  # Vermillion
  "CC"   = "#0072B2",  # Blue
  "KEGG" = "#E69F00"   # Orange
)

# 3. DEFINE PARAMETERS AND FILE LIST ----
integration_dir <- here(
  "results", "Script2_integration",
  "ATAC_padj0.05_RNA_FC1.25_padj0.05"
)
atac_padj_threshold <- 0.05
out_base  <- here("results",
                  "Script4_integration_visualisation")
cache_dir <- here("cache", "STRINGdb_cache")

comparisons <- list(
  list(name        = "AM133_vs_AMScr",
       display     = "AM133 vs AMScr",
       group       = "Group1_vs_Scramble",
       integration = here(integration_dir,
                          "per_comparison",
                          "integration_AM133_vs_AMScr.tsv"),
       annotated   = here("results",
                          "Script1_annotation",
                          "annotated",
                          "AM133_vs_AMScr_annotated.tsv")),
  list(name        = "AMAll_vs_AMScr",
       display     = "AMAll vs AMScr",
       group       = "Group1_vs_Scramble",
       integration = here(integration_dir,
                          "per_comparison",
                          "integration_AMAll_vs_AMScr.tsv"),
       annotated   = here("results",
                          "Script1_annotation",
                          "annotated",
                          "AMAll_vs_AMScr_annotated.tsv")),
  list(name        = "AM1_206_vs_AMScr",
       display     = "AM1/206 vs AMScr",
       group       = "Group1_vs_Scramble",
       integration = here(integration_dir,
                          "per_comparison",
                          "integration_AM1_206_vs_AMScr.tsv"),
       annotated   = here("results",
                          "Script1_annotation",
                          "annotated",
                          "AM1_206_vs_AMScr_annotated.tsv")),
  list(name        = "AM133_vs_AM1_206",
       display     = "AM133 vs AM1/206",
       group       = "Group2_vs_Each_Other",
       integration = here(integration_dir,
                          "per_comparison",
                          "integration_AM133_vs_AM1_206.tsv"),
       annotated   = here("results",
                          "Script1_annotation",
                          "annotated",
                          "AM133_vs_AM1_206_annotated.tsv")),
  list(name        = "AM133_vs_AMAll",
       display     = "AM133 vs AMAll",
       group       = "Group2_vs_Each_Other",
       integration = here(integration_dir,
                          "per_comparison",
                          "integration_AM133_vs_AMAll.tsv"),
       annotated   = here("results",
                          "Script1_annotation",
                          "annotated",
                          "AM133_vs_AMAll_annotated.tsv")),
  list(name        = "AM1_206_vs_AMAll",
       display     = "AM1/206 vs AMAll",
       group       = "Group2_vs_Each_Other",
       integration = here(integration_dir,
                          "per_comparison",
                          "integration_AM1_206_vs_AMAll.tsv"),
       annotated   = here("results",
                          "Script1_annotation",
                          "annotated",
                          "AM1_206_vs_AMAll_annotated.tsv"))
)

# 4. CREATE OUTPUT FOLDER STRUCTURE ----
for (sub in c("heatmap", "scatter", "volcano", "venn",
              "distance_histogram",
              "go_kegg/no_background",
              "go_kegg/with_background",
              "network")) {
  dir.create(here(out_base, sub),
             recursive = TRUE, showWarnings = FALSE)
}
dir.create(cache_dir,
           recursive = TRUE, showWarnings = FALSE)
dir.create(here("reproducibility",
                "Script4_integration_visualisation",
                "sessionInfo"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(here("reproducibility",
                "Script4_integration_visualisation",
                "md5"),
           recursive = TRUE, showWarnings = FALSE)
cat("Output folders created\n")

# 5. HELPER FUNCTIONS ----

# ----------------------------------------------------------
# save_plot(): saves a ggplot object as PNG and PDF.
# All plots saved at 300 DPI for publication quality.
# ----------------------------------------------------------
save_plot <- function(p, base_path, width, height) {
  ggsave(paste0(base_path, ".png"), p,
         width = width, height = height,
         dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(paste0(base_path, ".pdf"), p,
         width = width, height = height,
         bg = "white", limitsize = FALSE)
  cat("Saved:", basename(base_path), ".png and .pdf\n")
}

# ----------------------------------------------------------
# draw_venn_2way(): pairwise Venn diagram.
# Shows unique gene counts — each gene counted once
# regardless of how many ATAC peaks it has.
# Blue for ATAC, Vermillion for RNA — consistent across
# all pairwise Venn diagrams.
# ----------------------------------------------------------
draw_venn_2way <- function(n_atac, n_rna, n_overlap,
                           label_atac, label_rna,
                           title, base_path) {
  n_atac_only <- n_atac - n_overlap
  n_rna_only  <- n_rna  - n_overlap
  
  circles <- data.frame(
    x0    = c(-0.5,  0.5),
    y0    = c( 0,    0),
    r     = c( 1,    1),
    group = c("ATAC", "RNA")
  )
  labels <- data.frame(
    x     = c(-1.1,  0,    1.1),
    y     = c( 0,    0,    0),
    label = c(as.character(n_atac_only),
              as.character(n_overlap),
              as.character(n_rna_only))
  )
  cond_labels <- data.frame(
    x     = c(-1.2,  1.2),
    y     = c( 1.3,  1.3),
    label = c(label_atac, label_rna)
  )
  p <- ggplot() +
    geom_circle(
      data      = circles,
      aes(x0 = x0, y0 = y0, r = r, fill = group),
      alpha     = 0.3,
      color     = "gray40",
      linewidth = 0.8
    ) +
    geom_text(data = labels,
              aes(x = x, y = y, label = label),
              size = 6, fontface = "bold") +
    geom_text(data = cond_labels,
              aes(x = x, y = y, label = label),
              size = 4.5, fontface = "bold",
              color = "gray20") +
    scale_fill_manual(
      values = c("ATAC" = "#0072B2",
                 "RNA"  = "#D55E00")
    ) +
    coord_fixed() +
    labs(title = title) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title      = element_text(hjust  = 0.5,
                                     face   = "bold",
                                     size   = 13,
                                     margin = margin(b = 15)),
      plot.margin     = margin(20, 20, 20, 20)
    )
  save_plot(p, base_path, w = 6, h = 6)
}

# ----------------------------------------------------------
# draw_venn_3way(): 3-way Venn diagram for group-level
# overlap visualisation. Takes three gene vectors and
# calculates all seven region counts. Uses three distinct
# Okabe-Ito colours for the three circles.
# ----------------------------------------------------------
draw_venn_3way <- function(genes1, genes2, genes3,
                           labels, title, base_path) {
  only1  <- length(setdiff(genes1, union(genes2, genes3)))
  only2  <- length(setdiff(genes2, union(genes1, genes3)))
  only3  <- length(setdiff(genes3, union(genes1, genes2)))
  only12 <- length(setdiff(intersect(genes1, genes2),
                           genes3))
  only13 <- length(setdiff(intersect(genes1, genes3),
                           genes2))
  only23 <- length(setdiff(intersect(genes2, genes3),
                           genes1))
  all3   <- length(intersect(intersect(genes1, genes2),
                             genes3))
  
  circles <- data.frame(
    x0    = c( 0,   -0.6,  0.6),
    y0    = c( 0.5, -0.3, -0.3),
    r     = c( 1,    1,    1),
    group = c("C1", "C2", "C3")
  )
  numbers <- data.frame(
    x     = c( 0,   -1.0,  1.0,
               -0.6,  0.6,  0,   0),
    y     = c( 1.2, -0.7, -0.7,
               0.4,  0.4, -0.6, 0.1),
    label = c(as.character(only1),
              as.character(only2),
              as.character(only3),
              as.character(only12),
              as.character(only13),
              as.character(only23),
              as.character(all3))
  )
  cond_labels <- data.frame(
    x     = c( 0,    -1.8,   1.8),
    y     = c( 2.2,  -1.6,  -1.6),
    label = labels
  )
  p <- ggplot() +
    geom_circle(
      data      = circles,
      aes(x0 = x0, y0 = y0, r = r, fill = group),
      alpha     = 0.2,
      color     = "gray40",
      linewidth = 0.8
    ) +
    geom_text(data = numbers,
              aes(x = x, y = y, label = label),
              size = 5, fontface = "bold") +
    geom_text(data = cond_labels,
              aes(x = x, y = y, label = label),
              size = 3.8, fontface = "bold",
              color = "gray20") +
    scale_fill_manual(
      values = c("C1" = "#D55E00",
                 "C2" = "#0072B2",
                 "C3" = "#009E73")
    ) +
    coord_fixed(xlim = c(-3, 3),
                ylim = c(-2.5, 3)) +
    labs(title = title) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title      = element_text(hjust  = 0.5,
                                     face   = "bold",
                                     size   = 13,
                                     margin = margin(b = 15)),
      plot.margin     = margin(20, 20, 20, 20)
    )
  save_plot(p, base_path, w = 7, h = 8)
}

# 6. PRE-PASS: FIND GLOBAL SCALE MAXIMA ----
# Before generating any heatmaps or networks, a single
# pass over all integration files finds the maximum
# absolute fold change values. These are used as fixed
# symmetric colour scale limits across ALL heatmaps so
# the same colour always represents the same fold change
# magnitude — enabling direct cross-comparison.
cat("\n========================================\n")
cat("Pre-pass: finding global scale maxima...\n")
cat("========================================\n")

atac_max_global <- 0
rna_max_global  <- 0

for (comp in comparisons) {
  if (!file.exists(comp$integration)) next
  dt <- fread(comp$integration, sep = "\t")
  if (nrow(dt) == 0) next
  atac_max_global <- max(atac_max_global,
                         quantile(abs(dt$atac_log2FC),
                                  0.95, na.rm = TRUE))
  rna_max_global  <- max(rna_max_global,
                         quantile(abs(dt$rna_log2FC),
                                  0.95, na.rm = TRUE))
  rm(dt); gc()
}

cat("Global ATAC max |log2FC|:", round(atac_max_global, 3),
    "\n")
cat("Global RNA max |log2FC|:", round(rna_max_global, 3),
    "\n")

# 7. DEFINE HEATMAP COLOUR SCALES ----
# Defined after pre-pass so limits are based on actual data.
# ATAC: Vermillion positive, white zero, Blue negative.
# RNA: Orange positive, white zero, Bluish Green negative.
# These colour pairs do not overlap with concordance colours
# so heatmap columns and concordance sidebar are visually
# distinct.
atac_col <- colorRamp2(
  c(-atac_max_global, 0, atac_max_global),
  c("#4575B4", "white", "#D73027")
)
rna_col <- colorRamp2(
  c(-rna_max_global, 0, rna_max_global),
  c("#228B22", "white", "#FFD700")
)

cat("Heatmap colour scales defined\n")

# 8. BUILD GO/KEGG BACKGROUND GENE SET ----
# Background = all unique genes near significant ATAC peaks
# across all six comparisons. This is the biologically
# appropriate background for enrichment analysis — it
# represents the set of genes in accessible chromatin
# regions in this experiment, correcting for the fact
# that the input genes come from an already-accessible
# chromatin context.
cat("\n========================================\n")
cat("Building GO/KEGG background gene set...\n")
cat("========================================\n")

all_bg_genes <- c()
for (comp in comparisons) {
  atac <- fread(comp$annotated, sep = "\t")
  sig  <- atac[padj <= atac_padj_threshold]
  genes <- sig$external_gene_name[
    !is.na(sig$external_gene_name) &
      sig$external_gene_name != ""]
  all_bg_genes <- c(all_bg_genes, genes)
  rm(atac, sig); gc()
}
all_bg_genes <- unique(all_bg_genes)
cat("Total unique background genes:", length(all_bg_genes),
    "\n")

# Convert background to Entrez IDs
background_entrez <- tryCatch({
  bitr(all_bg_genes,
       fromType = "SYMBOL",
       toType   = "ENTREZID",
       OrgDb    = org.Gg.eg.db)
}, error = function(e) NULL)

if (!is.null(background_entrez)) {
  cat("Background Entrez IDs:",
      nrow(background_entrez), "\n")
} else {
  cat("Warning: background conversion failed\n")
}

# 9. SET UP STRINGDB ----
# STRINGdb version 11.5, Gallus gallus (species 9031).
# Score threshold 400 = medium confidence interactions —
# the standard threshold used in published network analyses.
# Cache directory stores downloaded interaction data
# locally to avoid repeated downloads on reruns.
cat("\n========================================\n")
cat("Setting up STRINGdb...\n")
cat("========================================\n")

options(timeout = 600)
string_db <- STRINGdb$new(
  version         = "11.5",
  species         = 9031,
  score_threshold = 400,
  input_directory = cache_dir
)
cat("STRINGdb ready\n")

# 10. MAIN PER-COMPARISON LOOP ----
cat("\n========================================\n")
cat("Starting per-comparison visualisation...\n")
cat("========================================\n")

# Named list to store network objects for all comparisons
# Populated during main loop, plotted after global
# maximum degree is determined in Block 11
network_list <- list()

# Storage for overlap genes — used in 3-way Venn diagrams
overlap_genes_by_comp <- list()

# Read summary file for Venn diagram counts
summary_dt <- fread(
  here(integration_dir, "summary",
       "integration_summary.tsv"),
  sep = "\t"
)

for (comp in comparisons) {
  
  cat("\n========================================\n")
  cat("Comparison:", comp$name, "\n")
  cat("========================================\n")
  
  # Check integration file exists
  if (!file.exists(comp$integration)) {
    cat("No integration file — skipping\n")
    overlap_genes_by_comp[[comp$name]] <- character(0)
    next
  }
  
  # Read integration file
  dt <- fread(comp$integration, sep = "\t")
  cat("Integrated rows loaded:", nrow(dt), "\n")
  
  # Store unique gene names for 3-way Venn
  overlap_genes_by_comp[[comp$name]] <-
    unique(dt$gene_name_atac)
  
  # Skip plots if no rows
  if (nrow(dt) == 0) {
    cat("No overlapping genes — skipping plots\n")
    next
  }
  
  # ----------------------------------------------------------
  # Step 10a: Heatmap
  # One row per gene — select most significant peak per gene
  # (lowest atac_padj, absolute atac_log2FC as tiebreaker).
  # Rows ordered by concordance then by descending absolute
  # atac_log2FC within each concordance category.
  # Fixed colour scales across all heatmaps.
  # ----------------------------------------------------------
  cat("Generating heatmap...\n")
  
  dt_gene <- dt[order(atac_padj, -abs(atac_log2FC))][
    !duplicated(gene_name_atac)]
  dt_gene <- dt_gene[order(
    concordance, -abs(atac_log2FC))]
  
  atac_mat <- matrix(
    dt_gene$atac_log2FC,
    nrow     = nrow(dt_gene),
    ncol     = 1,
    dimnames = list(dt_gene$gene_name_atac, "ATAC-seq")
  )
  rna_mat <- matrix(
    dt_gene$rna_log2FC,
    nrow     = nrow(dt_gene),
    ncol     = 1,
    dimnames = list(dt_gene$gene_name_atac, "RNA-seq")
  )
  
  # Concordance annotation sidebar
  conc_present <- concordance_colours[
    names(concordance_colours) %in%
      unique(dt_gene$concordance)]
  
  ha_right <- rowAnnotation(
    Concordance = dt_gene$concordance,
    col         = list(
      Concordance = conc_present),
    annotation_legend_param = list(
      Concordance = list(title = "Concordance"))
  )
  
  ht <- Heatmap(
    atac_mat,
    name            = "ATAC\nlog2FC",
    col             = atac_col,
    cluster_rows    = FALSE,
    cluster_columns = FALSE,
    show_row_names  = TRUE,
    row_names_side  = "left",
    row_names_gp    = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 10,
                           fontface = "bold"),
    width           = unit(1.5, "cm")
  ) +
    Heatmap(
      rna_mat,
      name            = "RNA\nlog2FC",
      col             = rna_col,
      cluster_rows    = FALSE,
      cluster_columns = FALSE,
      show_row_names  = FALSE,
      column_names_gp = gpar(fontsize = 10,
                             fontface = "bold"),
      width           = unit(1.5, "cm")
    ) +
    ha_right
  
  ht_height <- max(800, nrow(dt_gene) * 75 + 400)
  
  png(here(out_base, "heatmap",
           paste0("heatmap_", comp$name, ".png")),
      width = 2400, height = ht_height, res = 300)
  draw(ht,
       column_title    = paste(
         "ATAC-seq & RNA-seq Integration:",
         comp$display),
       column_title_gp = gpar(fontsize = 12,
                              fontface = "bold"),
       padding         = unit(c(5, 5, 5, 5), "mm"))
  dev.off()
  
  pdf(here(out_base, "heatmap",
           paste0("heatmap_", comp$name, ".pdf")),
      width  = 6,
      height = max(4, nrow(dt_gene) * 0.25 + 2))
  draw(ht,
       column_title    = paste(
         "ATAC-seq & RNA-seq Integration:",
         comp$display),
       column_title_gp = gpar(fontsize = 12,
                              fontface = "bold"),
       padding         = unit(c(5, 5, 5, 5), "mm"))
  dev.off()
  
  cat("Saved heatmap:", comp$name, "\n")
  
  # ----------------------------------------------------------
  # Step 10b: Scatter plot
  # ATAC log2FC (X) vs RNA log2FC (Y) coloured by
  # concordance. Point size proportional to joint
  # significance. Linear regression line shows overall
  # relationship between chromatin and expression changes.
  # Two versions: unlabelled and labelled with top genes
  # by absolute weighted score (one label per gene).
  # ----------------------------------------------------------
  cat("Generating scatter plots...\n")
  
  p_scatter_base <- ggplot(
    dt,
    aes(x     = atac_log2FC,
        y     = rna_log2FC,
        color = concordance,
        size  = neg_log10_pval)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "gray60", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "gray60", linewidth = 0.5) +
    geom_smooth(method = "lm", se = FALSE,
                color = "gray40", linewidth = 0.8,
                linetype = "dashed") +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = concordance_colours,
                       name   = "Concordance") +
    scale_size_continuous(
      range = c(2, 8),
      name  = "-log10(ATAC padj\nx RNA padj)") +
    labs(
      title = paste("ATAC-seq & RNA-seq Integration:",
                    comp$display),
      x     = "ATAC-seq log2FoldChange",
      y     = "RNA-seq log2FoldChange"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title       = element_text(hjust = 0.5,
                                      face  = "bold",
                                      size  = 12),
      legend.position  = "right",
      panel.grid.minor = element_blank()
    )
  
  save_plot(
    p_scatter_base,
    here(out_base, "scatter",
         paste0("scatter_", comp$name, "_nolabels")),
    width = 8, height = 6
  )
  
  # Labelled version — top genes by weighted score
  label_dt <- dt[order(-abs(weighted_score))][
    !duplicated(gene_name_atac)]
  
  p_scatter_labelled <- p_scatter_base +
    geom_text_repel(
      data        = label_dt,
      aes(label   = gene_name_atac),
      size        = 3,
      max.overlaps = 30,
      box.padding = 0.4,
      fontface    = "bold",
      color       = "black",
      show.legend = FALSE
    )
  
  save_plot(
    p_scatter_labelled,
    here(out_base, "scatter",
         paste0("scatter_", comp$name, "_labelled")),
    width = 8, height = 6
  )
  
  # ----------------------------------------------------------
  # Step 10c: Volcano plot
  # Combined score (X) vs joint -log10 p-value (Y).
  # Combined score = atac_log2FC x rna_log2FC — positive
  # values indicate concordant, negative indicate discordant.
  # Two versions: unlabelled and labelled.
  # ----------------------------------------------------------
  cat("Generating volcano plots...\n")
  
  p_volcano_base <- ggplot(
    dt,
    aes(x     = combined_score,
        y     = neg_log10_pval,
        color = concordance)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "gray50", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "gray50", linewidth = 0.5) +
    geom_point(alpha = 0.8, size = 2.5) +
    scale_color_manual(values = concordance_colours,
                       name   = "Concordance") +
    labs(
      title = paste("ATAC-seq & RNA-seq Integration:",
                    comp$display),
      x     = "Combined Score (ATAC x RNA log2FoldChange)",
      y     = "-log10(ATAC padj x RNA padj)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title       = element_text(hjust = 0.5,
                                      face  = "bold",
                                      size  = 12),
      legend.position  = "right",
      panel.grid.minor = element_blank()
    )
  
  save_plot(
    p_volcano_base,
    here(out_base, "volcano",
         paste0("volcano_", comp$name, "_nolabels")),
    width = 8, height = 6
  )
  
  p_volcano_labelled <- p_volcano_base +
    geom_text_repel(
      data        = label_dt,
      aes(label   = gene_name_atac),
      size        = 3,
      max.overlaps = 30,
      box.padding = 0.4,
      fontface    = "bold",
      color       = "black",
      show.legend = FALSE
    )
  
  save_plot(
    p_volcano_labelled,
    here(out_base, "volcano",
         paste0("volcano_", comp$name, "_labelled")),
    width = 8, height = 6
  )
  
  # ----------------------------------------------------------
  # Step 10d: Distance histogram
  # Shows distance from peak to nearest TSS for all
  # integrated peak-gene pairs, stratified by concordance.
  # All rows used (not deduplicated) — each peak-gene
  # pair contributes independently to the distribution.
  # X axis in kilobases, continuous, 30 bins.
  # Free Y scales allow each concordance category to be
  # assessed independently of differing group sizes.
  # ----------------------------------------------------------
  cat("Generating distance histograms...\n")
  
  p_dist <- ggplot(
    dt,
    aes(x    = abs(distanceToTSS) / 1000,
        fill = concordance)
  ) +
    geom_histogram(bins = 30, color = "white",
                   linewidth = 0.2) +
    facet_wrap(~concordance, scales = "free_y",
               nrow = 1) +
    scale_fill_manual(values = concordance_colours) +
    labs(
      title    = paste("Peak Distance to TSS:",
                       comp$display),
      subtitle = "All integrated peak-gene pairs",
      x        = "Distance to TSS (kb)",
      y        = "Number of peaks"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(hjust = 0.5,
                                      face  = "bold"),
      plot.subtitle    = element_text(hjust = 0.5,
                                      size  = 9,
                                      color = "gray40"),
      legend.position  = "none",
      strip.text       = element_text(face = "bold",
                                      size = 9),
      panel.grid.minor = element_blank()
    )
  
  save_plot(
    p_dist,
    here(out_base, "distance_histogram",
         paste0("dist_hist_", comp$name)),
    width = 12, height = 6
  )
  
  # ----------------------------------------------------------
  # Step 10e: Pairwise Venn diagram
  # Shows unique gene counts from the summary table —
  # each gene counted once regardless of peak count.
  # Blue = ATAC-seq genes, Vermillion = RNA-seq genes.
  # ----------------------------------------------------------
  cat("Generating pairwise Venn diagram...\n")
  
  comp_summary <- summary_dt[comparison == comp$name]
  
  if (nrow(comp_summary) > 0) {
    draw_venn_2way(
      n_atac     = comp_summary$atac_genes,
      n_rna      = comp_summary$rna_genes,
      n_overlap  = comp_summary$overlap,
      label_atac = "ATAC-seq",
      label_rna  = "RNA-seq",
      title      = paste("Gene Overlap:", comp$display),
      base_path  = here(out_base, "venn",
                        paste0("venn_", comp$name))
    )
  }
  
  # ----------------------------------------------------------
  # Step 10f: GO/KEGG enrichment — two versions
  # Version A: no background (whole genome default)
  # Version B: custom background (all genes near significant
  #   ATAC peaks across all comparisons)
  # Both versions: Count >= 2, pvalue <= 0.05
  # Gene symbols used as input — confirmed identical between
  # ATAC and RNA annotations for all overlapping genes.
  # ----------------------------------------------------------
  cat("Running GO/KEGG enrichment...\n")
  
  gene_symbols <- unique(dt$gene_name_atac)
  gene_symbols <- gene_symbols[
    !is.na(gene_symbols) & gene_symbols != ""]
  
  cat("  Genes for enrichment:", length(gene_symbols), "\n")
  
  if (length(gene_symbols) >= 3) {
    
    entrez_ids <- tryCatch({
      bitr(gene_symbols,
           fromType = "SYMBOL",
           toType   = "ENTREZID",
           OrgDb    = org.Gg.eg.db)
    }, error = function(e) NULL)
    
    if (!is.null(entrez_ids) && nrow(entrez_ids) > 0) {
      
      # --------------------------------------------------
      # go_kegg_plot(): runs enrichment for one version,
      # combines ontologies, filters, and saves plot.
      # universe: NULL = whole genome, character vector =
      #   custom background Entrez IDs.
      # subfolder: "no_background" or "with_background"
      # subtitle_bg: describes background used
      # --------------------------------------------------
      go_kegg_plot <- function(entrez_vec, universe,
                               subfolder, subtitle_bg) {
        
        results_list <- list()
        
        for (ont in c("BP", "MF", "CC")) {
          res <- tryCatch({
            enrichGO(
              gene          = entrez_vec,
              universe      = universe,
              OrgDb         = org.Gg.eg.db,
              ont           = ont,
              pAdjustMethod = "BH",
              pvalueCutoff  = 1.0,
              qvalueCutoff  = 1.0,
              readable      = TRUE
            )
          }, error = function(e) NULL)
          
          if (!is.null(res)) {
            df <- as.data.frame(res)
            df <- df[df$pvalue <= 0.05 &
                       df$Count >= 2, ]
            if (nrow(df) > 0) {
              df$ontology    <- ont
              df$neg_log10_p <- -log10(df$pvalue)
              results_list[[ont]] <- df[,
                                        c("Description", "pvalue",
                                          "neg_log10_p", "ontology", "Count")]
            }
          }
        }
        
        # KEGG — no universe parameter supported
        res_kegg <- tryCatch({
          enrichKEGG(
            gene          = entrez_vec,
            organism      = "gga",
            pAdjustMethod = "BH",
            pvalueCutoff  = 1.0,
            qvalueCutoff  = 1.0
          )
        }, error = function(e) NULL)
        
        if (!is.null(res_kegg)) {
          df_k <- as.data.frame(res_kegg)
          df_k <- df_k[df_k$pvalue <= 0.05 &
                         df_k$Count >= 2, ]
          if (nrow(df_k) > 0) {
            df_k$ontology    <- "KEGG"
            df_k$neg_log10_p <- -log10(df_k$pvalue)
            results_list[["KEGG"]] <- df_k[,
                                           c("Description", "pvalue",
                                             "neg_log10_p", "ontology", "Count")]
          }
        }
        
        if (length(results_list) == 0) {
          cat("  No terms found for:", subfolder, "\n")
          return(invisible(NULL))
        }
        
        combined_df <- do.call(rbind, results_list)
        combined_df$ontology <- factor(
          combined_df$ontology,
          levels = c("BP", "MF", "CC", "KEGG"))
        combined_df <- combined_df[
          !duplicated(combined_df$Description), ]
        combined_df <- combined_df[
          order(combined_df$ontology,
                combined_df$neg_log10_p), ]
        combined_df$Description <- factor(
          combined_df$Description,
          levels = combined_df$Description)
        
        n_terms <- nrow(combined_df)
        plot_h  <- max(5, n_terms * 0.35 + 2)
        
        p_go <- ggplot(
          combined_df,
          aes(x    = neg_log10_p,
              y    = Description,
              fill = ontology)
        ) +
          geom_bar(stat = "identity", width = 0.7) +
          geom_text(
            aes(label = paste0("n=", Count)),
            hjust = -0.1, size = 3, color = "gray30"
          ) +
          scale_fill_manual(
            values = ontology_colours,
            name   = "Ontology",
            labels = c(
              "BP"   = "Biological Process",
              "MF"   = "Molecular Function",
              "CC"   = "Cellular Component",
              "KEGG" = "KEGG Pathway"
            )
          ) +
          scale_x_continuous(
            expand = expansion(mult = c(0, 0.2))
          ) +
          labs(
            title    = paste("GO/KEGG Enrichment:",
                             comp$display),
            subtitle = paste0(subtitle_bg,
                              " | p \u2264 0.05",
                              " | count \u2265 2"),
            x        = expression(-log[10](p-value)),
            y        = NULL
          ) +
          coord_cartesian(clip = "off") +
          theme_bw(base_size = 11) +
          theme(
            plot.title    = element_text(hjust = 0.5,
                                         face  = "bold",
                                         size  = 12),
            plot.subtitle = element_text(hjust = 0.5,
                                         size  = 8,
                                         color = "gray40"),
            axis.text.y   = element_text(size  = 9),
            legend.position    = "right",
            panel.grid.minor   = element_blank(),
            panel.grid.major.y = element_blank(),
            plot.margin        = margin(5, 40, 5, 5)
          )
        
        ggsave(
          here(out_base, "go_kegg", subfolder,
               paste0("go_kegg_", comp$name,
                      "_", subfolder, ".png")),
          p_go, width = 12, height = plot_h,
          dpi = 300, bg = "white", limitsize = FALSE
        )
        ggsave(
          here(out_base, "go_kegg", subfolder,
               paste0("go_kegg_", comp$name,
                      "_", subfolder, ".pdf")),
          p_go, width = 12, height = plot_h,
          bg = "white", limitsize = FALSE
        )
        cat("  Saved GO/KEGG:", subfolder, "\n")
      }
      
      # Version A: no background (whole genome)
      go_kegg_plot(
        entrez_vec  = entrez_ids$ENTREZID,
        universe    = NULL,
        subfolder   = "no_background",
        subtitle_bg = "Background: whole genome (default)"
      )
      
      # Version B: custom background
      if (!is.null(background_entrez)) {
        go_kegg_plot(
          entrez_vec  = entrez_ids$ENTREZID,
          universe    = background_entrez$ENTREZID,
          subfolder   = "with_background",
          subtitle_bg = paste0(
            "Background: all genes near significant ",
            "ATAC peaks (padj \u2264 0.05)")
        )
      }
      
    } else {
      cat("  No Entrez IDs found — skipping GO/KEGG\n")
    }
  } else {
    cat("  Too few genes for GO/KEGG — skipping\n")
  }
  
  # ----------------------------------------------------------
  # Step 10g: Protein interaction network
  # Maps gene symbols to STRING IDs, retrieves interactions
  # at medium confidence (score >= 400), builds undirected
  # graph, removes isolated nodes, calculates node metrics.
  # Network objects stored in a named list for plotting
  # after global maximum degree is determined.
  # ----------------------------------------------------------
  cat("Generating network plot...\n")
  
  genes_for_network <- unique(dt$gene_name_atac)
  genes_for_network <- genes_for_network[
    !is.na(genes_for_network) &
      genes_for_network != ""]
  
  if (length(genes_for_network) >= 3) {
    
    gene_df <- data.frame(gene = genes_for_network,
                          stringsAsFactors = FALSE)
    mapped <- tryCatch({
      string_db$map(gene_df, "gene",
                    removeUnmappedRows = TRUE)
    }, error = function(e) NULL)
    
    if (!is.null(mapped) && nrow(mapped) > 0) {
      
      interactions <- tryCatch({
        string_db$get_interactions(mapped$STRING_id)
      }, error = function(e) NULL)
      
      if (!is.null(interactions) &&
          nrow(interactions) > 0) {
        
        g <- graph_from_data_frame(
          d        = interactions[, c("from", "to")],
          directed = FALSE,
          vertices = mapped$STRING_id
        )
        
        id_to_gene <- setNames(mapped$gene,
                               mapped$STRING_id)
        V(g)$name  <- id_to_gene[V(g)$name]
        g          <- delete_vertices(
          g, which(igraph::degree(g) == 0))
        
        if (vcount(g) > 0) {
          
          node_metrics <- data.frame(
            gene_name   = V(g)$name,
            degree      = as.numeric(igraph::degree(g)),
            betweenness = as.numeric(
              igraph::betweenness(g)),
            stringsAsFactors = FALSE
          )
          
          cat("  Network nodes:", vcount(g), "\n")
          cat("  Network edges:", ecount(g), "\n")
          
          V(g)$degree      <- node_metrics$degree
          V(g)$betweenness <- node_metrics$betweenness
          
          # Store in named list for Block 11
          network_list[[comp$name]] <- list(
            graph   = g,
            metrics = node_metrics
          )
        }
      }
    }
  }
  
  rm(dt, label_dt); gc()
}
# 11. NETWORK PLOTS WITH FIXED GLOBAL SCALE ----
# Find global maximum degree across all comparisons.
# Use as fixed maximum for scale_size_continuous so
# node sizes are directly comparable and the same
# legend applies to all network plots.
cat("\n========================================\n")
cat("Generating network plots with fixed scale...\n")
cat("========================================\n")

# Find global maximum degree from stored networks
all_degrees <- unlist(lapply(network_list,
                             function(x) x$metrics$degree))
max_degree_global <- if (length(all_degrees) > 0)
  max(all_degrees) else 10

cat("Global maximum degree:", max_degree_global, "\n")

# Generate one network plot per comparison
for (comp in comparisons) {
  
  if (is.null(network_list[[comp$name]])) {
    cat("No network for:", comp$name, "— skipping\n")
    next
  }
  
  g_obj       <- network_list[[comp$name]]$graph
  metrics_obj <- network_list[[comp$name]]$metrics
  
  # Hub score = degree normalised by global maximum
  # Fixed 0-1 scale ensures same legend across all plots
  metrics_obj$hub_score <- metrics_obj$degree /
    max_degree_global
  
  V(g_obj)$hub_score <- metrics_obj$hub_score
  
  tg <- as_tbl_graph(g_obj)
  
  p_network <- ggraph(tg, layout = "fr") +
    geom_edge_link(alpha = 0.4, color = "gray60",
                   width = 0.8) +
    geom_node_point(
      aes(size = degree, color = hub_score),
      alpha = 0.9
    ) +
    geom_node_text(
      aes(label = name),
      repel        = TRUE,
      fontface     = "bold",
      size         = 3,
      max.overlaps = 30
    ) +
    scale_color_gradient(
      low    = "#4575B4",
      high   = "#D73027",
      name   = "Hub Score\n(normalised\ndegree)",
      limits = c(0, 1)
    ) +
    scale_size_continuous(
      range  = c(3, 12),
      limits = c(1, max_degree_global),
      name   = "Connections"
    ) +
    labs(
      title    = paste("Protein Interaction Network:",
                       comp$display),
      subtitle = paste("Node size = connections |",
                       "Node colour = hub importance |",
                       "STRING score >= 400")
    ) +
    theme_graph(base_family = "sans") +
    theme(
      plot.title      = element_text(hjust = 0.5,
                                     face  = "bold",
                                     size  = 12),
      plot.subtitle   = element_text(hjust = 0.5,
                                     size  = 8,
                                     color = "gray40"),
      plot.background = element_rect(fill  = "white",
                                     color = NA),
      legend.position = "right"
    )
  
  save_plot(
    p_network,
    here(out_base, "network",
         paste0("network_", comp$name)),
    width = 12, height = 10
  )
  
  write.xlsx(
    metrics_obj[order(-metrics_obj$hub_score), ],
    here(out_base, "network",
         paste0("hub_genes_", comp$name, ".xlsx"))
  )
  
  cat("Saved network:", comp$name, "\n")
}
# 12. 3-WAY VENN DIAGRAMS ----
cat("\n========================================\n")
cat("Generating 3-way Venn diagrams...\n")
cat("========================================\n")

# Group 1 — antagomirs vs scramble
g1_genes <- list(
  overlap_genes_by_comp[["AM133_vs_AMScr"]],
  overlap_genes_by_comp[["AMAll_vs_AMScr"]],
  overlap_genes_by_comp[["AM1_206_vs_AMScr"]]
)

if (sum(sapply(g1_genes, length)) > 0) {
  draw_venn_3way(
    genes1    = g1_genes[[1]],
    genes2    = g1_genes[[2]],
    genes3    = g1_genes[[3]],
    labels    = c("AM133\nvs AMScr",
                  "AMAll\nvs AMScr",
                  "AM1/206\nvs AMScr"),
    title     = "Group 1: Antagomirs vs Scramble\nOverlapping Genes",
    base_path = here(out_base, "venn",
                     "venn_3way_Group1")
  )
  cat("Saved Group 1 3-way Venn\n")
}

# Group 2 — antagomirs vs each other
g2_genes <- list(
  overlap_genes_by_comp[["AM133_vs_AM1_206"]],
  overlap_genes_by_comp[["AM133_vs_AMAll"]],
  overlap_genes_by_comp[["AM1_206_vs_AMAll"]]
)

if (sum(sapply(g2_genes, length)) > 0) {
  draw_venn_3way(
    genes1    = g2_genes[[1]],
    genes2    = g2_genes[[2]],
    genes3    = g2_genes[[3]],
    labels    = c("AM133\nvs AM1/206",
                  "AM133\nvs AMAll",
                  "AM1/206\nvs AMAll"),
    title     = "Group 2: Antagomirs vs Each Other\nOverlapping Genes",
    base_path = here(out_base, "venn",
                     "venn_3way_Group2")
  )
  cat("Saved Group 2 3-way Venn\n")
}

# 13. REPRODUCIBILITY RECORDS ----
cat("\nSaving reproducibility records...\n")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

session_info <- capture.output(sessionInfo())
session_file <- here(
  "reproducibility",
  "Script4_integration_visualisation",
  "sessionInfo",
  paste0("sessionInfo_Script4_", timestamp, ".txt")
)
writeLines(session_info, session_file)
cat("Session info saved:", basename(session_file), "\n")

input_files <- c(
  sapply(comparisons, function(x) x$integration),
  sapply(comparisons, function(x) x$annotated),
  here(integration_dir, "summary",
       "integration_summary.tsv")
)
output_files   <- list.files(out_base,
                             recursive  = TRUE,
                             full.names = TRUE)
all_files      <- c(input_files, output_files)
existing_files <- all_files[file.exists(all_files)]

md5_records <- data.table(
  file   = basename(existing_files),
  path   = existing_files,
  md5sum = sapply(existing_files, tools::md5sum,
                  USE.NAMES = FALSE),
  type   = ifelse(existing_files %in% input_files,
                  "input", "output"),
  script = "Script4_integration_visualisation",
  run_at = format(Sys.time())
)
md5_file <- here(
  "reproducibility",
  "Script4_integration_visualisation",
  "md5",
  paste0("md5_Script4_", timestamp, ".csv")
)
fwrite(md5_records, md5_file)
cat("MD5 checksums saved:", basename(md5_file), "\n")

script_end <- Sys.time()
cat("\nScript completed:", format(script_end), "\n")
cat("Total runtime:",
    format(round(script_end - script_start, 2)), "\n")

save.image(here("Script4_integration_visualisation.RData"))
cat("Workspace saved\n")
cat("Script 4 complete!\n")
