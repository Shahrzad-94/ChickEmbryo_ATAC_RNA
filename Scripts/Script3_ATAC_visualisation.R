# ============================================================
# Script: Script3_ATAC_visualisation.R
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
#   Generates all ATAC-seq visualisations from annotated
#   peak files produced by Script 1. Plots include:
#   genomic feature pie charts, combined annotation bar
#   charts, distance-to-TSS bidirectional bar charts,
#   Manhattan plots with Excel output, and circos plots.
#   No re-annotation is performed — all annotation columns
#   are read directly from Script 1 output files.
#
# Input:
#   results/Script1_annotation/annotated/[name]_annotated.tsv
#
# Output:
#   results/Script3_ATAC_visualisation/
#     genomic_features/pie/
#     genomic_features/annobar/
#     genomic_features/distance_to_tss/
#     manhattan/
#     circos/
#   reproducibility/Script3_ATAC_visualisation/
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
library(ggpubr)
library(openxlsx)
library(dplyr)

cat("Working directory:", here(), "\n")

# Verify package versions match those used when script
# was written — warns but does not stop if mismatch found
expected_versions <- list(
  "here"       = "1.0.2",
  "data.table" = "1.18.4",
  "ggplot2"    = "4.0.3",
  "ggrepel"    = "0.9.8",
  "ggpubr"     = "0.6.3",
  "openxlsx"   = "4.2.8.1",
  "dplyr"      = "1.2.1"
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
# All colours defined once here and referenced by name
# throughout the script — guarantees consistency across
# all plots of the same type.
# All colours are from the Okabe-Ito colourblind-safe palette.

# Genomic feature categories — same colour for same feature
# across ALL pie charts and ALL bar charts
feature_colours <- c(
  "Promoter"          = "#0072B2",
  "5' UTR"            = "#56B4E9",
  "3' UTR"            = "#009E73",
  "Exon"              = "#E69F00",
  "Intron"            = "#CC79A7",
  "Distal Intergenic" = "#D55E00",
  "Downstream"        = "#F0E442"
)

# Manhattan plot colours — direction of chromatin change
manhattan_colours <- c(
  "Opening" = "#D55E00",
  "Closing" = "#0072B2"
)

# Circos plot colours — one colour per comparison
# Same colours used across all three circos versions
circos_colours <- c(
  "AM133 vs AMScr"   = "#D55E00",
  "AMAll vs AMScr"   = "#0072B2",
  "AM1/206 vs AMScr" = "#009E73",
  "AM1/206 vs AM133" = "#E69F00",
  "AMAll vs AM133"   = "#CC79A7",
  "AMAll vs AM1/206" = "#56B4E9"
)

# 3. DEFINE PARAMETERS AND FILE LIST ----
# ATAC significance threshold
atac_padj_threshold <- 0.05

# Six comparisons — each references its Script 1 output file
comparisons <- list(
  list(name      = "AM133_vs_AMScr",
       display   = "AM133 vs AMScr",
       group     = "Group1_vs_Scramble",
       annotated = here("results", "Script1_annotation",
                        "annotated",
                        "AM133_vs_AMScr_annotated.tsv")),
  list(name      = "AMAll_vs_AMScr",
       display   = "AMAll vs AMScr",
       group     = "Group1_vs_Scramble",
       annotated = here("results", "Script1_annotation",
                        "annotated",
                        "AMAll_vs_AMScr_annotated.tsv")),
  list(name      = "AM1_206_vs_AMScr",
       display   = "AM1/206 vs AMScr",
       group     = "Group1_vs_Scramble",
       annotated = here("results", "Script1_annotation",
                        "annotated",
                        "AM1_206_vs_AMScr_annotated.tsv")),
  list(name      = "AM133_vs_AM1_206",
       display   = "AM133 vs AM1/206",
       group     = "Group2_vs_Each_Other",
       annotated = here("results", "Script1_annotation",
                        "annotated",
                        "AM133_vs_AM1_206_annotated.tsv")),
  list(name      = "AM133_vs_AMAll",
       display   = "AM133 vs AMAll",
       group     = "Group2_vs_Each_Other",
       annotated = here("results", "Script1_annotation",
                        "annotated",
                        "AM133_vs_AMAll_annotated.tsv")),
  list(name      = "AM1_206_vs_AMAll",
       display   = "AM1/206 vs AMAll",
       group     = "Group2_vs_Each_Other",
       annotated = here("results", "Script1_annotation",
                        "annotated",
                        "AM1_206_vs_AMAll_annotated.tsv"))
)

# 4. CREATE OUTPUT FOLDER STRUCTURE ----
out_base <- here("results", "Script3_ATAC_visualisation")
dir.create(here(out_base, "genomic_features", "pie"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(here(out_base, "genomic_features", "annobar"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(here(out_base, "genomic_features",
                "distance_to_tss"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(here(out_base, "manhattan"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(here(out_base, "circos"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(here("reproducibility",
                "Script3_ATAC_visualisation",
                "sessionInfo"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(here("reproducibility",
                "Script3_ATAC_visualisation",
                "md5"),
           recursive = TRUE, showWarnings = FALSE)
cat("Output folders created\n")

# 5. HELPER FUNCTIONS ----

# ----------------------------------------------------------
# simplify_annotation(): converts detailed ChIPseeker
# annotation labels to clean category names for plotting.
# Uses grepl pattern matching on the annotation string.
# Applied consistently to all files so categories are
# identical across all plots.
# ----------------------------------------------------------
simplify_annotation <- function(annotation) {
  dplyr::case_when(
    grepl("^Promoter",          annotation) ~ "Promoter",
    grepl("^Exon",              annotation) ~ "Exon",
    grepl("^Intron",            annotation) ~ "Intron",
    grepl("^Downstream",        annotation) ~ "Downstream",
    grepl("^Distal Intergenic", annotation) ~ "Distal Intergenic",
    grepl("^5' UTR",            annotation) ~ "5' UTR",
    grepl("^3' UTR",            annotation) ~ "3' UTR",
    TRUE                                    ~ "Other"
  )
}

# ----------------------------------------------------------
# save_plot(): saves a ggplot object as both PNG and PDF.
# All plots saved at 300 DPI for publication quality.
# limitsize = FALSE allows tall plots to render correctly.
# ----------------------------------------------------------
save_plot <- function(p, base_path, width, height) {
  ggsave(paste0(base_path, ".png"), p,
         width     = width,
         height    = height,
         dpi       = 300,
         bg        = "white",
         limitsize = FALSE)
  ggsave(paste0(base_path, ".pdf"), p,
         width     = width,
         height    = height,
         bg        = "white",
         limitsize = FALSE)
  cat("Saved:", basename(base_path), ".png and .pdf\n")
}

# 6. GENOMIC FEATURE PIE CHARTS ----
# One pie chart per comparison.
# Reads the annotation column from Script 1 output directly.
# Filters to significant peaks (padj <= threshold).
# Simplifies annotation labels to clean categories.
# All feature categories are included in every legend
# regardless of whether that category is present in a
# given comparison, ensuring a consistent legend
# across all six charts. All charts use feature_colours.
cat("\n========================================\n")
cat("Generating genomic feature pie charts...\n")
cat("========================================\n")

for (comp in comparisons) {
  
  cat("Pie chart:", comp$name, "\n")
  
  # Read annotated file from Script 1
  dt <- fread(comp$annotated, sep = "\t")
  
  # Filter to significant peaks only
  sig <- dt[padj <= atac_padj_threshold]
  
  if (nrow(sig) == 0) {
    cat("No significant peaks — skipping\n")
    rm(dt, sig); gc(); next
  }
  
  # Simplify annotation labels
  sig[, feature := simplify_annotation(annotation)]
  
  # Count peaks per feature category and calculate percentage
  feat_counts <- sig[, .N, by = feature]
  feat_counts[, pct := N / sum(N) * 100]
  feat_counts[, label := paste0(feature, "\n",
                                round(pct, 1), "%")]
  
  # Enforce all feature categories in the legend,
  # including those absent from this comparison,
  # so the legend is identical across all pie charts
  all_levels <- names(feature_colours)
  present     <- feat_counts$feature
  missing_lvls <- setdiff(all_levels, present)
  
  if (length(missing_lvls) > 0) {
    missing_dt <- data.table(
      feature = missing_lvls,
      N       = 0L,
      pct     = 0
    )
    feat_counts <- rbindlist(list(feat_counts, missing_dt),
                             fill = TRUE)
  }
  
  feat_counts[, feature := factor(feature,
                                  levels = all_levels)]
  feat_counts[, label := paste0(feature, "\n",
                                round(pct, 1), "%")]
  
  p_pie <- ggplot(feat_counts,
                  aes(x    = "",
                      y    = pct,
                      fill = feature)) +
    geom_bar(stat     = "identity",
             width    = 1,
             color    = "white",
             linewidth = 0.3) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = feature_colours,
                      name   = "Genomic Feature",
                      drop   = FALSE) +
    labs(title = paste("Genomic Feature Distribution:",
                       comp$display),
         subtitle = paste0("Significant peaks (padj \u2264 ",
                           atac_padj_threshold, ")  |  n = ",
                           nrow(sig))) +
    theme_void() +
    theme(
      plot.title    = element_text(hjust  = 0.5,
                                   face   = "bold",
                                   size   = 13),
      plot.subtitle = element_text(hjust  = 0.5,
                                   size   = 9,
                                   color  = "gray40"),
      legend.position = "right",
      legend.title    = element_text(face = "bold",
                                     size = 10),
      legend.text     = element_text(size = 9),
      plot.margin     = margin(10, 10, 10, 10)
    )
  
  save_plot(
    p_pie,
    here(out_base, "genomic_features", "pie",
         paste0("pie_", comp$name)),
    width = 8, height = 6
  )
  
  rm(dt, sig, feat_counts, p_pie); gc()
}

# 7. COMBINED ANNOTATION BAR CHART (ANNOBAR) ----
# Three versions: all 6 comparisons, Group 1 only,
# Group 2 only. Each chart shows one 100% stacked bar
# per comparison, with comparisons on the X axis.
# All bars use feature_colours — same colours as pie
# charts for direct comparability. Annotation categories
# simplified identically to pie charts.
cat("\n========================================\n")
cat("Generating combined annotation bar charts...\n")
cat("========================================\n")

# Build combined feature table across all comparisons
all_feat_list <- list()

for (comp in comparisons) {
  
  dt  <- fread(comp$annotated, sep = "\t")
  sig <- dt[padj <= atac_padj_threshold]
  
  if (nrow(sig) == 0) {
    rm(dt, sig); gc(); next
  }
  
  sig[, feature := simplify_annotation(annotation)]
  feat_counts   <- sig[, .N, by = feature]
  feat_counts[, pct        := N / sum(N) * 100]
  feat_counts[, comparison := comp$name]
  feat_counts[, display    := comp$display]
  feat_counts[, group      := comp$group]
  
  all_feat_list[[comp$name]] <- feat_counts
  rm(dt, sig, feat_counts); gc()
}

all_feat_dt <- rbindlist(all_feat_list, fill = TRUE)

# Set factor levels for consistent ordering
all_feat_dt[, feature := factor(
  feature, levels = rev(names(feature_colours)))]
all_feat_dt[, display := factor(
  display,
  levels = sapply(comparisons, function(x) x$display))]

# ----------------------------------------------------------
# make_annobar(): builds one 100% stacked bar chart
# from a subset of the combined feature table.
# Each bar represents one comparison.
# Vertical bars with angled X axis labels matching
# standard annotation bar chart orientation in the field.
# ----------------------------------------------------------
make_annobar <- function(feat_dt, title, subtitle) {
  
  p <- ggplot(feat_dt,
              aes(x    = display,
                  y    = pct,
                  fill = feature)) +
    geom_bar(stat      = "identity",
             position  = "stack",
             color     = "white",
             linewidth = 0.3,
             width     = 0.4) +
    scale_fill_manual(values = feature_colours,
                      name   = "Genomic Feature",
                      drop   = FALSE) +
    scale_y_continuous(labels = function(x)
      paste0(x, "%"),
      expand = c(0, 0),
      limits = c(0, 101)) +
    labs(title    = title,
         subtitle = subtitle,
         x        = NULL,
         y        = "Percentage of peaks") +
    theme_bw(base_size = 12) +
    theme(
      plot.title         = element_text(hjust = 0.5,
                                        face  = "bold",
                                        size  = 13),
      plot.subtitle      = element_text(hjust = 0.5,
                                        size  = 9,
                                        color = "gray40"),
      axis.text.x        = element_text(angle = 35,
                                        hjust = 1,
                                        face  = "bold",
                                        size  = 10),
      legend.position    = "right",
      legend.title       = element_text(face  = "bold",
                                        size  = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )
  return(p)
}

# Version 1: All 6 comparisons
p_bar_all <- make_annobar(
  all_feat_dt,
  title    = "Genomic Feature Distribution — All Comparisons",
  subtitle = paste0("Significant peaks (padj \u2264 ",
                    atac_padj_threshold, ")")
)
save_plot(p_bar_all,
          here(out_base, "genomic_features", "annobar",
               "annobar_all_comparisons"),
          width = 10, height = 6)

# Version 2: Group 1 only (vs Scramble)
p_bar_g1 <- make_annobar(
  all_feat_dt[group == "Group1_vs_Scramble"],
  title    = "Genomic Feature Distribution — Group 1 (vs Scramble)",
  subtitle = paste0("Significant peaks (padj \u2264 ",
                    atac_padj_threshold, ")")
)
save_plot(p_bar_g1,
          here(out_base, "genomic_features", "annobar",
               "annobar_Group1_vs_Scramble"),
          width = 5, height = 6)

# Version 3: Group 2 only (vs Each Other)
p_bar_g2 <- make_annobar(
  all_feat_dt[group == "Group2_vs_Each_Other"],
  title    = "Genomic Feature Distribution — Group 2 (vs Each Other)",
  subtitle = paste0("Significant peaks (padj \u2264 ",
                    atac_padj_threshold, ")")
)
save_plot(p_bar_g2,
          here(out_base, "genomic_features", "annobar",
               "annobar_Group2_vs_Each_Other"),
          width = 5, height = 6)
# 8. DISTANCE TO TSS BIDIRECTIONAL BAR CHART ----
# One chart per comparison.
# Uses the distanceToTSS column from Script 1 output directly.
# Bins peaks into 12 distance bands (in kb) matching the
# reference style: upstream bands shown as negative values,
# downstream bands as positive values.
# Percentage of total significant peaks per band is shown.
# Y axis is fixed at 40% across all comparisons for
# direct visual comparability.
# Upstream bands coloured Blue, downstream with Vermillion
# to visually emphasise the bidirectional nature.
cat("\n========================================\n")
cat("Generating distance to TSS bar charts...\n")
cat("========================================\n")

# Define 12 distance bands in kb
# Negative = upstream of TSS, positive = downstream
tss_breaks <- c(-Inf, -100, -10, -5, -3, -1,
                0, 1, 3, 5, 10, 100, Inf)
tss_labels <- c("< -100", "-100 to -10", "-10 to -5",
                "-5 to -3", "-3 to -1", "-1 to 0",
                "0 to 1", "1 to 3", "3 to 5",
                "5 to 10", "10 to 100", "> 100")

# Upstream bands (negative distance) coloured Blue
# Downstream bands (positive distance) coloured Vermillion
tss_colours <- c(
  "< -100"       = "#0072B2",
  "-100 to -10"  = "#0072B2",
  "-10 to -5"    = "#0072B2",
  "-5 to -3"     = "#0072B2",
  "-3 to -1"     = "#0072B2",
  "-1 to 0"      = "#0072B2",
  "0 to 1"       = "#D55E00",
  "1 to 3"       = "#D55E00",
  "3 to 5"       = "#D55E00",
  "5 to 10"      = "#D55E00",
  "10 to 100"    = "#D55E00",
  "> 100"        = "#D55E00"
)

for (comp in comparisons) {
  
  cat("Distance to TSS:", comp$name, "\n")
  
  dt  <- fread(comp$annotated, sep = "\t")
  sig <- dt[padj <= atac_padj_threshold]
  
  if (nrow(sig) == 0) {
    cat("No significant peaks — skipping\n")
    rm(dt, sig); gc(); next
  }
  
  # Convert distanceToTSS from bp to kb
  sig[, dist_kb := distanceToTSS / 1000]
  
  # Bin into 12 distance bands
  sig[, tss_band := cut(dist_kb,
                        breaks = tss_breaks,
                        labels = tss_labels,
                        right  = TRUE,
                        include.lowest = TRUE)]
  
  # Count peaks per band and calculate percentage
  band_counts <- sig[, .N, by = tss_band]
  band_counts[, pct := N / sum(N) * 100]
  
  # Ensure all bands present even if count is zero
  all_bands <- data.table(
    tss_band = factor(tss_labels, levels = tss_labels))
  band_counts <- merge(all_bands, band_counts,
                       by = "tss_band", all.x = TRUE)
  band_counts[is.na(pct), pct := 0]
  band_counts[is.na(N),   N   := 0]
  
  p_tss <- ggplot(band_counts,
                  aes(x    = tss_band,
                      y    = pct,
                      fill = tss_band)) +
    geom_bar(stat  = "identity",
             color = "white",
             linewidth = 0.3) +
    scale_fill_manual(values = tss_colours) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0),
      limits = c(0, 40)
    ) +
    labs(
      title    = paste("Peak Distance to TSS:",
                       comp$display),
      subtitle = paste0("Significant peaks (padj \u2264 ",
                        atac_padj_threshold,
                        ")  |  n = ", nrow(sig)),
      x        = "Distance from TSS (kb)",
      y        = "Binding sites (%)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title         = element_text(hjust = 0.5,
                                        face  = "bold",
                                        size  = 13),
      plot.subtitle      = element_text(hjust = 0.5,
                                        size  = 9,
                                        color = "gray40"),
      axis.text.x        = element_text(angle = 45,
                                        hjust = 1,
                                        size  = 9),
      legend.position    = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )
  
  save_plot(
    p_tss,
    here(out_base, "genomic_features", "distance_to_tss",
         paste0("distToTSS_", comp$name)),
    width = 9, height = 6
  )
  
  rm(dt, sig, band_counts, p_tss); gc()
}

# 9. MANHATTAN PLOTS ----
# One Manhattan plot per comparison plus one combined 2x3 panel.
# All peaks shown — significant peaks (padj <= threshold)
# coloured by chromatin change direction.
# Non-significant peaks shown as small grey points.
# Y axis = -log10(padj), no cap — maximum across all files
# is 13.6 which is well within a readable range.
# Combined panel uses fixed Y axis limit = 15 so all
# 6 comparisons are directly comparable.
# Excel file saved per comparison with all peaks and a
# significant (TRUE/FALSE) column for dot identification.
cat("\n========================================\n")
cat("Generating Manhattan plots...\n")
cat("========================================\n")

# GalGal6 chromosome sizes in base pairs
# Source: UCSC Genome Browser GalGal6 assembly
chrom_sizes <- data.table(
  chrom = paste0("chr", c(1:28, "Z", "W")),
  size  = c(
    197608386, 150435657, 112617285, 92735134,
    59543457,  35381781,  37306461,  30030098,
    24191541,  20016774,  19807102,  21459726,
    18220052,  16022694,  13786070,  13269384,
    11017023,  10677681,  9955776,   13750054,
    6834642,   4664182,   6150831,   6585680,
    2041316,   5011791,   3742660,   4519568,
    73182983,  3030853
  )
)

# Build cumulative genome positions for X axis placement
# Chromosomes ordered chr1 to chr28, chrZ, chrW
chrom_sizes[, chrom_num := as.numeric(
  gsub("chrZ", "29", gsub("chrW", "30",
                          gsub("chr", "", chrom))))]
chrom_sizes <- chrom_sizes[order(chrom_num)]
chrom_sizes[, cum_start := c(0,
                             cumsum(size)[-.N])]
chrom_sizes[, cum_end   := cumsum(size)]
chrom_sizes[, mid       := cum_start + size / 2]

total_genome_size <- sum(chrom_sizes$size)

# Alternating grey backgrounds for chromosome bands
chrom_sizes[, bg    := ifelse(
  seq_len(.N) %% 2 == 0, "#F0F0F0", "#FAFAFA")]
chrom_sizes[, label := gsub("chr", "", chrom)]

# ----------------------------------------------------------
# make_manhattan(): builds one Manhattan plot from an
# annotated ATAC-seq file. Calculates cumulative genome
# position for each peak midpoint. Returns the ggplot
# object so it can be used in the combined panel.
# fixed_ymax: if provided, fixes Y axis limit for panel
# ----------------------------------------------------------
make_manhattan <- function(comp, fixed_ymax = NULL) {
  
  cat("Manhattan:", comp$name, "\n")
  
  dt <- fread(comp$annotated, sep = "\t")
  dt <- dt[Chr %in% chrom_sizes$chrom]
  dt <- dt[!is.na(padj) & padj > 0]
  
  # Calculate cumulative genome position for each peak
  dt <- merge(dt,
              chrom_sizes[, .(chrom, cum_start)],
              by.x  = "Chr",
              by.y  = "chrom",
              all.x = TRUE)
  dt[, peak_mid  := (Start + End) / 2]
  dt[, cum_pos   := cum_start + peak_mid]
  dt[, neg_log10 := -log10(padj)]
  dt[, direction := fifelse(log2FoldChange > 0,
                            "Opening", "Closing")]
  dt[, significant := padj <= atac_padj_threshold]
  
  y_max    <- if (!is.null(fixed_ymax)) fixed_ymax else
    max(dt$neg_log10, na.rm = TRUE) + 0.5
  sig_line <- -log10(atac_padj_threshold)
  
  cat("  Total peaks:", nrow(dt), "\n")
  cat("  Significant:", sum(dt$significant), "\n")
  cat("  Max -log10(padj):",
      round(max(dt$neg_log10, na.rm = TRUE), 2), "\n")
  
  # Save Excel file with all peaks and significant column
  # Allows identification of individual dots on the plot
  excel_out <- here(out_base, "manhattan",
                    paste0("manhattan_peaks_",
                           comp$name, ".xlsx"))
  write.xlsx(
    dt[, .(Chr, Start, End, log2FoldChange, padj,
           neg_log10_padj = neg_log10,
           direction, significant,
           ensembl_gene_id, external_gene_name,
           annotation, distanceToTSS)],
    excel_out
  )
  cat("  Saved Excel:", basename(excel_out), "\n")
  
  p <- ggplot() +
    
    # Alternating chromosome background rectangles
    geom_rect(
      data  = chrom_sizes,
      aes(xmin = cum_start, xmax = cum_end,
          ymin = 0,         ymax = y_max),
      fill  = chrom_sizes$bg,
      color = NA
    ) +
    
    # Non-significant peaks — small grey dots
    geom_point(
      data  = dt[significant == FALSE],
      aes(x = cum_pos, y = neg_log10),
      color = "gray80", size = 0.15, alpha = 0.4
    ) +
    
    # Significant peaks — coloured by direction
    geom_point(
      data  = dt[significant == TRUE],
      aes(x = cum_pos, y = neg_log10,
          color = direction),
      size  = 0.8, alpha = 0.8
    ) +
    
    # Significance threshold line
    geom_hline(
      yintercept = sig_line,
      linetype   = "dashed",
      color      = "black",
      linewidth  = 0.5
    ) +
    
    # Label for threshold line — placed outside plot on right
    annotate(
      "text",
      x     = total_genome_size * 1.01,
      y     = sig_line,
      label = paste0("padj \u2264 ", atac_padj_threshold),
      size  = 3,
      color = "black",
      hjust = 0,
      vjust = 0.5
    ) +
    
    scale_color_manual(values = manhattan_colours,
                       name   = "Chromatin\nChange") +
    scale_x_continuous(
      breaks = chrom_sizes$mid,
      labels = chrom_sizes$label,
      expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
      expand = c(0.02, 0),
      limits = c(0, y_max),
      breaks = pretty(c(0, y_max), n = 5)
    ) +
    labs(
      title = paste("ATAC-seq Differential Peaks:",
                    comp$display),
      x     = "Chromosome",
      y     = expression(-log[10](padj))
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 11) +
    theme(
      plot.background    = element_rect(fill  = "white",
                                        color = NA),
      panel.background   = element_rect(fill  = "white",
                                        color = NA),
      plot.title         = element_text(hjust = 0.5,
                                        face  = "bold",
                                        size  = 12),
      axis.text.x        = element_text(angle = 45,
                                        hjust = 1,
                                        size  = 7),
      axis.text.y        = element_text(size  = 9),
      axis.title         = element_text(size  = 10,
                                        face  = "bold"),
      legend.position    = "right",
      panel.grid         = element_blank(),
      panel.border       = element_rect(color     = "gray50",
                                        fill      = NA,
                                        linewidth = 0.5),
      plot.margin        = margin(10, 60, 10, 10)
    )
  
  # Save individual Manhattan plot
  save_plot(
    p,
    here(out_base, "manhattan",
         paste0("manhattan_", comp$name)),
    width = 18, height = 6
  )
  
  rm(dt); gc()
  return(p)
}

# Generate individual Manhattan plots
# Fixed Y axis at 15 for all (max across all files = 13.6)
manhattan_plots <- list()
for (comp in comparisons) {
  manhattan_plots[[comp$name]] <- make_manhattan(
    comp, fixed_ymax = 15
  )
}

# Combined 2x3 panel — fixed Y axis ensures direct
# comparability across all 6 comparisons
cat("Generating combined Manhattan panel...\n")

combined_manhattan <- ggarrange(
  plotlist      = manhattan_plots,
  nrow          = 2,
  ncol          = 3,
  common.legend = TRUE,
  legend        = "right",
  labels        = c("A", "B", "C", "D", "E", "F"),
  font.label    = list(size = 11, face = "bold")
)
combined_manhattan <- annotate_figure(
  combined_manhattan,
  top = text_grob(
    "Genome-wide ATAC-seq Differential Peak Distribution",
    face = "bold", size = 13
  )
)
ggsave(
  here(out_base, "manhattan",
       "manhattan_combined_all.png"),
  combined_manhattan,
  width = 22, height = 10, dpi = 300, bg = "white"
)
ggsave(
  here(out_base, "manhattan",
       "manhattan_combined_all.pdf"),
  combined_manhattan,
  width = 22, height = 10, bg = "white"
)
cat("Saved combined Manhattan panel\n")

# 10. CIRCOS PLOTS ----
# Three circular peak density plots:
#   - All 6 comparisons
#   - Group 1 only (vs Scramble)
#   - Group 2 only (vs Each Other)
# Peak density calculated per 1 Mb window per chromosome.
# Density normalised 0-1 within each chromosome so all
# chromosomes are visually comparable regardless of size.
# Uses circos_colours — Okabe-Ito palette as defined above.
cat("\n========================================\n")
cat("Generating circos plots...\n")
cat("========================================\n")

# Circos comparison file mapping — display name to file path
circos_file_map <- list(
  "AM133 vs AMScr"   = here("results", "Script1_annotation",
                            "annotated",
                            "AM133_vs_AMScr_annotated.tsv"),
  "AMAll vs AMScr"   = here("results", "Script1_annotation",
                            "annotated",
                            "AMAll_vs_AMScr_annotated.tsv"),
  "AM1/206 vs AMScr" = here("results", "Script1_annotation",
                            "annotated",
                            "AM1_206_vs_AMScr_annotated.tsv"),
  "AM1/206 vs AM133" = here("results", "Script1_annotation",
                            "annotated",
                            "AM133_vs_AM1_206_annotated.tsv"),
  "AMAll vs AM133"   = here("results", "Script1_annotation",
                            "annotated",
                            "AM133_vs_AMAll_annotated.tsv"),
  "AMAll vs AM1/206" = here("results", "Script1_annotation",
                            "annotated",
                            "AM1_206_vs_AMAll_annotated.tsv")
)

inner_radius <- 3.0
track_width  <- 0.7

# ----------------------------------------------------------
# build_density_tracks(): calculates normalised peak density
# per 1 Mb window per chromosome per comparison.
# Density normalised 0-1 within each chromosome so
# chromosomes of different sizes are visually comparable.
# Returns a data.table with genome position, normalised
# density, and track position for each window.
# ----------------------------------------------------------
build_density_tracks <- function(comp_names) {
  
  all_tracks <- list()
  
  for (i in seq_along(comp_names)) {
    comp_name <- comp_names[i]
    cat("Density:", comp_name, "\n")
    
    dt  <- fread(circos_file_map[[comp_name]], sep = "\t")
    dt  <- dt[padj <= atac_padj_threshold]
    dt  <- dt[Chr %in% chrom_sizes$chrom]
    
    base_r <- inner_radius + (i - 1) * track_width
    
    track_rows <- list()
    for (chr in chrom_sizes$chrom) {
      chr_row  <- chrom_sizes[chrom == chr]
      chr_data <- dt[Chr == chr]
      chr_size <- chr_row$size
      windows  <- seq(0, chr_size, by = 1000000)
      if (length(windows) < 2)
        windows <- c(0, chr_size)
      counts <- numeric(length(windows) - 1)
      if (nrow(chr_data) > 0) {
        mids <- (chr_data$Start + chr_data$End) / 2
        for (j in seq_along(counts)) {
          counts[j] <- sum(mids >= windows[j] &
                             mids <  windows[j + 1])
        }
      }
      max_c      <- max(counts, 1)
      norm       <- counts / max_c
      win_mids   <- (windows[-length(windows)] +
                       windows[-1]) / 2
      genome_pos <- chr_row$cum_start + win_mids
      
      track_rows[[chr]] <- data.table(
        genome_pos = genome_pos,
        norm       = norm,
        comparison = comp_name,
        base       = base_r,
        top        = base_r + norm * track_width
      )
    }
    all_tracks[[comp_name]] <- rbindlist(track_rows)
    rm(dt); gc()
  }
  rbindlist(all_tracks)
}

# ----------------------------------------------------------
# make_circos(): builds one ggplot2 polar circos plot.
# comp_names: vector of comparison display names to include.
# colors:     named colour vector for fill.
# title:      plot title.
# filename:   output file base name.
# ----------------------------------------------------------
make_circos <- function(comp_names, colors,
                        title, filename) {
  
  cat("Circos:", title, "\n")
  
  tracks   <- build_density_tracks(comp_names)
  n_tracks <- length(comp_names)
  outer_r  <- inner_radius + n_tracks * track_width
  label_r  <- outer_r + 0.4
  
  p <- ggplot() +
    
    # Chromosome background segments
    geom_rect(
      data      = chrom_sizes,
      aes(xmin  = cum_start,
          xmax  = cum_end,
          ymin  = inner_radius - 0.15,
          ymax  = outer_r + 0.15),
      fill      = "gray94",
      color     = "white",
      linewidth = 0.5
    ) +
    
    # Peak density bars per window per comparison
    geom_rect(
      data      = tracks[norm > 0],
      aes(xmin  = genome_pos - 400000,
          xmax  = genome_pos + 400000,
          ymin  = base,
          ymax  = top,
          fill  = comparison),
      alpha     = 0.9
    ) +
    
    scale_fill_manual(values = colors,
                      name   = "Comparison",
                      breaks = comp_names) +
    
    # Chromosome labels outside outermost track
    geom_text(
      data      = chrom_sizes,
      aes(x     = mid,
          y     = label_r,
          label = label),
      size      = 2.2,
      fontface  = "bold",
      color     = "gray20"
    ) +
    
    coord_polar(theta = "x", start = 0, direction = 1) +
    scale_x_continuous(limits = c(0, total_genome_size),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, label_r + 0.8)) +
    labs(title = title) +
    theme_void() +
    theme(
      plot.background = element_rect(fill  = "white",
                                     color = NA),
      plot.title      = element_text(hjust  = 0.5,
                                     face   = "bold",
                                     size   = 15,
                                     margin = margin(b = 10,
                                                     t = 15)),
      legend.position = c(0.88, 0.15),
      legend.title    = element_text(size = 9,
                                     face = "bold"),
      legend.text     = element_text(size = 8),
      legend.key.size = unit(0.4, "cm"),
      plot.margin     = margin(20, 20, 20, 20)
    )
  
  ggsave(here(out_base, "circos",
              paste0("circos_", filename, ".png")),
         p, width = 14, height = 14,
         dpi = 300, bg = "white")
  ggsave(here(out_base, "circos",
              paste0("circos_", filename, ".pdf")),
         p, width = 14, height = 14, bg = "white")
  cat("Saved circos:", filename, "\n")
}

# All 6 comparisons
make_circos(
  comp_names = names(circos_file_map),
  colors     = circos_colours,
  title      = paste("Genome-wide ATAC-seq Peak",
                     "Distribution \u2014 All Comparisons"),
  filename   = "all_comparisons"
)

# Group 1 only
make_circos(
  comp_names = c("AM133 vs AMScr",
                 "AMAll vs AMScr",
                 "AM1/206 vs AMScr"),
  colors     = circos_colours,
  title      = "Group 1: Antagomirs vs Scramble Control",
  filename   = "Group1_vs_Scramble"
)

# Group 2 only
make_circos(
  comp_names = c("AM1/206 vs AM133",
                 "AMAll vs AM133",
                 "AMAll vs AM1/206"),
  colors     = circos_colours,
  title      = "Group 2: Antagomirs vs Each Other",
  filename   = "Group2_vs_Each_Other"
)

# 11. REPRODUCIBILITY RECORDS ----
cat("\nSaving reproducibility records...\n")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

session_info <- capture.output(sessionInfo())
session_file <- here(
  "reproducibility", "Script3_ATAC_visualisation",
  "sessionInfo",
  paste0("sessionInfo_Script3_", timestamp, ".txt")
)
writeLines(session_info, session_file)
cat("Session info saved:", basename(session_file), "\n")

input_files  <- sapply(comparisons,
                       function(x) x$annotated)
output_files <- list.files(out_base,
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
  script = "Script3_ATAC_visualisation",
  run_at = format(Sys.time())
)
md5_file <- here(
  "reproducibility", "Script3_ATAC_visualisation",
  "md5",
  paste0("md5_Script3_", timestamp, ".csv")
)
fwrite(md5_records, md5_file)
cat("MD5 checksums saved:", basename(md5_file), "\n")

script_end <- Sys.time()
cat("\nScript completed:", format(script_end), "\n")
cat("Total runtime:",
    format(round(script_end - script_start, 2)), "\n")

save.image(here("Script3_ATAC_visualisation.RData"))
cat("Workspace saved\n")
cat("Script 3 complete!\n")