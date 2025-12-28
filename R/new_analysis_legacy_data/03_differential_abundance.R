#!/usr/bin/env Rscript
# =============================================================================
# 03_differential_abundance.R
# Differential abundance analysis using ALDEx2 and MaAsLin2
# =============================================================================

library(tidyverse)
library(vegan)

# Ensure dplyr::select is used (MASS masks it)
select <- dplyr::select

# Check for required packages
if (!requireNamespace("ALDEx2", quietly = TRUE)) {
  cat("Installing ALDEx2...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("ALDEx2")
}

library(ALDEx2)

# Note: MaAsLin2 not available for R 4.5, skipping that analysis
# library(Maaslin2)

# Set paths
project_dir <- "/home/david/projects/abx_metagenomic_correlations"
results_dir <- file.path(project_dir, "results/new_analysis_legacy_data")

# Load prepared data
load(file.path(project_dir, "data/prepared_data.RData"))

sample_metadata <- prepared_data$sample_metadata
species_matrix <- prepared_data$species_matrix
genus_matrix <- prepared_data$genus_matrix
high_confidence_groups <- prepared_data$high_confidence_groups

# =============================================================================
# 1. Prepare Data for Analysis
# =============================================================================

cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("DIFFERENTIAL ABUNDANCE ANALYSIS\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")

# Use high-confidence samples only
sample_hc <- sample_metadata %>% filter(high_confidence)
hc_samples <- sample_hc$sample_id

# Filter matrices to high-confidence samples
species_hc <- species_matrix[hc_samples, ]
genus_hc <- genus_matrix[intersect(hc_samples, rownames(genus_matrix)), ]

# Filter to taxa present in at least 10% of samples
prevalence_threshold <- 0.10
min_samples <- ceiling(nrow(genus_hc) * prevalence_threshold)

genus_prevalence <- colSums(genus_hc > 0)
genus_filtered <- genus_hc[, genus_prevalence >= min_samples]

cat("High-confidence samples:", nrow(sample_hc), "\n")
cat("Genera after prevalence filter:", ncol(genus_filtered), "\n\n")

# =============================================================================
# 2. ALDEx2 Analysis - Any Antibiotic Exposure
# =============================================================================

cat("=== ALDEx2: Any Antibiotic Exposure (7-day) ===\n\n")

# Prepare data for ALDEx2 (taxa as rows, samples as columns)
aldex_counts <- t(genus_filtered)

# Create conditions vector
aldex_meta <- sample_hc[match(colnames(aldex_counts), sample_hc$sample_id), ]
conditions_any <- ifelse(aldex_meta$abx_any_7d, "exposed", "unexposed")

# Run ALDEx2
set.seed(42)
aldex_any <- aldex(
  reads = aldex_counts,
  conditions = conditions_any,
  mc.samples = 128,
  test = "t",
  effect = TRUE,
  include.sample.summary = FALSE,
  verbose = TRUE
)

# Get significant results
aldex_any_sig <- aldex_any %>%
  rownames_to_column("taxon") %>%
  filter(wi.eBH < 0.1) %>%
  arrange(wi.eBH)

cat("\nSignificant taxa (BH < 0.1):", nrow(aldex_any_sig), "\n")
if (nrow(aldex_any_sig) > 0) {
  cat("\nTop results:\n")
  print(aldex_any_sig %>%
          select(taxon, effect, wi.ep, wi.eBH) %>%
          head(20))
}

# Save full results
aldex_any_full <- aldex_any %>%
  rownames_to_column("taxon") %>%
  arrange(wi.eBH)

write_csv(aldex_any_full,
          file.path(results_dir, "differential_abundance/aldex2_any_abx_7d.csv"))
cat("\nSaved: differential_abundance/aldex2_any_abx_7d.csv\n")

# =============================================================================
# 3. ALDEx2 Analysis - Anti-Anaerobic Antibiotic Exposure
# =============================================================================

cat("\n=== ALDEx2: Anti-Anaerobic Antibiotic Exposure ===\n\n")

conditions_anaerobic <- ifelse(aldex_meta$abx_anaerobic_7d, "exposed", "unexposed")

# Only run if we have both exposed and unexposed
if (length(unique(conditions_anaerobic)) == 2) {

  set.seed(42)
  aldex_anaerobic <- aldex(
    reads = aldex_counts,
    conditions = conditions_anaerobic,
    mc.samples = 128,
    test = "t",
    effect = TRUE,
    include.sample.summary = FALSE,
    verbose = TRUE
  )

  aldex_anaerobic_sig <- aldex_anaerobic %>%
    rownames_to_column("taxon") %>%
    filter(wi.eBH < 0.1) %>%
    arrange(wi.eBH)

  cat("\nSignificant taxa (BH < 0.1):", nrow(aldex_anaerobic_sig), "\n")
  if (nrow(aldex_anaerobic_sig) > 0) {
    cat("\nTop results:\n")
    print(aldex_anaerobic_sig %>%
            select(taxon, effect, wi.ep, wi.eBH) %>%
            head(20))
  }

  aldex_anaerobic_full <- aldex_anaerobic %>%
    rownames_to_column("taxon") %>%
    arrange(wi.eBH)

  write_csv(aldex_anaerobic_full,
            file.path(results_dir, "differential_abundance/aldex2_anti_anaerobic.csv"))
  cat("\nSaved: differential_abundance/aldex2_anti_anaerobic.csv\n")
}

# =============================================================================
# 4. ALDEx2 Analysis - Broad-Spectrum Antibiotic Exposure
# =============================================================================

cat("\n=== ALDEx2: Broad-Spectrum Antibiotic Exposure ===\n\n")

conditions_broad <- ifelse(aldex_meta$abx_broad_7d, "exposed", "unexposed")

if (length(unique(conditions_broad)) == 2) {

  set.seed(42)
  aldex_broad <- aldex(
    reads = aldex_counts,
    conditions = conditions_broad,
    mc.samples = 128,
    test = "t",
    effect = TRUE,
    include.sample.summary = FALSE,
    verbose = TRUE
  )

  aldex_broad_sig <- aldex_broad %>%
    rownames_to_column("taxon") %>%
    filter(wi.eBH < 0.1) %>%
    arrange(wi.eBH)

  cat("\nSignificant taxa (BH < 0.1):", nrow(aldex_broad_sig), "\n")
  if (nrow(aldex_broad_sig) > 0) {
    cat("\nTop results:\n")
    print(aldex_broad_sig %>%
            select(taxon, effect, wi.ep, wi.eBH) %>%
            head(20))
  }

  aldex_broad_full <- aldex_broad %>%
    rownames_to_column("taxon") %>%
    arrange(wi.eBH)

  write_csv(aldex_broad_full,
            file.path(results_dir, "differential_abundance/aldex2_broad_spectrum.csv"))
  cat("\nSaved: differential_abundance/aldex2_broad_spectrum.csv\n")
}

# =============================================================================
# 5. MaAsLin2 Analysis - SKIPPED (not available for R 4.5)
# =============================================================================

cat("\n=== MaAsLin2: Skipped (package not available for R 4.5) ===\n\n")

# =============================================================================
# 6. Generate Volcano Plot
# =============================================================================

cat("\n=== Generating Volcano Plot ===\n\n")

# Use ALDEx2 results for any abx
volcano_data <- aldex_any_full %>%
  mutate(
    significant = wi.eBH < 0.1,
    label = ifelse(significant, taxon, NA)
  )

p_volcano <- ggplot(volcano_data, aes(x = effect, y = -log10(wi.ep))) +
  geom_point(aes(color = significant), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "#d73027")) +
  ggrepel::geom_text_repel(
    aes(label = label),
    size = 3,
    max.overlaps = 15,
    na.rm = TRUE
  ) +
  labs(
    title = "Differential Abundance: Any Antibiotic Exposure (7-day)",
    subtitle = "ALDEx2 analysis, high-confidence patient groups",
    x = "Effect Size",
    y = "-log10(p-value)",
    color = "Significant (BH < 0.1)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(results_dir, "figures/volcano_any_abx.pdf"), p_volcano,
       width = 8, height = 7)
cat("Saved: figures/volcano_any_abx.pdf\n")

# =============================================================================
# 7. Summary of Key Findings
# =============================================================================

cat("\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("SUMMARY OF KEY FINDINGS\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")

# Focus on biologically relevant taxa
key_genera <- c("Bacteroides", "Parabacteroides", "Prevotella", "Clostridium",
                "Faecalibacterium", "Blautia", "Roseburia", "Ruminococcus",
                "Escherichia", "Klebsiella", "Enterobacter", "Enterococcus",
                "Lactobacillus", "Streptococcus")

cat("Key genera results from ALDEx2 (Any Abx):\n")
aldex_any_full %>%
  filter(grepl(paste(key_genera, collapse = "|"), taxon)) %>%
  select(taxon, effect, wi.ep, wi.eBH) %>%
  mutate(
    direction = ifelse(effect > 0, "Increased", "Decreased"),
    across(where(is.numeric), ~ round(.x, 4))
  ) %>%
  arrange(wi.ep) %>%
  head(20) %>%
  print()

if (exists("aldex_anaerobic")) {
  cat("\nKey genera results from ALDEx2 (Anti-anaerobic Abx):\n")
  aldex_anaerobic_full %>%
    filter(grepl(paste(key_genera, collapse = "|"), taxon)) %>%
    select(taxon, effect, wi.ep, wi.eBH) %>%
    mutate(
      direction = ifelse(effect > 0, "Increased", "Decreased"),
      across(where(is.numeric), ~ round(.x, 4))
    ) %>%
    arrange(wi.ep) %>%
    head(20) %>%
    print()
}

cat("\nDifferential abundance analysis complete!\n")
