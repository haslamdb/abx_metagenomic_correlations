#!/usr/bin/env Rscript
# =============================================================================
# 04_paired_sample_analysis.R
# Paired sample analysis - within-patient changes with antibiotic exposure
# =============================================================================

library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
library(broom.mixed)

# Set paths
project_dir <- "/home/david/projects/abx_metagenomic_correlations"
results_dir <- file.path(project_dir, "results/new_analysis_legacy_data")
figures_dir <- file.path(results_dir, "figures")

# Load data
load(file.path(project_dir, "data/bracken_count_matrices.RData"))  # Raw Bracken counts
load(file.path(project_dir, "data/prepared_data.RData"))  # Antibiotic metadata + paired samples

sample_metadata <- prepared_data$sample_metadata
species_matrix <- bracken_data$species_matrix  # Use raw Bracken counts
genus_matrix <- bracken_data$genus_matrix      # Use raw Bracken counts
paired_samples <- prepared_data$paired_samples
high_confidence_groups <- prepared_data$high_confidence_groups

# Filter to overlapping samples
common_samples <- intersect(sample_metadata$sample_id, rownames(species_matrix))
sample_metadata <- sample_metadata %>% filter(sample_id %in% common_samples)
species_matrix <- species_matrix[common_samples, ]
genus_matrix <- genus_matrix[intersect(common_samples, rownames(genus_matrix)), ]

# Recalculate diversity metrics on raw Bracken data
sample_metadata$shannon <- vegan::diversity(species_matrix, index = "shannon")
sample_metadata$richness <- rowSums(species_matrix > 0)

cat("Using raw Bracken counts:\n")
cat("  Species:", nrow(species_matrix), "samples x", ncol(species_matrix), "species\n")
cat("  Genus:", nrow(genus_matrix), "samples x", ncol(genus_matrix), "genera\n\n")

# =============================================================================
# 1. Calculate Paired Sample Metrics
# =============================================================================

cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("PAIRED SAMPLE ANALYSIS\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")

# Filter to high-confidence pairs
paired_hc <- paired_samples %>% filter(high_confidence)
cat("High-confidence sample pairs:", nrow(paired_hc), "\n")
cat("Pairs with Abx between samples:", sum(paired_hc$abx_between_any), "\n\n")

# Get diversity and functional group data
diversity_data <- sample_metadata %>%
  select(sample_id, shannon, richness, anaerobes_rel, enterobact_rel, enterococcus_rel)

# Calculate metrics for each pair
calculate_pair_metrics <- function(pair_row, species_mat, diversity_df) {

  sample1 <- pair_row$sample1
  sample2 <- pair_row$sample2

  # Check if samples exist in matrix
  if (!sample1 %in% rownames(species_mat) || !sample2 %in% rownames(species_mat)) {
    return(NULL)
  }

  # Get diversity data
  div1 <- diversity_df %>% filter(sample_id == sample1)
  div2 <- diversity_df %>% filter(sample_id == sample2)

  if (nrow(div1) == 0 || nrow(div2) == 0) return(NULL)

  # Calculate Bray-Curtis distance between the pair
  pair_matrix <- species_mat[c(sample1, sample2), ]
  pair_rel <- pair_matrix / rowSums(pair_matrix)
  bray_dist <- vegdist(pair_rel, method = "bray")[1]

  # Calculate changes (sample2 - sample1, since sample2 is later)
  data.frame(
    delta_shannon = div2$shannon - div1$shannon,
    delta_richness = div2$richness - div1$richness,
    delta_anaerobes = log2((div2$anaerobes_rel + 1e-6) / (div1$anaerobes_rel + 1e-6)),
    delta_enterobact = log2((div2$enterobact_rel + 1e-6) / (div1$enterobact_rel + 1e-6)),
    delta_enterococcus = log2((div2$enterococcus_rel + 1e-6) / (div1$enterococcus_rel + 1e-6)),
    bray_distance = bray_dist
  )
}

cat("Calculating paired sample metrics...\n")

# Calculate for all pairs
pair_metrics_list <- list()
for (i in 1:nrow(paired_hc)) {
  metrics <- calculate_pair_metrics(paired_hc[i, ], species_matrix, diversity_data)
  if (!is.null(metrics)) {
    pair_metrics_list[[i]] <- cbind(paired_hc[i, ], metrics)
  }
}

paired_results <- bind_rows(pair_metrics_list)
cat("Pairs with complete data:", nrow(paired_results), "\n\n")

# =============================================================================
# 2. Hypothesis Tests
# =============================================================================

cat("=== Statistical Tests ===\n\n")

# -----------------------------------------------------------------------------
# H1: Antibiotics reduce diversity (Shannon)
# -----------------------------------------------------------------------------

cat("H1: Antibiotic exposure reduces diversity change\n")
cat("   (More negative delta_shannon with Abx exposure)\n\n")

# Wilcoxon test
h1_wilcox <- wilcox.test(
  delta_shannon ~ abx_between_any,
  data = paired_results,
  alternative = "less"
)

cat("Wilcoxon test: delta_shannon by Abx exposure\n")
cat("  Median (no Abx):", median(paired_results$delta_shannon[!paired_results$abx_between_any]), "\n")
cat("  Median (Abx):", median(paired_results$delta_shannon[paired_results$abx_between_any]), "\n")
cat("  p-value:", format.pval(h1_wilcox$p.value), "\n\n")

# Mixed-effects model
h1_model <- lmer(
  delta_shannon ~ abx_between_any + interval_days + (1|MRN),
  data = paired_results
)
cat("Mixed-effects model: delta_shannon ~ Abx + interval + (1|patient)\n")
print(tidy(h1_model, effects = "fixed", conf.int = TRUE) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))))

# -----------------------------------------------------------------------------
# H2: Anti-anaerobic antibiotics reduce anaerobe abundance
# -----------------------------------------------------------------------------

cat("\n\nH2: Anti-anaerobic antibiotics reduce anaerobe abundance\n")
cat("   (More negative delta_anaerobes with anti-anaerobic Abx)\n\n")

h2_wilcox <- wilcox.test(
  delta_anaerobes ~ abx_between_anaerobic,
  data = paired_results,
  alternative = "less"
)

cat("Wilcoxon test: delta_anaerobes (log2FC) by anti-anaerobic Abx\n")
cat("  Median (no anti-anaerobic):", round(median(paired_results$delta_anaerobes[!paired_results$abx_between_anaerobic]), 3), "\n")
cat("  Median (anti-anaerobic):", round(median(paired_results$delta_anaerobes[paired_results$abx_between_anaerobic]), 3), "\n")
cat("  p-value:", format.pval(h2_wilcox$p.value), "\n\n")

h2_model <- lmer(
  delta_anaerobes ~ abx_between_anaerobic + interval_days + (1|MRN),
  data = paired_results
)
cat("Mixed-effects model: delta_anaerobes ~ anti-anaerobic Abx + interval + (1|patient)\n")
print(tidy(h2_model, effects = "fixed", conf.int = TRUE) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))))

# -----------------------------------------------------------------------------
# H3: Broad-spectrum antibiotics increase Enterobacteriaceae
# -----------------------------------------------------------------------------

cat("\n\nH3: Broad-spectrum antibiotics increase Enterobacteriaceae\n")
cat("   (More positive delta_enterobact with broad-spectrum Abx)\n\n")

h3_wilcox <- wilcox.test(
  delta_enterobact ~ abx_between_broad,
  data = paired_results,
  alternative = "greater"
)

cat("Wilcoxon test: delta_enterobact (log2FC) by broad-spectrum Abx\n")
cat("  Median (no broad):", round(median(paired_results$delta_enterobact[!paired_results$abx_between_broad]), 3), "\n")
cat("  Median (broad):", round(median(paired_results$delta_enterobact[paired_results$abx_between_broad]), 3), "\n")
cat("  p-value:", format.pval(h3_wilcox$p.value), "\n\n")

h3_model <- lmer(
  delta_enterobact ~ abx_between_broad + interval_days + (1|MRN),
  data = paired_results
)
cat("Mixed-effects model: delta_enterobact ~ broad-spectrum Abx + interval + (1|patient)\n")
print(tidy(h3_model, effects = "fixed", conf.int = TRUE) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))))

# -----------------------------------------------------------------------------
# H4: Antibiotics increase community instability (Bray-Curtis distance)
# -----------------------------------------------------------------------------

cat("\n\nH4: Antibiotic exposure increases community instability\n")
cat("   (Higher Bray-Curtis distance with Abx exposure)\n\n")

h4_wilcox <- wilcox.test(
  bray_distance ~ abx_between_any,
  data = paired_results,
  alternative = "greater"
)

cat("Wilcoxon test: Bray-Curtis distance by Abx exposure\n")
cat("  Median (no Abx):", round(median(paired_results$bray_distance[!paired_results$abx_between_any]), 3), "\n")
cat("  Median (Abx):", round(median(paired_results$bray_distance[paired_results$abx_between_any]), 3), "\n")
cat("  p-value:", format.pval(h4_wilcox$p.value), "\n\n")

h4_model <- lmer(
  bray_distance ~ abx_between_any + interval_days + (1|MRN),
  data = paired_results
)
cat("Mixed-effects model: Bray-Curtis ~ Abx + interval + (1|patient)\n")
print(tidy(h4_model, effects = "fixed", conf.int = TRUE) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))))

# -----------------------------------------------------------------------------
# H5: Dose-response relationship
# -----------------------------------------------------------------------------

cat("\n\nH5: Dose-response relationship (antibiotic-days)\n\n")

# Only for pairs with antibiotic exposure
exposed_pairs <- paired_results %>% filter(abx_between_any)

h5_cor_shannon <- cor.test(
  exposed_pairs$abx_between_days,
  exposed_pairs$delta_shannon,
  method = "spearman"
)

cat("Spearman correlation: Abx-days vs delta_shannon (exposed pairs only)\n")
cat("  rho:", round(h5_cor_shannon$estimate, 3), "\n")
cat("  p-value:", format.pval(h5_cor_shannon$p.value), "\n\n")

h5_cor_bray <- cor.test(
  exposed_pairs$abx_between_days,
  exposed_pairs$bray_distance,
  method = "spearman"
)

cat("Spearman correlation: Abx-days vs Bray-Curtis distance\n")
cat("  rho:", round(h5_cor_bray$estimate, 3), "\n")
cat("  p-value:", format.pval(h5_cor_bray$p.value), "\n")

# =============================================================================
# 3. Generate Plots
# =============================================================================

cat("\n=== Generating Plots ===\n\n")

# -----------------------------------------------------------------------------
# 3.1 Paired Trajectory Plots
# -----------------------------------------------------------------------------

# Shannon diversity change
p_shannon <- ggplot(paired_results, aes(x = abx_between_any, y = delta_shannon, fill = abx_between_any)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  scale_fill_manual(
    values = c("FALSE" = "#4575b4", "TRUE" = "#d73027"),
    labels = c("FALSE" = "No Abx", "TRUE" = "Abx between samples")
  ) +
  labs(
    title = "Change in Shannon Diversity Between Paired Samples",
    subtitle = paste("n =", nrow(paired_results), "pairs from high-confidence patient groups"),
    x = "Antibiotic Exposure Between Samples",
    y = "Change in Shannon Diversity (Sample 2 - Sample 1)",
    fill = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(figures_dir, "paired_delta_shannon.pdf"), p_shannon,
       width = 6, height = 5)
cat("Saved: figures/paired_delta_shannon.pdf\n")

# Anaerobe change by anti-anaerobic exposure
p_anaerobes <- ggplot(paired_results, aes(x = abx_between_anaerobic, y = delta_anaerobes, fill = abx_between_anaerobic)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  scale_fill_manual(
    values = c("FALSE" = "#4575b4", "TRUE" = "#d73027"),
    labels = c("FALSE" = "No anti-anaerobic", "TRUE" = "Anti-anaerobic Abx")
  ) +
  labs(
    title = "Change in Obligate Anaerobe Abundance",
    subtitle = "By anti-anaerobic antibiotic exposure between samples",
    x = "Anti-anaerobic Antibiotic Exposure",
    y = "Log2 Fold Change in Anaerobe Abundance",
    fill = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(figures_dir, "paired_delta_anaerobes.pdf"), p_anaerobes,
       width = 6, height = 5)
cat("Saved: figures/paired_delta_anaerobes.pdf\n")

# Bray-Curtis distance
p_bray <- ggplot(paired_results, aes(x = abx_between_any, y = bray_distance, fill = abx_between_any)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  scale_fill_manual(
    values = c("FALSE" = "#4575b4", "TRUE" = "#d73027"),
    labels = c("FALSE" = "No Abx", "TRUE" = "Abx between samples")
  ) +
  labs(
    title = "Community Stability: Bray-Curtis Distance Between Paired Samples",
    subtitle = "Higher values indicate more change between samples",
    x = "Antibiotic Exposure Between Samples",
    y = "Bray-Curtis Distance",
    fill = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(figures_dir, "paired_bray_curtis.pdf"), p_bray,
       width = 6, height = 5)
cat("Saved: figures/paired_bray_curtis.pdf\n")

# -----------------------------------------------------------------------------
# 3.2 Dose-Response Plot
# -----------------------------------------------------------------------------

p_dose <- ggplot(exposed_pairs, aes(x = abx_between_days, y = bray_distance)) +
  geom_point(alpha = 0.6, size = 2, color = "#d73027") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  labs(
    title = "Dose-Response: Antibiotic Days vs Community Change",
    subtitle = paste("Spearman rho =", round(h5_cor_bray$estimate, 3),
                     ", p =", format.pval(h5_cor_bray$p.value)),
    x = "Antibiotic Days Between Samples",
    y = "Bray-Curtis Distance"
  ) +
  theme_bw()

ggsave(file.path(figures_dir, "dose_response_bray.pdf"), p_dose,
       width = 6, height = 5)
cat("Saved: figures/dose_response_bray.pdf\n")

# -----------------------------------------------------------------------------
# 3.3 Combined Panel Plot
# -----------------------------------------------------------------------------

# Faceted plot of all metrics
metrics_long <- paired_results %>%
  select(PairNumber, abx_between_any, delta_shannon, delta_anaerobes, delta_enterobact, bray_distance) %>%
  pivot_longer(
    cols = c(delta_shannon, delta_anaerobes, delta_enterobact, bray_distance),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = case_when(
      metric == "delta_shannon" ~ "Shannon Diversity Change",
      metric == "delta_anaerobes" ~ "Anaerobes (log2FC)",
      metric == "delta_enterobact" ~ "Enterobacteriaceae (log2FC)",
      metric == "bray_distance" ~ "Bray-Curtis Distance"
    ),
    metric = factor(metric, levels = c(
      "Shannon Diversity Change", "Anaerobes (log2FC)",
      "Enterobacteriaceae (log2FC)", "Bray-Curtis Distance"
    ))
  )

p_panel <- ggplot(metrics_long, aes(x = abx_between_any, y = value, fill = abx_between_any)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(
    values = c("FALSE" = "#4575b4", "TRUE" = "#d73027"),
    labels = c("FALSE" = "No Abx", "TRUE" = "Abx")
  ) +
  labs(
    title = "Paired Sample Analysis: Microbiome Changes by Antibiotic Exposure",
    subtitle = paste("n =", nrow(paired_results), "high-confidence sample pairs"),
    x = "Antibiotic Exposure Between Samples",
    y = "Metric Value",
    fill = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(file.path(figures_dir, "paired_analysis_panel.pdf"), p_panel,
       width = 12, height = 4)
cat("Saved: figures/paired_analysis_panel.pdf\n")

# =============================================================================
# 4. Analysis by Patient Group
# =============================================================================

cat("\n=== Analysis by Patient Group ===\n\n")

group_summary <- paired_results %>%
  group_by(PatientGroup, abx_between_any) %>%
  summarise(
    n_pairs = n(),
    median_delta_shannon = median(delta_shannon),
    median_delta_anaerobes = median(delta_anaerobes),
    median_bray = median(bray_distance),
    .groups = "drop"
  )

print(group_summary)

# =============================================================================
# 5. Save Results
# =============================================================================

cat("\n=== Saving Results ===\n\n")

# Save paired sample metrics
write_csv(paired_results, file.path(results_dir, "paired_analysis/paired_sample_metrics.csv"))
cat("Saved: paired_analysis/paired_sample_metrics.csv\n")

# Compile hypothesis test results
hypothesis_tests <- data.frame(
  hypothesis = c(
    "H1: Abx reduces diversity",
    "H2: Anti-anaerobic reduces anaerobes",
    "H3: Broad-spectrum increases Enterobact",
    "H4: Abx increases instability",
    "H5a: Dose-response (Shannon)",
    "H5b: Dose-response (Bray-Curtis)"
  ),
  test = c(
    "Wilcoxon (one-sided)",
    "Wilcoxon (one-sided)",
    "Wilcoxon (one-sided)",
    "Wilcoxon (one-sided)",
    "Spearman correlation",
    "Spearman correlation"
  ),
  statistic = c(
    h1_wilcox$statistic,
    h2_wilcox$statistic,
    h3_wilcox$statistic,
    h4_wilcox$statistic,
    h5_cor_shannon$estimate,
    h5_cor_bray$estimate
  ),
  p_value = c(
    h1_wilcox$p.value,
    h2_wilcox$p.value,
    h3_wilcox$p.value,
    h4_wilcox$p.value,
    h5_cor_shannon$p.value,
    h5_cor_bray$p.value
  )
)

write_csv(hypothesis_tests, file.path(results_dir, "paired_analysis/hypothesis_tests.csv"))
cat("Saved: paired_analysis/hypothesis_tests.csv\n")

# =============================================================================
# 6. Summary
# =============================================================================

cat("\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("SUMMARY OF PAIRED ANALYSIS\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")

cat("Sample pairs analyzed:", nrow(paired_results), "\n")
cat("  - With Abx exposure:", sum(paired_results$abx_between_any), "\n")
cat("  - Without Abx:", sum(!paired_results$abx_between_any), "\n\n")

cat("Key Findings:\n")
cat("-" %>% rep(40) %>% paste(collapse = ""), "\n")

for (i in 1:nrow(hypothesis_tests)) {
  sig <- ifelse(hypothesis_tests$p_value[i] < 0.05, "**",
                ifelse(hypothesis_tests$p_value[i] < 0.1, "*", ""))
  cat(sprintf("%-45s p = %.4f %s\n",
              hypothesis_tests$hypothesis[i],
              hypothesis_tests$p_value[i],
              sig))
}

cat("\n* p < 0.10, ** p < 0.05\n")
cat("\nPaired sample analysis complete!\n")
