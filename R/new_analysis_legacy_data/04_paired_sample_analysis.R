#!/usr/bin/env Rscript
# =============================================================================
# 04_paired_sample_analysis.R
# Paired sample analysis - within-patient changes with INDIVIDUAL antibiotic exposure
# Uses LMM with covariate adjustment for other antibiotics
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
paired_dir <- file.path(results_dir, "paired_analysis")

dir.create(paired_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
load(file.path(project_dir, "data/bracken_count_matrices.RData"))  # Raw Bracken counts
load(file.path(project_dir, "data/prepared_data.RData"))  # Antibiotic metadata + paired samples

sample_metadata <- prepared_data$sample_metadata
species_matrix <- bracken_data$species_matrix  # Use raw Bracken counts
genus_matrix <- bracken_data$genus_matrix      # Use raw Bracken counts
paired_samples <- prepared_data$paired_samples
high_confidence_groups <- prepared_data$high_confidence_groups

# Define individual antibiotics to analyze (matching script 05)
individual_antibiotics <- c(
  "Pip_Tazo", "Meropenem", "Cefepime", "Ceftriaxone", "Ciprofloxacin",
  "Metronidazole", "Clindamycin", "TMP_SMX", "Vancomycin_IV", "Vancomycin_PO"
)

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
cat("PAIRED SAMPLE ANALYSIS - INDIVIDUAL ANTIBIOTICS\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")

# Filter to high-confidence pairs
paired_hc <- paired_samples %>% filter(high_confidence)
cat("High-confidence sample pairs:", nrow(paired_hc), "\n")
cat("Pairs with any Abx between samples:", sum(paired_hc$abx_between_any), "\n\n")

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
# 2. Summary of Individual Antibiotic Exposures
# =============================================================================

cat("=== Individual Antibiotic Exposure Summary ===\n\n")

exposure_summary <- data.frame(
  antibiotic = individual_antibiotics,
  n_pairs_exposed = sapply(individual_antibiotics, function(abx) {
    col <- paste0(abx, "_between")
    if (col %in% colnames(paired_results)) {
      sum(paired_results[[col]] > 0)
    } else {
      NA
    }
  }),
  n_pairs_unexposed = sapply(individual_antibiotics, function(abx) {
    col <- paste0(abx, "_between")
    if (col %in% colnames(paired_results)) {
      sum(paired_results[[col]] == 0)
    } else {
      NA
    }
  })
)
exposure_summary$pct_exposed <- round(100 * exposure_summary$n_pairs_exposed /
                                        nrow(paired_results), 1)
print(exposure_summary)

# =============================================================================
# 3. LMM Analysis for Each Individual Antibiotic
# =============================================================================

cat("\n=== LMM Analysis: Individual Antibiotic Effects ===\n")
cat("Model: outcome ~ target_abx + other_abx_covariates + interval_days + (1|MRN)\n\n")

# Define outcomes to test
outcomes <- c("delta_shannon", "delta_anaerobes", "delta_enterobact",
              "delta_enterococcus", "bray_distance")

outcome_labels <- c(
  delta_shannon = "Shannon Diversity Change",
  delta_anaerobes = "Anaerobes (log2FC)",
  delta_enterobact = "Enterobacteriaceae (log2FC)",
  delta_enterococcus = "Enterococcus (log2FC)",
  bray_distance = "Bray-Curtis Distance"
)

# Store all results
all_lmm_results <- list()

for (abx in individual_antibiotics) {

  target_col <- paste0(abx, "_between")

  # Check if column exists and has variation
  if (!target_col %in% colnames(paired_results)) {
    cat("Skipping", abx, "- column not found\n")
    next
  }

  n_exposed <- sum(paired_results[[target_col]] > 0)
  n_unexposed <- sum(paired_results[[target_col]] == 0)

  cat("=== ", abx, " ===\n", sep = "")
  cat("  Exposed pairs:", n_exposed, ", Unexposed:", n_unexposed, "\n")

  # Skip if too few exposed
  if (n_exposed < 5) {
    cat("  Skipping - too few exposed pairs\n\n")
    next
  }

  # Define covariates: other antibiotics + interval
  other_abx <- setdiff(individual_antibiotics, abx)
  other_cols <- paste0(other_abx, "_between")
  other_cols <- other_cols[other_cols %in% colnames(paired_results)]

  # Create binary exposure variable for target (exposed vs not)
  paired_results[[paste0(abx, "_exposed")]] <- paired_results[[target_col]] > 0

  # Test each outcome
  abx_results <- list()

  for (outcome in outcomes) {

    # Build formula with covariates
    # Use binary exposure for target, continuous (days) for others
    covariate_terms <- paste(other_cols, collapse = " + ")
    formula_str <- paste0(outcome, " ~ ", abx, "_exposed + ", covariate_terms,
                          " + interval_days + (1|MRN)")

    tryCatch({
      model <- lmer(as.formula(formula_str), data = paired_results,
                    control = lmerControl(optimizer = "bobyqa"))

      # Extract results for target antibiotic
      model_summary <- tidy(model, effects = "fixed", conf.int = TRUE) %>%
        filter(grepl(paste0(abx, "_exposed"), term))

      if (nrow(model_summary) > 0) {
        abx_results[[outcome]] <- model_summary %>%
          mutate(
            antibiotic = abx,
            outcome = outcome,
            outcome_label = outcome_labels[outcome],
            n_exposed = n_exposed,
            n_unexposed = n_unexposed
          )
      }
    }, error = function(e) {
      cat("    ", outcome, ": Model failed -", conditionMessage(e), "\n")
    })
  }

  if (length(abx_results) > 0) {
    all_lmm_results[[abx]] <- bind_rows(abx_results)

    # Print summary for this antibiotic
    cat("\n  Results (estimate = effect of exposure on outcome):\n")
    summary_df <- bind_rows(abx_results) %>%
      select(outcome_label, estimate, std.error, p.value) %>%
      mutate(
        sig = case_when(
          p.value < 0.01 ~ "**",
          p.value < 0.05 ~ "*",
          p.value < 0.1 ~ ".",
          TRUE ~ ""
        )
      )
    for (i in 1:nrow(summary_df)) {
      cat(sprintf("    %-25s: est = %7.3f, SE = %.3f, p = %.4f %s\n",
                  summary_df$outcome_label[i],
                  summary_df$estimate[i],
                  summary_df$std.error[i],
                  summary_df$p.value[i],
                  summary_df$sig[i]))
    }
  }

  cat("\n")
}

# Combine all results
combined_results <- bind_rows(all_lmm_results)

# =============================================================================
# 4. Summary Table: Significant Findings
# =============================================================================

cat("\n=== Summary: Significant Findings (p < 0.05) ===\n\n")

significant_results <- combined_results %>%
  filter(p.value < 0.05) %>%
  arrange(p.value) %>%
  select(antibiotic, outcome_label, estimate, std.error, conf.low, conf.high, p.value)

if (nrow(significant_results) > 0) {
  print(as.data.frame(significant_results))
} else {
  cat("No significant findings at p < 0.05\n")
}

cat("\n=== Marginal Findings (0.05 <= p < 0.1) ===\n\n")

marginal_results <- combined_results %>%
  filter(p.value >= 0.05 & p.value < 0.1) %>%
  arrange(p.value) %>%
  select(antibiotic, outcome_label, estimate, std.error, p.value)

if (nrow(marginal_results) > 0) {
  print(as.data.frame(marginal_results))
} else {
  cat("No marginal findings\n")
}

# =============================================================================
# 5. Generate Plots
# =============================================================================

cat("\n=== Generating Plots ===\n\n")

# 5.1 Forest plot of antibiotic effects on Shannon diversity
shannon_results <- combined_results %>%
  filter(outcome == "delta_shannon") %>%
  mutate(
    antibiotic = factor(antibiotic, levels = rev(individual_antibiotics)),
    significant = p.value < 0.05
  )

if (nrow(shannon_results) > 0) {
  p_forest_shannon <- ggplot(shannon_results, aes(x = estimate, y = antibiotic)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    geom_point(aes(color = significant), size = 3) +
    scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "#d73027"),
                       labels = c("FALSE" = "p >= 0.05", "TRUE" = "p < 0.05")) +
    labs(
      title = "Effect of Individual Antibiotics on Shannon Diversity Change",
      subtitle = "LMM with covariate adjustment for other antibiotics",
      x = "Estimated Effect on Delta Shannon (95% CI)",
      y = NULL,
      color = "Significance"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(file.path(figures_dir, "paired_individual_abx_shannon.pdf"), p_forest_shannon,
         width = 8, height = 6)
  cat("Saved: figures/paired_individual_abx_shannon.pdf\n")
}

# 5.2 Forest plot for Enterobacteriaceae
enterobact_results <- combined_results %>%
  filter(outcome == "delta_enterobact") %>%
  mutate(
    antibiotic = factor(antibiotic, levels = rev(individual_antibiotics)),
    significant = p.value < 0.05
  )

if (nrow(enterobact_results) > 0) {
  p_forest_enterobact <- ggplot(enterobact_results, aes(x = estimate, y = antibiotic)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    geom_point(aes(color = significant), size = 3) +
    scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "#d73027"),
                       labels = c("FALSE" = "p >= 0.05", "TRUE" = "p < 0.05")) +
    labs(
      title = "Effect of Individual Antibiotics on Enterobacteriaceae Change",
      subtitle = "LMM with covariate adjustment for other antibiotics",
      x = "Estimated Effect on Enterobacteriaceae log2FC (95% CI)",
      y = NULL,
      color = "Significance"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(file.path(figures_dir, "paired_individual_abx_enterobact.pdf"), p_forest_enterobact,
         width = 8, height = 6)
  cat("Saved: figures/paired_individual_abx_enterobact.pdf\n")
}

# 5.3 Heatmap of all effects
if (nrow(combined_results) > 0) {
  heatmap_data <- combined_results %>%
    mutate(
      sig_label = case_when(
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        p.value < 0.1 ~ ".",
        TRUE ~ ""
      )
    ) %>%
    select(antibiotic, outcome_label, estimate, sig_label)

  p_heatmap <- ggplot(heatmap_data, aes(x = outcome_label, y = antibiotic, fill = estimate)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sig_label), size = 4, color = "black") +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0,
      name = "Effect\nEstimate"
    ) +
    labs(
      title = "Individual Antibiotic Effects on Microbiome Outcomes",
      subtitle = "Paired sample analysis with LMM (* p<0.05, ** p<0.01, . p<0.1)",
      x = NULL,
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  ggsave(file.path(figures_dir, "paired_individual_abx_heatmap.pdf"), p_heatmap,
         width = 10, height = 7)
  cat("Saved: figures/paired_individual_abx_heatmap.pdf\n")
}

# =============================================================================
# 6. Save Results
# =============================================================================

cat("\n=== Saving Results ===\n\n")

# Save full results
write_csv(combined_results, file.path(paired_dir, "individual_abx_lmm_results.csv"))
cat("Saved: paired_analysis/individual_abx_lmm_results.csv\n")

# Save significant results
write_csv(significant_results, file.path(paired_dir, "individual_abx_significant.csv"))
cat("Saved: paired_analysis/individual_abx_significant.csv\n")

# Save paired sample metrics with individual exposures
write_csv(paired_results, file.path(paired_dir, "paired_sample_metrics.csv"))
cat("Saved: paired_analysis/paired_sample_metrics.csv\n")

# Save exposure summary
write_csv(exposure_summary, file.path(paired_dir, "exposure_summary.csv"))
cat("Saved: paired_analysis/exposure_summary.csv\n")

# =============================================================================
# 7. Summary
# =============================================================================

cat("\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("SUMMARY OF PAIRED ANALYSIS - INDIVIDUAL ANTIBIOTICS\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")

cat("Sample pairs analyzed:", nrow(paired_results), "\n")
cat("Individual antibiotics tested:", length(individual_antibiotics), "\n")
cat("Outcomes tested:", length(outcomes), "\n\n")

cat("Significant findings (p < 0.05):", nrow(significant_results), "\n")
cat("Marginal findings (0.05 <= p < 0.1):", nrow(marginal_results), "\n\n")

if (nrow(significant_results) > 0) {
  cat("Key significant associations:\n")
  cat("-" %>% rep(40) %>% paste(collapse = ""), "\n")
  for (i in 1:min(10, nrow(significant_results))) {
    direction <- ifelse(significant_results$estimate[i] > 0, "increases", "decreases")
    cat(sprintf("  %s %s %s (est=%.3f, p=%.4f)\n",
                significant_results$antibiotic[i],
                direction,
                significant_results$outcome_label[i],
                significant_results$estimate[i],
                significant_results$p.value[i]))
  }
}

cat("\nPaired sample analysis with individual antibiotics complete!\n")
