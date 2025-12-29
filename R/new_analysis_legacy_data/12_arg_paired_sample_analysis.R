#!/usr/bin/env Rscript
# =============================================================================
# 12_arg_paired_sample_analysis.R
# Paired sample analysis for ARGs - within-patient changes with INDIVIDUAL
# antibiotic exposure
# Uses LMM with covariate adjustment for other antibiotics
#
# Outcomes analyzed:
#   - Total ARG abundance (normalized)
#   - ARG diversity (Shannon)
#   - ARG richness
#   - Key resistance class changes (tetracycline, beta-lactam, etc.)
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
paired_dir <- file.path(results_dir, "arg_paired_analysis")

dir.create(paired_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Load Data
# =============================================================================

cat("=============================================================================\n")
cat("ARG Paired Sample Analysis\n")
cat("=============================================================================\n\n")

# Load ARG count matrices
load(file.path(project_dir, "data/arg_count_matrices.RData"))

# Load antibiotic metadata + paired samples
load(file.path(project_dir, "data/prepared_data.RData"))

sample_metadata <- prepared_data$sample_metadata
paired_samples <- prepared_data$paired_samples
high_confidence_groups <- prepared_data$high_confidence_groups

# ARG data
arg_count_matrix <- arg_data$arg_count_matrix
arg_rpm_matrix <- arg_data$arg_rpm_matrix

cat("ARG count matrix:", nrow(arg_count_matrix), "samples x", ncol(arg_count_matrix), "ARGs\n")
cat("ARG RPM matrix:", nrow(arg_rpm_matrix), "samples x", ncol(arg_rpm_matrix), "ARGs\n")
cat("Sample metadata:", nrow(sample_metadata), "samples\n")
cat("Paired samples:", nrow(paired_samples), "pairs\n\n")

# Define individual antibiotics to analyze
individual_antibiotics <- c(
  "Pip_Tazo", "Meropenem", "Cefepime", "Ceftriaxone", "Ciprofloxacin",
  "Metronidazole", "Clindamycin", "TMP_SMX", "Vancomycin_IV", "Vancomycin_PO"
)

# Filter to overlapping samples
common_samples <- intersect(sample_metadata$sample_id, rownames(arg_count_matrix))
sample_metadata <- sample_metadata %>% filter(sample_id %in% common_samples)
arg_count_matrix <- arg_count_matrix[common_samples, ]
arg_rpm_matrix <- arg_rpm_matrix[intersect(common_samples, rownames(arg_rpm_matrix)), ]

# =============================================================================
# 2. Calculate ARG Summary Metrics Per Sample
# =============================================================================

cat("=== Calculating ARG Summary Metrics ===\n\n")

# Calculate diversity metrics on ARG data
arg_diversity <- data.frame(
  sample_id = rownames(arg_count_matrix),
  arg_richness = rowSums(arg_count_matrix > 0),
  arg_shannon = vegan::diversity(arg_count_matrix, index = "shannon"),
  arg_total_reads = rowSums(arg_count_matrix)
)

# Calculate RPM-based metrics where available
if (nrow(arg_rpm_matrix) > 0) {
  arg_diversity$arg_total_rpm <- rowSums(arg_rpm_matrix[arg_diversity$sample_id, ], na.rm = TRUE)
}

# Define ARG resistance classes
resistance_classes <- list(
  tetracycline = "^Tet\\(",
  macrolide = "Mef|Msr|Erm|Ere",
  betalactam = "CTX|TEM|SHV|OXA|KPC|NDM|VIM|IMP|AmpC",
  vancomycin = "^Van[A-Z]",
  aminoglycoside = "Aac|Aph|Ant|Aad",
  fluoroquinolone = "Qnr|GyrA|ParC",
  sulfonamide = "Sul[0-9]|Dfr"
)

# Calculate relative abundance of each resistance class
for (class_name in names(resistance_classes)) {
  pattern <- resistance_classes[[class_name]]
  matching_cols <- grep(pattern, colnames(arg_count_matrix), ignore.case = TRUE)

  if (length(matching_cols) > 0) {
    class_counts <- rowSums(arg_count_matrix[, matching_cols, drop = FALSE])
    total_counts <- rowSums(arg_count_matrix)
    arg_diversity[[paste0(class_name, "_rel")]] <- class_counts / (total_counts + 1)
  } else {
    arg_diversity[[paste0(class_name, "_rel")]] <- 0
  }
}

cat("ARG diversity metrics calculated for", nrow(arg_diversity), "samples\n")
cat("Resistance classes quantified:", paste(names(resistance_classes), collapse = ", "), "\n\n")

# Merge with sample metadata
sample_metadata <- sample_metadata %>%
  left_join(arg_diversity, by = "sample_id")

# =============================================================================
# 3. Prepare Paired Sample Data
# =============================================================================

cat("=== Preparing Paired Sample Data ===\n\n")

# Filter to high-confidence pairs
paired_hc <- paired_samples %>% filter(high_confidence)
cat("High-confidence sample pairs:", nrow(paired_hc), "\n")
cat("Pairs with any Abx between samples:", sum(paired_hc$abx_between_any), "\n\n")

# Calculate metrics for each pair
calculate_pair_metrics <- function(pair_row, arg_diversity_df) {

  sample1 <- pair_row$sample1
  sample2 <- pair_row$sample2

  # Get ARG diversity data
  div1 <- arg_diversity_df %>% filter(sample_id == sample1)
  div2 <- arg_diversity_df %>% filter(sample_id == sample2)

  if (nrow(div1) == 0 || nrow(div2) == 0) return(NULL)

  # Calculate changes (sample2 - sample1)
  result <- data.frame(
    delta_arg_richness = div2$arg_richness - div1$arg_richness,
    delta_arg_shannon = div2$arg_shannon - div1$arg_shannon,
    delta_arg_total = log2((div2$arg_total_reads + 1) / (div1$arg_total_reads + 1))
  )

  # Add resistance class changes (log2 fold change)
  for (class_name in names(resistance_classes)) {
    col <- paste0(class_name, "_rel")
    if (col %in% colnames(div1) && col %in% colnames(div2)) {
      result[[paste0("delta_", class_name)]] <- log2(
        (div2[[col]] + 1e-6) / (div1[[col]] + 1e-6)
      )
    }
  }

  result
}

cat("Calculating paired sample metrics...\n")

# Calculate for all pairs
pair_metrics_list <- list()
for (i in 1:nrow(paired_hc)) {
  metrics <- calculate_pair_metrics(paired_hc[i, ], sample_metadata)
  if (!is.null(metrics)) {
    pair_metrics_list[[i]] <- cbind(paired_hc[i, ], metrics)
  }
}

paired_results <- bind_rows(pair_metrics_list)
cat("Pairs with complete ARG data:", nrow(paired_results), "\n\n")

# =============================================================================
# 4. Summary of Individual Antibiotic Exposures
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
# 5. LMM Analysis for Each Individual Antibiotic
# =============================================================================

cat("\n=== LMM Analysis: Individual Antibiotic Effects on ARGs ===\n")
cat("Model: outcome ~ target_abx + other_abx_covariates + interval_days + (1|MRN)\n\n")

# Define outcomes to test
outcomes <- c("delta_arg_richness", "delta_arg_shannon", "delta_arg_total",
              "delta_tetracycline", "delta_macrolide", "delta_betalactam",
              "delta_vancomycin", "delta_aminoglycoside")

# Filter to outcomes that exist in data
outcomes <- outcomes[outcomes %in% colnames(paired_results)]

outcome_labels <- c(
  delta_arg_richness = "ARG Richness Change",
  delta_arg_shannon = "ARG Shannon Diversity Change",
  delta_arg_total = "Total ARG (log2FC)",
  delta_tetracycline = "Tetracycline Resistance (log2FC)",
  delta_macrolide = "Macrolide Resistance (log2FC)",
  delta_betalactam = "Beta-lactam Resistance (log2FC)",
  delta_vancomycin = "Vancomycin Resistance (log2FC)",
  delta_aminoglycoside = "Aminoglycoside Resistance (log2FC)"
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
      cat(sprintf("    %-30s: est = %7.3f, SE = %.3f, p = %.4f %s\n",
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
# 6. Summary Table: Significant Findings
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
# 7. Generate Plots
# =============================================================================

cat("\n=== Generating Plots ===\n\n")

# 7.1 Forest plot of antibiotic effects on ARG richness
richness_results <- combined_results %>%
  filter(outcome == "delta_arg_richness") %>%
  mutate(
    antibiotic = factor(antibiotic, levels = rev(individual_antibiotics)),
    significant = p.value < 0.05
  )

if (nrow(richness_results) > 0) {
  p_forest_richness <- ggplot(richness_results, aes(x = estimate, y = antibiotic)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    geom_point(aes(color = significant), size = 3) +
    scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "#d73027"),
                       labels = c("FALSE" = "p >= 0.05", "TRUE" = "p < 0.05")) +
    labs(
      title = "Effect of Individual Antibiotics on ARG Richness Change",
      subtitle = "LMM with covariate adjustment for other antibiotics",
      x = "Estimated Effect on Delta ARG Richness (95% CI)",
      y = NULL,
      color = "Significance"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(file.path(figures_dir, "arg_paired_individual_abx_richness.pdf"), p_forest_richness,
         width = 8, height = 6)
  cat("Saved: figures/arg_paired_individual_abx_richness.pdf\n")
}

# 7.2 Forest plot for tetracycline resistance
tet_results <- combined_results %>%
  filter(outcome == "delta_tetracycline") %>%
  mutate(
    antibiotic = factor(antibiotic, levels = rev(individual_antibiotics)),
    significant = p.value < 0.05
  )

if (nrow(tet_results) > 0) {
  p_forest_tet <- ggplot(tet_results, aes(x = estimate, y = antibiotic)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    geom_point(aes(color = significant), size = 3) +
    scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "#d73027"),
                       labels = c("FALSE" = "p >= 0.05", "TRUE" = "p < 0.05")) +
    labs(
      title = "Effect of Individual Antibiotics on Tetracycline Resistance",
      subtitle = "LMM with covariate adjustment for other antibiotics",
      x = "Estimated Effect on Tetracycline Resistance log2FC (95% CI)",
      y = NULL,
      color = "Significance"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(file.path(figures_dir, "arg_paired_individual_abx_tetracycline.pdf"), p_forest_tet,
         width = 8, height = 6)
  cat("Saved: figures/arg_paired_individual_abx_tetracycline.pdf\n")
}

# 7.3 Heatmap of all effects
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
      title = "Individual Antibiotic Effects on ARG Outcomes",
      subtitle = "Paired sample analysis with LMM (* p<0.05, ** p<0.01, . p<0.1)",
      x = NULL,
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  ggsave(file.path(figures_dir, "arg_paired_individual_abx_heatmap.pdf"), p_heatmap,
         width = 12, height = 7)
  cat("Saved: figures/arg_paired_individual_abx_heatmap.pdf\n")
}

# =============================================================================
# 8. Vancomycin Route Comparison
# =============================================================================

cat("\n=== Vancomycin Route Comparison (IV vs PO) ===\n\n")

vanco_iv_results <- combined_results %>% filter(antibiotic == "Vancomycin_IV")
vanco_po_results <- combined_results %>% filter(antibiotic == "Vancomycin_PO")

if (nrow(vanco_iv_results) > 0 && nrow(vanco_po_results) > 0) {
  comparison <- vanco_iv_results %>%
    select(outcome_label, iv_estimate = estimate, iv_pvalue = p.value) %>%
    left_join(
      vanco_po_results %>%
        select(outcome_label, po_estimate = estimate, po_pvalue = p.value),
      by = "outcome_label"
    )

  cat("Vancomycin IV vs PO effects:\n")
  print(as.data.frame(comparison))
}

# =============================================================================
# 9. Save Results
# =============================================================================

cat("\n=== Saving Results ===\n\n")

# Save full results
write_csv(combined_results, file.path(paired_dir, "individual_abx_lmm_results.csv"))
cat("Saved: arg_paired_analysis/individual_abx_lmm_results.csv\n")

# Save significant results
if (nrow(significant_results) > 0) {
  write_csv(significant_results, file.path(paired_dir, "individual_abx_significant.csv"))
  cat("Saved: arg_paired_analysis/individual_abx_significant.csv\n")
}

# Save paired sample metrics
write_csv(paired_results, file.path(paired_dir, "paired_sample_metrics.csv"))
cat("Saved: arg_paired_analysis/paired_sample_metrics.csv\n")

# Save exposure summary
write_csv(exposure_summary, file.path(paired_dir, "exposure_summary.csv"))
cat("Saved: arg_paired_analysis/exposure_summary.csv\n")

# Save ARG diversity metrics
write_csv(arg_diversity, file.path(paired_dir, "sample_arg_metrics.csv"))
cat("Saved: arg_paired_analysis/sample_arg_metrics.csv\n")

# =============================================================================
# 10. Summary
# =============================================================================

cat("\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("SUMMARY OF ARG PAIRED ANALYSIS - INDIVIDUAL ANTIBIOTICS\n")
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

cat("\nARG paired sample analysis complete!\n")
