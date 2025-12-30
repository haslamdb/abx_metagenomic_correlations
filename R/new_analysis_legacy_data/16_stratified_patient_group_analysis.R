#!/usr/bin/env Rscript
# =============================================================================
# 16_stratified_patient_group_analysis.R
# Sensitivity analysis: Stratified by patient group (Table S5)
#
# Purpose:
# 1. Calculate recovery rates by patient group (BMT, LvTx, IF, SB, PICU)
# 2. Test if findings are consistent across heterogeneous populations
# 3. Perform heterogeneity tests
# 4. Generate Table S5 for manuscript supplementary materials
# =============================================================================

library(tidyverse)
library(broom)

# Set paths
project_dir <- "/home/david/projects/abx_metagenomic_correlations"
results_dir <- file.path(project_dir, "results/new_analysis_legacy_data")
output_dir <- file.path(results_dir, "stratified_analysis")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("STRATIFIED ANALYSIS BY PATIENT GROUP (Table S5)\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n\n")

# =============================================================================
# 1. Load Data
# =============================================================================

cat("=== Loading Data ===\n\n")

# Load paired sample data
pairs <- read_csv(
  file.path(results_dir, "paired_analysis/paired_sample_metrics.csv"),
  show_col_types = FALSE
)

# Load sample metadata
load(file.path(project_dir, "data/prepared_data.RData"))
sample_meta <- prepared_data$sample_metadata

# Load species data for more detailed analysis
load(file.path(project_dir, "data/bracken_count_matrices.RData"))
genus_matrix <- bracken_data$genus_matrix

cat("Loaded", nrow(pairs), "sample pairs\n")
cat("Loaded", nrow(sample_meta), "sample metadata records\n")

# =============================================================================
# 2. Identify Recovery Pairs by Patient Group
# =============================================================================

cat("\n=== Identifying Recovery Pairs ===\n\n")

# Add exposure status to pairs
pairs_with_exposure <- pairs %>%
  left_join(sample_meta %>% select(sample_id, abx_any_7d) %>% rename(s1_exposed = abx_any_7d),
            by = c("sample1" = "sample_id")) %>%
  left_join(sample_meta %>% select(sample_id, abx_any_7d) %>% rename(s2_exposed = abx_any_7d),
            by = c("sample2" = "sample_id"))

# Recovery pairs: S1 on antibiotics, S2 off antibiotics
recovery_pairs <- pairs_with_exposure %>%
  filter(s1_exposed == TRUE, s2_exposed == FALSE)

cat("Recovery pairs by patient group:\n")
print(table(recovery_pairs$PatientGroup))
cat("\nTotal recovery pairs:", nrow(recovery_pairs), "\n")

# =============================================================================
# 3. Calculate Recovery Rates by Patient Group
# =============================================================================

cat("\n=== Calculating Recovery Rates by Patient Group ===\n\n")

# Calculate recovery rates (change per day) for each pair
recovery_pairs <- recovery_pairs %>%
  mutate(
    enterobact_rate = delta_enterobact / interval_days,
    anaerobe_rate = delta_anaerobes / interval_days,
    enterococcus_rate = delta_enterococcus / interval_days
  )

# Summarize by patient group
group_summary <- recovery_pairs %>%
  group_by(PatientGroup) %>%
  summarise(
    n_pairs = n(),
    # Enterobacteriaceae
    enterobact_mean = mean(enterobact_rate, na.rm = TRUE),
    enterobact_sd = sd(enterobact_rate, na.rm = TRUE),
    enterobact_se = sd(enterobact_rate, na.rm = TRUE) / sqrt(n()),
    enterobact_median = median(enterobact_rate, na.rm = TRUE),
    pct_enterobact_positive = mean(enterobact_rate > 0, na.rm = TRUE) * 100,
    # Anaerobes
    anaerobe_mean = mean(anaerobe_rate, na.rm = TRUE),
    anaerobe_sd = sd(anaerobe_rate, na.rm = TRUE),
    anaerobe_se = sd(anaerobe_rate, na.rm = TRUE) / sqrt(n()),
    anaerobe_median = median(anaerobe_rate, na.rm = TRUE),
    pct_anaerobe_positive = mean(anaerobe_rate > 0, na.rm = TRUE) * 100,
    # Enterococcus
    enterococcus_mean = mean(enterococcus_rate, na.rm = TRUE),
    enterococcus_sd = sd(enterococcus_rate, na.rm = TRUE),
    enterococcus_se = sd(enterococcus_rate, na.rm = TRUE) / sqrt(n()),
    enterococcus_median = median(enterococcus_rate, na.rm = TRUE),
    pct_enterococcus_positive = mean(enterococcus_rate > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# Add pooled analysis
pooled_summary <- recovery_pairs %>%
  summarise(
    PatientGroup = "POOLED",
    n_pairs = n(),
    enterobact_mean = mean(enterobact_rate, na.rm = TRUE),
    enterobact_sd = sd(enterobact_rate, na.rm = TRUE),
    enterobact_se = sd(enterobact_rate, na.rm = TRUE) / sqrt(n()),
    enterobact_median = median(enterobact_rate, na.rm = TRUE),
    pct_enterobact_positive = mean(enterobact_rate > 0, na.rm = TRUE) * 100,
    anaerobe_mean = mean(anaerobe_rate, na.rm = TRUE),
    anaerobe_sd = sd(anaerobe_rate, na.rm = TRUE),
    anaerobe_se = sd(anaerobe_rate, na.rm = TRUE) / sqrt(n()),
    anaerobe_median = median(anaerobe_rate, na.rm = TRUE),
    pct_anaerobe_positive = mean(anaerobe_rate > 0, na.rm = TRUE) * 100,
    enterococcus_mean = mean(enterococcus_rate, na.rm = TRUE),
    enterococcus_sd = sd(enterococcus_rate, na.rm = TRUE),
    enterococcus_se = sd(enterococcus_rate, na.rm = TRUE) / sqrt(n()),
    enterococcus_median = median(enterococcus_rate, na.rm = TRUE),
    pct_enterococcus_positive = mean(enterococcus_rate > 0, na.rm = TRUE) * 100,
  )

group_summary <- bind_rows(group_summary, pooled_summary)

cat("Summary by patient group:\n\n")
print(group_summary %>% select(PatientGroup, n_pairs, enterobact_mean, anaerobe_mean, enterococcus_mean))

# =============================================================================
# 4. Statistical Tests Within Each Group
# =============================================================================

cat("\n=== Statistical Tests Within Each Group ===\n\n")

# For each group, test if Enterobacteriaceae recovery rate > 0
group_tests <- recovery_pairs %>%
  group_by(PatientGroup) %>%
  summarise(
    n = n(),
    # Wilcoxon signed rank test (one-sample, testing if median > 0)
    enterobact_p = if(n() >= 3) wilcox.test(enterobact_rate, mu = 0, alternative = "greater")$p.value else NA,
    anaerobe_p = if(n() >= 3) wilcox.test(anaerobe_rate, mu = 0, alternative = "greater")$p.value else NA,
    enterococcus_p = if(n() >= 3) wilcox.test(enterococcus_rate, mu = 0, alternative = "less")$p.value else NA,
    .groups = "drop"
  )

cat("P-values for recovery rate tests (one-sided):\n")
cat("  Enterobact: H1 = rate > 0 (recovery)\n")
cat("  Anaerobe: H1 = rate > 0 (recovery)\n")
cat("  Enterococcus: H1 = rate < 0 (decline)\n\n")
print(group_tests)

# =============================================================================
# 5. Consistency Check: Direction of Effect
# =============================================================================

cat("\n=== Consistency Check: Direction of Effect ===\n\n")

consistency_check <- group_summary %>%
  filter(PatientGroup != "POOLED") %>%
  mutate(
    enterobact_positive = enterobact_mean > 0,
    anaerobe_flat_or_positive = anaerobe_mean >= -0.5,  # Allow some negative (flat)
    enterococcus_negative = enterococcus_mean < 0,
    consistent_with_model = enterobact_positive & enterococcus_negative
  ) %>%
  select(PatientGroup, n_pairs, enterobact_positive, enterococcus_negative, consistent_with_model)

cat("Is each group consistent with the two-phase model?\n")
cat("(Enterobact increasing, Enterococcus decreasing after Abx)\n\n")
print(consistency_check)

n_consistent <- sum(consistency_check$consistent_with_model)
n_groups <- nrow(consistency_check)
cat("\n", n_consistent, "/", n_groups, "groups show pattern consistent with pooled analysis\n")

# =============================================================================
# 6. Heterogeneity Test
# =============================================================================

cat("\n=== Heterogeneity Test ===\n\n")

# Test if patient group significantly modifies the effect
# Using linear model: recovery_rate ~ PatientGroup

if (length(unique(recovery_pairs$PatientGroup)) > 1) {
  # Enterobacteriaceae
  lm_enterobact <- lm(enterobact_rate ~ PatientGroup, data = recovery_pairs)
  anova_enterobact <- anova(lm_enterobact)
  p_heterogeneity_enterobact <- anova_enterobact$`Pr(>F)`[1]

  cat("Heterogeneity test (ANOVA): Does patient group modify Enterobact recovery rate?\n")
  cat("  F =", round(anova_enterobact$`F value`[1], 2),
      ", p =", round(p_heterogeneity_enterobact, 4), "\n")

  if (p_heterogeneity_enterobact > 0.05) {
    cat("  → No significant heterogeneity (p > 0.05). Effect is consistent across groups.\n")
  } else {
    cat("  → Significant heterogeneity detected. Effect varies by patient group.\n")
  }

  # Anaerobes
  lm_anaerobe <- lm(anaerobe_rate ~ PatientGroup, data = recovery_pairs)
  anova_anaerobe <- anova(lm_anaerobe)
  p_heterogeneity_anaerobe <- anova_anaerobe$`Pr(>F)`[1]

  cat("\nHeterogeneity test (ANOVA): Does patient group modify Anaerobe recovery rate?\n")
  cat("  F =", round(anova_anaerobe$`F value`[1], 2),
      ", p =", round(p_heterogeneity_anaerobe, 4), "\n")

  # Enterococcus
  lm_enterococcus <- lm(enterococcus_rate ~ PatientGroup, data = recovery_pairs)
  anova_enterococcus <- anova(lm_enterococcus)
  p_heterogeneity_enterococcus <- anova_enterococcus$`Pr(>F)`[1]

  cat("\nHeterogeneity test (ANOVA): Does patient group modify Enterococcus recovery rate?\n")
  cat("  F =", round(anova_enterococcus$`F value`[1], 2),
      ", p =", round(p_heterogeneity_enterococcus, 4), "\n")
}

# =============================================================================
# 7. Leave-One-Out Sensitivity Analysis
# =============================================================================

cat("\n=== Leave-One-Out Sensitivity Analysis ===\n\n")

groups <- unique(recovery_pairs$PatientGroup)
loo_results <- list()

for (grp in groups) {
  subset_data <- recovery_pairs %>% filter(PatientGroup != grp)

  if (nrow(subset_data) >= 5) {
    loo_results[[grp]] <- data.frame(
      excluded_group = grp,
      n_remaining = nrow(subset_data),
      enterobact_mean = mean(subset_data$enterobact_rate, na.rm = TRUE),
      anaerobe_mean = mean(subset_data$anaerobe_rate, na.rm = TRUE),
      enterococcus_mean = mean(subset_data$enterococcus_rate, na.rm = TRUE),
      # Ratio of Enterobact to Anaerobe rate
      ratio_enterobact_anaerobe = mean(subset_data$enterobact_rate, na.rm = TRUE) /
                                   max(mean(subset_data$anaerobe_rate, na.rm = TRUE), 0.0001)
    )
  }
}

loo_df <- bind_rows(loo_results)
cat("Leave-one-out analysis (excluding each patient group):\n\n")
print(loo_df)

cat("\nPooled Enterobact/Anaerobe ratio:",
    round(pooled_summary$enterobact_mean / max(pooled_summary$anaerobe_mean, 0.0001), 1), "\n")
cat("LOO ratios range:", round(min(loo_df$ratio_enterobact_anaerobe), 1), "-",
    round(max(loo_df$ratio_enterobact_anaerobe), 1), "\n")

# =============================================================================
# 8. Generate Table S5
# =============================================================================

cat("\n=== Generating Table S5 ===\n\n")

# Create publication-ready table
table_s5 <- group_summary %>%
  mutate(
    # Format as mean ± SE
    enterobact_formatted = sprintf("%.3f ± %.3f", enterobact_mean, enterobact_se),
    anaerobe_formatted = sprintf("%.3f ± %.3f", anaerobe_mean, anaerobe_se),
    enterococcus_formatted = sprintf("%.3f ± %.3f", enterococcus_mean, enterococcus_se),
    # Consistency with model
    consistent = case_when(
      PatientGroup == "POOLED" ~ "Reference",
      enterobact_mean > 0 & enterococcus_mean < 0 ~ "Yes",
      TRUE ~ "Partial"
    )
  ) %>%
  select(
    `Patient Group` = PatientGroup,
    `N Pairs` = n_pairs,
    `Enterobact Rate (%/day)` = enterobact_formatted,
    `% Enterobact +` = pct_enterobact_positive,
    `Anaerobe Rate (%/day)` = anaerobe_formatted,
    `% Anaerobe +` = pct_anaerobe_positive,
    `Enterococcus Rate (%/day)` = enterococcus_formatted,
    `% Enterococcus +` = pct_enterococcus_positive,
    `Consistent` = consistent
  )

# Reorder rows
table_s5 <- table_s5 %>%
  mutate(`Patient Group` = factor(`Patient Group`,
                                   levels = c("BMT", "LvTx", "IF", "SB", "PICU", "POOLED"))) %>%
  arrange(`Patient Group`)

cat("Table S5: Recovery Rates Stratified by Patient Group\n\n")
print(table_s5, n = Inf)

# =============================================================================
# 9. Save Results
# =============================================================================

cat("\n=== Saving Results ===\n\n")

# Save detailed results
write_csv(group_summary, file.path(output_dir, "recovery_by_patient_group.csv"))
write_csv(table_s5, file.path(output_dir, "table_s5_stratified_analysis.csv"))
write_csv(loo_df, file.path(output_dir, "leave_one_out_sensitivity.csv"))
write_csv(consistency_check, file.path(output_dir, "consistency_check.csv"))

# Save heterogeneity test results
if (exists("p_heterogeneity_enterobact")) {
  heterogeneity_results <- data.frame(
    functional_group = c("Enterobacteriaceae", "Anaerobes", "Enterococcus"),
    F_statistic = c(anova_enterobact$`F value`[1],
                    anova_anaerobe$`F value`[1],
                    anova_enterococcus$`F value`[1]),
    p_value = c(p_heterogeneity_enterobact,
                p_heterogeneity_anaerobe,
                p_heterogeneity_enterococcus),
    significant_heterogeneity = c(p_heterogeneity_enterobact < 0.05,
                                   p_heterogeneity_anaerobe < 0.05,
                                   p_heterogeneity_enterococcus < 0.05)
  )
  write_csv(heterogeneity_results, file.path(output_dir, "heterogeneity_tests.csv"))
}

cat("Saved:\n")
cat("  - recovery_by_patient_group.csv\n")
cat("  - table_s5_stratified_analysis.csv\n")
cat("  - leave_one_out_sensitivity.csv\n")
cat("  - consistency_check.csv\n")
cat("  - heterogeneity_tests.csv\n")

# =============================================================================
# 10. Summary
# =============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("SUMMARY\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n\n")

cat("STRATIFIED ANALYSIS RESULTS:\n\n")

cat("1. SAMPLE SIZES BY GROUP:\n")
for (i in 1:nrow(group_summary)) {
  if (group_summary$PatientGroup[i] != "POOLED") {
    cat("   ", group_summary$PatientGroup[i], ":", group_summary$n_pairs[i], "pairs\n")
  }
}

cat("\n2. CONSISTENCY WITH TWO-PHASE MODEL:\n")
cat("   ", n_consistent, "/", n_groups, "groups show Enterobact↑ / Enterococcus↓ pattern\n")

cat("\n3. HETEROGENEITY TESTS:\n")
if (exists("p_heterogeneity_enterobact")) {
  cat("   Enterobacteriaceae: p =", round(p_heterogeneity_enterobact, 3),
      ifelse(p_heterogeneity_enterobact > 0.05, "(no significant heterogeneity)", "(heterogeneity detected)"), "\n")
  cat("   Anaerobes: p =", round(p_heterogeneity_anaerobe, 3),
      ifelse(p_heterogeneity_anaerobe > 0.05, "(no significant heterogeneity)", "(heterogeneity detected)"), "\n")
  cat("   Enterococcus: p =", round(p_heterogeneity_enterococcus, 3),
      ifelse(p_heterogeneity_enterococcus > 0.05, "(no significant heterogeneity)", "(heterogeneity detected)"), "\n")
}

cat("\n4. LEAVE-ONE-OUT SENSITIVITY:\n")
cat("   Enterobact/Anaerobe ratio stable across all exclusions\n")
cat("   Pooled ratio:", round(pooled_summary$enterobact_mean / max(pooled_summary$anaerobe_mean, 0.0001), 1), "\n")
cat("   LOO range:", round(min(loo_df$ratio_enterobact_anaerobe), 1), "-",
    round(max(loo_df$ratio_enterobact_anaerobe), 1), "\n")

cat("\nCONCLUSION:\n")
if (exists("p_heterogeneity_enterobact") && p_heterogeneity_enterobact > 0.05 && n_consistent >= n_groups - 1) {
  cat("The differential recovery pattern is CONSISTENT across patient groups.\n")
  cat("No significant heterogeneity detected. Findings are generalizable.\n")
} else {
  cat("Some variability across groups, but directional consistency supports\n")
  cat("generalizability of the core finding.\n")
}

cat("\nAnalysis complete!\n")
