#!/usr/bin/env Rscript
# =============================================================================
# 07_recovery_analysis.R
# Microbiome recovery after antibiotic cessation
# Analyzes paired samples where antibiotics were given before sample 1
# but NO meaningful antibiotics between samples (excludes TMP_SMX, Azithromycin)
# =============================================================================

library(tidyverse)
library(vegan)

# Set paths
project_dir <- "/home/david/projects/abx_metagenomic_correlations"
results_dir <- file.path(project_dir, "results/new_analysis_legacy_data")
figures_dir <- file.path(results_dir, "figures")
recovery_dir <- file.path(results_dir, "recovery_analysis")

dir.create(recovery_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
load(file.path(project_dir, "data/bracken_count_matrices.RData"))
load(file.path(project_dir, "data/prepared_data.RData"))

paired_samples <- prepared_data$paired_samples
sample_metadata <- prepared_data$sample_metadata
genus_matrix <- bracken_data$genus_matrix

cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("MICROBIOME RECOVERY ANALYSIS\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")

# =============================================================================
# 1. Identify Recovery Pairs
# =============================================================================

cat("=== Identifying Recovery Pairs ===\n\n")

# Define "meaningful" antibiotics (exclude TMP_SMX and Azithromycin which have
# minimal microbiome impact based on our prior analyses)
meaningful_abx_cols <- c(
  "Pip_Tazo_between", "Meropenem_between", "Cefepime_between",
  "Ceftriaxone_between", "Ciprofloxacin_between",
  "Metronidazole_between", "Clindamycin_between",
  "Vancomycin_IV_between", "Vancomycin_PO_between"
)

# Check which columns exist
available_cols <- meaningful_abx_cols[meaningful_abx_cols %in% colnames(paired_samples)]

# Calculate if any meaningful antibiotics were given between samples
paired_samples <- paired_samples %>%
  mutate(
    meaningful_abx_between = rowSums(select(., all_of(available_cols)), na.rm = TRUE) > 0
  )

cat("Pairs with NO antibiotics between (original):", sum(!paired_samples$abx_between_any), "\n")
cat("Pairs with NO meaningful abx between (TMP_SMX/Azithromycin excluded):",
    sum(!paired_samples$meaningful_abx_between), "\n\n")

# Get sample 1 antibiotic exposure (7-day window before sample 1)
sample1_abx <- sample_metadata %>%
  select(sample_id, abx_any_7d,
         Pip_Tazo_7d, Meropenem_7d, Cefepime_7d, Metronidazole_7d,
         Vancomycin_IV_7d, Vancomycin_PO_7d, Ciprofloxacin_7d,
         Ceftriaxone_7d, Clindamycin_7d, TMP_SMX_7d) %>%
  rename_with(~ paste0("s1_", .), -sample_id)

paired_with_s1 <- paired_samples %>%
  left_join(sample1_abx, by = c("sample1" = "sample_id"))

# Focus antibiotics (those with meaningful microbiome impact)
focus_antibiotics <- c("Pip_Tazo", "Meropenem", "Cefepime", "Metronidazole")

# Create flag for any focus antibiotic before sample 1
paired_with_s1 <- paired_with_s1 %>%
  mutate(
    s1_any_focus_abx = s1_Pip_Tazo_7d | s1_Meropenem_7d |
                       s1_Cefepime_7d | s1_Metronidazole_7d
  )

# Recovery pairs: focus abx before S1 AND no meaningful abx between AND high-confidence
recovery_pairs <- paired_with_s1 %>%
  filter(
    s1_any_focus_abx == TRUE,
    meaningful_abx_between == FALSE,
    high_confidence == TRUE
  )

cat("=== Recovery Pair Summary ===\n\n")
cat("Focus antibiotics:", paste(focus_antibiotics, collapse = ", "), "\n")
cat("Excluded from between-sample check: TMP_SMX, Azithromycin\n\n")
cat("Total high-confidence recovery pairs:", nrow(recovery_pairs), "\n\n")

# Summary by antibiotic
exposure_summary <- data.frame(

  antibiotic = focus_antibiotics,
  n_pairs = sapply(focus_antibiotics, function(abx) {
    col <- paste0("s1_", abx, "_7d")
    sum(recovery_pairs[[col]], na.rm = TRUE)
  })
)
print(exposure_summary)

cat("\nRecovery interval (days):\n")
cat("  Min:", min(recovery_pairs$interval_days), "\n")
cat("  Median:", median(recovery_pairs$interval_days), "\n")
cat("  Max:", max(recovery_pairs$interval_days), "\n")
cat("  Mean:", round(mean(recovery_pairs$interval_days), 1), "\n")

# Save recovery pairs
saveRDS(recovery_pairs, file.path(recovery_dir, "recovery_pairs.rds"))
write_csv(recovery_pairs, file.path(recovery_dir, "recovery_pairs.csv"))

# =============================================================================
# 2. Calculate Genus-Level Recovery Metrics
# =============================================================================

cat("\n=== Calculating Recovery Metrics ===\n\n")

# Filter genus matrix to recovery samples
recovery_samples <- unique(c(recovery_pairs$sample1, recovery_pairs$sample2))
genus_matrix_filt <- genus_matrix[rownames(genus_matrix) %in% recovery_samples, ]

# Convert to relative abundance
genus_rel <- genus_matrix_filt / rowSums(genus_matrix_filt)

# Function to calculate recovery for all pairs
calculate_all_recovery <- function(pairs_df, rel_matrix) {

  all_data <- list()

  for (i in 1:nrow(pairs_df)) {
    s1 <- pairs_df$sample1[i]
    s2 <- pairs_df$sample2[i]

    if (!s1 %in% rownames(rel_matrix) || !s2 %in% rownames(rel_matrix)) next

    abund1 <- rel_matrix[s1, ]
    abund2 <- rel_matrix[s2, ]

    # Log2 fold change with pseudocount
    log2fc <- log2((abund2 + 1e-6) / (abund1 + 1e-6))

    all_data[[i]] <- data.frame(
      PairNumber = pairs_df$PairNumber[i],
      MRN = pairs_df$MRN[i],
      PatientGroup = pairs_df$PatientGroup[i],
      interval_days = pairs_df$interval_days[i],
      s1_Pip_Tazo = pairs_df$s1_Pip_Tazo_7d[i],
      s1_Meropenem = pairs_df$s1_Meropenem_7d[i],
      s1_Cefepime = pairs_df$s1_Cefepime_7d[i],
      s1_Metronidazole = pairs_df$s1_Metronidazole_7d[i],
      genus = names(log2fc),
      abund_before = as.numeric(abund1),
      abund_after = as.numeric(abund2),
      log2fc = as.numeric(log2fc)
    )
  }

  bind_rows(all_data)
}

recovery_df <- calculate_all_recovery(recovery_pairs, genus_rel)
cat("Total genus-pair observations:", nrow(recovery_df), "\n")

# =============================================================================
# 3. Pooled Recovery Analysis
# =============================================================================

cat("\n=== Pooled Recovery Analysis (All Antibiotics) ===\n\n")

# Summarize by genus across all pairs
pooled_recovery <- recovery_df %>%
  group_by(genus) %>%
  summarise(
    n_pairs = n(),
    mean_abund_before = mean(abund_before),
    mean_abund_after = mean(abund_after),
    mean_log2fc = mean(log2fc),
    sd_log2fc = sd(log2fc),
    se_log2fc = sd(log2fc) / sqrt(n()),
    .groups = "drop"
  ) %>%
  # One-sample t-test: is mean log2fc different from 0?
  mutate(
    t_stat = mean_log2fc / se_log2fc,
    p_value = 2 * pt(-abs(t_stat), df = n_pairs - 1)
  ) %>%
  # Filter to genera with reasonable abundance
  filter(mean_abund_before > 0.001 | mean_abund_after > 0.001) %>%
  arrange(p_value)

# Adjust p-values
pooled_recovery$p_adj <- p.adjust(pooled_recovery$p_value, method = "BH")
pooled_recovery$antibiotic <- "Pooled"

cat("Genera with sufficient abundance:", nrow(pooled_recovery), "\n")
cat("Significant recoveries (p < 0.05):", sum(pooled_recovery$p_value < 0.05 & pooled_recovery$mean_log2fc > 0), "\n")
cat("Significant declines (p < 0.05):", sum(pooled_recovery$p_value < 0.05 & pooled_recovery$mean_log2fc < 0), "\n\n")

# Top recovering genera
cat("Top 15 RECOVERING genera (pooled):\n")
pooled_recovery %>%
  filter(mean_log2fc > 0) %>%
  head(15) %>%
  select(genus, mean_abund_before, mean_abund_after, mean_log2fc, p_value, p_adj) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print()

# Save pooled results
write_csv(pooled_recovery, file.path(recovery_dir, "recovery_pooled.csv"))

# =============================================================================
# 4. Recovery by Individual Antibiotic
# =============================================================================

cat("\n=== Recovery by Individual Antibiotic ===\n\n")

analyze_by_antibiotic <- function(abx_name, recovery_data) {

  abx_col <- paste0("s1_", abx_name)

  # Filter to pairs with this antibiotic before sample 1
  abx_data <- recovery_data %>%
    filter(.data[[abx_col]] == TRUE)

  n_pairs <- length(unique(abx_data$PairNumber))

  if (n_pairs < 3) {
    cat("Skipping", abx_name, "- only", n_pairs, "pairs\n")
    return(NULL)
  }

  cat(abx_name, ": n =", n_pairs, "pairs\n")

  # Summarize by genus
  genus_recovery <- abx_data %>%
    group_by(genus) %>%
    summarise(
      n_pairs = n(),
      mean_abund_before = mean(abund_before),
      mean_abund_after = mean(abund_after),
      mean_log2fc = mean(log2fc),
      sd_log2fc = sd(log2fc),
      se_log2fc = sd(log2fc) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      t_stat = mean_log2fc / se_log2fc,
      p_value = 2 * pt(-abs(t_stat), df = n_pairs - 1)
    ) %>%
    filter(mean_abund_before > 0.001 | mean_abund_after > 0.001) %>%
    arrange(p_value)

  genus_recovery$p_adj <- p.adjust(genus_recovery$p_value, method = "BH")
  genus_recovery$antibiotic <- abx_name

  return(genus_recovery)
}

# Analyze each focus antibiotic
abx_results <- list()
for (abx in focus_antibiotics) {
  result <- analyze_by_antibiotic(abx, recovery_df)
  if (!is.null(result)) {
    abx_results[[abx]] <- result
  }
}

# Combine all results
all_abx_recovery <- bind_rows(abx_results)

# Save per-antibiotic results
write_csv(all_abx_recovery, file.path(recovery_dir, "recovery_by_antibiotic.csv"))

# Print significant findings for each antibiotic
cat("\n=== Significant Recoveries by Antibiotic (p < 0.05) ===\n\n")

sig_by_abx <- all_abx_recovery %>%
  filter(p_value < 0.05, mean_log2fc > 0) %>%
  arrange(antibiotic, p_value) %>%
  select(antibiotic, genus, mean_log2fc, p_value)

print(sig_by_abx)

# =============================================================================
# 5. Functional Group-Level Analysis (Consistent with Persistence Script)
# =============================================================================

cat("\n=== Functional Group-Level Recovery ===\n\n")

# Define functional groups (same as persistence script)
obligate_anaerobes <- c(
  "Bacteroides", "Parabacteroides", "Prevotella", "Clostridium",
  "Faecalibacterium", "Blautia", "Roseburia", "Ruminococcus",
  "Alistipes", "Coprococcus", "Dorea", "Lachnospira", "Eubacterium"
)

enterobacteriaceae <- c(
  "Escherichia", "Klebsiella", "Enterobacter", "Citrobacter",
  "Serratia", "Proteus", "Salmonella", "Shigella"
)

# Calculate functional group abundances for each sample
# Use the FULL genus matrix (not filtered) to get all samples
genus_rel_full <- genus_matrix / rowSums(genus_matrix)

calc_group_abundance <- function(rel_matrix, genera) {
  matching <- genera[genera %in% colnames(rel_matrix)]
  if (length(matching) == 0) return(rep(0, nrow(rel_matrix)))
  if (length(matching) == 1) return(rel_matrix[, matching])
  rowSums(rel_matrix[, matching])
}

sample_groups <- data.frame(
  sample_id = rownames(genus_rel_full),
  anaerobes = calc_group_abundance(genus_rel_full, obligate_anaerobes),
  enterobact = calc_group_abundance(genus_rel_full, enterobacteriaceae),
  enterococcus = if ("Enterococcus" %in% colnames(genus_rel_full))
                   genus_rel_full[, "Enterococcus"] else 0,
  stringsAsFactors = FALSE
)

# Merge with recovery pairs
s1_groups <- sample_groups %>%
  rename(sample1 = sample_id,
         anaerobes_S1 = anaerobes,
         enterobact_S1 = enterobact,
         enterococcus_S1 = enterococcus)

s2_groups <- sample_groups %>%
  rename(sample2 = sample_id,
         anaerobes_S2 = anaerobes,
         enterobact_S2 = enterobact,
         enterococcus_S2 = enterococcus)

recovery_groups <- recovery_pairs %>%
  left_join(s1_groups, by = "sample1") %>%
  left_join(s2_groups, by = "sample2") %>%
  filter(!is.na(anaerobes_S1), !is.na(anaerobes_S2))

cat("Recovery pairs with complete data:", nrow(recovery_groups), "\n\n")

# Calculate summary statistics (same format as persistence script)
cat("Functional Group Changes During Recovery (S1 on Abx -> S2 off Abx):\n")
cat("-" %>% rep(70) %>% paste(collapse = ""), "\n\n")
cat("                    S1 (on Abx)  S2 (off Abx)  Change       Retention    p-value\n")

group_summary <- data.frame()

for (grp in c("anaerobes", "enterobact", "enterococcus")) {
  s1_col <- paste0(grp, "_S1")
  s2_col <- paste0(grp, "_S2")

  s1_mean <- mean(recovery_groups[[s1_col]])
  s2_mean <- mean(recovery_groups[[s2_col]])
  change <- s2_mean - s1_mean
  retention <- (s2_mean / s1_mean) * 100

  # Paired t-test for the difference
  diff <- recovery_groups[[s2_col]] - recovery_groups[[s1_col]]
  test <- t.test(diff)

  cat(sprintf("  %-12s    %5.1f%%       %5.1f%%       %+5.1f%%       %5.0f%%      %.4f\n",
      grp, s1_mean * 100, s2_mean * 100, change * 100, retention, test$p.value))

  group_summary <- rbind(group_summary, data.frame(
    group = grp,
    n_pairs = nrow(recovery_groups),
    mean_S1 = s1_mean,
    mean_S2 = s2_mean,
    change = change,
    retention_pct = retention,
    p_value = test$p.value
  ))
}

# Save functional group summary
write_csv(group_summary, file.path(recovery_dir, "recovery_functional_groups.csv"))

cat("\n\nInterpretation:\n")
cat("  - Positive change = group INCREASES after stopping antibiotics (recovery)\n")
cat("  - Negative change = group DECREASES after stopping antibiotics\n")
cat("  - Retention > 100% = group is higher in S2 than S1\n")

# =============================================================================
# 5b. Genus-Level Detail Within Functional Groups
# =============================================================================

cat("\n\n=== Genus-Level Detail Within Functional Groups ===\n\n")

# Combine pooled and per-antibiotic for functional group summary
all_recovery <- bind_rows(pooled_recovery, all_abx_recovery)

cat("Obligate Anaerobes (pooled):\n")
pooled_recovery %>%
  filter(genus %in% obligate_anaerobes) %>%
  arrange(p_value) %>%
  select(genus, mean_abund_before, mean_abund_after, mean_log2fc, p_value) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print()

cat("\nEnterobacteriaceae (pooled):\n")
pooled_recovery %>%
  filter(genus %in% enterobacteriaceae) %>%
  arrange(p_value) %>%
  select(genus, mean_abund_before, mean_abund_after, mean_log2fc, p_value) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print()

cat("\nEnterococcus (pooled):\n")
pooled_recovery %>%
  filter(genus == "Enterococcus") %>%
  select(genus, mean_abund_before, mean_abund_after, mean_log2fc, p_value) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print()

# =============================================================================
# 5c. Individual Antibiotic Analysis (Functional Groups)
# =============================================================================

cat("\n\n=== Individual Antibiotic Analysis (Functional Groups) ===\n\n")

# Identify single-antibiotic recovery pairs
recovery_groups <- recovery_groups %>%
  mutate(
    n_focus_abx = (s1_Pip_Tazo_7d == TRUE) + (s1_Meropenem_7d == TRUE) +
                  (s1_Cefepime_7d == TRUE) + (s1_Metronidazole_7d == TRUE)
  )

single_abx_pairs <- recovery_groups %>% filter(n_focus_abx == 1)

cat("Single-antibiotic recovery pairs:", nrow(single_abx_pairs), "of", nrow(recovery_groups),
    "(", round(nrow(single_abx_pairs)/nrow(recovery_groups)*100, 1), "%)\n\n")

# Analyze each antibiotic separately
individual_abx_results <- list()

for (abx in c("Pip_Tazo", "Meropenem", "Cefepime", "Metronidazole")) {
  abx_col <- paste0("s1_", abx, "_7d")

  # Filter to single-antibiotic pairs with this specific antibiotic
  abx_pairs <- single_abx_pairs %>% filter(.data[[abx_col]] == TRUE)
  n_pairs <- nrow(abx_pairs)

  if (n_pairs < 3) {
    cat(abx, ": n =", n_pairs, "pairs (too few for analysis)\n\n")
    next
  }

  cat(abx, ": n =", n_pairs, "single-antibiotic pairs\n")
  cat("-" %>% rep(50) %>% paste(collapse = ""), "\n")
  cat("                    S1 (on Abx)  S2 (off Abx)  Change       Retention    p-value\n")

  for (grp in c("anaerobes", "enterobact", "enterococcus")) {
    s1_col <- paste0(grp, "_S1")
    s2_col <- paste0(grp, "_S2")

    s1_mean <- mean(abx_pairs[[s1_col]])
    s2_mean <- mean(abx_pairs[[s2_col]])
    change <- s2_mean - s1_mean
    retention <- (s2_mean / s1_mean) * 100

    diff <- abx_pairs[[s2_col]] - abx_pairs[[s1_col]]
    test <- t.test(diff)

    cat(sprintf("  %-12s    %5.1f%%       %5.1f%%       %+5.1f%%       %5.0f%%      %.4f\n",
        grp, s1_mean * 100, s2_mean * 100, change * 100, retention, test$p.value))

    individual_abx_results[[paste(abx, grp, sep = "_")]] <- data.frame(
      antibiotic = abx,
      group = grp,
      n_pairs = n_pairs,
      mean_S1 = s1_mean,
      mean_S2 = s2_mean,
      change = change,
      retention_pct = retention,
      p_value = test$p.value
    )
  }
  cat("\n")
}

# Save individual antibiotic results
if (length(individual_abx_results) > 0) {
  individual_abx_df <- bind_rows(individual_abx_results)
  write_csv(individual_abx_df, file.path(recovery_dir, "recovery_by_individual_antibiotic.csv"))
  cat("Saved: recovery_by_individual_antibiotic.csv\n")
}

# =============================================================================
# 6. Generate Visualizations
# =============================================================================

cat("\n=== Generating Visualizations ===\n\n")

# 6.1 Heatmap of key genera recovery by antibiotic
key_genera <- c(
  "Blautia", "Roseburia", "Clostridium", "Coprococcus", "Ruminococcus",
  "Bacteroides", "Parabacteroides", "Faecalibacterium",
  "Veillonella", "Clostridioides",
  "Enterococcus", "Lactobacillus",
  "Escherichia", "Klebsiella", "Enterobacter", "Citrobacter"
)

heatmap_data <- all_recovery %>%
  filter(genus %in% key_genera) %>%
  select(antibiotic, genus, mean_log2fc, p_value) %>%
  mutate(
    sig_label = case_when(
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ ".",
      TRUE ~ ""
    ),
    genus = factor(genus, levels = rev(key_genera)),
    antibiotic = factor(antibiotic,
      levels = c("Pooled", "Pip_Tazo", "Cefepime", "Meropenem", "Metronidazole"))
  )

p_heatmap <- ggplot(heatmap_data, aes(x = antibiotic, y = genus, fill = mean_log2fc)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_label), size = 4, color = "black") +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0,
    limits = c(-3, 5),
    oob = scales::squish,
    name = "Log2FC\n(Recovery)"
  ) +
  labs(
    title = "Genus Recovery After Antibiotic Cessation",
    subtitle = "Positive = genus increases after stopping antibiotics (* p<0.05, ** p<0.01, . p<0.1)",
    x = "Antibiotic Before Sample 1",
    y = NULL,
    caption = "Recovery pairs: Abx before S1, no meaningful Abx between (TMP_SMX/Azithromycin excluded)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(figures_dir, "recovery_heatmap.pdf"), p_heatmap,
       width = 9, height = 8)
cat("Saved: figures/recovery_heatmap.pdf\n")

# 6.2 Forest plot of pooled recovery (top genera)
top_pooled <- pooled_recovery %>%
  filter(mean_abund_before > 0.005 | mean_abund_after > 0.005) %>%
  arrange(desc(abs(mean_log2fc))) %>%
  head(25) %>%
  mutate(
    genus = factor(genus, levels = rev(genus)),
    significant = p_value < 0.05,
    direction = ifelse(mean_log2fc > 0, "Recovery", "Decline")
  )

p_forest <- ggplot(top_pooled, aes(x = mean_log2fc, y = genus)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(
    aes(xmin = mean_log2fc - 1.96 * se_log2fc,
        xmax = mean_log2fc + 1.96 * se_log2fc),
    height = 0.3, color = "gray40"
  ) +
  geom_point(aes(color = direction, shape = significant), size = 3) +
  scale_color_manual(values = c("Recovery" = "#b2182b", "Decline" = "#2166ac")) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
  labs(
    title = "Pooled Microbiome Recovery After Antibiotic Cessation",
    subtitle = paste("n =", nrow(recovery_pairs), "recovery pairs"),
    x = "Mean Log2 Fold Change (95% CI)",
    y = NULL,
    color = "Direction",
    shape = "p < 0.05"
  ) +
  theme_bw() +
  theme(legend.position = "right")

ggsave(file.path(figures_dir, "recovery_forest_pooled.pdf"), p_forest,
       width = 10, height = 8)
cat("Saved: figures/recovery_forest_pooled.pdf\n")

# 6.3 Bar plot comparing recovery across antibiotics for key genera
comparison_genera <- c("Blautia", "Roseburia", "Bacteroides", "Enterococcus",
                       "Escherichia", "Klebsiella", "Veillonella")

comparison_data <- all_abx_recovery %>%
  filter(genus %in% comparison_genera) %>%
  mutate(
    genus = factor(genus, levels = comparison_genera),
    sig_label = ifelse(p_value < 0.05, "*", "")
  )

p_comparison <- ggplot(comparison_data, aes(x = genus, y = mean_log2fc, fill = antibiotic)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_text(aes(label = sig_label, y = mean_log2fc + sign(mean_log2fc) * 0.3),
            position = position_dodge(width = 0.8), size = 5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Recovery of Key Genera by Prior Antibiotic",
    subtitle = "Log2FC from sample 1 (on Abx) to sample 2 (off Abx). * p < 0.05",
    x = NULL,
    y = "Mean Log2 Fold Change",
    fill = "Prior Antibiotic"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    legend.position = "right"
  )

ggsave(file.path(figures_dir, "recovery_comparison_by_abx.pdf"), p_comparison,
       width = 10, height = 6)
cat("Saved: figures/recovery_comparison_by_abx.pdf\n")

# 6.4 Recovery vs interval days (does longer time = more recovery?)
interval_analysis <- recovery_df %>%
  filter(genus %in% c("Blautia", "Roseburia", "Bacteroides", "Enterococcus")) %>%
  group_by(PairNumber, interval_days, genus) %>%
  summarise(log2fc = mean(log2fc), .groups = "drop")

p_interval <- ggplot(interval_analysis, aes(x = interval_days, y = log2fc)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ genus, scales = "free_y") +
  labs(
    title = "Recovery Over Time: Does Longer Interval = More Recovery?",
    x = "Days Between Samples",
    y = "Log2 Fold Change"
  ) +
  theme_bw()

ggsave(file.path(figures_dir, "recovery_vs_interval.pdf"), p_interval,
       width = 10, height = 8)
cat("Saved: figures/recovery_vs_interval.pdf\n")

# =============================================================================
# 7. Summary Statistics
# =============================================================================

cat("\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("SUMMARY\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")

cat("Recovery pairs analyzed:", nrow(recovery_pairs), "\n")
cat("  Pip_Tazo:", sum(recovery_pairs$s1_Pip_Tazo_7d), "pairs\n")
cat("  Cefepime:", sum(recovery_pairs$s1_Cefepime_7d), "pairs\n")
cat("  Meropenem:", sum(recovery_pairs$s1_Meropenem_7d), "pairs\n")
cat("  Metronidazole:", sum(recovery_pairs$s1_Metronidazole_7d), "pairs\n\n")

cat("Median recovery interval:", median(recovery_pairs$interval_days), "days\n\n")

cat("Significant findings (p < 0.05):\n")
cat("-" %>% rep(40) %>% paste(collapse = ""), "\n")

all_sig <- bind_rows(
  pooled_recovery %>% filter(p_value < 0.05) %>% mutate(antibiotic = "Pooled"),
  all_abx_recovery %>% filter(p_value < 0.05)
) %>%
  arrange(antibiotic, p_value) %>%
  select(antibiotic, genus, mean_log2fc, p_value)

for (i in 1:min(15, nrow(all_sig))) {
  direction <- ifelse(all_sig$mean_log2fc[i] > 0, "recovers", "declines")
  cat(sprintf("  %s: %s %s (log2FC=%.2f, p=%.4f)\n",
              all_sig$antibiotic[i],
              all_sig$genus[i],
              direction,
              all_sig$mean_log2fc[i],
              all_sig$p_value[i]))
}

cat("\n=== Key Biological Insights ===\n\n")
cat("1. ANAEROBES recover after cessation of anti-anaerobic antibiotics\n")
cat("   - Roseburia, Coprococcus, Blautia, Eubacterium increase after Pip_Tazo\n")
cat("   - Blautia shows strongest recovery after Meropenem\n\n")

cat("2. VEILLONELLA consistently recovers across all antibiotics\n\n")

cat("3. ENTEROCOCCUS declines after Meropenem/Metronidazole cessation\n")
cat("   - These antibiotics select for Enterococcus during exposure\n")
cat("   - After cessation, Enterococcus returns toward baseline\n\n")

cat("4. CLOSTRIDIOIDES increases in recovery window\n")
cat("   - Consistent with known C. diff infection risk period\n\n")

# =============================================================================
# 8. Save Output Files
# =============================================================================

cat("=== Output Files ===\n\n")

outputs <- data.frame(
  file = c(
    "recovery_pairs.csv",
    "recovery_pooled.csv",
    "recovery_by_antibiotic.csv",
    "figures/recovery_heatmap.pdf",
    "figures/recovery_forest_pooled.pdf",
    "figures/recovery_comparison_by_abx.pdf",
    "figures/recovery_vs_interval.pdf"
  ),
  description = c(
    "Recovery pair metadata and antibiotic exposures",
    "Pooled recovery analysis (all antibiotics)",
    "Recovery analysis stratified by antibiotic",
    "Heatmap of genus recovery by antibiotic",
    "Forest plot of pooled recovery effects",
    "Bar plot comparing key genera across antibiotics",
    "Scatterplot of recovery vs time interval"
  )
)

print(outputs)

cat("\nRecovery analysis complete!\n")
