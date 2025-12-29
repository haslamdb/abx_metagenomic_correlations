#!/usr/bin/env Rscript
# =============================================================================
# 08_persistence_analysis.R
# Analysis of microbial persistence during antibiotic exposure
# Tests the hypothesis that Enterobacteriaceae/Enterococcus persist better
# than obligate anaerobes during antibiotic courses
# =============================================================================

library(tidyverse)
library(vegan)

# Set paths
project_dir <- "/home/david/projects/abx_metagenomic_correlations"
results_dir <- file.path(project_dir, "results/new_analysis_legacy_data")
figures_dir <- file.path(results_dir, "figures")
persist_dir <- file.path(results_dir, "persistence_analysis")

dir.create(persist_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
load(file.path(project_dir, "data/bracken_count_matrices.RData"))
load(file.path(project_dir, "data/prepared_data.RData"))

paired_samples <- prepared_data$paired_samples
sample_metadata <- prepared_data$sample_metadata
genus_matrix <- bracken_data$genus_matrix
species_matrix <- bracken_data$species_matrix

cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("PERSISTENCE ANALYSIS: Why Enterobacteriaceae/Enterococcus Dominate\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n\n")

# =============================================================================
# 1. Define Functional Groups with Individual Genera
# =============================================================================

cat("=== Functional Group Definitions ===\n\n")

# Obligate anaerobes - require anaerobic conditions, slow-growing
obligate_anaerobes <- c(
  "Bacteroides", "Parabacteroides", "Prevotella", "Alistipes",       # Bacteroidetes
  "Faecalibacterium", "Blautia", "Roseburia", "Ruminococcus",        # Firmicutes - Clostridia
  "Coprococcus", "Dorea", "Lachnospira", "Eubacterium", "Anaerostipes",
  "Clostridium", "Lachnoclostridium", "Clostridioides"
)

# Enterobacteriaceae - facultative anaerobes, stress-resistant
enterobacteriaceae <- c(
  "Escherichia", "Klebsiella", "Enterobacter", "Citrobacter",
  "Serratia", "Proteus", "Salmonella", "Shigella", "Cronobacter",
  "Raoultella", "Hafnia"
)

# Enterococcus - intrinsically resistant to many antibiotics
enterococcus <- c("Enterococcus")

# Other potential opportunists
other_opportunists <- c(
  "Staphylococcus", "Streptococcus", "Lactobacillus",
  "Pseudomonas", "Acinetobacter", "Stenotrophomonas"
)

# Print which genera we have in each group
cat("Checking available genera in each group:\n")
available_genera <- colnames(genus_matrix)

for (grp_name in c("obligate_anaerobes", "enterobacteriaceae", "enterococcus")) {
  grp <- get(grp_name)
  matches <- grp[grp %in% available_genera]
  cat(sprintf("  %s: %d/%d available\n", grp_name, length(matches), length(grp)))
}

# =============================================================================
# 2. Calculate Group and Genus Abundances for Each Sample
# =============================================================================

cat("\n=== Calculating Abundances ===\n\n")

# Convert to relative abundance
genus_rel <- genus_matrix / rowSums(genus_matrix)

# Function to calculate group abundance
calc_group_abundance <- function(rel_matrix, genera) {
  matching <- genera[genera %in% colnames(rel_matrix)]
  if (length(matching) == 0) return(rep(0, nrow(rel_matrix)))
  if (length(matching) == 1) return(rel_matrix[, matching])
  rowSums(rel_matrix[, matching])
}

# Calculate abundances for all samples
sample_abundances <- data.frame(
  sample_id = rownames(genus_rel),
  anaerobes = calc_group_abundance(genus_rel, obligate_anaerobes),
  enterobact = calc_group_abundance(genus_rel, enterobacteriaceae),
  enterococcus = calc_group_abundance(genus_rel, enterococcus),
  stringsAsFactors = FALSE
)

# Add individual genera abundances
key_genera <- c(
  # Obligate anaerobes
  "Bacteroides", "Parabacteroides", "Blautia", "Roseburia", "Faecalibacterium",
  "Ruminococcus", "Coprococcus", "Clostridium", "Alistipes",
  # Enterobacteriaceae
  "Escherichia", "Klebsiella", "Enterobacter", "Citrobacter",
  # Enterococcus
  "Enterococcus"
)

for (genus in key_genera) {
  if (genus %in% colnames(genus_rel)) {
    sample_abundances[[genus]] <- genus_rel[, genus]
  }
}

# Merge with metadata
sample_abundances <- sample_abundances %>%
  left_join(
    sample_metadata %>%
      select(sample_id, MRN, PatientGroup, abx_any_7d, abx_any_14d),
    by = "sample_id"
  )

# =============================================================================
# 3. Sample Exposure Status (for reference only)
# =============================================================================

# Define exposure status for sample-level data (used in visualizations)
sample_abundances <- sample_abundances %>%
  mutate(
    exposure_status = case_when(
      abx_any_14d == FALSE ~ "No_Abx_14d",
      abx_any_7d == TRUE ~ "Abx_7d",
      TRUE ~ "Abx_8_14d"
    )
  )

cat("=== Sample Distribution ===\n\n")
cat("Sample distribution by exposure status:\n")
print(table(sample_abundances$exposure_status))
cat("\nNote: Primary analysis uses PAIRED samples (Section 4), not cross-sectional.\n")

# =============================================================================
# 4. Paired Analysis: Persistence During Antibiotic Exposure
# =============================================================================

cat("\n=== Paired Analysis: Persistence During Antibiotics ===\n\n")

# Get pairs where antibiotics were given BETWEEN samples
abx_pairs <- paired_samples %>%
  filter(abx_between_any == TRUE, high_confidence == TRUE)

cat("Pairs with antibiotics between samples:", nrow(abx_pairs), "\n")

# Calculate abundances for paired samples
calc_paired_abundances <- function(pairs_df, abundances_df) {
  results <- list()

  for (i in 1:nrow(pairs_df)) {
    s1 <- pairs_df$sample1[i]
    s2 <- pairs_df$sample2[i]

    a1 <- abundances_df %>% filter(sample_id == s1)
    a2 <- abundances_df %>% filter(sample_id == s2)

    if (nrow(a1) == 0 || nrow(a2) == 0) next

    results[[i]] <- data.frame(
      PairNumber = pairs_df$PairNumber[i],
      MRN = pairs_df$MRN[i],
      interval_days = pairs_df$interval_days[i],
      # Group abundances before/after
      anaerobes_before = a1$anaerobes,
      anaerobes_after = a2$anaerobes,
      enterobact_before = a1$enterobact,
      enterobact_after = a2$enterobact,
      enterococcus_before = a1$enterococcus,
      enterococcus_after = a2$enterococcus,
      stringsAsFactors = FALSE
    )
  }

  bind_rows(results)
}

paired_abund <- calc_paired_abundances(abx_pairs, sample_abundances)
cat("Pairs with complete data:", nrow(paired_abund), "\n")

# Calculate persistence metrics
paired_abund <- paired_abund %>%
  mutate(
    # Log2 fold change during antibiotics
    anaerobes_log2fc = log2((anaerobes_after + 1e-6) / (anaerobes_before + 1e-6)),
    enterobact_log2fc = log2((enterobact_after + 1e-6) / (enterobact_before + 1e-6)),
    enterococcus_log2fc = log2((enterococcus_after + 1e-6) / (enterococcus_before + 1e-6)),
    # Absolute change
    anaerobes_delta = anaerobes_after - anaerobes_before,
    enterobact_delta = enterobact_after - enterobact_before,
    enterococcus_delta = enterococcus_after - enterococcus_before,
    # Persistence: still >1% of original abundance?
    anaerobes_persists = (anaerobes_after / (anaerobes_before + 1e-6)) > 0.01,
    enterobact_persists = (enterobact_after / (enterobact_before + 1e-6)) > 0.01,
    enterococcus_persists = (enterococcus_after / (enterococcus_before + 1e-6)) > 0.01,
    # Maintained >50% of original?
    anaerobes_maintained = (anaerobes_after / (anaerobes_before + 1e-6)) > 0.5,
    enterobact_maintained = (enterobact_after / (enterobact_before + 1e-6)) > 0.5,
    enterococcus_maintained = (enterococcus_after / (enterococcus_before + 1e-6)) > 0.5
  )

# Summarize persistence rates
cat("\n=== Persistence During Antibiotic Exposure ===\n")
cat("-" %>% rep(50) %>% paste(collapse = ""), "\n\n")

# Calculate summary statistics in consistent format with recovery script
cat("Functional Group Changes During Antibiotics (S1 pre-Abx -> S2 post-Abx):\n")
cat("-" %>% rep(70) %>% paste(collapse = ""), "\n\n")
cat("                    S1 (pre-Abx)  S2 (post-Abx)  Change       Retention    p-value\n")

persistence_summary <- data.frame()

for (grp in c("anaerobes", "enterobact", "enterococcus")) {
  before_col <- paste0(grp, "_before")
  after_col <- paste0(grp, "_after")

  s1_mean <- mean(paired_abund[[before_col]])
  s2_mean <- mean(paired_abund[[after_col]])
  change <- s2_mean - s1_mean
  retention <- (s2_mean / s1_mean) * 100

  # Paired t-test for the difference (consistent with recovery script)
  diff <- paired_abund[[after_col]] - paired_abund[[before_col]]
  test <- t.test(diff)

  cat(sprintf("  %-12s    %5.1f%%       %5.1f%%       %+5.1f%%       %5.0f%%      %.4f\n",
      grp, s1_mean * 100, s2_mean * 100, change * 100, retention, test$p.value))

  persistence_summary <- rbind(persistence_summary, data.frame(
    group = grp,
    n_pairs = nrow(paired_abund),
    mean_S1 = s1_mean,
    mean_S2 = s2_mean,
    change = change,
    retention_pct = retention,
    p_value = test$p.value
  ))
}

cat("\n\nInterpretation:\n")
cat("  - Positive change = group INCREASES during antibiotics\n")
cat("  - Negative change = group DECREASES during antibiotics (depletion)\n")
cat("  - Retention > 100% = group is higher after antibiotics than before\n")

# Save paired results
write_csv(paired_abund, file.path(persist_dir, "paired_persistence_data.csv"))

# =============================================================================
# 4b. Individual Antibiotic Analysis (Functional Groups)
# =============================================================================

cat("\n\n=== Individual Antibiotic Analysis (Functional Groups) ===\n\n")

# Get antibiotic exposure for each pair
abx_between_cols <- c("Pip_Tazo_between", "Meropenem_between", "Cefepime_between",
                      "Ceftriaxone_between", "Metronidazole_between", "Vancomycin_IV_between")
avail_abx_cols <- abx_between_cols[abx_between_cols %in% colnames(abx_pairs)]

# Count antibiotics per pair
abx_pairs <- abx_pairs %>%
  mutate(n_abx_between = rowSums(select(., all_of(avail_abx_cols)), na.rm = TRUE))

cat("Antibiotic exposure distribution:\n")
cat("  Single antibiotic:", sum(abx_pairs$n_abx_between == 1), "pairs (",
    round(sum(abx_pairs$n_abx_between == 1)/nrow(abx_pairs)*100, 1), "%)\n")
cat("  Multiple antibiotics:", sum(abx_pairs$n_abx_between > 1), "pairs (",
    round(sum(abx_pairs$n_abx_between > 1)/nrow(abx_pairs)*100, 1), "%)\n\n")

cat("Note: Most pairs have multiple antibiotics. Analysis below shows pairs\n")
cat("      that INCLUDE each antibiotic (may overlap with other antibiotics).\n\n")

# Merge functional group abundances with abx_pairs
s1_groups <- data.frame(
  sample1 = rownames(genus_rel),
  anaerobes_before = calc_group_abundance(genus_rel, obligate_anaerobes),
  enterobact_before = calc_group_abundance(genus_rel, enterobacteriaceae),
  enterococcus_before = calc_group_abundance(genus_rel, enterococcus),
  stringsAsFactors = FALSE
)

s2_groups <- data.frame(
  sample2 = rownames(genus_rel),
  anaerobes_after = calc_group_abundance(genus_rel, obligate_anaerobes),
  enterobact_after = calc_group_abundance(genus_rel, enterobacteriaceae),
  enterococcus_after = calc_group_abundance(genus_rel, enterococcus),
  stringsAsFactors = FALSE
)

abx_pairs_groups <- abx_pairs %>%
  left_join(s1_groups, by = "sample1") %>%
  left_join(s2_groups, by = "sample2") %>%
  filter(!is.na(anaerobes_before), !is.na(anaerobes_after))

# Analyze each antibiotic (includes analysis)
includes_abx_results <- list()

focus_antibiotics <- c("Pip_Tazo", "Meropenem", "Cefepime", "Ceftriaxone",
                       "Metronidazole", "Vancomycin_IV")

for (abx in focus_antibiotics) {
  abx_col <- paste0(abx, "_between")
  if (!abx_col %in% colnames(abx_pairs_groups)) next

  # Filter to pairs that include this antibiotic
  abx_subset <- abx_pairs_groups %>% filter(.data[[abx_col]] == TRUE)
  n_pairs <- nrow(abx_subset)

  if (n_pairs < 5) {
    cat(abx, ": n =", n_pairs, "pairs (too few for analysis)\n\n")
    next
  }

  # Calculate % that are single vs multiple
  n_single <- sum(abx_subset$n_abx_between == 1)
  pct_single <- round(n_single / n_pairs * 100, 0)

  cat(abx, ": n =", n_pairs, "pairs (", pct_single, "% monotherapy)\n")
  cat("-" %>% rep(55) %>% paste(collapse = ""), "\n")
  cat("                    S1 (pre-Abx)  S2 (post-Abx)  Change       Retention    p-value\n")

  for (grp in c("anaerobes", "enterobact", "enterococcus")) {
    before_col <- paste0(grp, "_before")
    after_col <- paste0(grp, "_after")

    s1_mean <- mean(abx_subset[[before_col]])
    s2_mean <- mean(abx_subset[[after_col]])
    change <- s2_mean - s1_mean
    retention <- (s2_mean / s1_mean) * 100

    diff <- abx_subset[[after_col]] - abx_subset[[before_col]]
    test <- t.test(diff)

    cat(sprintf("  %-12s    %5.1f%%       %5.1f%%       %+5.1f%%       %5.0f%%      %.4f\n",
        grp, s1_mean * 100, s2_mean * 100, change * 100, retention, test$p.value))

    includes_abx_results[[paste(abx, grp, sep = "_")]] <- data.frame(
      antibiotic = abx,
      analysis_type = "includes",
      group = grp,
      n_pairs = n_pairs,
      pct_monotherapy = pct_single,
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
if (length(includes_abx_results) > 0) {
  includes_abx_df <- bind_rows(includes_abx_results)
  write_csv(includes_abx_df, file.path(persist_dir, "persistence_by_antibiotic.csv"))
  cat("Saved: persistence_by_antibiotic.csv\n")
}

# =============================================================================
# 5. Genus-Level Breakdown of Persistence
# =============================================================================

cat("\n=== Genus-Level Persistence Breakdown ===\n\n")

# Calculate genus-level paired abundances
calc_genus_paired <- function(pairs_df, genus_rel_matrix) {
  results <- list()

  for (i in 1:nrow(pairs_df)) {
    s1 <- pairs_df$sample1[i]
    s2 <- pairs_df$sample2[i]

    if (!s1 %in% rownames(genus_rel_matrix) || !s2 %in% rownames(genus_rel_matrix)) next

    abund1 <- genus_rel_matrix[s1, ]
    abund2 <- genus_rel_matrix[s2, ]

    for (genus in key_genera) {
      if (!genus %in% names(abund1)) next

      results[[paste(i, genus)]] <- data.frame(
        PairNumber = pairs_df$PairNumber[i],
        genus = genus,
        abund_before = abund1[[genus]],
        abund_after = abund2[[genus]],
        stringsAsFactors = FALSE
      )
    }
  }

  bind_rows(results)
}

genus_paired <- calc_genus_paired(abx_pairs, genus_rel)

# Calculate log2FC
genus_paired <- genus_paired %>%
  mutate(
    log2fc = log2((abund_after + 1e-6) / (abund_before + 1e-6)),
    retention = (abund_after + 1e-6) / (abund_before + 1e-6)
  )

# Summarize by genus
genus_persistence <- genus_paired %>%
  group_by(genus) %>%
  summarise(
    n_pairs = n(),
    mean_before = mean(abund_before),
    mean_after = mean(abund_after),
    mean_log2fc = mean(log2fc),
    se_log2fc = sd(log2fc) / sqrt(n()),
    mean_retention = mean(retention),
    .groups = "drop"
  ) %>%
  mutate(
    t_stat = mean_log2fc / se_log2fc,
    p_value = 2 * pt(-abs(t_stat), df = n_pairs - 1),
    functional_group = case_when(
      genus %in% obligate_anaerobes ~ "Obligate Anaerobe",
      genus %in% enterobacteriaceae ~ "Enterobacteriaceae",
      genus == "Enterococcus" ~ "Enterococcus",
      TRUE ~ "Other"
    )
  ) %>%
  filter(mean_before > 0.001 | mean_after > 0.001) %>%
  arrange(p_value)

genus_persistence$p_adj <- p.adjust(genus_persistence$p_value, method = "BH")

cat("Genus-level persistence during antibiotics (sorted by p-value):\n")
genus_persistence %>%
  select(genus, functional_group, mean_before, mean_after, mean_log2fc, p_value) %>%
  mutate(
    mean_before = round(mean_before * 100, 2),
    mean_after = round(mean_after * 100, 2),
    mean_log2fc = round(mean_log2fc, 2),
    p_value = round(p_value, 4)
  ) %>%
  print(n = 20)

# Save genus-level results
write_csv(genus_persistence, file.path(persist_dir, "genus_persistence.csv"))

# =============================================================================
# 6. Dominance Analysis
# =============================================================================

cat("\n=== Dominance Patterns ===\n\n")

# Calculate dominance (Enterobact + Enterococcus > 30%)
sample_abundances <- sample_abundances %>%
  mutate(
    opportunist_total = enterobact + enterococcus,
    dominated = opportunist_total > 0.30
  )

dominance_by_exposure <- sample_abundances %>%
  filter(exposure_status %in% c("No_Abx_14d", "Abx_7d")) %>%
  group_by(exposure_status) %>%
  summarise(
    n = n(),
    pct_dominated = mean(dominated, na.rm = TRUE) * 100,
    mean_anaerobes = mean(anaerobes, na.rm = TRUE) * 100,
    mean_enterobact = mean(enterobact, na.rm = TRUE) * 100,
    mean_enterococcus = mean(enterococcus, na.rm = TRUE) * 100,
    mean_opportunist = mean(opportunist_total, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("Dominance by Enterobact/Enterococcus (>30% combined):\n")
print(dominance_by_exposure)

# Chi-squared test for dominance
dom_table <- table(
  sample_abundances$exposure_status,
  sample_abundances$dominated
)[c("No_Abx_14d", "Abx_7d"), ]
chi_test <- chisq.test(dom_table)
cat("\nChi-squared test for dominance pattern:\n")
print(chi_test)

# Note: Recovery analysis is now done in 07_recovery_analysis.R using consistent
# paired sample methodology. See recovery_analysis/recovery_functional_groups.csv
# for results directly comparable to this persistence analysis.

# =============================================================================
# 8. Visualizations
# =============================================================================

cat("\n=== Generating Visualizations ===\n\n")

# 7.1 Stacked bar: Healthy vs Exposed composition
composition_data <- sample_abundances %>%
  filter(exposure_status %in% c("No_Abx_14d", "Abx_7d")) %>%
  group_by(exposure_status) %>%
  summarise(
    `Obligate Anaerobes` = mean(anaerobes),
    `Enterobacteriaceae` = mean(enterobact),
    `Enterococcus` = mean(enterococcus),
    Other = 1 - mean(anaerobes) - mean(enterobact) - mean(enterococcus),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(`Obligate Anaerobes`, Enterobacteriaceae, Enterococcus, Other),
    names_to = "group",
    values_to = "abundance"
  ) %>%
  mutate(
    group = factor(group, levels = c("Other", "Enterococcus", "Enterobacteriaceae", "Obligate Anaerobes")),
    exposure_status = factor(exposure_status,
      levels = c("No_Abx_14d", "Abx_7d"),
      labels = c("No Antibiotics\n(14 days)", "Antibiotics\n(7 days)"))
  )

p_composition <- ggplot(composition_data, aes(x = exposure_status, y = abundance, fill = group)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(
    values = c(
      "Obligate Anaerobes" = "#2166ac",
      "Enterobacteriaceae" = "#b2182b",
      "Enterococcus" = "#d6604d",
      "Other" = "#cccccc"
    )
  ) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    title = "Microbiome Composition: Healthy-like vs Antibiotic-Exposed",
    x = NULL,
    y = "Mean Relative Abundance",
    fill = "Functional Group"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 11)
  )

ggsave(file.path(figures_dir, "persistence_composition.pdf"), p_composition,
       width = 7, height = 6)
cat("Saved: figures/persistence_composition.pdf\n")

# 7.2 Persistence during antibiotics (paired log2FC boxplot)
persistence_plot_data <- paired_abund %>%
  select(PairNumber, anaerobes_log2fc, enterobact_log2fc, enterococcus_log2fc) %>%
  pivot_longer(
    cols = ends_with("_log2fc"),
    names_to = "group",
    values_to = "log2fc"
  ) %>%
  mutate(
    group = case_when(
      group == "anaerobes_log2fc" ~ "Obligate\nAnaerobes",
      group == "enterobact_log2fc" ~ "Enterobacteriaceae",
      group == "enterococcus_log2fc" ~ "Enterococcus"
    ),
    group = factor(group, levels = c("Obligate\nAnaerobes", "Enterobacteriaceae", "Enterococcus"))
  )

p_persistence <- ggplot(persistence_plot_data, aes(x = group, y = log2fc, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  scale_fill_manual(
    values = c(
      "Obligate\nAnaerobes" = "#2166ac",
      "Enterobacteriaceae" = "#b2182b",
      "Enterococcus" = "#d6604d"
    )
  ) +
  labs(
    title = "Persistence During Antibiotic Exposure",
    subtitle = paste("n =", nrow(paired_abund), "paired samples with antibiotics between"),
    x = NULL,
    y = "Log2 Fold Change",
    caption = "Diamond = mean. Dashed line = no change."
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(figures_dir, "persistence_boxplot.pdf"), p_persistence,
       width = 7, height = 6)
cat("Saved: figures/persistence_boxplot.pdf\n")

# 7.3 Genus-level persistence heatmap
genus_plot_data <- genus_persistence %>%
  filter(genus %in% key_genera) %>%
  mutate(
    genus = factor(genus, levels = rev(key_genera)),
    sig_label = case_when(
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ ".",
      TRUE ~ ""
    )
  )

p_genus_persist <- ggplot(genus_plot_data, aes(x = 1, y = genus, fill = mean_log2fc)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_label), size = 4) +
  geom_text(aes(x = 2, label = round(mean_retention * 100, 0)), hjust = 0, size = 3) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0,
    limits = c(-5, 5),
    oob = scales::squish,
    name = "Log2FC"
  ) +
  facet_grid(functional_group ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Genus-Level Persistence During Antibiotics",
    subtitle = "* p<0.05, ** p<0.01. Numbers = % retention",
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

ggsave(file.path(figures_dir, "genus_persistence_heatmap.pdf"), p_genus_persist,
       width = 6, height = 8)
cat("Saved: figures/genus_persistence_heatmap.pdf\n")

# 7.4 Before/after scatter for each group
scatter_data <- paired_abund %>%
  select(PairNumber, anaerobes_before, anaerobes_after,
         enterobact_before, enterobact_after,
         enterococcus_before, enterococcus_after) %>%
  pivot_longer(
    cols = -PairNumber,
    names_to = c("group", "timepoint"),
    names_pattern = "(.+)_(before|after)",
    values_to = "abundance"
  ) %>%
  pivot_wider(names_from = timepoint, values_from = abundance) %>%
  mutate(
    group = case_when(
      group == "anaerobes" ~ "Obligate Anaerobes",
      group == "enterobact" ~ "Enterobacteriaceae",
      group == "enterococcus" ~ "Enterococcus"
    )
  )

p_scatter <- ggplot(scatter_data, aes(x = before, y = after)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  facet_wrap(~ group, scales = "free") +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Abundance Before vs After Antibiotic Exposure",
    subtitle = "Dashed line = no change. Points above = increase, below = decrease.",
    x = "Abundance Before Antibiotics",
    y = "Abundance After Antibiotics"
  ) +
  theme_bw()

ggsave(file.path(figures_dir, "persistence_scatter.pdf"), p_scatter,
       width = 10, height = 4)
cat("Saved: figures/persistence_scatter.pdf\n")

# 7.5 Forest plot of genus persistence
forest_data <- genus_persistence %>%
  filter(genus %in% key_genera) %>%
  arrange(functional_group, desc(mean_log2fc)) %>%
  mutate(
    genus = factor(genus, levels = rev(genus)),
    significant = p_value < 0.05,
    direction = ifelse(mean_log2fc > 0, "Increase", "Decrease")
  )

p_forest <- ggplot(forest_data, aes(x = mean_log2fc, y = genus)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(
    aes(xmin = mean_log2fc - 1.96 * se_log2fc,
        xmax = mean_log2fc + 1.96 * se_log2fc),
    height = 0.3, color = "gray40"
  ) +
  geom_point(aes(color = functional_group, shape = significant), size = 3) +
  scale_color_manual(
    values = c(
      "Obligate Anaerobe" = "#2166ac",
      "Enterobacteriaceae" = "#b2182b",
      "Enterococcus" = "#d6604d"
    )
  ) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
  labs(
    title = "Genus-Level Persistence During Antibiotic Exposure",
    subtitle = paste("n =", nrow(paired_abund), "paired samples"),
    x = "Mean Log2 Fold Change (95% CI)",
    y = NULL,
    color = "Functional Group",
    shape = "p < 0.05"
  ) +
  theme_bw() +
  theme(legend.position = "right")

ggsave(file.path(figures_dir, "genus_persistence_forest.pdf"), p_forest,
       width = 10, height = 7)
cat("Saved: figures/genus_persistence_forest.pdf\n")

# =============================================================================
# 7. Summary and Conclusions
# =============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("SUMMARY: PERSISTENCE DURING ANTIBIOTICS\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n\n")

cat("1. PAIRED ANALYSIS (n =", nrow(paired_abund), "pairs with Abx between S1 and S2):\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")
for (i in 1:nrow(persistence_summary)) {
  grp <- persistence_summary$group[i]
  grp_name <- case_when(
    grp == "anaerobes" ~ "Obligate Anaerobes",
    grp == "enterobact" ~ "Enterobacteriaceae",
    grp == "enterococcus" ~ "Enterococcus"
  )
  cat(sprintf("   %s: %.1f%% → %.1f%% (%.0f%% retention, p=%.4f)\n",
              grp_name,
              persistence_summary$mean_S1[i] * 100,
              persistence_summary$mean_S2[i] * 100,
              persistence_summary$retention_pct[i],
              persistence_summary$p_value[i]))
}

cat("\n2. INTERPRETATION:\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("   - Anaerobes: Decrease ~2 percentage points during antibiotics\n")
cat("   - Enterobact: Roughly stable during antibiotics\n")
cat("   - Enterococcus: Increases during antibiotics (selected for)\n")

cat("\n3. GENUS-LEVEL INSIGHTS:\n")
cat("-" %>% rep(50) %>% paste(collapse = ""), "\n")
sig_genera <- genus_persistence %>%
  filter(p_value < 0.05) %>%
  arrange(mean_log2fc)
for (i in 1:min(10, nrow(sig_genera))) {
  direction <- ifelse(sig_genera$mean_log2fc[i] > 0, "INCREASES", "DECREASES")
  cat(sprintf("   %s (%s) %s during Abx (log2FC=%.2f, p=%.4f)\n",
              sig_genera$genus[i],
              sig_genera$functional_group[i],
              direction,
              sig_genera$mean_log2fc[i],
              sig_genera$p_value[i]))
}

# =============================================================================
# 9. Save Output Files
# =============================================================================

cat("\n=== Output Files ===\n\n")

# Save summary statistics
write_csv(persistence_summary, file.path(persist_dir, "group_persistence_summary.csv"))
write_csv(dominance_by_exposure, file.path(persist_dir, "dominance_by_exposure.csv"))

# Save sample-level data
sample_export <- sample_abundances %>%
  select(sample_id, MRN, PatientGroup, exposure_status,
         anaerobes, enterobact, enterococcus, opportunist_total, dominated,
         all_of(key_genera[key_genera %in% colnames(sample_abundances)]))
write_csv(sample_export, file.path(persist_dir, "sample_abundances.csv"))

outputs <- data.frame(
  file = c(
    "paired_persistence_data.csv",
    "genus_persistence.csv",
    "group_persistence_summary.csv",
    "dominance_by_exposure.csv",
    "sample_abundances.csv",
    "figures/persistence_composition.pdf",
    "figures/persistence_boxplot.pdf",
    "figures/genus_persistence_heatmap.pdf",
    "figures/persistence_scatter.pdf",
    "figures/genus_persistence_forest.pdf"
  ),
  description = c(
    "Paired sample data with group abundances before/after antibiotics",
    "Genus-level persistence statistics",
    "Summary of group-level persistence (functional groups)",
    "Dominance patterns by exposure status",
    "Sample-level abundances for all key genera",
    "Stacked bar of composition by exposure",
    "Boxplot of log2FC during antibiotics",
    "Heatmap of genus-level persistence",
    "Scatter plot of before vs after",
    "Forest plot of genus persistence"
  )
)

print(outputs)

cat("\nPersistence analysis complete!\n")
