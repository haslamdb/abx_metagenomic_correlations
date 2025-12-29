#!/usr/bin/env Rscript
# 10_enterobact_temporal_dynamics.R
# Temporal dynamics of Enterobacteriaceae expansion and recovery
#
# Key questions:
# 1. When does Enterobacteriaceae expansion occur relative to antibiotics?
# 2. Do Enterobacteriaceae recover faster than anaerobes after antibiotics stop?
# 3. Do pathogenic Enterobacteriaceae (E. coli, K. pneumoniae) behave differently?

library(tidyverse)
library(lme4)
library(lmerTest)
library(broom.mixed)

# Set paths
base_dir <- "/home/david/projects/abx_metagenomic_correlations"
data_dir <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results/new_analysis_legacy_data/enterobact_temporal")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

cat("=== Enterobacteriaceae Temporal Dynamics Analysis ===\n\n")

# Load prepared data
load(file.path(data_dir, "prepared_data.RData"))
sample_metadata <- prepared_data$sample_metadata
paired_samples <- prepared_data$paired_samples

# Load species and genus matrices
load(file.path(data_dir, "bracken_count_matrices.RData"))
species_matrix <- bracken_data$species_matrix
genus_matrix <- bracken_data$genus_matrix

# Define obligate anaerobe genera
obligate_anaerobe_genera <- c(
  "Bacteroides", "Parabacteroides", "Prevotella", "Porphyromonas",
  "Fusobacterium", "Clostridium", "Clostridioides", "Peptostreptococcus",
  "Peptococcus", "Veillonella", "Megasphaera", "Blautia", "Coprococcus",
  "Roseburia", "Faecalibacterium", "Ruminococcus", "Eubacterium",
  "Anaerostipes", "Dorea", "Subdoligranulum", "Dialister",
  "Megamonas", "Phascolarctobacterium", "Acidaminococcus", "Alistipes"
)

# =============================================================================
# SECTION 1: Define taxonomic groups at species level
# =============================================================================

cat("Section 1: Defining taxonomic groups at species level...\n")

# Calculate relative abundances at species level
species_rel <- sweep(species_matrix, 1, rowSums(species_matrix), "/") * 100

# Define pathogenic vs commensal Enterobacteriaceae
pathogenic_enterobact <- c(
  "Escherichia coli",
  "Klebsiella pneumoniae",
  "Klebsiella oxytoca",
  "Enterobacter cloacae",
  "Citrobacter freundii",
  "Serratia marcescens"
)

# Get all Enterobacteriaceae species
enterobact_genera <- c("Escherichia", "Klebsiella", "Enterobacter", "Citrobacter",
                       "Serratia", "Salmonella", "Shigella", "Proteus",
                       "Morganella", "Providencia", "Hafnia", "Raoultella")

all_enterobact_species <- colnames(species_rel)[
  sapply(colnames(species_rel), function(x) {
    any(sapply(enterobact_genera, function(g) grepl(paste0("^", g), x)))
  })
]

commensal_enterobact <- setdiff(all_enterobact_species, pathogenic_enterobact)

cat(sprintf("  Pathogenic Enterobacteriaceae species: %d\n",
            sum(pathogenic_enterobact %in% colnames(species_rel))))
cat(sprintf("  Commensal Enterobacteriaceae species: %d\n",
            sum(commensal_enterobact %in% colnames(species_rel))))

# Calculate abundances for each group
sample_taxa <- data.frame(
  sample_id = rownames(species_rel)
) %>%
  mutate(
    # Pathogenic Enterobacteriaceae
    pathogenic_enterobact = rowSums(species_rel[, intersect(pathogenic_enterobact, colnames(species_rel)), drop = FALSE]),
    # Commensal Enterobacteriaceae
    commensal_enterobact = rowSums(species_rel[, intersect(commensal_enterobact, colnames(species_rel)), drop = FALSE]),
    # Total Enterobacteriaceae
    total_enterobact = pathogenic_enterobact + commensal_enterobact,
    # Individual pathogenic species
    e_coli = if("Escherichia coli" %in% colnames(species_rel)) species_rel[, "Escherichia coli"] else 0,
    k_pneumoniae = if("Klebsiella pneumoniae" %in% colnames(species_rel)) species_rel[, "Klebsiella pneumoniae"] else 0,
    e_cloacae = if("Enterobacter cloacae" %in% colnames(species_rel)) species_rel[, "Enterobacter cloacae"] else 0
  )

# =============================================================================
# SECTION 2: Merge with metadata and calculate time-since-antibiotics
# =============================================================================

cat("\nSection 2: Calculating time since last antibiotic exposure...\n")

# Load antibiotic exposure data
abx_data <- read_csv(file.path(data_dir, "sample_metadata_with_abx.csv"), show_col_types = FALSE)

# Merge sample taxa with metadata
analysis_data <- sample_taxa %>%
  left_join(abx_data, by = "sample_id") %>%
  filter(!is.na(MRN))  # Keep only samples with metadata

# Calculate anaerobe and Enterococcus abundances from genus matrix
genus_rel <- sweep(genus_matrix, 1, rowSums(genus_matrix), "/") * 100
genus_rel[is.nan(genus_rel)] <- 0

genus_abundances <- data.frame(
  sample_id = rownames(genus_rel)
) %>%
  mutate(
    anaerobes = rowSums(genus_rel[, intersect(obligate_anaerobe_genera, colnames(genus_rel)), drop = FALSE]),
    enterococcus = if("Enterococcus" %in% colnames(genus_rel)) genus_rel[, "Enterococcus"] else 0
  )

# Join with analysis data
analysis_data <- analysis_data %>%
  left_join(genus_abundances, by = "sample_id")

# Calculate days since last antibiotic (any antibiotic)
# This requires looking at the raw drug data - for now use the 7d/14d windows
# If abx_any_7d == 1, they had abx in last 7 days
# If abx_any_7d == 0 and abx_any_14d == 1, they had abx 8-14 days ago
# If both == 0, more than 14 days since abx (or never)

analysis_data <- analysis_data %>%
  mutate(
    abx_recency = case_when(
      abx_any_7d == 1 ~ "0-7 days",
      abx_any_14d == 1 & abx_any_7d == 0 ~ "8-14 days",
      TRUE ~ ">14 days or never"
    ),
    abx_recency = factor(abx_recency, levels = c("0-7 days", "8-14 days", ">14 days or never"))
  )

cat(sprintf("  Samples with abx 0-7 days ago: %d\n", sum(analysis_data$abx_recency == "0-7 days", na.rm = TRUE)))
cat(sprintf("  Samples with abx 8-14 days ago: %d\n", sum(analysis_data$abx_recency == "8-14 days", na.rm = TRUE)))
cat(sprintf("  Samples with abx >14 days ago or never: %d\n", sum(analysis_data$abx_recency == ">14 days or never", na.rm = TRUE)))

# =============================================================================
# SECTION 3: Cross-sectional analysis by time-since-antibiotics
# =============================================================================

cat("\nSection 3: Cross-sectional analysis by time since antibiotics...\n")

# Summary statistics by recency group
cross_sectional_summary <- analysis_data %>%
  filter(!is.na(abx_recency)) %>%
  group_by(abx_recency) %>%
  summarise(
    n = n(),
    # Enterobacteriaceae
    mean_total_enterobact = mean(total_enterobact, na.rm = TRUE),
    se_total_enterobact = sd(total_enterobact, na.rm = TRUE) / sqrt(n()),
    mean_pathogenic = mean(pathogenic_enterobact, na.rm = TRUE),
    se_pathogenic = sd(pathogenic_enterobact, na.rm = TRUE) / sqrt(n()),
    mean_commensal = mean(commensal_enterobact, na.rm = TRUE),
    se_commensal = sd(commensal_enterobact, na.rm = TRUE) / sqrt(n()),
    # Specific pathogens
    mean_e_coli = mean(e_coli, na.rm = TRUE),
    mean_k_pneumoniae = mean(k_pneumoniae, na.rm = TRUE),
    # Anaerobes and Enterococcus
    mean_anaerobes = mean(anaerobes, na.rm = TRUE),
    se_anaerobes = sd(anaerobes, na.rm = TRUE) / sqrt(n()),
    mean_enterococcus = mean(enterococcus, na.rm = TRUE),
    se_enterococcus = sd(enterococcus, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

print(cross_sectional_summary)

write_csv(cross_sectional_summary, file.path(results_dir, "cross_sectional_by_abx_recency.csv"))

# Statistical tests (Kruskal-Wallis for non-parametric)
cat("\nKruskal-Wallis tests for abundance by abx recency:\n")

taxa_to_test <- c("total_enterobact", "pathogenic_enterobact", "commensal_enterobact",
                  "e_coli", "k_pneumoniae", "anaerobes", "enterococcus")

kw_results <- map_df(taxa_to_test, function(taxon) {
  test_data <- analysis_data %>% filter(!is.na(abx_recency), !is.na(!!sym(taxon)))
  if(nrow(test_data) < 10) return(NULL)

  kw <- kruskal.test(as.formula(paste(taxon, "~ abx_recency")), data = test_data)

  tibble(
    taxon = taxon,
    statistic = kw$statistic,
    p_value = kw$p.value,
    df = kw$parameter
  )
})

print(kw_results)
write_csv(kw_results, file.path(results_dir, "kruskal_wallis_by_recency.csv"))

# =============================================================================
# SECTION 4: Paired sample analysis - Recovery rates
# =============================================================================

cat("\nSection 4: Paired sample recovery rate analysis...\n")

# Get paired sample data from prepared_data
# Add S1 antibiotic exposure from abx_data
paired_data <- paired_samples %>%
  filter(high_confidence) %>%
  # Add S1 antibiotic exposure
  left_join(
    abx_data %>% select(sample_id, abx_any_7d) %>%
      rename(abx_before_s1 = abx_any_7d),
    by = c("sample1" = "sample_id")
  ) %>%
  # Add species-level taxa to paired samples
  left_join(
    sample_taxa %>% select(sample_id, pathogenic_enterobact, commensal_enterobact,
                           total_enterobact, e_coli, k_pneumoniae, e_cloacae) %>%
      rename_with(~paste0(.x, "_s1"), -sample_id),
    by = c("sample1" = "sample_id")
  ) %>%
  left_join(
    sample_taxa %>% select(sample_id, pathogenic_enterobact, commensal_enterobact,
                           total_enterobact, e_coli, k_pneumoniae, e_cloacae) %>%
      rename_with(~paste0(.x, "_s2"), -sample_id),
    by = c("sample2" = "sample_id")
  ) %>%
  # Add anaerobes and enterococcus
  left_join(
    genus_abundances %>% rename(anaerobes_s1 = anaerobes, enterococcus_s1 = enterococcus),
    by = c("sample1" = "sample_id")
  ) %>%
  left_join(
    genus_abundances %>% rename(anaerobes_s2 = anaerobes, enterococcus_s2 = enterococcus),
    by = c("sample2" = "sample_id")
  )

# Identify recovery pairs (abx before S1, no meaningful abx between)
# Using same definition as recovery_analysis.R
meaningful_abx <- c("Pip_Tazo", "Meropenem", "Cefepime", "Ceftriaxone",
                    "Ciprofloxacin", "Metronidazole", "Clindamycin",
                    "Vancomycin_IV", "Vancomycin_PO")

recovery_pairs <- paired_data %>%
  filter(abx_before_s1 == 1) %>%
  # Check for no meaningful antibiotics between samples
  mutate(
    meaningful_abx_between = rowSums(across(any_of(paste0(meaningful_abx, "_between")),
                                            ~.x > 0), na.rm = TRUE) > 0
  ) %>%
  filter(!meaningful_abx_between)

# Calculate deltas and rates for recovery pairs
recovery_pairs <- recovery_pairs %>%
  mutate(
    # Deltas (S2 - S1)
    delta_pathogenic = pathogenic_enterobact_s2 - pathogenic_enterobact_s1,
    delta_commensal = commensal_enterobact_s2 - commensal_enterobact_s1,
    delta_total_enterobact = total_enterobact_s2 - total_enterobact_s1,
    delta_e_coli = e_coli_s2 - e_coli_s1,
    delta_k_pneumoniae = k_pneumoniae_s2 - k_pneumoniae_s1,
    delta_anaerobes = anaerobes_s2 - anaerobes_s1,
    delta_enterococcus = enterococcus_s2 - enterococcus_s1,
    # Rate of change (delta per day)
    rate_pathogenic = delta_pathogenic / interval_days,
    rate_commensal = delta_commensal / interval_days,
    rate_total_enterobact = delta_total_enterobact / interval_days,
    rate_e_coli = delta_e_coli / interval_days,
    rate_k_pneumoniae = delta_k_pneumoniae / interval_days,
    rate_anaerobes = delta_anaerobes / interval_days,
    rate_enterococcus = delta_enterococcus / interval_days
  )

cat(sprintf("  Recovery pairs identified: %d\n", nrow(recovery_pairs)))

# Summary of recovery rates
recovery_rate_summary <- recovery_pairs %>%
  summarise(
    n = n(),
    # Mean delta (total change)
    mean_delta_pathogenic = mean(delta_pathogenic, na.rm = TRUE),
    mean_delta_commensal = mean(delta_commensal, na.rm = TRUE),
    mean_delta_anaerobes = mean(delta_anaerobes, na.rm = TRUE),
    mean_delta_enterococcus = mean(delta_enterococcus, na.rm = TRUE),
    # Mean rate (change per day)
    mean_rate_pathogenic = mean(rate_pathogenic, na.rm = TRUE),
    mean_rate_commensal = mean(rate_commensal, na.rm = TRUE),
    mean_rate_anaerobes = mean(rate_anaerobes, na.rm = TRUE),
    mean_rate_enterococcus = mean(rate_enterococcus, na.rm = TRUE),
    # P-values (paired t-test on delta)
    p_pathogenic = t.test(delta_pathogenic)$p.value,
    p_commensal = t.test(delta_commensal)$p.value,
    p_anaerobes = t.test(delta_anaerobes)$p.value,
    p_enterococcus = t.test(delta_enterococcus)$p.value
  )

cat("\nRecovery Summary (change from S1 on-abx to S2 off-abx):\n")
print(recovery_rate_summary)

# Compare recovery rates directly
rate_comparison <- recovery_pairs %>%
  select(PairNumber, MRN, interval_days, starts_with("rate_")) %>%
  pivot_longer(cols = starts_with("rate_"),
               names_to = "taxon",
               values_to = "rate",
               names_prefix = "rate_") %>%
  group_by(taxon) %>%
  summarise(
    n = sum(!is.na(rate)),
    mean_rate = mean(rate, na.rm = TRUE),
    se_rate = sd(rate, na.rm = TRUE) / sqrt(n()),
    median_rate = median(rate, na.rm = TRUE),
    p_vs_zero = if(sum(!is.na(rate)) > 2) t.test(rate)$p.value else NA,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_rate))

cat("\nRecovery RATE comparison (% change per day):\n")
print(rate_comparison)

write_csv(rate_comparison, file.path(results_dir, "recovery_rate_comparison.csv"))

# =============================================================================
# SECTION 5: Persistence analysis - During antibiotics
# =============================================================================

cat("\nSection 5: Persistence analysis (during antibiotic exposure)...\n")

# Identify persistence pairs (antibiotics BETWEEN S1 and S2)
persistence_pairs <- paired_data %>%
  filter(abx_between_any == 1) %>%
  mutate(
    # Deltas (S2 - S1)
    delta_pathogenic = pathogenic_enterobact_s2 - pathogenic_enterobact_s1,
    delta_commensal = commensal_enterobact_s2 - commensal_enterobact_s1,
    delta_total_enterobact = total_enterobact_s2 - total_enterobact_s1,
    delta_e_coli = e_coli_s2 - e_coli_s1,
    delta_k_pneumoniae = k_pneumoniae_s2 - k_pneumoniae_s1,
    delta_anaerobes = anaerobes_s2 - anaerobes_s1,
    delta_enterococcus = enterococcus_s2 - enterococcus_s1,
    # Rate of change
    rate_pathogenic = delta_pathogenic / interval_days,
    rate_commensal = delta_commensal / interval_days,
    rate_total_enterobact = delta_total_enterobact / interval_days,
    rate_anaerobes = delta_anaerobes / interval_days,
    rate_enterococcus = delta_enterococcus / interval_days
  )

cat(sprintf("  Persistence pairs (abx between): %d\n", nrow(persistence_pairs)))

# Persistence rate comparison
persistence_rate_comparison <- persistence_pairs %>%
  select(PairNumber, MRN, interval_days, starts_with("rate_")) %>%
  pivot_longer(cols = starts_with("rate_"),
               names_to = "taxon",
               values_to = "rate",
               names_prefix = "rate_") %>%
  group_by(taxon) %>%
  summarise(
    n = sum(!is.na(rate)),
    mean_rate = mean(rate, na.rm = TRUE),
    se_rate = sd(rate, na.rm = TRUE) / sqrt(n()),
    median_rate = median(rate, na.rm = TRUE),
    p_vs_zero = if(sum(!is.na(rate)) > 2) t.test(rate)$p.value else NA,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_rate))

cat("\nPersistence RATE comparison (% change per day during abx):\n")
print(persistence_rate_comparison)

write_csv(persistence_rate_comparison, file.path(results_dir, "persistence_rate_comparison.csv"))

# =============================================================================
# SECTION 6: Compare pathogenic vs commensal Enterobacteriaceae dynamics
# =============================================================================

cat("\nSection 6: Pathogenic vs commensal Enterobacteriaceae comparison...\n")

# Recovery: pathogenic vs commensal
pathogen_vs_commensal_recovery <- recovery_pairs %>%
  summarise(
    n = n(),
    # Pathogenic
    pathogenic_s1 = mean(pathogenic_enterobact_s1, na.rm = TRUE),
    pathogenic_s2 = mean(pathogenic_enterobact_s2, na.rm = TRUE),
    pathogenic_change = mean(delta_pathogenic, na.rm = TRUE),
    pathogenic_rate = mean(rate_pathogenic, na.rm = TRUE),
    pathogenic_p = t.test(delta_pathogenic)$p.value,
    # Commensal
    commensal_s1 = mean(commensal_enterobact_s1, na.rm = TRUE),
    commensal_s2 = mean(commensal_enterobact_s2, na.rm = TRUE),
    commensal_change = mean(delta_commensal, na.rm = TRUE),
    commensal_rate = mean(rate_commensal, na.rm = TRUE),
    commensal_p = t.test(delta_commensal)$p.value
  )

cat("\nRecovery: Pathogenic vs Commensal Enterobacteriaceae:\n")
print(pathogen_vs_commensal_recovery)

# Persistence: pathogenic vs commensal
pathogen_vs_commensal_persistence <- persistence_pairs %>%
  summarise(
    n = n(),
    # Pathogenic
    pathogenic_s1 = mean(pathogenic_enterobact_s1, na.rm = TRUE),
    pathogenic_s2 = mean(pathogenic_enterobact_s2, na.rm = TRUE),
    pathogenic_change = mean(delta_pathogenic, na.rm = TRUE),
    pathogenic_rate = mean(rate_pathogenic, na.rm = TRUE),
    pathogenic_p = t.test(delta_pathogenic)$p.value,
    # Commensal
    commensal_s1 = mean(commensal_enterobact_s1, na.rm = TRUE),
    commensal_s2 = mean(commensal_enterobact_s2, na.rm = TRUE),
    commensal_change = mean(delta_commensal, na.rm = TRUE),
    commensal_rate = mean(rate_commensal, na.rm = TRUE),
    commensal_p = t.test(delta_commensal)$p.value
  )

cat("\nPersistence: Pathogenic vs Commensal Enterobacteriaceae:\n")
print(pathogen_vs_commensal_persistence)

# Combine into summary table
dynamics_comparison <- bind_rows(
  pathogen_vs_commensal_recovery %>% mutate(phase = "Recovery (after abx)"),
  pathogen_vs_commensal_persistence %>% mutate(phase = "Persistence (during abx)")
)

write_csv(dynamics_comparison, file.path(results_dir, "pathogenic_vs_commensal_dynamics.csv"))

# =============================================================================
# SECTION 7: Stratify by specific antibiotic type
# =============================================================================

cat("\nSection 7: Recovery dynamics by prior antibiotic type...\n")

# For recovery pairs, stratify by which antibiotic was before S1
# This requires knowing which specific antibiotics were given

# Check which antibiotics were given before S1
abx_before_cols <- grep("_7d$", names(abx_data), value = TRUE)
abx_before_cols <- abx_before_cols[!grepl("any_abx|anti_|broad_", abx_before_cols)]

# Add pre-S1 antibiotic exposure to recovery pairs
recovery_by_prior_abx <- recovery_pairs %>%
  left_join(
    abx_data %>%
      select(sample_id, all_of(abx_before_cols)) %>%
      rename_with(~paste0("s1_", .x), -sample_id),
    by = c("sample1" = "sample_id")
  )

# Analyze by prior gram-negative vs non-gram-negative antibiotic
gram_neg_abx <- c("Cefepime_7d", "Meropenem_7d", "Ciprofloxacin_7d", "Pip_Tazo_7d",
                  "Ceftriaxone_7d", "TMP_SMX_7d")
non_gram_neg_abx <- c("Vancomycin_7d", "Metronidazole_7d", "Clindamycin_7d")

recovery_by_prior_abx <- recovery_by_prior_abx %>%
  mutate(
    had_gram_neg_abx = rowSums(across(paste0("s1_", gram_neg_abx), ~.x > 0), na.rm = TRUE) > 0,
    had_non_gram_neg_abx = rowSums(across(paste0("s1_", non_gram_neg_abx), ~.x > 0), na.rm = TRUE) > 0,
    prior_abx_type = case_when(
      had_gram_neg_abx & !had_non_gram_neg_abx ~ "Gram-neg active only",
      !had_gram_neg_abx & had_non_gram_neg_abx ~ "Non-gram-neg only",
      had_gram_neg_abx & had_non_gram_neg_abx ~ "Both types",
      TRUE ~ "Unknown"
    )
  )

# Summary by prior antibiotic type
recovery_by_abx_type <- recovery_by_prior_abx %>%
  filter(prior_abx_type != "Unknown") %>%
  group_by(prior_abx_type) %>%
  summarise(
    n = n(),
    # Pathogenic Enterobact
    pathogenic_s1 = mean(pathogenic_enterobact_s1, na.rm = TRUE),
    pathogenic_s2 = mean(pathogenic_enterobact_s2, na.rm = TRUE),
    pathogenic_change = mean(delta_pathogenic, na.rm = TRUE),
    pathogenic_rate = mean(rate_pathogenic, na.rm = TRUE),
    # Anaerobes
    anaerobes_s1 = mean(anaerobes_s1, na.rm = TRUE),
    anaerobes_s2 = mean(anaerobes_s2, na.rm = TRUE),
    anaerobes_change = mean(delta_anaerobes, na.rm = TRUE),
    anaerobes_rate = mean(rate_anaerobes, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nRecovery by prior antibiotic type:\n")
print(recovery_by_abx_type)

write_csv(recovery_by_abx_type, file.path(results_dir, "recovery_by_prior_abx_type.csv"))

# =============================================================================
# SECTION 7b: Recovery by INDIVIDUAL antibiotic (Enterobact vs Anaerobe rates)
# =============================================================================

cat("\nSection 7b: Recovery rates by individual prior antibiotic...\n")

# Individual antibiotics to analyze
individual_abx <- c("Pip_Tazo", "Meropenem", "Cefepime", "Metronidazole",
                    "Ciprofloxacin", "Ceftriaxone", "Vancomycin", "Clindamycin")

# For each antibiotic, calculate recovery rates for Enterobact vs Anaerobes
recovery_by_individual_abx <- map_df(individual_abx, function(abx) {
  col_name <- paste0("s1_", abx, "_7d")

  # Check if column exists

if (!(col_name %in% names(recovery_by_prior_abx))) {
    return(NULL)
  }

  # Filter to pairs where this antibiotic was given before S1
  abx_pairs <- recovery_by_prior_abx %>%
    filter(!!sym(col_name) > 0)

  if (nrow(abx_pairs) < 3) {
    return(tibble(
      antibiotic = abx,
      n_pairs = nrow(abx_pairs),
      enterobact_rate = NA_real_,
      anaerobe_rate = NA_real_,
      enterococcus_rate = NA_real_,
      ratio_enterobact_anaerobe = NA_real_,
      enterobact_p = NA_real_,
      anaerobe_p = NA_real_
    ))
  }

  # Calculate rates
  enterobact_rate <- mean(abx_pairs$rate_total_enterobact, na.rm = TRUE)
  anaerobe_rate <- mean(abx_pairs$rate_anaerobes, na.rm = TRUE)
  enterococcus_rate <- mean(abx_pairs$rate_enterococcus, na.rm = TRUE)

  # P-values (is rate significantly different from 0?)
  enterobact_p <- tryCatch(
    t.test(abx_pairs$rate_total_enterobact)$p.value,
    error = function(e) NA_real_
  )
  anaerobe_p <- tryCatch(
    t.test(abx_pairs$rate_anaerobes)$p.value,
    error = function(e) NA_real_
  )

  # Ratio (handle division by zero/near-zero)
  ratio <- ifelse(abs(anaerobe_rate) > 0.001,
                  enterobact_rate / anaerobe_rate,
                  ifelse(enterobact_rate > 0, Inf, NA_real_))

  tibble(
    antibiotic = abx,
    n_pairs = nrow(abx_pairs),
    enterobact_rate = enterobact_rate,
    anaerobe_rate = anaerobe_rate,
    enterococcus_rate = enterococcus_rate,
    ratio_enterobact_anaerobe = ratio,
    enterobact_p = enterobact_p,
    anaerobe_p = anaerobe_p
  )
})

# Also add pathogenic vs commensal breakdown by antibiotic
recovery_pathogenic_by_abx <- map_df(individual_abx, function(abx) {
  col_name <- paste0("s1_", abx, "_7d")

  if (!(col_name %in% names(recovery_by_prior_abx))) {
    return(NULL)
  }

  abx_pairs <- recovery_by_prior_abx %>%
    filter(!!sym(col_name) > 0)

  if (nrow(abx_pairs) < 3) {
    return(tibble(
      antibiotic = abx,
      n_pairs = nrow(abx_pairs),
      pathogenic_rate = NA_real_,
      commensal_rate = NA_real_,
      ratio_path_commensal = NA_real_
    ))
  }

  pathogenic_rate <- mean(abx_pairs$rate_pathogenic, na.rm = TRUE)
  commensal_rate <- mean(abx_pairs$rate_commensal, na.rm = TRUE)

  ratio <- ifelse(abs(commensal_rate) > 0.001,
                  pathogenic_rate / commensal_rate,
                  ifelse(pathogenic_rate > 0, Inf, NA_real_))

  tibble(
    antibiotic = abx,
    n_pairs = nrow(abx_pairs),
    pathogenic_rate = pathogenic_rate,
    commensal_rate = commensal_rate,
    ratio_path_commensal = ratio
  )
})

# Merge the two analyses
recovery_by_individual_abx <- recovery_by_individual_abx %>%
  left_join(
    recovery_pathogenic_by_abx %>% select(antibiotic, pathogenic_rate, commensal_rate, ratio_path_commensal),
    by = "antibiotic"
  ) %>%
  arrange(desc(n_pairs))

cat("\nRecovery rates by individual prior antibiotic (paired samples):\n")
print(recovery_by_individual_abx %>%
        select(antibiotic, n_pairs, enterobact_rate, anaerobe_rate, ratio_enterobact_anaerobe,
               pathogenic_rate, commensal_rate) %>%
        mutate(across(where(is.numeric) & !matches("n_pairs"), ~round(.x, 3))))

write_csv(recovery_by_individual_abx, file.path(results_dir, "recovery_by_individual_antibiotic.csv"))

# Summary table for display
cat("\n\nKey finding - Enterobact vs Anaerobe recovery by antibiotic:\n")
cat("(Positive rate = increasing after abx stops, negative = decreasing)\n\n")

recovery_by_individual_abx %>%
  filter(n_pairs >= 5) %>%
  select(antibiotic, n_pairs, enterobact_rate, anaerobe_rate, ratio_enterobact_anaerobe) %>%
  mutate(
    enterobact_rate = sprintf("%+.2f%%/day", enterobact_rate),
    anaerobe_rate = sprintf("%+.2f%%/day", anaerobe_rate),
    ratio = ifelse(is.infinite(ratio_enterobact_anaerobe), "Inf",
                   ifelse(is.na(ratio_enterobact_anaerobe), "NA",
                          sprintf("%.1fx", ratio_enterobact_anaerobe)))
  ) %>%
  select(-ratio_enterobact_anaerobe) %>%
  rename(`Enterobact Rate` = enterobact_rate,
         `Anaerobe Rate` = anaerobe_rate,
         `Ratio (E/A)` = ratio) %>%
  print()

# =============================================================================
# SECTION 8: Visualizations
# =============================================================================

cat("\nSection 8: Creating visualizations...\n")

# Plot 1: Cross-sectional abundance by time since abx
plot_data <- analysis_data %>%
  filter(!is.na(abx_recency)) %>%
  select(sample_id, abx_recency, pathogenic_enterobact, commensal_enterobact,
         anaerobes, enterococcus) %>%
  pivot_longer(cols = c(pathogenic_enterobact, commensal_enterobact, anaerobes, enterococcus),
               names_to = "taxon", values_to = "abundance") %>%
  mutate(
    taxon = case_when(
      taxon == "pathogenic_enterobact" ~ "Pathogenic Enterobact",
      taxon == "commensal_enterobact" ~ "Commensal Enterobact",
      taxon == "anaerobes" ~ "Obligate Anaerobes",
      taxon == "enterococcus" ~ "Enterococcus"
    ),
    taxon = factor(taxon, levels = c("Pathogenic Enterobact", "Commensal Enterobact",
                                      "Obligate Anaerobes", "Enterococcus"))
  )

p1 <- ggplot(plot_data, aes(x = abx_recency, y = abundance, fill = taxon)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~taxon, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Pathogenic Enterobact" = "#E41A1C",
                               "Commensal Enterobact" = "#377EB8",
                               "Obligate Anaerobes" = "#4DAF4A",
                               "Enterococcus" = "#984EA3")) +
  labs(title = "Microbial Abundance by Time Since Antibiotic Exposure",
       subtitle = "Cross-sectional analysis of all samples",
       x = "Time since last antibiotic",
       y = "Relative abundance (%)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(results_dir, "figures", "cross_sectional_by_abx_recency.pdf"),
       p1, width = 12, height = 4)

# Plot 2: Recovery rates comparison
rate_plot_data <- rate_comparison %>%
  mutate(
    taxon_label = case_when(
      taxon == "pathogenic" ~ "Pathogenic\nEnterobact",
      taxon == "commensal" ~ "Commensal\nEnterobact",
      taxon == "total_enterobact" ~ "Total\nEnterobact",
      taxon == "anaerobes" ~ "Obligate\nAnaerobes",
      taxon == "enterococcus" ~ "Enterococcus",
      taxon == "e_coli" ~ "E. coli",
      taxon == "k_pneumoniae" ~ "K. pneumoniae"
    ),
    significant = p_vs_zero < 0.05
  )

p2 <- ggplot(rate_plot_data, aes(x = reorder(taxon_label, mean_rate), y = mean_rate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_col(aes(fill = mean_rate > 0), width = 0.7) +
  geom_errorbar(aes(ymin = mean_rate - se_rate, ymax = mean_rate + se_rate),
                width = 0.2) +
  geom_text(aes(label = ifelse(significant, "*", ""),
                y = mean_rate + sign(mean_rate) * se_rate + 0.05),
            size = 6) +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                    guide = "none") +
  labs(title = "Recovery Rate After Antibiotic Cessation",
       subtitle = sprintf("n = %d paired samples (abx before S1, no abx between S1-S2)",
                         nrow(recovery_pairs)),
       x = "",
       y = "Rate of change (% per day)") +
  theme_bw() +
  coord_flip()

ggsave(file.path(results_dir, "figures", "recovery_rates.pdf"),
       p2, width = 8, height = 6)

# Plot 3: Persistence vs Recovery comparison
combined_rates <- bind_rows(
  rate_comparison %>% mutate(phase = "Recovery\n(after abx)"),
  persistence_rate_comparison %>% mutate(phase = "Persistence\n(during abx)")
) %>%
  filter(taxon %in% c("pathogenic", "commensal", "anaerobes", "enterococcus")) %>%
  mutate(
    taxon_label = case_when(
      taxon == "pathogenic" ~ "Pathogenic Enterobact",
      taxon == "commensal" ~ "Commensal Enterobact",
      taxon == "anaerobes" ~ "Obligate Anaerobes",
      taxon == "enterococcus" ~ "Enterococcus"
    )
  )

p3 <- ggplot(combined_rates, aes(x = taxon_label, y = mean_rate, fill = phase)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_rate - se_rate, ymax = mean_rate + se_rate),
                position = position_dodge(0.8), width = 0.2) +
  scale_fill_manual(values = c("Recovery\n(after abx)" = "#4DAF4A",
                               "Persistence\n(during abx)" = "#E41A1C")) +
  labs(title = "Rate of Change: During vs After Antibiotic Exposure",
       x = "",
       y = "Rate of change (% per day)",
       fill = "Phase") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(results_dir, "figures", "persistence_vs_recovery_rates.pdf"),
       p3, width = 10, height = 6)

# Plot 4: Pathogenic vs commensal scatter
scatter_data <- recovery_pairs %>%
  select(PairNumber, delta_pathogenic, delta_commensal, interval_days) %>%
  filter(!is.na(delta_pathogenic) & !is.na(delta_commensal))

p4 <- ggplot(scatter_data, aes(x = delta_commensal, y = delta_pathogenic)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_point(aes(size = interval_days), alpha = 0.6, color = "#E41A1C") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  labs(title = "Pathogenic vs Commensal Enterobacteriaceae Recovery",
       subtitle = sprintf("n = %d recovery pairs", nrow(scatter_data)),
       x = "Change in commensal Enterobacteriaceae (%)",
       y = "Change in pathogenic Enterobacteriaceae (%)",
       size = "Days between\nsamples") +
  theme_bw()

# Add correlation
cor_test <- cor.test(scatter_data$delta_commensal, scatter_data$delta_pathogenic)
p4 <- p4 + annotate("text", x = Inf, y = Inf,
                    label = sprintf("r = %.2f, p = %.3f", cor_test$estimate, cor_test$p.value),
                    hjust = 1.1, vjust = 1.5)

ggsave(file.path(results_dir, "figures", "pathogenic_vs_commensal_recovery.pdf"),
       p4, width = 8, height = 6)

# =============================================================================
# SECTION 9: Summary statistics
# =============================================================================

cat("\nSection 9: Final summary...\n")

summary_stats <- list(
  cross_sectional = cross_sectional_summary,
  recovery_rates = rate_comparison,
  persistence_rates = persistence_rate_comparison,
  pathogenic_vs_commensal = dynamics_comparison,
  recovery_by_abx_type = recovery_by_abx_type
)

saveRDS(summary_stats, file.path(results_dir, "summary_statistics.rds"))

# Create text summary
sink(file.path(results_dir, "analysis_summary.txt"))
cat("=======================================================\n")
cat("ENTEROBACTERIACEAE TEMPORAL DYNAMICS ANALYSIS SUMMARY\n")
cat("=======================================================\n\n")

cat("1. CROSS-SECTIONAL ANALYSIS (all samples - NOT paired)\n")
cat("-------------------------------------------------------\n")
cat("NOTE: This compares different samples at different time points.\n")
cat("The '>14 days or never' group may include patients who never had abx.\n\n")
cat(sprintf("Total samples analyzed: %d\n", nrow(analysis_data)))
cat("\nMean abundance by time since antibiotics:\n")
print(cross_sectional_summary %>%
        select(abx_recency, n, mean_pathogenic, mean_commensal, mean_anaerobes, mean_enterococcus))

cat("\n\n2. RECOVERY ANALYSIS - PAIRED SAMPLES (same patient, S1->S2)\n")
cat("------------------------------------------------------------\n")
cat("Abx before S1, no meaningful abx between S1 and S2.\n")
cat("This is the KEY analysis - tracks actual within-patient recovery.\n\n")
cat(sprintf("Recovery pairs: %d\n", nrow(recovery_pairs)))
cat("\nRecovery rates (% change per day):\n")
print(rate_comparison)

cat("\n\n3. PERSISTENCE ANALYSIS (paired samples, abx between S1-S2)\n")
cat("------------------------------------------------------------\n")
cat(sprintf("Persistence pairs: %d\n", nrow(persistence_pairs)))
cat("\nPersistence rates (% change per day during abx):\n")
print(persistence_rate_comparison)

cat("\n\n4. PATHOGENIC vs COMMENSAL ENTEROBACTERIACEAE\n")
cat("----------------------------------------------\n")
cat("\nRecovery phase:\n")
cat(sprintf("  Pathogenic: S1=%.1f%% -> S2=%.1f%% (change=%.2f%%, rate=%.3f%%/day)\n",
            pathogen_vs_commensal_recovery$pathogenic_s1,
            pathogen_vs_commensal_recovery$pathogenic_s2,
            pathogen_vs_commensal_recovery$pathogenic_change,
            pathogen_vs_commensal_recovery$pathogenic_rate))
cat(sprintf("  Commensal:  S1=%.1f%% -> S2=%.1f%% (change=%.2f%%, rate=%.3f%%/day)\n",
            pathogen_vs_commensal_recovery$commensal_s1,
            pathogen_vs_commensal_recovery$commensal_s2,
            pathogen_vs_commensal_recovery$commensal_change,
            pathogen_vs_commensal_recovery$commensal_rate))

cat("\nPersistence phase:\n")
cat(sprintf("  Pathogenic: S1=%.1f%% -> S2=%.1f%% (change=%.2f%%, rate=%.3f%%/day)\n",
            pathogen_vs_commensal_persistence$pathogenic_s1,
            pathogen_vs_commensal_persistence$pathogenic_s2,
            pathogen_vs_commensal_persistence$pathogenic_change,
            pathogen_vs_commensal_persistence$pathogenic_rate))
cat(sprintf("  Commensal:  S1=%.1f%% -> S2=%.1f%% (change=%.2f%%, rate=%.3f%%/day)\n",
            pathogen_vs_commensal_persistence$commensal_s1,
            pathogen_vs_commensal_persistence$commensal_s2,
            pathogen_vs_commensal_persistence$commensal_change,
            pathogen_vs_commensal_persistence$commensal_rate))

cat("\n\n5. RECOVERY BY INDIVIDUAL ANTIBIOTIC (paired samples)\n")
cat("-----------------------------------------------------\n")
if (exists("recovery_by_individual_abx") && nrow(recovery_by_individual_abx) > 0) {
  cat("\nEnterobact vs Anaerobe recovery rates by prior antibiotic:\n")
  print(recovery_by_individual_abx %>%
          filter(n_pairs >= 3) %>%
          select(antibiotic, n_pairs, enterobact_rate, anaerobe_rate, ratio_enterobact_anaerobe) %>%
          mutate(across(where(is.numeric) & !matches("n_pairs|ratio"), ~round(.x, 3))))
}

cat("\n\n6. KEY FINDINGS\n")
cat("---------------\n")

# Determine key findings
recovery_enterobact_rate <- rate_comparison %>% filter(taxon == "total_enterobact") %>% pull(mean_rate)
recovery_anaerobe_rate <- rate_comparison %>% filter(taxon == "anaerobes") %>% pull(mean_rate)
recovery_pathogenic_rate <- rate_comparison %>% filter(taxon == "pathogenic") %>% pull(mean_rate)
recovery_commensal_rate <- rate_comparison %>% filter(taxon == "commensal") %>% pull(mean_rate)

cat(sprintf("- Total Enterobacteriaceae recovery rate: %.3f%%/day\n", recovery_enterobact_rate))
cat(sprintf("- Obligate anaerobe recovery rate: %.3f%%/day\n", recovery_anaerobe_rate))
cat(sprintf("- Ratio (Enterobact/Anaerobe): %.2fx\n",
            ifelse(recovery_anaerobe_rate != 0, recovery_enterobact_rate/recovery_anaerobe_rate, NA)))
cat(sprintf("\n- Pathogenic Enterobact recovery rate: %.3f%%/day\n", recovery_pathogenic_rate))
cat(sprintf("- Commensal Enterobact recovery rate: %.3f%%/day\n", recovery_commensal_rate))
cat(sprintf("- Ratio (Pathogenic/Commensal): %.2fx\n",
            ifelse(recovery_commensal_rate != 0, recovery_pathogenic_rate/recovery_commensal_rate, NA)))

sink()

cat("\nAnalysis complete! Results saved to:\n")
cat(sprintf("  %s\n", results_dir))

cat("\nKey output files:\n")
cat("  - cross_sectional_by_abx_recency.csv\n")
cat("  - recovery_rate_comparison.csv\n")
cat("  - persistence_rate_comparison.csv\n")
cat("  - pathogenic_vs_commensal_dynamics.csv\n")
cat("  - analysis_summary.txt\n")
cat("  - figures/*.pdf\n")
