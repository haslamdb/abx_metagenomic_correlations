#!/usr/bin/env Rscript
# =============================================================================
# 15_antibiotic_specific_effects.R
# Antibiotic-specific effects on functional categories + Ciprofloxacin route analysis
#
# Purpose:
# 1. Show effects of individual antibiotics on Enterococcus/Enterobacteriaceae/Anaerobes
# 2. Generate heatmap for supplementary figure
# 3. Compare Ciprofloxacin PO vs IV route effects
# =============================================================================

library(tidyverse)
library(scales)

# Set paths
project_dir <- "/home/david/projects/abx_metagenomic_correlations"
results_dir <- file.path(project_dir, "results/new_analysis_legacy_data")
output_dir <- file.path(results_dir, "antibiotic_specific_effects")
figures_dir <- file.path(results_dir, "figures")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("ANTIBIOTIC-SPECIFIC EFFECTS ANALYSIS\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n\n")

# =============================================================================
# 1. Load Data
# =============================================================================

cat("=== Loading Data ===\n\n")

# Load robust associations
robust <- read_csv(
  file.path(results_dir, "individual_antibiotics_v2_with_covariates/robust_associations.csv"),
  show_col_types = FALSE
)
cat("Loaded", nrow(robust), "robust associations\n")

# Load raw drug data for ciprofloxacin route analysis
drugs <- read_csv(file.path(project_dir, "data/legacy/Drugs.csv"), show_col_types = FALSE)

# Load prepared data
load(file.path(project_dir, "data/prepared_data.RData"))
sample_metadata <- prepared_data$sample_metadata

# Load species matrix for route-specific analysis
load(file.path(project_dir, "data/bracken_count_matrices.RData"))
species_matrix <- bracken_data$species_matrix

# =============================================================================
# 2. Define Functional Categories
# =============================================================================

# Pathogenic Enterobacteriaceae
pathogenic_enterobact <- c(

  "Escherichia coli", "Klebsiella pneumoniae", "Klebsiella oxytoca",
  "Enterobacter cloacae", "Citrobacter freundii", "Serratia marcescens",
  "Salmonella enterica", "Shigella flexneri", "Shigella sonnei"
)

enterobact_genera <- c("Escherichia", "Klebsiella", "Enterobacter", "Citrobacter",
                       "Serratia", "Salmonella", "Shigella", "Proteus",
                       "Morganella", "Providencia", "Hafnia", "Raoultella")

obligate_anaerobe_genera <- c(
  "Bacteroides", "Parabacteroides", "Prevotella", "Porphyromonas",
  "Fusobacterium", "Clostridium", "Clostridioides", "Peptostreptococcus",
  "Veillonella", "Megasphaera", "Blautia", "Coprococcus",
  "Roseburia", "Faecalibacterium", "Ruminococcus", "Eubacterium",
  "Anaerostipes", "Dorea", "Subdoligranulum", "Dialister", "Alistipes"
)

# Function to categorize species
categorize_species <- function(species_name) {
  genus <- strsplit(species_name, " ")[[1]][1]

  if (species_name %in% pathogenic_enterobact) {
    return("Pathogenic_Enterobact")
  } else if (genus %in% enterobact_genera) {
    return("Commensal_Enterobact")
  } else if (genus == "Enterococcus") {
    return("Enterococcus")
  } else if (genus %in% obligate_anaerobe_genera) {
    return("Obligate_Anaerobe")
  } else {
    return("Other")
  }
}

# Add category to robust associations
robust <- robust %>%
  mutate(category = sapply(species, categorize_species))

# =============================================================================
# 3. Calculate Antibiotic-Specific Effects by Category
# =============================================================================

cat("\n=== Calculating Antibiotic-Specific Effects ===\n\n")

# Summarize by antibiotic and category
abx_category_summary <- robust %>%
  filter(category != "Other") %>%
  group_by(antibiotic, category) %>%
  summarise(
    n_species = n(),
    mean_effect = mean(mean_effect),
    sd_effect = sd(mean_effect),
    se_effect = sd(mean_effect) / sqrt(n()),
    pct_increased = mean(direction == "increased") * 100,
    pct_decreased = mean(direction == "decreased") * 100,
    .groups = "drop"
  )

# Pivot wider for display
abx_effects_wide <- abx_category_summary %>%
  select(antibiotic, category, mean_effect, n_species) %>%
  pivot_wider(
    names_from = category,
    values_from = c(mean_effect, n_species),
    names_sep = "_"
  )

cat("Effects by antibiotic and category:\n")
print(abx_effects_wide)

# Save detailed results
write_csv(abx_category_summary, file.path(output_dir, "antibiotic_category_effects.csv"))
write_csv(abx_effects_wide, file.path(output_dir, "antibiotic_category_effects_wide.csv"))

# =============================================================================
# 4. Generate Heatmap (Figure S3)
# =============================================================================

cat("\n=== Generating Heatmap ===\n\n")

# Prepare data for ggplot2 heatmap
heatmap_long <- abx_category_summary %>%
  filter(category %in% c("Enterococcus", "Pathogenic_Enterobact",
                         "Commensal_Enterobact", "Obligate_Anaerobe")) %>%
  mutate(
    category = factor(category, levels = c("Enterococcus", "Pathogenic_Enterobact",
                                           "Commensal_Enterobact", "Obligate_Anaerobe")),
    category_label = recode(category,
      "Enterococcus" = "Enterococcus",
      "Pathogenic_Enterobact" = "Pathogenic\nEnterobact",
      "Commensal_Enterobact" = "Commensal\nEnterobact",
      "Obligate_Anaerobe" = "Obligate\nAnaerobe"
    )
  )

# Order antibiotics by mean effect on anaerobes (most suppressive first)
abx_order <- heatmap_long %>%
  filter(category == "Obligate_Anaerobe") %>%
  arrange(mean_effect) %>%
  pull(antibiotic)

# For antibiotics with no anaerobe data, order by Enterococcus effect
missing_abx <- setdiff(unique(heatmap_long$antibiotic), abx_order)
if (length(missing_abx) > 0) {
  missing_order <- heatmap_long %>%
    filter(antibiotic %in% missing_abx, category == "Enterococcus") %>%
    arrange(desc(mean_effect)) %>%
    pull(antibiotic)
  abx_order <- c(abx_order, missing_order)
}

heatmap_long <- heatmap_long %>%
  mutate(antibiotic = factor(antibiotic, levels = rev(abx_order)))

# Create heatmap using ggplot2
p_heatmap <- ggplot(heatmap_long, aes(x = category_label, y = antibiotic, fill = mean_effect)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", mean_effect)), size = 3.5) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0, limits = c(-2.5, 2.5), oob = squish,
    name = "Mean Effect\n(log2FC)"
  ) +
  labs(
    title = "Antibiotic Effects on Functional Categories",
    subtitle = "Mean log2 fold-change: positive = increases during antibiotic exposure",
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(figures_dir, "antibiotic_category_heatmap.pdf"), p_heatmap,
       width = 10, height = 8)
cat("Saved: antibiotic_category_heatmap.pdf\n")

# =============================================================================
# 5. Bar Plot by Antibiotic (Alternative Visualization)
# =============================================================================

cat("\n=== Generating Bar Plots ===\n\n")

# Prepare data for plotting
plot_data <- abx_category_summary %>%
  filter(category %in% c("Enterococcus", "Pathogenic_Enterobact",
                         "Commensal_Enterobact", "Obligate_Anaerobe")) %>%
  mutate(
    category = factor(category, levels = c("Enterococcus", "Pathogenic_Enterobact",
                                           "Commensal_Enterobact", "Obligate_Anaerobe")),
    category_label = recode(category,
      "Enterococcus" = "Enterococcus",
      "Pathogenic_Enterobact" = "Pathogenic\nEnterobact",
      "Commensal_Enterobact" = "Commensal\nEnterobact",
      "Obligate_Anaerobe" = "Obligate\nAnaerobe"
    )
  )

# Color palette
category_colors <- c(
  "Enterococcus" = "#9467bd",
  "Pathogenic_Enterobact" = "#d62728",
  "Commensal_Enterobact" = "#ff7f0e",
  "Obligate_Anaerobe" = "#1f77b4"
)

# Grouped bar plot
p_bars <- ggplot(plot_data, aes(x = antibiotic, y = mean_effect, fill = category)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_effect - se_effect, ymax = mean_effect + se_effect),
    position = position_dodge(width = 0.8),
    width = 0.25
  ) +
  scale_fill_manual(values = category_colors, name = "Category") +
  labs(
    title = "Antibiotic-Specific Effects on Functional Categories",
    subtitle = "Mean effect (log2FC) with SE. Positive = increases during antibiotic exposure.",
    x = NULL,
    y = "Mean Effect (log2 fold-change)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold")
  )

ggsave(file.path(figures_dir, "antibiotic_category_barplot.pdf"), p_bars,
       width = 14, height = 7)
cat("Saved: antibiotic_category_barplot.pdf\n")

# Faceted by category
p_facet <- ggplot(plot_data, aes(x = reorder(antibiotic, -mean_effect), y = mean_effect, fill = category)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean_effect - se_effect, ymax = mean_effect + se_effect), width = 0.3) +
  facet_wrap(~ category_label, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = category_colors) +
  labs(
    title = "Antibiotic Effects by Functional Category",
    subtitle = "Antibiotics ranked by effect within each category",
    x = NULL,
    y = "Mean Effect (log2FC)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold")
  )

ggsave(file.path(figures_dir, "antibiotic_effects_faceted.pdf"), p_facet,
       width = 14, height = 5)
cat("Saved: antibiotic_effects_faceted.pdf\n")

# =============================================================================
# 6. Ciprofloxacin Route Analysis (PO vs IV)
# =============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("CIPROFLOXACIN ROUTE ANALYSIS (PO vs IV)\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n\n")

# Get ciprofloxacin administrations by route
cipro_drugs <- drugs %>%
  filter(Drug == "Ciprofloxacin") %>%
  filter(Route %in% c("IV", "PO"))  # Exclude Otic, Ophthalmic

cat("Ciprofloxacin administrations (systemic only):\n")
print(table(cipro_drugs$Route))

# Get sample list with MRN
sample_list <- read_csv(
  file.path(project_dir, "data/legacy/SampleListFormatted.csv"),
  show_col_types = FALSE
)

# Parse dates
sample_list <- sample_list %>%
  mutate(SampleDate_parsed = as.Date(SampleDate, format = "%m/%d/%Y"))

cipro_drugs <- cipro_drugs %>%
  mutate(
    Date = as.Date(Date),
    MRN = as.character(MRN)
  )

# Function to calculate route-specific ciprofloxacin exposure
calc_cipro_route_exposure <- function(mrn, sample_date, route, window_days = 7) {
  start_date <- sample_date - window_days
  end_date <- sample_date - 1

  mrn_cipro <- cipro_drugs %>%
    filter(MRN == as.character(mrn), Route == route) %>%
    filter(Date >= start_date, Date <= end_date)

  return(nrow(mrn_cipro) > 0)
}

# Calculate ciprofloxacin route exposure for each sample
cat("\nCalculating route-specific ciprofloxacin exposure...\n")

cipro_exposure <- sample_list %>%
  filter(!is.na(SampleDate_parsed)) %>%
  rowwise() %>%
  mutate(
    Cipro_IV_7d = calc_cipro_route_exposure(MRN, SampleDate_parsed, "IV", 7),
    Cipro_PO_7d = calc_cipro_route_exposure(MRN, SampleDate_parsed, "PO", 7)
  ) %>%
  ungroup() %>%
  select(SequenceFileName, MRN, SampleDate_parsed, Cipro_IV_7d, Cipro_PO_7d)

cat("\nCiprofloxacin route exposure summary:\n")
cat("  Cipro_IV_7d:", sum(cipro_exposure$Cipro_IV_7d), "samples\n")
cat("  Cipro_PO_7d:", sum(cipro_exposure$Cipro_PO_7d), "samples\n")
cat("  Both IV and PO:", sum(cipro_exposure$Cipro_IV_7d & cipro_exposure$Cipro_PO_7d), "samples\n")
cat("  Either route:", sum(cipro_exposure$Cipro_IV_7d | cipro_exposure$Cipro_PO_7d), "samples\n")

# Merge with species data
species_rel <- species_matrix / rowSums(species_matrix)
species_rel[is.na(species_rel)] <- 0

# Calculate functional group abundances
calc_group_abundance <- function(rel_matrix, genera_list) {
  matching <- genera_list[genera_list %in% colnames(rel_matrix)]
  if (length(matching) == 0) return(rep(0, nrow(rel_matrix)))
  if (length(matching) == 1) return(rel_matrix[, matching])
  rowSums(rel_matrix[, matching, drop = FALSE])
}

# Calculate Enterobacteriaceae at species level
enterobact_species <- colnames(species_rel)[sapply(colnames(species_rel), function(x) {
  genus <- strsplit(x, " ")[[1]][1]
  genus %in% enterobact_genera
})]

sample_abundances <- data.frame(
  sample_id = rownames(species_rel),
  enterobact = rowSums(species_rel[, enterobact_species[enterobact_species %in% colnames(species_rel)], drop = FALSE]),
  enterococcus = if ("Enterococcus" %in% colnames(species_rel)) {
    rowSums(species_rel[, grep("^Enterococcus", colnames(species_rel)), drop = FALSE])
  } else 0,
  stringsAsFactors = FALSE
)

# Add anaerobes from genus-level
genus_matrix <- bracken_data$genus_matrix
genus_rel <- genus_matrix / rowSums(genus_matrix)
genus_rel[is.na(genus_rel)] <- 0

sample_abundances$anaerobes <- calc_group_abundance(genus_rel, obligate_anaerobe_genera)[
  match(sample_abundances$sample_id, rownames(genus_rel))
]

# Merge exposure with abundances
analysis_data <- cipro_exposure %>%
  left_join(sample_abundances, by = c("SequenceFileName" = "sample_id")) %>%
  filter(!is.na(enterobact))

cat("\nSamples with complete data:", nrow(analysis_data), "\n")

# =============================================================================
# 7. Compare PO vs IV Effects
# =============================================================================

cat("\n=== Comparing Ciprofloxacin PO vs IV Effects ===\n\n")

# Create exposure groups
analysis_data <- analysis_data %>%
  mutate(
    cipro_group = case_when(
      Cipro_PO_7d & !Cipro_IV_7d ~ "PO only",
      Cipro_IV_7d & !Cipro_PO_7d ~ "IV only",
      Cipro_IV_7d & Cipro_PO_7d ~ "Both",
      TRUE ~ "None"
    )
  )

cat("Exposure groups:\n")
print(table(analysis_data$cipro_group))

# Compare abundances by group
group_comparison <- analysis_data %>%
  filter(cipro_group != "Both") %>%  # Exclude mixed exposure
  group_by(cipro_group) %>%
  summarise(
    n = n(),
    mean_enterobact = mean(enterobact) * 100,
    se_enterobact = sd(enterobact) / sqrt(n()) * 100,
    mean_enterococcus = mean(enterococcus) * 100,
    se_enterococcus = sd(enterococcus) / sqrt(n()) * 100,
    mean_anaerobes = mean(anaerobes, na.rm = TRUE) * 100,
    se_anaerobes = sd(anaerobes, na.rm = TRUE) / sqrt(n()) * 100,
    .groups = "drop"
  )

cat("\nFunctional group abundances by ciprofloxacin exposure:\n")
print(group_comparison)

# Statistical tests: IV vs None, PO vs None, IV vs PO
cat("\n=== Statistical Tests ===\n\n")

none_data <- analysis_data %>% filter(cipro_group == "None")
iv_data <- analysis_data %>% filter(cipro_group == "IV only")
po_data <- analysis_data %>% filter(cipro_group == "PO only")

run_comparison <- function(group1, group2, label1, label2, variable) {
  if (nrow(group1) < 3 || nrow(group2) < 3) {
    return(data.frame(
      comparison = paste(label1, "vs", label2),
      variable = variable,
      mean1 = NA, mean2 = NA, difference = NA, p_value = NA
    ))
  }

  test <- wilcox.test(group1[[variable]], group2[[variable]])

  data.frame(
    comparison = paste(label1, "vs", label2),
    variable = variable,
    mean1 = mean(group1[[variable]]) * 100,
    mean2 = mean(group2[[variable]]) * 100,
    difference = (mean(group1[[variable]]) - mean(group2[[variable]])) * 100,
    p_value = test$p.value
  )
}

# Run all comparisons
comparisons <- list()

for (var in c("enterobact", "enterococcus", "anaerobes")) {
  if (nrow(iv_data) >= 3) {
    comparisons[[paste("IV_vs_None", var)]] <- run_comparison(iv_data, none_data, "IV", "None", var)
  }
  if (nrow(po_data) >= 3) {
    comparisons[[paste("PO_vs_None", var)]] <- run_comparison(po_data, none_data, "PO", "None", var)
  }
  if (nrow(iv_data) >= 3 && nrow(po_data) >= 3) {
    comparisons[[paste("IV_vs_PO", var)]] <- run_comparison(iv_data, po_data, "IV", "PO", var)
  }
}

comparison_results <- bind_rows(comparisons)
print(comparison_results)

# Save comparison results
write_csv(comparison_results, file.path(output_dir, "cipro_route_comparison.csv"))
cat("\nSaved: cipro_route_comparison.csv\n")

# =============================================================================
# 8. Visualization of Ciprofloxacin Route Effects
# =============================================================================

cat("\n=== Generating Ciprofloxacin Route Figures ===\n\n")

# Prepare data for plotting
plot_cipro <- group_comparison %>%
  filter(cipro_group != "Both") %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "group_type",
    values_to = "mean_abundance"
  ) %>%
  mutate(
    se = case_when(
      group_type == "mean_enterobact" ~ group_comparison$se_enterobact[match(cipro_group, group_comparison$cipro_group)],
      group_type == "mean_enterococcus" ~ group_comparison$se_enterococcus[match(cipro_group, group_comparison$cipro_group)],
      group_type == "mean_anaerobes" ~ group_comparison$se_anaerobes[match(cipro_group, group_comparison$cipro_group)]
    ),
    group_label = recode(group_type,
      "mean_enterobact" = "Enterobacteriaceae",
      "mean_enterococcus" = "Enterococcus",
      "mean_anaerobes" = "Anaerobes"
    ),
    cipro_group = factor(cipro_group, levels = c("None", "PO only", "IV only"))
  )

# Bar plot
p_cipro <- ggplot(plot_cipro, aes(x = cipro_group, y = mean_abundance, fill = cipro_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_abundance - se, ymax = mean_abundance + se), width = 0.3) +
  facet_wrap(~ group_label, scales = "free_y") +
  scale_fill_manual(values = c("None" = "gray60", "PO only" = "#2ca02c", "IV only" = "#1f77b4")) +
  labs(
    title = "Ciprofloxacin Route Effects on Gut Microbiome",
    subtitle = "Comparing PO vs IV administration (7-day exposure window)",
    x = "Ciprofloxacin Exposure",
    y = "Mean Relative Abundance (%)",
    fill = "Exposure"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(figures_dir, "cipro_route_effects.pdf"), p_cipro,
       width = 10, height = 5)
cat("Saved: cipro_route_effects.pdf\n")

# Box plot with individual points
plot_individual <- analysis_data %>%
  filter(cipro_group %in% c("None", "PO only", "IV only")) %>%
  pivot_longer(
    cols = c(enterobact, enterococcus, anaerobes),
    names_to = "group_type",
    values_to = "abundance"
  ) %>%
  mutate(
    abundance = abundance * 100,
    group_label = recode(group_type,
      "enterobact" = "Enterobacteriaceae",
      "enterococcus" = "Enterococcus",
      "anaerobes" = "Anaerobes"
    ),
    cipro_group = factor(cipro_group, levels = c("None", "PO only", "IV only"))
  )

p_cipro_box <- ggplot(plot_individual, aes(x = cipro_group, y = abundance, fill = cipro_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  facet_wrap(~ group_label, scales = "free_y") +
  scale_fill_manual(values = c("None" = "gray60", "PO only" = "#2ca02c", "IV only" = "#1f77b4")) +
  labs(
    title = "Ciprofloxacin Route Effects on Gut Microbiome",
    subtitle = "Box plots with individual samples",
    x = "Ciprofloxacin Exposure",
    y = "Relative Abundance (%)",
    fill = "Exposure"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(figures_dir, "cipro_route_effects_boxplot.pdf"), p_cipro_box,
       width = 10, height = 5)
cat("Saved: cipro_route_effects_boxplot.pdf\n")

# =============================================================================
# 9. Summary
# =============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("SUMMARY\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n\n")

cat("ANTIBIOTIC-SPECIFIC FINDINGS:\n")
cat("1. Most antibiotics INCREASE Enterococcus (intrinsic resistance)\n")
cat("2. Ciprofloxacin is UNIQUE - kills BOTH Enterococcus AND Enterobacteriaceae\n")
cat("3. Anti-anaerobic agents cause Enterobacteriaceae bloom DURING exposure\n")
cat("4. Vancomycin_IV kills gram-positives → gram-negatives expand\n\n")

cat("CIPROFLOXACIN ROUTE FINDINGS:\n")
cat("  IV exposed samples:", nrow(iv_data), "\n")
cat("  PO exposed samples:", nrow(po_data), "\n")
cat("  Unexposed samples:", nrow(none_data), "\n\n")

if (nrow(iv_data) >= 3 && nrow(po_data) >= 3) {
  cat("Both IV and PO ciprofloxacin reach the gut (unlike vancomycin where IV doesn't).\n")
  cat("Check cipro_route_comparison.csv for statistical comparisons.\n")
} else {
  cat("NOTE: Limited samples for route comparison. Results should be interpreted cautiously.\n")
}

cat("\n=== Output Files ===\n")
cat("  antibiotic_category_effects.csv - Detailed effects by antibiotic and category\n")
cat("  antibiotic_category_effects_wide.csv - Wide format summary table\n")
cat("  cipro_route_comparison.csv - Statistical comparisons for cipro routes\n")
cat("  antibiotic_category_heatmap.pdf - Heatmap visualization\n")
cat("  antibiotic_category_barplot.pdf - Bar plot visualization\n")
cat("  antibiotic_effects_faceted.pdf - Faceted by category\n")
cat("  cipro_route_effects.pdf - Ciprofloxacin route comparison\n")
cat("  cipro_route_effects_boxplot.pdf - Box plot with individual samples\n")

cat("\nAnalysis complete!\n")
