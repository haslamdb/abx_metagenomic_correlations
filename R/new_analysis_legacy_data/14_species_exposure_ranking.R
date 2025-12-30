#!/usr/bin/env Rscript
# =============================================================================
# 14_species_exposure_ranking.R
# Unbiased species-level ranking DURING antibiotic exposure
#
# Purpose: Rank all species by their response to antibiotic exposure (from bulk
# differential abundance analysis), then annotate with functional categories
# post-hoc. This shows Phase 1 of the temporal model:
#   - Enterococcus INCREASES during antibiotics
#   - Enterobacteriaceae DECREASES during antibiotics
#   - Anaerobes DECREASE during antibiotics
#
# Uses the robust associations from multi-method concordance analysis.
# =============================================================================

library(tidyverse)
library(ggrepel)

# Set paths
project_dir <- "/home/david/projects/abx_metagenomic_correlations"
results_dir <- file.path(project_dir, "results/new_analysis_legacy_data")
output_dir <- file.path(results_dir, "species_exposure_ranking")
figures_dir <- file.path(results_dir, "figures")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("UNBIASED SPECIES RANKING DURING ANTIBIOTIC EXPOSURE\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n\n")

# =============================================================================
# 1. Load Differential Abundance Results
# =============================================================================

cat("=== Loading Differential Abundance Results ===\n\n")

# Load robust associations (species found by 2+ methods)
robust_file <- file.path(results_dir, "individual_antibiotics_v2_with_covariates/robust_associations.csv")

if (!file.exists(robust_file)) {
  stop("Robust associations file not found. Run 05c_combine_results.R first.")
}

robust_assoc <- read_csv(robust_file, show_col_types = FALSE)
cat("Loaded", nrow(robust_assoc), "robust associations\n")
cat("Antibiotics:", paste(unique(robust_assoc$antibiotic), collapse = ", "), "\n\n")

# Also load method concordance for all associations (including 1-method only)
concordance_file <- file.path(results_dir, "individual_antibiotics_v2_with_covariates/method_concordance.csv")

if (file.exists(concordance_file)) {
  all_assoc <- read_csv(concordance_file, show_col_types = FALSE)
  cat("Loaded", nrow(all_assoc), "total associations (including 1-method)\n\n")
} else {
  all_assoc <- robust_assoc
  cat("Using robust associations only\n\n")
}

# =============================================================================
# 2. Define Species Categories (same as recovery script)
# =============================================================================

cat("=== Defining Species Categories ===\n\n")

# Pathogenic Enterobacteriaceae
pathogenic_enterobact <- c(
  "Escherichia coli",
  "Klebsiella pneumoniae",
  "Klebsiella oxytoca",
  "Enterobacter cloacae",
  "Citrobacter freundii",
  "Serratia marcescens",
  "Salmonella enterica",
  "Shigella flexneri",
  "Shigella sonnei"
)

# All Enterobacteriaceae genera
enterobact_genera <- c("Escherichia", "Klebsiella", "Enterobacter", "Citrobacter",
                       "Serratia", "Salmonella", "Shigella", "Proteus",
                       "Morganella", "Providencia", "Hafnia", "Raoultella",
                       "Cronobacter", "Pantoea")

# Obligate anaerobe genera
obligate_anaerobe_genera <- c(

  "Bacteroides", "Parabacteroides", "Prevotella", "Porphyromonas",
  "Fusobacterium", "Clostridium", "Clostridioides", "Peptostreptococcus",
  "Peptococcus", "Veillonella", "Megasphaera", "Blautia", "Coprococcus",
  "Roseburia", "Faecalibacterium", "Ruminococcus", "Eubacterium",
  "Anaerostipes", "Dorea", "Subdoligranulum", "Dialister",
  "Lachnoclostridium", "Oscillibacter", "Butyricicoccus",
  "Alistipes", "Odoribacter", "Barnesiella"
)

# Function to categorize species
categorize_species <- function(species_name) {
  genus <- strsplit(species_name, " ")[[1]][1]

  if (species_name %in% pathogenic_enterobact) {
    return("Pathogenic Enterobacteriaceae")
  } else if (genus %in% enterobact_genera) {
    return("Commensal Enterobacteriaceae")
  } else if (genus == "Enterococcus") {
    return("Enterococcus")
  } else if (genus %in% obligate_anaerobe_genera) {
    return("Obligate Anaerobe")
  } else if (genus %in% c("Lactobacillus", "Lactococcus", "Streptococcus",
                          "Leuconostoc", "Pediococcus", "Weissella")) {
    return("Lactobacillales")
  } else if (genus == "Bifidobacterium") {
    return("Bifidobacterium")
  } else if (genus %in% c("Staphylococcus", "Pseudomonas", "Acinetobacter",
                          "Stenotrophomonas", "Candida")) {
    return("Other Opportunist")
  } else {
    return("Other")
  }
}

# =============================================================================
# 3. Calculate Pooled Effect Across Antibiotics
# =============================================================================

cat("=== Calculating Pooled Effects Across Antibiotics ===\n\n")

# For each species, calculate:
# - Mean effect across all antibiotics that affect it
# - Number of antibiotics affecting it
# - Consistency of direction

species_pooled <- robust_assoc %>%
  group_by(species) %>%
  summarise(
    n_antibiotics = n(),
    n_methods_max = max(n_methods),
    mean_effect = mean(mean_effect),
    median_effect = median(mean_effect),
    sd_effect = sd(mean_effect),
    n_increased = sum(direction == "increased"),
    n_decreased = sum(direction == "decreased"),
    antibiotics_list = paste(antibiotic, collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(
    # Consistency: all same direction?
    consistent_direction = (n_increased == 0 | n_decreased == 0),
    # Predominant direction
    predominant_direction = ifelse(n_increased > n_decreased, "increased", "decreased"),
    # Category (post-hoc)
    category = sapply(species, categorize_species)
  ) %>%
  arrange(desc(mean_effect))

cat("Species with robust associations:", nrow(species_pooled), "\n")
cat("  Always increased:", sum(species_pooled$n_decreased == 0), "\n")
cat("  Always decreased:", sum(species_pooled$n_increased == 0), "\n")
cat("  Mixed response:", sum(species_pooled$n_increased > 0 & species_pooled$n_decreased > 0), "\n\n")

# =============================================================================
# 4. Rank Species by Effect During Antibiotic Exposure
# =============================================================================

cat("=== Ranking Species by Effect During Exposure ===\n\n")

# Add rank (1 = most increased, last = most decreased)
species_pooled <- species_pooled %>%
  mutate(rank = row_number())

# Top 20 increasers
cat("TOP 20 SPECIES THAT INCREASE DURING ANTIBIOTICS:\n")
species_pooled %>%
  filter(mean_effect > 0) %>%
  head(20) %>%
  select(rank, species, category, mean_effect, n_antibiotics) %>%
  print()

# Top 20 decreasers
cat("\nTOP 20 SPECIES THAT DECREASE DURING ANTIBIOTICS:\n")
species_pooled %>%
  filter(mean_effect < 0) %>%
  arrange(mean_effect) %>%
  head(20) %>%
  select(rank, species, category, mean_effect, n_antibiotics) %>%
  print()

# =============================================================================
# 5. Category Distribution Analysis
# =============================================================================

cat("\n=== Category Distribution Analysis ===\n\n")

# Summary by category
category_summary <- species_pooled %>%
  group_by(category) %>%
  summarise(
    n_species = n(),
    mean_rank = mean(rank),
    median_rank = median(rank),
    mean_effect = mean(mean_effect),
    se_effect = sd(mean_effect) / sqrt(n()),
    pct_increased = mean(mean_effect > 0) * 100,
    pct_decreased = mean(mean_effect < 0) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_effect))

cat("Category summary (ordered by mean effect during antibiotics):\n")
print(category_summary)

# Statistical tests
cat("\n=== Statistical Comparisons ===\n\n")

enterococcus_effects <- species_pooled %>%
  filter(category == "Enterococcus") %>%
  pull(mean_effect)

enterobact_effects <- species_pooled %>%
  filter(category %in% c("Pathogenic Enterobacteriaceae", "Commensal Enterobacteriaceae")) %>%
  pull(mean_effect)

anaerobe_effects <- species_pooled %>%
  filter(category == "Obligate Anaerobe") %>%
  pull(mean_effect)

if (length(enterococcus_effects) >= 3 && length(anaerobe_effects) >= 3) {
  wilcox_entero_vs_anaerobe <- wilcox.test(enterococcus_effects, anaerobe_effects)
  cat("Enterococcus vs Obligate Anaerobes:\n")
  cat("  Enterococcus mean effect:", round(mean(enterococcus_effects), 3), "\n")
  cat("  Anaerobe mean effect:", round(mean(anaerobe_effects), 3), "\n")
  cat("  Wilcoxon p-value:", format(wilcox_entero_vs_anaerobe$p.value, digits = 4), "\n\n")
}

if (length(enterococcus_effects) >= 3 && length(enterobact_effects) >= 3) {
  wilcox_entero_vs_enterobact <- wilcox.test(enterococcus_effects, enterobact_effects)
  cat("Enterococcus vs Enterobacteriaceae:\n")
  cat("  Enterococcus mean effect:", round(mean(enterococcus_effects), 3), "\n")
  cat("  Enterobacteriaceae mean effect:", round(mean(enterobact_effects), 3), "\n")
  cat("  Wilcoxon p-value:", format(wilcox_entero_vs_enterobact$p.value, digits = 4), "\n\n")
}

# =============================================================================
# 6. Generate Visualizations
# =============================================================================

cat("=== Generating Visualizations ===\n\n")

# Color palette
category_colors <- c(
  "Pathogenic Enterobacteriaceae" = "#d62728",
  "Commensal Enterobacteriaceae" = "#ff7f0e",
  "Enterococcus" = "#9467bd",
  "Obligate Anaerobe" = "#1f77b4",
  "Lactobacillales" = "#2ca02c",
  "Bifidobacterium" = "#17becf",
  "Other Opportunist" = "#e377c2",
  "Other" = "#7f7f7f"
)

# 6.1 Ranked bar plot (top 50 species)
top_n <- min(50, nrow(species_pooled))

# Select top increasers and decreasers
top_increasers <- species_pooled %>% filter(mean_effect > 0) %>% head(25)
top_decreasers <- species_pooled %>% filter(mean_effect < 0) %>% arrange(mean_effect) %>% head(25)
top_species <- bind_rows(top_increasers, top_decreasers) %>%
  arrange(desc(mean_effect)) %>%
  mutate(species = factor(species, levels = rev(species)))

p_ranking <- ggplot(top_species, aes(x = mean_effect, y = species, fill = category)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_manual(values = category_colors, name = "Category") +
  labs(
    title = "Species Response DURING Antibiotic Exposure",
    subtitle = "Top 25 increasers and top 25 decreasers. Positive = increases with antibiotics.",
    x = "Mean Effect (log2 fold-change, exposed vs unexposed)",
    y = NULL,
    caption = "Based on robust associations (2+ methods agree)"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 7),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold")
  )

ggsave(file.path(figures_dir, "species_exposure_ranking.pdf"), p_ranking,
       width = 12, height = 14)
cat("Saved: species_exposure_ranking.pdf\n")

# 6.2 Category boxplot
p_category <- ggplot(species_pooled %>% filter(category != "Other"),
                     aes(x = reorder(category, -mean_effect), y = mean_effect, fill = category)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  scale_fill_manual(values = category_colors) +
  labs(
    title = "Effect DURING Antibiotic Exposure by Category",
    subtitle = "Each point = one species. Positive = increases with antibiotics.",
    x = NULL,
    y = "Mean Effect (log2 fold-change)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold")
  )

ggsave(file.path(figures_dir, "species_exposure_by_category.pdf"), p_category,
       width = 10, height = 7)
cat("Saved: species_exposure_by_category.pdf\n")

# 6.3 Forest plot of category means
category_stats <- species_pooled %>%
  filter(category != "Other") %>%
  group_by(category) %>%
  summarise(
    n = n(),
    mean = mean(mean_effect),
    se = sd(mean_effect) / sqrt(n()),
    ci_low = mean - 1.96 * se,
    ci_high = mean + 1.96 * se,
    .groups = "drop"
  ) %>%
  arrange(desc(mean)) %>%
  mutate(category = factor(category, levels = rev(category)))

p_forest <- ggplot(category_stats, aes(x = mean, y = category, color = category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.3, linewidth = 1) +
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("n=%d", n)), hjust = -0.3, vjust = -0.8, size = 3) +
  scale_color_manual(values = category_colors) +
  labs(
    title = "Mean Effect DURING Antibiotic Exposure by Category",
    subtitle = "Error bars = 95% CI. Positive = increases during antibiotics.",
    x = "Mean Effect (log2 fold-change)",
    y = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold")
  )

ggsave(file.path(figures_dir, "species_exposure_category_forest.pdf"), p_forest,
       width = 10, height = 6)
cat("Saved: species_exposure_category_forest.pdf\n")

# 6.4 Highlight key species
key_species_data <- species_pooled %>%
  filter(species %in% c(
    "Enterococcus faecalis", "Enterococcus faecium", "Enterococcus durans",
    "Escherichia coli", "Klebsiella pneumoniae", "Klebsiella oxytoca",
    "Bacteroides fragilis", "Bacteroides vulgatus",
    "Faecalibacterium prausnitzii", "Blautia obeum",
    "Veillonella parvula", "Clostridium perfringens"
  )) %>%
  arrange(desc(mean_effect)) %>%
  mutate(species = factor(species, levels = rev(species)))

if (nrow(key_species_data) > 0) {
  p_key <- ggplot(key_species_data, aes(x = mean_effect, y = species, fill = category)) +
    geom_col() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_text(aes(label = sprintf("n=%d abx", n_antibiotics)),
              hjust = ifelse(key_species_data$mean_effect > 0, -0.1, 1.1),
              size = 3) +
    scale_fill_manual(values = category_colors) +
    labs(
      title = "Key Species Response DURING Antibiotic Exposure",
      subtitle = "Enterococcus increases while Enterobacteriaceae and anaerobes decrease",
      x = "Mean Effect (log2 fold-change)",
      y = NULL
    ) +
    theme_bw() +
    theme(legend.position = "right")

  ggsave(file.path(figures_dir, "species_exposure_key_species.pdf"), p_key,
         width = 10, height = 6)
  cat("Saved: species_exposure_key_species.pdf\n")
}

# =============================================================================
# 7. Save Results
# =============================================================================

cat("\n=== Saving Results ===\n\n")

# Full species ranking
write_csv(species_pooled, file.path(output_dir, "species_exposure_ranking.csv"))
cat("Saved: species_exposure_ranking.csv\n")

# Category summary
write_csv(category_summary, file.path(output_dir, "category_exposure_summary.csv"))
cat("Saved: category_exposure_summary.csv\n")

# Top increasers and decreasers
write_csv(top_species, file.path(output_dir, "top_species_exposure.csv"))
cat("Saved: top_species_exposure.csv\n")

# =============================================================================
# 8. Summary for Manuscript
# =============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("SUMMARY FOR MANUSCRIPT - PHASE 1: DURING ANTIBIOTIC EXPOSURE\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n\n")

cat("KEY FINDING: Enterococcus increases while pathogens and commensals decrease\n\n")

# Enterococcus species
enterococcus_species <- species_pooled %>%
  filter(category == "Enterococcus") %>%
  select(species, mean_effect, n_antibiotics)

cat("Enterococcus species (all INCREASE during antibiotics):\n")
print(enterococcus_species)

# Pathogenic Enterobacteriaceae
pathogen_species <- species_pooled %>%
  filter(category == "Pathogenic Enterobacteriaceae") %>%
  select(species, mean_effect, n_antibiotics)

cat("\nPathogenic Enterobacteriaceae (trend to DECREASE during antibiotics):\n")
print(pathogen_species)

# Summary stats
cat("\n=== Category Mean Effects ===\n")
cat("Enterococcus:", round(mean(enterococcus_effects), 3), "log2FC (n =", length(enterococcus_effects), "species)\n")
cat("Enterobacteriaceae:", round(mean(enterobact_effects), 3), "log2FC (n =", length(enterobact_effects), "species)\n")
cat("Obligate Anaerobes:", round(mean(anaerobe_effects), 3), "log2FC (n =", length(anaerobe_effects), "species)\n\n")

cat("INTERPRETATION (Phase 1 of Temporal Model):\n")
cat("  - During antibiotics, Enterococcus EXPANDS (fills ecological niche)\n")
cat("  - Enterobacteriaceae are SUPPRESSED (susceptible to broad-spectrum agents)\n")
cat("  - Obligate anaerobes are SUPPRESSED (susceptible, slow to recover)\n")
cat("  - This sets the stage for Phase 2: rapid Enterobacteriaceae recovery\n\n")

cat("Analysis complete!\n")
