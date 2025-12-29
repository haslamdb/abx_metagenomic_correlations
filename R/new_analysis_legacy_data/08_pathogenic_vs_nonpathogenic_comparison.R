#!/usr/bin/env Rscript
# =============================================================================
# 08_pathogenic_vs_nonpathogenic_comparison.R
# Compare antibiotic susceptibility between pathogenic and non-pathogenic
# species within the same genus (Escherichia, Klebsiella)
#
# Hypothesis: Pathogenic species may be more robust to antibiotic-induced
# microbiome perturbation than their non-pathogenic relatives
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

select <- dplyr::select
filter <- dplyr::filter

project_dir <- "/home/david/projects/abx_metagenomic_correlations"
results_dir <- file.path(project_dir, "results/new_analysis_legacy_data")

# Output directory
compare_dir <- file.path(results_dir, "pathogenic_vs_nonpathogenic")
dir.create(compare_dir, showWarnings = FALSE, recursive = TRUE)

cat("=============================================================================\n")
cat("Pathogenic vs Non-Pathogenic Species Comparison\n")
cat("=============================================================================\n\n")

# =============================================================================
# 1. Load Data
# =============================================================================

cat("Loading data...\n")

# Differential abundance results
aldex <- read_csv(file.path(results_dir, "individual_antibiotics_v2_with_covariates/aldex2/all_antibiotics_combined.csv"),
                  show_col_types = FALSE)
maaslin <- read_csv(file.path(results_dir, "individual_antibiotics_v2_with_covariates/maaslin3/all_antibiotics_combined.csv"),
                    show_col_types = FALSE)
ancombc <- read_csv(file.path(results_dir, "individual_antibiotics_v2_with_covariates/ancombc2/all_antibiotics_combined.csv"),
                    show_col_types = FALSE)

# Abundance data
load(file.path(project_dir, "data/bracken_count_matrices.RData"))
species_matrix <- bracken_data$species_matrix

cat("Loaded", ncol(species_matrix), "species from", nrow(species_matrix), "samples\n\n")

# =============================================================================
# 2. Define Analysis Function
# =============================================================================

analyze_species <- function(species_name, aldex, maaslin, ancombc, species_matrix) {

  # Check if species exists in data
  if (!species_name %in% colnames(species_matrix)) {
    return(NULL)
  }

  # Calculate prevalence and abundance
  counts <- species_matrix[, species_name]
  prevalence <- sum(counts > 0) / length(counts) * 100
  total_reads <- rowSums(species_matrix)
  rel_abund <- counts / total_reads * 100
  mean_abund <- if (sum(counts > 0) > 0) mean(rel_abund[counts > 0]) else 0

  # Get significant results from each method
  aldex_sig <- aldex %>% filter(species == species_name, pval_BH < 0.1)
  maaslin_sig <- maaslin %>% filter(species == species_name, pval_BH < 0.1)
  ancombc_sig <- ancombc %>% filter(species == species_name, diff_abundant == TRUE)

  # Count antibiotics with at least 1 method significant
  all_abx_sig <- unique(c(aldex_sig$antibiotic, maaslin_sig$antibiotic, ancombc_sig$antibiotic))

  # Count antibiotics with 2+ methods (robust)
  abx_counts <- table(c(aldex_sig$antibiotic, maaslin_sig$antibiotic, ancombc_sig$antibiotic))
  robust_abx <- names(abx_counts[abx_counts >= 2])

  # Get all ANCOMBC results for mean effect size
  ancombc_all <- ancombc %>% filter(species == species_name)
  mean_abs_effect <- if (nrow(ancombc_all) > 0) mean(abs(ancombc_all$lfc), na.rm = TRUE) else NA

  # Which antibiotics affect this species?
  abx_list <- if (length(all_abx_sig) > 0) paste(all_abx_sig, collapse = ", ") else "none"

  data.frame(
    species = species_name,
    prevalence = round(prevalence, 1),
    mean_abundance = round(mean_abund, 4),
    n_abx_any_sig = length(all_abx_sig),
    n_abx_robust = length(robust_abx),
    mean_abs_lfc = round(mean_abs_effect, 3),
    antibiotics_affected = abx_list,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 3. ESCHERICHIA Comparison
# =============================================================================

cat("=== ESCHERICHIA GENUS ===\n\n")

# Get all Escherichia species (exclude phages/viruses)
all_escherichia <- colnames(species_matrix)[grepl("^Escherichia", colnames(species_matrix))]
escherichia_bacteria <- all_escherichia[!grepl("phage|virus", all_escherichia, ignore.case = TRUE)]

cat("Found", length(escherichia_bacteria), "Escherichia bacterial species/strains\n")

# Analyze each species
esch_results <- lapply(escherichia_bacteria, function(sp) {
  analyze_species(sp, aldex, maaslin, ancombc, species_matrix)
})

esch_df <- bind_rows(esch_results) %>%
  filter(!is.na(species)) %>%
  mutate(
    pathogen_status = case_when(
      species == "Escherichia coli" ~ "PRIMARY_PATHOGEN",
      species == "Escherichia albertii" ~ "emerging_pathogen",
      species == "Escherichia fergusonii" ~ "opportunistic",
      TRUE ~ "commensal"
    )
  ) %>%
  arrange(desc(prevalence))

# Print results (filter to >5% prevalence for readability)
cat("\nEscherichia species with >5% prevalence:\n\n")
esch_prevalent <- esch_df %>% filter(prevalence > 5)
print(as.data.frame(esch_prevalent %>%
  select(species, pathogen_status, prevalence, mean_abundance,
         n_abx_any_sig, n_abx_robust, mean_abs_lfc)))

# Statistical comparison
cat("\n--- ESCHERICHIA STATISTICAL SUMMARY ---\n\n")

ecoli <- esch_df %>% filter(species == "Escherichia coli")
e_fergusonii <- esch_df %>% filter(species == "Escherichia fergusonii")
e_albertii <- esch_df %>% filter(species == "Escherichia albertii")
other_esch <- esch_df %>% filter(pathogen_status == "commensal", prevalence > 10)

cat("E. coli (primary pathogen, prevalence", ecoli$prevalence, "%):\n")
cat("  Antibiotics with ANY effect:", ecoli$n_abx_any_sig, "/9\n")
cat("  Antibiotics with ROBUST effect (2+ methods):", ecoli$n_abx_robust, "/9\n")
cat("  Mean |log fold change|:", ecoli$mean_abs_lfc, "\n")
cat("  Affected by:", ecoli$antibiotics_affected, "\n\n")

if (nrow(e_fergusonii) > 0) {
  cat("E. fergusonii (opportunistic, prevalence", e_fergusonii$prevalence, "%):\n")
  cat("  Antibiotics with ANY effect:", e_fergusonii$n_abx_any_sig, "/9\n")
  cat("  Antibiotics with ROBUST effect:", e_fergusonii$n_abx_robust, "/9\n")
  cat("  Mean |log fold change|:", e_fergusonii$mean_abs_lfc, "\n")
  cat("  Affected by:", e_fergusonii$antibiotics_affected, "\n\n")
}

if (nrow(e_albertii) > 0) {
  cat("E. albertii (emerging pathogen, prevalence", e_albertii$prevalence, "%):\n")
  cat("  Antibiotics with ANY effect:", e_albertii$n_abx_any_sig, "/9\n")
  cat("  Antibiotics with ROBUST effect:", e_albertii$n_abx_robust, "/9\n")
  cat("  Mean |log fold change|:", e_albertii$mean_abs_lfc, "\n")
  cat("  Affected by:", e_albertii$antibiotics_affected, "\n\n")
}

if (nrow(other_esch) > 0) {
  cat("Other commensal Escherichia (>10% prevalence, n =", nrow(other_esch), "):\n")
  cat("  Mean antibiotics with ANY effect:", round(mean(other_esch$n_abx_any_sig), 2), "/9\n")
  cat("  Mean antibiotics with ROBUST effect:", round(mean(other_esch$n_abx_robust), 2), "/9\n")
  cat("  Mean |log fold change|:", round(mean(other_esch$mean_abs_lfc, na.rm = TRUE), 3), "\n\n")

  # Compare E. coli to commensals
  cat("COMPARISON: Is E. coli more resistant to antibiotic perturbation?\n")
  cat("  E. coli affected by", ecoli$n_abx_any_sig, "antibiotics\n")
  cat("  Commensals affected by", round(mean(other_esch$n_abx_any_sig), 1), "antibiotics (mean)\n")
  cat("  E. coli percentile (lower = more resistant):",
      round(sum(ecoli$n_abx_any_sig <= other_esch$n_abx_any_sig) / nrow(other_esch) * 100, 1), "%\n")
}

# =============================================================================
# 4. KLEBSIELLA Comparison
# =============================================================================

cat("\n\n=== KLEBSIELLA GENUS ===\n\n")

# Get all Klebsiella species (exclude phages)
all_klebsiella <- colnames(species_matrix)[grepl("^Klebsiella", colnames(species_matrix))]
klebsiella_bacteria <- all_klebsiella[!grepl("phage|virus", all_klebsiella, ignore.case = TRUE)]

cat("Found", length(klebsiella_bacteria), "Klebsiella bacterial species/strains\n")

# Analyze each
kleb_results <- lapply(klebsiella_bacteria, function(sp) {
  analyze_species(sp, aldex, maaslin, ancombc, species_matrix)
})

kleb_df <- bind_rows(kleb_results) %>%
  filter(!is.na(species)) %>%
  mutate(
    pathogen_status = case_when(
      species == "Klebsiella pneumoniae" ~ "PRIMARY_PATHOGEN",
      species == "Klebsiella oxytoca" ~ "PRIMARY_PATHOGEN",
      species == "Klebsiella aerogenes" ~ "opportunistic",
      species == "Klebsiella variicola" ~ "opportunistic",
      species == "Klebsiella quasipneumoniae" ~ "opportunistic",
      TRUE ~ "commensal"
    )
  ) %>%
  arrange(desc(prevalence))

# Print results
cat("\nKlebsiella species with >5% prevalence:\n\n")
kleb_prevalent <- kleb_df %>% filter(prevalence > 5)
print(as.data.frame(kleb_prevalent %>%
  select(species, pathogen_status, prevalence, mean_abundance,
         n_abx_any_sig, n_abx_robust, mean_abs_lfc)))

# Statistical comparison
cat("\n--- KLEBSIELLA STATISTICAL SUMMARY ---\n\n")

kpn <- kleb_df %>% filter(species == "Klebsiella pneumoniae")
kox <- kleb_df %>% filter(species == "Klebsiella oxytoca")
k_aerogenes <- kleb_df %>% filter(species == "Klebsiella aerogenes")
other_kleb <- kleb_df %>% filter(pathogen_status == "commensal", prevalence > 5)

cat("K. pneumoniae (primary pathogen, prevalence", kpn$prevalence, "%):\n")
cat("  Antibiotics with ANY effect:", kpn$n_abx_any_sig, "/9\n")
cat("  Antibiotics with ROBUST effect:", kpn$n_abx_robust, "/9\n")
cat("  Mean |log fold change|:", kpn$mean_abs_lfc, "\n")
cat("  Affected by:", kpn$antibiotics_affected, "\n\n")

cat("K. oxytoca (primary pathogen, prevalence", kox$prevalence, "%):\n")
cat("  Antibiotics with ANY effect:", kox$n_abx_any_sig, "/9\n")
cat("  Antibiotics with ROBUST effect:", kox$n_abx_robust, "/9\n")
cat("  Mean |log fold change|:", kox$mean_abs_lfc, "\n")
cat("  Affected by:", kox$antibiotics_affected, "\n\n")

if (nrow(k_aerogenes) > 0) {
  cat("K. aerogenes (opportunistic, prevalence", k_aerogenes$prevalence, "%):\n")
  cat("  Antibiotics with ANY effect:", k_aerogenes$n_abx_any_sig, "/9\n")
  cat("  Antibiotics with ROBUST effect:", k_aerogenes$n_abx_robust, "/9\n")
  cat("  Mean |log fold change|:", k_aerogenes$mean_abs_lfc, "\n\n")
}

if (nrow(other_kleb) > 0) {
  cat("Other commensal Klebsiella (>5% prevalence, n =", nrow(other_kleb), "):\n")
  cat("  Mean antibiotics with ANY effect:", round(mean(other_kleb$n_abx_any_sig), 2), "/9\n")
  cat("  Mean antibiotics with ROBUST effect:", round(mean(other_kleb$n_abx_robust), 2), "/9\n")
  cat("  Mean |log fold change|:", round(mean(other_kleb$mean_abs_lfc, na.rm = TRUE), 3), "\n\n")
}

# =============================================================================
# 5. Create Comparison Plot
# =============================================================================

cat("Creating comparison plots...\n")

# Combine data for plotting
plot_data <- bind_rows(
  esch_df %>% mutate(genus = "Escherichia") %>% filter(prevalence > 10),
  kleb_df %>% mutate(genus = "Klebsiella") %>% filter(prevalence > 5)
) %>%
  mutate(
    is_pathogen = pathogen_status %in% c("PRIMARY_PATHOGEN", "emerging_pathogen"),
    pathogen_label = ifelse(is_pathogen, "Pathogenic", "Non-pathogenic")
  )

# Plot: Number of antibiotics affecting each species
p1 <- ggplot(plot_data, aes(x = reorder(species, -n_abx_any_sig), y = n_abx_any_sig,
                             fill = pathogen_label)) +
  geom_col() +
  facet_wrap(~ genus, scales = "free_x") +
  scale_fill_manual(values = c("Pathogenic" = "firebrick", "Non-pathogenic" = "steelblue")) +
  labs(
    title = "Number of Antibiotics Affecting Each Species",
    subtitle = "Lower = more resistant to antibiotic-induced perturbation",
    x = NULL,
    y = "Number of antibiotics with significant effect",
    fill = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom"
  )

ggsave(file.path(compare_dir, "antibiotics_per_species.pdf"), p1,
       width = 12, height = 6)

# Plot: Effect size comparison
p2 <- ggplot(plot_data %>% filter(!is.na(mean_abs_lfc)),
             aes(x = pathogen_label, y = mean_abs_lfc, fill = pathogen_label)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ genus) +
  scale_fill_manual(values = c("Pathogenic" = "firebrick", "Non-pathogenic" = "steelblue")) +
  labs(
    title = "Mean Absolute Effect Size by Pathogen Status",
    subtitle = "Lower = less affected by antibiotics overall",
    x = NULL,
    y = "Mean |log fold change|"
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(compare_dir, "effect_size_by_pathogen_status.pdf"), p2,
       width = 8, height = 5)

# Plot: Scatter of prevalence vs antibiotics affected
p3 <- ggplot(plot_data, aes(x = prevalence, y = n_abx_any_sig,
                             color = pathogen_label, size = mean_abundance)) +
  geom_point(alpha = 0.7) +
  geom_text(aes(label = ifelse(is_pathogen, species, "")),
            hjust = -0.1, vjust = 0, size = 3, show.legend = FALSE) +
  facet_wrap(~ genus, scales = "free") +
  scale_color_manual(values = c("Pathogenic" = "firebrick", "Non-pathogenic" = "steelblue")) +
  labs(
    title = "Prevalence vs Antibiotic Susceptibility",
    x = "Prevalence (%)",
    y = "Number of antibiotics with effect",
    color = NULL,
    size = "Mean abundance (%)"
  ) +
  theme_bw() +
  theme(legend.position = "right")

ggsave(file.path(compare_dir, "prevalence_vs_susceptibility.pdf"), p3,
       width = 12, height = 5)

# =============================================================================
# 6. Summary Statistics and Output
# =============================================================================

cat("\n=== OVERALL SUMMARY ===\n\n")

# Combine for summary
all_species <- bind_rows(
  esch_df %>% mutate(genus = "Escherichia"),
  kleb_df %>% mutate(genus = "Klebsiella")
) %>%
  filter(prevalence > 5)

pathogens <- all_species %>% filter(pathogen_status == "PRIMARY_PATHOGEN")
nonpathogens <- all_species %>% filter(pathogen_status == "commensal")

cat("Primary pathogens (n =", nrow(pathogens), "):\n")
cat("  Mean antibiotics with effect:", round(mean(pathogens$n_abx_any_sig), 2), "\n")
cat("  Mean antibiotics with robust effect:", round(mean(pathogens$n_abx_robust), 2), "\n")
cat("  Mean |effect size|:", round(mean(pathogens$mean_abs_lfc, na.rm = TRUE), 3), "\n\n")

cat("Commensal species (n =", nrow(nonpathogens), "):\n")
cat("  Mean antibiotics with effect:", round(mean(nonpathogens$n_abx_any_sig), 2), "\n")
cat("  Mean antibiotics with robust effect:", round(mean(nonpathogens$n_abx_robust), 2), "\n")
cat("  Mean |effect size|:", round(mean(nonpathogens$mean_abs_lfc, na.rm = TRUE), 3), "\n\n")

# Statistical test
if (nrow(pathogens) >= 2 && nrow(nonpathogens) >= 3) {
  wilcox_n <- wilcox.test(pathogens$n_abx_any_sig, nonpathogens$n_abx_any_sig)
  wilcox_effect <- wilcox.test(pathogens$mean_abs_lfc, nonpathogens$mean_abs_lfc)

  cat("Statistical tests (Wilcoxon rank-sum):\n")
  cat("  N antibiotics affected: p =", round(wilcox_n$p.value, 4), "\n")
  cat("  Mean effect size: p =", round(wilcox_effect$p.value, 4), "\n")
}

# Save results
write_csv(esch_df, file.path(compare_dir, "escherichia_comparison.csv"))
write_csv(kleb_df, file.path(compare_dir, "klebsiella_comparison.csv"))

summary_df <- data.frame(
  category = c("E. coli", "Other Escherichia (>10% prev)",
               "K. pneumoniae", "K. oxytoca", "Other Klebsiella (>5% prev)"),
  n_species = c(1, nrow(other_esch), 1, 1, nrow(other_kleb)),
  mean_prevalence = c(ecoli$prevalence, mean(other_esch$prevalence),
                      kpn$prevalence, kox$prevalence, mean(other_kleb$prevalence)),
  mean_abx_affected = c(ecoli$n_abx_any_sig, mean(other_esch$n_abx_any_sig),
                        kpn$n_abx_any_sig, kox$n_abx_any_sig, mean(other_kleb$n_abx_any_sig)),
  mean_abx_robust = c(ecoli$n_abx_robust, mean(other_esch$n_abx_robust),
                      kpn$n_abx_robust, kox$n_abx_robust, mean(other_kleb$n_abx_robust)),
  mean_effect_size = c(ecoli$mean_abs_lfc, mean(other_esch$mean_abs_lfc, na.rm = TRUE),
                       kpn$mean_abs_lfc, kox$mean_abs_lfc, mean(other_kleb$mean_abs_lfc, na.rm = TRUE))
)

print(summary_df)
write_csv(summary_df, file.path(compare_dir, "summary_comparison.csv"))

cat("\n=============================================================================\n")
cat("Analysis complete! Results saved to:", compare_dir, "\n")
cat("=============================================================================\n")
