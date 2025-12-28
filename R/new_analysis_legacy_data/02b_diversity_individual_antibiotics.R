#!/usr/bin/env Rscript
# =============================================================================
# 02b_diversity_individual_antibiotics.R
# Alpha/beta diversity analysis for individual antibiotics with covariate adjustment
# Models each antibiotic's effect while controlling for other antibiotics
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(tibble)
library(stringr)
library(vegan)
library(lme4)
library(lmerTest)
library(broom.mixed)

# =============================================================================
# Setup
# =============================================================================

project_dir <- "/home/david/projects/abx_metagenomic_correlations"
results_dir <- file.path(project_dir, "results/new_analysis_legacy_data")
figures_dir <- file.path(results_dir, "figures")
diversity_dir <- file.path(results_dir, "diversity")

dir.create(diversity_dir, showWarnings = FALSE, recursive = TRUE)

# Load prepared data
load(file.path(project_dir, "data/prepared_data.RData"))
load(file.path(project_dir, "data/legacy/AbxEffectData20191014"))

sample_metadata <- prepared_data$sample_metadata
species_matrix <- prepared_data$species_matrix

# =============================================================================
# 1. Calculate Individual Antibiotic Exposures (same as script 05)
# =============================================================================

cat("=== Calculating Individual Antibiotic Exposures ===\n\n")

top_antibiotics <- c("Pip/Tazo", "TMP/SMX", "Vancomycin", "Cefepime",
                     "Meropenem", "Ciprofloxacin", "Metronidazole",
                     "Azithromycin", "Cefazolin", "Ceftriaxone")

# Filter drug table to systemic routes
systemic_routes <- c("IV", "PO", "IM")
DrugTable_systemic <- DrugTable %>%
  filter(Route %in% systemic_routes) %>%
  mutate(Date = as.Date(Date))

# Function to calculate exposure
calc_abx_exposure <- function(mrn, sample_date, drug_name, drug_data, window = 7) {
  start_date <- sample_date - window
  end_date <- sample_date - 1

  exposed <- drug_data %>%
    filter(MRN == mrn,
           Drug == drug_name,
           Date >= start_date,
           Date <= end_date) %>%
    nrow() > 0

  return(exposed)
}

# Calculate exposure for each antibiotic
for (abx in top_antibiotics) {
  col_name <- paste0(gsub("/", "_", gsub("-", "_", abx)), "_7d")
  cat("  Calculating", abx, "exposure...\n")

  sample_metadata[[col_name]] <- sapply(1:nrow(sample_metadata), function(i) {
    calc_abx_exposure(
      sample_metadata$MRN[i],
      sample_metadata$SampleDate[i],
      abx,
      DrugTable_systemic
    )
  })
}

# Create clean column names mapping
abx_cols <- sapply(top_antibiotics, function(abx) {
  paste0(gsub("/", "_", gsub("-", "_", abx)), "_7d")
})
names(abx_cols) <- top_antibiotics

# Print exposure summary
cat("\n=== Antibiotic Exposure Summary ===\n\n")
exposure_summary <- data.frame(
  Antibiotic = top_antibiotics,
  Exposed = sapply(top_antibiotics, function(abx) sum(sample_metadata[[abx_cols[abx]]])),
  Total = nrow(sample_metadata)
) %>%
  mutate(Pct = round(100 * Exposed / Total, 1))
print(exposure_summary)

# =============================================================================
# 2. Prepare High-Confidence Samples
# =============================================================================

cat("\n=== Preparing High-Confidence Samples ===\n\n")

sample_hc <- sample_metadata %>% filter(high_confidence)
cat("High-confidence samples:", nrow(sample_hc), "\n")
cat("Unique patients:", n_distinct(sample_hc$MRN), "\n")

# Ensure patient_group exists
if (!"patient_group" %in% colnames(sample_hc)) {
  sample_hc$patient_group <- sample_hc$PatientGroup
}

cat("\nPatient group distribution:\n")
print(table(sample_hc$patient_group))

# =============================================================================
# 3. Alpha Diversity Models - Individual Antibiotics with Covariates
# =============================================================================

cat("\n=== Alpha Diversity: Individual Antibiotic Effects ===\n")
cat("Adjusting for: other antibiotics + patient group + random patient effect\n\n")

alpha_results <- list()

for (abx in top_antibiotics) {
  target_col <- abx_cols[abx]

  # Check sample sizes
  n_exposed <- sum(sample_hc[[target_col]])
  n_unexposed <- sum(!sample_hc[[target_col]])

  cat("Processing:", abx, "(n_exposed =", n_exposed, ", n_unexposed =", n_unexposed, ")\n")

  if (n_exposed < 10 || n_unexposed < 10) {
    cat("  Skipping - insufficient samples\n")
    next
  }

  # Build formula with all antibiotics as covariates
  other_abx_cols <- abx_cols[names(abx_cols) != abx]

  # Full model: target + other antibiotics + patient group + (1|patient)
  formula_str <- paste0(
    "shannon ~ ", target_col, " + patient_group + ",
    paste(other_abx_cols, collapse = " + "),
    " + (1|MRN)"
  )

  tryCatch({
    model <- lmer(as.formula(formula_str), data = sample_hc)

    # Extract target antibiotic effect
    coefs <- tidy(model, effects = "fixed", conf.int = TRUE)
    target_result <- coefs %>%
      filter(term == paste0(target_col, "TRUE")) %>%
      mutate(antibiotic = abx)

    if (nrow(target_result) > 0) {
      alpha_results[[abx]] <- target_result
      cat("  Effect:", round(target_result$estimate, 3),
          "(95% CI:", round(target_result$conf.low, 3), "-",
          round(target_result$conf.high, 3), ")",
          "p =", format(target_result$p.value, digits = 3), "\n")
    }
  }, error = function(e) {
    cat("  Error:", conditionMessage(e), "\n")
  })
}

# Combine results
alpha_combined <- bind_rows(alpha_results) %>%
  arrange(estimate) %>%
  mutate(
    significant = p.value < 0.05,
    direction = ifelse(estimate > 0, "Increases diversity", "Decreases diversity")
  )

cat("\n=== Shannon Diversity: Summary of Individual Antibiotic Effects ===\n\n")
print(alpha_combined %>%
        select(antibiotic, estimate, conf.low, conf.high, p.value, significant) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))))

# =============================================================================
# 4. Functional Group Models - Anaerobes, Enterobacteriaceae, Enterococcus
# =============================================================================

cat("\n=== Functional Groups: Individual Antibiotic Effects ===\n\n")

functional_groups <- list(
  anaerobes = list(col = "anaerobes_rel", name = "Obligate Anaerobes"),
  enterobact = list(col = "enterobact_rel", name = "Enterobacteriaceae"),
  enterococcus = list(col = "enterococcus_rel", name = "Enterococcus")
)

functional_results <- list()

for (fg_key in names(functional_groups)) {
  fg <- functional_groups[[fg_key]]
  cat("\n--- Modeling:", fg$name, "---\n")

  fg_results <- list()

  for (abx in top_antibiotics) {
    target_col <- abx_cols[abx]

    n_exposed <- sum(sample_hc[[target_col]])
    n_unexposed <- sum(!sample_hc[[target_col]])

    if (n_exposed < 10 || n_unexposed < 10) next

    other_abx_cols <- abx_cols[names(abx_cols) != abx]

    # Model log-transformed relative abundance
    formula_str <- paste0(
      "log10(", fg$col, " + 1e-6) ~ ", target_col, " + patient_group + ",
      paste(other_abx_cols, collapse = " + "),
      " + (1|MRN)"
    )

    tryCatch({
      model <- lmer(as.formula(formula_str), data = sample_hc)

      coefs <- tidy(model, effects = "fixed", conf.int = TRUE)
      target_result <- coefs %>%
        filter(term == paste0(target_col, "TRUE")) %>%
        mutate(
          antibiotic = abx,
          functional_group = fg$name
        )

      if (nrow(target_result) > 0) {
        fg_results[[abx]] <- target_result
      }
    }, error = function(e) {
      # Silent skip
    })
  }

  if (length(fg_results) > 0) {
    functional_results[[fg_key]] <- bind_rows(fg_results)

    cat("\nResults for", fg$name, ":\n")
    print(functional_results[[fg_key]] %>%
            select(antibiotic, estimate, conf.low, conf.high, p.value) %>%
            mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
            arrange(estimate))
  }
}

# Combine all functional group results
functional_combined <- bind_rows(functional_results) %>%
  mutate(significant = p.value < 0.05)

# =============================================================================
# 5. Beta Diversity - PERMANOVA with Individual Antibiotics
# =============================================================================

cat("\n=== Beta Diversity: PERMANOVA for Individual Antibiotics ===\n")
cat("Adjusting for: other antibiotics + patient group\n\n")

# Prepare species data
hc_samples <- sample_hc$sample_id
species_hc <- species_matrix[intersect(rownames(species_matrix), hc_samples), ]
species_rel_hc <- species_hc / rowSums(species_hc)

# Align metadata
meta_aligned <- sample_hc[match(rownames(species_rel_hc), sample_hc$sample_id), ]
rownames(meta_aligned) <- meta_aligned$sample_id

# Calculate Bray-Curtis
bray_dist_hc <- vegdist(species_rel_hc, method = "bray")

permanova_results <- list()

for (abx in top_antibiotics) {
  target_col <- abx_cols[abx]

  n_exposed <- sum(meta_aligned[[target_col]])
  n_unexposed <- sum(!meta_aligned[[target_col]])

  cat("Processing:", abx, "(n_exposed =", n_exposed, ")\n")

  if (n_exposed < 10 || n_unexposed < 10) {
    cat("  Skipping - insufficient samples\n")
    next
  }

  other_abx_cols <- abx_cols[names(abx_cols) != abx]

  # Build PERMANOVA formula
  formula_str <- paste0(
    "bray_dist_hc ~ ", target_col, " + patient_group + ",
    paste(other_abx_cols, collapse = " + ")
  )

  tryCatch({
    set.seed(42)
    perm_result <- adonis2(
      as.formula(formula_str),
      data = meta_aligned,
      permutations = 999,
      by = "margin"
    )

    # Extract target antibiotic result
    target_row <- perm_result[target_col, ]

    permanova_results[[abx]] <- data.frame(
      antibiotic = abx,
      R2 = target_row$R2,
      F_stat = target_row$F,
      p_value = target_row$`Pr(>F)`
    )

    cat("  R² =", round(target_row$R2, 4),
        ", F =", round(target_row$F, 2),
        ", p =", format(target_row$`Pr(>F)`, digits = 3), "\n")

  }, error = function(e) {
    cat("  Error:", conditionMessage(e), "\n")
  })
}

permanova_combined <- bind_rows(permanova_results) %>%
  arrange(desc(R2)) %>%
  mutate(significant = p_value < 0.05)

cat("\n=== PERMANOVA Summary ===\n\n")
print(permanova_combined %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))))

# =============================================================================
# 6. Generate Forest Plots
# =============================================================================

cat("\n=== Generating Forest Plots ===\n\n")

# --- Forest plot for Shannon diversity ---

alpha_plot_data <- alpha_combined %>%
  mutate(antibiotic = factor(antibiotic, levels = antibiotic[order(estimate)]))

p_alpha_forest <- ggplot(alpha_plot_data, aes(x = estimate, y = antibiotic)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, color = "gray30") +
  geom_point(aes(color = significant, size = abs(estimate)), shape = 16) +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "#E74C3C"),
                     labels = c("FALSE" = "p ≥ 0.05", "TRUE" = "p < 0.05")) +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  labs(
    title = "Effect of Individual Antibiotics on Shannon Diversity",
    subtitle = "Adjusted for other antibiotics, patient group, and patient random effect",
    x = "Effect on Shannon Diversity (95% CI)",
    y = NULL,
    color = "Significance"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10)
  )

ggsave(file.path(figures_dir, "shannon_by_individual_abx_forest.pdf"),
       p_alpha_forest, width = 8, height = 6)
cat("Saved: figures/shannon_by_individual_abx_forest.pdf\n")

# --- Forest plot for functional groups ---

if (nrow(functional_combined) > 0) {
  functional_plot_data <- functional_combined %>%
    mutate(
      label = paste0(antibiotic, " → ", functional_group),
      label = factor(label, levels = label[order(functional_group, estimate)])
    )

  p_functional_forest <- ggplot(functional_plot_data,
                                 aes(x = estimate, y = label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0.2, color = "gray30") +
    geom_point(aes(color = significant), size = 3, shape = 16) +
    facet_wrap(~ functional_group, scales = "free_y", ncol = 1) +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "#E74C3C")) +
    labs(
      title = "Effect of Individual Antibiotics on Functional Groups",
      subtitle = "Adjusted for other antibiotics + patient group (log10 relative abundance)",
      x = "Effect (log10 scale, 95% CI)",
      y = NULL,
      color = "p < 0.05"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "gray90")
    )

  ggsave(file.path(figures_dir, "functional_groups_by_individual_abx_forest.pdf"),
         p_functional_forest, width = 10, height = 12)
  cat("Saved: figures/functional_groups_by_individual_abx_forest.pdf\n")
}

# --- Bar plot for PERMANOVA R² ---

p_permanova <- ggplot(permanova_combined,
                       aes(x = reorder(antibiotic, R2), y = R2, fill = significant)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = ifelse(significant, "*", "")),
            vjust = -0.5, size = 6) +
  scale_fill_manual(values = c("FALSE" = "gray60", "TRUE" = "#3498DB"),
                    labels = c("FALSE" = "p ≥ 0.05", "TRUE" = "p < 0.05")) +
  labs(
    title = "Variance in Microbiome Composition Explained by Individual Antibiotics",
    subtitle = "PERMANOVA R² (adjusted for other antibiotics + patient group)",
    x = NULL,
    y = "R² (variance explained)",
    fill = "Significance"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(file.path(figures_dir, "permanova_by_individual_abx.pdf"),
       p_permanova, width = 8, height = 5)
cat("Saved: figures/permanova_by_individual_abx.pdf\n")

# =============================================================================
# 7. Save Results
# =============================================================================

cat("\n=== Saving Results ===\n\n")

# Alpha diversity results
write_csv(alpha_combined,
          file.path(diversity_dir, "alpha_diversity_individual_abx.csv"))
cat("Saved: diversity/alpha_diversity_individual_abx.csv\n")

# Functional group results
write_csv(functional_combined,
          file.path(diversity_dir, "functional_groups_individual_abx.csv"))
cat("Saved: diversity/functional_groups_individual_abx.csv\n")

# PERMANOVA results
write_csv(permanova_combined,
          file.path(diversity_dir, "permanova_individual_abx.csv"))
cat("Saved: diversity/permanova_individual_abx.csv\n")

# Combined summary table
summary_table <- alpha_combined %>%
  select(antibiotic, shannon_effect = estimate, shannon_p = p.value) %>%
  left_join(
    permanova_combined %>%
      select(antibiotic, permanova_R2 = R2, permanova_p = p_value),
    by = "antibiotic"
  ) %>%
  left_join(
    functional_combined %>%
      filter(functional_group == "Obligate Anaerobes") %>%
      select(antibiotic, anaerobes_effect = estimate, anaerobes_p = p.value),
    by = "antibiotic"
  ) %>%
  left_join(
    functional_combined %>%
      filter(functional_group == "Enterobacteriaceae") %>%
      select(antibiotic, enterobact_effect = estimate, enterobact_p = p.value),
    by = "antibiotic"
  ) %>%
  left_join(
    functional_combined %>%
      filter(functional_group == "Enterococcus") %>%
      select(antibiotic, enterococcus_effect = estimate, enterococcus_p = p.value),
    by = "antibiotic"
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

write_csv(summary_table,
          file.path(diversity_dir, "individual_abx_effects_summary.csv"))
cat("Saved: diversity/individual_abx_effects_summary.csv\n")

cat("\n=== Analysis Complete ===\n")
cat("\nKey outputs:\n")
cat("  - alpha_diversity_individual_abx.csv: Shannon diversity effects\n")
cat("  - functional_groups_individual_abx.csv: Anaerobes/Enterobact/Enterococcus effects\n")
cat("  - permanova_individual_abx.csv: Beta diversity (composition) effects\n")
cat("  - individual_abx_effects_summary.csv: Combined summary table\n")
cat("  - Forest plots in figures/\n")
