#!/usr/bin/env Rscript
# =============================================================================
# 05b_maaslin3_patch.R
# MaAsLin3-only patch run to fill in results from main analysis
# Run concurrently with 05_individual_antibiotic_effects.R
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(stringr)
library(purrr)
library(maaslin3)

# Ensure dplyr functions aren't masked
select <- dplyr::select
filter <- dplyr::filter

# Set paths
project_dir <- "/home/david/projects/abx_metagenomic_correlations"
results_dir <- file.path(project_dir, "results/new_analysis_legacy_data")

# Load data
load(file.path(project_dir, "data/prepared_data.RData"))
load(file.path(project_dir, "data/legacy/AbxEffectData20191014"))

sample_metadata <- prepared_data$sample_metadata
species_matrix <- prepared_data$species_matrix

# =============================================================================
# 1. Calculate Individual Antibiotic Exposures (same as main script)
# =============================================================================

cat("=== MaAsLin3 Patch Run ===\n\n")
cat("Calculating individual antibiotic exposures...\n")

top_antibiotics <- c("Pip/Tazo", "TMP/SMX", "Vancomycin", "Cefepime",
                     "Meropenem", "Ciprofloxacin", "Metronidazole",
                     "Azithromycin", "Cefazolin", "Ceftriaxone")

systemic_routes <- c("IV", "PO", "IM")
DrugTable_systemic <- DrugTable %>%
  dplyr::filter(Route %in% systemic_routes) %>%
  mutate(Date = as.Date(Date))

calc_abx_exposure <- function(mrn, sample_date, drug_name, drug_data, window = 7) {
  start_date <- sample_date - window
  end_date <- sample_date - 1

  exposed <- drug_data %>%
    dplyr::filter(MRN == mrn,
           Drug == drug_name,
           Date >= start_date,
           Date <= end_date) %>%
    nrow() > 0

  return(exposed)
}

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

abx_cols <- sapply(top_antibiotics, function(abx) {
  paste0(gsub("/", "_", gsub("-", "_", abx)), "_7d")
})
names(abx_cols) <- top_antibiotics

# =============================================================================
# 2. Prepare Species Data (same as main script)
# =============================================================================

cat("\nPreparing species data...\n")

sample_hc <- sample_metadata %>% dplyr::filter(high_confidence)
hc_samples <- sample_hc$sample_id

species_hc <- species_matrix[intersect(hc_samples, rownames(species_matrix)), ]
species_rel <- species_hc / rowSums(species_hc)

min_prevalence <- 0.10
min_mean_abundance <- 0.0001

prevalence <- colSums(species_hc > 0) / nrow(species_hc)
mean_abundance <- colMeans(species_rel)

keep_species <- (prevalence >= min_prevalence) & (mean_abundance >= min_mean_abundance)
species_filtered <- species_hc[, keep_species]

cat("Species filtering:\n")
cat("  Starting species:", ncol(species_hc), "\n")
cat("  After filters:", ncol(species_filtered), "\n")
cat("  Samples:", nrow(species_filtered), "\n")

species_counts <- t(species_filtered)

meta_aligned <- sample_hc[match(colnames(species_counts), sample_hc$sample_id), ]
rownames(meta_aligned) <- meta_aligned$sample_id

if (!"patient_group" %in% colnames(meta_aligned)) {
  group_cols <- grep("group|cohort|study", colnames(meta_aligned),
                     ignore.case = TRUE, value = TRUE)
  if (length(group_cols) > 0) {
    meta_aligned$patient_group <- meta_aligned[[group_cols[1]]]
  } else {
    meta_aligned$patient_group <- "Unknown"
  }
}

# =============================================================================
# 3. Output Directory (same location as main script)
# =============================================================================

indiv_dir <- file.path(results_dir, "individual_antibiotics_v2_with_covariates")
maaslin_dir <- file.path(indiv_dir, "maaslin3")
dir.create(maaslin_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 4. Run MaAsLin3 for Each Antibiotic
# =============================================================================

cat("\n=== Running MaAsLin3 Analysis ===\n")
cat("Adjusting for: other antibiotic exposures + patient group\n\n")

all_maaslin_results <- list()

for (abx in top_antibiotics) {
  target_col <- abx_cols[abx]

  exposure <- meta_aligned[[target_col]]
  n_exposed <- sum(exposure)
  n_unexposed <- sum(!exposure)

  cat("=== Processing:", abx, "===\n")
  cat("  Exposed:", n_exposed, ", Unexposed:", n_unexposed, "\n")

  if (n_exposed < 10 || n_unexposed < 10) {
    cat("  Skipping - insufficient samples in one group\n\n")
    next
  }

  # Define covariates
  other_abx <- setdiff(names(abx_cols), abx)
  other_abx_cols <- abx_cols[other_abx]

  all_abx_cols <- c(target_col, other_abx_cols)
  cols_to_select <- c("sample_id", all_abx_cols, "patient_group")

  covariates_df <- meta_aligned[, cols_to_select, drop = FALSE]

  for (col in all_abx_cols) {
    covariates_df[[col]] <- as.numeric(covariates_df[[col]])
  }

  rownames(covariates_df) <- covariates_df$sample_id

  # Run MaAsLin3
  cat("  Running MaAsLin3...\n")

  tryCatch({
    maaslin_input <- as.data.frame(t(species_counts))

    maaslin_meta <- covariates_df %>%
      dplyr::select(-sample_id) %>%
      mutate(patient_group = as.factor(patient_group))
    rownames(maaslin_meta) <- covariates_df$sample_id

    fixed_effects <- c(target_col, other_abx_cols, "patient_group")

    maaslin_out_dir <- file.path(maaslin_dir, gsub("/", "_", gsub("-", "_", abx)))

    maaslin_result <- maaslin3(
      input_data = maaslin_input,
      input_metadata = maaslin_meta,
      output = maaslin_out_dir,
      fixed_effects = fixed_effects,
      normalization = "TSS",
      transform = "LOG",
      min_prevalence = 0.05,
      cores = 4,
      plot_summary_plot = FALSE,
      plot_associations = FALSE
    )

    # Extract results - FIXED: use $results from the list
    if (!is.null(maaslin_result$fit_data_abundance$results)) {
      maaslin_df <- maaslin_result$fit_data_abundance$results %>%
        dplyr::filter(grepl(target_col, metadata, fixed = TRUE)) %>%
        mutate(antibiotic = abx, method = "MaAsLin3") %>%
        rename(species = feature, pval_BH = qval_individual) %>%
        arrange(pval_BH)

      all_maaslin_results[[abx]] <- maaslin_df

      file_name <- paste0("maaslin3_", gsub("/", "_", gsub("-", "_", abx)), ".csv")
      write_csv(maaslin_df, file.path(maaslin_dir, file_name))

      n_sig <- sum(maaslin_df$pval_BH < 0.1, na.rm = TRUE)
      cat("    Significant (q < 0.1):", n_sig, "\n")
    }
  }, error = function(e) {
    cat("    MaAsLin3 error:", conditionMessage(e), "\n")
  })

  cat("\n")
}

# =============================================================================
# 5. Combine Results
# =============================================================================

cat("=== Combining MaAsLin3 Results ===\n\n")

if (length(all_maaslin_results) > 0) {
  combined_maaslin <- bind_rows(all_maaslin_results)
  write_csv(combined_maaslin, file.path(maaslin_dir, "all_antibiotics_combined.csv"))
  cat("MaAsLin3: Saved", nrow(combined_maaslin), "results\n")

  # Summary by antibiotic
  summary_df <- combined_maaslin %>%
    group_by(antibiotic) %>%
    summarise(
      total = n(),
      significant = sum(pval_BH < 0.1, na.rm = TRUE),
      .groups = "drop"
    )

  cat("\nSummary by antibiotic:\n")
  print(as.data.frame(summary_df))
}

cat("\n=== MaAsLin3 Patch Run Complete ===\n")
