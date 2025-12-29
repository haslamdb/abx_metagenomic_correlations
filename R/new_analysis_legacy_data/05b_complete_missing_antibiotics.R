#!/usr/bin/env Rscript
# =============================================================================
# 05b_complete_missing_antibiotics.R
# Complete Ceftriaxone (MaAsLin3 + ANCOM-BC2) and Clindamycin (all 3 methods)
# Using same covariate structure as second run (Vancomycin_IV/PO split)
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(stringr)
library(ALDEx2)
library(maaslin3)
library(ANCOMBC)
library(phyloseq)

select <- dplyr::select
filter <- dplyr::filter

project_dir <- "/home/david/projects/abx_metagenomic_correlations"
results_dir <- file.path(project_dir, "results/new_analysis_legacy_data")

# Load data
load(file.path(project_dir, "data/bracken_count_matrices.RData"))
load(file.path(project_dir, "data/prepared_data.RData"))

sample_metadata <- prepared_data$sample_metadata
species_matrix <- bracken_data$species_matrix

cat("=== Completing Missing Antibiotic Analyses ===\n\n")

# Filter to common samples
common_samples <- intersect(sample_metadata$sample_id, rownames(species_matrix))
sample_metadata <- sample_metadata %>% filter(sample_id %in% common_samples)
species_matrix <- species_matrix[common_samples, ]

cat("Samples:", nrow(species_matrix), "\n")
cat("Species:", ncol(species_matrix), "\n\n")

# SAME antibiotic list as second run (with Vancomycin split)
all_antibiotics <- c("Pip_Tazo", "TMP_SMX", "Vancomycin_IV", "Vancomycin_PO", "Cefepime",
                     "Meropenem", "Ciprofloxacin", "Metronidazole",
                     "Ceftriaxone", "Clindamycin")

abx_cols <- paste0(all_antibiotics, "_7d")
names(abx_cols) <- all_antibiotics

# Prepare species data (same filtering as main script)
species_rel <- species_matrix / rowSums(species_matrix)
min_prevalence <- 0.10
min_mean_abundance <- 0.0001

prevalence <- colSums(species_matrix > 0) / nrow(species_matrix)
mean_abundance <- colMeans(species_rel)
keep_species <- (prevalence >= min_prevalence) & (mean_abundance >= min_mean_abundance)
species_filtered <- species_matrix[, keep_species]

cat("Species after filtering:", ncol(species_filtered), "\n\n")

species_counts <- t(species_filtered)

meta_aligned <- sample_metadata[match(colnames(species_counts), sample_metadata$sample_id), ]
rownames(meta_aligned) <- meta_aligned$sample_id
meta_aligned$patient_group <- meta_aligned$PatientGroup
meta_aligned$patient_group[is.na(meta_aligned$patient_group)] <- "Unknown"

# Output directories
indiv_dir <- file.path(results_dir, "individual_antibiotics_v2_with_covariates")
aldex_dir <- file.path(indiv_dir, "aldex2")
maaslin_dir <- file.path(indiv_dir, "maaslin3")
ancombc_dir <- file.path(indiv_dir, "ancombc2")

# =============================================================================
# Define what to run
# =============================================================================

# Ceftriaxone: only MaAsLin3 + ANCOM-BC2 (ALDEx2 already done)
# Clindamycin: all 3 methods

run_configs <- list(
  list(abx = "Ceftriaxone", run_aldex = FALSE, run_maaslin = TRUE, run_ancombc = TRUE),
  list(abx = "Clindamycin", run_aldex = TRUE, run_maaslin = TRUE, run_ancombc = TRUE)
)

for (config in run_configs) {
  abx <- config$abx
  target_col <- abx_cols[abx]

  exposure <- meta_aligned[[target_col]] > 0
  n_exposed <- sum(exposure, na.rm = TRUE)
  n_unexposed <- sum(!exposure, na.rm = TRUE)

  cat("=== Processing:", abx, "===\n")
  cat("  Exposed:", n_exposed, ", Unexposed:", n_unexposed, "\n")

  if (n_exposed < 10 || n_unexposed < 10) {
    cat("  Skipping - insufficient samples\n\n")
    next
  }

  # Build covariates (same as second run - includes Vancomycin_IV/PO split)
  other_abx <- setdiff(names(abx_cols), abx)
  other_abx_cols <- abx_cols[other_abx]

  all_abx_cols <- c(target_col, other_abx_cols)
  cols_to_select <- c("sample_id", all_abx_cols, "patient_group")

  covariates_df <- meta_aligned[, cols_to_select, drop = FALSE]
  for (col in all_abx_cols) {
    covariates_df[[col]] <- as.numeric(covariates_df[[col]])
  }
  rownames(covariates_df) <- covariates_df$sample_id

  # ---------------------------------------------------------------------------
  # ALDEx2
  # ---------------------------------------------------------------------------
  if (config$run_aldex) {
    cat("  Running ALDEx2...\n")

    tryCatch({
      set.seed(42)

      formula_str <- paste0("~ ", target_col, " + patient_group + ",
                            paste(other_abx_cols, collapse = " + "))
      mm <- model.matrix(as.formula(formula_str), data = covariates_df)

      aldex_clr <- aldex.clr(
        reads = species_counts,
        conds = mm,
        mc.samples = 128,
        verbose = FALSE
      )

      aldex_glm <- aldex.glm(aldex_clr, mm)

      target_pval_col <- grep(paste0("^", target_col, "(TRUE)?:pval\\.padj$"),
                              colnames(aldex_glm), value = TRUE)[1]
      target_est_col <- grep(paste0("^", target_col, "(TRUE)?:Est$"),
                             colnames(aldex_glm), value = TRUE)[1]

      if (!is.na(target_pval_col) && !is.na(target_est_col)) {
        aldex_df <- data.frame(
          species = rownames(aldex_glm),
          estimate = aldex_glm[[target_est_col]],
          pval_BH = aldex_glm[[target_pval_col]],
          antibiotic = abx,
          method = "ALDEx2_glm"
        ) %>% arrange(pval_BH)

        file_name <- paste0("aldex2_", abx, ".csv")
        write_csv(aldex_df, file.path(aldex_dir, file_name))

        n_sig <- sum(aldex_df$pval_BH < 0.1, na.rm = TRUE)
        cat("    Significant (BH < 0.1):", n_sig, "\n")
      }
    }, error = function(e) {
      cat("    ALDEx2 error:", conditionMessage(e), "\n")
    })
  } else {
    cat("  Skipping ALDEx2 (already complete)\n")
  }

  # ---------------------------------------------------------------------------
  # MaAsLin3
  # ---------------------------------------------------------------------------
  if (config$run_maaslin) {
    cat("  Running MaAsLin3...\n")

    tryCatch({
      maaslin_input <- as.data.frame(t(species_counts))

      maaslin_meta <- covariates_df %>%
        select(-sample_id) %>%
        mutate(patient_group = as.factor(patient_group))
      rownames(maaslin_meta) <- covariates_df$sample_id

      fixed_effects <- c(target_col, other_abx_cols, "patient_group")

      maaslin_out_dir <- file.path(maaslin_dir, abx)

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

      if (!is.null(maaslin_result$fit_data_abundance$results)) {
        maaslin_df <- maaslin_result$fit_data_abundance$results %>%
          filter(grepl(target_col, metadata, fixed = TRUE)) %>%
          mutate(antibiotic = abx, method = "MaAsLin3") %>%
          rename(species = feature, pval_BH = qval_individual) %>%
          arrange(pval_BH)

        file_name <- paste0("maaslin3_", abx, ".csv")
        write_csv(maaslin_df, file.path(maaslin_dir, file_name))

        n_sig <- sum(maaslin_df$pval_BH < 0.1, na.rm = TRUE)
        cat("    Significant (q < 0.1):", n_sig, "\n")
      }
    }, error = function(e) {
      cat("    MaAsLin3 error:", conditionMessage(e), "\n")
    })
  }

  # ---------------------------------------------------------------------------
  # ANCOM-BC2
  # ---------------------------------------------------------------------------
  if (config$run_ancombc) {
    cat("  Running ANCOM-BC2...\n")

    tryCatch({
      otu_table_ps <- otu_table(species_counts, taxa_are_rows = TRUE)
      sample_data_ps <- sample_data(covariates_df)
      sample_names(sample_data_ps) <- covariates_df$sample_id

      ps <- phyloseq(otu_table_ps, sample_data_ps)

      fix_formula <- paste0(target_col, " + patient_group + ",
                            paste(other_abx_cols, collapse = " + "))

      ancombc_result <- ancombc2(
        data = ps,
        fix_formula = fix_formula,
        p_adj_method = "BH",
        prv_cut = 0.05,
        lib_cut = 1000,
        s0_perc = 0.05,
        group = NULL,
        struc_zero = FALSE,
        neg_lb = FALSE,
        alpha = 0.05,
        global = FALSE,
        pairwise = FALSE,
        dunnet = FALSE,
        trend = FALSE,
        iter_control = list(tol = 1e-5, max_iter = 20, verbose = FALSE),
        em_control = list(tol = 1e-5, max_iter = 100),
        mdfdr_control = NULL,
        n_cl = 4,
        verbose = FALSE
      )

      res_df <- ancombc_result$res

      lfc_col <- paste0("lfc_", target_col)
      q_col <- paste0("q_", target_col)
      diff_col <- paste0("diff_", target_col)

      if (lfc_col %in% colnames(res_df)) {
        ancombc_df <- data.frame(
          species = res_df$taxon,
          lfc = res_df[[lfc_col]],
          pval_BH = res_df[[q_col]],
          diff_abundant = res_df[[diff_col]],
          antibiotic = abx,
          method = "ANCOMBC2"
        ) %>% arrange(pval_BH)

        file_name <- paste0("ancombc2_", abx, ".csv")
        write_csv(ancombc_df, file.path(ancombc_dir, file_name))

        n_sig <- sum(ancombc_df$diff_abundant, na.rm = TRUE)
        cat("    Significant (diff_abundant):", n_sig, "\n")
      }
    }, error = function(e) {
      cat("    ANCOM-BC2 error:", conditionMessage(e), "\n")
    })
  }

  cat("\n")
}

cat("=== Missing Antibiotics Complete ===\n")
cat("Now regenerate combined files with the updated script.\n")
