#!/usr/bin/env Rscript
# =============================================================================
# 11_arg_differential_abundance.R
# Differential abundance analysis for ARGs vs individual antibiotics
# Using ALDEx2, MaAsLin3, and ANCOM-BC2 with covariate adjustment
#
# Note: With ~3000 ARGs, FDR correction may be challenging. Consider:
#   - Consolidating to gene families/drug classes in a follow-up analysis
#   - Using more stringent prevalence filters
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(tibble)
library(stringr)
library(purrr)
library(ALDEx2)
library(maaslin3)
library(ANCOMBC)
library(phyloseq)

# Ensure dplyr functions aren't masked by MASS (loaded by ALDEx2)
select <- dplyr::select
filter <- dplyr::filter

# Set paths
project_dir <- "/home/david/projects/abx_metagenomic_correlations"
results_dir <- file.path(project_dir, "results/new_analysis_legacy_data")

# =============================================================================
# 1. Load ARG Data
# =============================================================================

cat("=============================================================================\n")
cat("ARG Differential Abundance Analysis\n")
cat("=============================================================================\n\n")

# Load ARG count matrices
load(file.path(project_dir, "data/arg_count_matrices.RData"))

# Load antibiotic exposure metadata
load(file.path(project_dir, "data/prepared_data.RData"))

sample_metadata <- prepared_data$sample_metadata
arg_count_matrix <- arg_data$arg_count_matrix

cat("ARG count matrix:", nrow(arg_count_matrix), "samples x", ncol(arg_count_matrix), "ARGs\n")
cat("Sample metadata:", nrow(sample_metadata), "samples\n\n")

# =============================================================================
# 2. Filter to samples with both ARG data and antibiotic metadata
# =============================================================================

common_samples <- intersect(sample_metadata$sample_id, rownames(arg_count_matrix))
cat("Samples with ARG data and antibiotic metadata:", length(common_samples), "\n")

# Filter to common samples
sample_metadata <- sample_metadata %>% dplyr::filter(sample_id %in% common_samples)
arg_count_matrix <- arg_count_matrix[common_samples, ]

cat("After filtering to common samples:\n")
cat("  ARG matrix:", nrow(arg_count_matrix), "samples x", ncol(arg_count_matrix), "ARGs\n")
cat("  Sample metadata:", nrow(sample_metadata), "samples\n\n")

# =============================================================================
# 3. Define Antibiotics to Analyze
# =============================================================================

# Top antibiotics (with Vancomycin split by route)
top_antibiotics <- c("Pip_Tazo", "TMP_SMX", "Vancomycin_IV", "Vancomycin_PO", "Cefepime",
                     "Meropenem", "Ciprofloxacin", "Metronidazole",
                     "Ceftriaxone", "Clindamycin")

# Create column name mapping
abx_cols <- paste0(top_antibiotics, "_7d")
names(abx_cols) <- top_antibiotics

# Summary of exposures
cat("=== Antibiotic Exposure Summary (7-day window) ===\n\n")
exposure_summary <- data.frame(
  Antibiotic = top_antibiotics,
  Exposed = sapply(abx_cols, function(col) {
    sum(sample_metadata[[col]] > 0, na.rm = TRUE)
  }),
  Unexposed = sapply(abx_cols, function(col) {
    sum(sample_metadata[[col]] == 0, na.rm = TRUE)
  })
)
print(exposure_summary)

# =============================================================================
# 4. Prepare ARG Data for Analysis
# =============================================================================

cat("\nPreparing ARG data for analysis...\n")

# Calculate relative abundance for filtering
arg_rel <- arg_count_matrix / rowSums(arg_count_matrix)

# Filter criteria
# Note: Using 5% prevalence to reduce multiple testing burden with ~3000 ARGs
min_prevalence <- 0.05  # Present in at least 5% of samples
min_mean_abundance <- 1e-5  # Low threshold since ARGs are sparse

prevalence <- colSums(arg_count_matrix > 0) / nrow(arg_count_matrix)
mean_abundance <- colMeans(arg_rel, na.rm = TRUE)

keep_args <- (prevalence >= min_prevalence)
arg_filtered <- arg_count_matrix[, keep_args]

cat("ARG filtering:\n")
cat("  Starting ARGs:", ncol(arg_count_matrix), "\n")
cat("  After 5% prevalence filter:", sum(prevalence >= min_prevalence), "\n")
cat("  Final ARGs:", ncol(arg_filtered), "\n")
cat("  Samples:", nrow(arg_filtered), "\n\n")

# Prepare for ALDEx2 (needs samples as columns)
arg_counts <- t(arg_filtered)

# Align metadata
meta_aligned <- sample_metadata[match(colnames(arg_counts), sample_metadata$sample_id), ]
rownames(meta_aligned) <- meta_aligned$sample_id

# Use PatientGroup column
meta_aligned$patient_group <- meta_aligned$PatientGroup
meta_aligned$patient_group[is.na(meta_aligned$patient_group)] <- "Unknown"

cat("Patient group distribution:\n")
print(table(meta_aligned$patient_group, useNA = "ifany"))

# =============================================================================
# 5. Create Output Directories
# =============================================================================

arg_dir <- file.path(results_dir, "arg_differential_abundance")
aldex_dir <- file.path(arg_dir, "aldex2")
maaslin_dir <- file.path(arg_dir, "maaslin3")
ancombc_dir <- file.path(arg_dir, "ancombc2")

dir.create(aldex_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(maaslin_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ancombc_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 6. Run Analyses for Each Antibiotic
# =============================================================================

cat("\n=== Running Differential Abundance Analyses (ARG Level) ===\n")
cat("Adjusting for: other antibiotic exposures + patient group\n\n")

all_aldex_results <- list()
all_maaslin_results <- list()
all_ancombc_results <- list()

for (abx in top_antibiotics) {
  target_col <- abx_cols[abx]

  # Get exposure status
  exposure <- meta_aligned[[target_col]] > 0
  n_exposed <- sum(exposure, na.rm = TRUE)
  n_unexposed <- sum(!exposure, na.rm = TRUE)

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

  # ---------------------------------------------------------------------------
  # 6a. ALDEx2 with GLM
  # ---------------------------------------------------------------------------

  cat("  Running ALDEx2...\n")

  tryCatch({
    set.seed(42)

    formula_str <- paste0("~ ", target_col, " + patient_group + ",
                          paste(other_abx_cols, collapse = " + "))

    mm <- model.matrix(as.formula(formula_str), data = covariates_df)

    aldex_clr <- aldex.clr(
      reads = arg_counts,
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
        arg = rownames(aldex_glm),
        estimate = aldex_glm[[target_est_col]],
        pval_BH = aldex_glm[[target_pval_col]],
        antibiotic = abx,
        method = "ALDEx2_glm"
      ) %>% arrange(pval_BH)

      all_aldex_results[[abx]] <- aldex_df

      file_name <- paste0("aldex2_", gsub("/", "_", gsub("-", "_", abx)), ".csv")
      write_csv(aldex_df, file.path(aldex_dir, file_name))

      n_sig <- sum(aldex_df$pval_BH < 0.1, na.rm = TRUE)
      cat("    Significant (BH < 0.1):", n_sig, "\n")
    }
  }, error = function(e) {
    cat("    ALDEx2 error:", conditionMessage(e), "\n")
  })

  # ---------------------------------------------------------------------------
  # 6b. MaAsLin3
  # ---------------------------------------------------------------------------

  cat("  Running MaAsLin3...\n")

  tryCatch({
    maaslin_input <- as.data.frame(t(arg_counts))

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
      min_prevalence = 0.03,  # Lower threshold for sparse ARG data
      cores = 4,
      plot_summary_plot = FALSE,
      plot_associations = FALSE
    )

    if (!is.null(maaslin_result$fit_data_abundance$results)) {
      maaslin_df <- maaslin_result$fit_data_abundance$results %>%
        dplyr::filter(grepl(target_col, metadata, fixed = TRUE)) %>%
        mutate(antibiotic = abx, method = "MaAsLin3") %>%
        rename(arg = feature, pval_BH = qval_individual) %>%
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

  # ---------------------------------------------------------------------------
  # 6c. ANCOM-BC2
  # ---------------------------------------------------------------------------

  cat("  Running ANCOM-BC2...\n")

  tryCatch({
    otu_table_ps <- otu_table(arg_counts, taxa_are_rows = TRUE)
    sample_data_ps <- sample_data(covariates_df)
    sample_names(sample_data_ps) <- covariates_df$sample_id

    ps <- phyloseq(otu_table_ps, sample_data_ps)

    fix_formula <- paste0(target_col, " + patient_group + ",
                          paste(other_abx_cols, collapse = " + "))

    ancombc_result <- ancombc2(
      data = ps,
      fix_formula = fix_formula,
      p_adj_method = "BH",
      prv_cut = 0.03,  # Lower threshold for sparse ARG data
      lib_cut = 100,   # Lower library cut for ARG data
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
        arg = res_df$taxon,
        lfc = res_df[[lfc_col]],
        pval_BH = res_df[[q_col]],
        diff_abundant = res_df[[diff_col]],
        antibiotic = abx,
        method = "ANCOMBC2"
      ) %>% arrange(pval_BH)

      all_ancombc_results[[abx]] <- ancombc_df

      file_name <- paste0("ancombc2_", gsub("/", "_", gsub("-", "_", abx)), ".csv")
      write_csv(ancombc_df, file.path(ancombc_dir, file_name))

      n_sig <- sum(ancombc_df$diff_abundant, na.rm = TRUE)
      cat("    Significant (diff_abundant):", n_sig, "\n")
    }
  }, error = function(e) {
    cat("    ANCOM-BC2 error:", conditionMessage(e), "\n")
  })

  cat("\n")
}

# =============================================================================
# 7. Combine Results
# =============================================================================

cat("=== Combining Results Across Methods ===\n\n")

if (length(all_aldex_results) > 0) {
  combined_aldex <- bind_rows(all_aldex_results)
  write_csv(combined_aldex, file.path(aldex_dir, "all_antibiotics_combined.csv"))
  cat("ALDEx2: Saved", nrow(combined_aldex), "results\n")
} else {
  combined_aldex <- data.frame()
}

if (length(all_maaslin_results) > 0) {
  combined_maaslin <- bind_rows(all_maaslin_results)
  write_csv(combined_maaslin, file.path(maaslin_dir, "all_antibiotics_combined.csv"))
  cat("MaAsLin3: Saved", nrow(combined_maaslin), "results\n")
} else {
  combined_maaslin <- data.frame()
}

if (length(all_ancombc_results) > 0) {
  combined_ancombc <- bind_rows(all_ancombc_results)
  write_csv(combined_ancombc, file.path(ancombc_dir, "all_antibiotics_combined.csv"))
  cat("ANCOM-BC2: Saved", nrow(combined_ancombc), "results\n")
} else {
  combined_ancombc <- data.frame()
}

# =============================================================================
# 8. Concordance Analysis
# =============================================================================

cat("\n=== Concordance Analysis: Associations Found by Multiple Methods ===\n\n")

get_significant <- function(df, method_name, pval_col = "pval_BH",
                            effect_col = NULL, threshold = 0.1) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  sig <- df %>%
    dplyr::filter(.data[[pval_col]] < threshold) %>%
    dplyr::select(antibiotic, arg, any_of(c(effect_col, "lfc", "estimate", "coef")))

  if (nrow(sig) > 0) {
    sig$method <- method_name
    effect_cols <- intersect(c("lfc", "estimate", "coef", effect_col), colnames(sig))
    if (length(effect_cols) > 0) {
      sig$effect <- sig[[effect_cols[1]]]
    }
  }
  return(sig)
}

sig_aldex <- if (nrow(combined_aldex) > 0) {
  get_significant(combined_aldex, "ALDEx2", "pval_BH", "estimate")
} else {
  data.frame(antibiotic = character(), arg = character(), effect = numeric(), method = character())
}

sig_maaslin <- if (nrow(combined_maaslin) > 0) {
  get_significant(combined_maaslin, "MaAsLin3", "pval_BH", "coef")
} else {
  data.frame(antibiotic = character(), arg = character(), effect = numeric(), method = character())
}

sig_ancombc <- if (nrow(combined_ancombc) > 0) {
  get_significant(combined_ancombc, "ANCOMBC2", "pval_BH", "lfc")
} else {
  data.frame(antibiotic = character(), arg = character(), effect = numeric(), method = character())
}

all_significant <- bind_rows(
  sig_aldex %>% dplyr::select(antibiotic, arg, effect, method),
  sig_maaslin %>% dplyr::select(antibiotic, arg, effect, method),
  sig_ancombc %>% dplyr::select(antibiotic, arg, effect, method)
)

if (nrow(all_significant) > 0) {
  concordance <- all_significant %>%
    group_by(antibiotic, arg) %>%
    summarise(
      n_methods = n(),
      methods = paste(sort(unique(method)), collapse = ", "),
      mean_effect = mean(effect, na.rm = TRUE),
      direction = ifelse(mean_effect > 0, "increased", "decreased"),
      .groups = "drop"
    ) %>%
    arrange(desc(n_methods), antibiotic, arg)

  robust_associations <- concordance %>%
    dplyr::filter(n_methods >= 2)

  cat("Associations found by 2+ methods:", nrow(robust_associations), "\n")
  cat("Associations found by all 3 methods:", sum(concordance$n_methods == 3), "\n\n")

  if (nrow(robust_associations) > 0) {
    cat("Top robust associations:\n")
    print(as.data.frame(head(robust_associations, 30)))
  }

  write_csv(concordance, file.path(arg_dir, "method_concordance.csv"))
  write_csv(robust_associations, file.path(arg_dir, "robust_associations.csv"))
} else {
  concordance <- data.frame()
  robust_associations <- data.frame()
  cat("No significant associations found.\n")
}

# =============================================================================
# 9. Summary Table
# =============================================================================

cat("\n=== Summary by Antibiotic (ARG Level) ===\n\n")

summary_by_abx <- data.frame(
  antibiotic = character(),
  aldex2_sig = integer(),
  maaslin3_sig = integer(),
  ancombc2_sig = integer(),
  robust_sig = integer()
)

for (abx in top_antibiotics) {
  aldex_n <- if (nrow(combined_aldex) > 0) sum(combined_aldex$antibiotic == abx & combined_aldex$pval_BH < 0.1, na.rm = TRUE) else 0
  maaslin_n <- if (nrow(combined_maaslin) > 0) sum(combined_maaslin$antibiotic == abx & combined_maaslin$pval_BH < 0.1, na.rm = TRUE) else 0
  ancombc_n <- if (nrow(combined_ancombc) > 0) sum(combined_ancombc$antibiotic == abx & combined_ancombc$diff_abundant, na.rm = TRUE) else 0
  robust_n <- if (nrow(robust_associations) > 0) sum(robust_associations$antibiotic == abx) else 0

  summary_by_abx <- rbind(summary_by_abx, data.frame(
    antibiotic = abx,
    aldex2_sig = aldex_n,
    maaslin3_sig = maaslin_n,
    ancombc2_sig = ancombc_n,
    robust_sig = robust_n
  ))
}

print(summary_by_abx)
write_csv(summary_by_abx, file.path(arg_dir, "summary_by_antibiotic.csv"))

# =============================================================================
# 10. Key ARG Classes Focus
# =============================================================================

cat("\n=== Key ARG Classes Results ===\n\n")

# Common resistance gene patterns
key_patterns <- c(
  "Tet" = "tetracycline",
  "Mef|Msr|Erm|Ere" = "macrolide",
  "CTX|TEM|SHV|OXA|KPC|NDM|VIM|IMP" = "beta-lactam",
  "Van" = "vancomycin",
  "Qnr|Aac.*6" = "fluoroquinolone/aminoglycoside",
  "Sul|Dfr" = "sulfonamide/trimethoprim",
  "Cat|Cml" = "chloramphenicol"
)

for (pattern_name in names(key_patterns)) {
  class_name <- key_patterns[pattern_name]

  class_hits <- all_significant %>%
    dplyr::filter(grepl(pattern_name, arg, ignore.case = TRUE))

  if (nrow(class_hits) > 0) {
    cat(class_name, "resistance:\n")
    class_summary <- class_hits %>%
      group_by(antibiotic) %>%
      summarise(
        n_args = n_distinct(arg),
        n_methods = n(),
        args = paste(unique(arg)[1:min(3, n_distinct(arg))], collapse = ", "),
        .groups = "drop"
      )
    print(as.data.frame(class_summary))
    cat("\n")
  }
}

# =============================================================================
# 11. Vancomycin Comparison (IV vs PO)
# =============================================================================

cat("\n=== Vancomycin Route Comparison (IV vs PO) ===\n\n")

if (nrow(all_significant) > 0) {
  vanco_iv <- all_significant %>% filter(antibiotic == "Vancomycin_IV")
  vanco_po <- all_significant %>% filter(antibiotic == "Vancomycin_PO")

  cat("Vancomycin IV associations:", n_distinct(vanco_iv$arg), "unique ARGs\n")
  cat("Vancomycin PO associations:", n_distinct(vanco_po$arg), "unique ARGs\n")

  # ARGs associated with both routes
  both_routes <- intersect(vanco_iv$arg, vanco_po$arg)
  cat("ARGs associated with both routes:", length(both_routes), "\n")

  if (length(both_routes) > 0) {
    cat("\nShared ARGs:\n")
    print(both_routes)
  }

  # Unique to each route
  iv_only <- setdiff(vanco_iv$arg, vanco_po$arg)
  po_only <- setdiff(vanco_po$arg, vanco_iv$arg)

  if (length(iv_only) > 0) {
    cat("\nARGs unique to IV:\n")
    print(iv_only[1:min(10, length(iv_only))])
  }
  if (length(po_only) > 0) {
    cat("\nARGs unique to PO:\n")
    print(po_only[1:min(10, length(po_only))])
  }
}

# =============================================================================
# 12. Save Metadata
# =============================================================================

analysis_info <- list(
  date = Sys.time(),
  level = "ARG",
  n_samples = ncol(arg_counts),
  n_args = nrow(arg_counts),
  antibiotics_analyzed = top_antibiotics,
  covariates = c("other_antibiotic_exposures", "patient_group"),
  methods = c("ALDEx2 (glm)", "MaAsLin3", "ANCOM-BC2"),
  prevalence_filter = 0.05,
  exposure_window = "7 days",
  note = "With ~3000 ARGs, FDR correction is challenging. Consider gene family consolidation."
)

saveRDS(analysis_info, file.path(arg_dir, "analysis_metadata.rds"))

cat("\n=== ARG Differential Abundance Analysis Complete ===\n")
cat("Results saved to:", arg_dir, "\n")
cat("\nKey outputs:\n")
cat("  - Individual antibiotic results in aldex2/, maaslin3/, ancombc2/\n")
cat("  - method_concordance.csv: All significant associations with method counts\n")
cat("  - robust_associations.csv: Associations found by 2+ methods\n")
cat("  - summary_by_antibiotic.csv: Overview of significant findings\n")
