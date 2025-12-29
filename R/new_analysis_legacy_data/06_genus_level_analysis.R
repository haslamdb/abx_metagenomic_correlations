#!/usr/bin/env Rscript
# =============================================================================
# 06_genus_level_analysis.R
# Differential abundance analysis for individual antibiotics at GENUS level
# Using ALDEx2, MaAsLin3, and ANCOM-BC2 with covariate adjustment
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

# Load data
load(file.path(project_dir, "data/bracken_count_matrices.RData"))  # Raw Bracken counts (BSI samples only)
load(file.path(project_dir, "data/prepared_data.RData"))  # Antibiotic exposure data

# Use antibiotic exposure from prepared_data
sample_metadata <- prepared_data$sample_metadata

# =============================================================================
# 1. Prepare Genus Data from Raw Bracken Counts
# =============================================================================

cat("=== Genus-Level Differential Abundance Analysis ===\n\n")
cat("Using raw Bracken genus counts (BSI samples only):\n")

genus_matrix <- bracken_data$genus_matrix  # Raw Bracken counts (no Human/Salinibacter)

cat("  Genus matrix:", nrow(genus_matrix), "samples x", ncol(genus_matrix), "genera\n")
cat("  Sample metadata (prepared_data):", nrow(sample_metadata), "samples\n")
cat("  Excluded: Homo, Salinibacter\n\n")

# =============================================================================
# 2. Filter to samples with both Bracken data and antibiotic metadata
# =============================================================================

# Get samples in BOTH genus matrix and prepared_data
common_samples <- intersect(sample_metadata$sample_id, rownames(genus_matrix))
cat("Samples with Bracken data and antibiotic metadata:", length(common_samples), "\n")

# Filter to common samples
sample_metadata <- sample_metadata %>% dplyr::filter(sample_id %in% common_samples)
genus_matrix <- genus_matrix[common_samples, ]

cat("After filtering to common samples:\n")
cat("  Genus matrix:", nrow(genus_matrix), "samples x", ncol(genus_matrix), "genera\n")
cat("  Sample metadata:", nrow(sample_metadata), "samples\n\n")

# Top antibiotics to analyze (using columns available in prepared_data)
# Note: Vancomycin split by route - IV (no gut penetration) vs PO (stays in gut lumen)
top_antibiotics <- c("Pip_Tazo", "TMP_SMX", "Vancomycin_IV", "Vancomycin_PO", "Cefepime",
                     "Meropenem", "Ciprofloxacin", "Metronidazole",
                     "Ceftriaxone", "Clindamycin")
# Note: Azithromycin excluded (too few exposures), Cefazolin not in prepared_data

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

# Calculate total antibiotic count for each sample (useful covariate)
sample_metadata$total_abx_count <- rowSums(
  sample_metadata[, abx_cols, drop = FALSE] > 0,
  na.rm = TRUE
)

# =============================================================================
# 3. Prepare Genus Data for Analysis
# =============================================================================

cat("\nPreparing genus data for analysis...\n")

# Use all samples (already filtered to common samples)
genus_use <- genus_matrix

# Calculate relative abundance for filtering
genus_rel <- genus_use / rowSums(genus_use)

# Filter criteria
min_prevalence <- 0.10  # Present in at least 10% of samples
min_mean_abundance <- 0.0001  # Mean relative abundance > 0.01%

prevalence <- colSums(genus_use > 0) / nrow(genus_use)
mean_abundance <- colMeans(genus_rel, na.rm = TRUE)

keep_genera <- (prevalence >= min_prevalence) & (mean_abundance >= min_mean_abundance)
genus_filtered <- genus_use[, keep_genera]

cat("Genus filtering:\n")
cat("  Starting genera:", ncol(genus_use), "\n")
cat("  After 10% prevalence filter:", sum(prevalence >= min_prevalence), "\n")
cat("  After 0.01% abundance filter:", sum(mean_abundance >= min_mean_abundance, na.rm = TRUE), "\n")
cat("  After BOTH filters:", ncol(genus_filtered), "\n")
cat("  Samples:", nrow(genus_filtered), "\n")

# Prepare for ALDEx2 (needs samples as columns)
genus_counts <- t(genus_filtered)

# Align metadata
meta_aligned <- sample_metadata[match(colnames(genus_counts), sample_metadata$sample_id), ]
rownames(meta_aligned) <- meta_aligned$sample_id

# Use PatientGroup column
meta_aligned$patient_group <- meta_aligned$PatientGroup
meta_aligned$patient_group[is.na(meta_aligned$patient_group)] <- "Unknown"

cat("\nPatient group distribution:\n")
print(table(meta_aligned$patient_group, useNA = "ifany"))

# =============================================================================
# 4. Create Output Directories
# =============================================================================

genus_dir <- file.path(results_dir, "individual_antibiotics_genus_level")
aldex_dir <- file.path(genus_dir, "aldex2")
maaslin_dir <- file.path(genus_dir, "maaslin3")
ancombc_dir <- file.path(genus_dir, "ancombc2")

dir.create(aldex_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(maaslin_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ancombc_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 5. Run Analyses for Each Antibiotic
# =============================================================================

cat("\n=== Running Differential Abundance Analyses (Genus Level) ===\n")
cat("Adjusting for: other antibiotic exposures + patient group\n\n")

all_aldex_results <- list()
all_maaslin_results <- list()
all_ancombc_results <- list()

for (abx in top_antibiotics) {
  target_col <- abx_cols[abx]

  # Get exposure status (numeric: days of exposure, convert to binary)
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
  # 5a. ALDEx2 with GLM
  # ---------------------------------------------------------------------------

  cat("  Running ALDEx2...\n")

  tryCatch({
    set.seed(42)

    formula_str <- paste0("~ ", target_col, " + patient_group + ",
                          paste(other_abx_cols, collapse = " + "))

    mm <- model.matrix(as.formula(formula_str), data = covariates_df)

    aldex_clr <- aldex.clr(
      reads = genus_counts,
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
        genus = rownames(aldex_glm),
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
  # 5b. MaAsLin3
  # ---------------------------------------------------------------------------

  cat("  Running MaAsLin3...\n")

  tryCatch({
    maaslin_input <- as.data.frame(t(genus_counts))

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

    # Extract results - use $results from the list
    if (!is.null(maaslin_result$fit_data_abundance$results)) {
      maaslin_df <- maaslin_result$fit_data_abundance$results %>%
        dplyr::filter(grepl(target_col, metadata, fixed = TRUE)) %>%
        mutate(antibiotic = abx, method = "MaAsLin3") %>%
        rename(genus = feature, pval_BH = qval_individual) %>%
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
  # 5c. ANCOM-BC2
  # ---------------------------------------------------------------------------

  cat("  Running ANCOM-BC2...\n")

  tryCatch({
    otu_table_ps <- otu_table(genus_counts, taxa_are_rows = TRUE)
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
        genus = res_df$taxon,
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
# 6. Combine Results
# =============================================================================

cat("=== Combining Results Across Methods ===\n\n")

if (length(all_aldex_results) > 0) {
  combined_aldex <- bind_rows(all_aldex_results)
  write_csv(combined_aldex, file.path(aldex_dir, "all_antibiotics_combined.csv"))
  cat("ALDEx2: Saved", nrow(combined_aldex), "results\n")
}

if (length(all_maaslin_results) > 0) {
  combined_maaslin <- bind_rows(all_maaslin_results)
  write_csv(combined_maaslin, file.path(maaslin_dir, "all_antibiotics_combined.csv"))
  cat("MaAsLin3: Saved", nrow(combined_maaslin), "results\n")
}

if (length(all_ancombc_results) > 0) {
  combined_ancombc <- bind_rows(all_ancombc_results)
  write_csv(combined_ancombc, file.path(ancombc_dir, "all_antibiotics_combined.csv"))
  cat("ANCOM-BC2: Saved", nrow(combined_ancombc), "results\n")
}

# =============================================================================
# 7. Concordance Analysis
# =============================================================================

cat("\n=== Concordance Analysis: Associations Found by Multiple Methods ===\n\n")

get_significant <- function(df, method_name, pval_col = "pval_BH",
                            effect_col = NULL, threshold = 0.1) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  sig <- df %>%
    dplyr::filter(.data[[pval_col]] < threshold) %>%
    dplyr::select(antibiotic, genus, any_of(c(effect_col, "lfc", "estimate", "coef")))

  if (nrow(sig) > 0) {
    sig$method <- method_name
    effect_cols <- intersect(c("lfc", "estimate", "coef", effect_col), colnames(sig))
    if (length(effect_cols) > 0) {
      sig$effect <- sig[[effect_cols[1]]]
    }
  }
  return(sig)
}

sig_aldex <- if (exists("combined_aldex") && nrow(combined_aldex) > 0) {
  get_significant(combined_aldex, "ALDEx2", "pval_BH", "estimate")
} else {
  data.frame(antibiotic = character(), genus = character(), effect = numeric(), method = character())
}

sig_maaslin <- if (exists("combined_maaslin") && nrow(combined_maaslin) > 0) {
  get_significant(combined_maaslin, "MaAsLin3", "pval_BH", "coef")
} else {
  data.frame(antibiotic = character(), genus = character(), effect = numeric(), method = character())
}

sig_ancombc <- if (exists("combined_ancombc") && nrow(combined_ancombc) > 0) {
  get_significant(combined_ancombc, "ANCOMBC2", "pval_BH", "lfc")
} else {
  data.frame(antibiotic = character(), genus = character(), effect = numeric(), method = character())
}

all_significant <- bind_rows(
  sig_aldex %>% dplyr::select(antibiotic, genus, effect, method),
  sig_maaslin %>% dplyr::select(antibiotic, genus, effect, method),
  sig_ancombc %>% dplyr::select(antibiotic, genus, effect, method)
)

if (nrow(all_significant) > 0) {
  concordance <- all_significant %>%
    group_by(antibiotic, genus) %>%
    summarise(
      n_methods = n(),
      methods = paste(sort(unique(method)), collapse = ", "),
      mean_effect = mean(effect, na.rm = TRUE),
      direction = ifelse(mean_effect > 0, "increased", "decreased"),
      .groups = "drop"
    ) %>%
    arrange(desc(n_methods), antibiotic, genus)

  robust_associations <- concordance %>%
    dplyr::filter(n_methods >= 2)

  cat("Associations found by 2+ methods:", nrow(robust_associations), "\n")
  cat("Associations found by all 3 methods:", sum(concordance$n_methods == 3), "\n\n")

  if (nrow(robust_associations) > 0) {
    cat("Top robust associations:\n")
    print(as.data.frame(head(robust_associations, 20)))
  }

  write_csv(concordance, file.path(genus_dir, "method_concordance.csv"))
  write_csv(robust_associations, file.path(genus_dir, "robust_associations.csv"))
}

# =============================================================================
# 8. Summary Table
# =============================================================================

cat("\n=== Summary by Antibiotic (Genus Level) ===\n\n")

summary_by_abx <- data.frame(antibiotic = character(),
                              aldex2_sig = integer(),
                              maaslin3_sig = integer(),
                              ancombc2_sig = integer(),
                              robust_sig = integer())

for (abx in top_antibiotics) {
  aldex_n <- if (exists("combined_aldex") && !is.null(combined_aldex)) sum(combined_aldex$antibiotic == abx & combined_aldex$pval_BH < 0.1, na.rm = TRUE) else 0
  maaslin_n <- if (exists("combined_maaslin") && !is.null(combined_maaslin)) sum(combined_maaslin$antibiotic == abx & combined_maaslin$pval_BH < 0.1, na.rm = TRUE) else 0
  ancombc_n <- if (exists("combined_ancombc") && !is.null(combined_ancombc)) sum(combined_ancombc$antibiotic == abx & combined_ancombc$diff_abundant, na.rm = TRUE) else 0
  robust_n <- if (exists("robust_associations")) sum(robust_associations$antibiotic == abx) else 0

  summary_by_abx <- rbind(summary_by_abx, data.frame(
    antibiotic = abx,
    aldex2_sig = aldex_n,
    maaslin3_sig = maaslin_n,
    ancombc2_sig = ancombc_n,
    robust_sig = robust_n
  ))
}

print(summary_by_abx)
write_csv(summary_by_abx, file.path(genus_dir, "summary_by_antibiotic.csv"))

# =============================================================================
# 9. Key Genera Focus
# =============================================================================

cat("\n=== Key Genera Results ===\n\n")

key_genera <- c("Enterococcus", "Escherichia", "Klebsiella", "Staphylococcus",
                "Clostridioides", "Bacteroides", "Prevotella")

for (genus_name in key_genera) {
  genus_hits <- all_significant %>%
    dplyr::filter(grepl(genus_name, genus, ignore.case = TRUE))

  if (nrow(genus_hits) > 0) {
    cat(genus_name, ":\n")
    genus_summary <- genus_hits %>%
      group_by(antibiotic) %>%
      summarise(
        n_methods = n(),
        methods = paste(sort(unique(method)), collapse = ", "),
        mean_effect = mean(effect, na.rm = TRUE),
        .groups = "drop"
      )
    print(as.data.frame(genus_summary))
    cat("\n")
  }
}

# =============================================================================
# 10. Save Metadata
# =============================================================================

analysis_info <- list(
  date = Sys.time(),
  level = "genus",
  n_samples = ncol(genus_counts),
  n_genera = nrow(genus_counts),
  antibiotics_analyzed = top_antibiotics,
  covariates = c("other_antibiotic_exposures", "patient_group"),
  methods = c("ALDEx2 (glm)", "MaAsLin3", "ANCOM-BC2"),
  prevalence_filter = 0.10,
  abundance_filter = 0.0001,
  exposure_window = "7 days"
)

saveRDS(analysis_info, file.path(genus_dir, "analysis_metadata.rds"))

cat("\n=== Genus-Level Analysis Complete ===\n")
cat("Results saved to:", genus_dir, "\n")
