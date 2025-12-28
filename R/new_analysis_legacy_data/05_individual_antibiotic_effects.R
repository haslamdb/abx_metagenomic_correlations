#!/usr/bin/env Rscript
# =============================================================================
# 05_individual_antibiotic_effects.R
# Differential abundance analysis for individual antibiotics at species level
# Using ALDEx2, MaAsLin3, and ANCOM-BC2 with covariate adjustment
# =============================================================================

# Load tidyverse components individually (tidyverse meta-package not installed)
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
load(file.path(project_dir, "data/prepared_data.RData"))
load(file.path(project_dir, "data/legacy/AbxEffectData20191014"))

sample_metadata <- prepared_data$sample_metadata
species_matrix <- prepared_data$species_matrix

# =============================================================================
# 1. Calculate Individual Antibiotic Exposures
# =============================================================================

cat("Calculating individual antibiotic exposures...\n")

# Top 10 antibiotics
top_antibiotics <- c("Pip/Tazo", "TMP/SMX", "Vancomycin", "Cefepime",
                     "Meropenem", "Ciprofloxacin", "Metronidazole",
                     "Azithromycin", "Cefazolin", "Ceftriaxone")

# Filter drug table to systemic routes
systemic_routes <- c("IV", "PO", "IM")
DrugTable_systemic <- DrugTable %>%
  dplyr::filter(Route %in% systemic_routes) %>%
  mutate(Date = as.Date(Date))

# Function to calculate exposure for a specific antibiotic
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

# Create clean column names for modeling
abx_cols <- sapply(top_antibiotics, function(abx) {
  paste0(gsub("/", "_", gsub("-", "_", abx)), "_7d")
})
names(abx_cols) <- top_antibiotics

# Summary of exposures
cat("\n=== Antibiotic Exposure Summary (7-day window) ===\n\n")
exposure_summary <- data.frame(
  Antibiotic = top_antibiotics,
  Exposed = sapply(top_antibiotics, function(abx) {
    sum(sample_metadata[[abx_cols[abx]]])
  }),
  Unexposed = sapply(top_antibiotics, function(abx) {
    sum(!sample_metadata[[abx_cols[abx]]])
  })
)
print(exposure_summary)

# Calculate total antibiotic count for each sample (useful covariate)
sample_metadata$total_abx_count <- rowSums(
  sample_metadata[, abx_cols, drop = FALSE],
  na.rm = TRUE
)

# =============================================================================
# 2. Prepare Species Data
# =============================================================================

cat("\nPreparing species data...\n")

# Use high-confidence samples
sample_hc <- sample_metadata %>% dplyr::filter(high_confidence)
hc_samples <- sample_hc$sample_id

species_hc <- species_matrix[intersect(hc_samples, rownames(species_matrix)), ]

# Calculate relative abundance for filtering
species_rel <- species_hc / rowSums(species_hc)

# Filter criteria (more stringent to reduce multiple testing burden)
min_prevalence <- 0.10  # Present in at least 10% of samples
min_mean_abundance <- 0.0001  # Mean relative abundance > 0.01%

prevalence <- colSums(species_hc > 0) / nrow(species_hc)
mean_abundance <- colMeans(species_rel)

# Apply both filters
keep_species <- (prevalence >= min_prevalence) & (mean_abundance >= min_mean_abundance)
species_filtered <- species_hc[, keep_species]

cat("Species filtering:\n")
cat("  Starting species:", ncol(species_hc), "\n")
cat("  After 10% prevalence filter:", sum(prevalence >= min_prevalence), "\n")
cat("  After 0.01% abundance filter:", sum(mean_abundance >= min_mean_abundance), "\n")
cat("  After BOTH filters:", ncol(species_filtered), "\n")
cat("Samples:", nrow(species_filtered), "\n")

# Prepare for ALDEx2 (needs samples as columns)
species_counts <- t(species_filtered)

# Align metadata with species data
meta_aligned <- sample_hc[match(colnames(species_counts), sample_hc$sample_id), ]
rownames(meta_aligned) <- meta_aligned$sample_id

# Check patient group availability
if (!"patient_group" %in% colnames(meta_aligned)) {
  # Try to find alternative column names
  group_cols <- grep("group|cohort|study", colnames(meta_aligned),
                     ignore.case = TRUE, value = TRUE)
  if (length(group_cols) > 0) {
    meta_aligned$patient_group <- meta_aligned[[group_cols[1]]]
    cat("Using", group_cols[1], "as patient_group\n")
  } else {
    cat("Warning: No patient group column found. Creating placeholder.\n")
    meta_aligned$patient_group <- "Unknown"
  }
}

cat("\nPatient group distribution:\n")
print(table(meta_aligned$patient_group))

# =============================================================================
# 3. Create Output Directories
# =============================================================================

# Use separate output directory to avoid overwriting legacy results
# Change to "individual_antibiotics" once legacy analysis is complete
indiv_dir <- file.path(results_dir, "individual_antibiotics_v2_with_covariates")
aldex_dir <- file.path(indiv_dir, "aldex2")
maaslin_dir <- file.path(indiv_dir, "maaslin3")
ancombc_dir <- file.path(indiv_dir, "ancombc2")

dir.create(aldex_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(maaslin_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ancombc_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 4. Run Analyses for Each Antibiotic (with covariate adjustment)
# =============================================================================

cat("\n=== Running Differential Abundance Analyses ===\n")
cat("Adjusting for: other antibiotic exposures + patient group\n\n")

all_aldex_results <- list()
all_maaslin_results <- list()
all_ancombc_results <- list()

for (abx in top_antibiotics) {
  target_col <- abx_cols[abx]

  # Get exposure status
  exposure <- meta_aligned[[target_col]]
  n_exposed <- sum(exposure)
  n_unexposed <- sum(!exposure)

  cat("=== Processing:", abx, "===\n")
  cat("  Exposed:", n_exposed, ", Unexposed:", n_unexposed, "\n")

  # Skip if too few in either group
  if (n_exposed < 10 || n_unexposed < 10) {
    cat("  Skipping - insufficient samples in one group\n\n")
    next
  }

  # ---------------------------------------------------------------------------
  # Define covariates: other antibiotics + patient group
  # ---------------------------------------------------------------------------

  other_abx <- setdiff(names(abx_cols), abx)
  other_abx_cols <- abx_cols[other_abx]

  # Build covariate dataframe - select columns explicitly
  all_abx_cols <- c(target_col, other_abx_cols)
  cols_to_select <- c("sample_id", all_abx_cols, "patient_group")

  covariates_df <- meta_aligned[, cols_to_select, drop = FALSE]

  # Convert antibiotic columns to numeric
  for (col in all_abx_cols) {
    covariates_df[[col]] <- as.numeric(covariates_df[[col]])
  }

  rownames(covariates_df) <- covariates_df$sample_id

  # ---------------------------------------------------------------------------
  # 4a. ALDEx2 with GLM (handles covariates)
  # ---------------------------------------------------------------------------

  cat("  Running ALDEx2...\n")

  tryCatch({
    set.seed(42)

    # Create model matrix with target antibiotic + other antibiotics + patient group
    # Note: ALDEx2 glm requires a model matrix
    formula_str <- paste0("~ ", target_col, " + patient_group + ",
                          paste(other_abx_cols, collapse = " + "))

    mm <- model.matrix(as.formula(formula_str), data = covariates_df)

    # Run ALDEx2 with CLR transformation
    aldex_clr <- aldex.clr(
      reads = species_counts,
      conds = mm,
      mc.samples = 128,
      verbose = FALSE
    )

    # Run GLM test
    aldex_glm <- aldex.glm(aldex_clr, mm)

    # Extract results for target antibiotic
    # Column names follow pattern: model.{covariate}Pr(>|t|).BH
    target_pval_col <- grep(paste0(target_col, ".*Pr.*BH"),
                            colnames(aldex_glm), value = TRUE)[1]
    target_est_col <- grep(paste0(target_col, ".*Estimate"),
                           colnames(aldex_glm), value = TRUE)[1]

    if (!is.na(target_pval_col) && !is.na(target_est_col)) {
      aldex_df <- data.frame(
        species = rownames(aldex_glm),
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
  # 4b. MaAsLin3 (designed for multivariable adjustment)
  # ---------------------------------------------------------------------------

  cat("  Running MaAsLin3...\n")

  tryCatch({
    # Prepare input for MaAsLin3
    # MaAsLin3 expects samples as rows, features as columns
    maaslin_input <- as.data.frame(t(species_counts))

    # Prepare metadata
    maaslin_meta <- covariates_df %>%
      dplyr::select(-sample_id) %>%
      mutate(patient_group = as.factor(patient_group))

    # Define fixed effects: target antibiotic + other antibiotics + patient group
    fixed_effects <- c(target_col, other_abx_cols, "patient_group")

    # Run MaAsLin3
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

    # Extract results for target antibiotic
    if (!is.null(maaslin_result$fit_data_abundance)) {
      maaslin_df <- maaslin_result$fit_data_abundance %>%
        dplyr::filter(grepl(target_col, metadata, fixed = TRUE)) %>%
        mutate(antibiotic = abx, method = "MaAsLin3") %>%
        rename(species = feature, pval_BH = qval) %>%
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
  # 4c. ANCOM-BC2 (with bias correction and sensitivity analysis)
  # ---------------------------------------------------------------------------

  cat("  Running ANCOM-BC2...\n")

  tryCatch({
    # Create phyloseq object
    otu_table_ps <- otu_table(species_counts, taxa_are_rows = TRUE)
    sample_data_ps <- sample_data(covariates_df)
    sample_names(sample_data_ps) <- covariates_df$sample_id

    ps <- phyloseq(otu_table_ps, sample_data_ps)

    # Build formula: target + other antibiotics + patient group
    fix_formula <- paste0(target_col, " + patient_group + ",
                          paste(other_abx_cols, collapse = " + "))

    # Run ANCOM-BC2
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

    # Extract results for target antibiotic
    res_df <- ancombc_result$res

    # Find columns for target antibiotic
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
# 5. Combine and Compare Results Across Methods
# =============================================================================

cat("=== Combining Results Across Methods ===\n\n")

# Combine ALDEx2 results
if (length(all_aldex_results) > 0) {
  combined_aldex <- bind_rows(all_aldex_results)
  write_csv(combined_aldex, file.path(aldex_dir, "all_antibiotics_combined.csv"))
  cat("ALDEx2: Saved", nrow(combined_aldex), "results\n")
}

# Combine MaAsLin3 results
if (length(all_maaslin_results) > 0) {
  combined_maaslin <- bind_rows(all_maaslin_results)
  write_csv(combined_maaslin, file.path(maaslin_dir, "all_antibiotics_combined.csv"))
  cat("MaAsLin3: Saved", nrow(combined_maaslin), "results\n")
}

# Combine ANCOM-BC2 results
if (length(all_ancombc_results) > 0) {
  combined_ancombc <- bind_rows(all_ancombc_results)
  write_csv(combined_ancombc, file.path(ancombc_dir, "all_antibiotics_combined.csv"))
  cat("ANCOM-BC2: Saved", nrow(combined_ancombc), "results\n")
}

# =============================================================================
# 6. Concordance Analysis - Find Robust Associations
# =============================================================================

cat("\n=== Concordance Analysis: Associations Found by Multiple Methods ===\n\n")

# Function to extract significant associations
get_significant <- function(df, method_name, pval_col = "pval_BH",
                            effect_col = NULL, threshold = 0.1) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  sig <- df %>%
    dplyr::filter(.data[[pval_col]] < threshold) %>%
    dplyr::select(antibiotic, species, any_of(c(effect_col, "lfc", "estimate", "coef")))

  if (nrow(sig) > 0) {
    sig$method <- method_name
    # Standardize effect column name
    effect_cols <- intersect(c("lfc", "estimate", "coef", effect_col), colnames(sig))
    if (length(effect_cols) > 0) {
      sig$effect <- sig[[effect_cols[1]]]
    }
  }
  return(sig)
}

# Get significant from each method
sig_aldex <- get_significant(combined_aldex, "ALDEx2", "pval_BH", "estimate")
sig_maaslin <- get_significant(combined_maaslin, "MaAsLin3", "pval_BH", "coef")
sig_ancombc <- get_significant(combined_ancombc, "ANCOMBC2", "pval_BH", "lfc")

# Combine
all_significant <- bind_rows(
  sig_aldex %>% dplyr::select(antibiotic, species, effect, method),
  sig_maaslin %>% dplyr::select(antibiotic, species, effect, method),
  sig_ancombc %>% dplyr::select(antibiotic, species, effect, method)
)

if (nrow(all_significant) > 0) {
  # Count methods per association
  concordance <- all_significant %>%
    group_by(antibiotic, species) %>%
    summarise(
      n_methods = n(),
      methods = paste(sort(unique(method)), collapse = ", "),
      mean_effect = mean(effect, na.rm = TRUE),
      direction = ifelse(mean_effect > 0, "increased", "decreased"),
      .groups = "drop"
    ) %>%
    arrange(desc(n_methods), antibiotic, species)

  # Filter to associations found by 2+ methods
  robust_associations <- concordance %>%
    dplyr::filter(n_methods >= 2)

  cat("Associations found by 2+ methods:", nrow(robust_associations), "\n")
  cat("Associations found by all 3 methods:", sum(concordance$n_methods == 3), "\n\n")

  if (nrow(robust_associations) > 0) {
    cat("Top robust associations:\n")
    print(as.data.frame(head(robust_associations, 20)))
  }

  write_csv(concordance, file.path(indiv_dir, "method_concordance.csv"))
  write_csv(robust_associations, file.path(indiv_dir, "robust_associations.csv"))
}

# =============================================================================
# 7. Summary Table by Antibiotic
# =============================================================================

cat("\n=== Summary by Antibiotic ===\n\n")

summary_by_abx <- data.frame(antibiotic = character(),
                              aldex2_sig = integer(),
                              maaslin3_sig = integer(),
                              ancombc2_sig = integer(),
                              robust_sig = integer())

for (abx in top_antibiotics) {
  aldex_n <- if (!is.null(combined_aldex)) sum(combined_aldex$antibiotic == abx & combined_aldex$pval_BH < 0.1, na.rm = TRUE) else 0
  maaslin_n <- if (!is.null(combined_maaslin)) sum(combined_maaslin$antibiotic == abx & combined_maaslin$pval_BH < 0.1, na.rm = TRUE) else 0
  ancombc_n <- if (!is.null(combined_ancombc)) sum(combined_ancombc$antibiotic == abx & combined_ancombc$diff_abundant, na.rm = TRUE) else 0
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
write_csv(summary_by_abx, file.path(indiv_dir, "summary_by_antibiotic.csv"))

# =============================================================================
# 8. Enterococcus Focus Analysis
# =============================================================================

cat("\n=== Enterococcus Results by Antibiotic (All Methods) ===\n\n")

enterococcus_all <- all_significant %>%
  dplyr::filter(grepl("Enterococcus", species, ignore.case = TRUE))

if (nrow(enterococcus_all) > 0) {
  enterococcus_summary <- enterococcus_all %>%
    group_by(antibiotic, species) %>%
    summarise(
      n_methods = n(),
      methods = paste(sort(unique(method)), collapse = ", "),
      mean_effect = mean(effect, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_methods), antibiotic)

  print(as.data.frame(enterococcus_summary))
  write_csv(enterococcus_summary, file.path(indiv_dir, "enterococcus_by_antibiotic.csv"))
} else {
  cat("No significant Enterococcus associations found.\n")
}

# =============================================================================
# 9. Save Analysis Metadata
# =============================================================================

analysis_info <- list(
  date = Sys.time(),
  n_samples = nrow(meta_aligned),
  n_species = nrow(species_counts),
  antibiotics_analyzed = top_antibiotics,
  covariates = c("other_antibiotic_exposures", "patient_group"),
  methods = c("ALDEx2 (glm)", "MaAsLin3", "ANCOM-BC2"),
  prevalence_filter = 0.05,
  exposure_window = "7 days",
  patient_groups = unique(meta_aligned$patient_group)
)

saveRDS(analysis_info, file.path(indiv_dir, "analysis_metadata.rds"))

cat("\n=== Analysis Complete ===\n")
cat("Results saved to:", indiv_dir, "\n")
cat("\nKey outputs:\n")
cat("  - Individual antibiotic results in aldex2/, maaslin3/, ancombc2/\n")
cat("  - method_concordance.csv: All significant associations with method counts\n")
cat("  - robust_associations.csv: Associations found by 2+ methods\n")
cat("  - summary_by_antibiotic.csv: Overview of significant findings\n")
