#!/usr/bin/env Rscript
# =============================================================================
# 00b_build_arg_count_matrices.R
# Build ARG (Antibiotic Resistance Gene) count and coverage matrices
# from legacy alignment data
#
# Outputs:
#   - arg_count_matrix: Raw ARG read counts per sample
#   - arg_coverage_matrix: Coverage breadth (0-1) per ARG per sample
#   - sample_bacterial_reads: Total bacterial reads from Bracken (for normalization)
#   - arg_rpm_matrix: ARG reads per million bacterial reads (normalized)
# =============================================================================

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(tibble)

# Set paths
project_dir <- "/home/david/projects/abx_metagenomic_correlations"
arg_dir <- file.path(project_dir, "data/arg_legacy")
bracken_dir <- file.path(project_dir, "data/kraken2_legacy")
output_dir <- file.path(project_dir, "data")

cat("=============================================================================\n")
cat("Building ARG Count and Coverage Matrices\n")
cat("=============================================================================\n\n")

# =============================================================================
# 1. Load target samples from BSI_Abx_Samples.csv
# =============================================================================

bsi_samples <- read.csv(file.path(project_dir, "data/legacy/BSI_Abx_Samples.csv"))
target_samples <- bsi_samples$sample_id
cat("Target samples from BSI_Abx_Samples.csv:", length(target_samples), "\n\n")

# Entries to exclude (QC failures, ambiguous mappings)
exclude_patterns <- c("^__", "too_low", "ambiguous")

# =============================================================================
# 2. Build ARG Count Matrix from .CDS.csv files
# =============================================================================

cat("=== Building ARG Count Matrix ===\n")

# Find all .CDS.csv files (excluding subdirectories)
cds_files <- list.files(
  arg_dir,
  pattern = "\\.CDS\\.csv$",
  full.names = TRUE
)

# Extract sample names
cds_sample_names <- gsub("\\.CDS\\.csv$", "", basename(cds_files))

# Filter to target samples
keep_cds <- cds_sample_names %in% target_samples
target_cds_files <- cds_files[keep_cds]
cat("Found", sum(keep_cds), "ARG files matching BSI samples\n")

# Read and combine all CDS files
arg_list <- lapply(seq_along(target_cds_files), function(i) {
  f <- target_cds_files[i]
  sample_name <- gsub("\\.CDS\\.csv$", "", basename(f))

  if (i %% 200 == 0) cat("  Processing sample", i, "of", length(target_cds_files), "\n")

  tryCatch({
    # Tab-separated, no header
    df <- read_tsv(f, col_names = c("gene", "count"),
                   col_types = cols(gene = col_character(), count = col_double()),
                   show_col_types = FALSE)

    df %>%
      mutate(sample = sample_name)
  }, error = function(e) {
    cat("    Error reading:", basename(f), "-", conditionMessage(e), "\n")
    return(NULL)
  })
})

# Combine into long format
arg_long <- bind_rows(arg_list)
cat("Total rows (long format):", nrow(arg_long), "\n")

# Filter out QC/ambiguous entries
arg_long_filtered <- arg_long %>%
  filter(!str_detect(gene, paste(exclude_patterns, collapse = "|")))

n_excluded <- nrow(arg_long) - nrow(arg_long_filtered)
cat("Excluded", n_excluded, "rows (QC failures, ambiguous)\n")

# Aggregate duplicates (sum counts for same sample-gene pairs)
arg_agg <- arg_long_filtered %>%
  group_by(sample, gene) %>%
  summarise(count = sum(count), .groups = "drop")

# Pivot to wide format (samples as rows, genes as columns)
arg_wide <- arg_agg %>%
  pivot_wider(names_from = gene, values_from = count, values_fill = 0)

# Convert to matrix
arg_count_matrix <- as.matrix(arg_wide[, -1])
rownames(arg_count_matrix) <- arg_wide$sample

cat("ARG count matrix:", nrow(arg_count_matrix), "samples x", ncol(arg_count_matrix), "genes\n\n")

# =============================================================================
# 3. Build ARG Coverage Matrix from .out files
# =============================================================================

cat("=== Building ARG Coverage Matrix ===\n")

# Find all .out files
out_files <- list.files(
  arg_dir,
  pattern = "\\.out$",
  full.names = TRUE
)

out_sample_names <- gsub("\\.out$", "", basename(out_files))

# Filter to target samples
keep_out <- out_sample_names %in% target_samples
target_out_files <- out_files[keep_out]
cat("Found", sum(keep_out), "coverage files matching BSI samples\n")

# Read and parse .out files
# Format: accession, source, feature, start, end, score, strand, phase, attributes, reads, bases_covered, length, coverage
coverage_list <- lapply(seq_along(target_out_files), function(i) {
  f <- target_out_files[i]
  sample_name <- gsub("\\.out$", "", basename(f))

  if (i %% 200 == 0) cat("  Processing sample", i, "of", length(target_out_files), "\n")

  tryCatch({
    # Read the GFF-like format
    df <- read_tsv(f, col_names = FALSE, col_types = cols(.default = col_character()),
                   show_col_types = FALSE)

    # Filter to CDS rows only (contains product annotation)
    cds_rows <- df %>%
      filter(X3 == "CDS") %>%
      select(attributes = X9, coverage = X13)

    if (nrow(cds_rows) == 0) return(NULL)

    # Extract gene product name from attributes
    cds_rows %>%
      mutate(
        # Extract product= field from attributes
        gene = str_extract(attributes, "product=[^;]+"),
        gene = str_replace(gene, "product=", ""),
        coverage = as.numeric(coverage),
        sample = sample_name
      ) %>%
      filter(!is.na(gene)) %>%
      select(sample, gene, coverage)

  }, error = function(e) {
    cat("    Error reading:", basename(f), "-", conditionMessage(e), "\n")
    return(NULL)
  })
})

# Combine coverage data
coverage_long <- bind_rows(coverage_list)
cat("Total coverage rows:", nrow(coverage_long), "\n")

# Filter out QC entries
coverage_filtered <- coverage_long %>%
  filter(!str_detect(gene, paste(exclude_patterns, collapse = "|")))

# For multiple entries of same gene, take max coverage
coverage_agg <- coverage_filtered %>%
  group_by(sample, gene) %>%
  summarise(coverage = max(coverage, na.rm = TRUE), .groups = "drop")

# Pivot to wide format
coverage_wide <- coverage_agg %>%
  pivot_wider(names_from = gene, values_from = coverage, values_fill = 0)

# Convert to matrix
arg_coverage_matrix <- as.matrix(coverage_wide[, -1])
rownames(arg_coverage_matrix) <- coverage_wide$sample

cat("ARG coverage matrix:", nrow(arg_coverage_matrix), "samples x", ncol(arg_coverage_matrix), "genes\n\n")

# =============================================================================
# 4. Calculate Total Bacterial Reads from Bracken (for normalization)
# =============================================================================

cat("=== Calculating Total Bacterial Reads ===\n")

bracken_files <- list.files(
  file.path(bracken_dir, "species"),
  pattern = "_species_abundance.txt$",
  full.names = TRUE
)

bracken_sample_names <- gsub("_species_abundance.txt", "", basename(bracken_files))

# Only process samples that are in our ARG matrix
samples_to_process <- rownames(arg_count_matrix)
keep_bracken <- bracken_sample_names %in% samples_to_process
target_bracken_files <- bracken_files[keep_bracken]

cat("Calculating bacterial reads for", length(target_bracken_files), "samples\n")

bacterial_reads <- sapply(target_bracken_files, function(f) {
  sample_name <- gsub("_species_abundance.txt", "", basename(f))

  tryCatch({
    df <- read_tsv(f, col_types = cols(.default = col_character()), show_col_types = FALSE)
    # Exclude Homo sapiens from total
    df_filtered <- df %>%
      filter(name != "Homo sapiens") %>%
      mutate(reads = as.numeric(new_est_reads))
    sum(df_filtered$reads, na.rm = TRUE)
  }, error = function(e) {
    NA
  })
})

names(bacterial_reads) <- gsub("_species_abundance.txt", "", basename(target_bracken_files))

# Create data frame
sample_bacterial_reads <- data.frame(
  sample = names(bacterial_reads),
  bacterial_reads = as.numeric(bacterial_reads)
)

cat("Bacterial reads calculated for", sum(!is.na(sample_bacterial_reads$bacterial_reads)), "samples\n")
cat("  Median bacterial reads:", median(sample_bacterial_reads$bacterial_reads, na.rm = TRUE), "\n")
cat("  Range:", min(sample_bacterial_reads$bacterial_reads, na.rm = TRUE), "-",
    max(sample_bacterial_reads$bacterial_reads, na.rm = TRUE), "\n\n")

# =============================================================================
# 5. Calculate Normalized ARG Matrix (Reads Per Million bacterial reads)
# =============================================================================

cat("=== Calculating Normalized ARG Matrix (RPM) ===\n")

# Match samples
common_samples <- intersect(rownames(arg_count_matrix), sample_bacterial_reads$sample)
cat("Samples with both ARG and Bracken data:", length(common_samples), "\n")

# Subset and align
arg_counts_aligned <- arg_count_matrix[common_samples, , drop = FALSE]
bac_reads_aligned <- sample_bacterial_reads$bacterial_reads[match(common_samples, sample_bacterial_reads$sample)]

# Calculate RPM: (ARG counts / bacterial reads) * 1,000,000
arg_rpm_matrix <- sweep(arg_counts_aligned, 1, bac_reads_aligned / 1e6, "/")

cat("ARG RPM matrix:", nrow(arg_rpm_matrix), "samples x", ncol(arg_rpm_matrix), "genes\n\n")

# =============================================================================
# 6. Quality Checks
# =============================================================================

cat("=== Quality Checks ===\n\n")

# ARG count summary
arg_totals <- rowSums(arg_count_matrix)
cat("Total ARG reads per sample:\n")
cat("  Min:", min(arg_totals), "\n")
cat("  Median:", median(arg_totals), "\n")
cat("  Max:", max(arg_totals), "\n")
cat("  Samples with 0 ARG reads:", sum(arg_totals == 0), "\n\n")

# Coverage summary
coverage_means <- rowMeans(arg_coverage_matrix, na.rm = TRUE)
cat("Mean ARG coverage per sample:\n")
cat("  Min:", round(min(coverage_means), 4), "\n")
cat("  Median:", round(median(coverage_means), 4), "\n")
cat("  Max:", round(max(coverage_means), 4), "\n\n")

# Top ARGs by prevalence
arg_prevalence <- colSums(arg_count_matrix > 0) / nrow(arg_count_matrix)
top_args <- sort(arg_prevalence, decreasing = TRUE)[1:20]
cat("Top 20 ARGs by prevalence:\n")
for (i in 1:min(20, length(top_args))) {
  cat(sprintf("  %2d. %s (%.1f%%)\n", i, names(top_args)[i], top_args[i] * 100))
}

# =============================================================================
# 7. Save Matrices
# =============================================================================

cat("\n=== Saving Matrices ===\n")

# Save as RData
arg_data <- list(
  # Raw counts (for differential abundance tools)
  arg_count_matrix = arg_count_matrix,
  # Coverage breadth (for filtering low-confidence detections)
  arg_coverage_matrix = arg_coverage_matrix,
  # Normalized (for exploratory analysis, visualization)
  arg_rpm_matrix = arg_rpm_matrix,
  # Library size info
  sample_bacterial_reads = sample_bacterial_reads,
  # Metadata
  n_samples = nrow(arg_count_matrix),
  n_args = ncol(arg_count_matrix),
  n_samples_with_coverage = nrow(arg_coverage_matrix),
  date_created = Sys.time()
)

save(arg_data, file = file.path(output_dir, "arg_count_matrices.RData"))
cat("Saved: data/arg_count_matrices.RData\n")

# Save CSVs for inspection
write.csv(arg_count_matrix, file.path(output_dir, "arg_counts_raw.csv"))
write.csv(arg_rpm_matrix, file.path(output_dir, "arg_counts_rpm.csv"))
write.csv(arg_coverage_matrix, file.path(output_dir, "arg_coverage.csv"))
write.csv(sample_bacterial_reads, file.path(output_dir, "sample_bacterial_reads.csv"), row.names = FALSE)

cat("Saved: data/arg_counts_raw.csv\n")
cat("Saved: data/arg_counts_rpm.csv\n")
cat("Saved: data/arg_coverage.csv\n")
cat("Saved: data/sample_bacterial_reads.csv\n")

cat("\n=== Done ===\n")
cat("ARG count matrix:", nrow(arg_count_matrix), "samples x", ncol(arg_count_matrix), "genes\n")
cat("ARG coverage matrix:", nrow(arg_coverage_matrix), "samples x", ncol(arg_coverage_matrix), "genes\n")
cat("ARG RPM matrix:", nrow(arg_rpm_matrix), "samples x", ncol(arg_rpm_matrix), "genes\n")
