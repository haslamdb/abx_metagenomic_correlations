#!/usr/bin/env Rscript
# =============================================================================
# 00_build_count_matrices_from_bracken.R
# Build genus and species count matrices from raw Bracken output files
# Filters out Human (Homo sapiens) and Salinibacter (common contaminant)
# =============================================================================

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(tibble)

# Set paths
project_dir <- "/home/david/projects/abx_metagenomic_correlations"
bracken_dir <- file.path(project_dir, "data/kraken2_legacy")
output_dir <- file.path(project_dir, "data")

# =============================================================================
# 1. Get sample IDs from BSI_Abx_Samples.csv - ONLY these samples
# =============================================================================

cat("=== Building Count Matrices from Bracken Files ===\n\n")

# Load BSI/Abx sample list - restrict to ONLY these samples
bsi_samples <- read.csv(file.path(project_dir, "data/legacy/BSI_Abx_Samples.csv"))
target_samples <- bsi_samples$sample_id
cat("Target samples from BSI_Abx_Samples.csv:", length(target_samples), "\n")

# Taxa to exclude (contaminants)
exclude_genera <- c("Homo", "Salinibacter")
exclude_species <- c("Homo sapiens", "Salinibacter ruber")

# =============================================================================
# 2. Build Genus Count Matrix
# =============================================================================

cat("Building genus count matrix...\n")

genus_files <- list.files(
  file.path(bracken_dir, "genus"),
  pattern = "_genus_abundance.txt$",
  full.names = TRUE
)

# Extract sample names from filenames
genus_sample_names <- gsub("_genus_abundance.txt", "", basename(genus_files))

# Filter to ONLY samples in BSI_Abx_Samples.csv
keep_genus <- genus_sample_names %in% target_samples

bsi_genus_files <- genus_files[keep_genus]
cat("  Found", sum(keep_genus), "files matching BSI_Abx_Samples\n")

# Read and combine all genus files
genus_list <- lapply(bsi_genus_files, function(f) {
  sample_name <- gsub("_genus_abundance.txt", "", basename(f))

  tryCatch({
    df <- read_tsv(f, col_types = cols(
      name = col_character(),
      taxonomy_id = col_double(),
      taxonomy_lvl = col_character(),
      kraken_assigned_reads = col_double(),
      added_reads = col_double(),
      new_est_reads = col_double(),
      fraction_total_reads = col_double()
    ), show_col_types = FALSE)

    # Use new_est_reads as the count
    df %>%
      select(name, count = new_est_reads) %>%
      mutate(sample = sample_name)
  }, error = function(e) {
    cat("    Error reading:", basename(f), "-", conditionMessage(e), "\n")
    return(NULL)
  })
})

# Combine into long format
genus_long <- bind_rows(genus_list)
cat("  Total rows (long format):", nrow(genus_long), "\n")

# Filter out contaminants
genus_long_filtered <- genus_long %>%
  filter(!name %in% exclude_genera)

n_excluded <- nrow(genus_long) - nrow(genus_long_filtered)
cat("  Excluded", n_excluded, "rows (Homo, Salinibacter)\n")

# Aggregate any duplicates (sum counts for same sample-genus pairs)
genus_agg <- genus_long_filtered %>%
  group_by(sample, name) %>%
  summarise(count = sum(count), .groups = "drop")

# Pivot to wide format (samples as rows, genera as columns)
genus_wide <- genus_agg %>%
  pivot_wider(names_from = name, values_from = count, values_fill = 0)

# Convert to matrix with rownames
genus_matrix <- as.matrix(genus_wide[, -1])
rownames(genus_matrix) <- genus_wide$sample

cat("  Final genus matrix:", nrow(genus_matrix), "samples x", ncol(genus_matrix), "genera\n")

# Verify no contaminants
if (any(colnames(genus_matrix) %in% exclude_genera)) {
  warning("Contaminant genera still present!")
} else {
  cat("  Verified: No Homo or Salinibacter in matrix\n")
}

# =============================================================================
# 3. Build Species Count Matrix
# =============================================================================

cat("\nBuilding species count matrix...\n")

species_files <- list.files(
  file.path(bracken_dir, "species"),
  pattern = "_species_abundance.txt$",
  full.names = TRUE
)

# Extract sample names from filenames
species_sample_names <- gsub("_species_abundance.txt", "", basename(species_files))

# Filter to ONLY samples in BSI_Abx_Samples.csv
keep_species <- species_sample_names %in% target_samples

bsi_species_files <- species_files[keep_species]
cat("  Found", sum(keep_species), "files matching BSI_Abx_Samples\n")

# Read and combine all species files
species_list <- lapply(bsi_species_files, function(f) {
  sample_name <- gsub("_species_abundance.txt", "", basename(f))

  tryCatch({
    df <- read_tsv(f, col_types = cols(
      name = col_character(),
      taxonomy_id = col_double(),
      taxonomy_lvl = col_character(),
      kraken_assigned_reads = col_double(),
      added_reads = col_double(),
      new_est_reads = col_double(),
      fraction_total_reads = col_double()
    ), show_col_types = FALSE)

    df %>%
      select(name, count = new_est_reads) %>%
      mutate(sample = sample_name)
  }, error = function(e) {
    cat("    Error reading:", basename(f), "-", conditionMessage(e), "\n")
    return(NULL)
  })
})

# Combine into long format
species_long <- bind_rows(species_list)
cat("  Total rows (long format):", nrow(species_long), "\n")

# Filter out contaminants
species_long_filtered <- species_long %>%
  filter(!name %in% exclude_species)

n_excluded_sp <- nrow(species_long) - nrow(species_long_filtered)
cat("  Excluded", n_excluded_sp, "rows (Homo sapiens, Salinibacter ruber)\n")

# Aggregate any duplicates
species_agg <- species_long_filtered %>%
  group_by(sample, name) %>%
  summarise(count = sum(count), .groups = "drop")

# Pivot to wide format
species_wide <- species_agg %>%
  pivot_wider(names_from = name, values_from = count, values_fill = 0)

# Convert to matrix with rownames
species_matrix <- as.matrix(species_wide[, -1])
rownames(species_matrix) <- species_wide$sample

cat("  Final species matrix:", nrow(species_matrix), "samples x", ncol(species_matrix), "species\n")

# Verify no contaminants
if (any(colnames(species_matrix) %in% exclude_species)) {
  warning("Contaminant species still present!")
} else {
  cat("  Verified: No Homo sapiens or Salinibacter ruber in matrix\n")
}

# =============================================================================
# 4. Quality Checks
# =============================================================================

cat("\n=== Quality Checks ===\n")

# Check row sums (library sizes)
genus_lib_sizes <- rowSums(genus_matrix)
species_lib_sizes <- rowSums(species_matrix)

cat("\nGenus library sizes:\n")
cat("  Min:", min(genus_lib_sizes), "\n")
cat("  Max:", max(genus_lib_sizes), "\n")
cat("  Median:", median(genus_lib_sizes), "\n")
cat("  Samples < 10,000 reads:", sum(genus_lib_sizes < 10000), "\n")

cat("\nSpecies library sizes:\n")
cat("  Min:", min(species_lib_sizes), "\n")
cat("  Max:", max(species_lib_sizes), "\n")
cat("  Median:", median(species_lib_sizes), "\n")
cat("  Samples < 10,000 reads:", sum(species_lib_sizes < 10000), "\n")

# Check for samples with very low counts (potential QC failures)
low_count_samples <- rownames(genus_matrix)[genus_lib_sizes < 10000]
if (length(low_count_samples) > 0) {
  cat("\nWarning:", length(low_count_samples), "samples have < 10,000 reads\n")
  cat("  First 10:", paste(head(low_count_samples, 10), collapse = ", "), "\n")
}

# =============================================================================
# 5. Save Matrices
# =============================================================================

cat("\n=== Saving Matrices ===\n")

# Save as RData
bracken_data <- list(
  genus_matrix = genus_matrix,
  species_matrix = species_matrix,
  excluded_genera = exclude_genera,
  excluded_species = exclude_species,
  n_samples = nrow(genus_matrix),
  n_genera = ncol(genus_matrix),
  n_species = ncol(species_matrix),
  date_created = Sys.time()
)

save(bracken_data, file = file.path(output_dir, "bracken_count_matrices.RData"))
cat("Saved: data/bracken_count_matrices.RData\n")

# Also save as CSV for inspection
write.csv(genus_matrix, file.path(output_dir, "genus_counts_raw.csv"))
write.csv(species_matrix, file.path(output_dir, "species_counts_raw.csv"))
cat("Saved: data/genus_counts_raw.csv\n")
cat("Saved: data/species_counts_raw.csv\n")

cat("\n=== Done ===\n")
cat("Genus matrix:", nrow(genus_matrix), "x", ncol(genus_matrix), "\n")
cat("Species matrix:", nrow(species_matrix), "x", ncol(species_matrix), "\n")
