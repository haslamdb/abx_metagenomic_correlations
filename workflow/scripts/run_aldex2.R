#!/usr/bin/env Rscript
#' ALDEx2 Differential Abundance Analysis
#'
#' Compositional-aware differential abundance testing comparing
#' samples with and without recent antibiotic exposure.

library(ALDEx2)
library(tidyverse)

# Load data from Snakemake
species_file <- snakemake@input[["species"]]
metadata_file <- snakemake@input[["metadata"]]
drug_file <- snakemake@input[["drug_exposure"]]
output_results <- snakemake@output[["results"]]
output_significant <- snakemake@output[["significant"]]

# Read data
species <- read.delim(species_file, row.names = 1, check.names = FALSE)
metadata <- read.delim(metadata_file, stringsAsFactors = FALSE)
drug_exposure <- read.delim(drug_file, stringsAsFactors = FALSE)

# Parse dates
metadata$sample_date <- as.Date(metadata$sample_date)
drug_exposure$date <- as.Date(drug_exposure$date)

# Calculate antibiotic exposure in 7 days before each sample
calculate_exposure <- function(sample_id, sample_date, patient_id, drug_data, window_days = 7) {
  start_date <- sample_date - window_days
  patient_drugs <- drug_data %>%
    filter(patient_id == !!patient_id,
           date >= start_date,
           date <= sample_date)
  nrow(patient_drugs) > 0
}

metadata$abx_7d <- mapply(
  calculate_exposure,
  metadata$sample_id,
  metadata$sample_date,
  metadata$patient_id,
  MoreArgs = list(drug_data = drug_exposure, window_days = 7)
)

message(sprintf("Samples with Abx in prior 7 days: %d/%d",
                sum(metadata$abx_7d), nrow(metadata)))

# Ensure samples are in same order
species <- species[metadata$sample_id, ]

# Transpose for ALDEx2 (taxa as rows, samples as columns)
species_t <- t(species)

# Convert counts to integers (ALDEx2 requirement)
# If using relative abundance, convert back to pseudo-counts
if (max(species_t) <= 1) {
  species_t <- round(species_t * 1e6)  # Convert to counts per million
}
species_t <- apply(species_t, c(1, 2), as.integer)

# Remove taxa with very low counts
min_count <- 10
taxa_keep <- rowSums(species_t > 0) >= (ncol(species_t) * 0.1)  # Present in 10% of samples
species_t <- species_t[taxa_keep, ]

message(sprintf("Analyzing %d taxa across %d samples", nrow(species_t), ncol(species_t)))

# Define conditions
conditions <- ifelse(metadata$abx_7d, "exposed", "unexposed")

# Run ALDEx2
message("Running ALDEx2 (this may take a while)...")

aldex_results <- aldex(
  reads = species_t,
  conditions = conditions,
  mc.samples = 128,
  test = "t",
  effect = TRUE,
  include.sample.summary = FALSE,
  verbose = TRUE
)

# Add taxon names
aldex_results$taxon <- rownames(aldex_results)

# Sort by effect size
aldex_results <- aldex_results %>%
  arrange(desc(abs(effect)))

# Filter significant results
significant <- aldex_results %>%
  filter(wi.eBH < 0.1 | we.eBH < 0.1)

message(sprintf("Found %d significant taxa (BH-adjusted p < 0.1)", nrow(significant)))

# Save results
write.table(aldex_results, output_results, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(significant, output_significant, sep = "\t", row.names = FALSE, quote = FALSE)

message("ALDEx2 analysis complete!")
