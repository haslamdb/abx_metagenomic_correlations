#!/usr/bin/env Rscript
# =============================================================================
# 01_rebuild_sample_metadata.R
# Rebuild sample metadata from SampleListFormatted.csv and compute antibiotic
# exposure from Drugs.csv for all samples with Bracken data
# =============================================================================

library(dplyr)
library(readr)
library(lubridate)
library(tidyr)
library(data.table)  # For faster operations

# Set paths
project_dir <- "/home/david/projects/abx_metagenomic_correlations"
setwd(project_dir)

cat("=== Rebuilding Sample Metadata ===\n\n")

# =============================================================================
# 1. Load Bracken count matrices to get sample names
# =============================================================================

load("data/bracken_count_matrices.RData")
bracken_samples <- rownames(bracken_data$genus_matrix)
cat("Bracken samples:", length(bracken_samples), "\n")

# =============================================================================
# 2. Load sample metadata from SampleListFormatted.csv
# =============================================================================

sample_list <- read_csv("data/legacy/SampleListFormatted.csv", show_col_types = FALSE) %>%
  filter(SequenceFileName %in% bracken_samples)

cat("Samples in SampleListFormatted with Bracken data:", nrow(sample_list), "\n")

# Parse dates - handle multiple formats
sample_list <- sample_list %>%
  mutate(SampleDate_parsed = mdy(SampleDate, quiet = TRUE)) %>%
  filter(!is.na(SampleDate_parsed))

cat("Samples with valid dates:", nrow(sample_list), "\n")

# =============================================================================
# 3. Load and pre-filter drug data
# =============================================================================

cat("\nLoading Drugs.csv...\n")

# Get the MRNs we need
sample_mrns <- unique(sample_list$MRN)

# Use data.table for fast reading and filtering
drugs_dt <- fread("data/legacy/Drugs.csv", select = c("Date", "MRN", "Drug", "Route"))
drugs_dt[, Date := as.Date(Date, origin = "1970-01-01")]
drugs_dt[, MRN := as.character(MRN)]

# Pre-filter to only MRNs we care about
drugs_dt <- drugs_dt[MRN %in% as.character(sample_mrns)]
setkey(drugs_dt, MRN, Date)

cat("Drug records for our samples:", nrow(drugs_dt), "\n")
cat("Unique MRNs with drug data:", length(unique(drugs_dt$MRN)), "\n")

# =============================================================================
# 4. Define antibiotics of interest
# =============================================================================

key_antibiotics <- c("Pip_Tazo", "Meropenem", "Cefepime", "Vancomycin",
                     "Metronidazole", "Ceftriaxone", "Ciprofloxacin",
                     "Clindamycin", "Azithromycin", "TMP_SMX")

# Drug name mapping (for matching in Drugs.csv)
drug_name_map <- c(
  "Pip_Tazo" = "Pip/Tazo",
  "Meropenem" = "Meropenem",
  "Cefepime" = "Cefepime",
  "Vancomycin" = "Vancomycin",
  "Metronidazole" = "Metronidazole",
  "Ceftriaxone" = "Ceftriaxone",
  "Ciprofloxacin" = "Ciprofloxacin",
  "Clindamycin" = "Clindamycin",
  "Azithromycin" = "Azithromycin",
  "TMP_SMX" = "TMP/SMX"
)

anaerobic_drugs <- c("Metronidazole", "Pip/Tazo", "Meropenem", "Clindamycin",
                     "Augmentin", "Ertapenem", "Amp/Sulb")
broad_spectrum <- c("Pip/Tazo", "Meropenem", "Cefepime", "Imipenem", "Ertapenem",
                    "Ciprofloxacin", "Levofloxacin")
gram_positive <- c("Vancomycin", "Daptomycin", "Linezolid", "Ceftaroline")
pseudomonal <- c("Pip/Tazo", "Meropenem", "Cefepime", "Ciprofloxacin",
                 "Tobramycin", "Gentamicin", "Amikacin", "Aztreonam")

# =============================================================================
# 5. Vectorized antibiotic exposure computation
# =============================================================================

cat("\nComputing antibiotic exposure...\n")

# Initialize result dataframe
sample_metadata <- sample_list %>%
  select(
    sample_id = SequenceFileName,
    MRN,
    PatientID = PMID,
    SampleDate = SampleDate_parsed,
    SampleType,
    PatientGroup,
    PatientSubGroup
  ) %>%
  mutate(MRN = as.character(MRN))

# Function to compute exposures for a window
compute_window_exposures <- function(samples_df, drugs_dt, window_days) {
  cat("  Computing", window_days, "-day window exposures...\n")

  n_samples <- nrow(samples_df)

  # Initialize result columns
  result <- data.frame(
    abx_any = logical(n_samples),
    abx_days = integer(n_samples),
    abx_anaerobic = logical(n_samples),
    abx_broad = logical(n_samples),
    abx_gram_pos = logical(n_samples),
    abx_pseudomonal = logical(n_samples),
    drugs = character(n_samples),
    stringsAsFactors = FALSE
  )

  # Add columns for individual antibiotics
  for (abx in key_antibiotics) {
    result[[abx]] <- integer(n_samples)
  }

  # Process in batches for speed
  batch_size <- 100
  n_batches <- ceiling(n_samples / batch_size)

  pb <- txtProgressBar(min = 0, max = n_batches, style = 3)

  for (b in 1:n_batches) {
    start_idx <- (b - 1) * batch_size + 1
    end_idx <- min(b * batch_size, n_samples)

    for (i in start_idx:end_idx) {
      mrn <- samples_df$MRN[i]
      sample_date <- samples_df$SampleDate[i]
      start_date <- sample_date - window_days

      # Get drugs for this patient in window
      patient_drugs <- drugs_dt[MRN == mrn & Date >= start_date & Date < sample_date]

      if (nrow(patient_drugs) == 0) {
        result$abx_any[i] <- FALSE
        result$abx_days[i] <- 0
        result$drugs[i] <- ""
      } else {
        unique_drugs <- unique(patient_drugs$Drug)
        unique_dates <- unique(patient_drugs$Date)

        result$abx_any[i] <- TRUE
        result$abx_days[i] <- length(unique_dates)
        result$abx_anaerobic[i] <- any(unique_drugs %in% anaerobic_drugs)
        result$abx_broad[i] <- any(unique_drugs %in% broad_spectrum)
        result$abx_gram_pos[i] <- any(unique_drugs %in% gram_positive)
        result$abx_pseudomonal[i] <- any(unique_drugs %in% pseudomonal)
        result$drugs[i] <- paste(unique_drugs, collapse = "; ")

        # Individual antibiotics
        for (abx in key_antibiotics) {
          drug_name <- drug_name_map[abx]
          abx_dates <- patient_drugs[Drug == drug_name]$Date
          result[[abx]][i] <- length(unique(abx_dates))
        }
      }
    }
    setTxtProgressBar(pb, b)
  }
  close(pb)

  return(result)
}

# Compute 7-day and 14-day windows
abx_7d <- compute_window_exposures(sample_metadata, drugs_dt, 7)
names(abx_7d) <- paste0(names(abx_7d), "_7d")

abx_14d <- compute_window_exposures(sample_metadata, drugs_dt, 14)
names(abx_14d) <- paste0(names(abx_14d), "_14d")

# Combine
sample_metadata <- bind_cols(sample_metadata, abx_7d, abx_14d)

# =============================================================================
# 6. Quality checks and summary
# =============================================================================

cat("\n=== Summary ===\n")
cat("Total samples with metadata:", nrow(sample_metadata), "\n")
cat("Samples with any abx (7d):", sum(sample_metadata$abx_any_7d), "\n")
cat("Samples with any abx (14d):", sum(sample_metadata$abx_any_14d), "\n")

cat("\nPatient groups:\n")
print(table(sample_metadata$PatientGroup, useNA = "ifany"))

cat("\n7-day antibiotic exposure counts:\n")
for (abx in key_antibiotics) {
  col_name <- paste0(gsub("/", "_", abx), "_7d")
  n_exposed <- sum(sample_metadata[[col_name]] > 0)
  cat(sprintf("  %s: %d samples\n", abx, n_exposed))
}

# =============================================================================
# 7. Save results
# =============================================================================

cat("\n=== Saving Results ===\n")

# Convert sample_id to rownames for compatibility
sample_metadata_matrix <- sample_metadata
rownames(sample_metadata_matrix) <- sample_metadata$sample_id

save(sample_metadata, file = "data/sample_metadata_rebuilt.RData")
cat("Saved: data/sample_metadata_rebuilt.RData\n")

write_csv(sample_metadata, "data/sample_metadata_rebuilt.csv")
cat("Saved: data/sample_metadata_rebuilt.csv\n")

# Also check overlap with Bracken matrices
cat("\n=== Bracken Matrix Overlap ===\n")
genus_overlap <- intersect(sample_metadata$sample_id, rownames(bracken_data$genus_matrix))
species_overlap <- intersect(sample_metadata$sample_id, rownames(bracken_data$species_matrix))
cat("Samples in both metadata and genus matrix:", length(genus_overlap), "\n")
cat("Samples in both metadata and species matrix:", length(species_overlap), "\n")

cat("\n=== Done ===\n")
