#!/usr/bin/env Rscript
# =============================================================================
# 01_load_and_prepare_data.R
# Load legacy RData and prepare enhanced data structures for analysis
# =============================================================================

library(tidyverse)
library(vegan)

# Set paths
project_dir <- "/home/david/projects/abx_metagenomic_correlations"
data_dir <- file.path(project_dir, "data/legacy")
results_dir <- file.path(project_dir, "results/new_analysis_legacy_data")

# =============================================================================
# 1. Load Legacy Data
# =============================================================================

cat("Loading legacy RData file...\n")
load(file.path(data_dir, "AbxEffectData20191014"))

# =============================================================================
# 2. Antibiotic Spectrum Classification
# =============================================================================

# Define antibiotic categories based on spectrum
abx_categories <- list(
  anti_anaerobic = c(
    "Metronidazole", "Pip/Tazo", "Meropenem", "Ertapenem", "Doripenem",
    "Clindamycin", "Augmentin", "Amp/Sulb", "Ticarcillin/clavulanate",
    "Cefoxitin", "Tigecycline"
  ),
  anti_gram_positive = c(
    "Vancomycin", "Daptomycin", "Linezolid", "Ceftaroline"
  ),
  broad_spectrum = c(
    "Meropenem", "Pip/Tazo", "Cefepime", "Doripenem", "Ertapenem",
    "Ceftazidime/avibactam", "Ceftolozane/tazobactam"
  ),
  anti_pseudomonal = c(
    "Cefepime", "Pip/Tazo", "Meropenem", "Ceftazidime", "Tobramycin",
    "Ciprofloxacin", "Aztreonam", "Ceftazidime/avibactam",
    "Ceftolozane/tazobactam", "Colistimethate", "Polymyxin B Sulfate"
  ),
  fluoroquinolones = c(
    "Ciprofloxacin", "Levofloxacin", "Moxifloxacin"
  ),
  aminoglycosides = c(
    "Gentamicin", "Tobramycin", "Amikacin"
  ),
  antifungals = c(
    "Fluconazole", "Voriconazole", "Posaconazole", "Caspofungin",
    "Micafungin", "Ambisome", "Amphotericin B", "Nystatin", "Clotrimazole"
  ),
  antivirals = c(
    "Acyclovir", "Valacyclovir", "Ganciclovir", "Valganciclovir",
    "Foscarnet", "Cidofovir", "Oseltamivir"
  )
)

# Key individual antibiotics to track
# Note: Vancomycin split by route - IV (no gut penetration) vs PO (stays in gut lumen)
key_antibiotics <- c(
  "Pip/Tazo", "Meropenem", "Cefepime", "Vancomycin_IV", "Vancomycin_PO",
  "Metronidazole", "Ceftriaxone", "Ciprofloxacin", "Clindamycin", "Azithromycin", "TMP/SMX"
)

# =============================================================================
# 3. Filter Drug Table to Systemic Antibiotics Only
# =============================================================================

# Exclude topical, ophthalmic, otic routes (not systemically absorbed)
systemic_routes <- c("IV", "PO", "IM")

DrugTable_systemic <- DrugTable %>%
  filter(Route %in% systemic_routes) %>%
  filter(!Drug %in% c(abx_categories$antifungals, abx_categories$antivirals)) %>%
  mutate(Date = as.Date(Date))

cat("Systemic antibiotic administrations:", nrow(DrugTable_systemic), "\n")
cat("Unique drugs:", length(unique(DrugTable_systemic$Drug)), "\n")

# =============================================================================
# 3b. Create Route-Specific Vancomycin Tables
# =============================================================================
# Vancomycin is special: IV does NOT reach gut lumen, PO stays IN gut lumen
# These have opposite effects on gut microbiome, so we split them

DrugTable_vanc_IV <- DrugTable %>%
  filter(Drug == "Vancomycin", Route == "IV") %>%
  mutate(Date = as.Date(Date))

DrugTable_vanc_PO <- DrugTable %>%
  filter(Drug == "Vancomycin", Route == "PO") %>%
  mutate(Date = as.Date(Date))

cat("\nVancomycin by route:\n")
cat("  IV administrations:", nrow(DrugTable_vanc_IV), "\n")
cat("  PO administrations:", nrow(DrugTable_vanc_PO), "\n")

# =============================================================================
# 4. Create Sample Metadata with Antibiotic Exposure Variables
# =============================================================================

# Start with AvailableDupSamples
sample_metadata <- AvailableDupSamples %>%
  mutate(
    SampleDate = as.Date(SampleDate),
    sample_id = SequenceFileName
  ) %>%
  select(sample_id, MRN, PatientID, SampleDate, SampleType, PatientGroup)

# Function to calculate antibiotic exposure for a sample
calculate_abx_exposure <- function(mrn, sample_date, drug_data, window_days = 7) {

  start_date <- sample_date - window_days
  end_date <- sample_date - 1  # Exclude day of sample

  # Get drugs in window
  drugs_in_window <- drug_data %>%
    filter(MRN == mrn, Date >= start_date, Date <= end_date)

  if (nrow(drugs_in_window) == 0) {
    return(list(
      abx_any = FALSE,
      abx_days = 0,
      abx_anaerobic = FALSE,
      abx_gram_pos = FALSE,
      abx_broad = FALSE,
      abx_pseudomonal = FALSE,
      drugs_list = ""
    ))
  }

  unique_drugs <- unique(drugs_in_window$Drug)
  abx_days <- length(unique(drugs_in_window$Date))

  list(
    abx_any = TRUE,
    abx_days = abx_days,
    abx_anaerobic = any(unique_drugs %in% abx_categories$anti_anaerobic),
    abx_gram_pos = any(unique_drugs %in% abx_categories$anti_gram_positive),
    abx_broad = any(unique_drugs %in% abx_categories$broad_spectrum),
    abx_pseudomonal = any(unique_drugs %in% abx_categories$anti_pseudomonal),
    drugs_list = paste(unique_drugs, collapse = "; ")
  )
}

cat("Calculating antibiotic exposure for each sample (7-day window)...\n")

# Calculate 7-day exposure for all samples
exposure_7d <- sample_metadata %>%
  rowwise() %>%
  mutate(
    exposure = list(calculate_abx_exposure(MRN, SampleDate, DrugTable_systemic, 7))
  ) %>%
  ungroup() %>%
  mutate(
    abx_any_7d = map_lgl(exposure, ~ .x$abx_any),
    abx_days_7d = map_dbl(exposure, ~ .x$abx_days),
    abx_anaerobic_7d = map_lgl(exposure, ~ .x$abx_anaerobic),
    abx_gram_pos_7d = map_lgl(exposure, ~ .x$abx_gram_pos),
    abx_broad_7d = map_lgl(exposure, ~ .x$abx_broad),
    abx_pseudomonal_7d = map_lgl(exposure, ~ .x$abx_pseudomonal),
    drugs_7d = map_chr(exposure, ~ .x$drugs_list)
  ) %>%
  select(-exposure)

cat("Calculating antibiotic exposure for each sample (14-day window)...\n")

# Calculate 14-day exposure
exposure_14d <- sample_metadata %>%
  rowwise() %>%
  mutate(
    exposure = list(calculate_abx_exposure(MRN, SampleDate, DrugTable_systemic, 14))
  ) %>%
  ungroup() %>%
  mutate(
    abx_any_14d = map_lgl(exposure, ~ .x$abx_any),
    abx_days_14d = map_dbl(exposure, ~ .x$abx_days),
    abx_anaerobic_14d = map_lgl(exposure, ~ .x$abx_anaerobic),
    abx_broad_14d = map_lgl(exposure, ~ .x$abx_broad)
  ) %>%
  select(sample_id, abx_any_14d, abx_days_14d, abx_anaerobic_14d, abx_broad_14d)

# Merge 7d and 14d exposure
sample_metadata <- exposure_7d %>%
  left_join(exposure_14d, by = "sample_id")

# Add individual key antibiotic exposures (7-day)
cat("Calculating individual antibiotic exposures...\n")

# Function to calculate route-specific vancomycin exposure
calc_vanc_route_exposure <- function(mrn, sample_date, drug_table, window_days = 7) {
  start_date <- sample_date - window_days
  end_date <- sample_date - 1

  drugs_in_window <- drug_table %>%
    filter(MRN == mrn, Date >= start_date, Date <= end_date)

  nrow(drugs_in_window) > 0
}

for (abx in key_antibiotics) {
  col_name <- paste0(gsub("/", "_", abx), "_7d")

  if (abx == "Vancomycin_IV") {
    # Calculate IV vancomycin exposure from route-specific table
    sample_metadata[[col_name]] <- mapply(
      calc_vanc_route_exposure,
      sample_metadata$MRN,
      sample_metadata$SampleDate,
      MoreArgs = list(drug_table = DrugTable_vanc_IV, window_days = 7)
    )
  } else if (abx == "Vancomycin_PO") {
    # Calculate PO vancomycin exposure from route-specific table
    sample_metadata[[col_name]] <- mapply(
      calc_vanc_route_exposure,
      sample_metadata$MRN,
      sample_metadata$SampleDate,
      MoreArgs = list(drug_table = DrugTable_vanc_PO, window_days = 7)
    )
  } else {
    # Standard grep-based approach for other antibiotics
    sample_metadata <- sample_metadata %>%
      mutate(
        !!col_name := map_lgl(drugs_7d, ~ grepl(gsub("_", "/", abx), .x, fixed = TRUE))
      )
  }
}

# Summary of vancomycin route-specific exposures
cat("\nVancomycin exposure summary:\n")
cat("  Vancomycin_IV_7d:", sum(sample_metadata$Vancomycin_IV_7d), "samples\n")
cat("  Vancomycin_PO_7d:", sum(sample_metadata$Vancomycin_PO_7d), "samples\n")
cat("  Both IV and PO:", sum(sample_metadata$Vancomycin_IV_7d & sample_metadata$Vancomycin_PO_7d), "samples\n")

# =============================================================================
# 5. Prepare Species and Genus Abundance Matrices
# =============================================================================

cat("Preparing abundance matrices...\n")

# Species matrix (remove the 'Sample' column, already in row names)
species_matrix <- SpeciesTableNR %>%
  select(-Sample) %>%
  as.matrix()
rownames(species_matrix) <- SpeciesTableNR$Sample

# Filter to samples in our metadata
common_samples <- intersect(rownames(species_matrix), sample_metadata$sample_id)
species_matrix <- species_matrix[common_samples, ]
sample_metadata <- sample_metadata %>% filter(sample_id %in% common_samples)

cat("Species matrix:", nrow(species_matrix), "samples x", ncol(species_matrix), "species\n")

# Genus matrix
genus_matrix <- GenusTableNR %>%
  select(-Sample) %>%
  as.matrix()
rownames(genus_matrix) <- GenusTableNR$Sample

# Filter to common samples
common_samples_genus <- intersect(rownames(genus_matrix), sample_metadata$sample_id)
genus_matrix_filtered <- genus_matrix[common_samples_genus, ]

cat("Genus matrix:", nrow(genus_matrix_filtered), "samples x", ncol(genus_matrix_filtered), "genera\n")

# =============================================================================
# 6. Create Functional Taxonomic Groups
# =============================================================================

cat("Creating functional taxonomic groups...\n")

# Define genus groups
obligate_anaerobes <- c(
  "Bacteroides", "Parabacteroides", "Prevotella", "Clostridium",
  "Faecalibacterium", "Blautia", "Roseburia", "Ruminococcus",
  "Alistipes", "Odoribacter", "Coprococcus", "Anaerostipes",
  "Eubacterium", "Dorea", "Lachnospira", "Oscillibacter"
)

enterobacteriaceae <- c(
  "Escherichia", "Klebsiella", "Enterobacter", "Citrobacter",
  "Serratia", "Proteus", "Morganella", "Salmonella", "Shigella"
)

enterococcus_genus <- c("Enterococcus")

lactobacillales <- c(
  "Lactobacillus", "Streptococcus", "Leuconostoc", "Pediococcus",
  "Lactococcus", "Weissella"
)

# Function to sum abundances for a genus group
sum_group_abundance <- function(matrix, genera_list) {
  # Find matching columns (genus names may have prefixes/suffixes)
  matching_cols <- sapply(genera_list, function(g) {
    grep(paste0("^", g, "$|^", g, "\\."), colnames(matrix), value = TRUE)
  })
  matching_cols <- unique(unlist(matching_cols))

  if (length(matching_cols) == 0) return(rep(0, nrow(matrix)))
  if (length(matching_cols) == 1) return(matrix[, matching_cols])
  rowSums(matrix[, matching_cols, drop = FALSE])
}

# Calculate group abundances (using genus matrix)
functional_groups <- data.frame(
  sample_id = rownames(genus_matrix_filtered),
  anaerobes_abs = sum_group_abundance(genus_matrix_filtered, obligate_anaerobes),
  enterobact_abs = sum_group_abundance(genus_matrix_filtered, enterobacteriaceae),
  enterococcus_abs = sum_group_abundance(genus_matrix_filtered, enterococcus_genus),
  lactobacillales_abs = sum_group_abundance(genus_matrix_filtered, lactobacillales)
)

# Calculate relative abundances
total_reads <- rowSums(genus_matrix_filtered)
functional_groups <- functional_groups %>%
  mutate(
    anaerobes_rel = anaerobes_abs / total_reads,
    enterobact_rel = enterobact_abs / total_reads,
    enterococcus_rel = enterococcus_abs / total_reads,
    lactobacillales_rel = lactobacillales_abs / total_reads
  )

# Merge with sample metadata
sample_metadata <- sample_metadata %>%
  left_join(functional_groups, by = "sample_id")

# =============================================================================
# 7. Calculate Alpha Diversity
# =============================================================================

cat("Calculating alpha diversity metrics...\n")

# Convert to relative abundance for diversity calculation
species_rel <- species_matrix / rowSums(species_matrix)

diversity_metrics <- data.frame(
  sample_id = rownames(species_matrix),
  shannon = diversity(species_matrix, index = "shannon"),
  simpson = diversity(species_matrix, index = "simpson"),
  richness = specnumber(species_matrix),
  evenness = diversity(species_matrix) / log(specnumber(species_matrix))
)

sample_metadata <- sample_metadata %>%
  left_join(diversity_metrics, by = "sample_id")

# =============================================================================
# 8. Data Quality Assessment by Patient Group
# =============================================================================

cat("\n=== Data Quality Assessment by Patient Group ===\n\n")

quality_summary <- sample_metadata %>%
  group_by(PatientGroup) %>%
  summarise(
    n_samples = n(),
    n_patients = n_distinct(MRN),
    samples_per_patient = round(n_samples / n_patients, 1),
    pct_abx_7d = round(100 * mean(abx_any_7d), 1),
    pct_abx_14d = round(100 * mean(abx_any_14d), 1),
    mean_shannon = round(mean(shannon, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(desc(n_samples))

print(quality_summary)

# Define high-confidence groups (inpatient, complete Abx data expected)
# Include SB patients but only Stool samples (exclude Fistula/Ostomy)
high_confidence_groups <- c("BMT", "IF", "LvTx", "PICU", "SB")

sample_metadata <- sample_metadata %>%
  mutate(
    high_confidence = (PatientGroup %in% high_confidence_groups) &
                      (SampleType == "Stool" | PatientGroup != "SB")
  )

cat("\n=== High-Confidence Subset ===\n")
cat("Groups:", paste(high_confidence_groups, collapse = ", "), "\n")
cat("Note: SB limited to Stool samples only (excluding Fistula/Ostomy)\n")
cat("Samples:", sum(sample_metadata$high_confidence), "\n")
cat("Patients:", n_distinct(sample_metadata$MRN[sample_metadata$high_confidence]), "\n")

# =============================================================================
# 9. Prepare Paired Sample Data
# =============================================================================

cat("\nPreparing paired sample data...\n")

# Use existing TwoWeekSamplePairs and enhance with exposure data
paired_samples <- TwoWeekSamplePairs %>%
  mutate(
    sample1 = SequenceFileName,
    sample2 = SequenceFileName.1,
    date1 = as.Date(SampleDate),
    date2 = as.Date(SampleDate2),
    interval_days = as.numeric(DateDiff)
  ) %>%
  select(PairNumber, MRN, PatientID, PatientGroup, sample1, sample2,
         date1, date2, interval_days)

# Calculate antibiotic exposure BETWEEN paired samples
# Now includes individual antibiotic exposures (days of exposure) for granular analysis
calculate_between_exposure <- function(mrn, date1, date2, drug_data,
                                        vanc_iv_data, vanc_po_data) {

  drugs_between <- drug_data %>%
    filter(MRN == mrn, Date > date1, Date < date2)

  # Route-specific vancomycin
  vanc_iv_between <- vanc_iv_data %>%
    filter(MRN == mrn, Date > date1, Date < date2)
  vanc_po_between <- vanc_po_data %>%
    filter(MRN == mrn, Date > date1, Date < date2)

  if (nrow(drugs_between) == 0) {
    result <- list(
      abx_between_any = FALSE,
      abx_between_days = 0,
      abx_between_anaerobic = FALSE,
      abx_between_broad = FALSE,
      drugs_between = "",
      # Individual antibiotics (days of exposure)
      Pip_Tazo_between = 0,
      Meropenem_between = 0,
      Cefepime_between = 0,
      Ceftriaxone_between = 0,
      Ciprofloxacin_between = 0,
      Metronidazole_between = 0,
      Clindamycin_between = 0,
      TMP_SMX_between = 0,
      Vancomycin_IV_between = 0,
      Vancomycin_PO_between = 0
    )
    return(result)
  }

  unique_drugs <- unique(drugs_between$Drug)
  abx_days <- length(unique(drugs_between$Date))

  # Calculate days of exposure for each individual antibiotic
  calc_drug_days <- function(drug_pattern) {
    matching <- drugs_between %>% filter(grepl(drug_pattern, Drug, fixed = TRUE))
    length(unique(matching$Date))
  }

  list(
    abx_between_any = TRUE,
    abx_between_days = abx_days,
    abx_between_anaerobic = any(unique_drugs %in% abx_categories$anti_anaerobic),
    abx_between_broad = any(unique_drugs %in% abx_categories$broad_spectrum),
    drugs_between = paste(unique_drugs, collapse = "; "),
    # Individual antibiotics (days of exposure between samples)
    Pip_Tazo_between = calc_drug_days("Pip/Tazo"),
    Meropenem_between = calc_drug_days("Meropenem"),
    Cefepime_between = calc_drug_days("Cefepime"),
    Ceftriaxone_between = calc_drug_days("Ceftriaxone"),
    Ciprofloxacin_between = calc_drug_days("Ciprofloxacin"),
    Metronidazole_between = calc_drug_days("Metronidazole"),
    Clindamycin_between = calc_drug_days("Clindamycin"),
    TMP_SMX_between = calc_drug_days("TMP/SMX"),
    Vancomycin_IV_between = length(unique(vanc_iv_between$Date)),
    Vancomycin_PO_between = length(unique(vanc_po_between$Date))
  )
}

cat("Calculating between-sample antibiotic exposure for pairs...\n")

paired_samples <- paired_samples %>%
  rowwise() %>%
  mutate(
    exposure = list(calculate_between_exposure(
      MRN, date1, date2, DrugTable_systemic, DrugTable_vanc_IV, DrugTable_vanc_PO
    ))
  ) %>%
  ungroup() %>%
  mutate(
    # Legacy category-based exposures (kept for backwards compatibility)
    abx_between_any = map_lgl(exposure, ~ .x$abx_between_any),
    abx_between_days = map_dbl(exposure, ~ .x$abx_between_days),
    abx_between_anaerobic = map_lgl(exposure, ~ .x$abx_between_anaerobic),
    abx_between_broad = map_lgl(exposure, ~ .x$abx_between_broad),
    drugs_between = map_chr(exposure, ~ .x$drugs_between),
    # Individual antibiotic exposures (days between samples)
    Pip_Tazo_between = map_dbl(exposure, ~ .x$Pip_Tazo_between),
    Meropenem_between = map_dbl(exposure, ~ .x$Meropenem_between),
    Cefepime_between = map_dbl(exposure, ~ .x$Cefepime_between),
    Ceftriaxone_between = map_dbl(exposure, ~ .x$Ceftriaxone_between),
    Ciprofloxacin_between = map_dbl(exposure, ~ .x$Ciprofloxacin_between),
    Metronidazole_between = map_dbl(exposure, ~ .x$Metronidazole_between),
    Clindamycin_between = map_dbl(exposure, ~ .x$Clindamycin_between),
    TMP_SMX_between = map_dbl(exposure, ~ .x$TMP_SMX_between),
    Vancomycin_IV_between = map_dbl(exposure, ~ .x$Vancomycin_IV_between),
    Vancomycin_PO_between = map_dbl(exposure, ~ .x$Vancomycin_PO_between)
  ) %>%
  select(-exposure)

# Print summary of individual antibiotic exposures between pairs
cat("\nIndividual antibiotic exposures between paired samples:\n")
indiv_abx <- c("Pip_Tazo", "Meropenem", "Cefepime", "Ceftriaxone", "Ciprofloxacin",
               "Metronidazole", "Clindamycin", "TMP_SMX", "Vancomycin_IV", "Vancomycin_PO")
for (abx in indiv_abx) {
  col <- paste0(abx, "_between")
  n_exposed <- sum(paired_samples[[col]] > 0)
  cat(sprintf("  %s: %d pairs (%.1f%%)\n", abx, n_exposed,
              100 * n_exposed / nrow(paired_samples)))
}

# Add high-confidence flag
paired_samples <- paired_samples %>%
  mutate(high_confidence = PatientGroup %in% high_confidence_groups)

cat("Total pairs:", nrow(paired_samples), "\n")
cat("Pairs with Abx exposure:", sum(paired_samples$abx_between_any), "\n")
cat("High-confidence pairs:", sum(paired_samples$high_confidence), "\n")

# =============================================================================
# 10. Save Prepared Data
# =============================================================================

cat("\nSaving prepared data...\n")

# Save quality summary
write_csv(quality_summary, file.path(results_dir, "data_quality/sample_counts_by_group.csv"))

# Save key objects for downstream analysis
prepared_data <- list(
  sample_metadata = sample_metadata,
  species_matrix = species_matrix,
  genus_matrix = genus_matrix_filtered,
  paired_samples = paired_samples,
  abx_categories = abx_categories,
  key_antibiotics = key_antibiotics,
  high_confidence_groups = high_confidence_groups
)

save(prepared_data, file = file.path(project_dir, "data/prepared_data.RData"))

cat("\nData preparation complete!\n")
cat("Saved to: data/prepared_data.RData\n")

# =============================================================================
# 11. Print Summary Statistics
# =============================================================================

cat("\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("SUMMARY\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("\nSamples:", nrow(sample_metadata), "\n")
cat("Unique patients:", n_distinct(sample_metadata$MRN), "\n")
cat("Sample pairs (2-week):", nrow(paired_samples), "\n")
cat("\nAntibiotic exposure (7-day window):\n")
cat("  Any Abx:", sum(sample_metadata$abx_any_7d),
    "(", round(100*mean(sample_metadata$abx_any_7d), 1), "%)\n")
cat("  Anti-anaerobic:", sum(sample_metadata$abx_anaerobic_7d),
    "(", round(100*mean(sample_metadata$abx_anaerobic_7d), 1), "%)\n")
cat("  Broad-spectrum:", sum(sample_metadata$abx_broad_7d),
    "(", round(100*mean(sample_metadata$abx_broad_7d), 1), "%)\n")
cat("\nDiversity (mean Shannon):", round(mean(sample_metadata$shannon), 2), "\n")
