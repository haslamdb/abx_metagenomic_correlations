#!/usr/bin/env Rscript
#' Paired Sample Analysis
#'
#' Analyzes within-patient microbiome changes associated with antibiotic exposure.
#' This is the core analysis that addresses the limitations of simple correlations.

library(vegan)
library(tidyverse)
library(lme4)
library(lmerTest)
library(yaml)

# Load configuration
config <- yaml::read_yaml("config/config.yaml")

# Load data from Snakemake
species_file <- snakemake@input[["species"]]
genus_file <- snakemake@input[["genus"]]
metadata_file <- snakemake@input[["metadata"]]
drug_file <- snakemake@input[["drug_exposure"]]
output_results <- snakemake@output[["results"]]
output_summary <- snakemake@output[["summary"]]

# Read data
species <- read.delim(species_file, row.names = 1, check.names = FALSE)
genus <- read.delim(genus_file, row.names = 1, check.names = FALSE)
metadata <- read.delim(metadata_file, stringsAsFactors = FALSE)
drug_exposure <- read.delim(drug_file, stringsAsFactors = FALSE)

# Parse dates
metadata$sample_date <- as.Date(metadata$sample_date)
drug_exposure$date <- as.Date(drug_exposure$date)

# Define taxonomic groups from config
obligate_anaerobes <- config$taxonomic_groups$obligate_anaerobes
enterobacteriaceae <- config$taxonomic_groups$enterobacteriaceae
enterococcus <- config$taxonomic_groups$enterococcus

# Define antibiotic classes from config
anti_anaerobic <- config$antibiotic_classes$anti_anaerobic
broad_spectrum <- config$antibiotic_classes$broad_spectrum
anti_gram_positive <- config$antibiotic_classes$anti_gram_positive

# Calculate taxonomic group abundances
calculate_group_abundance <- function(abundance_matrix, group_genera) {
  # Find matching columns (genus names)
  matching_cols <- colnames(abundance_matrix)[
    sapply(colnames(abundance_matrix), function(x) {
      any(sapply(group_genera, function(g) grepl(paste0("^", g), x, ignore.case = TRUE)))
    })
  ]

  if (length(matching_cols) == 0) {
    return(rep(0, nrow(abundance_matrix)))
  }

  rowSums(abundance_matrix[, matching_cols, drop = FALSE])
}

# Calculate alpha diversity
shannon <- diversity(species, index = "shannon")
richness <- specnumber(species)

# Calculate group abundances at genus level
anaerobe_abundance <- calculate_group_abundance(genus, obligate_anaerobes)
enterobact_abundance <- calculate_group_abundance(genus, enterobacteriaceae)
enterococcus_abundance <- calculate_group_abundance(genus, enterococcus)

# Create sample summary
sample_summary <- data.frame(
  sample_id = rownames(species),
  shannon = shannon,
  richness = richness,
  anaerobe_abundance = anaerobe_abundance,
  enterobact_abundance = enterobact_abundance,
  enterococcus_abundance = enterococcus_abundance,
  stringsAsFactors = FALSE
)

sample_summary <- merge(sample_summary, metadata, by = "sample_id")

# Function to get antibiotic exposure between two dates for a patient
get_exposure <- function(patient, date1, date2, drug_data, abx_classes) {
  patient_drugs <- drug_data %>%
    filter(patient_id == patient,
           date >= date1,
           date <= date2)

  list(
    any_abx = nrow(patient_drugs) > 0,
    n_abx_days = n_distinct(patient_drugs$date),
    anti_anaerobic = any(patient_drugs$drug %in% abx_classes$anti_anaerobic),
    broad_spectrum = any(patient_drugs$drug %in% abx_classes$broad_spectrum),
    anti_gram_positive = any(patient_drugs$drug %in% abx_classes$anti_gram_positive),
    drugs = paste(unique(patient_drugs$drug), collapse = "; ")
  )
}

# Find all valid sample pairs
message("Finding sample pairs...")

# Sort samples by patient and date
sample_summary <- sample_summary %>%
  arrange(patient_id, sample_date)

# Create pairs
pairs_list <- list()
pair_idx <- 1

for (patient in unique(sample_summary$patient_id)) {
  patient_samples <- sample_summary %>%
    filter(patient_id == patient) %>%
    arrange(sample_date)

  if (nrow(patient_samples) < 2) next

  for (i in 1:(nrow(patient_samples) - 1)) {
    for (j in (i + 1):nrow(patient_samples)) {
      days_between <- as.numeric(patient_samples$sample_date[j] - patient_samples$sample_date[i])

      # Check if within valid range
      if (days_between >= config$min_days_between_samples &&
          days_between <= config$max_days_between_samples) {

        sample1 <- patient_samples$sample_id[i]
        sample2 <- patient_samples$sample_id[j]

        # Calculate microbiome changes
        delta_shannon <- patient_samples$shannon[j] - patient_samples$shannon[i]
        delta_richness <- patient_samples$richness[j] - patient_samples$richness[i]

        # Log-fold change for abundances (add pseudocount)
        pseudo <- 1e-6
        delta_anaerobes <- log2((patient_samples$anaerobe_abundance[j] + pseudo) /
                                  (patient_samples$anaerobe_abundance[i] + pseudo))
        delta_enterobact <- log2((patient_samples$enterobact_abundance[j] + pseudo) /
                                   (patient_samples$enterobact_abundance[i] + pseudo))
        delta_enterococcus <- log2((patient_samples$enterococcus_abundance[j] + pseudo) /
                                     (patient_samples$enterococcus_abundance[i] + pseudo))

        # Calculate Bray-Curtis distance between paired samples
        bray_distance <- vegdist(
          rbind(species[sample1, ], species[sample2, ]),
          method = "bray"
        )[1]

        # Get antibiotic exposure
        abx_classes <- list(
          anti_anaerobic = anti_anaerobic,
          broad_spectrum = broad_spectrum,
          anti_gram_positive = anti_gram_positive
        )
        exposure <- get_exposure(
          patient,
          patient_samples$sample_date[i],
          patient_samples$sample_date[j],
          drug_exposure,
          abx_classes
        )

        pairs_list[[pair_idx]] <- data.frame(
          pair_id = pair_idx,
          patient_id = patient,
          sample1 = sample1,
          sample2 = sample2,
          date1 = patient_samples$sample_date[i],
          date2 = patient_samples$sample_date[j],
          days_between = days_between,
          patient_group = patient_samples$patient_group[i],
          delta_shannon = delta_shannon,
          delta_richness = delta_richness,
          delta_anaerobes = delta_anaerobes,
          delta_enterobact = delta_enterobact,
          delta_enterococcus = delta_enterococcus,
          bray_distance = bray_distance,
          any_abx = exposure$any_abx,
          n_abx_days = exposure$n_abx_days,
          anti_anaerobic = exposure$anti_anaerobic,
          broad_spectrum = exposure$broad_spectrum,
          anti_gram_positive = exposure$anti_gram_positive,
          drugs = exposure$drugs,
          stringsAsFactors = FALSE
        )
        pair_idx <- pair_idx + 1
      }
    }
  }
}

paired_results <- do.call(rbind, pairs_list)
message(sprintf("Found %d valid sample pairs from %d patients",
                nrow(paired_results), n_distinct(paired_results$patient_id)))

# Statistical tests
message("Running statistical tests...")

summary_results <- list()

# Test 1: Any antibiotic vs Bray-Curtis distance
if (sum(paired_results$any_abx) >= 5 && sum(!paired_results$any_abx) >= 5) {
  test1 <- wilcox.test(bray_distance ~ any_abx, data = paired_results)
  summary_results$abx_vs_distance <- data.frame(
    comparison = "Any antibiotic vs Bray-Curtis distance",
    group1_n = sum(!paired_results$any_abx),
    group2_n = sum(paired_results$any_abx),
    group1_median = median(paired_results$bray_distance[!paired_results$any_abx]),
    group2_median = median(paired_results$bray_distance[paired_results$any_abx]),
    p_value = test1$p.value,
    test = "Wilcoxon"
  )
}

# Test 2: Anti-anaerobic antibiotics vs anaerobe change
if (sum(paired_results$anti_anaerobic) >= 5 && sum(!paired_results$anti_anaerobic) >= 5) {
  test2 <- wilcox.test(delta_anaerobes ~ anti_anaerobic, data = paired_results)
  summary_results$antianaerobic_vs_anaerobes <- data.frame(
    comparison = "Anti-anaerobic Abx vs Anaerobe change",
    group1_n = sum(!paired_results$anti_anaerobic),
    group2_n = sum(paired_results$anti_anaerobic),
    group1_median = median(paired_results$delta_anaerobes[!paired_results$anti_anaerobic]),
    group2_median = median(paired_results$delta_anaerobes[paired_results$anti_anaerobic]),
    p_value = test2$p.value,
    test = "Wilcoxon"
  )
}

# Test 3: Broad-spectrum antibiotics vs Shannon diversity change
if (sum(paired_results$broad_spectrum) >= 5 && sum(!paired_results$broad_spectrum) >= 5) {
  test3 <- wilcox.test(delta_shannon ~ broad_spectrum, data = paired_results)
  summary_results$broadspectrum_vs_diversity <- data.frame(
    comparison = "Broad-spectrum Abx vs Shannon change",
    group1_n = sum(!paired_results$broad_spectrum),
    group2_n = sum(paired_results$broad_spectrum),
    group1_median = median(paired_results$delta_shannon[!paired_results$broad_spectrum]),
    group2_median = median(paired_results$delta_shannon[paired_results$broad_spectrum]),
    p_value = test3$p.value,
    test = "Wilcoxon"
  )
}

# Test 4: Mixed-effects model for anaerobe change
if (nrow(paired_results) >= 20) {
  tryCatch({
    model_anaerobes <- lmer(
      delta_anaerobes ~ anti_anaerobic + days_between + (1 | patient_id),
      data = paired_results
    )
    model_summary <- summary(model_anaerobes)
    coef_table <- as.data.frame(coef(model_summary))

    summary_results$lmer_anaerobes <- data.frame(
      comparison = "Mixed-effects: Anti-anaerobic Abx effect on anaerobes",
      estimate = coef_table["anti_anaerobicTRUE", "Estimate"],
      std_error = coef_table["anti_anaerobicTRUE", "Std. Error"],
      t_value = coef_table["anti_anaerobicTRUE", "t value"],
      p_value = coef_table["anti_anaerobicTRUE", "Pr(>|t|)"],
      test = "lmer"
    )
  }, error = function(e) {
    message("Mixed-effects model failed: ", e$message)
  })
}

# Combine summary results
summary_df <- do.call(rbind, lapply(summary_results, function(x) {
  # Ensure consistent columns
  x$estimate <- if ("estimate" %in% names(x)) x$estimate else NA
  x$std_error <- if ("std_error" %in% names(x)) x$std_error else NA
  x$t_value <- if ("t_value" %in% names(x)) x$t_value else NA
  x$group1_n <- if ("group1_n" %in% names(x)) x$group1_n else NA
  x$group2_n <- if ("group2_n" %in% names(x)) x$group2_n else NA
  x$group1_median <- if ("group1_median" %in% names(x)) x$group1_median else NA
  x$group2_median <- if ("group2_median" %in% names(x)) x$group2_median else NA
  x[, c("comparison", "group1_n", "group2_n", "group1_median", "group2_median",
        "estimate", "std_error", "t_value", "p_value", "test")]
}))

# Save results
write.table(paired_results, output_results, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(summary_df, output_summary, sep = "\t", row.names = FALSE, quote = FALSE)

message("Paired analysis complete!")
