#!/usr/bin/env Rscript
#' Calculate Alpha Diversity Metrics
#'
#' Calculates Shannon, Simpson, richness, and evenness for all samples.

library(vegan)
library(tidyverse)

# Load data from Snakemake
species_file <- snakemake@input[["species"]]
metadata_file <- snakemake@input[["metadata"]]
output_file <- snakemake@output[["alpha"]]

# Read species abundance matrix
species <- read.delim(species_file, row.names = 1, check.names = FALSE)

# Read metadata
metadata <- read.delim(metadata_file, stringsAsFactors = FALSE)

# Calculate diversity metrics
alpha_diversity <- data.frame(
  sample_id = rownames(species),
  shannon = diversity(species, index = "shannon"),
  simpson = diversity(species, index = "simpson"),
  invsimpson = diversity(species, index = "invsimpson"),
  richness = specnumber(species),
  evenness = diversity(species) / log(specnumber(species)),
  total_reads = rowSums(species)
)

# Handle edge cases (samples with 0 or 1 species)
alpha_diversity$evenness[is.nan(alpha_diversity$evenness)] <- NA

# Merge with metadata
alpha_diversity <- merge(alpha_diversity, metadata, by = "sample_id", all.x = TRUE)

# Save results
write.table(alpha_diversity, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

message(sprintf("Calculated diversity for %d samples", nrow(alpha_diversity)))
