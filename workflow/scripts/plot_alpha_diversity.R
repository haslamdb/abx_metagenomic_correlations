#!/usr/bin/env Rscript
#' Alpha Diversity Visualization
#'
#' Generate boxplots comparing diversity by antibiotic exposure status.

library(tidyverse)
library(ggplot2)

# Load data from Snakemake
alpha_file <- snakemake@input[["alpha"]]
metadata_file <- snakemake@input[["metadata"]]
drug_file <- snakemake@input[["drug_exposure"]]
output_file <- snakemake@output[[1]]

# Read data
alpha <- read.delim(alpha_file, stringsAsFactors = FALSE)
drug_exposure <- read.delim(drug_file, stringsAsFactors = FALSE)

# Parse dates
alpha$sample_date <- as.Date(alpha$sample_date)
drug_exposure$date <- as.Date(drug_exposure$date)

# Calculate 7-day exposure
alpha$abx_7d <- sapply(1:nrow(alpha), function(i) {
  start_date <- alpha$sample_date[i] - 7
  patient_drugs <- drug_exposure %>%
    filter(patient_id == alpha$patient_id[i],
           date >= start_date,
           date <= alpha$sample_date[i])
  nrow(patient_drugs) > 0
})

alpha$abx_status <- ifelse(alpha$abx_7d, "Antibiotics\n(7d prior)", "No antibiotics")

# Create multi-panel figure
p1 <- ggplot(alpha, aes(x = abx_status, y = shannon, fill = abx_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("Antibiotics\n(7d prior)" = "#E64B35",
                                "No antibiotics" = "#4DBBD5")) +
  labs(x = NULL, y = "Shannon Diversity Index") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10))

p2 <- ggplot(alpha, aes(x = abx_status, y = richness, fill = abx_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("Antibiotics\n(7d prior)" = "#E64B35",
                                "No antibiotics" = "#4DBBD5")) +
  labs(x = NULL, y = "Species Richness") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10))

# Statistical tests
test_shannon <- wilcox.test(shannon ~ abx_status, data = alpha)
test_richness <- wilcox.test(richness ~ abx_status, data = alpha)

# Add p-values to plots
p1 <- p1 + annotate("text", x = 1.5, y = max(alpha$shannon) * 0.95,
                    label = sprintf("p = %.3f", test_shannon$p.value))
p2 <- p2 + annotate("text", x = 1.5, y = max(alpha$richness) * 0.95,
                    label = sprintf("p = %.3f", test_richness$p.value))

# Combine plots
library(patchwork)
combined <- p1 + p2 +
  plot_annotation(title = "Alpha Diversity by Antibiotic Exposure Status")

# Save
ggsave(output_file, combined, width = 10, height = 5)

message("Alpha diversity plot saved!")
