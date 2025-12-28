#!/usr/bin/env Rscript
# =============================================================================
# 02_diversity_analysis.R
# Alpha and beta diversity analysis with antibiotic exposure
# =============================================================================

library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
library(broom.mixed)

# Set paths
project_dir <- "/home/david/projects/abx_metagenomic_correlations"
results_dir <- file.path(project_dir, "results/new_analysis_legacy_data")
figures_dir <- file.path(results_dir, "figures")

# Load prepared data
load(file.path(project_dir, "data/prepared_data.RData"))

sample_metadata <- prepared_data$sample_metadata
species_matrix <- prepared_data$species_matrix
high_confidence_groups <- prepared_data$high_confidence_groups

# =============================================================================
# 1. Alpha Diversity Analysis
# =============================================================================

cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("ALPHA DIVERSITY ANALYSIS\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")

# -----------------------------------------------------------------------------
# 1.1 Descriptive Statistics
# -----------------------------------------------------------------------------

cat("=== Descriptive Statistics ===\n\n")

# By antibiotic exposure (7-day)
alpha_by_abx <- sample_metadata %>%
  group_by(abx_any_7d) %>%
  summarise(
    n = n(),
    mean_shannon = mean(shannon),
    sd_shannon = sd(shannon),
    median_shannon = median(shannon),
    mean_richness = mean(richness),
    .groups = "drop"
  )

cat("Shannon diversity by 7-day Abx exposure:\n")
print(alpha_by_abx)

# By antibiotic spectrum
alpha_by_spectrum <- sample_metadata %>%
  summarise(
    `No Abx` = mean(shannon[!abx_any_7d]),
    `Any Abx` = mean(shannon[abx_any_7d]),
    `Anti-anaerobic` = mean(shannon[abx_anaerobic_7d]),
    `Broad-spectrum` = mean(shannon[abx_broad_7d]),
    `Gram-positive` = mean(shannon[abx_gram_pos_7d])
  )

cat("\nMean Shannon by antibiotic type:\n")
print(t(alpha_by_spectrum))

# -----------------------------------------------------------------------------
# 1.2 Mixed-Effects Models - All Samples
# -----------------------------------------------------------------------------

cat("\n=== Mixed-Effects Models (All Samples) ===\n\n")

# Model 1: Any antibiotic exposure
model_any <- lmer(
  shannon ~ abx_any_7d + PatientGroup + (1|MRN),
  data = sample_metadata
)

cat("Model: Shannon ~ Abx (any, 7d) + PatientGroup + (1|Patient)\n")
print(summary(model_any))

# Extract coefficients
coef_any <- tidy(model_any, effects = "fixed", conf.int = TRUE)
cat("\nFixed effects:\n")
print(coef_any %>% filter(term != "(Intercept)") %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))))

# Model 2: Antibiotic spectrum
model_spectrum <- lmer(
  shannon ~ abx_anaerobic_7d + abx_broad_7d + abx_gram_pos_7d + PatientGroup + (1|MRN),
  data = sample_metadata
)

cat("\n\nModel: Shannon ~ Abx spectrum + PatientGroup + (1|Patient)\n")
coef_spectrum <- tidy(model_spectrum, effects = "fixed", conf.int = TRUE)
print(coef_spectrum %>% filter(term != "(Intercept)") %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))))

# -----------------------------------------------------------------------------
# 1.3 Mixed-Effects Models - High-Confidence Groups Only
# -----------------------------------------------------------------------------

cat("\n=== Mixed-Effects Models (High-Confidence Groups Only) ===\n\n")
cat("Groups:", paste(high_confidence_groups, collapse = ", "), "\n\n")

sample_hc <- sample_metadata %>% filter(high_confidence)

model_hc <- lmer(
  shannon ~ abx_any_7d + PatientGroup + (1|MRN),
  data = sample_hc
)

cat("Model: Shannon ~ Abx (any, 7d) + PatientGroup + (1|Patient)\n")
coef_hc <- tidy(model_hc, effects = "fixed", conf.int = TRUE)
print(coef_hc %>% filter(term != "(Intercept)") %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))))

model_hc_spectrum <- lmer(
  shannon ~ abx_anaerobic_7d + abx_broad_7d + PatientGroup + (1|MRN),
  data = sample_hc
)

cat("\nModel: Shannon ~ Abx spectrum + PatientGroup + (1|Patient)\n")
coef_hc_spectrum <- tidy(model_hc_spectrum, effects = "fixed", conf.int = TRUE)
print(coef_hc_spectrum %>% filter(term != "(Intercept)") %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))))

# -----------------------------------------------------------------------------
# 1.4 Taxonomic Group Analysis
# -----------------------------------------------------------------------------

cat("\n=== Functional Taxonomic Groups ===\n\n")

# Anaerobe relative abundance by antibiotic exposure
cat("Anaerobe relative abundance:\n")
sample_metadata %>%
  group_by(abx_anaerobic_7d) %>%
  summarise(
    n = n(),
    mean_anaerobes = mean(anaerobes_rel, na.rm = TRUE),
    median_anaerobes = median(anaerobes_rel, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# Model for anaerobe abundance
model_anaerobes <- lmer(
  log10(anaerobes_rel + 1e-6) ~ abx_anaerobic_7d + PatientGroup + (1|MRN),
  data = sample_hc
)

cat("\nModel: log10(Anaerobes) ~ Anti-anaerobic Abx + PatientGroup + (1|Patient)\n")
coef_anaerobes <- tidy(model_anaerobes, effects = "fixed", conf.int = TRUE)
print(coef_anaerobes %>% filter(term != "(Intercept)") %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))))

# Enterobacteriaceae by broad-spectrum exposure
cat("\nEnterobacteriaceae relative abundance:\n")
sample_metadata %>%
  group_by(abx_broad_7d) %>%
  summarise(
    n = n(),
    mean_enterobact = mean(enterobact_rel, na.rm = TRUE),
    median_enterobact = median(enterobact_rel, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

model_enterobact <- lmer(
  log10(enterobact_rel + 1e-6) ~ abx_broad_7d + PatientGroup + (1|MRN),
  data = sample_hc
)

cat("\nModel: log10(Enterobacteriaceae) ~ Broad-spectrum Abx + PatientGroup + (1|Patient)\n")
coef_enterobact <- tidy(model_enterobact, effects = "fixed", conf.int = TRUE)
print(coef_enterobact %>% filter(term != "(Intercept)") %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))))

# =============================================================================
# 2. Beta Diversity Analysis
# =============================================================================

cat("\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("BETA DIVERSITY ANALYSIS\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")

# -----------------------------------------------------------------------------
# 2.1 Calculate Bray-Curtis Distance
# -----------------------------------------------------------------------------

# Relative abundance
species_rel <- species_matrix / rowSums(species_matrix)

# Bray-Curtis distance
bray_dist <- vegdist(species_rel, method = "bray")

# -----------------------------------------------------------------------------
# 2.2 PERMANOVA - All Samples
# -----------------------------------------------------------------------------

cat("=== PERMANOVA (All Samples) ===\n\n")

# Align metadata with distance matrix
meta_aligned <- sample_metadata[match(rownames(species_rel), sample_metadata$sample_id), ]

set.seed(42)
permanova_all <- adonis2(
  bray_dist ~ abx_any_7d + PatientGroup,
  data = meta_aligned,
  permutations = 999,
  by = "margin"
)

cat("Model: Bray-Curtis ~ Abx (any, 7d) + PatientGroup\n")
print(permanova_all)

# With spectrum
permanova_spectrum <- adonis2(
  bray_dist ~ abx_anaerobic_7d + abx_broad_7d + PatientGroup,
  data = meta_aligned,
  permutations = 999,
  by = "margin"
)

cat("\nModel: Bray-Curtis ~ Abx spectrum + PatientGroup\n")
print(permanova_spectrum)

# -----------------------------------------------------------------------------
# 2.3 PERMANOVA - High-Confidence Groups
# -----------------------------------------------------------------------------

cat("\n=== PERMANOVA (High-Confidence Groups) ===\n\n")

# Subset to high-confidence samples
hc_samples <- sample_metadata$sample_id[sample_metadata$high_confidence]
hc_idx <- which(rownames(species_rel) %in% hc_samples)

species_rel_hc <- species_rel[hc_idx, ]
bray_dist_hc <- vegdist(species_rel_hc, method = "bray")
meta_hc_aligned <- sample_metadata[match(rownames(species_rel_hc), sample_metadata$sample_id), ]

set.seed(42)
permanova_hc <- adonis2(
  bray_dist_hc ~ abx_any_7d + PatientGroup,
  data = meta_hc_aligned,
  permutations = 999,
  by = "margin"
)

cat("Model: Bray-Curtis ~ Abx (any, 7d) + PatientGroup\n")
print(permanova_hc)

permanova_hc_spectrum <- adonis2(
  bray_dist_hc ~ abx_anaerobic_7d + abx_broad_7d + PatientGroup,
  data = meta_hc_aligned,
  permutations = 999,
  by = "margin"
)

cat("\nModel: Bray-Curtis ~ Abx spectrum + PatientGroup\n")
print(permanova_hc_spectrum)

# -----------------------------------------------------------------------------
# 2.4 PCoA Ordination
# -----------------------------------------------------------------------------

cat("\n=== PCoA Ordination ===\n\n")

pcoa <- cmdscale(bray_dist_hc, k = 2, eig = TRUE)
variance_explained <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 1)

cat("Variance explained: PC1 =", variance_explained[1], "%, PC2 =", variance_explained[2], "%\n")

pcoa_df <- data.frame(
  sample_id = rownames(pcoa$points),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
) %>%
  left_join(sample_metadata, by = "sample_id")

# =============================================================================
# 3. Generate Plots
# =============================================================================

cat("\n=== Generating Plots ===\n\n")

# -----------------------------------------------------------------------------
# 3.1 Alpha Diversity Boxplot
# -----------------------------------------------------------------------------

p_alpha <- ggplot(sample_hc, aes(x = abx_any_7d, y = shannon, fill = abx_any_7d)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  facet_wrap(~ PatientGroup, nrow = 1) +
  scale_fill_manual(
    values = c("FALSE" = "#4575b4", "TRUE" = "#d73027"),
    labels = c("FALSE" = "No Abx", "TRUE" = "Abx (7d)")
  ) +
  labs(
    title = "Shannon Diversity by Antibiotic Exposure",
    subtitle = "High-confidence patient groups (BMT, IF, LvTx, PICU)",
    x = "Recent Antibiotic Exposure (7 days)",
    y = "Shannon Diversity Index",
    fill = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90")
  )

ggsave(file.path(figures_dir, "alpha_diversity_by_abx.pdf"), p_alpha,
       width = 10, height = 5)
cat("Saved: figures/alpha_diversity_by_abx.pdf\n")

# -----------------------------------------------------------------------------
# 3.2 Alpha Diversity by Spectrum
# -----------------------------------------------------------------------------

spectrum_long <- sample_hc %>%
  select(sample_id, shannon, abx_anaerobic_7d, abx_broad_7d, abx_gram_pos_7d) %>%
  pivot_longer(
    cols = c(abx_anaerobic_7d, abx_broad_7d, abx_gram_pos_7d),
    names_to = "abx_type",
    values_to = "exposed"
  ) %>%
  mutate(
    abx_type = case_when(
      abx_type == "abx_anaerobic_7d" ~ "Anti-anaerobic",
      abx_type == "abx_broad_7d" ~ "Broad-spectrum",
      abx_type == "abx_gram_pos_7d" ~ "Anti-gram-positive"
    )
  )

p_spectrum <- ggplot(spectrum_long, aes(x = exposed, y = shannon, fill = exposed)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.5) +
  facet_wrap(~ abx_type) +
  scale_fill_manual(
    values = c("FALSE" = "#4575b4", "TRUE" = "#d73027"),
    labels = c("FALSE" = "Not exposed", "TRUE" = "Exposed")
  ) +
  labs(
    title = "Shannon Diversity by Antibiotic Spectrum",
    x = "Exposed (7-day window)",
    y = "Shannon Diversity Index",
    fill = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(figures_dir, "alpha_diversity_by_spectrum.pdf"), p_spectrum,
       width = 9, height = 4)
cat("Saved: figures/alpha_diversity_by_spectrum.pdf\n")

# -----------------------------------------------------------------------------
# 3.3 PCoA Plot
# -----------------------------------------------------------------------------

p_pcoa <- ggplot(pcoa_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = abx_any_7d, shape = PatientGroup), alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c("FALSE" = "#4575b4", "TRUE" = "#d73027"),
    labels = c("FALSE" = "No Abx", "TRUE" = "Abx (7d)")
  ) +
  labs(
    title = "PCoA of Microbiome Composition (Bray-Curtis)",
    subtitle = "High-confidence patient groups",
    x = paste0("PC1 (", variance_explained[1], "%)"),
    y = paste0("PC2 (", variance_explained[2], "%)"),
    color = "Abx Exposure",
    shape = "Patient Group"
  ) +
  theme_bw() +
  theme(legend.position = "right")

ggsave(file.path(figures_dir, "pcoa_by_abx.pdf"), p_pcoa,
       width = 8, height = 6)
cat("Saved: figures/pcoa_by_abx.pdf\n")

# -----------------------------------------------------------------------------
# 3.4 Functional Groups Boxplot
# -----------------------------------------------------------------------------

functional_long <- sample_hc %>%
  select(sample_id, abx_anaerobic_7d, abx_broad_7d,
         anaerobes_rel, enterobact_rel, enterococcus_rel) %>%
  pivot_longer(
    cols = c(anaerobes_rel, enterobact_rel, enterococcus_rel),
    names_to = "group",
    values_to = "rel_abundance"
  ) %>%
  mutate(
    group = case_when(
      group == "anaerobes_rel" ~ "Obligate Anaerobes",
      group == "enterobact_rel" ~ "Enterobacteriaceae",
      group == "enterococcus_rel" ~ "Enterococcus"
    )
  )

p_functional <- ggplot(
  functional_long %>% filter(group == "Obligate Anaerobes"),
  aes(x = abx_anaerobic_7d, y = rel_abundance + 1e-6, fill = abx_anaerobic_7d)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_y_log10() +
  scale_fill_manual(
    values = c("FALSE" = "#4575b4", "TRUE" = "#d73027"),
    labels = c("FALSE" = "No anti-anaerobic Abx", "TRUE" = "Anti-anaerobic Abx")
  ) +
  labs(
    title = "Obligate Anaerobe Abundance by Anti-Anaerobic Antibiotic Exposure",
    x = "Anti-anaerobic Abx (7-day window)",
    y = "Relative Abundance (log scale)",
    fill = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(figures_dir, "anaerobes_by_anti_anaerobic_abx.pdf"), p_functional,
       width = 6, height = 5)
cat("Saved: figures/anaerobes_by_anti_anaerobic_abx.pdf\n")

# =============================================================================
# 4. Save Results
# =============================================================================

cat("\n=== Saving Results ===\n\n")

# Compile alpha diversity results
alpha_results <- list(
  model_any_all = tidy(model_any, effects = "fixed", conf.int = TRUE),
  model_spectrum_all = tidy(model_spectrum, effects = "fixed", conf.int = TRUE),
  model_any_hc = tidy(model_hc, effects = "fixed", conf.int = TRUE),
  model_spectrum_hc = tidy(model_hc_spectrum, effects = "fixed", conf.int = TRUE),
  model_anaerobes = tidy(model_anaerobes, effects = "fixed", conf.int = TRUE),
  model_enterobact = tidy(model_enterobact, effects = "fixed", conf.int = TRUE)
)

# Save to CSV
alpha_summary <- bind_rows(
  alpha_results$model_any_hc %>% mutate(model = "Shannon ~ Any Abx (HC)"),
  alpha_results$model_spectrum_hc %>% mutate(model = "Shannon ~ Spectrum (HC)"),
  alpha_results$model_anaerobes %>% mutate(model = "Anaerobes ~ Anti-anaerobic (HC)"),
  alpha_results$model_enterobact %>% mutate(model = "Enterobact ~ Broad-spectrum (HC)")
) %>%
  filter(!grepl("Intercept", term)) %>%
  select(model, term, estimate, std.error, conf.low, conf.high, p.value) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

write_csv(alpha_summary, file.path(results_dir, "diversity/alpha_diversity_results.csv"))
cat("Saved: diversity/alpha_diversity_results.csv\n")

# Save PERMANOVA results
permanova_results <- bind_rows(
  as.data.frame(permanova_hc) %>%
    rownames_to_column("term") %>%
    mutate(model = "Any Abx (HC)"),
  as.data.frame(permanova_hc_spectrum) %>%
    rownames_to_column("term") %>%
    mutate(model = "Spectrum (HC)")
)

write_csv(permanova_results, file.path(results_dir, "diversity/beta_diversity_permanova.csv"))
cat("Saved: diversity/beta_diversity_permanova.csv\n")

cat("\nDiversity analysis complete!\n")
