#!/usr/bin/env Rscript
# =============================================================================
# 07_cefepime_vs_piptazo_comparison.R
# Compare differential abundance effects of Cefepime vs Piperacillin-Tazobactam
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggVennDiagram)

select <- dplyr::select
filter <- dplyr::filter

project_dir <- "/home/david/projects/abx_metagenomic_correlations"
results_dir <- file.path(project_dir, "results/new_analysis_legacy_data")
indiv_dir <- file.path(results_dir, "individual_antibiotics_v2_with_covariates")
genus_dir <- file.path(results_dir, "individual_antibiotics_genus_level")

# Output directory for comparison
compare_dir <- file.path(results_dir, "cefepime_vs_piptazo")
dir.create(compare_dir, showWarnings = FALSE, recursive = TRUE)

cat("=============================================================================\n")
cat("Cefepime vs Piperacillin-Tazobactam Comparison\n")
cat("=============================================================================\n\n")

# =============================================================================
# 1. Load Combined Results
# =============================================================================

cat("Loading combined results...\n\n")

# Species level
aldex_species <- read_csv(file.path(indiv_dir, "aldex2/all_antibiotics_combined.csv"),
                          show_col_types = FALSE)
maaslin_species <- read_csv(file.path(indiv_dir, "maaslin3/all_antibiotics_combined.csv"),
                            show_col_types = FALSE)
ancombc_species <- read_csv(file.path(indiv_dir, "ancombc2/all_antibiotics_combined.csv"),
                            show_col_types = FALSE)

# Genus level
aldex_genus <- read_csv(file.path(genus_dir, "aldex2/all_antibiotics_combined.csv"),
                        show_col_types = FALSE)
maaslin_genus <- read_csv(file.path(genus_dir, "maaslin3/all_antibiotics_combined.csv"),
                          show_col_types = FALSE)
ancombc_genus <- read_csv(file.path(genus_dir, "ancombc2/all_antibiotics_combined.csv"),
                          show_col_types = FALSE)

# =============================================================================
# 2. Define Helper Functions
# =============================================================================

# Get significant taxa for an antibiotic (requiring 2+ methods)
get_robust_significant <- function(aldex_df, maaslin_df, ancombc_df,
                                   antibiotic, taxon_col = "species",
                                   threshold = 0.1) {

  # ALDEx2 significant
  sig_aldex <- aldex_df %>%
    filter(antibiotic == !!antibiotic, pval_BH < threshold) %>%
    select(taxon = all_of(taxon_col), effect_aldex = estimate)


  # MaAsLin3 significant
  sig_maaslin <- maaslin_df %>%
    filter(antibiotic == !!antibiotic, pval_BH < threshold) %>%
    select(taxon = all_of(taxon_col), effect_maaslin = coef)

  # ANCOM-BC2 significant
  if ("diff_abundant" %in% colnames(ancombc_df)) {
    sig_ancombc <- ancombc_df %>%
      filter(antibiotic == !!antibiotic, diff_abundant == TRUE) %>%
      select(taxon = all_of(taxon_col), effect_ancombc = lfc)
  } else {
    sig_ancombc <- ancombc_df %>%
      filter(antibiotic == !!antibiotic, pval_BH < threshold) %>%
      select(taxon = all_of(taxon_col), effect_ancombc = lfc)
  }

  # Combine and count methods
  all_sig <- bind_rows(
    sig_aldex %>% mutate(method = "ALDEx2"),
    sig_maaslin %>% mutate(method = "MaAsLin3"),
    sig_ancombc %>% mutate(method = "ANCOMBC2")
  )

  # Standardize effect column
  all_sig <- all_sig %>%
    mutate(effect = coalesce(effect_aldex, effect_maaslin, effect_ancombc))

  # Summarize by taxon
  summary_df <- all_sig %>%
    group_by(taxon) %>%
    summarise(
      n_methods = n(),
      methods = paste(sort(unique(method)), collapse = ", "),
      mean_effect = mean(effect, na.rm = TRUE),
      direction = ifelse(mean_effect > 0, "increased", "decreased"),
      .groups = "drop"
    )

  return(summary_df)
}

# Get any significant (1+ methods) for broader comparison
get_any_significant <- function(aldex_df, maaslin_df, ancombc_df,
                                antibiotic, taxon_col = "species",
                                threshold = 0.1) {

  sig_aldex <- aldex_df %>%
    filter(antibiotic == !!antibiotic, pval_BH < threshold) %>%
    pull(all_of(taxon_col))

  sig_maaslin <- maaslin_df %>%
    filter(antibiotic == !!antibiotic, pval_BH < threshold) %>%
    pull(all_of(taxon_col))

  if ("diff_abundant" %in% colnames(ancombc_df)) {
    sig_ancombc <- ancombc_df %>%
      filter(antibiotic == !!antibiotic, diff_abundant == TRUE) %>%
      pull(all_of(taxon_col))
  } else {
    sig_ancombc <- ancombc_df %>%
      filter(antibiotic == !!antibiotic, pval_BH < threshold) %>%
      pull(all_of(taxon_col))
  }

  return(list(
    aldex = sig_aldex,
    maaslin = sig_maaslin,
    ancombc = sig_ancombc,
    any = unique(c(sig_aldex, sig_maaslin, sig_ancombc)),
    robust = intersect(sig_aldex, sig_maaslin) %>%
      union(intersect(sig_aldex, sig_ancombc)) %>%
      union(intersect(sig_maaslin, sig_ancombc))
  ))
}

# =============================================================================
# 3. Species-Level Comparison
# =============================================================================

cat("=== SPECIES-LEVEL COMPARISON ===\n\n")

cef_species <- get_robust_significant(aldex_species, maaslin_species, ancombc_species,
                                      "Cefepime", "species")
pip_species <- get_robust_significant(aldex_species, maaslin_species, ancombc_species,
                                      "Pip_Tazo", "species")

# Get lists for Venn diagram (robust = 2+ methods)
cef_robust <- cef_species %>% filter(n_methods >= 2) %>% pull(taxon)
pip_robust <- pip_species %>% filter(n_methods >= 2) %>% pull(taxon)

cat("Cefepime: ", length(cef_robust), " robust species associations (2+ methods)\n")
cat("Pip_Tazo: ", length(pip_robust), " robust species associations (2+ methods)\n")
cat("Overlap:  ", length(intersect(cef_robust, pip_robust)), " species\n")
cat("Cefepime only: ", length(setdiff(cef_robust, pip_robust)), " species\n")
cat("Pip_Tazo only: ", length(setdiff(pip_robust, cef_robust)), " species\n\n")

# Create species Venn diagram
species_venn <- list(
  Cefepime = cef_robust,
  `Pip/Tazo` = pip_robust
)

venn_species <- ggVennDiagram(species_venn, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Species with Robust Associations (2+ methods)",
       subtitle = "Cefepime vs Piperacillin-Tazobactam") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(compare_dir, "venn_species_robust.pdf"), venn_species,
       width = 8, height = 6)

# Detailed comparison table
shared_species <- intersect(cef_robust, pip_robust)
if (length(shared_species) > 0) {
  shared_comparison <- cef_species %>%
    filter(taxon %in% shared_species) %>%
    select(species = taxon, cef_effect = mean_effect, cef_methods = n_methods) %>%
    left_join(
      pip_species %>%
        filter(taxon %in% shared_species) %>%
        select(species = taxon, pip_effect = mean_effect, pip_methods = n_methods),
      by = "species"
    ) %>%
    mutate(
      same_direction = sign(cef_effect) == sign(pip_effect),
      effect_diff = cef_effect - pip_effect
    ) %>%
    arrange(desc(abs(effect_diff)))

  cat("=== Shared Species (both antibiotics, 2+ methods) ===\n")
  cat("Same direction:", sum(shared_comparison$same_direction), "/",
      nrow(shared_comparison), "\n\n")
  print(as.data.frame(shared_comparison))

  write_csv(shared_comparison, file.path(compare_dir, "shared_species_comparison.csv"))
}

# =============================================================================
# 4. Genus-Level Comparison
# =============================================================================

cat("\n=== GENUS-LEVEL COMPARISON ===\n\n")

cef_genus <- get_robust_significant(aldex_genus, maaslin_genus, ancombc_genus,
                                    "Cefepime", "genus")
pip_genus <- get_robust_significant(aldex_genus, maaslin_genus, ancombc_genus,
                                    "Pip_Tazo", "genus")

cef_robust_g <- cef_genus %>% filter(n_methods >= 2) %>% pull(taxon)
pip_robust_g <- pip_genus %>% filter(n_methods >= 2) %>% pull(taxon)

cat("Cefepime: ", length(cef_robust_g), " robust genus associations\n")
cat("Pip_Tazo: ", length(pip_robust_g), " robust genus associations\n")
cat("Overlap:  ", length(intersect(cef_robust_g, pip_robust_g)), " genera\n\n")

# Create genus Venn diagram
genus_venn <- list(
  Cefepime = cef_robust_g,
  `Pip/Tazo` = pip_robust_g
)

venn_genus <- ggVennDiagram(genus_venn, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  labs(title = "Genera with Robust Associations (2+ methods)",
       subtitle = "Cefepime vs Piperacillin-Tazobactam") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(compare_dir, "venn_genus_robust.pdf"), venn_genus,
       width = 8, height = 6)

# Genus comparison table
shared_genera <- intersect(cef_robust_g, pip_robust_g)
if (length(shared_genera) > 0) {
  genus_comparison <- cef_genus %>%
    filter(taxon %in% shared_genera) %>%
    select(genus = taxon, cef_effect = mean_effect, cef_methods = n_methods) %>%
    left_join(
      pip_genus %>%
        filter(taxon %in% shared_genera) %>%
        select(genus = taxon, pip_effect = mean_effect, pip_methods = n_methods),
      by = "genus"
    ) %>%
    mutate(
      same_direction = sign(cef_effect) == sign(pip_effect),
      effect_diff = cef_effect - pip_effect
    ) %>%
    arrange(desc(abs(effect_diff)))

  cat("=== Shared Genera ===\n")
  print(as.data.frame(genus_comparison))

  write_csv(genus_comparison, file.path(compare_dir, "shared_genera_comparison.csv"))
}

# =============================================================================
# 5. Enterococcus Comparison
# =============================================================================

cat("\n=== ENTEROCOCCUS COMPARISON ===\n\n")

enterococcus_pattern <- "Enterococcus"

cef_entero <- cef_species %>%
  filter(grepl(enterococcus_pattern, taxon, ignore.case = TRUE)) %>%
  mutate(antibiotic = "Cefepime")

pip_entero <- pip_species %>%
  filter(grepl(enterococcus_pattern, taxon, ignore.case = TRUE)) %>%
  mutate(antibiotic = "Pip_Tazo")

cat("Cefepime Enterococcus associations:", nrow(cef_entero), "\n")
cat("  Robust (2+ methods):", sum(cef_entero$n_methods >= 2), "\n")
cat("Pip_Tazo Enterococcus associations:", nrow(pip_entero), "\n")
cat("  Robust (2+ methods):", sum(pip_entero$n_methods >= 2), "\n\n")

# Create Enterococcus Venn
cef_entero_robust <- cef_entero %>% filter(n_methods >= 2) %>% pull(taxon)
pip_entero_robust <- pip_entero %>% filter(n_methods >= 2) %>% pull(taxon)

if (length(cef_entero_robust) > 0 || length(pip_entero_robust) > 0) {
  entero_venn <- list(
    Cefepime = cef_entero_robust,
    `Pip/Tazo` = pip_entero_robust
  )

  venn_entero <- ggVennDiagram(entero_venn, label = "count", label_alpha = 0) +
    scale_fill_gradient(low = "white", high = "firebrick") +
    labs(title = "Enterococcus Species with Robust Associations",
         subtitle = "Cefepime vs Piperacillin-Tazobactam") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  ggsave(file.path(compare_dir, "venn_enterococcus.pdf"), venn_entero,
         width = 8, height = 6)
}

# Combined Enterococcus table
entero_combined <- bind_rows(cef_entero, pip_entero) %>%
  select(species = taxon, antibiotic, n_methods, mean_effect, direction) %>%
  pivot_wider(
    names_from = antibiotic,
    values_from = c(n_methods, mean_effect, direction),
    names_sep = "_"
  ) %>%
  mutate(
    shared = !is.na(n_methods_Cefepime) & !is.na(n_methods_Pip_Tazo)
  ) %>%
  arrange(desc(shared), species)

cat("Enterococcus species comparison:\n")
print(as.data.frame(entero_combined %>% head(20)))

write_csv(entero_combined, file.path(compare_dir, "enterococcus_comparison.csv"))

# =============================================================================
# 6. Enterobacteriaceae Comparison
# =============================================================================

cat("\n=== ENTEROBACTERIACEAE COMPARISON ===\n\n")

enterobact_genera <- c("Escherichia", "Klebsiella", "Enterobacter", "Citrobacter",
                       "Serratia", "Proteus", "Salmonella", "Shigella", "Morganella",
                       "Providencia", "Hafnia", "Edwardsiella", "Yersinia", "Cronobacter",
                       "Pantoea", "Raoultella", "Kosakonia", "Leclercia", "Pluralibacter")
enterobact_pattern <- paste(enterobact_genera, collapse = "|")

cef_enterobact <- cef_species %>%
  filter(grepl(enterobact_pattern, taxon, ignore.case = TRUE)) %>%
  mutate(antibiotic = "Cefepime")

pip_enterobact <- pip_species %>%
  filter(grepl(enterobact_pattern, taxon, ignore.case = TRUE)) %>%
  mutate(antibiotic = "Pip_Tazo")

cat("Cefepime Enterobacteriaceae associations:", nrow(cef_enterobact), "\n")
cat("  Robust (2+ methods):", sum(cef_enterobact$n_methods >= 2), "\n")
cat("  Increased:", sum(cef_enterobact$direction == "increased"), "\n")
cat("  Decreased:", sum(cef_enterobact$direction == "decreased"), "\n")
cat("Pip_Tazo Enterobacteriaceae associations:", nrow(pip_enterobact), "\n")
cat("  Robust (2+ methods):", sum(pip_enterobact$n_methods >= 2), "\n")
cat("  Increased:", sum(pip_enterobact$direction == "increased"), "\n")
cat("  Decreased:", sum(pip_enterobact$direction == "decreased"), "\n\n")

# Venn for Enterobacteriaceae
cef_enterobact_robust <- cef_enterobact %>% filter(n_methods >= 2) %>% pull(taxon)
pip_enterobact_robust <- pip_enterobact %>% filter(n_methods >= 2) %>% pull(taxon)

if (length(cef_enterobact_robust) > 0 || length(pip_enterobact_robust) > 0) {
  enterobact_venn <- list(
    Cefepime = cef_enterobact_robust,
    `Pip/Tazo` = pip_enterobact_robust
  )

  venn_enterobact <- ggVennDiagram(enterobact_venn, label = "count", label_alpha = 0) +
    scale_fill_gradient(low = "white", high = "purple") +
    labs(title = "Enterobacteriaceae with Robust Associations",
         subtitle = "Cefepime vs Piperacillin-Tazobactam") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  ggsave(file.path(compare_dir, "venn_enterobacteriaceae.pdf"), venn_enterobact,
         width = 8, height = 6)
}

# Combined table
enterobact_combined <- bind_rows(cef_enterobact, pip_enterobact) %>%
  select(species = taxon, antibiotic, n_methods, mean_effect, direction) %>%
  pivot_wider(
    names_from = antibiotic,
    values_from = c(n_methods, mean_effect, direction),
    names_sep = "_"
  ) %>%
  mutate(shared = !is.na(n_methods_Cefepime) & !is.na(n_methods_Pip_Tazo)) %>%
  arrange(desc(shared), species)

cat("Enterobacteriaceae comparison (top 20):\n")
print(as.data.frame(enterobact_combined %>% head(20)))

write_csv(enterobact_combined, file.path(compare_dir, "enterobacteriaceae_comparison.csv"))

# =============================================================================
# 7. Anaerobe Comparison
# =============================================================================

cat("\n=== ANAEROBE COMPARISON ===\n\n")

anaerobe_genera <- c("Bacteroides", "Prevotella", "Porphyromonas", "Fusobacterium",
                     "Clostridium", "Clostridioides", "Peptostreptococcus", "Peptoniphilus",
                     "Finegoldia", "Veillonella", "Parabacteroides", "Alistipes",
                     "Bilophila", "Desulfovibrio", "Sutterella", "Dialister",
                     "Megasphaera", "Acidaminococcus", "Anaerococcus", "Parvimonas",
                     "Eggerthella", "Slackia", "Atopobium", "Collinsella",
                     "Blautia", "Coprococcus", "Dorea", "Roseburia", "Ruminococcus",
                     "Faecalibacterium", "Eubacterium", "Lachnospira", "Butyrivibrio",
                     "Oscillibacter", "Subdoligranulum", "Megamonas", "Phascolarctobacterium")
anaerobe_pattern <- paste(anaerobe_genera, collapse = "|")

cef_anaerobe <- cef_species %>%
  filter(grepl(anaerobe_pattern, taxon, ignore.case = TRUE)) %>%
  mutate(antibiotic = "Cefepime")

pip_anaerobe <- pip_species %>%
  filter(grepl(anaerobe_pattern, taxon, ignore.case = TRUE)) %>%
  mutate(antibiotic = "Pip_Tazo")

cat("Cefepime Anaerobe associations:", nrow(cef_anaerobe), "\n")
cat("  Robust (2+ methods):", sum(cef_anaerobe$n_methods >= 2), "\n")
cat("  Increased:", sum(cef_anaerobe$direction == "increased"), "\n")
cat("  Decreased:", sum(cef_anaerobe$direction == "decreased"), "\n")
cat("Pip_Tazo Anaerobe associations:", nrow(pip_anaerobe), "\n")
cat("  Robust (2+ methods):", sum(pip_anaerobe$n_methods >= 2), "\n")
cat("  Increased:", sum(pip_anaerobe$direction == "increased"), "\n")
cat("  Decreased:", sum(pip_anaerobe$direction == "decreased"), "\n\n")

# Venn for anaerobes
cef_anaerobe_robust <- cef_anaerobe %>% filter(n_methods >= 2) %>% pull(taxon)
pip_anaerobe_robust <- pip_anaerobe %>% filter(n_methods >= 2) %>% pull(taxon)

anaerobe_venn <- list(
  Cefepime = cef_anaerobe_robust,
  `Pip/Tazo` = pip_anaerobe_robust
)

venn_anaerobe <- ggVennDiagram(anaerobe_venn, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "darkorange") +
  labs(title = "Anaerobes with Robust Associations",
       subtitle = "Cefepime vs Piperacillin-Tazobactam") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(compare_dir, "venn_anaerobes.pdf"), venn_anaerobe,
       width = 8, height = 6)

# Shared anaerobes with direction comparison
shared_anaerobes <- intersect(cef_anaerobe_robust, pip_anaerobe_robust)
if (length(shared_anaerobes) > 0) {
  anaerobe_shared <- cef_anaerobe %>%
    filter(taxon %in% shared_anaerobes) %>%
    select(species = taxon, cef_effect = mean_effect, cef_dir = direction) %>%
    left_join(
      pip_anaerobe %>%
        filter(taxon %in% shared_anaerobes) %>%
        select(species = taxon, pip_effect = mean_effect, pip_dir = direction),
      by = "species"
    ) %>%
    mutate(same_direction = cef_dir == pip_dir) %>%
    arrange(desc(same_direction), species)

  cat("Shared anaerobes (both antibiotics robust):\n")
  cat("Same direction:", sum(anaerobe_shared$same_direction), "/",
      nrow(anaerobe_shared), "\n\n")
  print(as.data.frame(anaerobe_shared))
}

# Full anaerobe comparison
anaerobe_combined <- bind_rows(cef_anaerobe, pip_anaerobe) %>%
  select(species = taxon, antibiotic, n_methods, mean_effect, direction) %>%
  pivot_wider(
    names_from = antibiotic,
    values_from = c(n_methods, mean_effect, direction),
    names_sep = "_"
  ) %>%
  mutate(shared = !is.na(n_methods_Cefepime) & !is.na(n_methods_Pip_Tazo)) %>%
  arrange(desc(shared), species)

write_csv(anaerobe_combined, file.path(compare_dir, "anaerobe_comparison.csv"))

# =============================================================================
# 8. Summary Statistics
# =============================================================================

cat("\n=============================================================================\n")
cat("SUMMARY: Cefepime vs Piperacillin-Tazobactam\n")
cat("=============================================================================\n\n")

summary_table <- data.frame(
  Category = c("All Species (robust)", "All Genera (robust)",
               "Enterococcus spp.", "Enterobacteriaceae", "Anaerobes"),
  Cefepime = c(length(cef_robust), length(cef_robust_g),
               sum(cef_entero$n_methods >= 2),
               sum(cef_enterobact$n_methods >= 2),
               sum(cef_anaerobe$n_methods >= 2)),
  Pip_Tazo = c(length(pip_robust), length(pip_robust_g),
               sum(pip_entero$n_methods >= 2),
               sum(pip_enterobact$n_methods >= 2),
               sum(pip_anaerobe$n_methods >= 2)),
  Overlap = c(length(intersect(cef_robust, pip_robust)),
              length(intersect(cef_robust_g, pip_robust_g)),
              length(intersect(cef_entero_robust, pip_entero_robust)),
              length(intersect(cef_enterobact_robust, pip_enterobact_robust)),
              length(intersect(cef_anaerobe_robust, pip_anaerobe_robust))),
  Cef_Only = c(length(setdiff(cef_robust, pip_robust)),
               length(setdiff(cef_robust_g, pip_robust_g)),
               length(setdiff(cef_entero_robust, pip_entero_robust)),
               length(setdiff(cef_enterobact_robust, pip_enterobact_robust)),
               length(setdiff(cef_anaerobe_robust, pip_anaerobe_robust))),
  Pip_Only = c(length(setdiff(pip_robust, cef_robust)),
               length(setdiff(pip_robust_g, cef_robust_g)),
               length(setdiff(pip_entero_robust, cef_entero_robust)),
               length(setdiff(pip_enterobact_robust, cef_enterobact_robust)),
               length(setdiff(pip_anaerobe_robust, cef_anaerobe_robust)))
)

print(summary_table)
write_csv(summary_table, file.path(compare_dir, "summary_comparison.csv"))

# =============================================================================
# 9. Create Combined Multi-Panel Figure
# =============================================================================

cat("\nCreating combined figure...\n")

# Combine Venns into one figure
library(patchwork)

combined_venn <- (venn_species + venn_genus) / (venn_entero + venn_anaerobe) +
  plot_annotation(
    title = "Cefepime vs Piperacillin-Tazobactam: Differential Abundance Comparison",
    subtitle = "Robust associations (2+ statistical methods)",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5)
    )
  )

ggsave(file.path(compare_dir, "combined_venn_diagrams.pdf"), combined_venn,
       width = 12, height = 10)

# =============================================================================
# 10. Effect Size Comparison Plot
# =============================================================================

cat("Creating effect size comparison plots...\n")

# For shared species, compare effect sizes
shared_all <- intersect(
  cef_species %>% filter(n_methods >= 1) %>% pull(taxon),
  pip_species %>% filter(n_methods >= 1) %>% pull(taxon)
)

if (length(shared_all) > 10) {
  effect_comparison <- cef_species %>%
    filter(taxon %in% shared_all) %>%
    select(taxon, cef_effect = mean_effect, cef_n = n_methods) %>%
    left_join(
      pip_species %>%
        filter(taxon %in% shared_all) %>%
        select(taxon, pip_effect = mean_effect, pip_n = n_methods),
      by = "taxon"
    ) %>%
    mutate(
      robust_both = cef_n >= 2 & pip_n >= 2,
      robust_either = cef_n >= 2 | pip_n >= 2
    )

  # Scatter plot of effect sizes
  effect_scatter <- ggplot(effect_comparison, aes(x = cef_effect, y = pip_effect)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue", alpha = 0.5) +
    geom_point(aes(color = robust_both, size = robust_either), alpha = 0.6) +
    scale_color_manual(values = c("gray50", "red"),
                       labels = c("One or neither robust", "Both robust")) +
    scale_size_manual(values = c(1.5, 3), guide = "none") +
    labs(
      x = "Cefepime Effect (log fold change)",
      y = "Pip/Tazo Effect (log fold change)",
      title = "Effect Size Correlation: Cefepime vs Pip/Tazo",
      subtitle = paste("Species significant in both (n =", length(shared_all), ")"),
      color = "Robustness"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    coord_fixed()

  ggsave(file.path(compare_dir, "effect_size_correlation.pdf"), effect_scatter,
         width = 8, height = 8)

  # Calculate correlation
  cor_test <- cor.test(effect_comparison$cef_effect, effect_comparison$pip_effect)
  cat("\nEffect size correlation (all shared species):\n")
  cat("  Pearson r =", round(cor_test$estimate, 3), "\n")
  cat("  p-value =", format(cor_test$p.value, digits = 3), "\n")

  # For robust only
  robust_both <- effect_comparison %>% filter(robust_both)
  if (nrow(robust_both) > 3) {
    cor_robust <- cor.test(robust_both$cef_effect, robust_both$pip_effect)
    cat("\nEffect size correlation (both robust):\n")
    cat("  Pearson r =", round(cor_robust$estimate, 3), "\n")
    cat("  p-value =", format(cor_robust$p.value, digits = 3), "\n")
  }

  write_csv(effect_comparison, file.path(compare_dir, "effect_size_comparison.csv"))
}

cat("\n=============================================================================\n")
cat("Analysis complete! Results saved to:", compare_dir, "\n")
cat("=============================================================================\n")
