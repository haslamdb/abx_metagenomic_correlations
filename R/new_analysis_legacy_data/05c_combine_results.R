#!/usr/bin/env Rscript
# =============================================================================
# 05c_combine_results.R
# Combine individual antibiotic results from ALDEx2, MaAsLin3, and ANCOM-BC2
# Run this after all individual analyses are complete
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(purrr)

select <- dplyr::select
filter <- dplyr::filter

project_dir <- "/home/david/projects/abx_metagenomic_correlations"
results_dir <- file.path(project_dir, "results/new_analysis_legacy_data")

# =============================================================================
# Configuration: Which antibiotics to include in combined results
# =============================================================================

# Current analysis antibiotics (with Vancomycin IV/PO split)
# Note: Vancomycin_PO may be excluded due to insufficient samples
target_antibiotics <- c("Pip_Tazo", "TMP_SMX", "Vancomycin_IV", "Vancomycin_PO", "Cefepime",
                        "Meropenem", "Ciprofloxacin", "Metronidazole",
                        "Ceftriaxone", "Clindamycin")

# =============================================================================
# Function to combine results from a directory
# =============================================================================

combine_method_results <- function(method_dir, method_name, file_prefix,
                                   target_abx, taxon_col = "species") {
  cat("\n=== Combining", method_name, "results ===\n")

  all_results <- list()

  for (abx in target_abx) {
    # Try common file naming patterns
    file_patterns <- c(
      paste0(file_prefix, "_", abx, ".csv"),
      paste0(file_prefix, "_", gsub("/", "_", gsub("-", "_", abx)), ".csv")
    )

    file_found <- FALSE
    for (pattern in file_patterns) {
      file_path <- file.path(method_dir, pattern)
      if (file.exists(file_path)) {
        df <- read_csv(file_path, show_col_types = FALSE)

        # Ensure antibiotic column is correct
        if (!"antibiotic" %in% colnames(df)) {
          df$antibiotic <- abx
        }

        all_results[[abx]] <- df
        cat("  Loaded:", abx, "(", nrow(df), "rows)\n")
        file_found <- TRUE
        break
      }
    }

    if (!file_found) {
      cat("  Missing:", abx, "\n")
    }
  }

  if (length(all_results) > 0) {
    combined <- bind_rows(all_results)
    output_file <- file.path(method_dir, "all_antibiotics_combined.csv")
    write_csv(combined, output_file)
    cat("  Saved combined file:", nrow(combined), "total rows\n")
    return(combined)
  } else {
    cat("  No results found!\n")
    return(NULL)
  }
}

# =============================================================================
# Combine Species-Level Results
# =============================================================================

cat("\n" , strrep("=", 60), "\n")
cat("Combining SPECIES-LEVEL Results\n")
cat(strrep("=", 60), "\n")

indiv_dir <- file.path(results_dir, "individual_antibiotics_v2_with_covariates")

# ALDEx2
combined_aldex <- combine_method_results(
  file.path(indiv_dir, "aldex2"),
  "ALDEx2", "aldex2", target_antibiotics, "species"
)

# MaAsLin3
combined_maaslin <- combine_method_results(
  file.path(indiv_dir, "maaslin3"),
  "MaAsLin3", "maaslin3", target_antibiotics, "species"
)

# ANCOM-BC2
combined_ancombc <- combine_method_results(
  file.path(indiv_dir, "ancombc2"),
  "ANCOM-BC2", "ancombc2", target_antibiotics, "species"
)

# =============================================================================
# Concordance Analysis - Species Level
# =============================================================================

cat("\n=== Concordance Analysis (Species Level) ===\n\n")

get_significant <- function(df, method_name, pval_col = "pval_BH",
                            effect_col = NULL, taxon_col = "species", threshold = 0.1) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  sig <- df %>%
    filter(.data[[pval_col]] < threshold)

  if (nrow(sig) == 0) return(NULL)

  # Standardize taxon column name
  if (taxon_col != "taxon" && taxon_col %in% colnames(sig)) {
    sig$taxon <- sig[[taxon_col]]
  }

  sig <- sig %>%
    select(antibiotic, taxon, any_of(c(effect_col, "lfc", "estimate", "coef")))

  sig$method <- method_name

  # Standardize effect column name
  effect_cols <- intersect(c("lfc", "estimate", "coef", effect_col), colnames(sig))
  if (length(effect_cols) > 0) {
    sig$effect <- sig[[effect_cols[1]]]
  } else {
    sig$effect <- NA_real_
  }

  return(sig)
}

# Get significant from each method
sig_aldex <- get_significant(combined_aldex, "ALDEx2", "pval_BH", "estimate", "species")
sig_maaslin <- get_significant(combined_maaslin, "MaAsLin3", "pval_BH", "coef", "species")
sig_ancombc <- get_significant(combined_ancombc, "ANCOMBC2", "pval_BH", "lfc", "species")

# Combine all significant
all_significant <- bind_rows(
  if (!is.null(sig_aldex)) select(sig_aldex, antibiotic, taxon, effect, method) else NULL,
  if (!is.null(sig_maaslin)) select(sig_maaslin, antibiotic, taxon, effect, method) else NULL,
  if (!is.null(sig_ancombc)) select(sig_ancombc, antibiotic, taxon, effect, method) else NULL
)

if (nrow(all_significant) > 0) {
  concordance <- all_significant %>%
    group_by(antibiotic, taxon) %>%
    summarise(
      n_methods = n(),
      methods = paste(sort(unique(method)), collapse = ", "),
      mean_effect = mean(effect, na.rm = TRUE),
      direction = ifelse(mean_effect > 0, "increased", "decreased"),
      .groups = "drop"
    ) %>%
    arrange(desc(n_methods), antibiotic, taxon)

  robust_associations <- concordance %>%
    filter(n_methods >= 2)

  cat("Associations found by 2+ methods:", nrow(robust_associations), "\n")
  cat("Associations found by all 3 methods:", sum(concordance$n_methods == 3), "\n\n")

  # Rename taxon back to species for output
  concordance <- concordance %>% rename(species = taxon)
  robust_associations <- robust_associations %>% rename(species = taxon)

  write_csv(concordance, file.path(indiv_dir, "method_concordance.csv"))
  write_csv(robust_associations, file.path(indiv_dir, "robust_associations.csv"))

  cat("Saved: method_concordance.csv, robust_associations.csv\n")
}

# =============================================================================
# Summary by Antibiotic - Species Level
# =============================================================================

cat("\n=== Summary by Antibiotic (Species Level) ===\n\n")

summary_by_abx <- data.frame(
  antibiotic = character(),
  aldex2_sig = integer(),
  maaslin3_sig = integer(),
  ancombc2_sig = integer(),
  robust_sig = integer(),
  stringsAsFactors = FALSE
)

for (abx in target_antibiotics) {
  aldex_n <- if (!is.null(combined_aldex)) {
    sum(combined_aldex$antibiotic == abx & combined_aldex$pval_BH < 0.1, na.rm = TRUE)
  } else 0

  maaslin_n <- if (!is.null(combined_maaslin)) {
    sum(combined_maaslin$antibiotic == abx & combined_maaslin$pval_BH < 0.1, na.rm = TRUE)
  } else 0

  ancombc_n <- if (!is.null(combined_ancombc)) {
    # ANCOM-BC2 uses diff_abundant column if available
    if ("diff_abundant" %in% colnames(combined_ancombc)) {
      sum(combined_ancombc$antibiotic == abx & combined_ancombc$diff_abundant, na.rm = TRUE)
    } else {
      sum(combined_ancombc$antibiotic == abx & combined_ancombc$pval_BH < 0.1, na.rm = TRUE)
    }
  } else 0

  robust_n <- if (exists("robust_associations") && nrow(robust_associations) > 0) {
    sum(robust_associations$species %in% combined_aldex$species[combined_aldex$antibiotic == abx] |
        robust_associations$species %in% combined_maaslin$species[combined_maaslin$antibiotic == abx] |
        robust_associations$species %in% combined_ancombc$species[combined_ancombc$antibiotic == abx])
  } else 0

  # Simpler robust count
  robust_n <- if (exists("robust_associations") && nrow(robust_associations) > 0) {
    sum(grepl(abx, robust_associations$species) |
        robust_associations$species %in% (all_significant %>% filter(antibiotic == abx) %>% pull(taxon)))
  } else 0

  # Actually just count from robust_associations directly
  if (exists("robust_associations") && "antibiotic" %in% colnames(robust_associations)) {
    robust_n <- sum(robust_associations$antibiotic == abx, na.rm = TRUE)
  } else {
    robust_n <- 0
  }

  summary_by_abx <- rbind(summary_by_abx, data.frame(
    antibiotic = abx,
    aldex2_sig = aldex_n,
    maaslin3_sig = maaslin_n,
    ancombc2_sig = ancombc_n,
    robust_sig = robust_n,
    stringsAsFactors = FALSE
  ))
}

print(summary_by_abx)
write_csv(summary_by_abx, file.path(indiv_dir, "summary_by_antibiotic.csv"))

# =============================================================================
# Enterococcus Focus
# =============================================================================

cat("\n=== Enterococcus Results (Species Level) ===\n\n")

enterococcus_all <- all_significant %>%
  filter(grepl("Enterococcus", taxon, ignore.case = TRUE))

if (nrow(enterococcus_all) > 0) {
  enterococcus_summary <- enterococcus_all %>%
    group_by(antibiotic, taxon) %>%
    summarise(
      n_methods = n(),
      methods = paste(sort(unique(method)), collapse = ", "),
      mean_effect = mean(effect, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(species = taxon) %>%
    arrange(desc(n_methods), antibiotic)

  print(as.data.frame(enterococcus_summary))
  write_csv(enterococcus_summary, file.path(indiv_dir, "enterococcus_by_antibiotic.csv"))
} else {
  cat("No significant Enterococcus associations found.\n")
}

# =============================================================================
# Enterobacteriaceae Focus (Species Level)
# =============================================================================

cat("\n=== Enterobacteriaceae Results (Species Level) ===\n\n")

# Enterobacteriaceae family genera
enterobact_genera <- c("Escherichia", "Klebsiella", "Enterobacter", "Citrobacter",
                       "Serratia", "Proteus", "Salmonella", "Shigella", "Morganella",
                       "Providencia", "Hafnia", "Edwardsiella", "Yersinia", "Cronobacter",
                       "Pantoea", "Raoultella", "Kosakonia", "Leclercia", "Pluralibacter")

enterobact_pattern <- paste(enterobact_genera, collapse = "|")

enterobact_all <- all_significant %>%
  filter(grepl(enterobact_pattern, taxon, ignore.case = TRUE))

if (nrow(enterobact_all) > 0) {
  enterobact_summary <- enterobact_all %>%
    group_by(antibiotic, taxon) %>%
    summarise(
      n_methods = n(),
      methods = paste(sort(unique(method)), collapse = ", "),
      mean_effect = mean(effect, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(species = taxon) %>%
    arrange(desc(n_methods), antibiotic)

  cat("Found", nrow(enterobact_summary), "significant Enterobacteriaceae associations\n")
  cat("By 3 methods:", sum(enterobact_summary$n_methods == 3), "\n")
  cat("By 2 methods:", sum(enterobact_summary$n_methods == 2), "\n\n")

  # Show top results
  top_enterobact <- enterobact_summary %>% filter(n_methods >= 2)
  if (nrow(top_enterobact) > 0) {
    cat("Robust associations (2+ methods):\n")
    print(as.data.frame(head(top_enterobact, 30)))
  }

  write_csv(enterobact_summary, file.path(indiv_dir, "enterobacteriaceae_by_antibiotic.csv"))
} else {
  cat("No significant Enterobacteriaceae associations found.\n")
}

# =============================================================================
# Anaerobe Focus (Species Level)
# =============================================================================

cat("\n=== Anaerobe Results (Species Level) ===\n\n")

# Common anaerobic genera (obligate and facultative anaerobes of clinical significance)
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

anaerobe_all <- all_significant %>%
  filter(grepl(anaerobe_pattern, taxon, ignore.case = TRUE))

if (nrow(anaerobe_all) > 0) {
  anaerobe_summary <- anaerobe_all %>%
    group_by(antibiotic, taxon) %>%
    summarise(
      n_methods = n(),
      methods = paste(sort(unique(method)), collapse = ", "),
      mean_effect = mean(effect, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(species = taxon) %>%
    arrange(desc(n_methods), antibiotic)

  cat("Found", nrow(anaerobe_summary), "significant anaerobe associations\n")
  cat("By 3 methods:", sum(anaerobe_summary$n_methods == 3), "\n")
  cat("By 2 methods:", sum(anaerobe_summary$n_methods == 2), "\n\n")

  # Summary by direction
  direction_summary <- anaerobe_summary %>%
    mutate(direction = ifelse(mean_effect > 0, "increased", "decreased")) %>%
    group_by(antibiotic, direction) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = direction, values_from = n, values_fill = 0)

  cat("Anaerobe changes by antibiotic (increased vs decreased):\n")
  print(as.data.frame(direction_summary))

  # Show top results
  top_anaerobe <- anaerobe_summary %>% filter(n_methods >= 2)
  if (nrow(top_anaerobe) > 0) {
    cat("\nRobust associations (2+ methods):\n")
    print(as.data.frame(head(top_anaerobe, 30)))
  }

  write_csv(anaerobe_summary, file.path(indiv_dir, "anaerobe_by_antibiotic.csv"))
} else {
  cat("No significant anaerobe associations found.\n")
}

# =============================================================================
# Combine Genus-Level Results
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("Combining GENUS-LEVEL Results\n")
cat(strrep("=", 60), "\n")

genus_dir <- file.path(results_dir, "individual_antibiotics_genus_level")

if (dir.exists(genus_dir)) {
  # ALDEx2
  combined_aldex_genus <- combine_method_results(
    file.path(genus_dir, "aldex2"),
    "ALDEx2", "aldex2", target_antibiotics, "genus"
  )

  # MaAsLin3
  combined_maaslin_genus <- combine_method_results(
    file.path(genus_dir, "maaslin3"),
    "MaAsLin3", "maaslin3", target_antibiotics, "genus"
  )

  # ANCOM-BC2
  combined_ancombc_genus <- combine_method_results(
    file.path(genus_dir, "ancombc2"),
    "ANCOM-BC2", "ancombc2", target_antibiotics, "genus"
  )

  # Concordance for genus level
  cat("\n=== Concordance Analysis (Genus Level) ===\n\n")

  sig_aldex_g <- get_significant(combined_aldex_genus, "ALDEx2", "pval_BH", "estimate", "genus")
  sig_maaslin_g <- get_significant(combined_maaslin_genus, "MaAsLin3", "pval_BH", "coef", "genus")
  sig_ancombc_g <- get_significant(combined_ancombc_genus, "ANCOMBC2", "pval_BH", "lfc", "genus")

  all_sig_genus <- bind_rows(
    if (!is.null(sig_aldex_g)) select(sig_aldex_g, antibiotic, taxon, effect, method) else NULL,
    if (!is.null(sig_maaslin_g)) select(sig_maaslin_g, antibiotic, taxon, effect, method) else NULL,
    if (!is.null(sig_ancombc_g)) select(sig_ancombc_g, antibiotic, taxon, effect, method) else NULL
  )

  if (nrow(all_sig_genus) > 0) {
    concordance_genus <- all_sig_genus %>%
      group_by(antibiotic, taxon) %>%
      summarise(
        n_methods = n(),
        methods = paste(sort(unique(method)), collapse = ", "),
        mean_effect = mean(effect, na.rm = TRUE),
        direction = ifelse(mean_effect > 0, "increased", "decreased"),
        .groups = "drop"
      ) %>%
      arrange(desc(n_methods), antibiotic, taxon) %>%
      rename(genus = taxon)

    robust_genus <- concordance_genus %>%
      filter(n_methods >= 2)

    cat("Associations found by 2+ methods:", nrow(robust_genus), "\n")
    cat("Associations found by all 3 methods:", sum(concordance_genus$n_methods == 3), "\n\n")

    write_csv(concordance_genus, file.path(genus_dir, "method_concordance.csv"))
    write_csv(robust_genus, file.path(genus_dir, "robust_associations.csv"))

    # Summary by antibiotic - genus level
    summary_genus <- data.frame(
      antibiotic = target_antibiotics,
      aldex2_sig = sapply(target_antibiotics, function(a) {
        if (!is.null(combined_aldex_genus)) sum(combined_aldex_genus$antibiotic == a & combined_aldex_genus$pval_BH < 0.1, na.rm = TRUE) else 0
      }),
      maaslin3_sig = sapply(target_antibiotics, function(a) {
        if (!is.null(combined_maaslin_genus)) sum(combined_maaslin_genus$antibiotic == a & combined_maaslin_genus$pval_BH < 0.1, na.rm = TRUE) else 0
      }),
      ancombc2_sig = sapply(target_antibiotics, function(a) {
        if (!is.null(combined_ancombc_genus)) {
          if ("diff_abundant" %in% colnames(combined_ancombc_genus)) {
            sum(combined_ancombc_genus$antibiotic == a & combined_ancombc_genus$diff_abundant, na.rm = TRUE)
          } else {
            sum(combined_ancombc_genus$antibiotic == a & combined_ancombc_genus$pval_BH < 0.1, na.rm = TRUE)
          }
        } else 0
      }),
      robust_sig = sapply(target_antibiotics, function(a) sum(robust_genus$antibiotic == a))
    )

    print(summary_genus)
    write_csv(summary_genus, file.path(genus_dir, "summary_by_antibiotic.csv"))
  }
} else {
  cat("Genus-level results directory not found.\n")
}

# =============================================================================
# Final Summary
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("COMBINING COMPLETE\n")
cat(strrep("=", 60), "\n\n")

cat("Species-level results:", indiv_dir, "\n")
if (dir.exists(genus_dir)) {
  cat("Genus-level results:", genus_dir, "\n")
}

cat("\nKey output files:\n")
cat("  - */all_antibiotics_combined.csv (per method)\n")
cat("  - method_concordance.csv\n")
cat("  - robust_associations.csv\n")
cat("  - summary_by_antibiotic.csv\n")
cat("  - enterococcus_by_antibiotic.csv\n")
cat("  - enterobacteriaceae_by_antibiotic.csv\n")
cat("  - anaerobe_by_antibiotic.csv\n")
