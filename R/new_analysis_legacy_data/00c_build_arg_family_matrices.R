#!/usr/bin/env Rscript
# =============================================================================
# 00c_build_arg_family_matrices.R
# Aggregate individual ARG variants to gene family level
#
# This reduces ~3000 individual gene variants to ~100-200 gene families,
# which improves:
#   - Statistical power (fewer tests, reduced FDR burden)
#   - Biological interpretability
#   - Robustness (less affected by annotation differences between variants)
# =============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# Set paths
project_dir <- "/home/david/projects/abx_metagenomic_correlations"
output_dir <- file.path(project_dir, "data")

cat("=============================================================================\n")
cat("Aggregating ARGs to Gene Family Level\n")
cat("=============================================================================\n\n")

# =============================================================================
# 1. Load ARG Data
# =============================================================================

load(file.path(output_dir, "arg_count_matrices.RData"))
arg_count_matrix <- arg_data$arg_count_matrix

cat("Input ARG matrix:", nrow(arg_count_matrix), "samples x",
    ncol(arg_count_matrix), "genes\n\n")

# =============================================================================
# 2. Define Gene Family Mapping Rules
# =============================================================================

# Function to map individual ARG names to gene families
map_to_family <- function(gene_name) {

  # ----- Beta-lactamases -----
  # TEM variants
  if (grepl("TEM-\\d+", gene_name)) {
    if (grepl("extended-spectrum", gene_name, ignore.case = TRUE)) {
      return("TEM_ESBL")
    } else if (grepl("inhibitor-resistant", gene_name, ignore.case = TRUE)) {
      return("TEM_inhibitor_resistant")
    } else {
      return("TEM_beta_lactamase")
    }
  }

  # SHV variants
  if (grepl("SHV-\\d+", gene_name)) {
    if (grepl("extended-spectrum", gene_name, ignore.case = TRUE)) {
      return("SHV_ESBL")
    } else if (grepl("inhibitor-resistant", gene_name, ignore.case = TRUE)) {
      return("SHV_inhibitor_resistant")
    } else {
      return("SHV_beta_lactamase")
    }
  }

  # CTX-M variants
  if (grepl("CTX-M", gene_name, ignore.case = TRUE)) {
    return("CTX_M_ESBL")
  }

  # OXA variants
  if (grepl("OXA-\\d+", gene_name)) {
    if (grepl("carbapenem", gene_name, ignore.case = TRUE)) {
      return("OXA_carbapenemase")
    } else {
      return("OXA_beta_lactamase")
    }
  }

  # KPC carbapenemases
  if (grepl("KPC", gene_name, ignore.case = TRUE)) {
    return("KPC_carbapenemase")
  }

  # NDM carbapenemases
  if (grepl("NDM", gene_name, ignore.case = TRUE)) {
    return("NDM_carbapenemase")
  }

  # VIM carbapenemases
  if (grepl("VIM-\\d+", gene_name)) {
    return("VIM_carbapenemase")
  }

  # IMP carbapenemases
  if (grepl("IMP-\\d+", gene_name)) {
    return("IMP_carbapenemase")
  }

  # CMY/AmpC beta-lactamases
  if (grepl("CMY|AmpC|MOX|DHA|ACT|MIR|FOX", gene_name, ignore.case = TRUE) &&
      grepl("beta-lactamase", gene_name, ignore.case = TRUE)) {
    return("AmpC_beta_lactamase")
  }

  # Other class C beta-lactamases
  if (grepl("class C.*beta-lactamase", gene_name, ignore.case = TRUE)) {
    return("Class_C_beta_lactamase")
  }

  # Other class A beta-lactamases
  if (grepl("class A.*beta-lactamase", gene_name, ignore.case = TRUE)) {
    return("Class_A_beta_lactamase_other")
  }

  # Metallo-beta-lactamases
  if (grepl("metallo-beta-lactamase", gene_name, ignore.case = TRUE)) {
    return("Metallo_beta_lactamase")
  }

  # ----- Aminoglycoside resistance -----
  # AAC acetyltransferases
  if (grepl("AAC\\(3\\)", gene_name)) {
    return("AAC3_acetyltransferase")
  }
  if (grepl("AAC\\(6", gene_name)) {
    return("AAC6_acetyltransferase")
  }
  if (grepl("AAC\\(2", gene_name)) {
    return("AAC2_acetyltransferase")
  }

  # APH phosphotransferases
  if (grepl("APH\\(3", gene_name)) {
    return("APH3_phosphotransferase")
  }
  if (grepl("APH\\(6", gene_name)) {
    return("APH6_phosphotransferase")
  }
  if (grepl("APH\\(2", gene_name)) {
    return("APH2_phosphotransferase")
  }
  if (grepl("APH\\(4", gene_name)) {
    return("APH4_phosphotransferase")
  }
  if (grepl("APH", gene_name) && grepl("phosphotransferase", gene_name, ignore.case = TRUE)) {
    return("APH_phosphotransferase_other")
  }

  # ANT nucleotidyltransferases
  if (grepl("ANT\\(|AadA|AadB|AadE", gene_name)) {
    return("ANT_nucleotidyltransferase")
  }

  # 16S rRNA methyltransferases (high-level aminoglycoside resistance)
  if (grepl("16S rRNA.*methyltransferase|Rmt|Arm", gene_name, ignore.case = TRUE)) {
    return("Aminoglycoside_16S_methyltransferase")
  }

  # ----- Vancomycin resistance -----
  if (grepl("VanA|VanH-A", gene_name) && !grepl("VanA-", gene_name)) {
    return("VanA_type")
  }
  if (grepl("VanB|VanH-B", gene_name)) {
    return("VanB_type")
  }
  if (grepl("VanC|VanT|VanXY-C", gene_name)) {
    return("VanC_type")
  }
  if (grepl("VanD|VanH-D", gene_name)) {
    return("VanD_type")
  }
  if (grepl("VanE|VanXY-E", gene_name)) {
    return("VanE_type")
  }
  if (grepl("VanG|VanXY-G", gene_name)) {
    return("VanG_type")
  }
  if (grepl("Van[XYRHSW]|VanZ", gene_name)) {
    return("Van_accessory")
  }

  # ----- Macrolide/Lincosamide/Streptogramin (MLS) resistance -----
  # Erm methyltransferases
  if (grepl("Erm\\(|Erm[A-Z]", gene_name)) {
    return("Erm_methyltransferase")
  }

  # Mph phosphotransferases
  if (grepl("Mph\\(|Mph[A-Z]", gene_name)) {
    return("Mph_phosphotransferase")
  }

  # Msr/Mef efflux
  if (grepl("Msr\\(|Msr[A-Z]", gene_name)) {
    return("Msr_ABC_efflux")
  }
  if (grepl("Mef\\(|Mef[A-Z]", gene_name)) {
    return("Mef_MFS_efflux")
  }

  # Ere esterases
  if (grepl("Ere\\(|Ere[A-Z]", gene_name)) {
    return("Ere_esterase")
  }

  # Lincosamide nucleotidyltransferases
  if (grepl("Lnu\\(|Lin[A-Z]", gene_name)) {
    return("Lnu_nucleotidyltransferase")
  }

  # Lsa ABC-F proteins
  if (grepl("Lsa\\(|Lsa[A-Z]", gene_name)) {
    return("Lsa_ABC_F")
  }

  # Vga ABC-F proteins
  if (grepl("Vga\\(|Vga[A-Z]", gene_name)) {
    return("Vga_ABC_F")
  }

  # ----- Tetracycline resistance -----
  # Ribosomal protection proteins
  if (grepl("Tet\\([MOSQWOB3246]|TetB\\(P\\)|tetracycline.*ribosomal protection", gene_name, ignore.case = TRUE)) {
    return("Tet_ribosomal_protection")
  }

  # Efflux pumps
  if (grepl("Tet\\([AKLCDEGHJ]|TetA\\(P\\)|tetracycline.*efflux", gene_name, ignore.case = TRUE)) {
    return("Tet_efflux")
  }

  # Enzymatic inactivation
  if (grepl("Tet\\(X|tetracycline.*monooxygenase", gene_name, ignore.case = TRUE)) {
    return("Tet_enzymatic")
  }

  # ----- Fluoroquinolone resistance -----
  # Qnr proteins
  if (grepl("Qnr[ABCDSE]|qnr", gene_name, ignore.case = TRUE)) {
    return("Qnr_quinolone_resistance")
  }

  # Aac(6')-Ib-cr (ciprofloxacin-modifying)
  if (grepl("Aac.*cr|ciprofloxacin", gene_name, ignore.case = TRUE)) {
    return("AAC6_Ib_cr_fluoroquinolone")
  }

  # ----- Sulfonamide resistance -----
  if (grepl("Sul[123]|dihydropteroate synthase.*Sul", gene_name, ignore.case = TRUE)) {
    return("Sul_sulfonamide")
  }

  # ----- Trimethoprim resistance -----
  if (grepl("Dfr[A-Z]|dihydrofolate reductase.*Dfr", gene_name, ignore.case = TRUE)) {
    return("Dfr_trimethoprim")
  }

  # ----- Chloramphenicol resistance -----
  if (grepl("Cat[ABC]|chloramphenicol acetyltransferase", gene_name, ignore.case = TRUE)) {
    return("Cat_chloramphenicol")
  }
  if (grepl("Cml|FloR|chloramphenicol.*efflux", gene_name, ignore.case = TRUE)) {
    return("Chloramphenicol_efflux")
  }

  # ----- Rifampin resistance -----
  if (grepl("Arr|rifampin", gene_name, ignore.case = TRUE)) {
    return("Arr_rifampin")
  }

  # ----- Multidrug efflux pumps -----
  # Oqx efflux
  if (grepl("Oqx[AB]", gene_name)) {
    return("Oqx_RND_efflux")
  }

  # Mex efflux
  if (grepl("Mex[A-Z]", gene_name)) {
    return("Mex_RND_efflux")
  }

  # AcrAB-TolC
  if (grepl("Acr[AB]|TolC", gene_name)) {
    return("AcrAB_RND_efflux")
  }

  # QacE (quaternary ammonium)
  if (grepl("QacE|qac", gene_name, ignore.case = TRUE)) {
    return("Qac_efflux")
  }

  # ----- ABC-F ribosomal protection -----
  if (grepl("ABC-F.*ribosomal protection|OptrA", gene_name, ignore.case = TRUE)) {
    return("ABC_F_ribosomal_protection")
  }

  # ----- Colistin resistance -----
  if (grepl("Mcr-|MCR-|phosphoethanolamine transferase", gene_name, ignore.case = TRUE)) {
    return("Mcr_colistin")
  }

  # ----- Fosfomycin resistance -----
  if (grepl("Fos[ABCX]|fosfomycin", gene_name, ignore.case = TRUE)) {
    return("Fos_fosfomycin")
  }

  # ----- Other/unclassified -----
  # Generic efflux
  if (grepl("efflux|transporter", gene_name, ignore.case = TRUE)) {
    return("Other_efflux")
  }

  # Generic beta-lactamase
  if (grepl("beta-lactamase", gene_name, ignore.case = TRUE)) {
    return("Other_beta_lactamase")
  }

  # Default: keep original name if no pattern matches
  return(gene_name)
}

# =============================================================================
# 3. Apply Mapping to All Genes
# =============================================================================

gene_names <- colnames(arg_count_matrix)
cat("Mapping", length(gene_names), "genes to families...\n")

# Create mapping dataframe
gene_mapping <- data.frame(
  gene = gene_names,
  family = sapply(gene_names, map_to_family),
  stringsAsFactors = FALSE
)

# Summary of mapping
family_counts <- table(gene_mapping$family)
n_families <- length(family_counts)
n_mapped <- sum(gene_mapping$family != gene_mapping$gene)

cat("  Genes mapped to families:", n_mapped, "\n")
cat("  Genes kept as-is:", sum(gene_mapping$family == gene_mapping$gene), "\n")
cat("  Total families:", n_families, "\n\n")

# Show top families by gene count
cat("Top 30 gene families by number of variants:\n")
top_families <- sort(family_counts, decreasing = TRUE)[1:30]
for (i in 1:30) {
  cat(sprintf("  %2d. %s: %d variants\n", i, names(top_families)[i], top_families[i]))
}

# =============================================================================
# 4. Aggregate Counts to Family Level
# =============================================================================

cat("\n\nAggregating counts to family level...\n")

# Convert count matrix to long format
arg_long <- as.data.frame(arg_count_matrix) %>%
  tibble::rownames_to_column("sample") %>%
  pivot_longer(cols = -sample, names_to = "gene", values_to = "count")

# Add family mapping
arg_long <- arg_long %>%
  left_join(gene_mapping, by = "gene")

# Aggregate by family (sum counts within each family)
family_agg <- arg_long %>%
  group_by(sample, family) %>%
  summarise(count = sum(count), .groups = "drop")

# Pivot back to wide format
family_wide <- family_agg %>%
  pivot_wider(names_from = family, values_from = count, values_fill = 0)

# Convert to matrix
arg_family_matrix <- as.matrix(family_wide[, -1])
rownames(arg_family_matrix) <- family_wide$sample

cat("Output family matrix:", nrow(arg_family_matrix), "samples x",
    ncol(arg_family_matrix), "families\n")

# =============================================================================
# 5. Also Create Coverage Matrix at Family Level
# =============================================================================

cat("\nAggregating coverage to family level...\n")

if (!is.null(arg_data$arg_coverage_matrix)) {
  arg_coverage_matrix <- arg_data$arg_coverage_matrix

  # For coverage, take the max within each family (most confident detection)
  coverage_long <- as.data.frame(arg_coverage_matrix) %>%
    tibble::rownames_to_column("sample") %>%
    pivot_longer(cols = -sample, names_to = "gene", values_to = "coverage")

  # Some genes in coverage may not be in count matrix - filter to common genes
  common_genes <- intersect(coverage_long$gene, gene_mapping$gene)
  coverage_long <- coverage_long %>%
    filter(gene %in% common_genes) %>%
    left_join(gene_mapping, by = "gene")

  # Max coverage within family
  coverage_family_agg <- coverage_long %>%
    group_by(sample, family) %>%
    summarise(coverage = max(coverage, na.rm = TRUE), .groups = "drop")

  coverage_family_wide <- coverage_family_agg %>%
    pivot_wider(names_from = family, values_from = coverage, values_fill = 0)

  arg_family_coverage_matrix <- as.matrix(coverage_family_wide[, -1])
  rownames(arg_family_coverage_matrix) <- coverage_family_wide$sample

  cat("Family coverage matrix:", nrow(arg_family_coverage_matrix), "samples x",
      ncol(arg_family_coverage_matrix), "families\n")
} else {
  arg_family_coverage_matrix <- NULL
}

# =============================================================================
# 6. Quality Checks
# =============================================================================

cat("\n=== Quality Checks ===\n\n")

# Total counts should match
original_totals <- rowSums(arg_count_matrix)
family_totals <- rowSums(arg_family_matrix)

# Match samples
common_samples <- intersect(names(original_totals), rownames(arg_family_matrix))
orig_subset <- original_totals[common_samples]
fam_subset <- family_totals[common_samples]

cat("Count preservation check:\n")
cat("  Original total counts:", sum(orig_subset), "\n")
cat("  Family total counts:", sum(fam_subset), "\n")
cat("  Match:", all.equal(sum(orig_subset), sum(fam_subset)), "\n\n")

# Family prevalence
family_prevalence <- colSums(arg_family_matrix > 0) / nrow(arg_family_matrix)
cat("Top 20 families by prevalence:\n")
top_prev <- sort(family_prevalence, decreasing = TRUE)[1:20]
for (i in 1:20) {
  cat(sprintf("  %2d. %s: %.1f%%\n", i, names(top_prev)[i], top_prev[i] * 100))
}

# =============================================================================
# 7. Save Results
# =============================================================================

cat("\n=== Saving Results ===\n")

# Update arg_data with family-level matrices
arg_family_data <- list(
  # Family-level matrices
  arg_family_matrix = arg_family_matrix,
  arg_family_coverage_matrix = arg_family_coverage_matrix,

  # Gene-to-family mapping
  gene_mapping = gene_mapping,

  # Original matrices (reference)
  arg_count_matrix = arg_data$arg_count_matrix,
  arg_coverage_matrix = arg_data$arg_coverage_matrix,
  arg_rpm_matrix = arg_data$arg_rpm_matrix,
  sample_bacterial_reads = arg_data$sample_bacterial_reads,

  # Metadata
  n_samples = nrow(arg_family_matrix),
  n_families = ncol(arg_family_matrix),
  n_original_genes = ncol(arg_data$arg_count_matrix),
  date_created = Sys.time()
)

save(arg_family_data, file = file.path(output_dir, "arg_family_matrices.RData"))
cat("Saved: data/arg_family_matrices.RData\n")

# Save CSVs
write.csv(arg_family_matrix, file.path(output_dir, "arg_family_counts.csv"))
write.csv(gene_mapping, file.path(output_dir, "arg_gene_to_family_mapping.csv"), row.names = FALSE)
cat("Saved: data/arg_family_counts.csv\n")
cat("Saved: data/arg_gene_to_family_mapping.csv\n")

cat("\n=== Done ===\n")
cat("Reduced from", ncol(arg_data$arg_count_matrix), "genes to",
    ncol(arg_family_matrix), "families\n")
cat("Reduction factor:", round(ncol(arg_data$arg_count_matrix) / ncol(arg_family_matrix), 1), "x\n")
