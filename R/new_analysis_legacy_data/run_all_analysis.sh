#!/bin/bash
# Run all analysis scripts in sequence
# Note: Vancomycin is now split by route (IV vs PO) in all analyses

cd /home/david/projects/abx_metagenomic_correlations

# Activate r-env conda environment
source ~/miniforge3/etc/profile.d/conda.sh
conda activate r-env

echo "=========================================="
echo "Starting full analysis pipeline..."
echo "Note: Vancomycin split by route (IV vs PO)"
echo "=========================================="

# Script 00 - Build count matrices from Bracken
echo ""
echo ">>> Running 00_build_count_matrices_from_bracken.R..."
echo ">>> Started at: $(date)"
Rscript R/new_analysis_legacy_data/00_build_count_matrices_from_bracken.R 2>&1 | tee results/new_analysis_legacy_data/00_build_count_matrices.log
echo ">>> Completed at: $(date)"

# Script 01 - Load and prepare data (creates prepared_data.RData with Vanc IV/PO split)
echo ""
echo ">>> Running 01_load_and_prepare_data.R..."
echo ">>> Started at: $(date)"
Rscript R/new_analysis_legacy_data/01_load_and_prepare_data.R 2>&1 | tee results/new_analysis_legacy_data/01_load_and_prepare_data.log
echo ">>> Completed at: $(date)"

# Script 02 - Diversity analysis
echo ""
echo ">>> Running 02_diversity_analysis.R..."
echo ">>> Started at: $(date)"
Rscript R/new_analysis_legacy_data/02_diversity_analysis.R 2>&1 | tee results/new_analysis_legacy_data/02_diversity_analysis.log
echo ">>> Completed at: $(date)"

# Script 02b - Diversity by individual antibiotic
echo ""
echo ">>> Running 02b_diversity_individual_antibiotics.R..."
echo ">>> Started at: $(date)"
Rscript R/new_analysis_legacy_data/02b_diversity_individual_antibiotics.R 2>&1 | tee results/new_analysis_legacy_data/02b_diversity_individual_antibiotics.log
echo ">>> Completed at: $(date)"

# Script 03 - Differential abundance
echo ""
echo ">>> Running 03_differential_abundance.R..."
echo ">>> Started at: $(date)"
Rscript R/new_analysis_legacy_data/03_differential_abundance.R 2>&1 | tee results/new_analysis_legacy_data/03_differential_abundance.log
echo ">>> Completed at: $(date)"

# Script 04 - Paired sample analysis
echo ""
echo ">>> Running 04_paired_sample_analysis.R..."
echo ">>> Started at: $(date)"
Rscript R/new_analysis_legacy_data/04_paired_sample_analysis.R 2>&1 | tee results/new_analysis_legacy_data/04_paired_sample_analysis.log
echo ">>> Completed at: $(date)"

# Script 05 - Species-level differential abundance (individual antibiotics)
echo ""
echo ">>> Running 05_individual_antibiotic_effects.R..."
echo ">>> Started at: $(date)"
Rscript R/new_analysis_legacy_data/05_individual_antibiotic_effects.R 2>&1 | tee results/new_analysis_legacy_data/05_species_analysis.log
echo ">>> Completed at: $(date)"

# Script 06 - Genus-level differential abundance
echo ""
echo ">>> Running 06_genus_level_analysis.R..."
echo ">>> Started at: $(date)"
Rscript R/new_analysis_legacy_data/06_genus_level_analysis.R 2>&1 | tee results/new_analysis_legacy_data/06_genus_analysis.log
echo ">>> Completed at: $(date)"

# Script 06 - Network visualization
echo ""
echo ">>> Running 06_network_visualization.R..."
echo ">>> Started at: $(date)"
Rscript R/new_analysis_legacy_data/06_network_visualization.R 2>&1 | tee results/new_analysis_legacy_data/06_network_visualization.log
echo ">>> Completed at: $(date)"

echo ""
echo "=========================================="
echo "All analyses completed!"
echo "=========================================="
