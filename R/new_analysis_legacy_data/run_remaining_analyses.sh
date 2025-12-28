#!/bin/bash
# Run remaining analysis scripts in sequence

cd /home/david/projects/abx_metagenomic_correlations

# Activate r-env conda environment
source ~/miniforge3/etc/profile.d/conda.sh
conda activate r-env

echo "=========================================="
echo "Starting remaining analyses..."
echo "=========================================="

# Script 06 - Genus-level analysis
echo ""
echo ">>> Running 06_genus_level_analysis.R..."
echo ">>> Started at: $(date)"
Rscript R/new_analysis_legacy_data/06_genus_level_analysis.R 2>&1 | tee results/genus_analysis.log
echo ">>> Completed at: $(date)"

# Script 02 - Diversity analysis
echo ""
echo ">>> Running 02_diversity_analysis.R..."
echo ">>> Started at: $(date)"
Rscript R/new_analysis_legacy_data/02_diversity_analysis.R 2>&1 | tee results/new_analysis_legacy_data/diversity_analysis.log
echo ">>> Completed at: $(date)"

# Script 02b - Diversity by individual antibiotic
echo ""
echo ">>> Running 02b_diversity_individual_antibiotics.R..."
echo ">>> Started at: $(date)"
Rscript R/new_analysis_legacy_data/02b_diversity_individual_antibiotics.R 2>&1 | tee results/new_analysis_legacy_data/diversity_individual_antibiotics.log
echo ">>> Completed at: $(date)"

# Script 03 - Differential abundance
echo ""
echo ">>> Running 03_differential_abundance.R..."
echo ">>> Started at: $(date)"
Rscript R/new_analysis_legacy_data/03_differential_abundance.R 2>&1 | tee results/new_analysis_legacy_data/differential_abundance.log
echo ">>> Completed at: $(date)"

# Script 04 - Paired sample analysis
echo ""
echo ">>> Running 04_paired_sample_analysis.R..."
echo ">>> Started at: $(date)"
Rscript R/new_analysis_legacy_data/04_paired_sample_analysis.R 2>&1 | tee results/new_analysis_legacy_data/paired_sample_analysis.log
echo ">>> Completed at: $(date)"

echo ""
echo "=========================================="
echo "All analyses completed!"
echo "=========================================="
