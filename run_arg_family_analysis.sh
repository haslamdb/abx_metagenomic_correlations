#!/bin/bash
# =============================================================================
# run_arg_family_analysis.sh
# Run ARG family-level differential abundance analysis
#
# This script:
#   1. Aggregates ~3000 individual ARGs to ~100-200 gene families
#   2. Runs differential abundance (ALDEx2, MaAsLin3, ANCOM-BC2)
#      with patient_group as covariate
#
# Usage: ./run_arg_family_analysis.sh
# =============================================================================

set -e  # Exit on error

# Activate conda environment with maaslin3
source ~/miniforge3/etc/profile.d/conda.sh
conda activate r-env

PROJECT_DIR="/home/david/projects/abx_metagenomic_correlations"
LOG_DIR="${PROJECT_DIR}/logs"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

mkdir -p "$LOG_DIR"

echo "=============================================="
echo "ARG Family-Level Analysis Pipeline"
echo "Started: $(date)"
echo "=============================================="

# Script 1: Build family-level matrices
echo ""
echo "[1/2] Aggregating ARGs to gene families..."
echo "      (TEM, SHV, CTX-M, etc.)"
echo ""

Rscript "${PROJECT_DIR}/R/new_analysis_legacy_data/00c_build_arg_family_matrices.R" 2>&1 | tee "${LOG_DIR}/00c_build_arg_family_matrices_${TIMESTAMP}.log"

echo ""
echo "[1/2] Complete!"
echo ""

# Script 2: Differential abundance analysis at family level
echo "[2/2] Running family-level differential abundance analysis..."
echo "      (ALDEx2, MaAsLin3, ANCOM-BC2 with patient_group covariate)"
echo ""

Rscript "${PROJECT_DIR}/R/new_analysis_legacy_data/11b_arg_family_differential_abundance.R" 2>&1 | tee "${LOG_DIR}/11b_arg_family_differential_abundance_${TIMESTAMP}.log"

echo ""
echo "[2/2] Complete!"
echo ""

echo "=============================================="
echo "ARG Family-Level Analysis Pipeline Complete"
echo "Finished: $(date)"
echo "=============================================="
echo ""
echo "Results saved to:"
echo "  - results/new_analysis_legacy_data/arg_family_differential_abundance/"
echo ""
echo "Key improvements over gene-level analysis:"
echo "  - Reduced from ~3000 genes to ~100-200 families"
echo "  - Lower multiple testing burden = more power"
echo "  - Patient group included as covariate"
echo "  - More interpretable results"
echo ""
echo "Logs saved to:"
echo "  - ${LOG_DIR}/00c_build_arg_family_matrices_${TIMESTAMP}.log"
echo "  - ${LOG_DIR}/11b_arg_family_differential_abundance_${TIMESTAMP}.log"
