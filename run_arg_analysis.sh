#!/bin/bash
# =============================================================================
# run_arg_analysis.sh
# Run ARG differential abundance and paired sample analyses in sequence
# Usage: ./run_arg_analysis.sh
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
echo "ARG Analysis Pipeline"
echo "Started: $(date)"
echo "=============================================="

# Script 1: Differential abundance analysis
echo ""
echo "[1/2] Running ARG differential abundance analysis..."
echo "      (ALDEx2, MaAsLin3, ANCOM-BC2)"
echo ""

Rscript "${PROJECT_DIR}/R/new_analysis_legacy_data/11_arg_differential_abundance.R" 2>&1 | tee "${LOG_DIR}/11_arg_differential_abundance_${TIMESTAMP}.log"

echo ""
echo "[1/2] Complete!"
echo ""

# Script 2: Paired sample analysis
echo "[2/2] Running ARG paired sample analysis..."
echo "      (LMM with individual antibiotic effects)"
echo ""

Rscript "${PROJECT_DIR}/R/new_analysis_legacy_data/12_arg_paired_sample_analysis.R" 2>&1 | tee "${LOG_DIR}/12_arg_paired_sample_analysis_${TIMESTAMP}.log"

echo ""
echo "[2/2] Complete!"
echo ""

echo "=============================================="
echo "ARG Analysis Pipeline Complete"
echo "Finished: $(date)"
echo "=============================================="
echo ""
echo "Results saved to:"
echo "  - results/new_analysis_legacy_data/arg_differential_abundance/"
echo "  - results/new_analysis_legacy_data/arg_paired_analysis/"
echo ""
echo "Logs saved to:"
echo "  - ${LOG_DIR}/11_arg_differential_abundance_${TIMESTAMP}.log"
echo "  - ${LOG_DIR}/12_arg_paired_sample_analysis_${TIMESTAMP}.log"
