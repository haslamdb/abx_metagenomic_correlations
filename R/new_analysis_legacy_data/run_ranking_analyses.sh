#!/bin/bash
# =============================================================================
# run_ranking_analyses.sh
# Run species ranking analyses (during exposure + recovery) in a screen session
#
# Usage:
#   ./run_ranking_analyses.sh          # Run in foreground
#   ./run_ranking_analyses.sh screen   # Run in detached screen session
#
# To reattach to screen session:
#   screen -r ranking_analysis
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
LOG_DIR="$PROJECT_DIR/logs"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

mkdir -p "$LOG_DIR"

run_analyses() {
    echo "=============================================="
    echo "SPECIES RANKING ANALYSES"
    echo "Started: $(date)"
    echo "Project: $PROJECT_DIR"
    echo "=============================================="
    echo ""

    cd "$PROJECT_DIR"

    # Script 1: During antibiotic exposure ranking
    echo "=== Running 14_species_exposure_ranking.R ==="
    echo "Started: $(date)"
    Rscript "$SCRIPT_DIR/14_species_exposure_ranking.R" 2>&1
    EXIT1=$?
    echo "Finished: $(date), Exit code: $EXIT1"
    echo ""

    # Script 2: Recovery ranking
    echo "=== Running 13_species_recovery_ranking.R ==="
    echo "Started: $(date)"
    Rscript "$SCRIPT_DIR/13_species_recovery_ranking.R" 2>&1
    EXIT2=$?
    echo "Finished: $(date), Exit code: $EXIT2"
    echo ""

    echo "=============================================="
    echo "ALL ANALYSES COMPLETE"
    echo "Finished: $(date)"
    echo "Exit codes: Script 14=$EXIT1, Script 13=$EXIT2"
    echo "=============================================="

    if [ $EXIT1 -eq 0 ] && [ $EXIT2 -eq 0 ]; then
        echo "SUCCESS: All scripts completed without errors"
    else
        echo "WARNING: One or more scripts had errors"
    fi
}

# Check if running in screen mode
if [ "$1" == "screen" ]; then
    LOG_FILE="$LOG_DIR/ranking_analyses_${TIMESTAMP}.log"
    echo "Starting analyses in detached screen session..."
    echo "Log file: $LOG_FILE"
    echo ""
    echo "To monitor progress:"
    echo "  screen -r ranking_analysis"
    echo ""
    echo "To check log:"
    echo "  tail -f $LOG_FILE"
    echo ""

    screen -dmS ranking_analysis bash -c "source ~/.bashrc; $(declare -f run_analyses); run_analyses 2>&1 | tee $LOG_FILE"

    echo "Screen session 'ranking_analysis' started."
else
    # Run in foreground
    run_analyses
fi
