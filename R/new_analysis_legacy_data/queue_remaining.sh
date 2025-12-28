#!/bin/bash
# Wait for script 05 to finish, then run remaining scripts

cd /home/david/projects/abx_metagenomic_correlations

# Wait for 05_individual_antibiotic_effects.R to finish
echo "Waiting for species-level analysis (05) to complete..."
while pgrep -f "05_individual_antibiotic_effects.R" > /dev/null; do
    sleep 60
done

echo "Species-level analysis complete at $(date)"
echo "Starting remaining analyses..."

# Run the remaining scripts in the existing screen session
screen -S abx_analysis -X stuff "./R/new_analysis_legacy_data/run_remaining_analyses.sh\n"
