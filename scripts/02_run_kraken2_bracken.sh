#!/bin/bash

# Run Kraken2 and Bracken on BSI/Abx samples
# Process each sample completely (Kraken2 + Bracken) before moving to next
# Uses host-filtered reads from 01_remove_host_reads.sh

# Configuration
KRAKEN2_DB="/shardpool/reference_data/Kraken2DB"
THREADS=24
PROJECT_DIR="/home/david/projects/abx_metagenomic_correlations"
OUTPUT_DIR="/fastpool/analysis/mss_pipeline/kraken2"
SAMPLE_FILE="${PROJECT_DIR}/samples/bsi_abx_samples_complete.csv"
LOG_FILE="${PROJECT_DIR}/logs/kraken2_bracken.log"

# Use host-filtered (nonhost) reads from shared location
NONHOST_READS_DIR="/fastpool/analysis/mss_pipeline/nonhost_reads"

# Create output directories
mkdir -p $OUTPUT_DIR
mkdir -p "${PROJECT_DIR}/logs"

# Initialize log
echo "Kraken2/Bracken Pipeline - BSI/Abx Samples" > $LOG_FILE
echo "Started $(date)" >> $LOG_FILE
echo "Kraken2 DB: $KRAKEN2_DB" >> $LOG_FILE
echo "================================" >> $LOG_FILE

# Count total samples (excluding header)
TOTAL_SAMPLES=$(tail -n +2 $SAMPLE_FILE | wc -l)

echo "Processing $TOTAL_SAMPLES samples with Kraken2 and Bracken..."
echo "Each sample will be fully processed (Kraken2 + Bracken) before moving to next"
echo "================================"
echo ""

CURRENT=0
tail -n +2 $SAMPLE_FILE | while IFS=',' read -r sample_name fastq_path R1_file R2_file naming_pattern; do
    CURRENT=$((CURRENT + 1))

    echo "========================================="
    echo "Sample $CURRENT/$TOTAL_SAMPLES: $sample_name"
    echo "========================================="
    echo "$(date): Processing $sample_name ($CURRENT/$TOTAL_SAMPLES)" >> $LOG_FILE

    # Use host-filtered (nonhost) reads
    R1_NONHOST="${NONHOST_READS_DIR}/${sample_name}_R1.nonhost.fastq.gz"
    R2_NONHOST="${NONHOST_READS_DIR}/${sample_name}_R2.nonhost.fastq.gz"

    # Skip if already processed
    if [[ -f "${OUTPUT_DIR}/${sample_name}_species.bracken" && -f "${OUTPUT_DIR}/${sample_name}_genus.bracken" ]]; then
        echo "Already processed, skipping..."
        echo "$(date): $sample_name - SKIPPED (already processed)" >> $LOG_FILE
        continue
    fi

    # Check if nonhost reads exist
    if [[ ! -f "$R1_NONHOST" || ! -f "$R2_NONHOST" ]]; then
        echo "ERROR: Nonhost reads not found for $sample_name"
        echo "  Expected R1: $R1_NONHOST"
        echo "  Expected R2: $R2_NONHOST"
        echo "  Run 01_remove_host_reads.sh first"
        echo "$(date): $sample_name - ERROR (nonhost reads not found)" >> $LOG_FILE
        continue
    fi

    START_TIME=$(date +%s)

    # Run Kraken2
    echo "Running Kraken2..."
    kraken2 --db $KRAKEN2_DB \
        --threads $THREADS \
        --paired \
        --output ${OUTPUT_DIR}/${sample_name}.kraken \
        --report ${OUTPUT_DIR}/${sample_name}.report \
        "$R1_NONHOST" \
        "$R2_NONHOST" \
        2>> $LOG_FILE

    # Run Bracken - Species level
    echo "Running Bracken (species level)..."
    bracken -d $KRAKEN2_DB \
        -i ${OUTPUT_DIR}/${sample_name}.report \
        -o ${OUTPUT_DIR}/${sample_name}_species.bracken \
        -w ${OUTPUT_DIR}/${sample_name}_species.report \
        -r 150 \
        -l S \
        -t 5 \
        2>> $LOG_FILE

    # Run Bracken - Genus level
    echo "Running Bracken (genus level)..."
    bracken -d $KRAKEN2_DB \
        -i ${OUTPUT_DIR}/${sample_name}.report \
        -o ${OUTPUT_DIR}/${sample_name}_genus.bracken \
        -w ${OUTPUT_DIR}/${sample_name}_genus.report \
        -r 150 \
        -l G \
        -t 5 \
        2>> $LOG_FILE

    # Remove large .kraken output file to save space (keep .report files)
    echo "Removing large .kraken file..."
    rm -f ${OUTPUT_DIR}/${sample_name}.kraken

    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))

    echo "Completed: $sample_name in ${ELAPSED}s"
    echo "$(date): $sample_name - COMPLETED in ${ELAPSED}s" >> $LOG_FILE
    echo ""
done

echo "================================"
echo "Kraken2 and Bracken analysis complete!"
echo "Results saved to $OUTPUT_DIR"
echo "================================"
echo "Pipeline finished $(date)" >> $LOG_FILE
