#!/bin/bash

# Remove host (human/mouse) reads from BSI/Abx samples
# Uses bowtie2 to align against host database and keeps unaligned reads
# Based on preterm_arg workflow

# Configuration
BOWTIE2="/home/david/environments/bowtie2/bin/bowtie2"
SAMTOOLS="/home/david/environments/bowtie2/bin/samtools"
BT2_INDEX="/bulkpool/reference_data/hostBT2/host_human_mouse"
THREADS=24
OUTPUT_DIR="/fastpool/analysis/mss_pipeline/nonhost_reads"
PROJECT_DIR="/home/david/projects/abx_metagenomic_correlations"
SAMPLE_FILE="${PROJECT_DIR}/samples/bsi_abx_samples_complete.csv"
LOG_FILE="${PROJECT_DIR}/logs/host_removal.log"
STATS_FILE="${OUTPUT_DIR}/host_removal_stats_bsi_abx.csv"

# Create output directories
mkdir -p $OUTPUT_DIR
mkdir -p "${PROJECT_DIR}/logs"

# Initialize log file
echo "Host Read Removal Pipeline - BSI/Abx Samples" > $LOG_FILE
echo "Started $(date)" >> $LOG_FILE
echo "Bowtie2 index: $BT2_INDEX" >> $LOG_FILE
echo "Threads: $THREADS" >> $LOG_FILE
echo "================================" >> $LOG_FILE

# Initialize stats CSV (only if it doesn't exist, to support resume)
if [[ ! -f "$STATS_FILE" ]]; then
    echo "sample_name,input_reads,nonhost_reads,host_reads,host_pct,elapsed_seconds" > $STATS_FILE
fi

# Count total samples (excluding header)
TOTAL_SAMPLES=$(tail -n +2 $SAMPLE_FILE | wc -l)

echo "Processing $TOTAL_SAMPLES samples for host read removal..."
echo "Output directory: $OUTPUT_DIR"
echo "================================"
echo ""

CURRENT=0
tail -n +2 $SAMPLE_FILE | while IFS=',' read -r sample_name fastq_path R1_file R2_file naming_pattern; do
    CURRENT=$((CURRENT + 1))

    echo "========================================="
    echo "Sample $CURRENT/$TOTAL_SAMPLES: $sample_name"
    echo "========================================="
    echo "$(date): Processing $sample_name ($CURRENT/$TOTAL_SAMPLES)" >> $LOG_FILE

    R1_INPUT="${fastq_path}/${R1_file}"
    R2_INPUT="${fastq_path}/${R2_file}"
    R1_OUTPUT="${OUTPUT_DIR}/${sample_name}_R1.nonhost.fastq.gz"
    R2_OUTPUT="${OUTPUT_DIR}/${sample_name}_R2.nonhost.fastq.gz"

    # Skip if output already exists
    if [[ -f "$R1_OUTPUT" && -f "$R2_OUTPUT" ]]; then
        echo "Output files already exist, skipping..."
        echo "$(date): $sample_name - SKIPPED (output exists)" >> $LOG_FILE
        continue
    fi

    # Check if input files exist
    if [[ ! -f "$R1_INPUT" || ! -f "$R2_INPUT" ]]; then
        echo "ERROR: Input files not found for $sample_name"
        echo "  R1: $R1_INPUT"
        echo "  R2: $R2_INPUT"
        echo "$(date): $sample_name - ERROR (input not found)" >> $LOG_FILE
        continue
    fi

    START_TIME=$(date +%s)

    # Align to host genome and extract non-host reads in one pipeline
    # -f 12: both reads unmapped
    # -F 256: not secondary alignment
    echo "Aligning to host genome and extracting non-host reads..."

    $BOWTIE2 -x $BT2_INDEX \
        -1 "$R1_INPUT" \
        -2 "$R2_INPUT" \
        --very-sensitive \
        -p $THREADS \
        --un-conc-gz "${OUTPUT_DIR}/${sample_name}_R%.nonhost.fastq.gz" \
        2>> $LOG_FILE \
        | $SAMTOOLS view -bS - > /dev/null

    # Check if output was created
    if [[ -f "${OUTPUT_DIR}/${sample_name}_R1.nonhost.fastq.gz" && -f "${OUTPUT_DIR}/${sample_name}_R2.nonhost.fastq.gz" ]]; then
        END_TIME=$(date +%s)
        ELAPSED=$((END_TIME - START_TIME))

        # Get read counts for reporting
        INPUT_READS=$(zcat "$R1_INPUT" 2>/dev/null | wc -l | awk '{print $1/4}')
        OUTPUT_READS=$(zcat "$R1_OUTPUT" 2>/dev/null | wc -l | awk '{print $1/4}')
        HOST_READS=$((INPUT_READS - OUTPUT_READS))
        HOST_PCT=$(echo "scale=2; $HOST_READS * 100 / $INPUT_READS" | bc 2>/dev/null || echo "N/A")

        echo "Completed in ${ELAPSED}s"
        echo "  Input reads: $INPUT_READS"
        echo "  Non-host reads: $OUTPUT_READS"
        echo "  Host reads removed: $HOST_READS ($HOST_PCT%)"
        echo "$(date): $sample_name - COMPLETED in ${ELAPSED}s (removed $HOST_PCT% host)" >> $LOG_FILE
        echo "${sample_name},${INPUT_READS},${OUTPUT_READS},${HOST_READS},${HOST_PCT},${ELAPSED}" >> $STATS_FILE
    else
        echo "ERROR: Failed to create output files for $sample_name"
        echo "$(date): $sample_name - ERROR (output not created)" >> $LOG_FILE
    fi

    echo ""
done

echo "================================"
echo "Host read removal complete!"
echo "Results saved to $OUTPUT_DIR"
echo "Log file: $LOG_FILE"
echo "================================"
echo "Pipeline finished $(date)" >> $LOG_FILE
