#!/bin/bash
# Transfer legacy Kraken/Bracken files from lambda-quad
# Source: lambda-quad (192.168.1.13)
# Destination: /home/david/projects/abx_metagenomic_correlations/data/kraken2_legacy/

SOURCE_HOST="david@192.168.1.13"
SOURCE_PORT="2424"
SOURCE_DIR="~/Documents/Alignments/KrakenAlignments"
DEST_DIR="/home/david/projects/abx_metagenomic_correlations/data/kraken2_legacy"

echo "=== Transferring Legacy Kraken/Bracken Files ==="
echo "Source: ${SOURCE_HOST}:${SOURCE_DIR}"
echo "Destination: ${DEST_DIR}"
echo ""

# Create destination directory
mkdir -p "$DEST_DIR"

# Transfer all relevant files
# - *_genus_abundance.txt (Bracken genus-level)
# - *_species_abundance.txt (Bracken species-level)
# - *.kraken.report (Kraken reports)
# - *.classic_bracken.report (Bracken reports)

SSH_CMD="ssh -p ${SOURCE_PORT}"

echo "Transferring Bracken genus abundance files..."
mkdir -p "${DEST_DIR}/genus"
rsync -avP -e "${SSH_CMD}" --include="*_genus_abundance.txt" --exclude="*" \
    "${SOURCE_HOST}:${SOURCE_DIR}/" "${DEST_DIR}/genus/"

echo ""
echo "Transferring Bracken species abundance files..."
mkdir -p "${DEST_DIR}/species"
rsync -avP -e "${SSH_CMD}" --include="*_species_abundance.txt" --exclude="*" \
    "${SOURCE_HOST}:${SOURCE_DIR}/" "${DEST_DIR}/species/"

echo ""
echo "Transferring Kraken reports..."
mkdir -p "${DEST_DIR}/reports"
rsync -avP -e "${SSH_CMD}" --include="*.kraken.report" --exclude="*" \
    "${SOURCE_HOST}:${SOURCE_DIR}/" "${DEST_DIR}/reports/"

echo ""
echo "=== Transfer Complete ==="
echo "Files saved to: ${DEST_DIR}"
ls -la "${DEST_DIR}"
