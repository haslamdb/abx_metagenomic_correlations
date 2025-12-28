#!/bin/bash
#
# Database Setup Script for Antibiotic-Microbiome Correlations Pipeline
# ======================================================================
#
# This script downloads and sets up all required databases:
# - Human genome (GRCh38) for host removal
# - Kraken2 Standard database for taxonomic classification
# - MetaPhlAn database for functional profiling
# - HUMAnN databases (ChocoPhlAn and UniRef90)
# - ABRicate databases (CARD, ResFinder, NCBI)
# - AMRFinderPlus database
#
# Usage: bash setup_databases.sh /path/to/database/directory
#
# WARNING: This requires ~300-500 GB of disk space and takes several hours!

set -euo pipefail

# Check arguments
if [ $# -ne 1 ]; then
    echo "Usage: $0 <database_directory>"
    echo "Example: $0 /data/databases"
    exit 1
fi

DB_DIR="$1"
THREADS=${THREADS:-8}

echo "=============================================="
echo "Database Setup for Metagenomics Pipeline"
echo "=============================================="
echo "Database directory: $DB_DIR"
echo "Threads: $THREADS"
echo "=============================================="

# Create directory structure
mkdir -p "$DB_DIR"/{kraken2,humann,metaphlan,host,abricate,amrfinder}

# -----------------------------------------------------------------------------
# 1. Human Genome for Host Removal
# -----------------------------------------------------------------------------
echo ""
echo "[1/6] Setting up human genome (GRCh38) for host removal..."

HOST_DIR="$DB_DIR/host"
cd "$HOST_DIR"

if [ ! -f "GRCh38.1.bt2" ]; then
    echo "Downloading GRCh38 genome..."
    wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

    echo "Decompressing..."
    gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GRCh38.fna

    echo "Building Bowtie2 index (this takes ~1 hour)..."
    bowtie2-build --threads $THREADS GRCh38.fna GRCh38

    echo "Host genome setup complete!"
else
    echo "Host genome already exists, skipping..."
fi

# -----------------------------------------------------------------------------
# 2. Kraken2 Database
# -----------------------------------------------------------------------------
echo ""
echo "[2/6] Setting up Kraken2 database..."

KRAKEN_DIR="$DB_DIR/kraken2"
cd "$KRAKEN_DIR"

if [ ! -f "hash.k2d" ]; then
    echo "Downloading Kraken2 Standard database (~70GB)..."
    echo "This may take several hours depending on your connection..."

    # Option 1: Pre-built database from Ben Langmead's collection
    wget -q https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240112.tar.gz

    echo "Extracting..."
    tar -xzf k2_standard_20240112.tar.gz
    rm k2_standard_20240112.tar.gz

    echo "Kraken2 database setup complete!"
else
    echo "Kraken2 database already exists, skipping..."
fi

# -----------------------------------------------------------------------------
# 3. MetaPhlAn Database
# -----------------------------------------------------------------------------
echo ""
echo "[3/6] Setting up MetaPhlAn database..."

METAPHLAN_DIR="$DB_DIR/metaphlan"

if [ ! -d "$METAPHLAN_DIR/mpa_vJun23_CHOCOPhlAnSGB_202307" ]; then
    echo "Downloading MetaPhlAn database..."

    # Create conda environment temporarily if needed
    metaphlan --install --bowtie2db "$METAPHLAN_DIR"

    echo "MetaPhlAn database setup complete!"
else
    echo "MetaPhlAn database already exists, skipping..."
fi

# -----------------------------------------------------------------------------
# 4. HUMAnN Databases
# -----------------------------------------------------------------------------
echo ""
echo "[4/6] Setting up HUMAnN databases..."

HUMANN_DIR="$DB_DIR/humann"

# ChocoPhlAn database (nucleotide)
if [ ! -d "$HUMANN_DIR/chocophlan" ]; then
    echo "Downloading HUMAnN ChocoPhlAn database (~15GB)..."
    humann_databases --download chocophlan full "$HUMANN_DIR"
    echo "ChocoPhlAn database setup complete!"
else
    echo "ChocoPhlAn database already exists, skipping..."
fi

# UniRef90 database (protein) - THIS IS LARGE (~30GB compressed, ~60GB uncompressed)
if [ ! -d "$HUMANN_DIR/uniref" ]; then
    echo "Downloading HUMAnN UniRef90 database (~30GB compressed)..."
    echo "This is the largest download and may take several hours..."
    humann_databases --download uniref uniref90_diamond "$HUMANN_DIR"
    echo "UniRef90 database setup complete!"
else
    echo "UniRef90 database already exists, skipping..."
fi

# -----------------------------------------------------------------------------
# 5. ABRicate Databases
# -----------------------------------------------------------------------------
echo ""
echo "[5/6] Setting up ABRicate databases..."

echo "Updating ABRicate databases..."
abricate --setupdb

# Update specific databases
for db in card resfinder ncbi argannot; do
    echo "Updating $db database..."
    abricate --db $db --check || abricate --db $db --setupdb
done

echo "ABRicate databases setup complete!"

# -----------------------------------------------------------------------------
# 6. AMRFinderPlus Database
# -----------------------------------------------------------------------------
echo ""
echo "[6/6] Setting up AMRFinderPlus database..."

AMRFINDER_DIR="$DB_DIR/amrfinder"

echo "Downloading AMRFinderPlus database..."
amrfinder --update --database "$AMRFINDER_DIR"

echo "AMRFinderPlus database setup complete!"

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "Database Setup Complete!"
echo "=============================================="
echo ""
echo "Database locations:"
echo "  Host genome:    $DB_DIR/host/GRCh38"
echo "  Kraken2:        $DB_DIR/kraken2"
echo "  MetaPhlAn:      $DB_DIR/metaphlan"
echo "  HUMAnN (nucl):  $DB_DIR/humann/chocophlan"
echo "  HUMAnN (prot):  $DB_DIR/humann/uniref"
echo "  AMRFinder:      $DB_DIR/amrfinder"
echo ""
echo "Update your config/config.yaml with these paths!"
echo ""

# Print disk usage
echo "Disk usage:"
du -sh "$DB_DIR"/*

echo ""
echo "Done!"
