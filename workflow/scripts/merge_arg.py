#!/usr/bin/env python3
"""
Merge ARG detection results into an abundance matrix.

Input: List of ABRicate and AMRFinderPlus output files
Output: Combined ARG abundance matrix (samples x ARG families)
"""

import pandas as pd
from pathlib import Path
import sys
from collections import defaultdict


def load_abricate_file(filepath):
    """Load ABRicate output and aggregate by gene family."""
    try:
        df = pd.read_csv(filepath, sep="\t")
        if len(df) == 0:
            return {}

        # Aggregate by GENE (ARG family)
        gene_counts = df.groupby("GENE").size().to_dict()
        return gene_counts
    except Exception as e:
        print(f"Warning: Could not load {filepath}: {e}", file=sys.stderr)
        return {}


def load_amrfinder_file(filepath):
    """Load AMRFinderPlus output and aggregate by gene symbol."""
    try:
        df = pd.read_csv(filepath, sep="\t")
        if len(df) == 0:
            return {}

        # Aggregate by Gene symbol
        gene_counts = df.groupby("Gene symbol").size().to_dict()
        return gene_counts
    except Exception as e:
        print(f"Warning: Could not load {filepath}: {e}", file=sys.stderr)
        return {}


def extract_sample_name(filepath):
    """Extract sample name from filepath."""
    name = Path(filepath).stem
    # Remove suffixes like _abricate_card, _amrfinder
    for suffix in ["_abricate_card", "_abricate_resfinder", "_abricate_ncbi", "_amrfinder"]:
        name = name.replace(suffix, "")
    return name


def merge_arg_files(card_files, amrfinder_files, output_file):
    """Merge ARG detection results into a matrix."""

    # Dictionary to store results: {sample: {gene: count}}
    results = defaultdict(lambda: defaultdict(int))

    # Process CARD results
    for filepath in card_files:
        sample = extract_sample_name(filepath)
        gene_counts = load_abricate_file(filepath)
        for gene, count in gene_counts.items():
            results[sample][f"CARD_{gene}"] += count

    # Process AMRFinder results
    for filepath in amrfinder_files:
        sample = extract_sample_name(filepath)
        gene_counts = load_amrfinder_file(filepath)
        for gene, count in gene_counts.items():
            results[sample][f"AMR_{gene}"] += count

    # Convert to DataFrame
    df = pd.DataFrame.from_dict(results, orient="index")
    df = df.fillna(0).astype(int)
    df = df.sort_index()

    # Save
    df.to_csv(output_file, sep="\t")

    print(f"Merged ARG results: {len(df)} samples, {len(df.columns)} genes")


if __name__ == "__main__":
    # When run by Snakemake
    merge_arg_files(
        card_files=snakemake.input.card,
        amrfinder_files=snakemake.input.amrfinder,
        output_file=snakemake.output[0]
    )
