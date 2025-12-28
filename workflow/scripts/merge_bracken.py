#!/usr/bin/env python3
"""
Merge Bracken abundance files into a single matrix.

Input: List of Bracken output files
Output: Combined abundance matrix (samples x taxa)
"""

import pandas as pd
from pathlib import Path
import sys


def load_bracken_file(filepath):
    """Load a single Bracken output file."""
    df = pd.read_csv(filepath, sep="\t")

    # Extract sample name from filename
    sample_name = Path(filepath).stem.replace("_species", "").replace("_genus", "").replace("_phylum", "")

    # Create series with taxon as index
    abundance = df.set_index("name")["fraction_total_reads"]
    abundance.name = sample_name

    return abundance


def merge_bracken_files(input_files, output_file):
    """Merge multiple Bracken files into a matrix."""

    # Load all files
    abundance_series = []
    for filepath in input_files:
        try:
            series = load_bracken_file(filepath)
            abundance_series.append(series)
        except Exception as e:
            print(f"Warning: Could not load {filepath}: {e}", file=sys.stderr)

    # Combine into DataFrame
    merged_df = pd.concat(abundance_series, axis=1)

    # Fill missing values with 0
    merged_df = merged_df.fillna(0)

    # Transpose so samples are rows
    merged_df = merged_df.T

    # Sort by sample name
    merged_df = merged_df.sort_index()

    # Save
    merged_df.to_csv(output_file, sep="\t")

    print(f"Merged {len(abundance_series)} samples, {len(merged_df.columns)} taxa")


if __name__ == "__main__":
    # When run by Snakemake
    merge_bracken_files(
        input_files=snakemake.input,
        output_file=snakemake.output[0]
    )
