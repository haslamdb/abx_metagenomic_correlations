#!/usr/bin/env python3
"""
Compute antibiotic exposure for all samples with Bracken data.
Much faster than R for large drug table operations.
"""

import pandas as pd
import numpy as np
from datetime import timedelta
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Paths
PROJECT_DIR = Path("/home/david/projects/abx_metagenomic_correlations")
DATA_DIR = PROJECT_DIR / "data"
LEGACY_DIR = DATA_DIR / "legacy"

print("=== Computing Antibiotic Exposure ===\n")

# =============================================================================
# 1. Load Bracken sample names
# =============================================================================
print("Loading Bracken sample names...")

# Read the genus counts to get sample names (first column is row names)
genus_counts = pd.read_csv(DATA_DIR / "genus_counts_raw.csv", usecols=[0])
bracken_samples = set(genus_counts.iloc[:, 0].tolist())
print(f"  Bracken samples: {len(bracken_samples)}")

# =============================================================================
# 2. Load sample metadata
# =============================================================================
print("Loading sample metadata...")

sample_list = pd.read_csv(LEGACY_DIR / "SampleListFormatted.csv")
sample_list = sample_list[sample_list['SequenceFileName'].isin(bracken_samples)].copy()
print(f"  Samples in SampleListFormatted with Bracken data: {len(sample_list)}")

# Parse dates
sample_list['SampleDate_parsed'] = pd.to_datetime(sample_list['SampleDate'], format='%m/%d/%Y', errors='coerce')
sample_list = sample_list[sample_list['SampleDate_parsed'].notna()].copy()
print(f"  Samples with valid dates: {len(sample_list)}")

# =============================================================================
# 3. Load drug data
# =============================================================================
print("\nLoading Drugs.csv...")

drugs = pd.read_csv(LEGACY_DIR / "Drugs.csv", usecols=['Date', 'MRN', 'Drug', 'Route'])
drugs['Date'] = pd.to_datetime(drugs['Date'])  # Already in YYYY-MM-DD format
drugs['MRN'] = drugs['MRN'].astype(str)
print(f"  Total drug records: {len(drugs):,}")

# Pre-filter to MRNs we care about
sample_mrns = set(sample_list['MRN'].astype(str).unique())
drugs = drugs[drugs['MRN'].isin(sample_mrns)].copy()
print(f"  Drug records for our samples: {len(drugs):,}")

# Create index for fast lookup
drugs = drugs.set_index(['MRN', 'Date']).sort_index()

# =============================================================================
# 4. Define antibiotics
# =============================================================================

KEY_ANTIBIOTICS = [
    'Pip/Tazo', 'Meropenem', 'Cefepime', 'Vancomycin', 'Metronidazole',
    'Ceftriaxone', 'Ciprofloxacin', 'Clindamycin', 'Azithromycin', 'TMP/SMX',
    'Cefazolin', 'Gentamicin', 'Ampicillin', 'Augmentin', 'Levofloxacin'
]

ANAEROBIC_DRUGS = {'Metronidazole', 'Pip/Tazo', 'Meropenem', 'Clindamycin',
                   'Augmentin', 'Ertapenem', 'Amp/Sulb'}
BROAD_SPECTRUM = {'Pip/Tazo', 'Meropenem', 'Cefepime', 'Imipenem', 'Ertapenem',
                  'Ciprofloxacin', 'Levofloxacin'}
GRAM_POSITIVE = {'Vancomycin', 'Daptomycin', 'Linezolid', 'Ceftaroline'}
PSEUDOMONAL = {'Pip/Tazo', 'Meropenem', 'Cefepime', 'Ciprofloxacin',
               'Tobramycin', 'Gentamicin', 'Amikacin', 'Aztreonam'}

# =============================================================================
# 5. Compute exposure for each sample
# =============================================================================
print("\nComputing antibiotic exposure...")

def compute_exposure(mrn, sample_date, window_days):
    """Compute antibiotic exposure for a sample."""
    start_date = sample_date - timedelta(days=window_days)
    end_date = sample_date - timedelta(days=1)  # Exclude sample date itself

    try:
        # Get all drugs for this MRN
        if mrn not in drugs.index.get_level_values(0):
            return None

        mrn_drugs = drugs.loc[mrn]

        # Filter to date range
        mask = (mrn_drugs.index >= start_date) & (mrn_drugs.index <= end_date)
        window_drugs = mrn_drugs[mask]

        if len(window_drugs) == 0:
            return None

        return window_drugs.reset_index()
    except:
        return None

# Process all samples
results = []
n_samples = len(sample_list)

for idx, (_, row) in enumerate(sample_list.iterrows()):
    if idx % 100 == 0:
        print(f"  Processing {idx}/{n_samples}...", end='\r')

    sample_id = row['SequenceFileName']
    mrn = str(row['MRN'])
    sample_date = row['SampleDate_parsed']

    result = {
        'sample_id': sample_id,
        'MRN': mrn,
        'PatientID': row.get('PMID', ''),
        'SampleDate': sample_date,
        'SampleType': row.get('SampleType', ''),
        'PatientGroup': row.get('PatientGroup', ''),
    }

    # Compute for 7-day and 14-day windows
    for window in [7, 14]:
        suffix = f'_{window}d'
        window_drugs = compute_exposure(mrn, sample_date, window)

        if window_drugs is None or len(window_drugs) == 0:
            result[f'abx_any{suffix}'] = False
            result[f'abx_days{suffix}'] = 0
            result[f'abx_anaerobic{suffix}'] = False
            result[f'abx_broad{suffix}'] = False
            result[f'abx_gram_pos{suffix}'] = False
            result[f'abx_pseudomonal{suffix}'] = False
            result[f'drugs{suffix}'] = ''
            for abx in KEY_ANTIBIOTICS:
                col = abx.replace('/', '_').replace('-', '_') + suffix
                result[col] = 0
        else:
            unique_drugs = set(window_drugs['Drug'].unique())
            unique_dates = window_drugs['Date'].nunique()

            result[f'abx_any{suffix}'] = True
            result[f'abx_days{suffix}'] = unique_dates
            result[f'abx_anaerobic{suffix}'] = bool(unique_drugs & ANAEROBIC_DRUGS)
            result[f'abx_broad{suffix}'] = bool(unique_drugs & BROAD_SPECTRUM)
            result[f'abx_gram_pos{suffix}'] = bool(unique_drugs & GRAM_POSITIVE)
            result[f'abx_pseudomonal{suffix}'] = bool(unique_drugs & PSEUDOMONAL)
            result[f'drugs{suffix}'] = '; '.join(sorted(unique_drugs))

            # Individual antibiotics (count days)
            for abx in KEY_ANTIBIOTICS:
                col = abx.replace('/', '_').replace('-', '_') + suffix
                abx_days = window_drugs[window_drugs['Drug'] == abx]['Date'].nunique()
                result[col] = abx_days

    results.append(result)

print(f"  Processed {n_samples}/{n_samples} samples")

# =============================================================================
# 6. Create output dataframe
# =============================================================================
print("\nCreating output dataframe...")

df = pd.DataFrame(results)

# Summary
print("\n=== Summary ===")
print(f"Total samples: {len(df)}")
print(f"Samples with any abx (7d): {df['abx_any_7d'].sum()}")
print(f"Samples with any abx (14d): {df['abx_any_14d'].sum()}")

print("\nIndividual antibiotic exposure (7-day window):")
for abx in KEY_ANTIBIOTICS:
    col = abx.replace('/', '_').replace('-', '_') + '_7d'
    if col in df.columns:
        n_exposed = (df[col] > 0).sum()
        print(f"  {abx}: {n_exposed} samples")

print("\nPatient groups:")
print(df['PatientGroup'].value_counts())

# =============================================================================
# 7. Save results
# =============================================================================
print("\n=== Saving Results ===")

# Save as CSV
output_csv = DATA_DIR / "sample_metadata_with_abx.csv"
df.to_csv(output_csv, index=False)
print(f"Saved: {output_csv}")

# Also save as pickle for faster R loading
output_pkl = DATA_DIR / "sample_metadata_with_abx.pkl"
df.to_pickle(output_pkl)
print(f"Saved: {output_pkl}")

print("\n=== Done ===")
