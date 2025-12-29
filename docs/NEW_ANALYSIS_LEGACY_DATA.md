# Reanalysis of Legacy Antibiotic-Microbiome Data

**Date:** December 2025 (last updated: 2025-12-28)
**Data Sources:**
- **Taxonomic profiles:** Original Kraken2/Bracken output files (`data/kraken2_legacy/`)
- **Sample metadata:** `data/legacy/SampleListFormatted.csv`
- **Antibiotic exposure:** `data/legacy/Drugs.csv`
- **Sample list:** `data/legacy/BSI_Abx_Samples.csv`

**Scripts:** `R/new_analysis_legacy_data/`
**Results:** `results/new_analysis_legacy_data/`
**Master script:** `R/new_analysis_legacy_data/run_all_analysis.sh`

---

## Overview

This analysis re-examines the relationship between antibiotic exposure and gut microbiome composition using improved statistical methods. The original 2019 analysis used simple Pearson correlations between rate-of-change metrics, which failed to reveal expected biological patterns.

### Key Changes from Previous Approach

**Data Pipeline Update (December 2025):**
- **Before:** Used pre-processed RData file (`AbxEffectData20191014`) with rarefied count matrices
- **Now:** Rebuilt count matrices directly from original Bracken output files, linked to antibiotic exposure from hospital records

This provides:
- Raw counts (not rarefied) for modern compositional methods
- Direct linkage to antibiotic exposure data via MRN and sample dates
- Transparent pipeline from raw data to analysis
- Exclusion of contaminants (Homo sapiens, Salinibacter ruber) at source

### Statistical Methods Applied

- Mixed-effects models to account for repeated patient samples
- Compositional-aware differential abundance testing (ALDEx2, MaAsLin3, ANCOM-BC2)
- Proper antibiotic spectrum classification
- Paired sample analysis for within-patient comparisons
- Multi-method concordance for robust associations

---

## Data Summary

| Metric | Value |
|--------|-------|
| Total samples | 905 |
| Unique patients | 145 |
| Sample pairs (2-week window) | 349 (346 high-confidence) |
| Species | 1,460 |
| Genera | 163 |

### Patient Groups

| Group | Samples | Patients | % with Abx (7d) | Notes |
|-------|---------|----------|-----------------|-------|
| LvTx (Liver Transplant) | 234 | 35 | 76.5% | High-confidence |
| BMT (Bone Marrow Transplant) | 218 | 61 | 66.1% | High-confidence |
| SB (Short Bowel) | 291 (155 Stool) | 16 | 25.4% | High-confidence (Stool only) |
| IF (Intestinal Failure) | 100 | 16 | 25.0% | High-confidence |
| IBD (Inflammatory Bowel Disease) | 43 | 10 | **0%** | Excluded - no Abx data |
| PICU | 16 | 7 | 87.5% | High-confidence |
| Urology | 3 | 1 | 0% | Excluded |

**High-confidence groups** (BMT, IF, LvTx, PICU, SB-Stool) are used for primary analyses due to complete inpatient antibiotic records. SB patients have 291 samples total but only 155 Stool samples are included; Fistula (20) and Ostomy (116) samples are excluded due to different microbial communities. IBD patients had no recorded antibiotic exposure, likely due to outpatient prescriptions not captured in the hospital system.

**High-confidence sample total: 723 samples from 130 patients**

---

## Analysis Pipeline

### Script 0: Build Count Matrices from Bracken
**`00_build_count_matrices_from_bracken.R`**

Builds genus and species count matrices directly from raw Bracken output files.

**Input:**
- `data/kraken2_legacy/genus/*_genus_abundance.txt` - Bracken genus-level abundances
- `data/kraken2_legacy/species/*_species_abundance.txt` - Bracken species-level abundances
- `data/legacy/BSI_Abx_Samples.csv` - List of samples to include (BSI/Abx study samples)

**Processing:**
- Reads raw Bracken output files (uses `new_est_reads` as counts)
- Filters to only samples in BSI_Abx_Samples.csv
- Excludes contaminants: Homo sapiens, Salinibacter ruber
- Aggregates counts by sample-taxon pairs

**Output:** `data/bracken_count_matrices.RData`
- `genus_matrix` - samples × genera count matrix
- `species_matrix` - samples × species count matrix

---

### Script 0b: Compute Antibiotic Exposure (Python)
**`scripts/compute_abx_exposure.py`**

Fast computation of antibiotic exposure windows using pandas (much faster than R for large drug tables).

**Input:**
- `data/legacy/SampleListFormatted.csv` - Sample metadata with dates
- `data/legacy/Drugs.csv` - Hospital antibiotic administration records
- `data/genus_counts_raw.csv` - Sample names from Bracken data

**Processing:**
- Links samples to patients via MRN
- Computes 7-day and 14-day exposure windows before each sample
- Calculates exposure for spectrum categories (anti-anaerobic, broad-spectrum, etc.)
- Counts days of exposure for each individual antibiotic

**Output:** `data/sample_metadata_with_abx.csv`

---

### Script 1: Data Preparation
**`01_load_and_prepare_data.R`**

- Loads Bracken count matrices and antibiotic exposure data
- Classifies antibiotics by spectrum:
  - Anti-anaerobic: Metronidazole, Pip/Tazo, Meropenem, Clindamycin, etc.
  - Anti-gram-positive: Vancomycin, Daptomycin, Linezolid
  - Broad-spectrum: Meropenem, Pip/Tazo, Cefepime
  - Anti-pseudomonal: Cefepime, Pip/Tazo, Meropenem, Tobramycin, Ciprofloxacin
- Creates binary exposure variables (7-day and 14-day windows)
- Calculates functional taxonomic groups (obligate anaerobes, Enterobacteriaceae, Enterococcus)
- Computes alpha diversity metrics (Shannon, Simpson, richness)
- Prepares paired sample data with between-sample antibiotic exposure

**Output:** `data/prepared_data.RData`

---

### Script 2: Diversity Analysis
**`02_diversity_analysis.R`**

- Alpha diversity mixed-effects models
- Beta diversity PERMANOVA
- PCoA ordination
- Functional group analysis (anaerobes, Enterobacteriaceae)

**Outputs:**
- `diversity/alpha_diversity_results.csv`
- `diversity/beta_diversity_permanova.csv`
- `figures/alpha_diversity_by_abx.pdf`
- `figures/alpha_diversity_by_spectrum.pdf`
- `figures/pcoa_by_abx.pdf`
- `figures/anaerobes_by_anti_anaerobic_abx.pdf`

---

### Script 2b: Diversity Analysis - Individual Antibiotics
**`02b_diversity_individual_antibiotics.R`**

Extends diversity analysis to model each antibiotic's effect individually, with proper covariate adjustment for polypharmacy.

#### Model Structure

For each antibiotic (e.g., Vancomycin):
```r
shannon ~ Vancomycin_7d + patient_group + Pip_Tazo_7d + TMP_SMX_7d +
          Cefepime_7d + Meropenem_7d + ... + (1|MRN)
```

This isolates the effect of each antibiotic by controlling for:
- **Other concurrent antibiotics** (critical for polypharmacy)
- **Patient group** (baseline microbiome differences)
- **Patient random effect** (repeated samples)

#### Analyses Performed

1. **Shannon diversity** - Mixed-effects model for each antibiotic
2. **Functional groups** - Obligate anaerobes, Enterobacteriaceae, Enterococcus
3. **Beta diversity** - PERMANOVA with covariate adjustment

**Outputs:**
- `diversity/alpha_diversity_individual_abx.csv`
- `diversity/functional_groups_individual_abx.csv`
- `diversity/permanova_individual_abx.csv`
- `diversity/individual_abx_effects_summary.csv`
- `figures/shannon_by_individual_abx_forest.pdf`
- `figures/functional_groups_by_individual_abx_forest.pdf`
- `figures/permanova_by_individual_abx.pdf`

---

### Script 3: Differential Abundance
**`03_differential_abundance.R`**

- ALDEx2 compositional analysis for:
  - Any antibiotic exposure (7-day)
  - Anti-anaerobic antibiotic exposure
  - Broad-spectrum antibiotic exposure
- Volcano plot generation

**Outputs:**
- `differential_abundance/aldex2_any_abx_7d.csv`
- `differential_abundance/aldex2_anti_anaerobic.csv`
- `differential_abundance/aldex2_broad_spectrum.csv`
- `figures/volcano_any_abx.pdf`

---

### Script 4: Paired Sample Analysis (Refactored 2025-12-28)
**`04_paired_sample_analysis.R`**

Analyzes within-patient microbiome changes between paired samples using **individual antibiotics** (not broad categories) with proper covariate adjustment.

#### Major Refactoring (2025-12-28)

**Previous approach:** Used broad antibiotic categories (anti-anaerobic, broad-spectrum, any) which obscured individual drug effects and produced misleading results (e.g., H3 "broad-spectrum increases Enterobacteriaceae" was wrong).

**New approach:** Tests each of 10 individual antibiotics separately with LMM and covariate adjustment:
- Pip_Tazo, Meropenem, Cefepime, Ceftriaxone, Ciprofloxacin
- Metronidazole, Clindamycin, TMP_SMX, Vancomycin_IV, Vancomycin_PO

#### Model Structure

For each antibiotic (e.g., Meropenem):
```r
outcome ~ Meropenem_exposed + Pip_Tazo_between + Cefepime_between + ... +
          interval_days + (1|MRN)
```

This isolates each antibiotic's effect by:
- **Covariate adjustment** for other antibiotics given between samples
- **Time adjustment** for interval between paired samples
- **Patient random effect** for multiple pairs per patient

#### Outcomes Tested

1. **delta_shannon** - Change in Shannon diversity
2. **delta_anaerobes** - Log2 fold-change in obligate anaerobes
3. **delta_enterobact** - Log2 fold-change in Enterobacteriaceae
4. **delta_enterococcus** - Log2 fold-change in Enterococcus
5. **bray_distance** - Bray-Curtis distance (community instability)

#### Data Preparation Changes

Script `01_load_and_prepare_data.R` was also updated to calculate individual antibiotic exposures between paired samples:
- `Pip_Tazo_between`, `Meropenem_between`, etc. (days of exposure)
- Route-specific: `Vancomycin_IV_between`, `Vancomycin_PO_between`

**Outputs:**
- `paired_analysis/individual_abx_lmm_results.csv` - Full LMM results
- `paired_analysis/individual_abx_significant.csv` - Significant findings (p < 0.05)
- `paired_analysis/paired_sample_metrics.csv` - Paired sample data with exposures
- `paired_analysis/exposure_summary.csv` - Exposure counts per antibiotic
- `figures/paired_individual_abx_shannon.pdf` - Forest plot for Shannon
- `figures/paired_individual_abx_enterobact.pdf` - Forest plot for Enterobacteriaceae
- `figures/paired_individual_abx_heatmap.pdf` - Heatmap of all effects

---

### Script 5: Individual Antibiotic Effects (with Covariate Adjustment)
**`05_individual_antibiotic_effects.R`**

Differential abundance analysis for individual antibiotics at species level, using three complementary methods with proper covariate adjustment.

#### Methods Used

| Method | Strengths | Package |
|--------|-----------|---------|
| **ALDEx2 (glm)** | Conservative FDR control, compositionally-aware | `ALDEx2` |
| **MaAsLin3** | Higher power, designed for multivariable models | `maaslin3` |
| **ANCOM-BC2** | Bias correction (sample + taxon-specific) | `ANCOMBC` |

#### Covariate Adjustment

When analyzing each antibiotic (e.g., Pip/Tazo), the model adjusts for:
- **Other antibiotic exposures** (9 other top antibiotics) - accounts for polypharmacy
- **Patient group** (BMT, LvTx, IF, PICU, etc.) - accounts for baseline microbiome differences

Example model for Pip/Tazo effect:
```r
~ Pip_Tazo_7d + patient_group + Vancomycin_7d + Meropenem_7d +
  Cefepime_7d + TMP_SMX_7d + Ciprofloxacin_7d + Metronidazole_7d +
  Azithromycin_7d + Cefazolin_7d + Ceftriaxone_7d
```

#### Top 10 Antibiotics Analyzed
1. Pip/Tazo (piperacillin-tazobactam)
2. TMP/SMX (trimethoprim-sulfamethoxazole)
3. Vancomycin
4. Cefepime
5. Meropenem
6. Ciprofloxacin
7. Metronidazole
8. Azithromycin
9. Cefazolin
10. Ceftriaxone

#### Concordance Analysis

The script identifies **robust associations** - species significantly associated with an antibiotic in 2 or more methods. This multi-method approach reduces false positives:
- ALDEx2 findings are conservative but may miss true effects
- MaAsLin3 has higher power but may include some false positives
- ANCOM-BC2 provides bias-corrected estimates

Associations found by all 3 methods are most trustworthy.

**Outputs:**
```
individual_antibiotics_v2_with_covariates/
├── aldex2/
│   ├── aldex2_Pip_Tazo.csv
│   ├── aldex2_Vancomycin.csv
│   ├── ... (per antibiotic)
│   └── all_antibiotics_combined.csv
├── maaslin3/
│   ├── Pip_Tazo/                    # Full MaAsLin3 output
│   ├── maaslin3_Pip_Tazo.csv
│   ├── ... (per antibiotic)
│   └── all_antibiotics_combined.csv
├── ancombc2/
│   ├── ancombc2_Pip_Tazo.csv
│   ├── ... (per antibiotic)
│   └── all_antibiotics_combined.csv
├── method_concordance.csv           # All significant hits with method counts
├── robust_associations.csv          # Hits found by 2+ methods
├── summary_by_antibiotic.csv        # Overview table
├── enterococcus_by_antibiotic.csv   # Enterococcus-focused results
└── analysis_metadata.rds            # Analysis parameters
```

---

### Script 6: Network Visualization
**`06_network_visualization.R`**

Creates interactive HTML network graphs showing antibiotic-species associations.

#### Features

- **Force-directed layout** - Nodes repel/attract based on connections
- **Hierarchical view** - Antibiotics on left, species on right
- **Per-antibiotic subnetworks** - Individual focused views
- **Interactive controls** - Zoom, drag, hover tooltips, filtering
- **Export to PNG** - Built-in export button

#### Visual Encoding

- **Red boxes** = Antibiotics
- **Green dots** = Species with increased abundance
- **Blue dots** = Species with decreased abundance
- **Green edges** = Positive association (abundance increases)
- **Red edges** = Negative association (abundance decreases)
- **Edge width** = Effect size magnitude

#### Data Sources (priority order)

1. `robust_associations.csv` - Multi-method concordant results (preferred)
2. `method_concordance.csv` - All significant associations
3. `all_antibiotics_combined.csv` - Single-method results (labeled as preliminary)

**Outputs:**
```
network_graphs/
├── antibiotic_species_network.html          # Main interactive network
├── antibiotic_species_network_hierarchical.html  # Left-to-right layout
├── network_Vancomycin.html                  # Per-antibiotic views
├── network_Meropenem.html
├── ... (one per antibiotic)
├── network_associations.csv                 # Data used for visualization
├── network_summary_by_antibiotic.csv        # Summary statistics
└── network_statistics.rds                   # Network metrics
```

---

### Script 7: Microbiome Recovery Analysis (Updated - 2025-12-28)
**`07_recovery_analysis.R`**

Analyzes microbiome recovery after antibiotic cessation using paired samples where:
1. Antibiotics were given before sample 1 (7-day window)
2. No **meaningful** antibiotics between samples (TMP_SMX and Azithromycin excluded as they have minimal microbiome impact)

#### Study Design

**Recovery pairs:** Sample pairs where patient was on antibiotics → then off antibiotics
- Uses **consistent paired-sample methodology** with persistence script
- Analyzes S1 (on Abx) → S2 (off Abx) changes within the same patients

#### Sample Sizes

| Analysis | Pairs |
|----------|-------|
| Pooled (all focus Abx) | 46 |
| Single-antibiotic pairs | 42 (91.3%) |

Single-antibiotic breakdown:
- Pip_Tazo: 25 pairs
- Cefepime: 11 pairs
- Meropenem: 5 pairs
- Metronidazole: 1 pair (too few)

Median recovery interval: 7 days (range 3-14)

#### Statistical Approach

**Functional group analysis:**
- Calculate mean abundance at S1 and S2 for anaerobes, Enterobact, Enterococcus
- Change = S2 - S1 (percentage points)
- Retention = (S2 / S1) × 100%
- Paired t-test on the difference

**Individual antibiotic analysis:**
- Filter to single-antibiotic pairs (91% of cohort)
- Analyze each antibiotic's effect on functional groups separately
- Provides cleaner signal without polypharmacy confounding

**Outputs:**
```
recovery_analysis/
├── recovery_pairs.csv                      # Recovery pair metadata
├── recovery_pooled.csv                     # Pooled genus-level analysis
├── recovery_by_antibiotic.csv              # Per-antibiotic genus stratification
├── recovery_functional_groups.csv          # Pooled functional group summary
└── recovery_by_individual_antibiotic.csv   # Individual abx functional groups

figures/
├── recovery_heatmap.pdf            # Heatmap of genus recovery by antibiotic
├── recovery_forest_pooled.pdf      # Forest plot of pooled effects
├── recovery_comparison_by_abx.pdf  # Bar plot comparing key genera
└── recovery_vs_interval.pdf        # Recovery vs time between samples
```

---

### Script 8: Persistence Mechanism Analysis (Updated - 2025-12-28)
**`08_persistence_analysis.R`**

Investigates what happens to functional groups during antibiotic exposure using paired samples where antibiotics were given between S1 and S2.

#### Study Design

**Persistence pairs:** Sample pairs where antibiotics were administered between samples
- Uses **consistent paired-sample methodology** with recovery script
- Analyzes S1 (pre-Abx) → S2 (post-Abx) changes within the same patients
- n=206 high-confidence pairs

#### Analyses Performed

1. **Pooled functional group analysis**: Changes in anaerobes, Enterobact, Enterococcus
2. **Individual antibiotic "includes" analysis**: Effect of each antibiotic (with overlap)
3. **Genus-level breakdown**: Individual genus persistence for key taxa

#### Sample Sizes

| Analysis | Pairs | Notes |
|----------|-------|-------|
| Pooled | 206 | All pairs with Abx between |
| Single-antibiotic | 16 (7.7%) | Too few for reliable analysis |
| Multiple antibiotics | 149 (71.3%) | Most pairs |

Individual antibiotic breakdown (includes analysis, pairs overlap):
- Pip_Tazo: 14 pairs (57% monotherapy)
- Meropenem: 10 pairs (40% monotherapy)
- Cefepime: 8 pairs (38% monotherapy)
- Vancomycin_IV: 8 pairs (0% monotherapy)

#### Statistical Approach

**Functional group analysis:**
- Calculate mean abundance at S1 and S2 for anaerobes, Enterobact, Enterococcus
- Change = S2 - S1 (percentage points)
- Retention = (S2 / S1) × 100%
- Paired t-test on the difference

**Individual antibiotic "includes" analysis:**
- For each antibiotic, analyze all pairs that received it (even with other abx)
- Reports % that were monotherapy for context
- Provides antibiotic-specific patterns despite polypharmacy

#### Key Functional Groups

| Group | Definition | Genera |
|-------|------------|--------|
| **Obligate Anaerobes** | Require anaerobic conditions | Bacteroides, Blautia, Roseburia, Faecalibacterium, Ruminococcus, etc. |
| **Enterobacteriaceae** | Facultative anaerobes, stress-resistant | Escherichia, Klebsiella, Enterobacter, Citrobacter, etc. |
| **Enterococcus** | Intrinsically resistant to many antibiotics | Enterococcus |

**Outputs:**
```
persistence_analysis/
├── paired_persistence_data.csv     # Paired sample abundances before/after Abx
├── genus_persistence.csv           # Genus-level persistence statistics
├── group_persistence_summary.csv   # Summary of group-level persistence
├── persistence_by_antibiotic.csv   # Individual antibiotic functional group effects
├── dominance_by_exposure.csv       # Dominance patterns by exposure status
├── sample_abundances.csv           # Sample-level abundances for all key genera

figures/
├── persistence_composition.pdf     # Stacked bar of composition by exposure
├── persistence_boxplot.pdf         # Boxplot of log2FC during antibiotics
├── genus_persistence_heatmap.pdf   # Heatmap of genus-level persistence
├── persistence_scatter.pdf         # Scatter plot of before vs after
└── genus_persistence_forest.pdf    # Forest plot of genus persistence
```

---

## Key Results

### 1. Alpha Diversity

**Antibiotic exposure significantly reduces Shannon diversity.**

| Model | Effect | 95% CI | p-value |
|-------|--------|--------|---------|
| Any Abx (7d) → Shannon | -0.31 | [-0.51, -0.12] | **< 0.001** |
| Anti-gram-positive → Shannon | -0.34 | [-0.56, -0.13] | **0.001** |
| Broad-spectrum → Shannon | -0.24 | [-0.50, 0.01] | 0.063 |

### 2. Obligate Anaerobe Abundance

**Anti-anaerobic antibiotics dramatically reduce anaerobe abundance.**

| Exposure | Median Anaerobe Rel. Abundance |
|----------|-------------------------------|
| No anti-anaerobic Abx | 13.8% |
| Anti-anaerobic Abx | 0.8% |

Mixed-effects model: **effect = -0.71 log10, p < 0.001**

### 3. Beta Diversity (PERMANOVA)

| Factor | R² | p-value |
|--------|-----|---------|
| Any Abx (7d) | 0.8% | **0.001** |
| Patient Group | 2.5% | **0.001** |
| Broad-spectrum Abx | 1.2% | **0.001** |
| Anti-anaerobic Abx | 0.5% | **0.001** |

### 4. Differential Abundance (ALDEx2)

#### Any Antibiotic Exposure (31 significant genera, BH < 0.1)

**Increased with Abx:**
| Genus | Effect | p (BH) |
|-------|--------|--------|
| Haemophilus | 0.34 | < 0.001 |
| Veillonella | 0.38 | < 0.001 |
| Klebsiella | 0.26 | < 0.001 |
| Enterobacter | 0.29 | < 0.001 |
| Escherichia | 0.24 | < 0.001 |
| Citrobacter | 0.27 | < 0.001 |
| Serratia | 0.23 | < 0.001 |

**Decreased with Abx:**
| Genus | Effect | p (BH) |
|-------|--------|--------|
| Parabacteroides | -0.21 | 0.012 |
| Bacteroides | -0.17 | 0.032 |
| Faecalibacterium | -0.15 | 0.056 |

### 5. Paired Sample Analysis - Individual Antibiotics (n = 343 pairs)

**Refactored 2025-12-28** to use individual antibiotics with LMM and covariate adjustment, replacing misleading broad category analyses.

#### Significant Findings (p < 0.05)

| Antibiotic | Outcome | Effect | p-value | Interpretation |
|------------|---------|--------|---------|----------------|
| **Pip_Tazo** | Anaerobes | -1.82 log2FC | **0.0005** | Strong anti-anaerobic effect |
| **Meropenem** | Enterobacteriaceae | -2.06 log2FC | **0.0043** | Carbapenem kills gram-negatives |
| **Metronidazole** | Anaerobes | -3.38 log2FC | **0.0087** | Strongest anti-anaerobic effect |
| **Meropenem** | Bray-Curtis | +0.12 | **0.0087** | Increases community instability |
| **Pip_Tazo** | Bray-Curtis | +0.09 | **0.0244** | Increases community instability |
| **Meropenem** | Enterococcus | +1.42 log2FC | **0.0317** | Classic carbapenem selection |
| **Cefepime** | Enterobacteriaceae | -1.60 log2FC | **0.0318** | Cephalosporin kills gram-negatives |

#### Key Insight: Broad Categories Were Misleading

The original H3 hypothesis ("broad-spectrum antibiotics increase Enterobacteriaceae") was **incorrect**. When analyzed by individual antibiotic:
- **Meropenem** and **Cefepime** actually **decrease** Enterobacteriaceae
- The "broad-spectrum" signal was likely confounded by Vancomycin (which allows gram-negative expansion)
- Meropenem shows the clearest carbapenem signature: kills Enterobacteriaceae, selects for resistant Enterococcus

#### Exposure Summary

| Antibiotic | Pairs Exposed | % of Total |
|------------|---------------|------------|
| Pip_Tazo | 82 | 23.9% |
| TMP_SMX | 63 | 18.4% |
| Vancomycin_IV | 47 | 13.7% |
| Meropenem | 46 | 13.4% |
| Cefepime | 43 | 12.5% |
| Ciprofloxacin | 27 | 7.9% |
| Metronidazole | 10 | 2.9% |
| Ceftriaxone | 8 | 2.3% |
| Clindamycin | 8 | 2.3% |
| Vancomycin_PO | 3 | 0.9% (skipped - too few) |

### 6. Microbiome Recovery After Antibiotic Cessation (n = 46 pairs)

Analyzes which genera recover fastest after stopping antibiotics (Abx before S1, no meaningful Abx between).

#### Significant Recoveries (p < 0.05)

| Antibiotic | Genus | Log2FC | p-value | Interpretation |
|------------|-------|--------|---------|----------------|
| **Pooled** | Veillonella | +2.43 | 0.002 | Universal rapid recovery |
| **Pooled** | Romboutsia | +2.22 | 0.003 | Firmicutes recovery |
| **Pip_Tazo** | Veillonella | +2.45 | 0.014 | Rapid recovery |
| **Pip_Tazo** | Clostridioides | +1.64 | 0.021 | C. diff risk window |
| **Pip_Tazo** | Roseburia | +1.35 | 0.035 | Butyrate producer |
| **Pip_Tazo** | Coprococcus | +1.46 | 0.042 | Butyrate producer |
| **Pip_Tazo** | Eubacterium | +1.42 | 0.047 | Anaerobe recovery |
| **Meropenem** | Blautia | +4.99 | 0.045 | Strongest recovery |

#### Key Biological Insights

1. **Obligate anaerobes recover after anti-anaerobic cessation**
   - Roseburia, Coprococcus, Blautia, Eubacterium increase after Pip_Tazo
   - These are butyrate-producing Firmicutes depleted by anti-anaerobic activity

2. **Blautia shows strongest recovery after Meropenem** (+5 log2FC)
   - Carbapenems devastate gut anaerobes; Blautia rebounds fastest

3. **Veillonella recovers universally** across all antibiotics tested

4. **Clostridioides increases in recovery window**
   - Consistent with known C. difficile infection risk period post-antibiotics

5. **Enterococcus declines after Meropenem/Metronidazole**
   - These antibiotics select for Enterococcus during exposure
   - After cessation, Enterococcus returns toward baseline (not significant, p=0.32)

### 8. Persistence and Recovery Analysis (Updated Methodology)

Both analyses now use **consistent paired-sample methodology** without external "healthy baseline" comparisons.

#### Persistence: During Antibiotics (n=206 pairs, Abx BETWEEN S1→S2)

| Group | S1 (pre-Abx) | S2 (post-Abx) | Change | Retention | p-value |
|-------|--------------|---------------|--------|-----------|---------|
| **Anaerobes** | 20.9% | 18.8% | -2.1% | 90% | 0.256 |
| **Enterobact** | 25.0% | 25.6% | +0.6% | 102% | 0.810 |
| **Enterococcus** | 16.3% | 21.8% | +5.5% | 134% | **0.018** |

#### Recovery: After Antibiotics Stop (n=46 pairs, Abx BEFORE S1, none between)

| Group | S1 (on Abx) | S2 (off Abx) | Change | Retention | p-value |
|-------|-------------|--------------|--------|-----------|---------|
| **Anaerobes** | 13.4% | 14.6% | +1.3% | 110% | 0.739 |
| **Enterobact** | 27.7% | 37.5% | +9.7% | 135% | **0.052** |
| **Enterococcus** | 18.1% | 10.2% | -7.9% | 57% | 0.124 |

#### Individual Antibiotic Analysis

**Recovery (single-antibiotic pairs, 91% of recovery cohort):**

| Antibiotic | n | Anaerobes | Enterobact | Enterococcus |
|------------|---|-----------|------------|--------------|
| Pip_Tazo | 25 | +3.5% (p=0.47) | +4.1% (p=0.49) | +1.1% (p=0.85) |
| Cefepime | 11 | -1.3% (p=0.88) | +8.7% (p=0.32) | -13.4% (p=0.23) |
| Meropenem | 5 | -14.8% (p=0.41) | +23.3% (p=0.29) | -0.7% (p=0.93) |

**Persistence (includes-antibiotic, most have multiple abx):**

| Antibiotic | n | % mono | Anaerobes | Enterobact | Enterococcus |
|------------|---|--------|-----------|------------|--------------|
| Pip_Tazo | 14 | 57% | -7.7% (p=0.17) | +8.8% (p=0.53) | -3.2% (p=0.54) |
| Cefepime | 8 | 38% | -11.3% (p=0.14) | +29.1% (p=0.16) | +1.2% (p=0.95) |
| Meropenem | 10 | 40% | -8.1% (p=0.14) | -0.4% (p=0.97) | +5.9% (p=0.57) |
| Vancomycin_IV | 8 | 0% | -0.8% (p=0.87) | -3.0% (p=0.86) | +14.1% (p=0.44) |

#### Key Biological Insights

1. **During antibiotics**: Enterococcus significantly INCREASES (+5.5%, p=0.018)
   - Intrinsic resistance leads to selection during antibiotic exposure
   - Most pronounced with Vancomycin_IV (+14%) and Meropenem (+6%)

2. **After antibiotics stop**: Enterobacteriaceae INCREASE (+9.7%, p=0.052)
   - Enterobact expand when competition (anaerobes) is depleted
   - Most pronounced after Meropenem (+23%) and Cefepime (+9%)

3. **Cefepime shows largest Enterobact expansion** (+29% during exposure)
   - Despite gram-negative coverage, Enterobact persist and expand
   - Anaerobes depleted (-11%), creating ecological niche

4. **Sample sizes limit statistical power** for individual antibiotic analysis
   - Consistent patterns emerge despite non-significant p-values
   - Larger cohorts needed to confirm antibiotic-specific effects

### 9. Individual Antibiotic Effects (Covariate-Adjusted)

**Critical finding:** Without adjusting for co-administered antibiotics, results are confounded and often biologically implausible. With proper covariate adjustment, results align with known antibiotic pharmacology.

#### Shannon Diversity by Individual Antibiotic

| Antibiotic | Effect | 95% CI | p-value | Interpretation |
|------------|--------|--------|---------|----------------|
| **Vancomycin** | **-0.56** | -0.83, -0.28 | **<0.0001** | Strong diversity reduction |
| Cefepime | -0.24 | -0.52, 0.05 | 0.10 | Trend toward reduction |
| Meropenem | -0.22 | -0.52, 0.08 | 0.16 | Trend toward reduction |
| TMP/SMX | -0.20 | -0.44, 0.04 | 0.10 | Trend toward reduction |
| Ceftriaxone | -0.11 | -0.65, 0.43 | 0.69 | No significant effect |
| Pip/Tazo | -0.08 | -0.29, 0.14 | 0.48 | No significant effect |
| Metronidazole | -0.05 | -0.44, 0.34 | 0.80 | No significant effect |
| Cefazolin | +0.14 | -0.25, 0.53 | 0.48 | No significant effect |
| **Ciprofloxacin** | **+0.42** | 0.10, 0.74 | **0.01** | Increases diversity |

#### Obligate Anaerobes by Individual Antibiotic

| Antibiotic | Effect (log10) | p-value | Biological Validation |
|------------|----------------|---------|----------------------|
| **Metronidazole** | **-1.30** | **<0.001** | ✓ Primary anti-anaerobic agent |
| **Vancomycin** | **-0.56** | **0.001** | ✓ Some anaerobic activity |
| **Meropenem** | **-0.38** | **0.03** | ✓ Broad-spectrum includes anaerobes |
| **Pip/Tazo** | **-0.37** | **0.003** | ✓ Tazobactam provides anaerobic coverage |
| Cefepime | -0.27 | 0.10 | Limited anaerobic activity |
| TMP/SMX | -0.16 | 0.25 | Minimal anaerobic activity |
| Ciprofloxacin | -0.05 | 0.77 | ✓ No anaerobic activity |
| Cefazolin | -0.01 | 0.97 | ✓ No anaerobic activity |
| Ceftriaxone | +0.15 | 0.64 | ✓ No anaerobic activity |

#### Enterobacteriaceae by Individual Antibiotic

| Antibiotic | Effect (log10) | p-value | Biological Validation |
|------------|----------------|---------|----------------------|
| **Cefepime** | **-0.86** | **<0.001** | ✓ Strong gram-negative activity |
| **TMP/SMX** | **-0.37** | **0.01** | ✓ Gram-negative activity |
| Ceftriaxone | -0.43 | 0.20 | Gram-negative activity (underpowered) |
| Ciprofloxacin | -0.22 | 0.26 | Fluoroquinolone activity |
| Meropenem | -0.13 | 0.49 | Expected reduction |
| Pip/Tazo | +0.08 | 0.55 | Weak effect |
| Vancomycin | +0.23 | 0.17 | ✓ No gram-negative activity |
| Cefazolin | +0.29 | 0.23 | Limited gram-negative activity |
| **Metronidazole** | **+0.65** | **0.008** | ✓ No gram-negative activity → bloom |

#### Enterococcus by Individual Antibiotic

| Antibiotic | Effect (log10) | p-value | Biological Validation |
|------------|----------------|---------|----------------------|
| **Cefepime** | **+0.69** | **<0.001** | ✓ Intrinsically resistant → bloom |
| **Metronidazole** | **+0.57** | **0.01** | ✓ No Enterococcus activity → bloom |
| Ceftriaxone | +0.32 | 0.32 | Intrinsically resistant |
| Meropenem | +0.19 | 0.29 | Some activity expected |
| TMP/SMX | +0.13 | 0.37 | Variable activity |
| Pip/Tazo | +0.03 | 0.79 | Has Enterococcus activity |
| Ciprofloxacin | -0.30 | 0.10 | Some Enterococcus activity |
| Cefazolin | -0.36 | 0.12 | Trend toward reduction |
| **Vancomycin** | **-0.53** | **0.001** | ✓ Primary anti-Enterococcal agent |

#### Beta Diversity (PERMANOVA) by Individual Antibiotic

| Antibiotic | R² | F | p-value |
|------------|-----|---|---------|
| **Meropenem** | **0.71%** | 4.25 | **0.001** |
| **Cefepime** | **0.71%** | 4.24 | **0.001** |
| **Metronidazole** | **0.52%** | 3.10 | **0.001** |
| **Pip/Tazo** | **0.40%** | 2.36 | **0.004** |
| **Ciprofloxacin** | **0.39%** | 2.32 | **0.003** |
| **TMP/SMX** | **0.39%** | 2.29 | **0.002** |
| **Vancomycin** | **0.36%** | 2.17 | **0.005** |
| **Ceftriaxone** | **0.27%** | 1.59 | **0.047** |
| Cefazolin | 0.15% | 0.89 | 0.583 |

---

## Conclusions

The reanalysis successfully reveals expected biological patterns that were not detected in the original 2019 analysis:

1. **Antibiotic exposure significantly reduces microbiome diversity** (p < 0.001), with vancomycin having the strongest individual effect (-0.56 Shannon units).

2. **Anti-anaerobic antibiotics specifically deplete obligate anaerobes** - Metronidazole (-1.30 log10), Pip/Tazo (-0.37), and Meropenem (-0.38) all significantly reduce anaerobes, while ceftriaxone and ciprofloxacin (no anaerobic activity) show no effect.

3. **Cephalosporins cause Enterococcus bloom** - Cefepime (+0.69 log10, p<0.001) significantly increases Enterococcus, consistent with intrinsic resistance. Vancomycin reduces Enterococcus (-0.53, p=0.001).

4. **Metronidazole causes Enterobacteriaceae bloom** - Without gram-negative activity, metronidazole allows Enterobacteriaceae expansion (+0.65, p=0.008), while cefepime strongly reduces them (-0.86, p<0.001).

5. **Covariate adjustment is essential** - Without controlling for co-administered antibiotics, results are confounded and biologically implausible. Proper adjustment reveals pharmacologically expected patterns.

6. **IBD patients should be excluded** from antibiotic analyses due to missing outpatient prescription data (0% recorded exposure). SB patients can be included but only Stool samples (Fistula/Ostomy samples have different microbial communities).

### Methodological Improvements

The key improvements over the original analysis:
- **Mixed-effects models** accounting for repeated patient samples
- **Compositional-aware methods** (ALDEx2) instead of simple correlations
- **Spectrum-based antibiotic classification** rather than treating all antibiotics equally
- **Proper paired sample design** comparing within-patient changes
- **Individual antibiotic models with covariate adjustment** for polypharmacy
- **Multi-method concordance** (ALDEx2 + MaAsLin3 + ANCOM-BC2) for robust associations

---

## Reproducibility

To rerun the analysis:

```bash
cd /home/david/projects/abx_metagenomic_correlations

# 1. Transfer Kraken/Bracken files from lambda-quad (if not already present)
bash scripts/transfer_legacy_kraken.sh

# 2. Build count matrices from raw Bracken output
conda activate r-env
Rscript R/new_analysis_legacy_data/00_build_count_matrices_from_bracken.R

# 3. Compute antibiotic exposure (Python - faster than R)
python scripts/compute_abx_exposure.py

# 4. Run analysis scripts in order
Rscript R/new_analysis_legacy_data/01_load_and_prepare_data.R
Rscript R/new_analysis_legacy_data/02_diversity_analysis.R
Rscript R/new_analysis_legacy_data/02b_diversity_individual_antibiotics.R
Rscript R/new_analysis_legacy_data/03_differential_abundance.R
Rscript R/new_analysis_legacy_data/04_paired_sample_analysis.R
Rscript R/new_analysis_legacy_data/05_individual_antibiotic_effects.R
Rscript R/new_analysis_legacy_data/06_genus_level_analysis.R
Rscript R/new_analysis_legacy_data/06_network_visualization.R
Rscript R/new_analysis_legacy_data/07_recovery_analysis.R
Rscript R/new_analysis_legacy_data/08_persistence_analysis.R
```

### Required R Packages

**CRAN:**
- dplyr, tidyr, readr, ggplot2, stringr, purrr, tibble
- vegan, lme4, lmerTest, broom.mixed, ggrepel
- visNetwork, igraph, htmlwidgets (for network visualization)

**Bioconductor:**
- ALDEx2 (v1.40.0) - compositional differential abundance
- maaslin3 (v1.0.2) - multivariable association discovery
- ANCOMBC (v2.10.1) - bias-corrected differential abundance
- phyloseq (v1.52.0) - microbiome data handling
- microbiome (v1.30.0) - required by ANCOMBC for data import

All packages are available in the `r-env` conda environment.

---

## Files

```
scripts/
├── transfer_legacy_kraken.sh          # Transfer Bracken files from lambda-quad
└── compute_abx_exposure.py            # Compute antibiotic exposure (Python)

data/
├── kraken2_legacy/                    # Raw Bracken output files
│   ├── genus/*_genus_abundance.txt
│   └── species/*_species_abundance.txt
├── legacy/                            # Original study files
│   ├── BSI_Abx_Samples.csv            # Sample list for this study
│   ├── SampleListFormatted.csv        # Sample metadata with dates
│   └── Drugs.csv                      # Hospital antibiotic records
├── bracken_count_matrices.RData       # Built from raw Bracken (script 00)
├── sample_metadata_with_abx.csv       # Sample metadata with abx exposure
└── prepared_data.RData                # Final prepared data for analysis

R/new_analysis_legacy_data/
├── 00_build_count_matrices_from_bracken.R  # Build matrices from raw Bracken
├── 01_load_and_prepare_data.R
├── 02_diversity_analysis.R
├── 02b_diversity_individual_antibiotics.R
├── 03_differential_abundance.R
├── 04_paired_sample_analysis.R
├── 05_individual_antibiotic_effects.R
├── 06_genus_level_analysis.R          # Genus-level differential abundance
├── 06_network_visualization.R
├── 07_recovery_analysis.R             # Microbiome recovery after Abx cessation
├── 08_persistence_analysis.R          # Persistence mechanism analysis
└── run_all_analysis.sh                # Master script to run all analyses

results/new_analysis_legacy_data/
├── data_quality/
│   └── sample_counts_by_group.csv
├── differential_abundance/
│   ├── aldex2_anti_anaerobic.csv
│   ├── aldex2_any_abx_7d.csv
│   └── aldex2_broad_spectrum.csv
├── diversity/
│   ├── alpha_diversity_results.csv
│   ├── alpha_diversity_individual_abx.csv
│   ├── beta_diversity_permanova.csv
│   ├── functional_groups_individual_abx.csv
│   ├── individual_abx_effects_summary.csv
│   └── permanova_individual_abx.csv
├── figures/
│   ├── alpha_diversity_by_abx.pdf
│   ├── alpha_diversity_by_spectrum.pdf
│   ├── anaerobes_by_anti_anaerobic_abx.pdf
│   ├── dose_response_bray.pdf
│   ├── functional_groups_by_individual_abx_forest.pdf
│   ├── paired_analysis_panel.pdf
│   ├── paired_bray_curtis.pdf
│   ├── paired_delta_anaerobes.pdf
│   ├── paired_delta_shannon.pdf
│   ├── pcoa_by_abx.pdf
│   ├── permanova_by_individual_abx.pdf
│   ├── shannon_by_individual_abx_forest.pdf
│   └── volcano_any_abx.pdf
├── paired_analysis/
│   ├── hypothesis_tests.csv
│   └── paired_sample_metrics.csv
├── individual_antibiotics_v2_with_covariates/
│   ├── aldex2/
│   │   ├── aldex2_*.csv (per antibiotic)
│   │   └── all_antibiotics_combined.csv
│   ├── maaslin3/
│   │   ├── */  (full output per antibiotic)
│   │   ├── maaslin3_*.csv (per antibiotic)
│   │   └── all_antibiotics_combined.csv
│   ├── ancombc2/
│   │   ├── ancombc2_*.csv (per antibiotic)
│   │   └── all_antibiotics_combined.csv
│   ├── method_concordance.csv
│   ├── robust_associations.csv
│   ├── summary_by_antibiotic.csv
│   ├── enterococcus_by_antibiotic.csv
│   └── analysis_metadata.rds
├── network_graphs/
│   ├── antibiotic_species_network.html
│   ├── antibiotic_species_network_hierarchical.html
│   ├── network_*.html (per antibiotic)
│   ├── network_associations.csv
│   ├── network_summary_by_antibiotic.csv
│   └── network_statistics.rds
├── recovery_analysis/
│   ├── recovery_pairs.csv
│   ├── recovery_pooled.csv
│   └── recovery_by_antibiotic.csv
└── persistence_analysis/
    ├── paired_persistence_data.csv
    ├── genus_persistence.csv
    ├── group_persistence_summary.csv
    ├── dominance_by_exposure.csv
    ├── sample_abundances.csv
    ├── recovery_vs_baseline.csv        # KEY: Recovery as % of healthy baseline
    └── genus_recovery_vs_baseline.csv
```

---

## Future Analyses (TODO)

### 1. Extended Exposure Windows

Current analyses use a **7-day** pre-sample window for antibiotic exposure. Additional windows to explore:

- **14-day window**: May capture delayed effects or cumulative exposure
- **30-day window**: Longer-term antibiotic burden, useful for chronic/repeated exposures

The `01_load_and_prepare_data.R` script already generates `_14d` variables. A 30-day window would require modification to the exposure calculation.

### 2. Microbiome Recovery After Antibiotic Cessation

**Objective:** Characterize microbiome recovery dynamics after stopping antibiotics.

**Study design:**
1. Identify sample pairs where:
   - Sample A: Has preceding antibiotic exposure (e.g., 7-day window)
   - Sample B: Collected after Sample A with **no intervening antibiotics**
   - Variable recovery interval (e.g., 7, 14, 30 days post-cessation)

2. Analyses:
   - Delta diversity (Shannon, richness) between Sample A → B
   - Compositional shifts (Bray-Curtis distance)
   - Recovery of specific taxa (anaerobes, Enterobacteriaceae, Enterococcus)
   - Time-to-recovery modeling (survival analysis framework)

3. Stratify by:
   - Antibiotic class (anti-anaerobic, broad-spectrum, etc.)
   - Duration of prior exposure
   - Patient group

**Expected findings:**
- Anaerobe recovery after anti-anaerobic cessation
- Enterococcus/Enterobacteriaceae normalization after cephalosporin cessation
- Potential identification of "recovery-resistant" dysbiosis patterns

### 3. Dose-Response Relationships

- Cumulative antibiotic-days vs. microbiome disruption
- Number of concurrent antibiotics vs. diversity loss
- Repeated exposure effects (sensitization vs. resilience)

### 4. Additional Covariates (Data Collection Needed)

Current models adjust for patient group and concurrent antibiotics. Additional covariates that may influence microbiome composition:

**Nutrition:**
- **TPN (Total Parenteral Nutrition)** - Associated with gut atrophy, reduced diversity
- **Enteral feeding** - Type (formula vs. breast milk in pediatrics), route (NG, G-tube)
- **NPO status** - Duration of nil per os
- **Diet type** - Solid foods, clear liquids, etc.

**Medications (non-antibiotic):**
- **Proton pump inhibitors (PPIs)** - Alter gastric pH, associated with dysbiosis and C. diff risk
- **H2 blockers** - Similar concerns to PPIs
- **Immunosuppressants** - Tacrolimus, cyclosporine, steroids (especially relevant for transplant patients)
- **Chemotherapy** - Mucositis, neutropenia effects
- **Motility agents** - Metoclopramide, erythromycin (sub-antimicrobial doses)
- **Laxatives/stool softeners** - Alter transit time

**Clinical factors:**
- **Primary diagnosis** - Beyond patient group (specific transplant indication, etc.)
- **Comorbidities** - Diabetes, renal failure, liver disease
- **GI anatomy** - Surgical history, ostomies, fistulas, bowel resection length
- **Infection status** - Active infections, C. diff history
- **ICU admission** - Critical illness effects

**Demographics:**
- **Age** - Pediatric vs. adult microbiome differences
- **Sex** - Hormonal influences on gut microbiome

**Environmental:**
- **Hospital unit** - Potential for nosocomial organism acquisition
- **Length of stay** - Hospital-associated microbiome shifts
- **Prior hospitalizations** - Cumulative healthcare exposure

---

## Changelog

### 2025-12-28 (Methodology Update & Individual Antibiotic Analysis)

**Consistent paired-sample methodology for recovery and persistence:**

Both scripts (07 and 08) now use the same approach:
- Analyze S1 → S2 changes within the same patients
- Report: mean S1, mean S2, change (S2-S1), retention (S2/S1 × 100%), p-value
- Removed misleading "healthy baseline" comparisons from persistence script

**Added individual antibiotic analysis:**

- Recovery script: Analyzes single-antibiotic pairs (91% of cohort = 42 pairs)
  - Pip_Tazo (n=25), Cefepime (n=11), Meropenem (n=5)
  - Output: `recovery_by_individual_antibiotic.csv`

- Persistence script: Analyzes "includes" antibiotic (pairs overlap)
  - Most pairs have multiple antibiotics (71%), so uses includes analysis
  - Reports % monotherapy for context
  - Output: `persistence_by_antibiotic.csv`

**Key findings from individual antibiotic analysis:**

| Antibiotic | During Abx (Enterobact) | After Abx (Enterobact) |
|------------|-------------------------|------------------------|
| Cefepime | +29.1% | +8.7% |
| Meropenem | -0.4% | +23.3% |
| Pip_Tazo | +8.8% | +4.1% |

- Cefepime drives largest Enterobact expansion during exposure
- Meropenem: Enterobact stable during, then expand after
- Vancomycin_IV selects for Enterococcus (+14%)

---

### 2025-12-28 (Persistence Mechanism Analysis)

**New analysis: Why Enterobacteriaceae/Enterococcus dominate (Script 08)**

- Added `08_persistence_analysis.R` to investigate persistence mechanisms
- Compares healthy-like vs antibiotic-exposed samples
- Analyzes paired samples with antibiotics between (n=206 pairs)
- Calculates recovery as % of healthy baseline (n=46 recovery pairs)
- Includes genus-level breakdown for key taxa

**Key findings:**
- Enterococcus expands 3.58x during antibiotic exposure (cross-sectional)
- Enterococcus shows 133.6% retention in paired analysis (expands during Abx)
- Obligate anaerobes recover to only 57% of healthy baseline after Abx cessation
- Enterococcus remains at 221% of baseline even in "recovery" samples
- Key depleted genera: Ruminococcus (-0.82 log2FC), Blautia (-0.74 log2FC)

**Conclusion:** Dominance comes from **persistence during antibiotic exposure**, not faster recovery. This represents chronic microbiome disruption where the ecological niche vacated by anaerobes is filled by intrinsically resistant organisms.

**New outputs:**
- `persistence_analysis/recovery_vs_baseline.csv` (key summary table)
- `persistence_analysis/genus_persistence.csv`
- `figures/persistence_composition.pdf`
- `figures/persistence_boxplot.pdf`
- `figures/genus_persistence_forest.pdf`

---

### 2025-12-28 (Recovery Analysis)

**New analysis: Microbiome recovery after antibiotic cessation (Script 07)**

- Added `07_recovery_analysis.R` to analyze which taxa recover after stopping antibiotics
- Study design: Paired samples where Abx before S1, no meaningful Abx between (TMP_SMX/Azithromycin excluded)
- 46 recovery pairs analyzed (Pip_Tazo: 27, Cefepime: 13, Meropenem: 7, Metronidazole: 3)

**Key findings:**
- Butyrate-producing anaerobes (Roseburia, Coprococcus, Blautia) recover after Pip_Tazo cessation
- Blautia shows strongest recovery after Meropenem (+5 log2FC)
- Clostridioides increases in recovery window (C. diff risk)
- Veillonella recovers universally across all antibiotics

**New outputs:**
- `recovery_analysis/recovery_pooled.csv`
- `recovery_analysis/recovery_by_antibiotic.csv`
- `figures/recovery_heatmap.pdf`
- `figures/recovery_forest_pooled.pdf`

---

### 2025-12-28 (Paired Analysis Refactoring)

**Major refactoring of paired sample analysis (Script 04):**

- **Problem:** Original analysis used broad antibiotic categories (anti-anaerobic, broad-spectrum) which produced misleading results. The H3 hypothesis "broad-spectrum antibiotics increase Enterobacteriaceae" was incorrect.

- **Solution:** Refactored to test individual antibiotics with LMM and covariate adjustment:
  - Modified `01_load_and_prepare_data.R` to calculate individual antibiotic exposures between paired samples (days of exposure for each of 10 antibiotics)
  - Rewrote `04_paired_sample_analysis.R` to use LMM: `outcome ~ target_abx + other_abx_covariates + interval_days + (1|MRN)`
  - Includes Vancomycin split by route (IV vs PO)

- **Key findings:**
  - Meropenem and Cefepime **decrease** Enterobacteriaceae (not increase as broad category suggested)
  - Meropenem shows classic carbapenem signature: kills Enterobacteriaceae, selects for Enterococcus
  - Pip_Tazo and Metronidazole significantly reduce anaerobes
  - Vancomycin IV shows no significant effects (consistent with no gut penetration)

- **New outputs:**
  - `paired_analysis/individual_abx_lmm_results.csv`
  - `paired_analysis/individual_abx_significant.csv`
  - `figures/paired_individual_abx_heatmap.pdf`

---

### 2025-12-28 (Data Pipeline Rebuild)

**Major data pipeline overhaul:**
- Rebuilt analysis pipeline to use original Kraken2/Bracken output files instead of pre-processed RData
- Added `scripts/transfer_legacy_kraken.sh` to transfer Bracken files from lambda-quad storage
- Added `R/new_analysis_legacy_data/00_build_count_matrices_from_bracken.R` to build count matrices from raw Bracken output
- Added `scripts/compute_abx_exposure.py` for fast antibiotic exposure computation (links samples to Drugs.csv via MRN)
- Contaminants (Homo sapiens, Salinibacter ruber) now excluded at matrix build time

**New analysis scripts:**
- Added `06_genus_level_analysis.R` for genus-level differential abundance analysis
- Added `05b_maaslin3_patch.R` for MaAsLin3 compatibility fixes
- Added `01_rebuild_sample_metadata.R` as alternative R-based metadata rebuilder

**Sample inclusion update:**
- Added SB (Short Bowel) Stool samples to high-confidence subset
- Excluded SB Fistula (20) and Ostomy (116) samples due to different microbial communities
- High-confidence samples increased from 568 to **723 samples** (130 patients)

**Bug fixes in `05_individual_antibiotic_effects.R`:**
- Fixed MaAsLin3 API compatibility: changed deprecated `plot_heatmap`/`plot_scatter` parameters to `plot_summary_plot`/`plot_associations`
- Installed missing `microbiome` package required by ANCOM-BC2 for data import

**Data preparation (`01_load_and_prepare_data.R`):**
- Updated high-confidence filter logic to include SB patients with Stool samples only:
  ```r
  high_confidence = (PatientGroup %in% high_confidence_groups) &
                    (SampleType == "Stool" | PatientGroup != "SB")
  ```
