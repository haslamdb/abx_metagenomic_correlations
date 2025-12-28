# Reanalysis of Legacy Antibiotic-Microbiome Data

**Date:** December 2025 (last updated: 2025-12-28)
**Data Source:** `data/legacy/AbxEffectData20191014` (RData from 2019 analysis)
**Scripts:** `R/new_analysis_legacy_data/`
**Results:** `results/new_analysis_legacy_data/`

---

## Overview

This analysis re-examines the relationship between antibiotic exposure and gut microbiome composition using improved statistical methods. The original 2019 analysis used simple Pearson correlations between rate-of-change metrics, which failed to reveal expected biological patterns. This reanalysis applies:

- Mixed-effects models to account for repeated patient samples
- Compositional-aware differential abundance testing (ALDEx2)
- Proper antibiotic spectrum classification
- Paired sample analysis for within-patient comparisons

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

### Script 1: Data Preparation
**`01_load_and_prepare_data.R`**

- Loads legacy RData file
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

### Script 4: Paired Sample Analysis
**`04_paired_sample_analysis.R`**

- Within-patient microbiome changes
- Hypothesis testing with Wilcoxon and mixed-effects models
- Dose-response analysis

**Outputs:**
- `paired_analysis/paired_sample_metrics.csv`
- `paired_analysis/hypothesis_tests.csv`
- `figures/paired_delta_shannon.pdf`
- `figures/paired_delta_anaerobes.pdf`
- `figures/paired_bray_curtis.pdf`
- `figures/paired_analysis_panel.pdf`
- `figures/dose_response_bray.pdf`

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

### 5. Paired Sample Analysis (n = 346 high-confidence pairs)

| Hypothesis | Test | Effect | p-value |
|------------|------|--------|---------|
| Anti-anaerobic reduces anaerobes | Mixed-effects | -1.61 log2FC | **0.005** |
| Broad-spectrum affects Enterobact | Mixed-effects | -1.28 log2FC | **0.04** |
| Dose-response (Abx-days vs Shannon) | Spearman | rho = -0.14 | 0.07 |
| Abx increases Bray-Curtis distance | Mixed-effects | +0.03 | 0.43 |

### 6. Individual Antibiotic Effects (Covariate-Adjusted)

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

# Activate R environment
conda activate r-env

# Run scripts in order
Rscript R/new_analysis_legacy_data/01_load_and_prepare_data.R
Rscript R/new_analysis_legacy_data/02_diversity_analysis.R
Rscript R/new_analysis_legacy_data/02b_diversity_individual_antibiotics.R
Rscript R/new_analysis_legacy_data/03_differential_abundance.R
Rscript R/new_analysis_legacy_data/04_paired_sample_analysis.R
Rscript R/new_analysis_legacy_data/05_individual_antibiotic_effects.R
Rscript R/new_analysis_legacy_data/06_network_visualization.R
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
R/new_analysis_legacy_data/
├── 01_load_and_prepare_data.R
├── 02_diversity_analysis.R
├── 02b_diversity_individual_antibiotics.R
├── 03_differential_abundance.R
├── 04_paired_sample_analysis.R
├── 05_individual_antibiotic_effects.R
└── 06_network_visualization.R

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
└── network_graphs/
    ├── antibiotic_species_network.html
    ├── antibiotic_species_network_hierarchical.html
    ├── network_*.html (per antibiotic)
    ├── network_associations.csv
    ├── network_summary_by_antibiotic.csv
    └── network_statistics.rds
```

---

## Changelog

### 2025-12-28

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
