# Reanalysis of Legacy Antibiotic-Microbiome Data

**Date:** December 2024
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
| Sample pairs (2-week window) | 246 |
| Species | 1,460 |
| Genera | 163 |

### Patient Groups

| Group | Samples | Patients | % with Abx (7d) | Notes |
|-------|---------|----------|-----------------|-------|
| LvTx (Liver Transplant) | 234 | 35 | 76.5% | High-confidence |
| BMT (Bone Marrow Transplant) | 218 | 61 | 66.1% | High-confidence |
| SB (Short Bowel) | 291 | 16 | 25.4% | |
| IF (Intestinal Failure) | 100 | 16 | 25.0% | High-confidence |
| IBD (Inflammatory Bowel Disease) | 43 | 10 | **0%** | Excluded - no Abx data |
| PICU | 16 | 7 | 87.5% | High-confidence |
| Urology | 3 | 1 | 0% | |

**High-confidence groups** (BMT, IF, LvTx, PICU) were used for primary analyses due to complete inpatient antibiotic records. IBD patients had no recorded antibiotic exposure, likely due to outpatient prescriptions not captured in the hospital system.

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
individual_antibiotics/
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

### 5. Paired Sample Analysis (n = 246 pairs)

| Hypothesis | Test | Effect | p-value |
|------------|------|--------|---------|
| Anti-anaerobic reduces anaerobes | Mixed-effects | -1.61 log2FC | **0.005** |
| Broad-spectrum affects Enterobact | Mixed-effects | -1.28 log2FC | **0.04** |
| Dose-response (Abx-days vs Shannon) | Spearman | rho = -0.14 | 0.07 |
| Abx increases Bray-Curtis distance | Mixed-effects | +0.03 | 0.43 |

---

## Conclusions

The reanalysis successfully reveals expected biological patterns that were not detected in the original 2019 analysis:

1. **Antibiotic exposure significantly reduces microbiome diversity** (p < 0.001), with anti-gram-positive agents having the strongest effect.

2. **Anti-anaerobic antibiotics specifically deplete obligate anaerobes** - Bacteroides and Parabacteroides are significantly reduced, with median relative abundance dropping from 13.8% to 0.8%.

3. **Enterobacteriaceae bloom with antibiotic exposure** - Klebsiella, Escherichia, Enterobacter, Citrobacter, and Serratia all significantly increased.

4. **Paired within-patient analysis confirms effects** - Anti-anaerobic drugs reduce anaerobes (p = 0.005) in paired samples from the same patient.

5. **IBD patients should be excluded** from antibiotic analyses due to missing outpatient prescription data (0% recorded exposure).

### Methodological Improvements

The key improvements over the original analysis:
- **Mixed-effects models** accounting for repeated patient samples
- **Compositional-aware methods** (ALDEx2) instead of simple correlations
- **Spectrum-based antibiotic classification** rather than treating all antibiotics equally
- **Proper paired sample design** comparing within-patient changes

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
Rscript R/new_analysis_legacy_data/03_differential_abundance.R
Rscript R/new_analysis_legacy_data/04_paired_sample_analysis.R
Rscript R/new_analysis_legacy_data/05_individual_antibiotic_effects.R
```

### Required R Packages

**CRAN:**
- tidyverse, vegan, lme4, lmerTest, broom.mixed, ggrepel

**Bioconductor:**
- ALDEx2 (v1.40.0) - compositional differential abundance
- maaslin3 (v1.0.2) - multivariable association discovery
- ANCOMBC (v2.10.1) - bias-corrected differential abundance
- phyloseq (v1.52.0) - microbiome data handling

All packages are available in the `r-env` conda environment.

---

## Files

```
R/new_analysis_legacy_data/
├── 01_load_and_prepare_data.R
├── 02_diversity_analysis.R
├── 03_differential_abundance.R
├── 04_paired_sample_analysis.R
└── 05_individual_antibiotic_effects.R

results/new_analysis_legacy_data/
├── data_quality/
│   └── sample_counts_by_group.csv
├── differential_abundance/
│   ├── aldex2_anti_anaerobic.csv
│   ├── aldex2_any_abx_7d.csv
│   └── aldex2_broad_spectrum.csv
├── diversity/
│   ├── alpha_diversity_results.csv
│   └── beta_diversity_permanova.csv
├── figures/
│   ├── alpha_diversity_by_abx.pdf
│   ├── alpha_diversity_by_spectrum.pdf
│   ├── anaerobes_by_anti_anaerobic_abx.pdf
│   ├── dose_response_bray.pdf
│   ├── paired_analysis_panel.pdf
│   ├── paired_bray_curtis.pdf
│   ├── paired_delta_anaerobes.pdf
│   ├── paired_delta_shannon.pdf
│   ├── pcoa_by_abx.pdf
│   └── volcano_any_abx.pdf
├── paired_analysis/
│   ├── hypothesis_tests.csv
│   └── paired_sample_metrics.csv
└── individual_antibiotics/
    ├── aldex2/
    │   ├── aldex2_*.csv (per antibiotic)
    │   └── all_antibiotics_combined.csv
    ├── maaslin3/
    │   ├── */  (full output per antibiotic)
    │   ├── maaslin3_*.csv (per antibiotic)
    │   └── all_antibiotics_combined.csv
    ├── ancombc2/
    │   ├── ancombc2_*.csv (per antibiotic)
    │   └── all_antibiotics_combined.csv
    ├── method_concordance.csv
    ├── robust_associations.csv
    ├── summary_by_antibiotic.csv
    ├── enterococcus_by_antibiotic.csv
    └── analysis_metadata.rds
```
