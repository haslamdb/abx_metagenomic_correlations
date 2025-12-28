# Antibiotic-Microbiome Correlations: Reanalysis Plan (v3)

## Overview

This plan focuses on re-analyzing existing Kraken/Bracken data using improved statistical methods. The previous analysis used simple Pearson correlations between rate-of-change metrics, which did not reveal expected biological patterns (e.g., anaerobe depletion with metronidazole).

**Data Source**: `data/legacy/AbxEffectData20191014` (RData file)

**Scope**: Microbiome composition and antibiotic-species correlations only (no ARG or functional pathway analysis in this phase)

---

## Phase 1: Data Loading and Validation

### 1.1 Load Existing Data Structures

Key objects expected from RData:
- `SpeciesTableNR` - Rarefied species abundance (~915 samples × ~1460 species)
- `GenusTableNR` - Rarefied genus abundance (~915 samples × ~303 genera)
- `DrugTable` - Antibiotic exposure data with dates
- `TwoWeekSamplePairs` / `ThreeMonthSamplePairs` - Paired sample metadata
- `DupStool` / `AvailableDupSamples` - Sample metadata

### 1.2 Data Quality Assessment

For each patient group, assess:
- Number of samples and unique patients
- Number of valid paired samples (2-week and 3-month windows)
- Antibiotic exposure data completeness
- Sequencing depth distribution

### 1.3 Patient Group Stratification

| Group | Description | Expected Abx Data Quality |
|-------|-------------|---------------------------|
| BMT | Bone marrow transplant | High (inpatient) |
| IF | Intestinal failure | High (inpatient) |
| LvTx | Liver transplant | High (inpatient) |
| PICU | Pediatric ICU | High (inpatient) |
| SB | Short bowel | Moderate |
| Urology | Urology patients | Moderate |
| IBD | Inflammatory bowel disease | Low (outpatient, may miss oral Abx) |

**Decision**: Run primary analyses on "high-confidence" groups (BMT, IF, LvTx, PICU), then validate on full dataset.

---

## Phase 2: Enhanced Antibiotic Exposure Variables

### 2.1 Binary Exposure Windows

For each sample, calculate:
```
abx_any_7d        - Any antibiotic in 7 days before sample
abx_any_14d       - Any antibiotic in 14 days before sample
abx_days_7d       - Number of antibiotic-days in prior 7 days
abx_days_14d      - Number of antibiotic-days in prior 14 days
```

### 2.2 Spectrum-Based Classification

| Category | Antibiotics |
|----------|-------------|
| Anti-anaerobic | Metronidazole, Piperacillin/Tazobactam, Meropenem, Ertapenem, Clindamycin, Amoxicillin/Clavulanate, Ampicillin/Sulbactam |
| Anti-gram-positive | Vancomycin, Daptomycin, Linezolid, Ceftaroline |
| Broad-spectrum | Meropenem, Piperacillin/Tazobactam, Cefepime, Imipenem |
| Anti-Pseudomonal | Cefepime, Piperacillin/Tazobactam, Meropenem, Tobramycin, Ciprofloxacin, Aztreonam |
| Narrow-spectrum | Penicillin, Ampicillin, Cephalexin, Nitrofurantoin |

Variables to create:
```
abx_anaerobic_7d    - Anti-anaerobic in prior 7 days
abx_gram_pos_7d     - Anti-gram-positive in prior 7 days
abx_broad_7d        - Broad-spectrum in prior 7 days
```

### 2.3 Specific Drug Exposures

For key antibiotics, create binary 7-day exposure flags:
- `pip_tazo_7d` - Piperacillin/Tazobactam
- `meropenem_7d` - Meropenem
- `cefepime_7d` - Cefepime
- `vancomycin_7d` - Vancomycin
- `metronidazole_7d` - Metronidazole

### 2.4 For Paired Samples

Calculate exposure between sample pairs:
```
abx_between_any     - Any antibiotic between samples
abx_between_days    - Total antibiotic-days between samples
abx_between_anaerobic - Anti-anaerobic between samples
abx_between_broad   - Broad-spectrum between samples
```

---

## Phase 3: Taxonomic Groupings

### 3.1 Functional Groups (Genus Level)

| Group | Genera | Biological Rationale |
|-------|--------|---------------------|
| Obligate anaerobes | Bacteroides, Parabacteroides, Prevotella, Clostridium, Faecalibacterium, Blautia, Roseburia, Ruminococcus, Alistipes, Odoribacter | Target of anti-anaerobic Abx |
| Enterobacteriaceae | Escherichia, Klebsiella, Enterobacter, Citrobacter, Serratia, Proteus | Often bloom after broad-spectrum Abx |
| Enterococcus | Enterococcus | Blooms with Abx pressure, VRE concern |
| Lactobacillales | Lactobacillus, Streptococcus, Leuconostoc, Pediococcus | Gram-positive fermenters |

### 3.2 Phylum-Level Summaries

- Total Firmicutes
- Total Bacteroidetes
- Total Proteobacteria
- Firmicutes/Bacteroidetes ratio

### 3.3 Diversity Metrics

- Shannon diversity
- Simpson diversity
- Observed richness (species count)
- Pielou's evenness

---

## Phase 4: Statistical Analyses

### 4.1 Cross-Sectional Analysis: Alpha Diversity

**Question**: Does recent antibiotic exposure reduce diversity?

```r
# Mixed-effects model accounting for repeated patients
library(lme4)
library(lmerTest)

model <- lmer(
  shannon ~ abx_any_7d + patient_group + (1|patient_id),
  data = sample_metadata
)

# Test specific antibiotic classes
model_spectrum <- lmer(
  shannon ~ abx_anaerobic_7d + abx_broad_7d + patient_group + (1|patient_id),
  data = sample_metadata
)
```

### 4.2 Cross-Sectional Analysis: Beta Diversity

**Question**: Does antibiotic exposure explain community composition differences?

```r
library(vegan)

# Bray-Curtis distance
bray_dist <- vegdist(species_matrix, method = "bray")

# PERMANOVA
adonis2(
  bray_dist ~ abx_any_7d + patient_group,
  data = metadata,
  permutations = 999,
  by = "margin"
)

# PCoA visualization colored by antibiotic exposure
pcoa <- cmdscale(bray_dist, k = 2, eig = TRUE)
```

### 4.3 Differential Abundance: ALDEx2

**Question**: Which taxa differ between Abx-exposed and unexposed samples?

```r
library(ALDEx2)

# Compositional-aware differential abundance
aldex_results <- aldex(
  reads = t(species_counts),  # taxa as rows, samples as columns
  conditions = factor(metadata$abx_any_7d),
  mc.samples = 128,
  test = "t",
  effect = TRUE
)

# Significant taxa (BH-corrected)
sig_taxa <- aldex_results[aldex_results$wi.eBH < 0.1, ]
```

### 4.4 Differential Abundance: MaAsLin2

**Question**: Which taxa are associated with antibiotic exposure after adjusting for covariates?

```r
library(Maaslin2)

maaslin_results <- Maaslin2(
  input_data = species_matrix,
  input_metadata = metadata,
  output = "maaslin2_output",
  fixed_effects = c("abx_anaerobic_7d", "abx_broad_7d", "patient_group"),
  random_effects = c("patient_id"),
  normalization = "TSS",
  transform = "LOG",
  min_abundance = 0.0001,
  min_prevalence = 0.1
)
```

### 4.5 Paired Sample Analysis (Primary Analysis)

**Question**: Within the same patient, how does the microbiome change with/without antibiotic exposure?

```r
# For each sample pair, calculate:
# 1. Microbiome change metrics
# 2. Whether antibiotics were given between samples

paired_analysis <- function(pairs_df, species_matrix, drug_data) {

  results <- pairs_df %>%
    rowwise() %>%
    mutate(
      # Diversity change
      delta_shannon = shannon[sample2] - shannon[sample1],

      # Taxonomic group changes (log2 fold change)
      delta_anaerobes = log2((anaerobe_sum[sample2] + 1) / (anaerobe_sum[sample1] + 1)),
      delta_enterobact = log2((enterobact_sum[sample2] + 1) / (enterobact_sum[sample1] + 1)),
      delta_enterococcus = log2((enterococcus[sample2] + 1) / (enterococcus[sample1] + 1)),

      # Community stability (Bray-Curtis)
      bray_distance = vegdist(rbind(species_matrix[sample1,], species_matrix[sample2,]))[1]
    )

  return(results)
}

# Statistical tests
# H1: Anti-anaerobic antibiotics decrease anaerobe abundance
wilcox.test(delta_anaerobes ~ abx_between_anaerobic, data = paired_results)

# H2: Antibiotics increase community disruption
wilcox.test(bray_distance ~ abx_between_any, data = paired_results)

# H3: Dose-response relationship
cor.test(paired_results$abx_between_days, paired_results$delta_shannon, method = "spearman")

# Mixed-effects for proper inference
model <- lmer(
  delta_anaerobes ~ abx_between_anaerobic + interval_days + (1|patient_id),
  data = paired_results
)
```

### 4.6 Patient Group Stratified Analysis

Run key analyses separately for:
1. **High-confidence groups**: BMT, IF, LvTx, PICU (complete Abx data expected)
2. **All groups**: Full dataset for comparison
3. **Sensitivity analysis**: Exclude IBD to assess impact of missing outpatient Abx

---

## Phase 5: Key Hypotheses

| # | Hypothesis | Analysis | Expected Direction |
|---|------------|----------|-------------------|
| H1 | Recent Abx exposure reduces alpha diversity | Mixed-effects model | Negative |
| H2 | Anti-anaerobic Abx reduce obligate anaerobe abundance | ALDEx2/paired analysis | Negative |
| H3 | Broad-spectrum Abx increase Enterobacteriaceae | ALDEx2/paired analysis | Positive |
| H4 | Broad-spectrum Abx increase Enterococcus | ALDEx2/paired analysis | Positive |
| H5 | Abx exposure increases within-patient instability | Paired Bray-Curtis | Positive |
| H6 | Metronidazole specifically depletes Bacteroides | MaAsLin2 | Negative |
| H7 | Vancomycin increases Enterococcus | MaAsLin2 | Positive |

---

## Phase 6: Visualizations

1. **Alpha diversity boxplots** - By Abx exposure status and patient group
2. **PCoA/NMDS ordination** - Colored by Abx exposure, shaped by patient group
3. **Volcano plots** - From ALDEx2/MaAsLin2 results
4. **Paired sample trajectories** - Before/after connected by lines, colored by Abx exposure
5. **Heatmap** - Top differentially abundant taxa by patient group
6. **Functional group barplots** - Anaerobes, Enterobacteriaceae proportions by Abx status

---

## Implementation Order

1. Load RData and assess data quality by patient group
2. Create enhanced antibiotic exposure variables
3. Calculate taxonomic groupings and diversity metrics
4. Run cross-sectional diversity analyses
5. Run differential abundance (ALDEx2, MaAsLin2)
6. Run paired sample analysis
7. Generate visualizations
8. Compile results and interpret

---

## Output Files

```
results/
├── data_quality/
│   ├── sample_counts_by_group.csv
│   └── abx_completeness_by_group.csv
├── diversity/
│   ├── alpha_diversity_results.csv
│   └── beta_diversity_permanova.csv
├── differential_abundance/
│   ├── aldex2_results.csv
│   └── maaslin2_results.csv
├── paired_analysis/
│   ├── paired_sample_metrics.csv
│   └── paired_hypothesis_tests.csv
└── figures/
    ├── alpha_diversity_boxplot.pdf
    ├── pcoa_by_abx.pdf
    ├── volcano_plot.pdf
    └── paired_trajectories.pdf
```
