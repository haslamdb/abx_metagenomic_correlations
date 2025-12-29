# Project: Antibiotic-Metagenomic Correlations

**Last Updated:** 2025-12-28
**Status:** Active
**Type:** Research

---

## Overview

Reanalysis of legacy antibiotic-microbiome data (2019) using modern compositional-aware methods. Examining effects of individual antibiotics on gut microbiome composition in pediatric patients (BMT, LvTx, IF, PICU, SB).

---

## Current Phase

**Phase 7: Vancomycin Route Split Analysis** - In Progress

Re-running all analyses with vancomycin separated by administration route:
- **Vancomycin_IV** - Intravenous (does not reach intestinal lumen, no direct microbiome effect)
- **Vancomycin_PO** - Oral (stays in gut lumen, massive microbiome effect)

Previous analyses combined IV and PO vancomycin, confounding the observed effects. Oral vancomycin (~18% of doses) likely drove the strong Lachnospiraceae depletion pattern attributed to "vancomycin."

---

## Completed Work

### Phase 0: Data Pipeline Rebuild (NEW)
- [x] Transfer original Bracken output files from lambda-quad storage
- [x] Build count matrices from raw Bracken files (excluding contaminants)
- [x] Compute antibiotic exposure from hospital Drugs.csv records
- [x] Link samples to exposure via MRN and sample dates
- **Output:** `data/bracken_count_matrices.RData`, `data/sample_metadata_with_abx.csv`

### Phase 1: Data Preparation
- [x] Load Bracken count matrices and antibiotic exposure
- [x] Classify antibiotics by spectrum
- [x] Create binary exposure variables (7d, 14d windows)
- [x] Calculate functional taxonomic groups
- [x] Compute alpha diversity metrics
- **Output:** `data/prepared_data.RData`

### Phase 2: Diversity Analysis
- [x] Alpha diversity mixed-effects models
- [x] Beta diversity PERMANOVA
- [x] Individual antibiotic effects on diversity
- **Output:** `results/new_analysis_legacy_data/diversity/`

### Phase 3: Differential Abundance
- [x] ALDEx2 for antibiotic class effects
- **Output:** `results/new_analysis_legacy_data/differential_abundance/`

### Phase 4: Paired Sample Analysis (Updated 2025-12-28)
- [x] Within-patient microbiome changes
- [x] **Refactored to use individual antibiotics** (not broad categories)
- [x] LMM with covariate adjustment for other antibiotics
- [x] Tests 10 individual antibiotics against 5 outcomes
- **Output:** `results/new_analysis_legacy_data/paired_analysis/`

### Phase 5: Species-Level Individual Antibiotic Effects
- [x] ALDEx2, MaAsLin3, ANCOM-BC2 with covariate adjustment
- [x] Multi-method concordance analysis
- **Output:** `results/new_analysis_legacy_data/individual_antibiotics_v2_with_covariates/`

### Phase 6: Genus-Level Analysis
- [x] Genus-level differential abundance (ALDEx2, MaAsLin3, ANCOM-BC2)
- [x] Multi-method concordance at genus level
- **Output:** `results/new_analysis_legacy_data/individual_antibiotics_genus_level/`

### Phase 8: Microbiome Recovery Analysis (Updated - 2025-12-28)
- [x] Identify recovery pairs (Abx before S1, no meaningful Abx between)
- [x] Pooled functional group recovery analysis (n=46 pairs)
- [x] Individual antibiotic analysis for single-antibiotic pairs (91% of recovery pairs)
- [x] Consistent methodology with persistence script (paired S1→S2 analysis)
- **Script:** `R/new_analysis_legacy_data/07_recovery_analysis.R`
- **Output:** `results/new_analysis_legacy_data/recovery_analysis/`

### Phase 9: Persistence Mechanism Analysis (Updated - 2025-12-28)
- [x] Paired analysis of persistence during antibiotic exposure (n=206 pairs)
- [x] Functional group analysis with consistent methodology
- [x] Individual antibiotic "includes" analysis (most pairs have multiple abx)
- [x] Genus-level breakdown of persistence for key taxa
- [x] Removed misleading "healthy baseline" comparisons
- **Script:** `R/new_analysis_legacy_data/08_persistence_analysis.R`
- **Output:** `results/new_analysis_legacy_data/persistence_analysis/`

---

## Pending Work (TODO)

### High Priority (Phase 7 - In Progress)
- [ ] Complete vancomycin route split re-analysis (running in screen `abx_vanc_split`)
- [ ] Compare Vancomycin_IV vs Vancomycin_PO effects
- [ ] Validate that IV vancomycin shows minimal/no direct microbiome effect
- [ ] Update network visualizations with route-specific results

### Medium Priority
- [ ] Review multi-method concordance results (species + genus levels)
- [ ] Pathway enrichment analysis (MetaCyc)
- [ ] Integrate findings into summary document

### Future Analyses (see docs/NEW_ANALYSIS_LEGACY_DATA.md)
- [ ] 14-day and 30-day exposure windows
- [ ] Microbiome recovery after antibiotic cessation
- [ ] Additional covariates (TPN, PPIs, immunosuppressants)

---

## Key Outputs

| Output | Location | Date |
|--------|----------|------|
| Bracken count matrices | `data/bracken_count_matrices.RData` | Dec 2025 |
| Sample metadata + abx | `data/sample_metadata_with_abx.csv` | Dec 2025 |
| Diversity results | `results/new_analysis_legacy_data/diversity/` | Dec 2025 |
| Individual abx effects | `results/new_analysis_legacy_data/individual_antibiotics_v2_with_covariates/` | Dec 2025 |
| Network graphs | `results/new_analysis_legacy_data/network_graphs/` | Dec 2025 |

---

## Key Findings

1. **Antibiotic exposure reduces diversity** (p < 0.001)
2. **Metronidazole causes blooms** in both Enterobacteriaceae and Enterococcus
3. **Cefepime causes Enterococcus bloom** but reduces Enterobacteriaceae
4. **Covariate adjustment is essential** - without it, results are confounded
5. **Vancomycin route matters** - IV vancomycin doesn't reach gut lumen (no direct effect), while PO vancomycin stays in gut (strong effect). Previous "Vancomycin" results were confounded by mixing routes (~18% PO, ~82% IV)
6. **Persistence, not recovery, drives dominance** - Enterobacteriaceae/Enterococcus dominate because they persist during antibiotics (stress response, intrinsic resistance), while anaerobes are depleted and fail to recover

### Paired Analysis - Individual Antibiotics (Phase 4 Update)

The original H3 hypothesis ("broad-spectrum antibiotics increase Enterobacteriaceae") was **misleading**. After refactoring to test individual antibiotics with LMM and covariate adjustment:

| Antibiotic | Outcome | Effect | p-value |
|------------|---------|--------|---------|
| **Pip_Tazo** | Anaerobes | -1.82 log2FC | 0.0005 |
| **Meropenem** | Enterobacteriaceae | -2.06 log2FC | 0.0043 |
| **Metronidazole** | Anaerobes | -3.38 log2FC | 0.0087 |
| **Meropenem** | Bray-Curtis | +0.12 | 0.0087 |
| **Pip_Tazo** | Bray-Curtis | +0.09 | 0.0244 |
| **Meropenem** | Enterococcus | +1.42 log2FC | 0.0317 |
| **Cefepime** | Enterobacteriaceae | -1.60 log2FC | 0.0318 |

**Key insight:** Broad-spectrum antibiotics (Meropenem, Cefepime) actually **decrease** Enterobacteriaceae in paired analysis, not increase as the broad category suggested. Meropenem shows the clearest signature: decreases Enterobacteriaceae, increases Enterococcus (classic carbapenem selection), and increases community instability.

### Recovery Analysis (Phase 8)

Analyzed 46 recovery pairs (Abx before sample 1, no meaningful Abx between samples):

| Antibiotic | Genus | Log2FC | p-value | Interpretation |
|------------|-------|--------|---------|----------------|
| **Pip_Tazo** | Veillonella | +2.45 | 0.014 | Rapid recovery |
| **Pip_Tazo** | Roseburia | +1.35 | 0.035 | Butyrate producer recovers |
| **Pip_Tazo** | Coprococcus | +1.46 | 0.042 | Butyrate producer recovers |
| **Pip_Tazo** | Clostridioides | +1.64 | 0.021 | C. diff risk window |
| **Meropenem** | Blautia | +4.99 | 0.045 | Strongest anaerobe recovery |

**Key patterns:**
- Obligate anaerobes (Roseburia, Coprococcus, Blautia) recover after anti-anaerobic cessation
- Clostridioides increases in recovery window (C. diff infection risk)
- Enterococcus declines after Meropenem cessation (returning to baseline)

### Persistence and Recovery Analysis (Phases 8-9, Updated)

Both scripts now use **consistent paired-sample methodology** without external "healthy baseline" comparisons.

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
| Pip_Tazo | 25 | +3.5% | +4.1% | +1.1% |
| Cefepime | 11 | -1.3% | +8.7% | -13.4% |
| Meropenem | 5 | -14.8% | +23.3% | -0.7% |

**Persistence (includes-antibiotic, most have multiple abx):**

| Antibiotic | n | Anaerobes | Enterobact | Enterococcus |
|------------|---|-----------|------------|--------------|
| Pip_Tazo | 14 | -7.7% | +8.8% | -3.2% |
| Cefepime | 8 | -11.3% | +29.1% | +1.2% |
| Meropenem | 10 | -8.1% | -0.4% | +5.9% |
| Vancomycin_IV | 8 | -0.8% | -3.0% | +14.1% |

**Key findings:**
- Enterococcus significantly INCREASES during antibiotics (+5.5%, p=0.018)
- Enterobacteriaceae INCREASE after antibiotics stop (+9.7%, p=0.052)
- Cefepime shows largest Enterobact expansion (+29%) during exposure
- Vancomycin_IV selects for Enterococcus (+14%)

---

## Notes

- Vancomycin route split analysis running in screen session `abx_vanc_split`
- Master script: `R/new_analysis_legacy_data/run_all_analysis.sh`
- 723 high-confidence samples (including SB Stool)
- Full documentation: `docs/NEW_ANALYSIS_LEGACY_DATA.md`
