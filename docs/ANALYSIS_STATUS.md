# Project: Antibiotic-Metagenomic Correlations

**Last Updated:** 2025-12-28
**Status:** Active
**Type:** Research

---

## Overview

Reanalysis of legacy antibiotic-microbiome data (2019) using modern compositional-aware methods. Examining effects of individual antibiotics on gut microbiome composition in pediatric patients (BMT, LvTx, IF, PICU, SB).

---

## Current Phase

**Phase 5: Individual Antibiotic Effects** - In Progress

Running species-level differential abundance analysis for 10 antibiotics using ALDEx2, MaAsLin3, and ANCOM-BC2 with covariate adjustment.

---

## Completed Work

### Phase 1: Data Preparation
- [x] Load legacy RData, classify antibiotics by spectrum
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

### Phase 4: Paired Sample Analysis
- [x] Within-patient microbiome changes
- **Output:** `results/new_analysis_legacy_data/paired_analysis/`

---

## Pending Work (TODO)

### High Priority
- [ ] Complete 05_individual_antibiotic_effects.R run (currently running)
- [ ] Review multi-method concordance results
- [ ] Validate Enterococcus findings (unexpected Pip/Tazo association)

### Medium Priority
- [ ] Pathway enrichment analysis (MetaCyc)
- [ ] Network visualization updates
- [ ] Integrate findings into summary document

### Future Analyses (see docs/NEW_ANALYSIS_LEGACY_DATA.md)
- [ ] 14-day and 30-day exposure windows
- [ ] Microbiome recovery after antibiotic cessation
- [ ] Additional covariates (TPN, PPIs, immunosuppressants)

---

## Key Outputs

| Output | Location | Date |
|--------|----------|------|
| Diversity results | `results/new_analysis_legacy_data/diversity/` | Dec 2025 |
| Individual abx effects | `results/new_analysis_legacy_data/individual_antibiotics_v2_with_covariates/` | In progress |
| Network graphs | `results/new_analysis_legacy_data/network_graphs/` | Dec 2025 |

---

## Key Findings

1. **Antibiotic exposure reduces diversity** (p < 0.001), Vancomycin strongest effect
2. **Metronidazole causes blooms** in both Enterobacteriaceae and Enterococcus
3. **Cefepime causes Enterococcus bloom** but reduces Enterobacteriaceae
4. **Covariate adjustment is essential** - without it, results are confounded

---

## Notes

- Analysis running in screen session `abx_analysis`
- 723 high-confidence samples (including SB Stool)
- Full documentation: `docs/NEW_ANALYSIS_LEGACY_DATA.md`
