# Project: Antibiotic-Metagenomic Correlations

**Last Updated:** 2025-12-29
**Status:** Active
**Type:** Research

---

## Overview

Reanalysis of legacy antibiotic-microbiome data (2019) using modern compositional-aware methods. Examining effects of individual antibiotics on gut microbiome composition in pediatric patients (BMT, LvTx, IF, PICU, SB).

---

## Current Phase

**Phase 7: Vancomycin Route Split Analysis** - **COMPLETE**

Re-ran all analyses with vancomycin separated by administration route:
- **Vancomycin_IV** - Intravenous (does not reach intestinal lumen, no direct microbiome effect) - **COMPLETE**
- **Vancomycin_PO** - Oral (stays in gut lumen, massive microbiome effect) - **SKIPPED** (only 8 exposed samples, need ≥10)

Previous analyses combined IV and PO vancomycin, confounding the observed effects. Oral vancomycin (~18% of doses) likely drove the strong Lachnospiraceae depletion pattern attributed to "vancomycin."

### Species-Level Analysis Status (Script 05)

| Antibiotic | ALDEx2 | MaAsLin3 | ANCOM-BC2 | Status |
|------------|--------|----------|-----------|--------|
| Pip_Tazo | ✓ | ✓ | ✓ | Complete |
| TMP_SMX | ✓ | ✓ | ✓ | Complete |
| Vancomycin_IV | ✓ | ✓ | ✓ | Complete |
| Cefepime | ✓ | ✓ | ✓ | Complete |
| Meropenem | ✓ | ✓ | ✓ | Complete |
| Ciprofloxacin | ✓ | ✓ | ✓ | Complete |
| Metronidazole | ✓ | ✓ | ✓ | Complete |
| Ceftriaxone | ✓ | ✓ | ✓ | Complete |
| Clindamycin | ✓ | ✓ | ✓ | Complete |

### Genus-Level Analysis Status (Script 06)
- **Status:** Complete
- Uses Vancomycin_IV/PO split covariates

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

### Phase 10: Combined Results and Comparative Analyses (2025-12-29)
- [x] Combined species-level and genus-level results from all three methods
- [x] Multi-method concordance analysis (577 robust species associations, 147 genus associations)
- [x] Enterobacteriaceae-focused analysis (499 associations, 2 by all 3 methods)
- [x] Anaerobe-focused analysis (760 associations, 22 by all 3 methods)
- [x] Cefepime vs Pip/Tazo comparative analysis with Venn diagrams
- [x] Pathogenic vs non-pathogenic species comparison
- **Scripts:** `05c_combine_results.R`, `07_cefepime_vs_piptazo_comparison.R`, `08_pathogenic_vs_nonpathogenic_comparison.R`
- **Output:** `results/new_analysis_legacy_data/combined_results/`, `results/new_analysis_legacy_data/cefepime_vs_piptazo/`

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

### Phase 11: Enterobacteriaceae Temporal Dynamics (2025-12-29)
- [x] Species-level analysis distinguishing pathogenic vs commensal Enterobacteriaceae
- [x] Recovery RATE analysis (% change per day) from paired samples
- [x] Cross-sectional analysis by time-since-antibiotics (all samples)
- [x] Individual antibiotic stratification of recovery rates
- [x] Comparison with 07_recovery_analysis.R (see below)
- **Script:** `R/new_analysis_legacy_data/10_enterobact_temporal_dynamics.R`
- **Output:** `results/new_analysis_legacy_data/enterobact_temporal/`

**How this differs from 07_recovery_analysis.R:**

| Aspect | 07_recovery_analysis.R | 10_enterobact_temporal_dynamics.R |
|--------|------------------------|-----------------------------------|
| Taxonomic level | Genus-level | Species-level (Enterobacteriaceae) |
| Enterobacteriaceae | Pooled together | **Split: pathogenic vs commensal** |
| Metric | Total change (Δ) | **Rate of change (%/day)** |
| Focus | Which genera recover after which antibiotics | **Do Enterobact recover faster than anaerobes?** |
| Additional | - | Cross-sectional validation |

---

## Pending Work (TODO)

### High Priority
- [x] ~~Complete vancomycin route split re-analysis~~ - **DONE** (all 9 antibiotics complete)
- [x] ~~Validate that IV vancomycin shows minimal/no direct microbiome effect~~ - **Confirmed** (only 16 robust associations vs 56-137 for other antibiotics)
- [ ] Update network visualizations with route-specific results

### Medium Priority
- [ ] Review multi-method concordance results (species + genus levels)
- [ ] Pathway enrichment analysis (MetaCyc)
- [ ] Integrate findings into summary document

### Future Analyses (see docs/NEW_ANALYSIS_LEGACY_DATA.md)
- [ ] **ARG abundance analysis** - Run same antibiotic-association analysis with antibiotic resistance gene (ARG) abundance results
- [ ] 14-day and 30-day exposure windows
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

### Multi-Method Differential Abundance Results (Species Level)

**Completed antibiotics (7/9) with Vancomycin_IV/PO split covariates:**

| Antibiotic | ALDEx2 sig | MaAsLin3 sig | ANCOM-BC2 sig | Robust (2+ methods) |
|------------|------------|--------------|---------------|---------------------|
| Pip_Tazo | 33 | 120 | 424 | 137 |
| Meropenem | 22 | 84 | 345 | 92 |
| Ciprofloxacin | 34 | 81 | 303 | 89 |
| Metronidazole | 15 | 98 | 87 | 86 |
| TMP_SMX | 11 | 72 | 310 | 70 |
| Cefepime | 11 | 63 | 104 | 56 |
| Vancomycin_IV | 0 | 16 | 307 | 16 |

**Total: 546 robust associations** (found by 2+ methods), **57 found by all 3 methods**

### Highest-Confidence Species Associations (All 3 Methods Agree)

#### Enterococcus Expansion (VRE Risk)
| Antibiotic | Species | Effect (log2FC) |
|------------|---------|-----------------|
| Cefepime | *E. faecalis* | +2.1 |
| Cefepime | *E. raffinosus* | +2.2 |
| Cefepime | *E. durans* | +1.8 |
| Meropenem | *E. faecalis* | +2.1 |
| Meropenem | Multiple *Enterococcus* spp. | +1.6 to +2.1 |

#### Opportunistic Pathogen Expansion
| Antibiotic | Species | Effect (log2FC) |
|------------|---------|-----------------|
| Meropenem | *Pseudomonas aeruginosa* | +2.2 |
| Meropenem | *Staphylococcus aureus* | +2.2 |
| Ciprofloxacin | *Candida tropicalis* | +2.2 |
| Ciprofloxacin | *[Candida] glabrata* | +2.9 |

#### Butyrate-Producer Depletion (Gut Barrier Risk)
| Antibiotic | Species | Effect (log2FC) |
|------------|---------|-----------------|
| Pip_Tazo | *Anaerostipes hadrus* | -1.1 |
| Pip_Tazo | *Roseburia inulinivorans* | -1.0 |
| Pip_Tazo | *Ruminococcus faecis* | -1.2 |
| Pip_Tazo | *[Eubacterium] hallii* | -1.2 |

#### Veillonella Depletion (Consistent Marker)
| Antibiotic | Species | Effect (log2FC) |
|------------|---------|-----------------|
| TMP_SMX | *V. parvula*, *V. tobetsuensis*, etc. (9 spp.) | -1.5 to -2.0 |
| Meropenem | *Veillonella* spp. (4 species) | -1.5 to -1.7 |
| Metronidazole | *V. parvula* | -1.8 |

#### Ciprofloxacin Paradox (Gram-negative Kill → Lactobacillus Bloom)
| Direction | Species | Effect (log2FC) |
|-----------|---------|-----------------|
| ↓ Decreased | *E. coli*, *Salmonella enterica* | -1.6 |
| ↑ Increased | *Lactobacillus* spp. (7+ species) | +1.8 to +3.0 |

#### Vancomycin_IV Results (Limited Direct Effect)
Only 16 robust associations found (vs 56-137 for other antibiotics):
- ↑ *Candida parapsilosis*, *Lactobacillus* spp., *Leuconostoc gelidum*
- ↓ *Bifidobacterium longum*, *Eggerthella* sp.
- Consistent with IV vancomycin not reaching gut lumen directly

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

### Cefepime vs Pip/Tazo Comparison (Phase 10)

Direct comparison of effects across taxonomic levels:

| Level | Cefepime Only | Shared | Pip/Tazo Only |
|-------|---------------|--------|---------------|
| Species | 55 | 56 | 143 |
| Genera | 20 | 18 | 37 |
| Enterococcus | 9 | 4 | 7 |
| Enterobacteriaceae | 16 | 9 | 68 |
| Anaerobes | 18 | 47 | 152 |

**Key findings:**
- Pip/Tazo has ~2.5x more species affected than Cefepime
- Both drugs increase Enterococcus (VRE risk)
- Pip/Tazo massively depletes anaerobes (136 decreased vs 13 increased)
- Cefepime more selective for Enterobacteriaceae
- Effect sizes highly correlated for shared species (r=0.925)

### Pathogenic vs Non-Pathogenic Comparison (Phase 10)

Tested hypothesis: Are pathogenic species more resistant to antibiotic perturbation?

**Result: Hypothesis NOT supported**

| Species Type | Mean Abx Affecting | Mean Robust Abx |
|--------------|--------------------|-----------------|
| E. coli (pathogenic) | 3.0 | 2.0 |
| Other Escherichia (commensal) | 1.7 | 1.4 |
| K. pneumoniae/oxytoca (pathogenic) | 1.5 | 0.5 |
| Other Klebsiella (commensal) | 1.9 | 1.4 |

Pathogens are affected by MORE antibiotics, not fewer. They have more robust associations (consistent detection across methods), possibly due to higher prevalence enabling better statistical power.

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

### Enterobacteriaceae Temporal Dynamics (Phase 11)

**Recovery RATES from paired samples (n=58 pairs):**

| Taxon | Rate (%/day) | p-value | Interpretation |
|-------|--------------|---------|----------------|
| **Total Enterobact** | **+1.45** | **0.036** | Rapid expansion after abx |
| **Pathogenic Enterobact** | **+1.23** | **0.043** | E. coli, K. pneumoniae expand faster |
| Commensal Enterobact | +0.22 | 0.405 | Slower expansion |
| Anaerobes | +0.001 | 0.999 | Essentially no recovery |
| **Enterococcus** | **-1.09** | **0.035** | Returns to baseline |

**Critical finding:** Enterobacteriaceae recover **2,424x faster** than anaerobes after antibiotics stop.

**Recovery rates by individual antibiotic (paired samples):**

| Antibiotic | n | Enterobact Rate | Anaerobe Rate | Ratio |
|------------|---|-----------------|---------------|-------|
| **Cefepime** | 12 | **+1.61%/day** | +0.18%/day | **8.7x** |
| **Pip_Tazo** | 8 | **+1.65%/day** | +1.27%/day | 1.3x |
| **Vancomycin** | 7 | **+1.15%/day** | -0.62%/day | -1.9x |

**Key insight:** After Cefepime exposure, Enterobacteriaceae expand **8.7x faster** than anaerobes recover. This explains why Enterobacteriaceae dominate the post-antibiotic gut - they rapidly fill the ecological niche vacated by depleted anaerobes.

---

## Notes

- All differential abundance analyses complete (9 antibiotics × 3 methods)
- Master script: `R/new_analysis_legacy_data/run_all_analysis.sh`
- 723 high-confidence samples (including SB Stool)
- Full documentation: `docs/NEW_ANALYSIS_LEGACY_DATA.md`
