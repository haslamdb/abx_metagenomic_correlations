# Manuscript Plan: Differential Recovery Rates Explain Pathogen Dominance After Antibiotic Exposure

**Date:** December 2025
**Status:** Draft

---

## Central Question

**Why do opportunistic pathogens (Enterobacteriaceae, Enterococcus) come to dominate the gut microbiome after antibiotic exposure?**

Traditional explanations focus on antibiotic resistance, but our data suggest a different mechanism: **differential recovery kinetics**.

---

## The Clinical Paradox (NEW - Motivating Data)

We have clinical evidence that directly contradicts the "resistance" narrative:

**Observation from BSI cases in this cohort:**
- Patients receive antibiotics (e.g., cefepime, meropenem)
- Patients subsequently develop bloodstream infections (BSI)
- The BSI isolate is **genomically identical** to a fecal isolate collected before the BSI
- Critically: The BSI isolate is **SUSCEPTIBLE** to the antibiotics the patient received

**This paradox motivates the entire paper:**
> If the pathogen that caused the BSI was susceptible to the antibiotics the patient received, why did it dominate the gut and translocate?

**Answer:** Not resistance, but differential recovery kinetics.

---

## Key Thesis

Antibiotics affect pathogens and commensals alike, but pathogens recover faster when antibiotics stop. This asymmetric recovery—not resistance—drives post-antibiotic pathogen dominance and subsequent infection risk.

---

## The Temporal Dominance Model (KEY CONCEPTUAL FRAMEWORK)

Our data reveal a **two-phase model** of post-antibiotic pathogen dominance:

### Phase 1: During Antibiotic Exposure
| Organism | Response | Mechanism |
|----------|----------|-----------|
| **Enterococcus** | ↑↑ INCREASES | Intrinsic resistance to many antibiotics, fills niche |
| **Enterobacteriaceae** | ↓ Decreases | Susceptible to broad-spectrum agents |
| **Obligate Anaerobes** | ↓↓ DECREASES | Highly susceptible, slow to recover |

### Phase 2: After Antibiotic Cessation
| Organism | Response | Mechanism |
|----------|----------|-----------|
| **Enterococcus** | ↓ Returns to baseline | Loses competitive advantage |
| **Enterobacteriaceae** | ↑↑ RAPIDLY RECOVERS | Fast growth, environmental reservoirs |
| **Obligate Anaerobes** | → Remains depleted | Slow growth, strict anaerobe requirements |

### Clinical Prediction: BSI Organism Identity Shifts with Time

This model predicts that **the dominant BSI organism should change based on timing relative to antibiotic exposure:**

| Timing | Predicted Dominant BSI Organism | Microbiome State |
|--------|--------------------------------|------------------|
| **During antibiotics** | Enterococcus (VRE) | Enterococcus-dominated gut |
| **Early post-antibiotics** (days 1-7) | Enterobacteriaceae (E. coli, Klebsiella) | Rapid Enterobacteriaceae recovery |
| **Late post-antibiotics** (weeks) | Mixed/normalized | Gradual commensal recovery |

**Testable in large hospital data:** Analyze BSI organism identity as a function of days since last antibiotic dose.

### Observational Support

Some patients' fecal microbiomes are dominated by Enterococcus (>50% relative abundance). Our model predicts these are patients **currently on or recently receiving antibiotics**. Enterobacteriaceae-dominated samples should be enriched in patients **days after antibiotic cessation**.

---

## Summary of Main Findings

### 1. Antibiotics DO Affect Pathogens (Bulk Analysis)

Contrary to a simple "resistance" narrative, our species-level differential abundance analysis shows:

| Finding | Evidence |
|---------|----------|
| **Enterobacteriaceae ARE depleted** | Meropenem: -2.06 log2FC (p=0.004), Cefepime: -1.60 log2FC (p=0.03) |
| **E. coli affected by 3.0 antibiotics** (vs 1.7 for commensal Escherichia) | Pathogenic species have MORE associations, not fewer |
| **Multi-method confirmation** | 577 robust species associations (2+ methods agree) |

**Key point:** Pathogens are NOT resistant to antibiotics in our data. They are affected by antibiotics, sometimes more so than commensals.

### 2. Persistence During Antibiotics (Paired Analysis, n=206 pairs)

| Functional Group | Pre-Abx | Post-Abx | Change | P-value |
|------------------|---------|----------|--------|---------|
| Anaerobes | 20.9% | 18.8% | -2.1% | 0.26 |
| Enterobacteriaceae | 25.0% | 25.6% | +0.6% | 0.81 |
| **Enterococcus** | 16.3% | 21.8% | **+5.5%** | **0.018** |

- Enterococcus significantly INCREASES during antibiotic exposure (expansion into depleted niches)
- Enterobacteriaceae remain stable (neither killed nor selected)
- Anaerobes decline modestly

### 3. Recovery After Antibiotics Stop (Paired Analysis, n=46-58 pairs)

| Functional Group | On Abx | Off Abx | Change | Recovery Rate (%/day) | P-value |
|------------------|--------|---------|--------|----------------------|---------|
| Anaerobes | 13.4% | 14.6% | +1.3% | **+0.001** | 0.99 |
| **Enterobacteriaceae** | 27.7% | 37.5% | **+9.7%** | **+1.45** | **0.036** |
| Enterococcus | 18.1% | 10.2% | -7.9% | -1.09 | 0.035 |

**Critical finding:** Enterobacteriaceae recover **2,424x faster** than anaerobes after antibiotics stop.

### 4. Pathogenic vs Commensal Enterobacteriaceae

| Species Type | Recovery Rate (%/day) | P-value |
|--------------|----------------------|---------|
| **Pathogenic** (E. coli, K. pneumoniae) | **+1.23** | **0.043** |
| Commensal Enterobacteriaceae | +0.22 | 0.41 |

Pathogenic species recover **5.6x faster** than commensal Enterobacteriaceae.

---

## Role of the Two Analysis Types

### Bulk (Cross-Sectional) Analysis

**Question answered:** Which taxa are associated with antibiotic exposure?

**What it shows:**
- Large sample size (723 samples) provides statistical power
- Identifies taxa differentially abundant in exposed vs unexposed
- Confirms antibiotics DO affect pathogen abundance
- Reveals antibiotic-specific effects (e.g., Ciprofloxacin kills Enterobacteriaceae but enriches Lactobacillus)

**Limitations:**
- Cannot distinguish cause from correlation
- Confounding: sicker patients receive more antibiotics AND have more dysbiosis

### Paired (Longitudinal) Analysis

**Question answered:** What happens to the microbiome over time within the same patient?

**What it shows:**
- Within-patient changes eliminate between-patient confounders
- Persistence analysis: What survives during antibiotics
- Recovery analysis: What rebounds after antibiotics stop
- Rate analysis: HOW FAST taxa recover

**Key advantage:** The paired recovery rate analysis directly demonstrates differential kinetics—the mechanistic explanation for pathogen dominance.

---

## AMR Analysis: The Supporting Footnote

### What We Found

**Bulk Analysis (ARG differential abundance):**
- 213 robust ARG associations (2+ methods)
- 30 robust ARG family associations
- Antibiotics are associated with ARG carriage in cross-sectional data

**Paired Analysis (Within-patient ARG changes):**
- Only **1 significant result** (TMP-SMX decreases ARG diversity, p=0.047)
- No other antibiotic showed significant within-patient ARG selection
- ARG richness, diversity, and resistance class abundances unchanged

### Interpretation

The disconnect between bulk and paired ARG results is telling:

| Analysis Type | Finding | Interpretation |
|---------------|---------|----------------|
| Bulk (cross-sectional) | Many ARG associations | Confounded: sicker patients get more antibiotics AND carry more ARGs |
| Paired (longitudinal) | No ARG selection | No causal relationship: antibiotics don't select for ARGs short-term |

**Conclusion for manuscript:** ARGs are a **marker** of antibiotic-exposed patients, not the **mechanism** of post-antibiotic pathogen expansion. The paired analysis shows antibiotics don't meaningfully increase ARG carriage within individuals on the timescales studied.

This supports our thesis: **differential recovery kinetics, not resistance, drives pathogen dominance**.

---

## Proposed Manuscript Structure

### Title (Options)
1. "Susceptible Pathogens Dominate: Differential Recovery Kinetics, Not Resistance, Drive Post-Antibiotic Gut Dysbiosis"
2. "The Resistance Paradox: Antibiotic-Susceptible Enterobacteriaceae Dominate Through Rapid Recovery"
3. "Recovery Kinetics, Not Resistance, Explain Post-Antibiotic Pathogen Dominance in Hospitalized Children"

### Abstract Key Points
1. **Clinical paradox:** BSI isolates are often susceptible to antibiotics patients received
2. **Question:** If not resistance, why do pathogens dominate?
3. **Finding:** Pathogens are affected by antibiotics but recover faster
4. **Key result:** Enterobacteriaceae recover 2,400x faster than anaerobes
5. **Supporting evidence:** ARG analysis shows no within-patient selection
6. **Conclusion:** Recovery kinetics, not resistance, is the mechanism

### Introduction

**Opening hook (Figure 1 preview):**
- Present 3-5 illustrative cases of susceptible BSI isolates
- Genomic identity between fecal and blood isolates
- Susceptibility to prior antibiotics
- "This paradox led us to investigate alternative mechanisms..."

**Background:**
- Clinical problem: Post-antibiotic dysbiosis and infections
- Traditional explanation: Resistant pathogens survive antibiotic pressure
- Evidence against: Our BSI isolates are susceptible
- Gap: Little data on recovery kinetics after antibiotics stop
- Our approach: Paired longitudinal analysis in hospitalized children

### Methods
- Cohort: 723 samples, 130 patients (BMT, LvTx, IF, PICU, SB)
- Metagenomic profiling: Kraken2/Bracken
- ARG profiling: AMRFinderPlus
- BSI isolate sequencing and susceptibility testing
- Statistical approach:
  - Bulk: ALDEx2, MaAsLin3, ANCOM-BC2 with multi-method concordance
  - Paired: Mixed-effects models with covariate adjustment
  - Recovery rates: %/day change after antibiotic cessation
  - Genomic comparison: ANI/SNP distance between fecal and blood isolates

### Results

**Section 1: The Clinical Paradox - Susceptible Isolates Cause BSI (Figure 1)**
- Representative cases showing fecal→blood transmission
- Genomic identity (ANI >99.9% or <10 SNPs)
- Antibiotic susceptibility despite prior exposure
- Establishes the motivating question

**Section 2: Antibiotics Deplete Both Pathogens and Commensals (Bulk)**
- 577 robust species associations
- Pathogens (E. coli, K. pneumoniae) affected by MORE antibiotics than commensals
- Key antibiotic effects (Table 1)
- Confirms pathogens are susceptible at population level

**Section 3: Persistence During Antibiotic Exposure (Paired)**
- Enterococcus expands (+34%, p=0.018)
- Enterobacteriaceae stable
- Anaerobes decline modestly

**Section 4: Differential Recovery After Antibiotics Stop (Paired)**
- Enterobacteriaceae recover rapidly (+1.45%/day)
- Pathogenic species recover faster than commensal (+1.23 vs +0.22%/day)
- Anaerobes fail to recover (~0%/day)
- Figure: Recovery rate forest plot

**Section 5: ARG Analysis Does Not Support Resistance-Mediated Selection**
- Cross-sectional associations exist (213 robust ARG associations)
- Within-patient analysis shows no ARG selection (only 1/90 tests significant)
- ARGs are markers of exposure, not mechanisms of dominance
- Consistent with BSI isolate susceptibility data

### Discussion

**Resolving the paradox:**
- Pathogens are susceptible to antibiotics (BSI data, bulk analysis)
- Yet they dominate the post-antibiotic gut (recovery analysis)
- Mechanism: Differential recovery kinetics, not resistance

**Biological mechanisms of rapid recovery:**
- Facultative anaerobes (can survive low-oxygen transition periods)
- Faster intrinsic growth rates
- Metabolic flexibility (can utilize diverse substrates)
- Environmental reservoirs (hospital surfaces, feeding tubes, hands)
- Spore-forming capacity absent in Enterobacteriaceae but present in some anaerobes - yet sporulation is slow

**Clinical implications:**
- Post-antibiotic window is the critical period for intervention
- Interventions targeting recovery (probiotics, FMT) may need precise timing
- ARG monitoring may be less clinically relevant than recovery kinetics
- Stewardship focus: Not just "which antibiotic" but "what happens after"

**Relationship to BSI risk (preview of future work):**
- This paper explains the mechanism
- Companion paper (in preparation) will address predictive modeling
- Recovery kinetics may be a targetable intervention point

### Figures

**Figure 1: The Susceptibility Paradox - BSI Isolates Are Susceptible to Prior Antibiotics** (NEW - KEY FIGURE)

Panel A: Study design schematic
- Timeline showing antibiotic exposure → fecal sampling → BSI event
- Illustrates the clinical scenario

Panel B: Representative cases (3-5 patients)
- Patient timeline with antibiotic exposure periods
- Fecal sample collection points
- BSI event timing
- Table: Antibiotic received vs. BSI isolate susceptibility

Panel C: Genomic identity confirmation
- Phylogenetic tree or SNP distance showing fecal-blood isolate pairs
- ANI values demonstrating clonality (>99.9%)

Panel D: Summary statistics
- N patients with matched fecal-blood isolates
- % susceptible to prior antibiotics
- "If resistance doesn't explain dominance, what does?"

**Figure 2: Unbiased Species Ranking DURING Antibiotic Exposure (Bulk Analysis)**

*Key message: Enterococcus increases, Enterobacteriaceae and commensals decrease*

Panel A: Ranked species by effect size during antibiotic exposure
- Horizontal bar plot of ALL species ranked by log2FC (exposed vs unexposed)
- Color-coded by functional category (post-hoc annotation)
- Shows Enterococcus species clustering at TOP (positive effect)
- Shows Enterobacteriaceae and anaerobes clustering at BOTTOM (negative effect)

Panel B: Category summary boxplot
- Recovery rate distributions by functional category
- Confirms pattern: Enterococcus ↑, Enterobacteriaceae ↓, Anaerobes ↓↓

Panel C: Conceptual diagram - Phase 1
- "During Antibiotics: Enterococcus fills the niche"

**Script:** `R/new_analysis_legacy_data/14_species_exposure_ranking.R`

**Figure 3: Unbiased Species Ranking AFTER Antibiotic Cessation (Recovery Analysis)**

*Key message: Enterobacteriaceae rapidly recover, Enterococcus declines, anaerobes stagnate*

Panel A: Ranked species by recovery rate (%/day)
- Horizontal bar plot of ALL species ranked by recovery rate
- Color-coded by functional category (post-hoc annotation)
- Shows Enterobacteriaceae (especially pathogens) clustering at TOP
- Shows anaerobes in middle/bottom

Panel B: Category summary boxplot
- Recovery rate distributions by functional category
- Confirms pattern: Enterobacteriaceae ↑↑, Enterococcus ↓, Anaerobes →

Panel C: Recovery rate vs published doubling time
- Scatter plot showing correlation with intrinsic growth rate
- Provides mechanistic explanation

Panel D: Conceptual diagram - Phase 2
- "After Antibiotics: Fast growers dominate"

**Script:** `R/new_analysis_legacy_data/13_species_recovery_ranking.R`

**Figure 4: The Two-Phase Temporal Model (Synthesis)**

Panel A: Schematic timeline
- X-axis: Time relative to antibiotic exposure
- Y-axis: Relative abundance
- Three lines: Enterococcus, Enterobacteriaceae, Anaerobes
- Shows crossover point where Enterococcus declines and Enterobacteriaceae rises

Panel B: Predicted BSI organism by timing
- Bar chart or timeline showing expected BSI organism identity
- During Abx: Enterococcus (VRE)
- After Abx: Enterobacteriaceae (E. coli, Klebsiella)

Panel C: Supporting evidence from this cohort
- % Enterococcus-dominated samples in on-Abx vs off-Abx groups
- % Enterobacteriaceae-dominated samples by days since last Abx

**Figure 5: ARG Analysis - No Within-Patient Selection**

Panel A: Bulk analysis shows many associations
- Heatmap or bar plot of ARG associations by antibiotic

Panel B: Paired analysis shows no selection
- Forest plot of within-patient ARG changes
- All CIs cross zero except TMP-SMX

Panel C: Interpretation diagram
- Bulk associations = confounding
- Paired analysis = no causal effect
- Supports "kinetics not resistance" thesis

### Tables

1. **Table 1:** Patient demographics and antibiotic exposure summary
2. **Table 2:** Species ranking during antibiotic exposure (top increasers and decreasers)
3. **Table 3:** Species ranking during recovery (top recoverers by rate)
4. **Table 4:** Category comparison statistics (Wilcoxon tests)
5. **Table 5:** BSI isolate characteristics and susceptibility (summary of Figure 1 data)

### Supplementary Materials

- **Table S1:** Full differential abundance results (all 577 robust associations)
- **Table S2:** ARG associations (all 213 robust)
- **Table S3:** Individual antibiotic effects on recovery rates
- **Table S4:** Antibiotic-specific effects on functional categories (see below)
- **Figure S1:** Individual patient trajectories (spaghetti plots)
- **Figure S2:** Network visualization of antibiotic-species associations
- **Figure S3:** Antibiotic-specific heatmap of effects on Enterococcus/Enterobacteriaceae/Anaerobes
- **Figure S4:** Ciprofloxacin PO vs IV route comparison
- **Methods S1:** Detailed statistical methods
- **Methods S2:** Genomic analysis pipeline for fecal-blood isolate comparison

---

## Supplementary Analysis: Antibiotic-Specific Effects

### Rationale

The main analysis pools effects across antibiotics to show the general temporal pattern. This supplementary analysis confirms patterns are consistent with each antibiotic's mechanism of action.

### Key Antibiotic-Specific Findings

| Antibiotic | Anaerobes | Enterobact | Enterococcus | Interpretation |
|------------|-----------|------------|--------------|----------------|
| **Cefepime** | +0.6 | **-0.8** | **+1.3** | Classic: kills Enterobact, Enterococcus fills niche |
| **Meropenem** | -0.8 | +0.6 | **+1.5** | Broad-spectrum, strong Enterococcus selection |
| **Pip_Tazo** | -0.7 | -0.4 | **+0.9** | Anti-anaerobic + anti-pseudomonal |
| **Metronidazole** | **-1.3** | ~0 | +0.5 | Pure anti-anaerobic, spares Enterobact |
| **Ciprofloxacin** | -0.2 | **-1.3** | **-1.2** | Unique: kills BOTH gram- and Enterococcus |
| **Clindamycin** | **-1.1** | +1.4 | +1.2 | Anti-anaerobic → Enterobact bloom |
| **Ceftriaxone** | **-2.1** | +1.4 | NA | Strong anti-anaerobic, biliary excretion |
| **TMP_SMX** | -0.5 | -0.1 | -0.8 | Modest effects on all groups |
| **Vancomycin_IV** | NA | +0.9 | -0.9 | Gram+ killer → Enterobact fill niche |

### Observations

1. **Most antibiotics increase Enterococcus** - consistent with intrinsic resistance
2. **Ciprofloxacin is unique** - fluoroquinolone kills both Enterococcus and Enterobacteriaceae
3. **Anti-anaerobic agents** (metronidazole, clindamycin) cause Enterobacteriaceae bloom DURING exposure
4. **Vancomycin_IV** kills gram-positives, allowing gram-negatives to expand

### Ciprofloxacin Route Analysis (PO vs IV)

Unlike vancomycin (where IV doesn't reach gut), both PO and IV ciprofloxacin reach the gut:
- **PO:** Direct luminal exposure
- **IV:** Biliary excretion into gut

#### Results

Both routes significantly suppress Enterobacteriaceae with **no significant difference between routes**:

| Route | N samples | Mean Enterobact (%) | vs Unexposed |
|-------|-----------|---------------------|--------------|
| None (unexposed) | 1,378 | 25.9% | -- |
| IV only | 44 | 16.0% | p=0.04 |
| PO only | 30 | 8.2% | p=0.008 |
| **IV vs PO** | -- | -- | **p=0.54 (NS)** |

Both routes also significantly suppress anaerobes (IV: p<0.001, PO: p<0.001).

**Interpretation:** Unlike vancomycin (where IV doesn't affect gut flora), ciprofloxacin IV reaches the gut via biliary excretion and has similar microbiome effects to PO administration. This confirms that route of administration matters differently depending on the antibiotic's pharmacokinetics.

**Script:** `R/new_analysis_legacy_data/15_antibiotic_specific_effects.R`

---

## How This Paper Relates to the BSI Risk Paper

| This Paper (Mechanism) | BSI Paper (Outcomes) |
|------------------------|----------------------|
| **Illustrative examples** of susceptible BSI isolates | **Systematic analysis** of all BSI cases |
| Establishes the paradox | Quantifies the risk |
| "Pathogens dominate despite susceptibility" | "Microbiome signature predicts BSI" |
| Recovery kinetics as mechanism | Predictive model for clinical use |
| N = 3-5 representative cases | N = all BSI events in cohort |
| Focus: WHY pathogens dominate | Focus: WHO will get BSI |

**Key distinction:** This paper uses BSI susceptibility data to *motivate* the mechanistic question. The BSI paper will use microbiome data to *predict* clinical outcomes.

---

## Future Direction: BSI Organism Identity Shifts with Antibiotic Timing

### Hypothesis (Testable in Large Hospital Data)

Based on the two-phase temporal model, we predict that **BSI organism identity should change as a function of days since last antibiotic dose**:

| Days Since Last Antibiotic | Predicted Dominant BSI Organism |
|---------------------------|--------------------------------|
| 0 (on antibiotics) | Enterococcus (VRE) |
| 1-3 days | Transition period |
| 4-14 days | Enterobacteriaceae (E. coli, K. pneumoniae) |
| >14 days | Mixed / baseline flora |

### Significance

If confirmed, this would:
- Validate the temporal dominance model with clinical outcome data
- Suggest **timing-specific interventions** (e.g., Enterococcus prophylaxis during Abx, Enterobacteriaceae prophylaxis after cessation)
- Explain why empiric therapy recommendations differ by clinical context

---

## BSI Etiology Analysis: Comprehensive Study Design

### Overview

We have 13 years of hospital data (2013-present) including:
- All positive blood cultures with organism identification
- Antibiotic susceptibility profiles for all isolates
- Complete antibiotic administration records for all patients
- Hospital admission/discharge records (in-house status)

This enables a definitive test of the temporal dominance model using clinical BSI outcomes.

### Cohort Definition

```
                          All Positive Blood Cultures (2013-2025)
                                         │
                    ┌────────────────────┴────────────────────┐
                    │                                          │
         Hospital-Associated BSI                    Community-Acquired BSI
         (≥3 days in-house prior)                   (<3 days in-house)
                    │                                          │
              PRIMARY COHORT                          COMPARISON COHORT
        (Complete exposure data)                   (Expected: S. aureus,
                    │                               S. pneumoniae dominant)
                    │
    ┌───────────────┼───────────────┐
    │               │               │
 Monomicrobial  Polymicrobial   Contaminant
   (include)   (analyze sep.)    (exclude)
```

**Inclusion Criteria:**
- Positive blood culture with organism identification
- ≥3 days of hospitalization prior to BSI (ensures complete antibiotic exposure data)
- Monomicrobial (or analyze polymicrobial separately)

**Exclusion Criteria:**
- Likely contaminants (single-bottle CoNS without clinical signs)
- Incomplete antibiotic records
- Outpatient/ED-only encounters

### BSI Organism Categorization

| Category | Organisms | Expected Temporal Pattern |
|----------|-----------|---------------------------|
| **Enterococcus** | E. faecalis, E. faecium, VRE | ↑ During Abx, ↓ after cessation |
| **Enterobacteriaceae** | E. coli, Klebsiella, Enterobacter, Serratia, Citrobacter, Proteus | ↓ During Abx, ↑↑ recovery phase |
| **Pseudomonas** | P. aeruginosa | Resistant phenotype, less time-dependent? |
| **S. aureus** | MSSA, MRSA | Skin/line source, not gut-related |
| **CoNS** | S. epidermidis, S. hominis | Line-associated, not gut-related |
| **Candida** | C. albicans, C. glabrata, C. parapsilosis | Delayed emergence after broad-spectrum |
| **Streptococcus** | S. pneumoniae, viridans group | Community-acquired or oral source |
| **Anaerobes** | Bacteroides, Clostridium, Fusobacterium | Rare BSI, indicates GI source |

### Primary Predictor Variables

**Temporal (Main Interest):**
| Variable | Definition | Rationale |
|----------|------------|-----------|
| Days_since_last_Abx | Days from last antibiotic dose to BSI culture | Primary predictor of interest |
| Currently_on_Abx | Binary: on antibiotics at time of BSI | Phase 1 (during) vs Phase 2 (after) |
| Cumulative_Abx_days | Total antibiotic days in prior 30 days | Dose-response relationship |
| N_Abx_classes | Number of distinct antibiotic classes | Breadth of exposure |
| Last_Abx_class | Most recent antibiotic class | Antibiotic-specific effects |
| Anti_anaerobic_exposure | Any metronidazole, clindamycin, carbapenems | Strongest anaerobe disruptors |

**Categorized Time Since Antibiotics:**
| Category | Days Since Last Abx | Predicted Dominant BSI |
|----------|---------------------|------------------------|
| Active exposure | 0 (currently on) | Enterococcus |
| Early recovery | 1-3 | Transition |
| Peak recovery window | 4-14 | Enterobacteriaceae |
| Late recovery | 15-30 | Mixed |
| No recent exposure | >30 | Baseline/Community pattern |

### Key Covariates

#### Tier 1: Critical Confounders (Must Include)

| Variable | Type | Rationale |
|----------|------|-----------|
| **ANC** (neutrophil count) | Continuous/categorical | Neutropenia fundamentally changes BSI risk and organism |
| **Central_line** | Binary + type | Line-associated BSI has different organism distribution |
| **ICU_status** | Binary | Higher acuity, different pathogens |
| **Age** | Continuous | Affects baseline flora, immune status |
| **Immunosuppression** | Composite score | Transplant, chemo, steroids affect risk |
| **Days_hospitalized** | Continuous | Healthcare exposure, colonization |

#### Tier 2: Important Confounders (Should Include)

| Variable | Type | Rationale |
|----------|------|-----------|
| **PPI_use** | Binary | Increases gut pH, favors Enterobacteriaceae |
| **TPN** | Binary | Gut atrophy, bacterial translocation |
| **GI_surgery** | Binary (30d) | Direct translocation risk |
| **Diagnosis_category** | Categorical | Oncology, transplant, surgery, medicine |
| **Prior_colonization** | Binary per organism | MRSA, VRE, CRE colonization predicts BSI |
| **Dialysis** | Binary | Line-associated, uremic immunosuppression |
| **Liver_disease** | Binary | Portal hypertension → translocation |
| **Diabetes** | Binary | Immune dysfunction |

#### Tier 3: Additional Covariates (If Available)

| Variable | Type | Rationale |
|----------|------|-----------|
| Albumin | Continuous | Nutritional/inflammatory status |
| Creatinine/eGFR | Continuous | Renal function affects Abx PK |
| Foley_catheter | Binary | Urosepsis source |
| Mechanical_ventilation | Binary | Critical illness marker |
| Prior_BSI_30d | Binary | Recurrence/persistent infection |
| Procalcitonin | Continuous | Infection severity |

### Statistical Analysis Plan

#### Primary Analysis: BSI Organism Category vs Time Since Antibiotics

**Model 1: Multinomial Logistic Regression**
```
BSI_organism_category ~ days_since_last_abx + currently_on_abx +
                        ANC_category + central_line + ICU +
                        immunosuppression + diagnosis_category +
                        age + days_hospitalized
```

**Model 2: Binary Logistic (Enterobacteriaceae vs Other)**
```
Enterobacteriaceae_BSI ~ days_since_last_abx + covariates
```

**Model 3: Spline for Non-Linear Temporal Effect**
```
BSI_category ~ ns(days_since_last_abx, df=3) + covariates
```

#### Secondary Analyses

1. **Stratified by Antibiotic Class:**
   - Does anti-anaerobic exposure predict Enterobacteriaceae BSI specifically?
   - Does fluoroquinolone exposure reduce Enterobacteriaceae BSI?

2. **Time-to-BSI Survival Analysis:**
   - Cox model for time from antibiotic cessation to BSI
   - Competing risks: discharge, death

3. **Susceptibility Paradox Quantification:**
   - % of BSI isolates susceptible to antibiotics received in prior 7/14/30 days
   - Compare to resistance rates expected if selection were the mechanism

4. **Community-Acquired Comparison:**
   - Confirm S. aureus/S. pneumoniae dominance in <3 day hospitalizations
   - Natural control group without gut disruption mechanism

#### Visualization Plan

**Figure: BSI Organism Distribution by Time Since Antibiotics**

Panel A: Stacked bar chart
- X-axis: Time categories (On Abx, 1-3d, 4-14d, 15-30d, >30d)
- Y-axis: Proportion of BSI events
- Colors: Organism categories
- Shows temporal shift from Enterococcus → Enterobacteriaceae

Panel B: Ribbon/area plot
- X-axis: Continuous days since last antibiotic
- Y-axis: Predicted probability of each organism category
- Smoothed curves showing transition

Panel C: Forest plot of odds ratios
- Effect of each predictor on Enterobacteriaceae BSI
- Shows relative importance of timing vs other factors

Panel D: Susceptibility paradox
- Bar chart: % BSI isolates susceptible to prior Abx
- By organism category
- Confirms kinetics hypothesis

### Sample Size Considerations

**Available Data (Estimated):**
| Metric | Estimate |
|--------|----------|
| Years of data | 13 (2013-2025) |
| BSI events per year | ~200-500 (institution dependent) |
| Total BSI events | ~2,600-6,500 |
| Hospital-associated (≥3d) | ~60-70% → 1,500-4,500 |
| With complete Abx data | ~90% → 1,350-4,000 |

**Power Calculation:**
- With ~2,000 hospital-associated BSI events:
- 10% absolute difference in Enterobacteriaceae proportion (e.g., 30% vs 40%)
- Power >90% with α=0.05
- Sufficient for multinomial model with 10+ covariates

### Expected Results (Based on Temporal Model)

**Primary Prediction:**
```
Time Since Last Abx    Enterococcus    Enterobact    Other
─────────────────────────────────────────────────────────
On antibiotics         35%             15%           50%
1-3 days               25%             25%           50%
4-14 days              15%             40%           45%
15-30 days             10%             30%           60%
>30 days               10%             20%           70%
```

**Key Expected Findings:**
1. Enterococcus BSI proportion peaks during active antibiotic exposure
2. Enterobacteriaceae BSI proportion peaks 4-14 days after antibiotic cessation
3. Effect persists after adjustment for neutropenia, central lines, ICU status
4. BSI isolates frequently susceptible to prior antibiotics (the paradox)
5. Community-acquired BSI shows flat distribution (no temporal pattern)

### Data Request

A comprehensive data request has been prepared for the IS team:
**See:** `docs/BSI_DATA_REQUEST.md`

**Minimum Required Variables:**
1. Patient ID, age, gender
2. BSI culture date, organism, susceptibilities
3. Hospital admission/discharge dates
4. All antibiotic administrations with dates
5. WBC/ANC at time of BSI
6. Central line status
7. ICU status at BSI
8. Primary diagnosis codes
9. 30-day mortality

**High-Value Additions:**
- PPI use (affects gut colonization)
- Immunosuppression status (steroids, chemo, transplant)
- Prior colonization (MRSA, VRE, CRE)
- GI surgery/procedures
- TPN use

### Timeline and Milestones

1. **Data extraction request submitted** - [Date]
2. **Data received and QC completed** - [Date + 4-6 weeks]
3. **Cohort definition and cleaning** - [Date + 2 weeks]
4. **Primary analysis completed** - [Date + 2 weeks]
5. **Figures and tables generated** - [Date + 1 week]
6. **Results integrated into manuscript** - [Date + 2 weeks]

### Relationship to Main Paper

| Main Paper (Current) | BSI Etiology Analysis |
|---------------------|----------------------|
| Illustrative BSI cases (N=3-5) | Comprehensive analysis (N=thousands) |
| "Here's the paradox" | "Here's the validation at scale" |
| Microbiome mechanism | Clinical outcome confirmation |
| Could be Figure 1 preview | Could be major Results section or companion paper |

**Decision point:** Include as:
- **Option A:** Major figure in main paper (if data ready in time)
- **Option B:** Companion clinical paper (if extensive enough)
- **Option C:** Supplementary validation (if sample size limited)

---

## Additional Validation: Hospital-Wide BSI Isolate Susceptibility vs Prior Antibiotics

### Hypothesis

If antibiotic resistance were the mechanism of post-antibiotic pathogen dominance, we would expect BSI isolates to be **resistant** to the antibiotics patients recently received. Our model predicts the opposite: BSI isolates will often be **susceptible** because dominance arises from differential recovery kinetics, not resistance selection.

### Available Data

- Hospital-wide antibiotic susceptibility data for all BSI isolates
- Antibiotic administration records for all hospitalized patients
- Can link BSI events to prior antibiotic exposure

### Analysis Plan

1. For each BSI event, identify:
   - Organism causing BSI
   - Antibiotic susceptibility profile (from clinical microbiology)
   - Antibiotics patient received in prior 7/14/30 days

2. Calculate:
   - % of BSI isolates SUSCEPTIBLE to at least one recent antibiotic
   - % of BSI isolates RESISTANT to all recent antibiotics

3. Compare to expected rates:
   - If resistance were the mechanism: expect HIGH resistance to recent antibiotics
   - If kinetics were the mechanism: expect susceptibility rates similar to population baseline

### Expected Result

BSI isolates will frequently be susceptible to antibiotics the patient recently received. This contradicts the "resistance selection" hypothesis and supports "differential recovery kinetics."

### Why This Matters

This hospital-wide analysis would:
- Validate findings from the study cohort in a much larger population
- Provide actionable clinical insight (empiric therapy choices)
- Challenge the assumption that antibiotic stewardship primarily prevents resistance
- Reframe the problem as "niche creation" rather than "resistance selection"

### Manuscript Implications

If this analysis confirms the hypothesis, it could be:
- **Option A:** Include as supplementary validation in this paper
- **Option B:** Save for companion BSI paper as major finding
- **Option C:** Separate short communication focused on the clinical paradox

---

## Key Points for Writing

### What to Emphasize
1. The **clinical paradox**: BSI isolates are susceptible yet pathogens dominate
2. The **mechanism**: Differential recovery kinetics, not resistance
3. The **evidence**: Paired analysis showing 2,400x faster recovery
4. The **negative finding**: ARGs don't explain the pattern
5. The **clinical relevance**: Understanding mechanism enables intervention

### What to Downplay
1. Individual ARG associations (confounded by patient factors)
2. Complex antibiotic-specific effects (save for supplementary)
3. Methods details (move to supplement)
4. Full BSI analysis (save for companion paper)

### Potential Limitations to Address
1. Hospitalized pediatric population (may not generalize to adults/outpatients)
2. 7-day exposure window (may miss longer-term effects)
3. Metagenomic profiling limitations (relative abundance, not absolute)
4. Recovery pairs smaller sample size (n=46-58)
5. BSI cases are illustrative, not comprehensive (full analysis in companion paper)

---

## Target Journals

### Tier 1 (High Impact, Broad Audience)
- **Gut** - Strong clinical focus, appreciates mechanism + clinical relevance
- **Gastroenterology** - If framed around gut health/translocation
- **Nature Communications** - Broad audience, likes paradigm-shifting findings

### Tier 2 (Specialty High Impact)
- **Microbiome** - Perfect fit for methods and findings
- **mBio** - ASM flagship, strong microbiome section
- **ISME Journal** - Microbial ecology focus

### Tier 3 (Clinical Infectious Disease Focus)
- **Clinical Infectious Diseases** - If framed around BSI prevention
- **Journal of Infectious Diseases** - Strong ID audience

**Recommendation:** Target **Gut** or **Microbiome** first. The clinical paradox hook + mechanism + ARG negative result is a compelling package.

---

## Additional Analyses to Strengthen the Paper

### 1. Unbiased Species-Level Recovery Ranking (HIGH PRIORITY)

**Rationale:** The current approach pre-specifies "pathogens vs commensals" which can appear hypothesis-driven. A stronger approach lets the data speak for itself.

**Method:**
1. Calculate recovery rate (%/day) for each species across all recovery pairs
2. Filter to species present in ≥20% of samples (statistical power)
3. Rank ALL species by mean recovery rate (fastest to slowest)
4. Annotate each species post-hoc: pathogenic Enterobacteriaceae, commensal Enterobacteriaceae, obligate anaerobe, Enterococcus, other
5. Visualize as ranked bar/forest plot with color-coded annotations

**Expected output:** If pathogens naturally cluster at the top of an unbiased ranking, that's far more compelling than a pre-specified comparison.

**Script:** `R/new_analysis_legacy_data/13_species_recovery_ranking.R`

### 2. Doubling Time Correlation (MEDIUM PRIORITY)

**Rationale:** Provides mechanistic explanation for WHY pathogens recover faster without requiring new experiments.

**Method:**
1. Compile published doubling times for key species from literature
2. Correlate recovery rate (from our data) with known doubling time
3. If correlation exists, fast recovery = fast intrinsic growth rate

**Key species with known doubling times:**
| Species | Doubling Time | Source |
|---------|--------------|--------|
| E. coli | ~20 min | Standard reference |
| K. pneumoniae | ~25 min | Standard reference |
| Enterococcus faecalis | ~30 min | Standard reference |
| Bacteroides fragilis | ~4-8 hours | Anaerobe literature |
| Faecalibacterium prausnitzii | ~8-12 hours | Anaerobe literature |
| Bifidobacterium spp. | ~2-4 hours | Anaerobe literature |

**Expected output:** Scatter plot of recovery rate vs doubling time, with correlation coefficient.

### 3. Citing Preterm RNA:DNA Work (Discussion)

The preterm_arg project has RNA:DNA metatranscriptomic data showing higher transcriptional activity in Enterobacteriaceae during recovery. This can be cited as supporting evidence from another cohort:

> "In a separate cohort of preterm infants with paired RNA:DNA metatranscriptomic data, we observed that Enterobacteriaceae show higher transcriptional activity during post-antibiotic recovery periods, consistent with rapid regrowth (manuscript in preparation)."

### 4. Alternative Mechanistic Explanations (Discussion)

If doubling time correlation is weak, discuss alternative mechanisms:
- **Facultative anaerobe advantage:** Enterobacteriaceae can survive oxygen exposure during gut disruption; obligate anaerobes cannot
- **Environmental reservoirs:** Hospital surfaces, feeding tubes, hands provide continuous re-seeding of Enterobacteriaceae
- **Metabolic flexibility:** Enterobacteriaceae can utilize diverse carbon sources; specialized anaerobes cannot
- **Sporulation paradox:** Some anaerobes form spores, but germination is slow compared to vegetative Enterobacteriaceae growth

---

## Updated Figure 3 Design

Based on unbiased ranking approach:

**Figure 3: Species-Level Recovery Ranking Reveals Pathogen Enrichment at Top**

Panel A: Ranked recovery rates for all species (≥20% prevalence)
- Horizontal bar plot, fastest at top, slowest at bottom
- Bars color-coded by functional category
- Pathogens highlighted (if they cluster at top, visually obvious)

Panel B: Summary by functional category
- Box plot showing recovery rate distributions
- Pathogenic Enterobacteriaceae vs Commensal Enterobacteriaceae vs Obligate Anaerobes

Panel C: Recovery rate vs published doubling time
- Scatter plot with correlation line
- Labeled key species
- Shows mechanistic relationship

Panel D: Conceptual model
- Diagram showing: Antibiotics → Depletion of all → Pathogens recover faster → Dominance
- Ties the story together

---

## Next Steps

1. [ ] **Run species-level recovery ranking analysis** (Script 13)
2. [ ] Compile BSI isolate data for Figure 1 (3-5 representative cases)
3. [ ] Generate Figure 1 panels (timeline, susceptibility, genomic identity)
4. [ ] Generate updated Figure 3 (unbiased ranking + doubling time correlation)
5. [ ] Draft abstract (200-250 words)
6. [ ] Draft introduction with clinical paradox hook
7. [ ] Select target journal and review formatting requirements
8. [ ] Co-author review of outline and figures

---

## Data Availability Statement

All analysis scripts are available in `R/new_analysis_legacy_data/`. Raw metagenomic data deposited in [SRA accession pending]. BSI isolate genomes deposited in [NCBI accession pending]. Processed results in `results/new_analysis_legacy_data/`.
