# Antibiotic Effects on Microbiome: Revised Analysis Plan

## Background

Previous analysis attempted to correlate antibiotic exposure with microbiome changes using
Pearson correlations between rate of drug administration and rate of species/genus change.
This approach did not reveal expected patterns (e.g., anaerobe depletion with metronidazole).

This revised plan addresses methodological limitations and leverages updated bioinformatics tools.

---

## Phase 1: Sample Re-processing with Updated Tools

### 1.1 Raw Data Inventory

**Action**: Catalog all raw FASTQ files
- Location: [Need to confirm path to raw sequences]
- Create manifest linking sequence files to sample metadata
- Verify paired-end status and read lengths

### 1.2 Quality Control

**Tools**: FastQC, MultiQC, fastp

```bash
# Per-sample QC
fastp \
  --in1 {sample}_R1.fastq.gz \
  --in2 {sample}_R2.fastq.gz \
  --out1 {sample}_R1.clean.fastq.gz \
  --out2 {sample}_R2.clean.fastq.gz \
  --html {sample}_fastp.html \
  --json {sample}_fastp.json \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --thread 8
```

### 1.3 Host Read Removal

**Tool**: Bowtie2 against human genome (GRCh38)

```bash
bowtie2 \
  -x /path/to/GRCh38/bowtie2_index \
  -1 {sample}_R1.clean.fastq.gz \
  -2 {sample}_R2.clean.fastq.gz \
  --un-conc-gz {sample}_host_removed_%.fastq.gz \
  -p 16 \
  -S /dev/null
```

### 1.4 Taxonomic Classification

**Tools**: Kraken2 (v2.1.3+) + Bracken (latest)

**Database**: Use updated database (Standard-16 or PlusPF for fungi/protozoa)

```bash
# Kraken2 classification
kraken2 \
  --db /path/to/kraken2_db \
  --paired \
  --gzip-compressed \
  --threads 16 \
  --report {sample}_kraken2_report.txt \
  --output {sample}_kraken2_output.txt \
  {sample}_host_removed_1.fastq.gz \
  {sample}_host_removed_2.fastq.gz

# Bracken abundance estimation (species level)
bracken \
  -d /path/to/kraken2_db \
  -i {sample}_kraken2_report.txt \
  -o {sample}_bracken_species.txt \
  -w {sample}_bracken_species_report.txt \
  -r 150 \
  -l S \
  -t 10

# Bracken abundance estimation (genus level)
bracken \
  -d /path/to/kraken2_db \
  -i {sample}_kraken2_report.txt \
  -o {sample}_bracken_genus.txt \
  -w {sample}_bracken_genus_report.txt \
  -r 150 \
  -l G \
  -t 10
```

### 1.5 Functional Profiling

**Tool**: HUMAnN 4 (latest version)

```bash
# Concatenate paired reads for HUMAnN
cat {sample}_host_removed_1.fastq.gz {sample}_host_removed_2.fastq.gz > {sample}_concat.fastq.gz

# Run HUMAnN 4
humann \
  --input {sample}_concat.fastq.gz \
  --output {sample}_humann_output \
  --threads 16 \
  --metaphlan-options "--bowtie2db /path/to/metaphlan_db"

# Outputs:
#   - {sample}_genefamilies.tsv (UniRef90 gene families)
#   - {sample}_pathabundance.tsv (MetaCyc pathways)
#   - {sample}_pathcoverage.tsv (pathway coverage)
```

**Post-processing**:
```bash
# Normalize to copies per million (CPM)
humann_renorm_table --input {sample}_genefamilies.tsv --output {sample}_genefamilies_cpm.tsv --units cpm
humann_renorm_table --input {sample}_pathabundance.tsv --output {sample}_pathabundance_cpm.tsv --units cpm

# Regroup to useful categories (e.g., KEGG orthologs, EC numbers)
humann_regroup_table --input {sample}_genefamilies_cpm.tsv --output {sample}_ko_cpm.tsv --groups uniref90_ko
```

### 1.6 Antibiotic Resistance Gene Analysis

**Tools**: ABRicate (CARD, ResFinder, NCBI AMR databases) or AMRFinderPlus

#### Option A: ABRicate on assembled contigs
```bash
# First assemble (if not already done)
megahit \
  -1 {sample}_host_removed_1.fastq.gz \
  -2 {sample}_host_removed_2.fastq.gz \
  -o {sample}_megahit_assembly \
  -t 16

# Run ABRicate against multiple databases
for db in card resfinder ncbi argannot; do
  abricate \
    --db $db \
    --minid 80 \
    --mincov 60 \
    {sample}_megahit_assembly/final.contigs.fa \
    > {sample}_abricate_${db}.tsv
done
```

#### Option B: Read-based ARG quantification with ShortBRED or CARD-RGI
```bash
# ShortBRED against CARD markers
shortbred_quantify.py \
  --markers CARD_shortbred_markers.fa \
  --wgs {sample}_concat.fastq.gz \
  --results {sample}_shortbred_card.tsv \
  --tmp {sample}_shortbred_tmp
```

#### Option C: AMRFinderPlus on assembled contigs
```bash
amrfinder \
  -n {sample}_megahit_assembly/final.contigs.fa \
  -o {sample}_amrfinder.tsv \
  --threads 8
```

---

## Phase 2: Data Integration and Preprocessing

### 2.1 Merge Sample Tables

Create combined abundance tables:
- `species_abundance_matrix.tsv` (samples × species)
- `genus_abundance_matrix.tsv` (samples × genera)
- `pathway_abundance_matrix.tsv` (samples × MetaCyc pathways)
- `ko_abundance_matrix.tsv` (samples × KEGG orthologs)
- `arg_abundance_matrix.tsv` (samples × ARG families)

### 2.2 Sample Metadata Integration

Merge with clinical data:
- Patient ID, MRN (de-identified)
- Sample date
- Patient group (BMT, IBD, IF, LvTx, PICU, SB, Urology)
- Age at sample
- **Antibiotic exposure** (detailed below)

### 2.3 Enhanced Antibiotic Exposure Data

For each sample, calculate:

| Variable | Description |
|----------|-------------|
| `abx_any_7d` | Any antibiotic in 7 days before sample |
| `abx_any_14d` | Any antibiotic in 14 days before sample |
| `abx_days_14d` | Number of antibiotic-days in prior 14 days |
| `abx_anaerobic_7d` | Anti-anaerobic antibiotic in 7 days |
| `abx_gram_pos_7d` | Anti-gram-positive antibiotic in 7 days |
| `abx_broad_7d` | Broad-spectrum antibiotic in 7 days |
| `pip_tazo_7d` | Pip/Tazo specifically in 7 days |
| `meropenem_7d` | Meropenem specifically in 7 days |
| `vancomycin_7d` | Vancomycin specifically in 7 days |
| `metronidazole_7d` | Metronidazole specifically in 7 days |

**Antibiotic Spectrum Classification**:
```
Anti-anaerobic: Metronidazole, Piperacillin/Tazobactam, Meropenem, Ertapenem,
                Clindamycin, Amoxicillin/Clavulanate, Ampicillin/Sulbactam

Anti-gram-positive: Vancomycin, Daptomycin, Linezolid, Ceftaroline

Broad-spectrum: Meropenem, Piperacillin/Tazobactam, Cefepime, Imipenem

Anti-Pseudomonal: Cefepime, Piperacillin/Tazobactam, Meropenem, Tobramycin,
                  Ciprofloxacin, Aztreonam
```

### 2.4 Create Functional Taxonomic Groups

Aggregate species into biologically meaningful groups:

| Group | Taxa | Rationale |
|-------|------|-----------|
| Obligate anaerobes | Bacteroides, Parabacteroides, Prevotella, Clostridium, Faecalibacterium, Blautia, Roseburia, Ruminococcus | Target of anti-anaerobic antibiotics |
| Enterobacteriaceae | Escherichia, Klebsiella, Enterobacter, Citrobacter, Serratia | Often expand after broad-spectrum Abx |
| Enterococcus | Enterococcus | Expands with vancomycin resistance |
| Lactobacillales | Lactobacillus, Streptococcus, Enterococcus | Gram-positive fermenters |
| Proteobacteria | All Proteobacteria | Dysbiosis marker |
| Firmicutes | All Firmicutes | Major phylum |
| Bacteroidetes | All Bacteroidetes | Major phylum |

---

## Phase 3: Statistical Analysis

### 3.1 Diversity Analyses

#### Alpha Diversity
```r
library(vegan)
library(phyloseq)

# Calculate per-sample metrics
alpha_metrics <- data.frame(
  sample = rownames(species_matrix),
  shannon = diversity(species_matrix, index = "shannon"),
  simpson = diversity(species_matrix, index = "simpson"),
  richness = specnumber(species_matrix),
  evenness = diversity(species_matrix) / log(specnumber(species_matrix))
)

# Mixed-effects model: Does recent antibiotic exposure reduce diversity?
library(lme4)
library(lmerTest)

model_shannon <- lmer(
  shannon ~ abx_any_7d + patient_group + age_years + (1|patient_id),
  data = merged_data
)
summary(model_shannon)

# Test specific antibiotic classes
model_anaerobic <- lmer(
  shannon ~ abx_anaerobic_7d + abx_broad_7d + patient_group + (1|patient_id),
  data = merged_data
)
```

#### Beta Diversity
```r
# Calculate Bray-Curtis distances
bray_dist <- vegdist(species_matrix, method = "bray")

# PERMANOVA: Does antibiotic exposure explain community variation?
adonis2(
  bray_dist ~ abx_any_7d + patient_group + age_years,
  data = metadata,
  permutations = 999,
  by = "margin"
)

# For paired samples: Calculate within-patient distance
# Compare distance when antibiotics given vs. not given
paired_distances <- calculate_paired_bray(sample_pairs, abx_status)
wilcox.test(distance ~ abx_between, data = paired_distances)
```

### 3.2 Differential Abundance Analysis

#### Using ALDEx2 (compositional-aware)
```r
library(ALDEx2)

# Compare samples with vs without recent antibiotic exposure
abx_status <- ifelse(metadata$abx_any_7d == TRUE, "exposed", "unexposed")

aldex_results <- aldex(
  reads = t(species_matrix),  # taxa as rows
  conditions = abx_status,
  mc.samples = 128,
  test = "t",
  effect = TRUE,
  include.sample.summary = FALSE
)

# Filter significant results
sig_taxa <- aldex_results[aldex_results$wi.eBH < 0.1, ]
```

#### Using MaAsLin2 (multivariate)
```r
library(Maaslin2)

# Full model with covariates
maaslin_results <- Maaslin2(
  input_data = species_matrix,
  input_metadata = metadata,
  output = "maaslin2_output",
  fixed_effects = c("abx_anaerobic_7d", "abx_broad_7d", "patient_group", "age_years"),
  random_effects = c("patient_id"),
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  min_abundance = 0.0001,
  min_prevalence = 0.1
)
```

### 3.3 Paired Sample Analysis (Primary Analysis)

This is the key improvement - use proper paired design:

```r
# For each patient with paired samples, calculate:
# 1. Change in diversity
# 2. Change in functional group abundances
# 3. Whether antibiotics were given between samples

paired_analysis <- function(sample_pairs, species_matrix, drug_exposure) {

  results <- data.frame()

  for (i in 1:nrow(sample_pairs)) {
    sample1 <- sample_pairs$sample1[i]
    sample2 <- sample_pairs$sample2[i]
    patient <- sample_pairs$patient_id[i]

    # Calculate microbiome changes
    delta_shannon <- shannon[sample2] - shannon[sample1]
    delta_anaerobes <- log2((anaerobe_sum[sample2] + 1) / (anaerobe_sum[sample1] + 1))
    delta_enterobact <- log2((enterobact_sum[sample2] + 1) / (enterobact_sum[sample1] + 1))
    bray_distance <- vegdist(rbind(species_matrix[sample1,], species_matrix[sample2,]))[1]

    # Get antibiotic exposure between samples
    date1 <- sample_pairs$date1[i]
    date2 <- sample_pairs$date2[i]
    abx_between <- get_abx_exposure(patient, date1, date2, drug_exposure)

    results <- rbind(results, data.frame(
      patient_id = patient,
      interval_days = as.numeric(date2 - date1),
      delta_shannon = delta_shannon,
      delta_anaerobes = delta_anaerobes,
      delta_enterobact = delta_enterobact,
      bray_distance = bray_distance,
      any_abx = abx_between$any,
      anaerobic_abx = abx_between$anaerobic,
      broad_spectrum = abx_between$broad,
      abx_days = abx_between$total_days
    ))
  }

  return(results)
}

# Statistical tests
# H1: Anti-anaerobic antibiotics decrease anaerobe abundance
wilcox.test(delta_anaerobes ~ anaerobic_abx, data = paired_results)

# H2: Antibiotics increase community disruption (Bray-Curtis distance)
wilcox.test(bray_distance ~ any_abx, data = paired_results)

# H3: Dose-response relationship
cor.test(paired_results$abx_days, paired_results$delta_shannon, method = "spearman")

# Mixed-effects for proper inference
library(lme4)
model <- lmer(delta_anaerobes ~ anaerobic_abx + interval_days + (1|patient_id),
              data = paired_results)
```

### 3.4 Functional Pathway Analysis

```r
# Test specific pathways of interest
# E.g., butyrate production, amino acid biosynthesis

# HUMAnN4 pathway abundance
pathway_results <- Maaslin2(
  input_data = pathway_matrix,
  input_metadata = metadata,
  output = "pathway_maaslin2",
  fixed_effects = c("abx_any_7d", "patient_group"),
  random_effects = c("patient_id"),
  normalization = "TSS",
  transform = "LOG"
)

# Focus on SCFA production pathways
scfa_pathways <- c(
  "PWY-5022",  # 4-aminobutanoate degradation
  "PWY-5676",  # acetyl-CoA fermentation to butanoate II
  "P163-PWY"   # lysine fermentation to acetate and butanoate
)
```

### 3.5 Antibiotic Resistance Gene Analysis

```r
# Does antibiotic exposure increase ARG abundance?

# Aggregate ARGs by class
arg_classes <- list(
  beta_lactamase = c("CTX-M", "TEM", "SHV", "OXA", "KPC", "NDM"),
  aminoglycoside = c("aac", "aph", "ant"),
  vancomycin = c("vanA", "vanB", "vanC"),
  fluoroquinolone = c("qnr", "gyrA", "parC")
)

# Test association
model_arg <- lmer(
  total_arg_rpkm ~ abx_any_14d + patient_group + (1|patient_id),
  data = merged_data
)

# Specific hypothesis: Vancomycin increases vancomycin resistance genes
model_van <- lmer(
  vancomycin_arg ~ vancomycin_14d + patient_group + (1|patient_id),
  data = merged_data
)
```

---

## Phase 4: Visualization and Reporting

### 4.1 Key Figures

1. **Alpha diversity by antibiotic exposure** (boxplots)
2. **PCoA/NMDS colored by antibiotic status** (ordination)
3. **Heatmap of differentially abundant taxa** (clustered heatmap)
4. **Volcano plots** from MaAsLin2/ALDEx2
5. **Paired sample plots** (before/after with lines connecting paired samples)
6. **Functional group trajectories** (anaerobes, Enterobacteriaceae over time)
7. **ARG abundance by exposure** (bar charts)

### 4.2 Summary Tables

1. Sample demographics and clinical characteristics
2. Antibiotic exposure summary
3. Differential abundance results (significant taxa)
4. Model coefficients and confidence intervals

---

## Phase 5: Pipeline Implementation

### 5.1 Suggested Workflow Manager

Use **Snakemake** or **Nextflow** for reproducibility:

```
workflow/
├── Snakefile
├── config.yaml
├── envs/
│   ├── kraken2.yaml
│   ├── humann.yaml
│   └── r_analysis.yaml
├── rules/
│   ├── qc.smk
│   ├── host_removal.smk
│   ├── taxonomy.smk
│   ├── functional.smk
│   └── arg.smk
└── scripts/
    ├── merge_bracken.py
    ├── calculate_exposure.R
    └── run_analysis.R
```

### 5.2 Computational Requirements

| Step | CPU cores | RAM | Time (per sample) |
|------|-----------|-----|-------------------|
| QC (fastp) | 8 | 4 GB | 5 min |
| Host removal | 16 | 8 GB | 15 min |
| Kraken2 | 16 | 64 GB | 10 min |
| Bracken | 1 | 4 GB | 1 min |
| HUMAnN4 | 16 | 32 GB | 2-4 hours |
| Assembly | 16 | 32 GB | 30-60 min |
| ARG detection | 8 | 16 GB | 10 min |

**Total estimate**: 4-6 hours per sample (dominated by HUMAnN4)

For ~900 samples with parallelization on a cluster: 1-2 weeks

---

## Key Hypotheses to Test

1. **H1**: Recent antibiotic exposure (7-14 days) reduces alpha diversity
2. **H2**: Anti-anaerobic antibiotics specifically reduce obligate anaerobe abundance
3. **H3**: Broad-spectrum antibiotics increase Enterobacteriaceae relative abundance
4. **H4**: Antibiotic exposure increases within-patient microbiome instability (Bray-Curtis distance between paired samples)
5. **H5**: Antibiotic exposure increases ARG abundance
6. **H6**: Vancomycin exposure specifically increases vancomycin resistance genes and Enterococcus
7. **H7**: Butyrate-producing pathway abundance decreases with antibiotic exposure

---

## Next Steps

1. [ ] Locate and inventory all raw FASTQ files
2. [ ] Update/install Kraken2, Bracken, HUMAnN4, and ARG tools
3. [ ] Download latest databases (Kraken2 Standard, MetaPhlAn, ChocoPhlAn, UniRef90, CARD)
4. [ ] Create Snakemake/Nextflow pipeline
5. [ ] Run QC and host removal on all samples
6. [ ] Run taxonomic and functional profiling
7. [ ] Merge results and integrate with clinical data
8. [ ] Execute statistical analysis plan
9. [ ] Generate figures and tables
10. [ ] Interpret and write up results
