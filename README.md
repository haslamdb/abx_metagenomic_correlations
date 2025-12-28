# Antibiotic-Microbiome Correlations Pipeline

A reproducible Snakemake pipeline for analyzing the effects of antibiotic exposure on gut microbiome composition using metagenomic sequencing data.

## Overview

This pipeline processes paired-end metagenomic sequencing data to:

1. **Quality control** and host read removal
2. **Taxonomic profiling** using Kraken2/Bracken
3. **Functional profiling** using HUMAnN 4
4. **Antibiotic resistance gene detection** using ABRicate/AMRFinderPlus
5. **Statistical analysis** correlating antibiotic exposure with microbiome changes

## Requirements

- Linux operating system (tested on Ubuntu 22.04)
- Conda/Mamba package manager
- At least 128 GB RAM (for Kraken2 with standard database)
- At least 500 GB disk space for databases
- Snakemake 7.0+

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/yourusername/abx_metagenomic_correlations.git
cd abx_metagenomic_correlations
```

### 2. Install Snakemake

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake mamba
conda activate snakemake
```

### 3. Download and set up databases

```bash
# This will download Kraken2, MetaPhlAn, HUMAnN, and CARD databases
# WARNING: This requires ~300 GB of disk space and takes several hours
bash resources/setup_databases.sh /path/to/database/directory
```

### 4. Configure the pipeline

Edit `config/config.yaml` to specify:
- Path to your sample sheet
- Path to databases
- Analysis parameters

### 5. Prepare your sample sheet

Create a tab-separated file with columns:
```
sample_id	fastq_1	fastq_2	patient_id	sample_date	patient_group
```

See `config/samples_example.tsv` for format.

### 6. Run the pipeline

```bash
# Dry run to check workflow
snakemake --dry-run

# Run with 16 cores
snakemake --cores 16 --use-conda

# Run on a cluster (SLURM example)
snakemake --profile config/slurm
```

## Pipeline Steps

```
Raw FASTQ files
       │
       ▼
┌──────────────────┐
│  Quality Control │  (fastp)
│  - Adapter trim  │
│  - Quality filter│
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│  Host Removal    │  (Bowtie2 vs GRCh38)
└────────┬─────────┘
         │
         ├─────────────────────────────────────┐
         │                                     │
         ▼                                     ▼
┌──────────────────┐                 ┌──────────────────┐
│  Taxonomic       │                 │  Functional      │
│  Classification  │                 │  Profiling       │
│  (Kraken2)       │                 │  (HUMAnN 4)      │
└────────┬─────────┘                 └────────┬─────────┘
         │                                     │
         ▼                                     ▼
┌──────────────────┐                 ┌──────────────────┐
│  Abundance       │                 │  - Gene families │
│  Estimation      │                 │  - Pathways      │
│  (Bracken)       │                 │  - Coverage      │
└────────┬─────────┘                 └────────┬─────────┘
         │                                     │
         └──────────────┬──────────────────────┘
                        │
                        ▼
              ┌──────────────────┐
              │  Assembly        │  (MEGAHIT)
              └────────┬─────────┘
                       │
                       ▼
              ┌──────────────────┐
              │  ARG Detection   │  (ABRicate/AMRFinderPlus)
              └────────┬─────────┘
                       │
                       ▼
              ┌──────────────────┐
              │  Statistical     │  (R: vegan, ALDEx2, MaAsLin2)
              │  Analysis        │
              └──────────────────┘
```

## Output Structure

```
results/
├── qc/
│   ├── {sample}_fastp.html
│   └── {sample}_fastp.json
├── host_removed/
│   └── {sample}_host_removed_{1,2}.fastq.gz
├── taxonomy/
│   ├── kraken2/
│   │   └── {sample}_kraken2_report.txt
│   └── bracken/
│       ├── {sample}_species.tsv
│       └── {sample}_genus.tsv
├── functional/
│   ├── {sample}_genefamilies.tsv
│   ├── {sample}_pathabundance.tsv
│   └── {sample}_pathcoverage.tsv
├── assembly/
│   └── {sample}/final.contigs.fa
├── arg/
│   ├── {sample}_abricate_card.tsv
│   └── {sample}_amrfinder.tsv
├── merged/
│   ├── species_abundance.tsv
│   ├── genus_abundance.tsv
│   ├── pathway_abundance.tsv
│   └── arg_abundance.tsv
├── analysis/
│   ├── alpha_diversity.tsv
│   ├── differential_abundance/
│   └── paired_analysis/
├── figures/
│   └── *.pdf
└── tables/
    └── *.tsv
```

## Clinical Data Integration

The pipeline expects antibiotic exposure data in the format:

```
patient_id	date	drug	route
```

Configure the path in `config/config.yaml`.

## Statistical Analyses

The R analysis scripts perform:

1. **Alpha diversity analysis** - Shannon, Simpson, richness by antibiotic exposure
2. **Beta diversity analysis** - PERMANOVA, Bray-Curtis ordination
3. **Differential abundance** - ALDEx2 and MaAsLin2
4. **Paired sample analysis** - Within-patient changes associated with antibiotics
5. **ARG association** - Antibiotic exposure vs resistance gene abundance

## Key Hypotheses Tested

1. Antibiotic exposure reduces alpha diversity
2. Anti-anaerobic antibiotics reduce obligate anaerobe abundance
3. Broad-spectrum antibiotics increase Enterobacteriaceae
4. Antibiotic exposure increases ARG abundance
5. Vancomycin increases vancomycin resistance genes

## Citation

If you use this pipeline, please cite:
- [Your publication here]
- Kraken2: Wood et al. (2019) Genome Biology
- HUMAnN: Beghini et al. (2021) eLife
- ALDEx2: Fernandes et al. (2014) Microbiome

## License

MIT License - see LICENSE file

## Contact

[Your contact information]
