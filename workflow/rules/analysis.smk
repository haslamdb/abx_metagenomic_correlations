"""
Statistical Analysis Rules
==========================
Diversity analysis, differential abundance, and paired sample comparisons.
"""


rule calculate_diversity:
    """Calculate alpha diversity metrics for all samples."""
    input:
        species=f"{RESULTS}/merged/species_abundance.tsv",
        metadata=config["samples"],
    output:
        alpha=f"{RESULTS}/analysis/alpha_diversity.tsv",
    log:
        f"{LOGS}/analysis/diversity.log",
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/calculate_diversity.R"


rule differential_abundance_aldex2:
    """Run ALDEx2 differential abundance analysis."""
    input:
        species=f"{RESULTS}/merged/species_abundance.tsv",
        metadata=config["samples"],
        drug_exposure=config["drug_exposure"],
    output:
        results=f"{RESULTS}/analysis/differential_abundance/aldex2_results.tsv",
        significant=f"{RESULTS}/analysis/differential_abundance/aldex2_significant.tsv",
    log:
        f"{LOGS}/analysis/aldex2.log",
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/run_aldex2.R"


rule differential_abundance_maaslin2:
    """Run MaAsLin2 differential abundance analysis."""
    input:
        species=f"{RESULTS}/merged/species_abundance.tsv",
        metadata=config["samples"],
        drug_exposure=config["drug_exposure"],
    output:
        dir=directory(f"{RESULTS}/analysis/differential_abundance/maaslin2"),
    log:
        f"{LOGS}/analysis/maaslin2.log",
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/run_maaslin2.R"


rule paired_analysis:
    """Analyze paired samples for within-patient changes."""
    input:
        species=f"{RESULTS}/merged/species_abundance.tsv",
        genus=f"{RESULTS}/merged/genus_abundance.tsv",
        metadata=config["samples"],
        drug_exposure=config["drug_exposure"],
    output:
        results=f"{RESULTS}/analysis/paired_analysis/paired_results.tsv",
        summary=f"{RESULTS}/analysis/paired_analysis/paired_summary.tsv",
    log:
        f"{LOGS}/analysis/paired.log",
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/paired_analysis.R"


rule arg_analysis:
    """Analyze ARG abundance in relation to antibiotic exposure."""
    input:
        arg=f"{RESULTS}/merged/arg_abundance.tsv",
        metadata=config["samples"],
        drug_exposure=config["drug_exposure"],
    output:
        results=f"{RESULTS}/analysis/arg_analysis/arg_exposure_association.tsv",
    log:
        f"{LOGS}/analysis/arg.log",
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/arg_analysis.R"


rule plot_alpha_diversity:
    """Generate alpha diversity boxplots."""
    input:
        alpha=f"{RESULTS}/analysis/alpha_diversity.tsv",
        metadata=config["samples"],
        drug_exposure=config["drug_exposure"],
    output:
        f"{RESULTS}/figures/alpha_diversity_boxplot.pdf",
    log:
        f"{LOGS}/plots/alpha_diversity.log",
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/plot_alpha_diversity.R"


rule plot_ordination:
    """Generate PCoA ordination plots."""
    input:
        species=f"{RESULTS}/merged/species_abundance.tsv",
        metadata=config["samples"],
        drug_exposure=config["drug_exposure"],
    output:
        f"{RESULTS}/figures/pcoa_bray_curtis.pdf",
    log:
        f"{LOGS}/plots/ordination.log",
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/plot_ordination.R"


rule plot_differential_abundance:
    """Generate volcano plots for differential abundance."""
    input:
        results=f"{RESULTS}/analysis/differential_abundance/aldex2_results.tsv",
    output:
        f"{RESULTS}/figures/differential_abundance_volcano.pdf",
    log:
        f"{LOGS}/plots/volcano.log",
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/plot_volcano.R"


rule generate_report:
    """Generate final analysis report."""
    input:
        alpha=f"{RESULTS}/analysis/alpha_diversity.tsv",
        aldex2=f"{RESULTS}/analysis/differential_abundance/aldex2_results.tsv",
        paired=f"{RESULTS}/analysis/paired_analysis/paired_results.tsv",
        figures=expand(
            f"{RESULTS}/figures/{{fig}}.pdf",
            fig=["alpha_diversity_boxplot", "pcoa_bray_curtis", "differential_abundance_volcano"]
        ),
    output:
        f"{RESULTS}/analysis/final_report.html",
    log:
        f"{LOGS}/analysis/report.log",
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/generate_report.Rmd"
