"""
Antibiotic Resistance Gene Detection Rules
===========================================
ABRicate and AMRFinderPlus for ARG identification.
"""


rule abricate_card:
    """Run ABRicate with CARD database."""
    input:
        contigs=f"{RESULTS}/assembly/{{sample}}/final.contigs.fa",
    output:
        f"{RESULTS}/arg/{{sample}}_abricate_card.tsv",
    log:
        f"{LOGS}/abricate/{{sample}}_card.log",
    threads: 8
    conda:
        "../envs/arg.yaml"
    params:
        minid=config.get("abricate_minid", 80),
        mincov=config.get("abricate_mincov", 60),
    shell:
        """
        abricate \
            --db card \
            --minid {params.minid} \
            --mincov {params.mincov} \
            --threads {threads} \
            {input.contigs} \
            > {output} \
            2> {log}
        """


rule abricate_resfinder:
    """Run ABRicate with ResFinder database."""
    input:
        contigs=f"{RESULTS}/assembly/{{sample}}/final.contigs.fa",
    output:
        f"{RESULTS}/arg/{{sample}}_abricate_resfinder.tsv",
    log:
        f"{LOGS}/abricate/{{sample}}_resfinder.log",
    threads: 8
    conda:
        "../envs/arg.yaml"
    params:
        minid=config.get("abricate_minid", 80),
        mincov=config.get("abricate_mincov", 60),
    shell:
        """
        abricate \
            --db resfinder \
            --minid {params.minid} \
            --mincov {params.mincov} \
            --threads {threads} \
            {input.contigs} \
            > {output} \
            2> {log}
        """


rule abricate_ncbi:
    """Run ABRicate with NCBI AMR database."""
    input:
        contigs=f"{RESULTS}/assembly/{{sample}}/final.contigs.fa",
    output:
        f"{RESULTS}/arg/{{sample}}_abricate_ncbi.tsv",
    log:
        f"{LOGS}/abricate/{{sample}}_ncbi.log",
    threads: 8
    conda:
        "../envs/arg.yaml"
    params:
        minid=config.get("abricate_minid", 80),
        mincov=config.get("abricate_mincov", 60),
    shell:
        """
        abricate \
            --db ncbi \
            --minid {params.minid} \
            --mincov {params.mincov} \
            --threads {threads} \
            {input.contigs} \
            > {output} \
            2> {log}
        """


rule amrfinderplus:
    """Run AMRFinderPlus for comprehensive ARG detection."""
    input:
        contigs=f"{RESULTS}/assembly/{{sample}}/final.contigs.fa",
    output:
        f"{RESULTS}/arg/{{sample}}_amrfinder.tsv",
    log:
        f"{LOGS}/amrfinder/{{sample}}.log",
    threads: 8
    conda:
        "../envs/arg.yaml"
    shell:
        """
        amrfinder \
            -n {input.contigs} \
            -o {output} \
            --threads {threads} \
            --plus \
            2> {log}
        """


rule summarize_arg:
    """Create per-sample ARG summary."""
    input:
        card=f"{RESULTS}/arg/{{sample}}_abricate_card.tsv",
        amrfinder=f"{RESULTS}/arg/{{sample}}_amrfinder.tsv",
    output:
        f"{RESULTS}/arg/{{sample}}_arg_summary.tsv",
    run:
        import pandas as pd

        # Load CARD results
        card_df = pd.read_csv(input.card, sep="\t")
        card_genes = card_df["GENE"].nunique() if len(card_df) > 0 else 0

        # Load AMRFinder results
        amr_df = pd.read_csv(input.amrfinder, sep="\t")
        amr_genes = len(amr_df) if len(amr_df) > 0 else 0

        # Write summary
        summary = pd.DataFrame({
            "sample": [wildcards.sample],
            "card_genes": [card_genes],
            "amrfinder_genes": [amr_genes],
        })
        summary.to_csv(output[0], sep="\t", index=False)
