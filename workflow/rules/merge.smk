"""
Merge Rules
===========
Combine per-sample outputs into abundance matrices.
"""


rule merge_species:
    """Merge Bracken species abundance across all samples."""
    input:
        expand(f"{RESULTS}/taxonomy/bracken/{{sample}}_species.tsv", sample=SAMPLES),
    output:
        f"{RESULTS}/merged/species_abundance.tsv",
    log:
        f"{LOGS}/merge/species.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge_bracken.py"


rule merge_genus:
    """Merge Bracken genus abundance across all samples."""
    input:
        expand(f"{RESULTS}/taxonomy/bracken/{{sample}}_genus.tsv", sample=SAMPLES),
    output:
        f"{RESULTS}/merged/genus_abundance.tsv",
    log:
        f"{LOGS}/merge/genus.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge_bracken.py"


rule merge_pathways:
    """Merge HUMAnN pathway abundances across all samples."""
    input:
        expand(f"{RESULTS}/functional/{{sample}}/{{sample}}_pathabundance_cpm.tsv", sample=SAMPLES),
    output:
        f"{RESULTS}/merged/pathway_abundance.tsv",
    log:
        f"{LOGS}/merge/pathways.log",
    conda:
        "../envs/functional.yaml"
    params:
        input_dir=f"{RESULTS}/functional",
    shell:
        """
        humann_join_tables \
            --input {params.input_dir} \
            --output {output} \
            --file_name pathabundance_cpm \
            2> {log}
        """


rule merge_genefamilies:
    """Merge HUMAnN gene families across all samples."""
    input:
        expand(f"{RESULTS}/functional/{{sample}}/{{sample}}_genefamilies_cpm.tsv", sample=SAMPLES),
    output:
        f"{RESULTS}/merged/genefamilies_abundance.tsv",
    log:
        f"{LOGS}/merge/genefamilies.log",
    conda:
        "../envs/functional.yaml"
    params:
        input_dir=f"{RESULTS}/functional",
    shell:
        """
        humann_join_tables \
            --input {params.input_dir} \
            --output {output} \
            --file_name genefamilies_cpm \
            2> {log}
        """


rule merge_ko:
    """Merge KEGG ortholog abundances across all samples."""
    input:
        expand(f"{RESULTS}/functional/{{sample}}/{{sample}}_ko_cpm.tsv", sample=SAMPLES),
    output:
        f"{RESULTS}/merged/ko_abundance.tsv",
    log:
        f"{LOGS}/merge/ko.log",
    conda:
        "../envs/functional.yaml"
    params:
        input_dir=f"{RESULTS}/functional",
    shell:
        """
        humann_join_tables \
            --input {params.input_dir} \
            --output {output} \
            --file_name ko_cpm \
            2> {log}
        """


rule merge_arg:
    """Merge ARG results across all samples."""
    input:
        card=expand(f"{RESULTS}/arg/{{sample}}_abricate_card.tsv", sample=SAMPLES),
        amrfinder=expand(f"{RESULTS}/arg/{{sample}}_amrfinder.tsv", sample=SAMPLES),
    output:
        f"{RESULTS}/merged/arg_abundance.tsv",
    log:
        f"{LOGS}/merge/arg.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge_arg.py"
