"""
Functional Profiling Rules
==========================
HUMAnN 4 for pathway and gene family quantification.
"""


rule concat_reads:
    """Concatenate paired reads for HUMAnN input."""
    input:
        r1=f"{RESULTS}/host_removed/{{sample}}_host_removed_1.fastq.gz",
        r2=f"{RESULTS}/host_removed/{{sample}}_host_removed_2.fastq.gz",
    output:
        temp(f"{RESULTS}/functional/{{sample}}_concat.fastq.gz"),
    shell:
        """
        cat {input.r1} {input.r2} > {output}
        """


rule humann:
    """Run HUMAnN 4 for functional profiling."""
    input:
        reads=f"{RESULTS}/functional/{{sample}}_concat.fastq.gz",
        nucleotide_db=config["databases"]["humann_nucleotide"],
        protein_db=config["databases"]["humann_protein"],
    output:
        genefamilies=f"{RESULTS}/functional/{{sample}}/{{sample}}_genefamilies.tsv",
        pathabundance=f"{RESULTS}/functional/{{sample}}/{{sample}}_pathabundance.tsv",
        pathcoverage=f"{RESULTS}/functional/{{sample}}/{{sample}}_pathcoverage.tsv",
    log:
        f"{LOGS}/humann/{{sample}}.log",
    threads: 16
    conda:
        "../envs/functional.yaml"
    params:
        output_dir=lambda wc: f"{RESULTS}/functional/{wc.sample}",
        metaphlan_db=config["databases"]["metaphlan"],
    shell:
        """
        humann \
            --input {input.reads} \
            --output {params.output_dir} \
            --nucleotide-database {input.nucleotide_db} \
            --protein-database {input.protein_db} \
            --metaphlan-options "--bowtie2db {params.metaphlan_db}" \
            --threads {threads} \
            --output-basename {wildcards.sample} \
            2> {log}
        """


rule humann_renorm:
    """Normalize HUMAnN output to copies per million (CPM)."""
    input:
        genefamilies=f"{RESULTS}/functional/{{sample}}/{{sample}}_genefamilies.tsv",
        pathabundance=f"{RESULTS}/functional/{{sample}}/{{sample}}_pathabundance.tsv",
    output:
        genefamilies_cpm=f"{RESULTS}/functional/{{sample}}/{{sample}}_genefamilies_cpm.tsv",
        pathabundance_cpm=f"{RESULTS}/functional/{{sample}}/{{sample}}_pathabundance_cpm.tsv",
    log:
        f"{LOGS}/humann/{{sample}}_renorm.log",
    conda:
        "../envs/functional.yaml"
    shell:
        """
        humann_renorm_table \
            --input {input.genefamilies} \
            --output {output.genefamilies_cpm} \
            --units cpm \
            2>> {log}

        humann_renorm_table \
            --input {input.pathabundance} \
            --output {output.pathabundance_cpm} \
            --units cpm \
            2>> {log}
        """


rule humann_regroup:
    """Regroup gene families to KEGG orthologs."""
    input:
        genefamilies=f"{RESULTS}/functional/{{sample}}/{{sample}}_genefamilies_cpm.tsv",
    output:
        ko=f"{RESULTS}/functional/{{sample}}/{{sample}}_ko_cpm.tsv",
        ec=f"{RESULTS}/functional/{{sample}}/{{sample}}_ec_cpm.tsv",
    log:
        f"{LOGS}/humann/{{sample}}_regroup.log",
    conda:
        "../envs/functional.yaml"
    shell:
        """
        humann_regroup_table \
            --input {input.genefamilies} \
            --output {output.ko} \
            --groups uniref90_ko \
            2>> {log}

        humann_regroup_table \
            --input {input.genefamilies} \
            --output {output.ec} \
            --groups uniref90_level4ec \
            2>> {log}
        """
