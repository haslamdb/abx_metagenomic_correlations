"""
Taxonomic Classification Rules
==============================
Kraken2 classification and Bracken abundance estimation.
"""


rule kraken2:
    """Classify reads using Kraken2."""
    input:
        r1=f"{RESULTS}/host_removed/{{sample}}_host_removed_1.fastq.gz",
        r2=f"{RESULTS}/host_removed/{{sample}}_host_removed_2.fastq.gz",
        db=config["databases"]["kraken2"],
    output:
        report=f"{RESULTS}/taxonomy/kraken2/{{sample}}_kraken2_report.txt",
        output=f"{RESULTS}/taxonomy/kraken2/{{sample}}_kraken2_output.txt",
    log:
        f"{LOGS}/kraken2/{{sample}}.log",
    threads: 16
    conda:
        "../envs/taxonomy.yaml"
    params:
        confidence=config.get("kraken2_confidence", 0.0),
    shell:
        """
        kraken2 \
            --db {input.db} \
            --paired \
            --gzip-compressed \
            --threads {threads} \
            --confidence {params.confidence} \
            --report {output.report} \
            --output {output.output} \
            {input.r1} {input.r2} \
            2> {log}
        """


rule bracken_species:
    """Estimate species-level abundances using Bracken."""
    input:
        report=f"{RESULTS}/taxonomy/kraken2/{{sample}}_kraken2_report.txt",
        db=config["databases"]["kraken2"],
    output:
        bracken=f"{RESULTS}/taxonomy/bracken/{{sample}}_species.tsv",
        report=f"{RESULTS}/taxonomy/bracken/{{sample}}_species_report.txt",
    log:
        f"{LOGS}/bracken/{{sample}}_species.log",
    conda:
        "../envs/taxonomy.yaml"
    params:
        read_len=config.get("bracken_read_length", 150),
        threshold=config.get("bracken_threshold", 10),
    shell:
        """
        bracken \
            -d {input.db} \
            -i {input.report} \
            -o {output.bracken} \
            -w {output.report} \
            -r {params.read_len} \
            -l S \
            -t {params.threshold} \
            2> {log}
        """


rule bracken_genus:
    """Estimate genus-level abundances using Bracken."""
    input:
        report=f"{RESULTS}/taxonomy/kraken2/{{sample}}_kraken2_report.txt",
        db=config["databases"]["kraken2"],
    output:
        bracken=f"{RESULTS}/taxonomy/bracken/{{sample}}_genus.tsv",
        report=f"{RESULTS}/taxonomy/bracken/{{sample}}_genus_report.txt",
    log:
        f"{LOGS}/bracken/{{sample}}_genus.log",
    conda:
        "../envs/taxonomy.yaml"
    params:
        read_len=config.get("bracken_read_length", 150),
        threshold=config.get("bracken_threshold", 10),
    shell:
        """
        bracken \
            -d {input.db} \
            -i {input.report} \
            -o {output.bracken} \
            -w {output.report} \
            -r {params.read_len} \
            -l G \
            -t {params.threshold} \
            2> {log}
        """


rule bracken_phylum:
    """Estimate phylum-level abundances using Bracken."""
    input:
        report=f"{RESULTS}/taxonomy/kraken2/{{sample}}_kraken2_report.txt",
        db=config["databases"]["kraken2"],
    output:
        bracken=f"{RESULTS}/taxonomy/bracken/{{sample}}_phylum.tsv",
        report=f"{RESULTS}/taxonomy/bracken/{{sample}}_phylum_report.txt",
    log:
        f"{LOGS}/bracken/{{sample}}_phylum.log",
    conda:
        "../envs/taxonomy.yaml"
    params:
        read_len=config.get("bracken_read_length", 150),
        threshold=config.get("bracken_threshold", 10),
    shell:
        """
        bracken \
            -d {input.db} \
            -i {input.report} \
            -o {output.bracken} \
            -w {output.report} \
            -r {params.read_len} \
            -l P \
            -t {params.threshold} \
            2> {log}
        """
