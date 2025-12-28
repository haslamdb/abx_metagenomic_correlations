"""
Assembly Rules
==============
Metagenomic assembly using MEGAHIT for downstream ARG analysis.
"""


rule megahit:
    """Assemble reads using MEGAHIT."""
    input:
        r1=f"{RESULTS}/host_removed/{{sample}}_host_removed_1.fastq.gz",
        r2=f"{RESULTS}/host_removed/{{sample}}_host_removed_2.fastq.gz",
    output:
        contigs=f"{RESULTS}/assembly/{{sample}}/final.contigs.fa",
    log:
        f"{LOGS}/megahit/{{sample}}.log",
    threads: 16
    conda:
        "../envs/assembly.yaml"
    params:
        out_dir=lambda wc: f"{RESULTS}/assembly/{wc.sample}",
        min_contig=config.get("megahit_min_contig", 500),
    shell:
        """
        # Remove output directory if exists (MEGAHIT requirement)
        rm -rf {params.out_dir}

        megahit \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.out_dir} \
            --min-contig-len {params.min_contig} \
            -t {threads} \
            2> {log}
        """


rule assembly_stats:
    """Calculate assembly statistics."""
    input:
        contigs=f"{RESULTS}/assembly/{{sample}}/final.contigs.fa",
    output:
        stats=f"{RESULTS}/assembly/{{sample}}/assembly_stats.txt",
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        # Calculate basic stats
        total_contigs=$(grep -c ">" {input.contigs})
        total_bp=$(grep -v ">" {input.contigs} | tr -d '\n' | wc -c)
        n50=$(seqkit stats -T {input.contigs} | tail -1 | cut -f7)

        echo -e "sample\ttotal_contigs\ttotal_bp\tN50" > {output.stats}
        echo -e "{wildcards.sample}\t$total_contigs\t$total_bp\t$n50" >> {output.stats}
        """
