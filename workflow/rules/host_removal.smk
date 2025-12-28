"""
Host Read Removal Rules
=======================
Remove human reads using Bowtie2 alignment to GRCh38.
"""


rule bowtie2_host_removal:
    """Align reads to human genome and extract unaligned (non-host) reads."""
    input:
        r1=f"{RESULTS}/qc/{{sample}}_R1.clean.fastq.gz",
        r2=f"{RESULTS}/qc/{{sample}}_R2.clean.fastq.gz",
        index=config["databases"]["host_index"],
    output:
        r1=f"{RESULTS}/host_removed/{{sample}}_host_removed_1.fastq.gz",
        r2=f"{RESULTS}/host_removed/{{sample}}_host_removed_2.fastq.gz",
    log:
        f"{LOGS}/bowtie2/{{sample}}.log",
    threads: 16
    conda:
        "../envs/alignment.yaml"
    params:
        index_base=lambda wc, input: input.index.replace(".1.bt2", ""),
    shell:
        """
        bowtie2 \
            -x {params.index_base} \
            -1 {input.r1} \
            -2 {input.r2} \
            --un-conc-gz {RESULTS}/host_removed/{wildcards.sample}_host_removed_%.fastq.gz \
            --threads {threads} \
            -S /dev/null \
            2> {log}
        """


rule count_host_reads:
    """Count host vs non-host reads for QC purposes."""
    input:
        original_r1=lambda wc: get_fastq(wc, 1),
        filtered_r1=f"{RESULTS}/host_removed/{{sample}}_host_removed_1.fastq.gz",
    output:
        f"{RESULTS}/host_removed/{{sample}}_host_stats.txt",
    conda:
        "../envs/qc.yaml"
    shell:
        """
        orig=$(zcat {input.original_r1} | wc -l | awk '{{print $1/4}}')
        filt=$(zcat {input.filtered_r1} | wc -l | awk '{{print $1/4}}')
        host=$((orig - filt))
        pct=$(echo "scale=2; $host * 100 / $orig" | bc)

        echo -e "sample\toriginal_reads\thost_reads\tnon_host_reads\tpercent_host" > {output}
        echo -e "{wildcards.sample}\t$orig\t$host\t$filt\t$pct" >> {output}
        """
