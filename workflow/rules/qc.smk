"""
Quality Control Rules
=====================
Adapter trimming and quality filtering using fastp.
"""


rule fastp:
    """Run fastp for adapter trimming and quality filtering."""
    input:
        r1=lambda wc: get_fastq(wc, 1),
        r2=lambda wc: get_fastq(wc, 2),
    output:
        r1=f"{RESULTS}/qc/{{sample}}_R1.clean.fastq.gz",
        r2=f"{RESULTS}/qc/{{sample}}_R2.clean.fastq.gz",
        html=f"{RESULTS}/qc/{{sample}}_fastp.html",
        json=f"{RESULTS}/qc/{{sample}}_fastp.json",
    log:
        f"{LOGS}/fastp/{{sample}}.log",
    threads: 8
    conda:
        "../envs/qc.yaml"
    params:
        quality=config.get("fastp_quality", 20),
        length=config.get("fastp_min_length", 50),
    shell:
        """
        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --html {output.html} \
            --json {output.json} \
            --detect_adapter_for_pe \
            --qualified_quality_phred {params.quality} \
            --length_required {params.length} \
            --thread {threads} \
            2> {log}
        """
