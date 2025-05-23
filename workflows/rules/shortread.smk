"""
Snakemake workflow for processing short-read sequencing data.

This workflow handles paired-end short reads by:
1. Ensuring all input files are gzipped
2. Merging multiple lanes if present
3. Quality control and trimming using fastp

Required config keys:
    fastq_dir: Directory containing input fastq files
    results_dir: Directory for output files

Input files should follow the naming convention:
    {sample}_R{1,2}_*.f*q*

Author: Austin Varela
"""

from pathlib import Path
import glob

# Config for input/output directories
config.setdefault("registry", "biohpc_aav4003")
config["fastq_dir"] = "raw_data"
config["results_dir"] = "results"
config["low_complexity_filter"] = True

# Get sample names from fastq files
def get_samples():
    fastqs = glob.glob(f"{config['fastq_dir']}/*_R1_*.f*q*")
    return list(set([Path(f).name.split("_R1_")[0] for f in fastqs]))


# Input functions
def get_fastq_r1(wildcards):
    return sorted(glob.glob(f"{config['fastq_dir']}/{wildcards.sample}_R1_*.f*q*"))


def get_fastq_r2(wildcards):
    return sorted(glob.glob(f"{config['fastq_dir']}/{wildcards.sample}_R2_*.f*q*"))


# Pipeline rules
rule all:
    input:
        expand(
            "{results}/fastp/{sample}_R{read}.fastq.gz",
            results=config["results_dir"],
            sample=get_samples(),
            read=[1, 2]
        )

rule ensure_gzipped:
    input:
        fastq="{dir}/{file}.fastq"
    output:
        gz="{dir}/{file}.fastq.gz"
    shell:
        "gzip -c {input.fastq} > {output.gz}"

rule merge_lanes:
    input:
        r1=get_fastq_r1,
        r2=get_fastq_r2
    output:
        r1=temp("{results}/merged/{sample}_R1.fastq.gz"),
        r2=temp("{results}/merged/{sample}_R2.fastq.gz")
    shell:
        """
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
        """

rule fastp:
    container: f"docker://${config['registry']}/multiqc:latest"
    input:
        r1="{results}/merged/{sample}_R1.fastq.gz",
        r2="{results}/merged/{sample}_R2.fastq.gz"
    output:
        r1="{results}/fastp/{sample}_R1.fastq.gz",
        r2="{results}/fastp/{sample}_R2.fastq.gz",
        json="{results}/fastp/{sample}.json",
        html="{results}/fastp/{sample}.html"
    log:
        "{results}/fastp/{sample}.log"
    threads: 4
    params:
        low_complexity = "--low_complexity_filter" if config["low_complexity_filter"] else ""
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              --json {output.json} --html {output.html} \
              --detect_adapter_for_pe \
              --thread {threads} {params.low_complexity} > {log} 2>&1
        """

wildcard_constraints:
    sample="[^/]+",
    read="[12]"
