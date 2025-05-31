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
import re

# Config for input/output directories
config.setdefault("registry", "biohpc_aav4003")
config.setdefault("fastq_dir", "raw_data")
config.setdefault("results_dir", "results")
config.setdefault("low_complexity_filter", True)

# Pipeline rules
# rule all:
#     input:
#         expand(
#             "{results}/fastp/{sample}_R{read}.fastq.gz",
#             results=config["results_dir"],
#             sample=get_samples(),
#             read=[1, 2]
#         )

SAMPLE_RE = r"[A-Za-z0-9._-]+"
READ_RE   = r"[12]"
LANE_RE   = r"\d{3}"

wildcard_constraints:
    sample = SAMPLE_RE,
    read   = READ_RE,
    lane   = LANE_RE

def get_samples():
    return sorted({
        re.sub(r"_L\d{3}(_R[12].*)?$", "", Path(f).stem)   # strip lane & read
                   for f in glob.glob(f"{config['fastq_dir']}/*_R1_*.f*q*")})


def get_fastp_dir():
    return os.path.abspath(f"{config['results_dir']}/fastp")


rule ensure_gzipped:
    """
    Copy *.fastq.gz as-is, or gzip *.fastq, so every lane/read file ends up
    in {results_dir}/gzipped/…
    """
    input:
        fq= lambda wc: glob.glob(f"{config['fastq_dir']}/{wc.sample}_L{wc.lane}_R{wc.read}_001.*"),
    output:
        gz=temp(f"{config['results_dir']}/gzipped/{{sample}}_L{{lane}}_R{{read}}_001.fastq.gz")
    threads: 2
    shell: r"""
        mkdir -p $(dirname {output.gz})
         if [[ ! -s {input.fq} ]]; then
            echo "Missing or empty input: {input.fq}" >&2
            exit 1
        fi
        if [[ "{input.fq}" == *.gz ]]; then
            cp {input.fq} {output.gz}
        else
            gzip -c {input.fq} > {output.gz}
        fi
    """


rule merge_lanes:
    input:
        r1=lambda wc: expand(
            f"{config['results_dir']}/gzipped/{{sample}}_L{{lane}}_R1_001.fastq.gz",
            sample=[wc.sample],
            lane=[Path(f).stem.split("_L")[1][:3] for f in glob.glob(f"{config['fastq_dir']}/{wc.sample}_L*_R1_*.f*q*")]
        ),
        r2=lambda wc: expand(
            f"{config['results_dir']}/gzipped/{{sample}}_L{{lane}}_R2_001.fastq.gz",
            sample=[wc.sample],
            lane=[Path(f).stem.split("_L")[1][:3] for f in glob.glob(f"{config['fastq_dir']}/{wc.sample}_L*_R2_*.f*q*")]
        )
    output:
        r1=temp(f"{config['results_dir']}/merged/{{sample}}_L001_R1_001.fastq.gz"),
        r2=temp(f"{config['results_dir']}/merged/{{sample}}_L001_R2_001.fastq.gz")
    shell:
        """
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
        """


rule fastp:
    container: f"{config['workflow_dir']}/docker/multiqc/multiqc_image.sif"
    shadow: "copy-minimal"  # Required when using .sif
    input:
        r1=f"{config['results_dir']}/merged/{{sample}}_L001_R1_001.fastq.gz",
        r2=f"{config['results_dir']}/merged/{{sample}}_L001_R2_001.fastq.gz",
    output:
        r1=f"{get_fastp_dir()}/{{sample}}_L001_R1_001.fastq.gz",
        r2=f"{get_fastp_dir()}/{{sample}}_L001_R2_001.fastq.gz",
        json=f"{get_fastp_dir()}/{{sample}}.json",
        html=f"{get_fastp_dir()}/{{sample}}.html",
    log:
        f"{get_fastp_dir()}/{{sample}}.log"
    threads: 8
    params:
        low_complexity = "--low_complexity_filter" if config["low_complexity_filter"] else ""
    shell:
        """
        touch {log}
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              --json {output.json} --html {output.html} \
              -R "{wildcards.sample} Report" \
              --detect_adapter_for_pe \
              --thread {threads} {params.low_complexity} > {log} 2>&1
        """


rule fastp_report_only:
    container: f"{config['workflow_dir']}/docker/multiqc/multiqc_image.sif"
    shadow: "copy-minimal"  # Required when using .sif
    input:
        r1=f"{config['results_dir']}/merged/{{sample}}_L001_R1_001.fastq.gz",
        r2=f"{config['results_dir']}/merged/{{sample}}_L001_R2_001.fastq.gz",
    output:
        json=f"{get_fastp_dir()}/{{sample}}_report.json",
        html=f"{get_fastp_dir()}/{{sample}}_report.html",
    log:
        f"{get_fastp_dir()}/{{sample}}.log"
    threads: 8
    shell:
        """
        touch {log}
        fastp -i {input.r1} -I {input.r2} \
              --json {output.json} --html {output.html} \
              -R "{wildcards.sample} Report" \
              --disable_trim_poly_g \
              --disable_adapter_trimming \
              --disable_length_filtering \
              --disable_quality_filtering \
              --detect_adapter_for_pe \
              --overrepresentation_analysis \
              --stdout /dev/null \
              --thread {threads} > {log} 2>&1
        """
