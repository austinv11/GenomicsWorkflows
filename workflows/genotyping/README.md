# Genotyping Workflow
A basic end-to-end pipeline for mapping WGS reads with to hg38 and calling mutations.

The pipeline will perform the following steps:
- Run fastp to trim sequencing and calculate QC metrics.
- Download and index the hg38 reference genome.
- Alignment with either BWA-MEM, Bowtie2, minimap2.
- Variant calling with either GATK Mutect2 or GLIMPSE2 (for low-coverage)
- Structural variant calling with Delly2.
- Run MultiQC to generate a report of the QC metrics from fastp and alignments.

Final outputs include: A MultiQC report, ...

```bash
snakemake \
  --snakefile GenomicsWorkflows/workflows/genotyping/Snakefile \
  --cores 32 --max-threads 32 \
  --configfile config.yaml \
  --software-deployment-method apptainer   # or --software-deployment-method conda
```

View an example config file [here](config.yaml).
