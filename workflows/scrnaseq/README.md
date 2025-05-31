# scRNA-seq Workflow
A basic end-to-end pipeline for mapping paired-end short reads with 10X Cellranger with support for 3' and fixed RNA
sequencing.

The pipeline will perform the following steps:
- Run fastp to get QC metrics.
- Run Cellranger to map reads to single-cells.
- Run MultiQC to generate a report of the QC metrics from fastp and Cellranger.
- Convert the Cellranger output to an AnnData h5ad file.
- Combine multiple samples into a single h5ad file, with a `sample` obs field.
- Run Scanpy to perform basic filtering and preprocessing of the h5ad file, including HVG calling, PCA, and clustering.

Final outputs include: A MultiQC report, a combined h5ad file, and Scanpy preprocessing reports.

```bash
snakemake \
  --snakefile GenomicsWorkflows/workflows/scrnaseq/Snakefile \
  --cores 32 --max-threads 32 \
  --configfile config.yaml \
  --software-deployment-method apptainer   # or --software-deployment-method conda
```

View an example config file [here](config.yaml).
