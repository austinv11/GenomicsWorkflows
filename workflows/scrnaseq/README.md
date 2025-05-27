# scRNA-seq Workflow
A basic end-to-end pipeline for mapping paired-end short reads with 10X Cellranger with support for 3' and fixed RNA
sequencing.

```bash
snakemake \
  --snakefile GenomicsWorkflows/workflows/scrnaseq/Snakefile \
  --cores 16 --max-threads 16 \
  --configfile config.yaml \
  --software-deployment-method apptainer \
  --software-deployment-method conda
```
