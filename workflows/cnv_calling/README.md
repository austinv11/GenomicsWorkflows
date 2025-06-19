# scRNA-seq CNV Calling Workflow
This workflow takes an annotated h5ad file (i.e. from the celltype_assignment workflow) and performs CNV/tumor cell calling using:
- inferCNV

```bash
snakemake \
  --snakefile GenomicsWorkflows/workflows/cnv_calling/Snakefile \
  --cores 32 --max-threads 32 \
  --configfile config.yaml \
  --software-deployment-method apptainer \   # or --software-deployment-method conda
  insert_{method}_assignments
```

View an example config file [here](config.yaml).
