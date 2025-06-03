# scRNA-seq Cell Type Assignment Workflow
This workflow takes an h5ad file (i.e. from the scrnaseq workflow) and performs cell type assignment using one of the following methods:
- CellAssign
- SingleR
- CellTypist
- ProjectTILs
- Azimuth
- Garnett

```bash
snakemake \
  --snakefile GenomicsWorkflows/workflows/celltype_assignment/Snakefile \
  --cores 32 --max-threads 32 \
  --configfile config.yaml \
  --software-deployment-method apptainer   # or --software-deployment-method conda
```

View an example config file [here](config.yaml).
