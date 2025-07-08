# DE-analysis Workflow
This workflow performs differential expression analysis using the following methods:
- `MAST`
- `DESeq2` (psuedo-bulk)
- TODO

```bash
snakemake \
  --snakefile GenomicsWorkflows/workflows/de_analysis/Snakefile \
  --cores 32 --max-threads 32 \
  --configfile config.yaml \
  --software-deployment-method apptainer \  # or --software-deployment-method conda
  compute_{method}_DEGs
```

View an example config file [here](config.yaml).
