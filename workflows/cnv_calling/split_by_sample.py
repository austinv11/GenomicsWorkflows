from pathlib import Path

import scanpy as sc

adata = sc.read_h5ad(snakemake.input.input_file)

# Split the AnnData object by sample
for sample_name in adata.obs[snakemake.params.sample_column].unique():
    sample_adata = adata[adata.obs[snakemake.params.sample_column] == sample_name].copy()
    output_file = Path(snakemake.output.h5ad_files) / f"{sample_name}.h5ad"
    sample_adata.write_h5ad(output_file)
