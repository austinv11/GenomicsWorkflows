import scanpy as sc

adata = sc.read_h5ad(snakemake.input.input_file)

# Split the AnnData object by sample
for sample_name in adata.obs[snakemake.params.sample_name].unique():
    sample_adata = adata[adata.obs[snakemake.params.sample_name] == sample_name].copy()
    output_file = snakemake.output.h5ad_files / f"{sample_name}.h5ad"
    sample_adata.write_h5ad(output_file)
