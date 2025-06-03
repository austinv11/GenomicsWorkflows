import scanpy as sc
from glob import glob
import os

tenx_outdir = snakemake.input.tenx_output
outfile = snakemake.output.output_file
sample_name = snakemake.input.sample_name

# Find the 10X output directory containing the `filtered_feature_bc_matrix` directory
# Recursively search for the directory
tenx_dirs = glob(os.path.join(tenx_outdir, "**/filtered_feature_bc_matrix"), recursive=True)

adata = sc.read_10x_mtx(
    tenx_dirs[0],
    gex_only=False
)

adata.obs['sample'] = sample_name

# Save the AnnData object to an h5ad file
adata.write_h5ad(outfile)
