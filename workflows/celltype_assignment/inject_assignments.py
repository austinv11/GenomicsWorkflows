import scanpy as sc
import pandas as pd



adata = sc.read_h5ad(snakemake.input.input_file)
assignments = pd.read_csv(snakemake.input.assignments, index_col=0)
adata.obs = adata.obs.join(assignments, how='left')
adata.write(snakemake.output.output_file, compression='gzip')
