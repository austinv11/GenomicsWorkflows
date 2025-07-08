# Snakemake script to load an anndata and save all unique groups to an output file
import scanpy as sc

adata = sc.read_h5ad(snakemake.input.input_file)
# Extract unique groups from the specified observation
unique_groups = adata.obs[snakemake.config['comparison_obs']].unique()
# Save the unique groups to the output file
with open(snakemake.output.groups, 'w') as f:
    for group in unique_groups:
        f.write(f"{group}\n")

exit(0)  # Exit successfully