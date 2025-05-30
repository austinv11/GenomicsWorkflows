import scanpy as sc

# Load all h5ad files into a list of AnnData objects
adatas = []
for h5ad_file in snakemake.input.h5ad_files:
    adata = sc.read_h5ad(h5ad_file)
    # adata.obs['sample'] = "_".join(h5ad_file.split("/")[-1].split("_")[:-1])  # Extract sample name from filename
    adatas.append(adata)

# Concatenate all AnnData objects
adatas = sc.concat(
    adatas,
    label='sample',  # Use 'sample' as the label for concatenation
    join='outer',  # Use outer join to keep all genes
    index_unique=None  # Avoid unique index to keep all cells
)
adatas.obs_names_make_unique()

# Save the concatenated AnnData object to an h5ad file
adatas.write_h5ad(snakemake.output.combined_h5ad)
