import scanpy as sc
import pandas as pd
import pyranges as pr

adata = sc.read_h5ad(snakemake.input.input_file)

state_map = {1: 'loss',
             2: 'loss',
             3: 'neutral',
             4: 'gain',
             5: 'gain',
             6: 'gain'}

cytoband = pr.read_bed(snakemake.input.cytoband)
cytoband = cytoband[~cytoband.Name.str.contains("_")]
cb_arms = (cytoband.assign(arm=cytoband.Name.str[0])  # 'p' or 'q'
              .df.groupby(["Chromosome", "arm"])[["Start", "End"]]
              .agg({"Start":"min", "End":"max"})
              .reset_index())
cb_arms = pr.PyRanges(cb_arms)

# Get all arm names
arm_names = cb_arms.df.apply(lambda x: f"{x['Chromosome']}{x['arm']}", axis=1)

for cell_grouping_file, gene_cnv_file in zip(
    snakemake.input.cell_grouping_files,
    snakemake.input.gene_cnv_files,
):
    cell_grouping = pd.read_csv(cell_grouping_file, sep="\t")
    gene_cnv = pd.read_csv(gene_cnv_file, sep="\t")

    gene_cnv['state'] = gene_cnv['state'].map(state_map)

    cnv_pr = pr.PyRanges(gene_cnv.rename(columns={"chr": "Chromosome",
                                             "start": "Start", "end": "End"}))
    annot = cnv_pr.join(cb_arms).df
    annot["arm_id"] = annot["Chromosome"] + annot["arm"]

    summary = (annot.groupby(["cell_group_name", "arm_id"])
                .state
                .agg(lambda x: x.mode()[0] if not x.mode().empty else "neutral")
                .unstack()).reset_index()
    # Expand summary_df to map to all cells according to the cell_grouping mapping
    summary_df = pd.merge(
        cell_grouping[["cell", "cell_group_name"]],
        summary,
        on="cell_group_name",
        how="left"
    ).set_index("cell").fillna("neutral")

    if 'infercnv' not in adata.obsm:
        adata.obsm['infercnv'] = pd.DataFrame(index=adata.obs_names,
                                               columns=arm_names,
                                               data="neutral")
        adata.obs['fraction_genes_altered'] = 0.0

    adata.obsm['infercnv'].loc[summary_df.index, summary_df.columns] = summary_df.values

    frac_genes_altered = (
        pd.merge(cell_grouping[["cell", "cell_group_name"]], gene_cnv, on="cell_group_name", how="left")
        .groupby("cell").state
        .apply(lambda x: (x != "neutral").mean())
    ).fillna(0).reset_index().rename(columns={"state": "fraction_genes_altered"})

    adata.obs.loc[frac_genes_altered['cell'], 'fraction_genes_altered'] = frac_genes_altered['fraction_genes_altered']


# Save the updated adata object
adata.write_h5ad(snakemake.output.output_file, compression='gzip')