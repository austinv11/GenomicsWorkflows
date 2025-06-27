# Based on https://www.10xgenomics.com/analysis-guides/creating-single-cell-references-for-xenium-custom-panel-design-from-seurat-or-anndata
import pandas as pd  # pandas
import scanpy as sc  # scanpy
import pybiomart
import scipy # scipy
import scipy.sparse as sparse
import zipfile
import tempfile
import os
import gzip

# Define function
def h5ad_to_10x(adata,
                gene_id_key="gene_id",
                gene_name_key="gene_name",
                cell_type_key="cell_type",
                output_path="matrix.zip",
                organism="hsapiens",
                barcode_key=None,
                subsample_rate=None):
    """
    Convert an AnnData h5ad file to 10x format,
    then produce a zipfile of the final matrix.

    Arguments:
        adata: A scanpy AnnData object.
        gene_id_key: The key that the gene IDs are under.
        gene_name_key: The key that the gene names are under.
        cell_type_key: The key that the cell types are under.
        organism: The organism to use for gene ID conversion.
        barcode_key: Optional key that the barcodes are under. If not set,
            will assume they are the index of `obs`.
        output_path: Path to write the zipfile to.
        subsample_rate: Optional argument for subsampling.
            If provided, should be between 0 and 1.
    """
    if subsample_rate:
        sc.pp.sample(adata, subsample_rate)

    var_df = adata.var.reset_index()
    if gene_id_key not in var_df.columns and gene_name_key not in var_df.columns:
        raise ValueError(f"Neither {gene_id_key} nor {gene_name_key} found in var DataFrame.")
    elif gene_id_key not in var_df.columns and gene_name_key in var_df.columns:
        # Query to fill in gene_id_key with gene_name_key
        print(f"Warning: {gene_id_key} not found in var DataFrame. Using {gene_name_key} to join with gene name.")
        annot = sc.queries.biomart_annotations(
            organism, ["ensembl_gene_id", "external_gene_name"]
        )
        # Join the var_df with the annotations, ensure that we do not get duplicates and order is preserved
        var_df = var_df.merge(annot.drop_duplicates(subset=["external_gene_name"]), left_on=gene_name_key, right_on="external_gene_name", how="left")
        # Rename the columns to match the expected keys
        var_df.rename(columns={"ensembl_gene_id": gene_id_key}, inplace=True)
    elif gene_id_key in var_df.columns and gene_name_key not in var_df.columns:
        # Query to fill in gene_name_key with gene_id_key
        print(f"Warning: {gene_name_key} not found in var DataFrame. Using {gene_id_key} to join with gene name.")
        annot = sc.queries.biomart_annotations(
            organism, ["ensembl_gene_id", "external_gene_name"]
        )
        # Join the var_df with the annotations, ensure that we do not get duplicates and order is preserved
        var_df = var_df.merge(annot, left_on=gene_id_key, right_on="ensembl_gene_id", how="left")
        # Rename the columns to match the expected keys
        var_df.rename(columns={"external_gene_name": gene_name_key}, inplace=True)

    genes = var_df[[gene_id_key, gene_name_key]]
    genes["feature_type"] = ["Gene Expression"] * len(genes)

    if barcode_key:
        barcodes = adata.obs[[barcode_key]]
    else:
        barcodes = pd.DataFrame(adata.obs.index)

    celltypes = adata.obs[[cell_type_key]].reset_index()
    celltypes.columns = ["barcode", "annotation"]

    # Ensure that the X matrix is an integer type sparse matrix
    adata.X = adata.X.astype(int)
    with tempfile.TemporaryDirectory() as tmp_dir:
        with gzip.open(os.path.join(tmp_dir, "matrix.mtx.gz"), "w") as handle:
            # Write in compressed sparse column matrix format
            scipy.io.mmwrite(handle, sparse.csc_matrix(adata.X.T), field="integer", precision=0)

        genes.to_csv(os.path.join(tmp_dir, "features.tsv.gz"), sep="\t", index=False, header=False, compression="gzip")
        barcodes.to_csv(os.path.join(tmp_dir, "barcodes.tsv.gz"), sep="\t", index=False, header=False, compression="gzip")
        celltypes.to_csv(os.path.join(tmp_dir, "celltypes.csv"), index=False)

        with zipfile.ZipFile(output_path, "w") as zip_handle:
            for file in ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz", "celltypes.csv"]:
                zip_handle.write(os.path.join(tmp_dir, file), arcname=file)

    # Check the final output size, the panel designer requires the file to be 2G or smaller
    output_size = os.path.getsize(output_path)
    gb_size = output_size / (1024 * 1024 * 1024)  # Convert to GB
    if gb_size > 2:
        # Delete the output file if it exceeds the size limit
        os.remove(output_path)
        raise ValueError(f"Output file {output_path} is larger than 2 GB (current size={gb_size} GB). "
                         "Please reduce the number of genes or cells, or use a greater subsampling rate.")
    else:
        print(f"Output file {output_path} created successfully with size {gb_size:.2f} GB.")
        print("You can now upload this file to the 10x Xenium Panel Designer.")


if __name__ == "__main__":
    import argparse
    # Note if var_names has no column name, you can use "index" as reset_index() will automatically create that column in our code
    parser = argparse.ArgumentParser(description="Convert h5ad to 10x Xenium Panel Designer format.")
    parser.add_argument("h5ad_file", type=str, help="Path to the input h5ad file.")
    parser.add_argument("--gene_id_key", type=str, default="gene_id", help="Key for gene IDs.")
    parser.add_argument("--gene_name_key", type=str, default="gene_name", help="Key for gene names.")
    parser.add_argument("--cell_type_key", type=str, default="cell_type", help="Key for cell types.")
    parser.add_argument("--barcode_key", type=str, help="Key for barcodes (optional).")
    parser.add_argument("--subsample_rate", type=float, default=None, help="Subsample rate (between 0 and 1).")
    parser.add_argument("--output_path", type=str, default="matrix.zip", help="Output path for the zipfile.")
    parser.add_argument("--organism", type=str, default="hsapiens", help="Organism for gene ID conversion.")

    args = parser.parse_args()

    adata = sc.read_h5ad(args.h5ad_file)
    h5ad_to_10x(adata, args.gene_id_key, args.gene_name_key, args.cell_type_key,
                args.output_path, args.organism, args.barcode_key, args.subsample_rate)
