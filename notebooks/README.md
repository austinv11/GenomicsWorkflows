# Jupyter Notebooks
This directory contains Jupyter notebooks that are used by the Snakemake workflows in this repository to
generate reports and analysis.

## Notebooks:
- `scrnaseq_preprocessing_pt1.ipynb`: Performs basic QC of a h5ad file, filtering low quality cells.
- `scrnaseq_preprocessing_pt2.ipynb`: Performs basic preprocessing of a h5ad file, including HVG calling, PCA, and clustering.
- `doublet_detection.ipynb`: Performs doublet detection on a h5ad file using scDblFinder.
- `cellassign_celltype_assignment.ipynb`: Performs cell type assignment using CellAssign on a h5ad file.
- `garnett_train_classifier.ipynb`: Trains a Garnett classifier on a h5ad file.
- `garnett_celltype_assignment.ipynb`: Performs cell type assignment using Garnett on a h5ad file.
- `singler_celltype_assignment.ipynb`: Performs cell type assignment using SingleR on a h5ad file with a specified reference.
- `infercnv_cnv_calling.ipynb`: Performs CNV calling using infercnv on a h5ad file.
