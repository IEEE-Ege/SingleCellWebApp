import scanpy as sc
import anndata
import os

def read_10x_matrix_folder(matrix_dir: str, sample_name: str = None) -> anndata.AnnData:
    """
    Loads a 10x Genomics dataset (already extracted folder) using Scanpy.

    Parameters:
        matrix_dir (str): Path to the directory containing matrix.mtx, barcodes.tsv, and features.tsv or genes.tsv
        sample_name (str, optional): If given, adds this as a 'source' column in adata.obs

    Returns:
        anndata.AnnData: The loaded AnnData object
    """
    if not os.path.isdir(matrix_dir):
        raise FileNotFoundError(f"Directory not found: {matrix_dir}")

    adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", cache=True)

    if sample_name:
        adata.obs["source"] = sample_name

    return adata
# for running the function
matrix_path = "/content/temp_extracted/filtered_gene_bc_matrices/hg19"

adata = read_10x_matrix_folder(matrix_path, sample_name="pbmc3k")
print(adata)
print("🧪 OBS columns:", adata.obs.columns.tolist())

