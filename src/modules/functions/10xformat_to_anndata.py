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

if __name__ == "__main__":
    test_matrix_dir = "/Users/eceyigit/Desktop/filtered_gene_bc_matrices 2/hg19"  # Replace with actual path
    test_sample_name = "test_sample"

    print(f"Attempting to load data from: {test_matrix_dir}")

    try:
        adata = read_10x_matrix_folder(test_matrix_dir, test_sample_name)
        print("\nSuccessfully loaded data!")
        print(f"AnnData object shape: {adata.shape}")
        print(f"Sample name: {test_sample_name}")
        print("\nFirst few observations:")
        print(adata.obs.head())
    except Exception as e:
        print(f"\nError occurred: {str(e)}")
        print("Please check:")
        print(f"- Does the directory exist? {os.path.exists(test_matrix_dir)}")
        print(f"- Does it contain matrix.mtx, barcodes.tsv, and features.tsv files?")