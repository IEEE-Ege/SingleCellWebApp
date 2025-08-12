import tarfile
import os
import pandas as pd
import scipy.io
import scanpy as sc
import anndata
from scipy import sparse
def extract_10x_tar_to_anndata(tar_path: str, extract_dir: str = "./temp") -> anndata.AnnData:
    os.makedirs(extract_dir, exist_ok=True)

    with tarfile.open(tar_path, "r") as tar:
        tar.extractall(path=extract_dir)

    matrix_file = None
    barcodes_file = None
    genes_file = None

    for root, _, files in os.walk(extract_dir):
        for file in files:
            path = os.path.join(root, file)
            if file.endswith("matrix.mtx"):
                matrix_file = path
            elif file.endswith("barcodes.tsv"): # More specific check for barcodes file
                barcodes_file = path
            elif file.endswith("genes.tsv") or file.endswith("features.tsv"): # Check for both genes.tsv and features.tsv
                genes_file = path

    if not all([matrix_file, barcodes_file, genes_file]): # checking for all the necessary files if not foun an error will be raised.

        # Construct the missing list with descriptive strings instead of None
        missing = []
        if matrix_file is None:
            missing.append(".mtx matrix file")
        if barcodes_file is None:
            missing.append(".tsv file")
        if genes_file is None:
            missing.append(".tsv file")

        raise FileNotFoundError(f"Required file(s) not found in tar file: {', '.join(missing)}")

    matrix = scipy.io.mmread(matrix_file)
    barcodes = pd.read_csv(barcodes_file, header=None)[0].tolist()

    genes_df = pd.read_csv(genes_file, header=None)

    if genes_df.shape[1] >= 2: # if there are 2 or more columns in genes file, this code takes the first column as gene id and the second as gene names
        genes = genes_df[1].tolist()
        gene_ids = genes_df[0].tolist()
    else:
        genes = genes_df[0].tolist()
        gene_ids = genes_df[0].tolist()

    matrix_transposed = matrix.transpose().tocsc()  # transposing the matrice because anndata format is obsxvar

    adata = anndata.AnnData(X=matrix_transposed)
    adata.obs_names = barcodes
    adata.var_names = genes
    if gene_ids and gene_ids != genes:

         adata.var['gene_id'] = gene_ids
# this code checks if there's a list called gene_ids and is it's different from the existing genes list if both are true it adds these gene_ids as a new column named 'gene_id' to the AnnData object's gene-specific information (adata.var).

    return adata
#testing
if __name__ == "__main__":
    test_file = "/pbmc3k_filtered_gene_bc_matrices.tar.gz"
    try:
        adata = extract_10x_tar_to_anndata(test_file)
        print("test is succesfull, AnnData :")
        print(adata)
    except FileNotFoundError:
        print(f" Error: {test_file} not found! Please try another file path.")



