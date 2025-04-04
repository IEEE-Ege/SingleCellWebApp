import scanpy as sc

def load_file(file_path):
    adata = sc.read_10x_mtx(file_path, var_names='gene_symbols', cache=True)
    return adata