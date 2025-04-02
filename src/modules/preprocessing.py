import scanpy as sc

def filter_data(adata, min_genes=200, min_cells=3):
    """
    Filters the dataset:
    - Ensures that cells have at least 'min_genes' genes.
    - Ensures that genes are present in at least 'min_cells' cells.
    """
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata

def calculate_qc_metrics(adata):
    """
    Computes quality control (QC) metrics and adds them to the dataset.
    """
    adata.var["mt"] = adata.var_names.str.startswith("MT-")  # Identify mitochondrial genes
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    return adata

def select_highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5):
    """
    Identifies highly variable genes in the dataset.
    """
    sc.pp.highly_variable_genes(adata, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp)
    return adata
