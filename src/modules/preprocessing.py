import scanpy as sc
import pandas as pd
from anndata import AnnData

def filter_data(
        adata: AnnData, 
        min_genes: int = 200, 
        min_cells: int = 3
) -> AnnData:
    """
    Filters the dataset.

    Args:
        adata (AnnData): The single-cell dataset.
        min_genes (int): Minimum number of genes expressed per cell. Default is 200.
        min_cells (int): Minimum number of cells expressing a gene. Default is 3.

    Returns:
        AnnData: Filtered dataset.
    """
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata

def calculate_qc_metrics(adata: AnnData) -> AnnData:
    """
    Computes quality control (QC) metrics and adds them to the dataset.
    Checks for mitochondrial genes with both uppercase and lowercase 'mt-' prefix.

    Args:
        adata (AnnData): The filtered dataset.

    Returns:
        AnnData: Dataset with QC metrics added.
    """
    if adata.var_names.str.startswith("MT-").any():
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
    elif adata.var_names.str.startswith("mt-").any():
        adata.var["mt"] = adata.var_names.str.startswith("mt-")
    else:
        # In case neither is found, check in a case-insensitive way just to be sure
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=["mt"], 
        percent_top=None, 
        log1p=False, 
        inplace=True
    )
    return adata

def select_highly_variable_genes(
        adata: AnnData, 
        min_mean: float = 0.0125, 
        max_mean: int = 3, 
        min_disp: float = 0.5
    ) -> AnnData:
    """
    Identifies highly variable genes in the dataset.

    Args:
        adata (AnnData): The dataset with QC metrics.
        min_mean (float): Minimum mean expression threshold. Default is 0.0125.
        max_mean (float): Maximum mean expression threshold. Default is 3.
        min_disp (float): Minimum dispersion threshold. Default is 0.5.

    Returns:
        AnnData: Dataset with highly variable gene information.
    """
    sc.pp.highly_variable_genes(
        adata, 
        min_mean=min_mean, 
        max_mean=max_mean, 
        min_disp=min_disp
    )
    return adata
