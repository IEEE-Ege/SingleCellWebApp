# Import libraries
import scanpy as sc
import matplotlib.pyplot as plt
import tempfile

# PCA plot function
def pca_plot(adata, log=True, color='CST3'):
    """_summary_

    Args:
        adata (_type_): _description_æ
        log (bool, optional): _description_. Defaults to True.
        color (str, optional): _description_. Defaults to 'CST3'.

    Raises:
        ValueError: _description_
    """
    if adata is None:
        raise ValueError("Input AnnData object (adata) is None.")
    sc.pl.pca_variance_ratio(adata, log=log)
    sc.pl.pca(adata, color=color)

# Highly variable genes plot
def highvarGen_pl(adata, min_mean=0.0125, max_mean=3, min_disp=0.5):
    """_summary_

    Args:
        adata (_type_): _description_
        min_mean (float, optional): _description_. Defaults to 0.0125.
        max_mean (int, optional): _description_. Defaults to 3.
        min_disp (float, optional): _description_. Defaults to 0.5.

    Raises:
        ValueError: _description_
        ValueError: _description_
    """
    if adata is None:
        raise ValueError("Input AnnData object (adata) is None.")
    if "Cell type" not in adata.obs:
        raise ValueError("'Cell type' not found in adata.obs. Please check your annotations.")
    
    sc.pl.umap(adata, color="Cell type", show=False)
    sc.pl.highly_variable_genes(adata)

# Violin plot function
def violin_pl(adata):
    """_summary_

    Args:
        adata (_type_): _description_

    Raises:
        ValueError: _description_
        ValueError: _description_
    """
    if adata is None:
        raise ValueError("Input AnnData object (adata) is None.")
    if "Cell type" not in adata.obs:
        raise ValueError("'Cell type' not found in adata.obs. Please check your annotations.")
    
    sc.pl.umap(adata, color="Cell type", show=False)
    sc.pl.violin(
        adata,
        ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
        jitter=0.4,
        multi_panel=True
    )

# Scatter plot function
def scatter_plotting(adata, x='total_counts', y='pct_counts_mt'):
    """_summary_

    Args:
        adata (_type_): _description_
        x (str, optional): _description_. Defaults to 'total_counts'.
        y (str, optional): _description_. Defaults to 'pct_counts_mt'.

    Raises:
        ValueError: _description_
        ValueError: _description_
    """
    if adata is None:
        raise ValueError("Input AnnData object (adata) is None.")
    if "Cell type" not in adata.obs:
        raise ValueError("'Cell type' not found in adata.obs. Please check your annotations.")
    
    sc.pl.umap(adata, color="Cell type", show=False)
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


