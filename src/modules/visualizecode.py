# Import libraries
import scanpy as sc
import matplotlib.pyplot as plt
import tempfile

# PCA plot function
def pca_plot(adata, log=True, color='CST3'):
    if adata is None:
        raise ValueError("Input AnnData object (adata) is None.")
    sc.pl.pca_variance_ratio(adata, log=log)
    sc.pl.pca(adata, color=color)

# Highly variable genes plot
def highvarGen_pl(adata, min_mean=0.0125, max_mean=3, min_disp=0.5):
    if adata is None:
        raise ValueError("Input AnnData object (adata) is None.")
    if "Cell type" not in adata.obs:
        raise ValueError("'Cell type' not found in adata.obs. Please check your annotations.")
    
sc.pl.umap(adata, color="Cell type", show=False)
sc.pl.highly_variable_genes(adata)

# Violin plot function
def violin_pl(adata):
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
    if adata is None:
        raise ValueError("Input AnnData object (adata) is None.")
    if "Cell type" not in adata.obs:
        raise ValueError("'Cell type' not found in adata.obs. Please check your annotations.")
    
sc.pl.umap(adata, color="Cell type", show=False)
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


def demo(adata,):
    return pca_plot
    return highvarGen_pl
    return violin_pl
    return scatter_plotting
    
