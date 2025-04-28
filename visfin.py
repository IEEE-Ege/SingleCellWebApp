# Import libraries
import scanpy as sc
import matplotlib.pyplot as plt
import tempfile

# violin_plotting(adata, n_components)
# pca_plotting(adata)
# scatter_plotting(adata, x, y)
# plotting_highly_variable_genes(adata, min_mean, max_mean, min_disp)


def pca_plot(adata,log=True,color= 'CST3'):
    sc.pp.pca(adata)
    sc.pl.pca_variance_ratio(adata, log=True)
    sc.pl.pca(adata, color) 


def highvarGen_pl(adata, min_mean=0.0125, max_mean= 3, min_disp= 0.5):  
    if "Cell type" in adata.obs:
        sc.pl.umap(adata, color="Cell type", show=False)
    else:
        print("Error: 'Cell type' variable couldn't find in adata.obs.")
    sc.pl.highly_variable_genes(adata)
    
def violin_plotting(adata):    
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'],jitter=0.4, multi_panel=True)

def scatter_plotting(adata,x,y):
    if "Cell type" in adata.obs:
        sc.pl.umap(adata, color="Cell type", show=False)
    else:
        print("Error: 'Cell type' variable couldn't find in adata.obs.")
        
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
