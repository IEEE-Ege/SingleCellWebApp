import scanpy as sc
adata = sc.datasets.pbmc3k()
sc.pp.neighbors(adata)
sc.tl.leiden(adata, resolution=1.0)
sc.settings.verbosity = 2  # reduce the verbosity

def get_rank_genes_groups_violin(adata, groups='0', n_genes=8):
   """
   -Plots a violin plot for the top ranked genes in specified groups (clusters).
   
    Parameters:
    - adata: AnnData object with results from rank_genes_groups
    - groups: List or string specifying the groups (clusters) to plot (default is all groups)
    - n_genes: Number of top genes to display in the plot (default is 8)

   """
   sc.pl.rank_genes_groups_violin(adata, groups=groups, n_genes=n_genes)
   return adata
