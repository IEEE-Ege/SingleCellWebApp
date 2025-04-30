import scanpy as sc
import pandas as pd
from get_marker_genes import get_marker_genes  # defined function in main code

# Demo dataset
adata = sc.datasets.pbmc3k()

#preprocessing and clustering (same with the main code)
sc.pp.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata)
sc.tl.leiden(adata)  # Leiden clustering yap

# ranking genes between groups
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')

# recalling the function
scores_and_pvals = get_scores_and_pvals(adata)
print(scores_and_pvals)