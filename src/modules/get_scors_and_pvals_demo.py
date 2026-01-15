from get_scors_and_pvals import get_scores_and_pvals
import scanpy as sc
import pandas as pd

# Demo dataset
adata = sc.datasets.pbmc3k() 

# ranking genes between groups
sc.pp.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata)
sc.tl.leiden(adata)
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')

# recalling the function
scores_and_pvals = get_scores_and_pvals(adata)
print(scores_and_pvals)
 