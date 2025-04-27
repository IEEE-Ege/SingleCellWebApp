import scanpy as sc
import pandas as pd

# Demo dataset
adata = sc.datasets.pbmc3k() 

# ranking genes between groups
sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')

# recalling the function
scores_and_pvals = get_scores_and_pvals(adata)
print(scores_and_pvals)
