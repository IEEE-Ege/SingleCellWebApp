from rank_genes_violin import get_rank_genes_groups_violin
import scanpy as sc

#sample dataset
adata = sc.datasets.pbmc3k()

# calculating neighborhoods
sc.pp.neighbors(adata)

sc.tl.leiden(adata, resolution=1.0)

# Differential expression analysis (with t-test)
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')


# calling the function and visualizing
get_rank_genes_groups_violin(adata, groups='0', n_genes=8)
