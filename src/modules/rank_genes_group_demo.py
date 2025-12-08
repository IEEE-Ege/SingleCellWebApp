from rankgenesgroup import rank_genes_groups
import scanpy as sc

adata = sc.datasets.pbmc3k()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)

# PCA and neighborhoods
sc.pp.pca(adata, svd_solver='arpack')  # Bu satır uyarıyı çözer
sc.pp.neighbors(adata, n_pcs=40)

# Clustering (Leiden)
sc.tl.leiden(adata, resolution=1.0)

sc.settings.verbosity = 2

#to save results
results_file = "pbmc3k_rank_genes.h5ad"
sc.tl.rank_genes_groups(adata, groupby='leiden' , method='t-test' , n_genes=25 , sharey=False)
adata.write(results_file)
