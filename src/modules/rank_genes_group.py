import scanpy as sc

# Örnek veri setini yükle
adata = sc.datasets.pbmc3k()

# Temel ön işleme
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)

# PCA ve komşuluk hesaplama
sc.pp.pca(adata, svd_solver='arpack')  # Bu satır uyarıyı çözer
sc.pp.neighbors(adata, n_pcs=40)

# Clustering (Leiden)
sc.tl.leiden(adata, resolution=1.0)

# Verbosity ayarı
sc.settings.verbosity = 2

# Sonuçları kaydetmek için bir dosya ismi belirle
results_file = "pbmc3k_rank_genes.h5ad"

# Gen sıralama ve görselleştirme fonksiyonu
def rank_genes_groups(adata, method='t-test', n_genes=25, sharey=False):
    """
    - Perform differential expression analysis and plot top marker genes.
    """
    sc.tl.rank_genes_groups(adata, groupby='leiden', method=method)
    sc.pl.rank_genes_groups(adata, n_genes=n_genes, sharey=sharey)
    adata.write(results_file)

# Fonksiyonu çağır
rank_genes_groups(adata)


