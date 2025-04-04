import scanpy as sc

def preprocess(adata, target_sum = 1e4, max_value = 10):
    data_normalization(adata, target_sum, max_value)
    #other preprocessing methods here

def data_normalization(adata, target_sum, max_value):
    sc.pp.normalize_total(adata, target_sum)
    sc.pp.log1p(adata)
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value)