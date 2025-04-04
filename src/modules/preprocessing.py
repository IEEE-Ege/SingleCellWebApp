import scanpy as sc

def normalize_and_scale(adata, target_sum, max_value):
    sc.pp.normalize_total(adata, target_sum)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value)

    vars_to_check = ["total_counts", "pct_counts_mt"]
    try:
        for var in vars_to_check: #checking for the var names if adata has them
            if var not in adata.obs.columns:
                raise ValueError(f"{var} is not found in adata.var_names")
            
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt']) #if adata has the parameters do regress out
    except Exception as e:
        print(f"An error has occurred: {e}") 


def preprocess(adata, target_sum = 1e4, max_value = 10):
    normalize_and_scale(adata, target_sum, max_value)
    #other preprocessing methods here
    return adata
