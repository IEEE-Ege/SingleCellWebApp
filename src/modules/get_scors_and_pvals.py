
import scanpy as sc
import pandas as pd

def get_scores_and_pvals(adata):
   result = adata.uns['rank_genes_groups']
   groups = result['names'].dtype.names
   df=pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(5)
   return df.head(5)