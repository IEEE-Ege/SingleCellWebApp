import scanpy as sc
import pandas as pd

adata="Hw3covid_Data_AllCells.h5ad"

sc.pp.neighbors(adata)
sc.tl.leiden(adata, resolution=1.0)
sc.settings.verbosity = 2  # reduce the verbosity


def get_marker_genes(adata):
   adata = sc.read(results_file)
   pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
   return df 



