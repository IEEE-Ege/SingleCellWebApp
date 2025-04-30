import scanpy as sc
import pandas as pd
from rank_genes_group import rank_genes_groups

# Load the dataset
adata = sc.datasets.pbmc3k()

# Apply PCA to reduce dimensionality (to avoid high-dimension warning)
sc.pp.pca(adata, svd_solver='arpack')

# Calculate neighbors on PCA-reduced data
sc.pp.neighbors(adata)

# Run Leiden algorithm (with future-proofing for 'igraph' backend)
sc.tl.leiden(adata, resolution=1.0, flavor='igraph', directed=False, n_iterations=2)

# Run differential expression analysis to rank genes
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')


# Reduce verbosity for cleaner output
sc.settings.verbosity = 2

# Define a function to get gene names and p-values from rank_genes_groups
def get_scores_and_pvals(adata):
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    df = pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
         for group in groups for key in ['names', 'pvals']}).head(5)
    return df

# Run the function to get the results
scores_and_pvals = get_scores_and_pvals(adata)
print(scores_and_pvals)