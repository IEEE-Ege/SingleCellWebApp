import numpy as np 
import pandas as pd
from anndata import AnnData
import scanpy as sc

# 1. Create dummy data
X = np.random.rand(100, 10)  # 100 cells, 10 genes
var_names = [f"Gene{i}" for i in range(10)]
obs = pd.DataFrame({
    # Create a 'leiden' column with categorical values ('0', '1', '2')
    'leiden': pd.Categorical(np.random.choice(['0', '1', '2'], size=100))
})
# Create the AnnData object
adata = AnnData(X=X, obs=obs)
adata.var_names = var_names # Set the gene names

print("Dummy AnnData object created with 'leiden' clusters.")
print(f" Original cluster categories:\n{adata.obs['leiden'].cat.categories}")


# 2. Define the new cluster names
# Ensure the number of new names matches the number of categories in 'leiden'
new_cluster_names = ['CD4 T cells', 'Monocytes', 'B cells']
print(f"\nDefined new cluster names: {new_cluster_names}")


# 3. Define the rename_clusters function
def rename_clusters(adata, cluster_algo: str, new_cluster_names: list):
    """
    Renames cluster labels in adata.obs using a chosen clustering algorithm key and new names.

    Parameters:
        adata (AnnData): The annotated data matrix.
        cluster_algo (str): The key in `adata.obs` representing clustering (e.g., 'leiden', 'louvain').
        new_cluster_names (list): A list of new cluster names.

    Returns:
        None: Updates adata in place.
    """
    print(f"Attempting to rename clusters for column '{cluster_algo}'...")

    # Check if the cluster key exists in adata.obs
    if cluster_algo not in adata.obs:
        raise ValueError(f"'{cluster_algo}' not found in adata.obs. Make sure clustering has been performed.")

    # Check if the column is categorical, as rename_categories works on categorical columns
    if not pd.api.types.is_categorical_dtype(adata.obs[cluster_algo]):
        raise TypeError(f"'{cluster_algo}' column must be categorical.")

    # Get the number of existing clusters (categories)
    n_clusters = len(adata.obs[cluster_algo].cat.categories)

    # Check if the number of new names matches the number of clusters
    if len(new_cluster_names) != n_clusters:
        raise ValueError(f"Expected {n_clusters} new names for '{cluster_algo}', but got {len(new_cluster_names)}.")

    # Perform the renaming
    adata.rename_categories(cluster_algo, new_cluster_names)
    print(f" Cluster names updated successfully in '{cluster_algo}'.")


# 4. Apply the function
print("\nApplying rename_clusters function...")
rename_clusters(adata, cluster_algo='leiden', new_cluster_names=new_cluster_names)


# 5. See the results
print("\n Updated cluster names:")
print(adata.obs['leiden'].cat.categories)