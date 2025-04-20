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

    if cluster_algo not in adata.obs:
        raise ValueError(f"'{cluster_algo}' not found in adata.obs. Make sure clustering has been run.")

    n_clusters = len(adata.obs[cluster_algo].cat.categories)

    if len(new_cluster_names) != n_clusters:
        raise ValueError(f"Number of new names ({len(new_cluster_names)}) does not match number of clusters ({n_clusters}).")

    adata.rename_categories(cluster_algo, new_cluster_names)
    print(f"Clusters renamed using '{cluster_algo}'.")


