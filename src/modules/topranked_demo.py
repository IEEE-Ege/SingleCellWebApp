import pandas as pd

def get_top_ranked_genes(adata, top_n=5):
    """
    Retrieves the top `top_n` marker genes for each cluster from
    `adata.uns['rank_genes_groups']`.

    Parameters:
    - adata (AnnData): Annotated data matrix after differential expression analysis (`sc.tl.rank_genes_groups`).
    - top_n (int): Number of top-ranked genes to retrieve for each cluster.

    Returns:
    - pd.DataFrame: Long-format DataFrame with 'cluster', 'rank', 'gene_name', 'pval', and 'score'.
    """
    # Check if 'rank_genes_groups' exists in adata.uns
    if 'rank_genes_groups' not in adata.uns:
        raise ValueError("No 'rank_genes_groups' found in adata.uns. Run `sc.tl.rank_genes_groups()` first.")

    # Validate the top_n parameter
    if not isinstance(top_n, int) or top_n <= 0:
        raise ValueError("`top_n` must be a positive integer.")

    # Get the results dictionary
    result = adata.uns['rank_genes_groups']
    # Extract cluster names from the result structure
    groups = result['names'].dtype.names

    data = [] # List to store data for the DataFrame
    # Iterate through each cluster
    for cluster in groups:
        # Extract top_n gene names for the current cluster
        gene_names = result['names'][cluster][:top_n]
        # Extract top_n p-values, using .get to handle potential absence of 'pvals'
        pvals = result.get('pvals', {}).get(cluster, [None]*top_n)
        # Extract top_n scores, using .get to handle potential absence of 'scores'
        scores = result.get('scores', {}).get(cluster, [None]*top_n)

        # Iterate through the top genes and their associated values for the cluster
        for rank, (gene, pval, score) in enumerate(zip(gene_names, pvals, scores), 1):
            # Append a dictionary for each gene to the data list
            data.append({
                'cluster': cluster,
                'rank': rank,
                'gene_name': gene,
                'pval': pval,
                'score': score
            })

    # Create and return a pandas DataFrame from the collected data
    return pd.DataFrame(data)