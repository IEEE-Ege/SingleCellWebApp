import pandas as pd

def extract_top_genes(adata, n_top=5):
    """
    Extracts top-ranked genes and p-values from adata.uns['rank_genes_groups'].

    Parameters:
    - adata: AnnData object after rank_genes_groups has been run.
    - n_top: Number of top genes to return (default is 5).

    Returns:
    - A tuple of two DataFrames:
        1. Top gene names per group.
        2. A combined DataFrame of top genes and p-values.
    """

    if 'rank_genes_groups' not in adata.uns:
        raise ValueError("No 'rank_genes_groups' results found in adata.uns. Did you run `sc.tl.rank_genes_groups`?")

    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names

    top_gene_names = pd.DataFrame(result['names']).head(n_top)

    combined_df = pd.DataFrame({
        f"{group}_n": result['names'][group][:n_top]
        for group in groups
    })
    combined_df_pvals = pd.DataFrame({
        f"{group}_p": result['pvals'][group][:n_top]
        for group in groups
    })

    combined = pd.concat([combined_df, combined_df_pvals], axis=1)

    return top_gene_names, combined


