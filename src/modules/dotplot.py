import scanpy as sc

def visualize_marker_genes(adata, marker_genes, cluster_key='leiden'):
    """
    Visualizes marker gene expression across annotated cell clusters.

    Parameters:
        adata (AnnData): Annotated data matrix containing single-cell RNA-seq data.
        marker_genes (list or dict): Marker genes to visualize.
        cluster_key (str): Key in `adata.obs` to group cells by (default: 'leiden').

    Returns:
        None: Displays dotplot and stacked violin plot of marker gene expression.
    """

    # Dotplot to show average expression and percent of expressing cells
    sc.pl.dotplot(adata, marker_genes, groupby=cluster_key)

    # Stacked violin plot to show gene expression distribution per cluster
    sc.pl.stacked_violin(adata, marker_genes, groupby=cluster_key, rotation=90)



