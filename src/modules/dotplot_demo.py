import scanpy as sc
import warnings
import pandas as pd

def visualize_marker_genes(adata, marker_genes, cluster_key='leiden'):
    """
    Visualizes marker gene expression across annotated cell clusters.

    Parameters:
        adata (AnnData): Annotated data matrix containing single-cell RNA-seq data.
        marker_genes (list or dict): Marker genes to visualize. Can be a list or a dict (for grouped markers).
        cluster_key (str): Key in `adata.obs` to group cells by (default: 'leiden').

    Returns:
        None: Displays dotplot and stacked violin plot of marker gene expression.
    """
    # Check if cluster_key exists in adata.obs columns
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"'{cluster_key}' column not found in adata.obs. Available columns: {list(adata.obs.columns)}")

    # Validate marker genes depending on type (list or dict)
    if isinstance(marker_genes, dict):
        # Flatten the list of genes from the dictionary values for validation
        flat_genes = [gene for genes in marker_genes.values() for gene in genes]
    elif isinstance(marker_genes, list):
        flat_genes = marker_genes
    else:
        # Raise an error if marker_genes is not a list or dictionary
        raise TypeError("marker_genes must be a list or a dictionary of gene lists.")

    # Check gene existence in adata.var_names
    valid_genes = []
    missing_genes = []
    for gene in flat_genes:
        if gene in adata.var_names:
            valid_genes.append(gene)
        else:
            missing_genes.append(gene)

    # Warn the user about any genes that were not found
    if missing_genes:
        warnings.warn(f"The following genes are not found in adata.var_names and will be ignored: {missing_genes}")

    # If no valid genes are left after checking, raise an error
    if not valid_genes:
        raise ValueError("None of the provided marker genes were found in adata.var_names.")

    # Filter marker_genes based on valid_genes
    if isinstance(marker_genes, dict):
        # Filter genes within each group in the dictionary
        marker_genes = {k: [g for g in v if g in valid_genes] for k, v in marker_genes.items()}
        # Remove any groups that become empty after filtering
        marker_genes = {k: v for k, v in marker_genes.items() if v}

    else:
        # If it was a list, just use the list of valid genes
        marker_genes = valid_genes

    # Generate the dotplot
    print("Generating dotplot...")
    sc.pl.dotplot(adata, marker_genes, groupby=cluster_key)
    print("Dotplot generated.")

    # Generate the stacked violin plot
    print("Generating stacked violin plot...")
    sc.pl.stacked_violin(adata, marker_genes, groupby=cluster_key, rotation=90)
    print("Stacked violin plot generated.")

# Note: These plotting functions typically display the plots automatically in interactive environments
# If running as a script, you might need to add:
# import matplotlib.pyplot as plt
# plt.show()