
import scanpy as sc 

sc.pp.neighbors(adata)
sc.tl.leiden(adata, resolution=1.0)
sc.settings.verbosity = 2  # reduce the verbosity

def rank_genes_groups(adata, method='t-test', n_genes=25, sharey=False):
    """
    -Perform differential expression analysis and plot top marker genes.

    Parameters:
    - adata: AnnData object
    - method: Method for differential testing ('t-test', 'wilcoxon', 'logreg', etc.)
    - n_genes: Number of top genes to plot
    - sharey: Whether to share the y-axis across plots
    """

    sc.tl.rank_genes_groups(adata, groupby='leiden', method=method)
    sc.pl.rank_genes_groups(adata, n_genes=n_genes, sharey=sharey)
    adata.write(results_file)
