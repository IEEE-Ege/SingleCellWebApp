import scanpy as sc
import pandas as pd

# Load the dataset
adata = sc.datasets.pbmc3k()

# Log-transform the raw count data (required for differential expression analysis)
sc.pp.log1p(adata)

# Apply PCA and calculate neighbors
sc.pp.pca(adata, svd_solver='arpack')  # Apply PCA
sc.pp.neighbors(adata, n_pcs=40)  # Calculate neighbors after PCA

# Run Leiden algorithm (with flavor='igraph' for future compatibility)
sc.tl.leiden(adata, resolution=1.0, flavor="igraph", directed=False, n_iterations=2)

# Run differential expression analysis to rank genes
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')

# Save the adata object with the rank_genes_groups results
adata.write("ranked_genes_results.h5ad")  # Save the results to a file

# Set verbosity level
sc.settings.verbosity = 2

# Define the function to get marker genes
def get_marker_genes(adata, results_file):
    # Read the results file
    adata = sc.read(results_file)
    
    # Check if 'rank_genes_groups' exists in adata.uns
    if 'rank_genes_groups' in adata.uns:
        # Extract the top 5 marker genes as a DataFrame
        df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
        return df
    else:
        print("Error: 'rank_genes_groups' not found in adata.uns.")
        return None

# Run the function to get marker genes
marker_genes = get_marker_genes(adata, results_file="ranked_genes_results.h5ad")
if marker_genes is not None:
    print(marker_genes)
