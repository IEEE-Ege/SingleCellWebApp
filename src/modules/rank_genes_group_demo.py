import scanpy as sc

# using sample data
adata = sc.datasets.pbmc3k()

# calculating neighborhoods
sc.pp.neighbors(adata)

# applying leiden 
sc.tl.leiden(adata, resolution=1.0)
sc.settings.verbosity = 2



# running the function
rank_genes_groups(adata, method='t-test', n_genes=10, sharey=True, results_file='ranked_genes_results.h5ad')

# to check if the file saved succesfully
print("Results saved to 'ranked_genes_results.h5ad'")
