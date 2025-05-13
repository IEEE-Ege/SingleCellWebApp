import scanpy as sc
import Neighborhood as nb  # Your functions are here!
import umap

#Load the dataset
print("Loading dataset...")
adata = sc.datasets.pbmc3k()
print(f"Dataset loaded! Number of cells: {adata.n_obs}, Number of genes: {adata.n_vars}\n")

#Compute the neighbor graph
print("Computing the neighborhood graph...")
nb.compute(adata, n_neighbors=10, n_pcs=40)

#Compute UMAP
print("Computing UMAP embedding...")
sc.tl.umap(adata)

#Visualize the embedding
print("Plotting UMAP...")
genes_of_interest = ['CD3D', 'MS4A1', 'GNLY']
nb.embed(adata, color=genes_of_interest)

#Perform clustering
print("Performing Leiden clustering...")
nb.cluster(adata, color=['leiden'] + genes_of_interest)

#Save the results
results_file = "pbmc3k_final_results.h5ad"
print(f"Saving results to {results_file}...")
nb.save(results_file, adata)
