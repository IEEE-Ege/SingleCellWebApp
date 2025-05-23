import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

# Step 1: Create random data
n_cells, n_genes = 100, 2000
np.random.seed(42)

data = np.random.rand(n_cells, n_genes)
adata = sc.AnnData(X=data)

adata.obs['cell_type'] = ['type1' if i < 50 else 'type2' for i in range(n_cells)]
adata.var['gene_id'] = [f"gene{i}" for i in range(n_genes)]

# Step 2: Apply PCA
def run_pca(adata, n_comps=50, svd_solver='arpack'):
    """Runs PCA and stores the computed components."""
    if 'X_pca' in adata.obsm:
        print("PCA already computed. Overwriting previous results...")

    try:
        print(f"Running PCA with {n_comps} components using {svd_solver} solver...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.scale(adata)
        sc.tl.pca(adata, n_comps=n_comps, svd_solver=svd_solver)
        print("PCA completed.")
    except Exception as e:
        raise Exception(f"PCA sırasında hata oluştu: {e}")

run_pca(adata, n_comps=10, svd_solver='arpack')

# Step 3: Plot PCA
def plot_pca(adata, color=None):
    """Plots PCA results, colored by a specified attribute (if provided)."""
    try:
        print(f"Plotting PCA, color by: {color or 'default'}")
        sc.pl.pca(adata, color=color, show=False)
        plt.title(f"PCA - Colored by {color if color else 'default'}")
        plt.show()
    except KeyError:
        raise Exception(f"'{color}' özelliği bulunamadı, renklemek için geçersiz.")
    except Exception as e:
        raise Exception(f"PCA görselleştirmesinde hata oluştu: {e}")

plot_pca(adata, color='cell_type')

# Step 4: Plot explained variance
def plot_variance(adata, log=True):
    """Plots the variance explained by PCA components."""
    try:
        print("Plotting explained variance...")
        sc.pl.pca_variance_ratio(adata, log=log, show=False)
        plt.title("PCA: Explained Variance")
        plt.show()
    except Exception as e:
        raise Exception(f"Varyans grafiğinde hata oluştu: {e}")

plot_variance(adata)

# Step 5: Save PCA results
def save_results(adata, results_file="pca_results.h5ad"):
    """Saves the PCA results to an H5AD file."""
    try:
        print(f"Saving results to {results_file}...")
        adata.write(results_file)
        print("Save successful.")
    except Exception as e:
        raise Exception(f"Sonuçları kaydederken hata oluştu: {e}")

save_results(adata, "pca_results.h5ad")

# Step 6: Retrieve processed AnnData object
def get_adata(adata):
    """Returns the processed AnnData object."""
    return adata

processed_adata = get_adata(adata)
