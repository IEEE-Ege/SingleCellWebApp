import scanpy as sc
import matplotlib.pyplot as plt

def run_pca(adata, n_comps=50, svd_solver='arpack'):
    """Runs PCA and stores the computed components."""
    try:
        if 'X_pca' in adata.obsm:
            print("PCA already computed. Overwriting previous results...")
        
        print(f"Running PCA with {n_comps} components using {svd_solver} solver...")
        sc.pp.normalize_total(adata, target_sum=1e4)  # Normalization
        sc.pp.log1p(adata)  # Log transformation
        sc.pp.scale(adata)  # Scaling
        sc.tl.pca(adata, n_comps=n_comps, svd_solver=svd_solver)
        
        print("PCA completed.")
    except Exception as e:
        print(f"Error during PCA: {e}")
        raise  

def plot_pca(adata, color=None):
    """Plots PCA results, colored by a specified attribute (if provided)."""
    try:
        print(f"Plotting PCA, color by: {color or 'default'}")
        sc.pl.pca(adata, color=color, show=False)
        plt.title(f"PCA - Colored by {color if color else 'default'}")
        plt.show()
    except KeyError:
        print(f"Warning: '{color}' not found. Using default coloring.")
        sc.pl.pca(adata, show=False)
        plt.title("PCA - Default Coloring")
        plt.show()
    except Exception as e:
        print(f"Error in PCA plot: {e}")
        raise  

def plot_variance(adata, log=True):
    """Plots the variance explained by PCA components."""
    try:
        print("Plotting explained variance...")
        sc.pl.pca_variance_ratio(adata, log=log, show=False)
        plt.title("PCA: Explained Variance")
        plt.show()
    except Exception as e:
        print(f"Error in variance plot: {e}")
        raise  

def save_results(adata, results_file="pca_results.h5ad"):
    """Saves the PCA results to an H5AD file."""
    try:
        print(f"Saving results to {results_file}...")
        adata.write(results_file)
        print("Save successful.")
    except Exception as e:
        print(f"Error saving results: {e}")
        raise  

def get_adata(adata):
    """Returns the processed AnnData object."""
    return adata
