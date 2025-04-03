import scanpy as sc
import matplotlib.pyplot as plt

class PCAAnalyzer:
    def __init__(self, adata, results_file="pca_results.h5ad"):
        """Class for performing PCA analysis on AnnData objects."""
        self.adata = adata
        self.results_file = results_file
    
    def run_pca(self, n_comps=50, svd_solver='arpack'):
        """Runs PCA and stores the computed components."""
        try:
            if 'X_pca' in self.adata.obsm:
                print("PCA already computed. Overwriting previous results...")
            
            print(f"Running PCA with {n_comps} components using {svd_solver} solver...")
            sc.pp.normalize_total(self.adata, target_sum=1e4)  # Normalization
            sc.pp.log1p(self.adata)  # Log transformation
            sc.pp.scale(self.adata)  # Scaling
            sc.tl.pca(self.adata, n_comps=n_comps, svd_solver=svd_solver)
            
            print("PCA completed.")
        except Exception as e:
            print(f"Error during PCA: {e}")
            raise  

    def plot_pca(self, color=None):
        """Plots PCA results, colored by a specified attribute (if provided)."""
        try:
            print(f"Plotting PCA, color by: {color or 'default'}")
            sc.pl.pca(self.adata, color=color, show=False)
            plt.title(f"PCA - Colored by {color if color else 'default'}")
            plt.show()
        except KeyError:
            print(f"Warning: '{color}' not found. Using default coloring.")
            sc.pl.pca(self.adata, show=False)
            plt.title("PCA - Default Coloring")
            plt.show()
        except Exception as e:
            print(f"Error in PCA plot: {e}")
            raise  
    
    def plot_variance(self, log=True):
        """Plots the variance explained by PCA components."""
        try:
            print("Plotting explained variance...")
            sc.pl.pca_variance_ratio(self.adata, log=log, show=False)
            plt.title("PCA: Explained Variance")
            plt.show()
        except Exception as e:
            print(f"Error in variance plot: {e}")
            raise  
    
    def save_results(self):
        """Saves the PCA results to an H5AD file."""
        try:
            print(f"Saving results to {self.results_file}...")
            self.adata.write(self.results_file)
            print("Save successful.")
        except Exception as e:
            print(f"Error saving results: {e}")
            raise  
    
    def get_adata(self):
        """Returns the processed AnnData object."""
        return self.adata