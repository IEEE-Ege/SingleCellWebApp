

import scanpy as sc
import matplotlib.pyplot as plt

class PCAAnalyzer:
    def __init__(self, adata, results_file="pca_results.h5ad"):
        """Class initializer for PCA analysis."""
        self.adata = adata
        self.results_file = results_file

   # n_comps: explains how many components are going to be calculated by pca
   # svd_solver: explains which algorithm is going to be used when calculating pca
    def run_pca(self, n_comps=50, svd_solver='arpack'): # calculated pca components
        """Runs PCA analysis and stores the computed components."""
        try:
            print(f"Running PCA (n_comps={n_comps}, solver={sed_solver})...")
            sc.tl.pca(self.adata, n_comps=n_comps, svd_solver=svd_solver)
            print("PCA completed successfully!")
            print("Generated fields:")
            print(f"- Cell PCs: adata.obsm['X_pca'] (shape: {self.adata.obsm['X_pca'].shape})")
            print(f"- Gene PCs: adata.varm['PCs'] (shape: {self.adata.varm['PCs'].shape})")
        except Exception as e:
            print(f"PCA error: {str(e)}")
            raise  # in case of error it helps by launching the error

    def plot_pca(self, color='CST3'): # pca scatter plot
        """Visualizes PCA results with a specified coloring.
        If the specified attribute is missing, defaults to standard coloring.
        """
        try:
            print(f"Plotting PCA (color: {color})...")
            sc.pl.pca(self.adata, color=color, show=False)
            plt.title(f"PCA - Colored by {color}")
            plt.show()
        except KeyError:
            print(f"Warning: '{color}' attribute not found in the dataset, using default coloring.")
            sc.pl.pca(self.adata, show=False)
            plt.title("PCA - Default Coloring")
            plt.show()
        except Exception as e:
            print(f"Visualization error: {str(e)}")
            raise  

    def plot_variance(self, log=True): # pca variance ratio bar graph
        """Plots the variance explained by PCA components."""
        try:
            print("Plotting variance ratio...")
            sc.pl.pca_variance_ratio(self.adata, log=log, show=False)
            plt.title("PCA: Explained Variance by Components")
            plt.show()
        except Exception as e:
            print(f"Variance plot error: {str(e)}")
            raise  

    def save_results(self):
        """Saves the analysis results to the specified file."""
        try:
            print(f"Saving results to '{self.results_file}'...")
            self.adata.write(self.results_file)
            print("Save successful!")
        except Exception as e:
            print(f"Save error: {str(e)}")
            raise  

    def get_adata(self):
        """Returns the processed AnnData object."""
        return self.adata