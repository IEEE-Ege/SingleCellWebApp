import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
import os

# === Plot helper ===
def save_and_show(filename, folder='plots') -> None:
    """
    Saves the current matplotlib plot to a file and displays it.

    Parameters:
    filename (str): Name of the output image file 
    folder (str): Folder where the image will be saved (default is "plots").
    """
    if not os.path.exists(folder):
        os.makedirs(folder)
    plt.savefig(f"{folder}/{filename}")
    plt.show()

def load_file_demo(file_path: str) -> sc.AnnData:
    """
    Loads a .h5ad file, applies basic quality control preprocessing,
    and plots the top expressed genes.

    Steps:
    - Reads file and saves a copy
    - Plots highest expressed genes
    - Filters cells and genes
    - Annotates mitochondrial genes
    - Calculates QC metrics
    - Normalizes and log-transforms data
    - Identifies highly variable genes

    Parameters:
    file_path (str): Path to the .h5ad file

    Returns:
    AnnData: Preprocessed AnnData object
    """
    adata = sc.read_h5ad(file_path)
    adata.write("adata_demo.h5ad")
    print("Data saved as adata_demo.h5ad")

    print("Plotting highest expressed genes...")
    sc.pl.highest_expr_genes(adata, n_top=20, show=False)
    save_and_show("highest_expr_genes.png")

    print("Filtering cells and genes...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    print("Annotating mitochondrial genes...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')

    print("Calculating QC metrics...")
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    print("Normalizing and log-transforming...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    print("Finding highly variable genes...")
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    print("Preprocessing done.")
    return adata

def pca_plot(adata: sc.AnnData, log: bool = True, color: str = 'CST3') -> None:
    """
    Performs PCA on the dataset and visualizes results.

    Parameters:
    adata (AnnData): Preprocessed AnnData object
    log (bool): Whether to log scale the variance ratio plot
    color (str): Gene or metadata field to color PCA by
    """
    sc.tl.pca(adata, svd_solver='arpack')

    sc.pl.pca_variance_ratio(adata, log=log, show=False)
    save_and_show("pca_variance_ratio.png")

    sc.pl.pca(adata, color=color, show=False)
    save_and_show("pca_plot.png")

def highvarGen_pl(adata: sc.AnnData) -> None:
    """
    Plots various visualizations related to highly variable genes and quality metrics.

    - UMAP by cell type (if available)
    - Highly variable genes plot
    - Violin plots for quality control stats
    """
    if "Cell type" in adata.obs:
        sc.pl.umap(adata, color="Cell type", show=False)
        save_and_show("umap_cell_type.png")
    else:
        print("Error: 'Cell type' variable couldn't find in adata.obs.")

    sc.pl.highly_variable_genes(adata, show=False)
    save_and_show("highly_variable_genes.png")

    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True, show=False)
    save_and_show("violin_qc.png")

def scatter_plotting(adata: sc.AnnData) -> None:
    """
    Generates scatter plots for common QC visualizations:

    - UMAP by cell type (if available)
    - Total counts vs. mitochondrial percentage
    - Total counts vs. number of genes
    """
    if "Cell type" in adata.obs:
        sc.pl.umap(adata, color="Cell type", show=False)
        save_and_show("umap_cell_type_scatter.png")
    else:
        print("Error: 'Cell type' variable couldn't find in adata.obs.")

    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', show=False)
    save_and_show("scatter_total_vs_pct_mt.png")

    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show=False)
    save_and_show("scatter_total_vs_n_genes.png")

# === MAIN RUN ===
file_path = r"C:\Users\gediz\Downloads\Hw3covid_Data_AllCells.h5ad"

adata = load_file_demo(file_path)

pca_plot(adata)
highvarGen_pl(adata)
scatter_plotting(adata)