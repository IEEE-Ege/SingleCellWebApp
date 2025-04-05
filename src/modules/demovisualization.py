import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

def load_file_demo(file_path):
    adata = sc.read_h5ad(file_path)
    adata.write("adata_demo.h5ad")
    print("Data saved as adata_demo.h5ad")

    print("Plotting highest expressed genes...")
    sc.pl.highest_expr_genes(adata, n_top=20, show=False)
    plt.show()

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

# PCA and plotting functions
def pca_plot(adata, log=True, color='CST3'):
    # Run PCA if not already done
    sc.tl.pca(adata, svd_solver='arpack')

    # Plot PCA variance ratio
    sc.pl.pca_variance_ratio(adata, log=log, show=False)
    plt.show()

    # Plot PCA
    sc.pl.pca(adata, color=color, show=False)
    plt.show()

def highvarGen_pl(adata):
    if "Cell type" in adata.obs:
        sc.pl.umap(adata, color="Cell type", show=False)
        plt.show()
    else:
        print("Error: 'Cell type' variable couldn't find in adata.obs.")
    
    sc.pl.highly_variable_genes(adata, show=False)
    plt.show()
    
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True, show=False)
    plt.show()

def scatter_plotting(adata):
    if "Cell type" in adata.obs:
        sc.pl.umap(adata, color="Cell type", show=False)
        plt.show()
    else:
        print("Error: 'Cell type' variable couldn't find in adata.obs.")
        
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', show=False)
    plt.show()

    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show=False)
    plt.show()

# === MAIN RUN ===
file_path = "Hw3covid_Data_AllCells.h5ad"
adata = load_file_demo(file_path)

pca_plot(adata)
highvarGen_pl(adata)
scatter_plotting(adata)
