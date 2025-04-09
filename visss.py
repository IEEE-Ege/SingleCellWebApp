import scanpy as sc
import matplotlib.pyplot as plt

def load_file_demo(file_path):
    adata = sc.read_h5ad(file_path)
    adata.write("adata_demo.h5ad")
    print("Data saved as adata_demo.h5ad")
    pca_plot(adata)
    highvarGen_pl(adata)
    scatter_plotting(adata)

    sc.pl.highest_expr_genes(adata, n_top=20)
    
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)

    print("Preprocessing done.")
    return adata

def pca_plot(adata,log=True,color= 'CST3'):
    sc.pp.pca(adata)
    sc.pl.pca_variance_ratio(adata, log=log)
    sc.pl.pca(adata, color) 
    plt.show()
def highvarGen_pl(adata):
    if "Cell type" in adata.obs:
        sc.pl.umap(adata, color="Cell type")
    else:
        print("Error: 'Cell type' variable couldn't find in adata.obs.")
    
    sc.pl.highly_variable_genes(adata, show=False)
    plt.show()
    
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True, show=False)
    plt.show()

def scatter_plotting(adata):
    if "Cell type" in adata.obs:
        sc.pl.umap(adata, color="Cell type", show=False)
    else:
        print("Error: 'Cell type' variable couldn't find in adata.obs.")
        
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
    plt.show()
# --- Call the function here ---
file_path = "Hw3covid_Data_AllCells.h5ad"
adata = load_file_demo(file_path)

# Then call the other plotting functions as needed:
# pca_plot(adata)
# highvarGen_pl(adata)
# scatter_plotting(adata)
