import scanpy as sc
import numpy as np
import pandas as pd
from preprocessing import filter_data, calculate_qc_metrics, select_highly_variable_genes

def load_file_demo(file_path):
    adata = sc.read_h5ad(file_path)
    adata.write("adata_demo.h5ad")
    print("Data saved as adata_demo.h5ad")

file_path = "C:/Users/pc/OneDrive/Masaüstü/bioinformatic/Hw3covid_Data_AllCells.h5ad"
adata = load_file_demo(file_path)

# Step 1: Create example data (random counts matrix)
np.random.seed(42)  # For reproducibility

n_cells = 1000  # Number of cells
n_genes = 5000  # Number of genes

# Create random gene expression data (Poisson distribution)
counts = np.random.poisson(lam=1.0, size=(n_cells, n_genes))  # Shape: (n_cells, n_genes)

# Create an AnnData object
adata = sc.AnnData(counts)

# Assign gene and cell names
adata.var_names = [f"Gene_{i}" for i in range(n_genes)]  # Gene names: Gene_0, Gene_1, ..., Gene_4999
adata.obs_names = [f"Cell_{i}" for i in range(n_cells)]  # Cell names: Cell_0, Cell_1, ..., Cell_999

# Before adding 'MT-', note that var_names is not a list but an Index object.
# So we need to convert var_names to a list before making changes.
var_names_list = adata.var_names.tolist()

# Add some mitochondrial genes for demonstration (every 10th gene is mitochondrial)
for i in range(0, n_genes, 10):
    var_names_list[i] = f"MT-{var_names_list[i]}"  # Rename every 10th gene to start with "MT-"

# We assign the modified list back as an Index object
adata.var_names = pd.Index(var_names_list)

# Step 2: Apply filter_data to filter cells and genes
adata_filtered = filter_data(adata, min_genes=200, min_cells=3)
print(f"Filtered AnnData: {adata_filtered.shape} (cells, genes)")

# Step 3: Calculate QC metrics
adata_qc = calculate_qc_metrics(adata_filtered)
print("\nQC Metrics (first few rows of obs):")
print(adata_qc.obs.head())  # View the QC metrics such as 'pct_counts_mt' (percentage of mitochondrial genes)

# Step 4: Select highly variable genes
adata_hvg = select_highly_variable_genes(adata_qc, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("\nHighly Variable Genes (first few rows of 'highly_variable' column):")
print(adata_hvg.var[['highly_variable']].head())  # Display which genes are highly variable

# Step 5: Count how many genes are highly variable
num_highly_variable_genes = adata_hvg.var['highly_variable'].sum()
print(f"\nNumber of highly variable genes selected: {num_highly_variable_genes}")
