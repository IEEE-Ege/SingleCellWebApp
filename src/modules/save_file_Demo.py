from load_file import load_file
from save import save_result
file_path = "/Users/azratuncay/Downloads/filtered_gene_bc_matrices/hg19"
adata = load_file(file_path) 
results_file = 'pbmc3k.h5ad' 
save_result(adata, results_file)
results_file

