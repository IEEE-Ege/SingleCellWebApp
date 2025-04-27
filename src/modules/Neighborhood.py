import shiny  #Only needed if you are planning to use it later
import scanpy as sc
import umap

adata = sc.read("C:/Users/pc/.ipython/HW3/Hw3covid_Data_AllCells.h5ad")

#computing
def compute(adata, n_neighbors=10, n_pcs=40):
    """_summary_

    Args:
        Hw3covid_Data_AllCells (_type_): _description_
        n_neighbors (int, optional): _description_. Defaults to 10.
        n_pcs (int, optional): _description_. Defaults to 40.
    """
    try:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        print("Neighbors computed successfully.")
    except Exception as e:
        print(f"Error in computing neighbors: {e}")

#embedding
def embed(adata, color=['CST3', 'NKG7', 'PPBP']):
    """_summary_

    Args:
        adata (_type_): _description_
        color (list, optional): _description_. Defaults to ['CST3', 'NKG7', 'PPBP'].
    """
    try:
        sc.pl.umap(adata, color=color)
        sc.pl.umap(adata, color=color)  # Note: This plots twice; check if you really want this
        print("Embedding and plotting successful.")
    except Exception as e:
        print(f"Error in embedding: {e}")

#clustering
def cluster(adata, color=['leiden', 'CST3', 'NKG7']):
    """_summary_

    Args:
        adata (_type_): _description_
        color (list, optional): _description_. Defaults to ['leiden', 'CST3', 'NKG7'].
    """
    try:
        sc.tl.leiden(adata)
        sc.pl.umap(adata, color=color)
        print("Clustering and plotting successful.")
    except Exception as e:
        print(f"Error in clustering: {e}")

#saving
def save(results_file, adata):
    """_summary_

    Args:
        results_file (_type_): _description_
        adata (_type_): _description_
    """
    try:
        adata.write(results_file)
        print(f"Data saved successfully to {results_file}.")
    except Exception as e:
        print(f"Error in saving data: {e}")
