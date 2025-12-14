from prefect import flow, task


from louv_yeni_son import (
    upload_dataset as backend_upload_dataset,
    store_dataset as backend_store_dataset,
    quality_control as backend_quality_control,
    normalization as backend_normalization,
    log_gc_errors as backend_log_gc_errors,
    dimensionality_reduction_pca as backend_dimred_pca,
    clustering as backend_clustering,
    save_results as backend_save_results,
    generate_umap_visualization as backend_generate_umap_visualization,
    render_interactive_umap as backend_render_interactive_umap,
    generate_output_report as backend_generate_output_report,
)

# === KONFİG  ===
DATA_PATH = "./COVID.h5ad"   # TODO: kendi dosya path'inle değiştir
MIN_GENES = 200
MAX_GENES = 5000
HVG_METHOD = "seurat_v3"
N_TOP_GENES = 2000
PERPLEXITY = 30.0
CLUSTERING_METHOD = "leiden"    # veya "louvain"
RESOLUTION = 0.8
OUTPUT_H5AD_PATH = "data/output_processed.h5ad"  # istersen None yap


@task
def upload_dataset(path: str):
    return backend_upload_dataset(path)

@task
def store_dataset(adata):
    return backend_store_dataset(adata)

@task(retries=3, retry_delay_seconds=5)
def quality_control(adata):
    return backend_quality_control(adata, MIN_GENES, MAX_GENES)

@task
def normalization(adata):
    return backend_normalization(adata, HVG_METHOD, N_TOP_GENES)

@task
def log_gc_errors():
    backend_log_gc_errors()

@task(retries=3, retry_delay_seconds=5)
def dimensionality_reduction_pca(adata):
    return backend_dimred_pca(adata, PERPLEXITY)

@task()
def clustering(adata):
    return backend_clustering(adata, CLUSTERING_METHOD, RESOLUTION)

@task
def save_results(adata):
    backend_save_results(adata, OUTPUT_H5AD_PATH)

@task
def generate_umap_visualization(adata):
    backend_generate_umap_visualization(adata)

@task
def render_interactive_umap(adata):
    backend_render_interactive_umap(adata)

@task
def generate_output_report():
    backend_generate_output_report()


@flow
def my_flow():
    # 1) Upload & store
    adata = upload_dataset(DATA_PATH)
    adata = store_dataset(adata)

    # 2) QC
    try:
        adata = quality_control(adata)
    except Exception as e:
        print(f"QC hata verdi: {e}")
        log_gc_errors()
        print("QC tekrar denendi ama başarısız oldu. Normalization ve downstream task'ler devam edecek.")

    # 3) Normalization + HVG
    adata = normalization(adata)

    # 4) PCA + UMAP + t-SNE + clustering
    try:
        adata = dimensionality_reduction_pca(adata)
        adata = clustering(adata)
    except Exception as e:
        print(f"PCA/clustering hata verdi: {e}")
        print("PCA tekrar denendi ama başarısız oldu, clustering atlandı.")

    # 5) Son adımlar
    # save_results(adata)
    generate_umap_visualization(adata)
    render_interactive_umap(adata)
    generate_output_report()


if __name__ == "__main__":
    my_flow()
