from prefect import flow, task
import random
import louv_yeni_son as app # Mevcut uygulamayı içe aktar

@task
def upload_dataset():
    app.upload_dataset()
    print("Uploading dataset...")

@task
def store_dataset():
    app.store_dataset()
    print("Storing dataset...")

@task(retries=3, retry_delay_seconds=5)
def quality_control() -> bool:
    app.quality_control()
    print("Performing quality control...")
    qc_passed = random.random() > 0.5
    if not qc_passed:
        print("QC failed!")
        raise ValueError("QC failed!")
    print("QC succeeded!")
    return True

@task
def normalization():
    app.normalization()
    print("Normalizing data...")

@task
def log_gc_errors():
    app.log_gc_errors()
    print("Logging GC errors...")

@task(retries=3, retry_delay_seconds=5)
def dimensionality_reduction_pca() -> bool:
    app.dimensionality_reduction_pca()
    print("Performing PCA for dimensionality reduction...")
    if random.random() < 0.5:
        print("PCA failed!")
        raise ValueError("PCA failed!")
    print("PCA succeeded!")
    return True

@task()
def clustering():
    app.clustering()
    print("Clustering data...")

@task
def save_results():
    app.save_results()
    print("Saving results...")

@task
def generate_umap_visualization():
    app.generate_umap_visualization()
    print("Generating UMAP visualization...")

@task
def render_interactive_umap():
    app.render_interactive_umap()
    print("Rendering interactive UMAP...")

@task
def generate_output_report():
    app.generate_output_report()
    print("Generating output report...")

@flow
def my_flow():
    upload_dataset()
    store_dataset()
    
    try:
        qc_success = quality_control().result()
    except Exception:
        log_gc_errors()
        print("QC tekrar denendi ama başarısız oldu. Normalization ve downstream task'ler devam edecek.")
    
    normalization()
    
    try:
        pca_success = dimensionality_reduction_pca().result()
        if pca_success:
            clustering()
    except Exception:
        print("PCA tekrar denendi ama başarısız oldu, clustering atlandı.")
    
    save_results()
    generate_umap_visualization()
    render_interactive_umap()
    generate_output_report()

if __name__ == "__main__":
    my_flow()

