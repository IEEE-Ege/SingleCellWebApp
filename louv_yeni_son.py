# -*- coding: utf-8 -*-
"""Scanpy Analysis App (Updated with New Marker Genes)"""

from shiny import App, ui, render, reactive
from shiny.types import FileInfo
import matplotlib.pyplot as plt
import scanpy as sc

# === ORTAK PIPELINE FONKSİYONLARI (Hem Prefect hem Shiny kullanacak) ===

def _load_adata(path: str):
    """H5AD dosyasını yükle ve var_names'i string'e çevir."""
    adata = sc.read_h5ad(path)
    adata.var_names = adata.var_names.astype(str)
    return adata

def _run_quality_control(adata, min_genes: int, max_genes: int):
    """QC: min_genes ve max_genes filtresi uygula."""
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs.n_genes < max_genes, :].copy()
    return adata

def _run_normalization_and_hvg(adata, hvg_method: str, n_top_genes: int):
    """Normalize, log-transform, HVG seçimi ve scale."""
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        flavor=hvg_method,
        n_top_genes=n_top_genes
    )
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    return adata

def _run_dimensionality_reduction(adata, perplexity: float):
    """PCA, komşular, UMAP ve t-SNE."""
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.tsne(adata, perplexity=perplexity)
    return adata

def _run_clustering(adata, method: str, resolution: float):
    """Leiden veya Louvain clustering."""
    if method == "louvain":
        sc.tl.louvain(adata, resolution=resolution, random_state=42)
    else:
        sc.tl.leiden(adata, resolution=resolution, random_state=42)
    return adata

def _save_results_to_disk(adata, out_path: str | None = None):
    """
    Sonuçları diske kaydetmek için basit bir örnek.
    Şimdilik opsiyonel; istersen path verebilirsin.
    """
    if out_path is not None:
        adata.write(out_path)


# === BACKEND FONKSİYONLAR (Prefect buradakileri çağıracak) ===
# Buradakiler GERÇEK pipeline adımlarını çağırıyor, sadece print değil.

def upload_dataset(path: str):
    print("[louv_yeni] upload_dataset backend fonksiyonu çalıştı.")
    adata = _load_adata(path)
    return adata

def store_dataset(adata):
    """
    Şimdilik sadece geçiş yapıyor.
    İstersen burada geçici dosyaya yazma vs. yapabilirsin.
    """
    print("[louv_yeni] store_dataset backend fonksiyonu çalıştı.")
    return adata

def quality_control(adata, min_genes: int, max_genes: int):
    print("[louv_yeni] quality_control backend fonksiyonu çalıştı.")
    adata = _run_quality_control(adata, min_genes, max_genes)
    # Burada istersen "hiç hücre kalmadı" kontrolü yapıp hata fırlatabilirsin
    if adata.n_obs == 0:
        raise ValueError("QC sonrası hiç hücre kalmadı.")
    return adata

def normalization(adata, hvg_method: str, n_top_genes: int):
    print("[louv_yeni] normalization backend fonksiyonu çalıştı.")
    adata = _run_normalization_and_hvg(adata, hvg_method, n_top_genes)
    return adata

def log_gc_errors():
    # QC veya başka bir adım hata verirse Prefect burayı çağıracak
    print("[louv_yeni] log_gc_errors backend fonksiyonu çalıştı. (Burada log dosyasına yazılabilir)")

def dimensionality_reduction_pca(adata, perplexity: float):
    print("[louv_yeni] dimensionality_reduction_pca backend fonksiyonu çalıştı.")
    adata = _run_dimensionality_reduction(adata, perplexity)
    return adata

def clustering(adata, method: str, resolution: float):
    print("[louv_yeni] clustering backend fonksiyonu çalıştı.")
    adata = _run_clustering(adata, method, resolution)
    return adata

def save_results(adata, out_path: str | None = None):
    print("[louv_yeni] save_results backend fonksiyonu çalıştı.")
    _save_results_to_disk(adata, out_path)

def generate_umap_visualization(adata):
    """
    Prefect için placeholder. İstersen burada figürü kaydedebilirsin.
    Shiny tarafı zaten ayrı plot ediyor.
    """
    print("[louv_yeni] generate_umap_visualization backend fonksiyonu çalıştı.")

def render_interactive_umap(adata):
    """
    Prefect için placeholder. Örn. bir HTML dosyası üretilebilir.
    """
    print("[louv_yeni] render_interactive_umap backend fonksiyonu çalıştı.")

def generate_output_report():
    """
    Prefect için placeholder. Örn. bir PDF/HTML raporu oluşturulabilir.
    """
    print("[louv_yeni] generate_output_report backend fonksiyonu çalıştı.")



# === SHINY UI KISMI ===

# UI
app_ui = ui.page_fluid(
    ui.layout_sidebar(
        ui.sidebar(
            ui.panel_title("SCA-Web"),
            ui.input_file("file", "Upload .h5ad file", accept=[".h5ad"]),
            ui.input_slider("min_genes", "Minimum genes per cell", min=0, max=500, value=200),
            ui.input_slider("max_genes", "Maximum genes per cell", min=500, max=10000, value=5000),
            ui.input_select("hvg_method", "Highly Variable Gene Selection Method",
                            choices=["seurat_v3", "seurat", "cell_ranger"],
                            selected="seurat_v3"),
            ui.input_slider("n_top_genes", "Top HVGs to keep", min=100, max=5000, value=2000),
            ui.input_select("reduction_method", "Dimensionality Reduction Method",
                            choices=["pca"], selected="pca"),
            ui.input_slider("perplexity", "t-SNE Perplexity", min=5, max=50, value=30),
            ui.input_action_button("run", "Run Analysis", class_="btn-success"),
            ui.panel_title("Clustering"),
            ui.input_select("clustering_method", "Clustering Algorithm",
                            choices=["leiden", "louvain"], selected="leiden"),
            ui.input_slider("resolution", "Resolution", min=0.1, max=2.0, value=0.8),
            ui.input_action_button("find_clusters", "Find Clusters", class_="btn-primary")
        ),
        [
            ui.card(
                ui.card_header("Dataset Summary"),
                ui.output_text_verbatim("summary"),
                id="summary_card"
            ),
            ui.card(
                ui.card_header("Highly Variable Genes"),
                ui.output_plot("hvg_plot")
            ),
            ui.card(
                ui.card_header("UMAP Plot"),
                ui.output_plot("umap_plot")
            ),
            ui.card(
                ui.card_header("t-SNE Plot"),
                ui.output_plot("tsne_plot")
            ),
            ui.card(
                ui.card_header("Clustered UMAP Plot"),
                ui.output_plot("clustered_umap_plot")
            ),
            ui.card(
                ui.card_header("Clustered t-SNE Plot"),
                ui.output_plot("clustered_tsne_plot")
            ),
            ui.card(
                ui.card_header("Marker Gene MatrixPlot"),
                ui.output_plot("marker_matrix_plot")
            )
        ]
    ),
    ui.TagList(
        ui.tags.style("""
            .card {
                background-color: #ffffff;
                box-shadow: 0 4px 8px rgba(0,0,0,0.1);
                border-radius: 8px;
                margin-bottom: 20px;
                padding: 20px;
            }
            .card-header {
                font-size: 1.2rem;
                font-weight: 600;
                background-color: #4CAF50;
                color: white;
                padding: 10px;
                border-top-left-radius: 8px;
                border-top-right-radius: 8px;
            }
            .shiny-input-slider {
                margin-bottom: 15px;
            }
            .card-body {
                padding: 10px;
                display: flex;
                justify-content: center;
                align-items: center;
            }
            .shiny-output-container {   
                margin: 0 auto;
            }
            #summary {
                justify-content: flex-end;
            }
            .shiny-image-output {
                max-width: 700px;
                max-height: 600px;
                    }
        """)
    )
)

# Server logic
def server(input, output, session):

    # NEW Marker Gene Dictionary
    marker_genes_dict = {
        'B-cell': ['CD79A', 'MS4A1'],
        'T-cell': ['CD3D'],
        'T-cell CD8+': ['CD8A', 'CD8B'],
        'NK': ['GNLY', 'NKG7'],
        'Myeloid': ['CST3', 'LYZ'],
        'Monocytes': ['FCGR3A'],
        'Dendritic': ['FCER1A']
    }

    @reactive.event(input.run)
    def process_file():
        """
        Eski haliyle bütün pipeline burada yazılıydı.
        Artık backend fonksiyonlarını kullanıyoruz ki Prefect ile ortak olsun.
        """
        file: list[FileInfo] | None = input.file()
        if not file:
            return None

        path = file[0]["datapath"]

        # === Prefect ile ortak backend fonksiyonları ===
        adata = upload_dataset(path)
        adata = store_dataset(adata)
        adata = quality_control(adata, input.min_genes(), input.max_genes())
        adata = normalization(adata, input.hvg_method(), input.n_top_genes())
        adata = dimensionality_reduction_pca(adata, input.perplexity())

        return adata

    @reactive.event(input.find_clusters)
    def find_clusters():
        adata = process_file()
        if adata is None:
            return None

        method = input.clustering_method()
        resolution = input.resolution()

        adata = clustering(adata, method, resolution)

        return adata

    @output
    @render.text
    def summary():
        adata = process_file()
        if adata is None:
            return "Upload a file and click Run."
        return f"{adata.shape[0]} cells, {adata.shape[1]} genes after filtering."

    @output
    @render.plot
    def hvg_plot():
        adata = process_file()
        if adata is None:
            return
        plt.figure()
        sc.pl.highly_variable_genes(adata, show=False)
        return plt.gcf()

    @output
    @render.plot
    def umap_plot():
        adata = process_file()
        if adata is None:
            return
        plt.figure()
        if "cell_type" in adata.obs.columns:
            sc.pl.umap(adata, color="cell_type", show=False)
        else:
            sc.pl.umap(adata, show=False)
        return plt.gcf()

    @output
    @render.plot
    def tsne_plot():
        adata = process_file()
        if adata is None:
            return
        plt.figure()
        if "cell_type" in adata.obs.columns:
            sc.pl.tsne(adata, color="cell_type", show=False)
        else:
            sc.pl.tsne(adata, show=False)
        return plt.gcf()

    @output
    @render.plot
    def clustered_umap_plot():
        adata = find_clusters()
        if adata is None:
            return
        plt.figure()
        clustering_key = "leiden" if "leiden" in adata.obs.columns else "louvain"
        sc.pl.umap(adata, color=clustering_key, show=False)
        return plt.gcf()

    @output
    @render.plot
    def clustered_tsne_plot():
        adata = find_clusters()
        if adata is None:
            return
        plt.figure()
        clustering_key = "leiden" if "leiden" in adata.obs.columns else "louvain"
        sc.pl.tsne(adata, color=clustering_key, show=False)
        return plt.gcf()

    @output
    @render.plot
    def marker_matrix_plot():
        adata = find_clusters()
        if adata is None:
            plt.figure()
            plt.text(0.5, 0.5, 'No data loaded.', ha='center', va='center')
            plt.axis('off')
            return plt.gcf()

        if adata.raw is None:
            plt.figure()
            plt.text(0.5, 0.5, 'No raw data available.', ha='center', va='center')
            plt.axis('off')
            return plt.gcf()

        plt.figure()
        clustering_key = "leiden" if "leiden" in adata.obs.columns else "louvain"

        fixed_marker_genes_dict = {
            group: [str(gene) for gene in genes]
            for group, genes in marker_genes_dict.items()
        }

        available_genes = set(adata.var_names)

        filtered_marker_genes_dict = {
            cluster: [gene for gene in genes if gene in available_genes]
            for cluster, genes in fixed_marker_genes_dict.items()
        }

        sc.pl.matrixplot(
            adata,
            filtered_marker_genes_dict,
            groupby=clustering_key,
            use_raw=False,
            show=False
        )
        return plt.gcf()


app = App(app_ui, server)
