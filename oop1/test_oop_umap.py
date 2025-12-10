# test_oop_umap.py

import scanpy as sc
import matplotlib.pyplot as plt

from core import VizParams
from service import VisualizationService
from backend import MatplotlibBackend
from save_strategy import NoSaveStrategy


# -------------------------------------------------
# 0) Opsiyonel: Özel backend (sağ panel için aks paylaşımı)
# -------------------------------------------------
class AxSharingBackend(MatplotlibBackend):
    """
    Bu backend, yeni bir figure/axes oluşturmak yerine
    dışarıdan verilen bir axes üzerine çizer.

    Amaç: Solda Scanpy, sağda OOP UMAP olacak şekilde
    aynı figürde 1x2 subplot içinde karşılaştırma yapmak.
    """

    def __init__(self, ax):
        super().__init__()
        self._ax = ax
        self._fig = ax.figure

    def init_canvas(self, width: float, height: float, dpi: int) -> None:
        """
        VisualizationService init_canvas çağırdığı için
        bu metot override edilip 'hiçbir şey yapmıyor'.
        Fig/Ax zaten dışarıdan geldi.
        """
        # Canvas zaten var, hiçbir şey yapmıyoruz.
        return


# -------------------------------------------------
# 1) Veriyi oku
# -------------------------------------------------
# Gerekirse bu path'i kendi dosya yoluna göre değiştir:
adata = sc.read_h5ad("/Users/birsen/Desktop/oop1/pbmc3k_raw.h5ad")

print(adata)
print(adata.obs.head())


# -------------------------------------------------
# 2) Preprocessing (Scanpy klasik PBMC3k pipeline)
# -------------------------------------------------

# Hücre ve gen filtreleme
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)

# log1p
sc.pp.log1p(adata)

# Highly variable genler
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
)

# Sadece HVG'ler
adata = adata[:, adata.var["highly_variable"]].copy()

# Ölçekleme
sc.pp.scale(adata, max_value=10)

# PCA
sc.tl.pca(adata, svd_solver="arpack")

# Komşuluk grafiği
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)


# -------------------------------------------------
# 3) Scanpy ile UMAP + Leiden cluster
# -------------------------------------------------
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# Scanpy'nin hesapladığı UMAP embedding
scanpy_umap = adata.obsm["X_umap"]
leiden = adata.obs["leiden"]


# -------------------------------------------------
# 4) Karşılaştırma figürü: 1x2 subplot
# -------------------------------------------------
fig, (ax_left, ax_right) = plt.subplots(
    1, 2,
    figsize=(12, 6),
    dpi=120
)

# --- Sol panel: Scanpy UMAP ---
sc_left = ax_left.scatter(
    scanpy_umap[:, 0],
    scanpy_umap[:, 1],
    c=leiden.astype("category").cat.codes,
    s=5,
    alpha=0.8,
    cmap="viridis",
)
ax_left.set_title("Scanpy UMAP (sc.tl.umap)")
ax_left.set_xlabel("UMAP-1")
ax_left.set_ylabel("UMAP-2")
fig.colorbar(sc_left, ax=ax_left, label="leiden (Scanpy)")


# --- Sağ panel: Senin OOP UMAP + VisualizationService ---
# Backend olarak AxSharingBackend kullanıyoruz ki sağ panelin ax'ine çizsin
backend = AxSharingBackend(ax_right)

# İstersen burada UMAP parametrelerini de oynayabilirsin:
params = VizParams(
    n_neighbors=15,
    min_dist=0.1,
    metric="euclidean",
    random_state=42,
    point_size=5.0,
    alpha=0.8,
    cmap="viridis",
    color_labels=True,
    title="OOP UMAP (VisualizationService + UMAPStrategy)",
    save_enabled=False,   # bu testte dosya kaydetme
)

# Varsayılan kayıt davranışı yerine hiçbir şey kaydetmeyen strategy
service = VisualizationService(
    backend=backend,
    save_strategy=NoSaveStrategy(),
)

ctx = service.run_umap(
    data=adata.obsm["X_pca"],   # 👈 UMAP PCA uzayında
    labels=leiden.cat.codes,
    params=params,
    show=False,                 # figürün tamamını en sonda plt.show ile açıyoruz
)

# Sağ panel başlığı (params.title zaten backend'de de kullanılabilir)
ax_right.set_title(params.title or "OOP UMAP")


# -------------------------------------------------
# 5) Sonuçları göster
# -------------------------------------------------
plt.tight_layout()
plt.show()

# Eğer istersen figürü dışarı kaydet:
# fig.savefig("compare_scanpy_vs_oop_umap.png", dpi=150, bbox_inches="tight")
