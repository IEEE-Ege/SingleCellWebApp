# core.py
# Veri akışı ve konfigürasyonun standartlaştırılması

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Optional, Sequence, Dict, List


@dataclass
class VizParams:
    """
    UMAP + plot + kayıt + canvas ayarlarını tutan config sınıfı.
    Strategy, backend, service hep bu obje üzerinden konuşur.
    """

    # ---------- UMAP ayarları ----------
    n_neighbors: int = 15
    min_dist: float = 0.1
    metric: str = "euclidean"
    random_state: Optional[int] = 42
    n_components: int = 2  # genelde 2D UMAP

    # ---------- Canvas / figür boyutları ----------
    width: float = 8.0    # inches (matplotlib)
    height: float = 6.0
    dpi: int = 120

    # ---------- Nokta / scatter ayarları ----------
    point_size: float = 5.0
    alpha: float = 0.8

    # ---------- Renk / etiket ayarları ----------
    cmap: str = "viridis"      # Eğer labels sayısal ise
    color_labels: bool = True  # Kategorik/etiket bazlı renklendirme
    label_name: Optional[str] = None  # Örn: "cluster", "louvain" vs.

    # ---------- Kayıt / output ayarları ----------
    save_dir: str = "figures"
    save_name: str = "umap"
    save_format: str = "png"
    save_enabled: bool = True  # SaveStrategy kullanırken işine yarar

    # ---------- Eksen / başlık ----------
    title: Optional[str] = "UMAP Projection"
    x_label: str = "UMAP-1"
    y_label: str = "UMAP-2"
    show_axis: bool = True

    def update(self, **kwargs: Any) -> None:
        """
        ****
        UMAP parametrelerini sonradan değiştirme imkanı.
        Parametreleri runtime'da değiştirmek için:
            ctx.params.update(n_neighbors=30, min_dist=0.05, title="Deneme UMAP")
        Yanlış isim verilirse AttributeError fırlatır ki hatayı erken görelim.
        """
        for key, value in kwargs.items():
            if not hasattr(self, key):
                raise AttributeError(f"VizParams has no attribute '{key}'")
            setattr(self, key, value)


@dataclass
class VizContext:
    """
    Farklı strategy'ler arasında ortak kullanılacak context yapısı.
    Data, parametreler, embedding, backend ve hataları burada tutar.
    """
    # Girdi tarafı
    data: Any
    labels: Optional[Sequence] = None
    params: VizParams = field(default_factory=VizParams)

    # Hesaplama sonrası doldurulacak
    embedding: Optional[Any] = None

    # Backend tarafının dolduracağı alanlar
    figure: Any = None
    ax: Any = None
    backend: Any = None  # Örn: MatplotlibBackend instance

    # Hata ve ek artefakt yönetimi
    errors: List[str] = field(default_factory=list)
    artifacts: Dict[str, Any] = field(default_factory=dict)

    def add_error(self, msg: str) -> None:
        """Context içinde hata mesajı biriktir."""
        self.errors.append(msg)

    @property
    def has_errors(self) -> bool:
        return len(self.errors) > 0

    def attach_backend(self, backend: Any) -> None:
        """Backend instance'ını context'e kaydet.
        service backend'i context ile ilişkilendirir.****
        """
        self.backend = backend

    def to_dict(self) -> Dict[str, Any]:
        """
        Debug / logging için context özetini sözlük formatında döner.
        AnnData veya büyük objeleri komple dump etmiyoruz, sadece boyut / key bilgisi.
        """
        return {
            "params": asdict(self.params),
            "has_errors": self.has_errors,
            "errors": list(self.errors),
            "embedding_shape": getattr(self.embedding, "shape", None),
            "artifacts": list(self.artifacts.keys()),
        }
