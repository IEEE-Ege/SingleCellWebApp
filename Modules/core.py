# core.py
#veri akışı standartlaştırma 
from dataclasses import dataclass, field
from typing import Any, Optional, Sequence
#parametreleri kolayca değiştirme ve tutma için 
#UMAP + plot + canvas ayarlarını tutan config sınıfı.


@dataclass #detaylı düşün??????
class VizParams: #tüm parametreleri tutar 
    """
    Sadece UMAP için basit parametre seti.
    Strategy, backend, service hep aynı objeyle konuşsun.
    Data + sonuçlar tek yerde olsun
    """
    n_neighbors: int = 15
    min_dist: float = 0.1
    metric: str = "euclidean"
    random_state: int = 42

    width: float = 8.0    # in inches (matplotlib)
    height: float = 6.0
    dpi: int = 120

    point_size: float = 5.0
    alpha: float = 0.8

    cmap: str = "viridis"  # Eğer labels numeric ise
    # Eğer kategorik etiketlerle renklendirme istersen:
    color_labels: bool = True

#farklı strategyler için ortak context yapısı
@dataclass
class VizContext:
    """
    UMAP çizimi için ihtiyaç olunan veriyi tutar.
    """
    data: Any
    labels: Optional[Sequence] = None
    params: VizParams = field(default_factory=VizParams)

    # prepare sonrası doldurulacak:
    embedding: Optional[Any] = None

    # backend tarafı dolduracak:
    figure: Any = None
    ax: Any = None
