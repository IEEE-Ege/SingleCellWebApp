# umap_strategy.py
from __future__ import annotations

from dataclasses import dataclass

import umap  # pip install umap-learn

from core import VizContext, VizParams
from backend import PlotBackend

#hesaplama ve çizim kısımları ayrıldı
"""umap yerine farklı başka algoritmalarla değişitirilebilir
    → PacMAP
    → TriMAP
    → PHATE
    → veya GPU UMAP gibi şeyler istersek burası değişir """
"""backend → çizim layer
core → parametre & context
service → pipeline
strategy → çizilecek metod
transform → ham hesaplama"""

@dataclass
class UMAPTransform:
    """
    UMAP dönüşümünü yapan sınıf (SRP – sadece embedding).
    model ayarlama 
    fit_transform() hesaplama 
    UMAP API çağrısı 
    
    """
    n_neighbors: int = 15
    min_dist: float = 0.1
    metric: str = "euclidean"
    random_state: int = 42

    def fit_transform(self, X):
        reducer = umap.UMAP(
            n_neighbors=self.n_neighbors,
            min_dist=self.min_dist,
            metric=self.metric,
            random_state=self.random_state,
        )
        embedding = reducer.fit_transform(X)
        return embedding


#    hangi görselleştirme yöntemi kullanılacak
class UMAPStrategy:
    """
    Sadece UMAP çizimi yapan strateji.
    OOP olsun diye Strategy pattern'in tek örneği gibi düşünebilirsin.
    """

    def __init__(self, transform: UMAPTransform | None = None) -> None:
        self.transform = transform

    def prepare(self, ctx: VizContext) -> None:
        """
        UMAP embedding'i hesapla ve ctx içine koy.
        Hata olursa raise ETMİYOR, sadece print yapıyor.
        """
        if self.transform is None:
            params: VizParams = ctx.params
            self.transform = UMAPTransform(
                n_neighbors=params.n_neighbors,
                min_dist=params.min_dist,
                metric=params.metric,
                random_state=params.random_state,
            )

        try:
            embedding = self.transform.fit_transform(ctx.data)
            ctx.embedding = embedding
        except Exception as e:
            # Burada kullanıcıyı kırmamak için sadece mesaj basıyoruz.
            print(f"[UMAPStrategy] UMAP hesaplanırken hata oluştu: {e}")
            ctx.embedding = None

    def render(self, ctx: VizContext, backend: PlotBackend) -> None:
        """
        Backend-agnostic render: sadece backend.draw_umap çağırıyoruz.
        """
        if ctx.embedding is None:
            print("[UMAPStrategy] Embedding yok, çizim yapılmadı.")
            return

        p = ctx.params
        backend.draw_umap(
            embedding=ctx.embedding,
            labels=ctx.labels,
            point_size=p.point_size,
            alpha=p.alpha,
            cmap=p.cmap,
            color_labels=p.color_labels,
        )
