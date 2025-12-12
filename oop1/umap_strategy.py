# umap_strategy.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import umap  # pip install umap-learn

from core import VizContext, VizParams
from backend import PlotBackend


"""
Katmanlar:
- transform  → sadece sayısal UMAP hesaplama (embedding)
- strategy   → embedding + backend ile çizim (UMAP'e özel)
core         → parametre & context
backend      → çizim layer (Matplotlib, Plotly vs.)
service      → pipeline / facade
"""


@dataclass
class UMAPTransform:
    """
    sadece matematik/algoritma***
    UMAP dönüşümünü yapan sınıf (SRP – sadece embedding).
    - model ayarlama
    - fit_transform() ile hesaplama
    """

    n_neighbors: int = 15
    min_dist: float = 0.1
    metric: str = "euclidean"
    random_state: int | None = 42
    n_components: int = 2

    def fit_transform(self, X: Any):
        reducer = umap.UMAP(
            n_neighbors=self.n_neighbors,
            min_dist=self.min_dist,
            metric=self.metric,
            random_state=self.random_state,
            n_components=self.n_components,
        )
        embedding = reducer.fit_transform(X)
        return embedding

    @classmethod
    def from_params(cls, params: VizParams) -> "UMAPTransform":
        """
        bu embedding ile ne yapıyoruz backend’le çizim, context yönetimi.
        VizParams içinden UMAPTransform objesi üretmek için yardımcı constructor.
        Parametreleri sonradan değiştirdiğinde (params.update(...)) buraya otomatik yansır.
        """
        return cls(
            n_neighbors=params.n_neighbors,
            min_dist=params.min_dist,
            metric=params.metric,
            random_state=params.random_state,
            n_components=params.n_components,
        )


class UMAPStrategy:
    """
    Sadece UMAP çizimi yapan strateji.
    İleride PacMAP / TriMAP / PHATE gibi başka transform'larla da değiştirilebilir.
    """

    def __init__(self, transform: UMAPTransform | None = None) -> None:
        self.transform = transform

    def prepare(self, ctx: VizContext) -> None:
        """
        UMAP embedding'i hesapla ve ctx.embedding içine koy.
        Hata olursa context'e error ekler, embedding None kalır.
        """
        # Transform objesi yoksa, her çağrıda güncel VizParams üzerinden oluştur.
        if self.transform is None:
            self.transform = UMAPTransform.from_params(ctx.params)
        else:
            # Dışarıdan inject edilen transform varsa ama params güncellendiyse,
            # senkronize etmek istiyorsan burada da update edebilirsin (opsiyonel).
            self.transform.n_neighbors = ctx.params.n_neighbors
            self.transform.min_dist = ctx.params.min_dist
            self.transform.metric = ctx.params.metric
            self.transform.random_state = ctx.params.random_state
            self.transform.n_components = ctx.params.n_components

        try:
            embedding = self.transform.fit_transform(ctx.data)
            ctx.embedding = embedding
            ctx.artifacts["umap_embedding"] = embedding
        except Exception as e:
            ctx.add_error(f"UMAP computation error: {e}")
            ctx.embedding = None

    def render(self, ctx: VizContext, backend: PlotBackend) -> None:
        """
    
        Backend-agnostic render: sadece backend.draw_umap çağırıyoruz.
        """
        if ctx.embedding is None:
            ctx.add_error("[UMAPStrategy] Embedding yok, çizim yapılmadı.")
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
