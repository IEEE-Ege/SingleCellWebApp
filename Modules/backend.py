# backend.py çizim nasıl yapılacak
#Adapter Pattern uygulaması
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Optional, Sequence

import matplotlib.pyplot as plt
#Matplotlib mi, Plotly mi, Bokeh mi burası biliyor geri kalan kısımlar bilmiyor 

class PlotBackend(ABC):
    """
    UMAP için minimal backend arayüzü.
    çizimi nası yapıcaz 
    """

    @abstractmethod
    def init_canvas(self, width: float, height: float, dpi: int) -> None:
        ...

    @abstractmethod
    def draw_umap(
        self,
        embedding,
        labels: Optional[Sequence] = None,
        point_size: float = 5.0,
        alpha: float = 0.8,
        cmap: str = "viridis",
        color_labels: bool = True,
    ) -> None:
        ...

    @abstractmethod
    def save(self, path: str) -> None:
        ...

    @abstractmethod
    def show(self) -> None:
        ...

    @abstractmethod
    def get_figure(self) -> Any:
        ...

    @abstractmethod
    def get_axes(self) -> Any:
        ...


class MatplotlibBackend(PlotBackend):
    """
    Matplotlib adapter'ı.
    """
    #canvas oluşturma 
    def __init__(self) -> None:
        self._fig: Any = None
        self._ax: Any = None

    def init_canvas(self, width: float, height: float, dpi: int) -> None:
        self._fig, self._ax = plt.subplots(figsize=(width, height), dpi=dpi)

    def draw_umap(
        self,
        embedding,
        labels: Optional[Sequence] = None,
        point_size: float = 5.0,
        alpha: float = 0.8,
        cmap: str = "viridis",
        color_labels: bool = True,
    ) -> None:
        if self._ax is None:
            self.init_canvas(8, 6, 120)

        x = embedding[:, 0]
        y = embedding[:, 1]

        if labels is None:
            sc = self._ax.scatter(x, y, s=point_size, alpha=alpha, cmap=cmap)
        else:
            if color_labels:
                sc = self._ax.scatter(
                    x,
                    y,
                    c=labels,
                    s=point_size,
                    alpha=alpha,
                    cmap=cmap,
                )
            else:
                sc = self._ax.scatter(
                    x,
                    y,
                    s=point_size,
                    alpha=alpha,
                )

        self._ax.set_xlabel("UMAP-1")
        self._ax.set_ylabel("UMAP-2")
        self._ax.set_title("UMAP Projection")

        # Etiketler sayısal ise colorbar güzel olur
        if labels is not None and color_labels:
            self._fig.colorbar(sc, ax=self._ax, label="labels")

    def save(self, path: str) -> None:
        if self._fig is not None:
            self._fig.savefig(path, bbox_inches="tight")

    def show(self) -> None:
        if self._fig is not None:
            plt.show()

    def get_figure(self) -> Any:
        return self._fig

    def get_axes(self) -> Any:
        return self._ax
