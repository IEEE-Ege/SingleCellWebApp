# save_strategy.py
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Optional
import os

from core import VizContext
from backend import PlotBackend


class SaveStrategy(ABC):
    """
    Kaydetme davranışını soyutlayan Strategy.
    Kaydetme politikası ayrı bir Strategy
    İster hiç kaydetme, ister otomatik dosya adı üret, ister özel path kullan.
    """

    @abstractmethod
    def save(
        self,
        ctx: VizContext,
        backend: PlotBackend,
        save_path: Optional[str] = None,
    ) -> None:
        ...


class NoSaveStrategy(SaveStrategy):
    """Hiçbir şey kaydetmez, sadece figürü context'te bırakır."""

    def save(
        self,
        ctx: VizContext,
        backend: PlotBackend,
        save_path: Optional[str] = None,
    ) -> None:
        # Sadece figür referansını artifacts içine atalım
        ctx.artifacts["fig"] = backend.get_figure()


class AutoPathSaveStrategy(SaveStrategy):
    """
    Kaydetme için otomatik path üretir:
    - Eğer save_path verilmişse onu kullanır
    - Verilmemişse ctx.params içindeki save_dir, save_name, save_format üzerinden üretir
    """

    def save(
        self,
        ctx: VizContext,
        backend: PlotBackend,
        save_path: Optional[str] = None,
    ) -> None:
        params = ctx.params

        # Kullanıcı kaydetmek istemiyorsa çık
        if not params.save_enabled:
            return

        # Dışarıdan özel path verilmişse onu kullan
        if save_path is None:
            os.makedirs(params.save_dir, exist_ok=True)
            filename = f"{params.save_name}.{params.save_format}"
            save_path = os.path.join(params.save_dir, filename)

        try:
            backend.save(save_path)
            ctx.artifacts["fig_path"] = save_path
        except Exception as e:
            ctx.add_error(f"Save error: {e}")
