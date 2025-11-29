# service.py
from __future__ import annotations

from typing import Any, Optional, Sequence

from core import VizContext, VizParams
from backend import MatplotlibBackend, PlotBackend
from umap_strategy import UMAPStrategy
"Kullanıcı için API basitleştirme"
"Parametre karmaşasını ayrı katmana taşır strategy context vs"
"Backend değiştirmeyi kolaylaştırır service "
"plot ekleme kolaylığı sağlar"

class VisualizationService: #içinde olanları dışarıya yansıtmaz
    #FACADE PATTERN kullanıcının gördüğü yüz 
    """
    Facade:
    Kullanıcı sadece run_umap(...) çağırır, gerisi içeride halledilir.
    modullerdik, genişletilebilrlik ve ayarları değiştirme kolaylığı sağlar.
    """

    def __init__(
        self,
        backend: Optional[PlotBackend] = None,
        strategy: Optional[UMAPStrategy] = None,
    ) -> None:
        self.backend = backend or MatplotlibBackend()
        self.strategy = strategy or UMAPStrategy()

    def run_umap(
        self,
        data: Any,
        labels: Optional[Sequence] = None,
        params: Optional[VizParams] = None,
        save_path: Optional[str] = None,
        show: bool = True,
    ) -> VizContext:
        """
        Uçtan uca:
        - Context oluştur
        - UMAPStrategy.prepare (embedding)
        - backend.init_canvas çizim alanı oluşturur 
        - UMAPStrategy.render
        - save/show
        """
        ctx = VizContext(
            data=data,
            labels=labels,
            params=params or VizParams(),
        )

        # 1) UMAP hazırla
        self.strategy.prepare(ctx)

        # 2) Canvas
        p = ctx.params
        self.backend.init_canvas(p.width, p.height, p.dpi)

        # 3) Çiz
        self.strategy.render(ctx, self.backend)

        # 4) Kaydet ve/veya göster
        if save_path:
            try:
                self.backend.save(save_path)
            except Exception as e:
                print(f"[VisualizationService] Kaydetme sırasında hata: {e}")

        if show:
            try:
                self.backend.show()
            except Exception as e:
                print(f"[VisualizationService] Gösterme sırasında hata: {e}")

        # Context'i geri döndürüyorum, embedding vs. istersen kullan.
        ctx.figure = self.backend.get_figure()
        ctx.ax = self.backend.get_axes()
        return ctx
