# service.py
from __future__ import annotations

from typing import Any, Optional, Sequence, Any as AnyType

from core import VizContext, VizParams
from backend import MatplotlibBackend, PlotBackend
from umap_strategy import UMAPStrategy
from save_strategy import SaveStrategy, AutoPathSaveStrategy, NoSaveStrategy


class VisualizationService:
    # FACADE PATTERN – kullanıcının gördüğü yüz
    """
    Facade:
    Kullanıcı sadece run_umap(...) çağırır, gerisi içeride halledilir.
    - Modülerlik
    - Genişletilebilirlik
    - Ayarları değiştirme kolaylığı
    """

    def __init__(
        self,
        backend: Optional[PlotBackend] = None,
        umap_strategy: Optional[UMAPStrategy] = None,
        save_strategy: Optional[SaveStrategy] = None,
    ) -> None:
        # Varsayılan backend ve stratejiler
        self.backend: PlotBackend = backend or MatplotlibBackend()
        self.umap_strategy: UMAPStrategy = umap_strategy or UMAPStrategy()
        # Varsayılan save strategy: parametrelere göre otomatik path
        self.save_strategy: SaveStrategy = save_strategy or AutoPathSaveStrategy()

    def _build_context(
        self,
        data: Any,
        labels: Optional[Sequence] = None,
        params: Optional[VizParams] = None,
        param_overrides: Optional[dict[str, AnyType]] = None,
    ) -> VizContext:
        """
        ****Context oluşturma ve parametre override işini ayrı fonksiyona aldık.
        run_umap çağrısında ek parametre verdiğinde (n_neighbors=30 gibi),
        onları tek tek VizParams’a elinle aktarmak zorunda değilsin → otomatik update.
        """
        ctx = VizContext(
            data=data,
            labels=labels,
            params=params or VizParams(),
        )

        # UMAP / plot parametrelerini runtime'da değiştirmek için
        if param_overrides:
            ctx.params.update(**param_overrides)

        return ctx

    def run_umap(
        self,
        data: Any,
        labels: Optional[Sequence] = None,
        params: Optional[VizParams] = None,
        save_path: Optional[str] = None,
        show: bool = True,
        save_strategy: Optional[SaveStrategy] = None,
        **param_overrides: AnyType,
    ) -> VizContext:
        """
        Uçtan uca UMAP akışı:
        - Context oluştur
        - Parametre override (n_neighbors, min_dist, title, save_name vs.)
        - UMAPStrategy.prepare → embedding hesapla
        - Backend.init_canvas → çizim alanı
        - UMAPStrategy.render → scatter çizimi
        - SaveStrategy.save → dosyaya kaydet / kaydetme / context'e sadece fig koy
        - İstenirse show()

        Örnek kullanım:
            service = VisualizationService()
            ctx = service.run_umap(
                data=X,
                labels=clusters,
                n_neighbors=30,
                min_dist=0.05,
                title="High-res UMAP",
                save_name="umap_highres"
            )
        """
        # 1) Context + parametre override
        ctx = self._build_context(
            data=data,
            labels=labels,
            params=params,
            param_overrides=param_overrides or None,
        )

        # Backend'i context'e iliştir (debug vs. için)
        ctx.attach_backend(self.backend)

        # 2) UMAP embedding hazırla
        self.umap_strategy.prepare(ctx)

        # Hesaplama sırasında hata olduysa çizime geçme
        if ctx.has_errors:
            # Burada istersen loglama vs. de yapabilirsin
            return ctx

        # 3) Canvas oluştur
        p = ctx.params
        self.backend.init_canvas(p.width, p.height, p.dpi)

        # 4) Çizim
        self.umap_strategy.render(ctx, self.backend)

        # 5) Kaydetme (SaveStrategy ile)
        effective_save_strategy = save_strategy or self.save_strategy
        effective_save_strategy.save(ctx, self.backend, save_path=save_path)

        # 6) Gösterme
        if show:
            try:
                self.backend.show()
            except Exception as e:
                ctx.add_error(f"Show error: {e}")

        # 7) Figure ve axes referanslarını context'e koy
        ctx.figure = self.backend.get_figure()
        ctx.ax = self.backend.get_axes()
        ctx.artifacts.setdefault("backend", self.backend)

        return ctx
