# flexvis_app.py
from __future__ import annotations

# --- Matplotlib'i non-interactive moda al ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import io, base64, inspect
import numpy as np
from typing import Optional, Dict, Any

# --- Dash ---
from dash import Dash, html, dcc, Input, Output

# --- CombinedPolicy/REG: iki farklı modülden biri olabilir ---
try:
    from .flexvis import CombinedPolicy, REG  # varsa bunu kullan
except ModuleNotFoundError:
    from .flexvis import CombinedPolicy, REG              # yoksa bunu dene

# =============================================================================
# Adapters (MPL -> Dash image figure)
# =============================================================================
def fig_to_base64(fig) -> str:
    """Matplotlib Figure -> base64 PNG (tarayıcıda gösterim için)"""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return "data:image/png;base64," + base64.b64encode(buf.read()).decode("utf-8")

def mpl_image_figure(img64: str) -> dict:
    """Base64 PNG'i Dash Graph'ın kabul ettiği figure=dict yapısına sar."""
    return {
        "data": [],
        "layout": {
            "images": [{
                "source": img64,
                "xref": "paper",
                "yref": "paper",
                "x": 0, "y": 1,
                "sizex": 1, "sizey": 1,
                "xanchor": "left", "yanchor": "top",
                "layer": "below"
            }],
            "xaxis": {"visible": False, "range": [0, 1]},
            "yaxis": {"visible": False, "range": [0, 1], "scaleanchor": "x"},
            "margin": {"l": 0, "r": 0, "t": 0, "b": 0},
        }
    }

def to_dash_figure(fig_obj):
    """
    Plotly fig ise direkt döndür, Matplotlib fig ise base64'e çevirip imaj olarak döndür.
    """
    # Plotly fig kontrolü: obje 'to_dict' metoduna sahipse güvenle geçilebilir
    if hasattr(fig_obj, "to_dict"):
        return fig_obj
    # Matplotlib Figure ise
    return mpl_image_figure(fig_to_base64(fig_obj))

# =============================================================================
# CombinedPolicy kurulumunu güvenli yapan yardımcı
# =============================================================================
def _supported_params() -> set:
    """CombinedPolicy'nin desteklediği __init__ parametrelerini oku."""
    return set(inspect.signature(CombinedPolicy).parameters.keys())

def build_combo(
    heat_backend: str,
    umap_backend: str,
    normalizer: str,
    clusterer: str,
    cmap_scale: str,
    pre_umap_pca: Optional[int],
) -> CombinedPolicy:
    """
    CombinedPolicy'yi kurar. CombinedPolicy imzasında bulunmayan anahtarları
    otomatik olarak filtreler (sürüm farklarına dayanıklı).
    """
    sup = _supported_params()

    # UMAP scatter backend tercih sırası
    if umap_backend in REG.scatter_plotters:
        umap_scatter = umap_backend
    elif "plotly_express" in REG.scatter_plotters:
        umap_scatter = "plotly_express"
    elif "plotly" in REG.scatter_plotters:
        umap_scatter = "plotly"
    else:
        umap_scatter = "matplotlib"

    # Iskelet konfig
    cfg: Dict[str, Any] = {
        # ortak
        "normalizer": normalizer,
        # heatmap
        "heat_clusterer": clusterer,
        "heat_cmap_scale": cmap_scale,
        "heat_plotter": heat_backend,
        "heat_plot_kwargs": dict(show=False),
        # umap
        "umap_neighbors": "sklearn",
        "umap_reducer": "umap",
        "umap_scatter": umap_scatter,
        "umap_neighbors_kwargs": dict(n_neighbors=15, metric="euclidean"),
        "umap_reducer_kwargs": dict(n_components=2, n_neighbors=15, min_dist=0.1,
                                    metric="euclidean", random_state=7),
        "umap_scatter_kwargs": dict(show=False),
        # opsiyonel: bazı sürümlerde yok
        "pre_umap_pca": pre_umap_pca,
    }

    # Sadece CombinedPolicy'nin bildiği anahtarları geçir
    safe_cfg = {k: v for k, v in cfg.items() if (k in sup and v is not None)}
    return CombinedPolicy(**safe_cfg)

# =============================================================================
# Demo veri (yerine kendi verini koyabilirsin)
# =============================================================================
rng = np.random.default_rng(42)
X = rng.normal(size=(40, 25)) * 2 + rng.normal(size=(40, 1)) * 3
rows = [f"Gene_{i}" for i in range(X.shape[0])]
cols = [f"Cell_{j}" for j in range(X.shape[1])]
point_labels = [f"C{i}" for i in range(X.shape[0])]

# =============================================================================
# Dash Uygulaması
# =============================================================================
app = Dash(__name__)
app.layout = html.Div([
    html.Div([
        html.Label("Heatmap backend"),
        dcc.Dropdown(["matplotlib", "plotly"], "plotly", id="heat-backend"),

        html.Label("UMAP backend"),
        dcc.Dropdown(["plotly_express", "matplotlib"], "plotly_express", id="umap-backend"),

        html.Label("Normalizer"),
        dcc.Dropdown(list(REG.normalizers.keys()), "zscore_rows", id="normalizer"),

        html.Label("Clustering"),
        dcc.Dropdown(list(REG.clusterers.keys()), "simple", id="clusterer"),

        html.Label("Colormap scale"),
        dcc.Dropdown(list(REG.cmap_scales.keys()), "center_zero", id="cmap_scale"),

        html.Label("Pre-UMAP PCA components"),
        dcc.Slider(min=0, max=X.shape[1], step=1, value=min(10, X.shape[1]),
                   marks={i: str(i) for i in range(0, X.shape[1] + 1)}, id="pca-components"),
    ], style={"width": "30%", "display": "inline-block", "verticalAlign": "top", "padding": "10px"}),

    html.Div([
        html.H4("Heatmap"),
        dcc.Graph(id="heatmap-graph"),

        html.H4("UMAP"),
        dcc.Graph(id="umap-graph"),
    ], style={"width": "65%", "display": "inline-block", "paddingLeft": "20px"})
])

# =============================================================================
# Callback
# =============================================================================
@app.callback(
    Output("heatmap-graph", "figure"),
    Output("umap-graph", "figure"),
    Input("heat-backend", "value"),
    Input("umap-backend", "value"),
    Input("normalizer", "value"),
    Input("clusterer", "value"),
    Input("cmap_scale", "value"),
    Input("pca-components", "value"),
)
def update_plots(heat_backend, umap_backend, normalizer, clusterer, cmap_scale, pca_components):
    try:
        pre_pca = int(pca_components) if (pca_components is not None and int(pca_components) > 0) else None

        combo = build_combo(
            heat_backend=heat_backend,
            umap_backend=umap_backend,
            normalizer=normalizer,
            clusterer=clusterer,
            cmap_scale=cmap_scale,
            pre_umap_pca=pre_pca,
        )

        # Pipeline’i çalıştır
        (heat_fig, umap_fig), meta = combo.apply(
            X, row_labels=rows, col_labels=cols, point_labels=point_labels,
            heat_title="Heatmap", umap_title="UMAP"
        )

        # Dash'e uygun figür nesneleri
        heat_out = to_dash_figure(heat_fig)
        umap_out = to_dash_figure(umap_fig)

        return heat_out, umap_out

    except Exception as e:
        # Hata durumunda ekrana yaz ve boş figür dön
        msg = f"Callback error: {type(e).__name__}: {e}"
        print(msg)
        fallback = {"data": [], "layout": {"title": msg}}
        return fallback, fallback

# =============================================================================
# Main
# =============================================================================
if __name__ == "__main__":
    # İstersen tarayıcıyı otomatik açabilirsin:
    # import webbrowser; webbrowser.open("http://127.0.0.1:8050/")
    app.run(debug=True)
