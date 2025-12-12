import matplotlib
matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import io, base64

from dash import Dash, html, dcc, Input, Output
from viz_flexible_interactive import CombinedPolicy, REG

# -------------------------
# Helper functions
# -------------------------
def fig_to_base64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return "data:image/png;base64," + base64.b64encode(buf.read()).decode("utf-8")

def plot_heat_mpl_dash(X, row_labels=None, col_labels=None):
    fig, ax = plt.subplots(figsize=(8,6))
    im = ax.imshow(X, aspect="auto", interpolation="nearest", cmap="viridis")
    if row_labels:
        ax.set_yticks(range(len(row_labels)))
        ax.set_yticklabels(row_labels, fontsize=8)
    if col_labels:
        ax.set_xticks(range(len(col_labels)))
        ax.set_xticklabels(col_labels, fontsize=8, rotation=90)
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    return fig_to_base64(fig)

def plot_umap_mpl_dash(emb, labels=None):
    fig, ax = plt.subplots(figsize=(8,6))
    if labels is None:
        labels = np.arange(len(emb))
    scatter = ax.scatter(emb[:,0], emb[:,1], c=labels, cmap="tab20")
    fig.colorbar(scatter, ax=ax, label="Label")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("UMAP")
    fig.tight_layout()
    return fig_to_base64(fig)

def mpl_to_dash_figure(base64_img):
    return {
        "data": [],
        "layout": {
            "images": [{
                "source": base64_img,
                "xref": "paper",
                "yref": "paper",
                "x": 0, "y": 1,
                "sizex": 1, "sizey": 1,
                "xanchor": "left", "yanchor": "top",
                "layer": "below"
            }],
            "xaxis": {"visible": False},
            "yaxis": {"visible": False},
            "margin": {"l":0, "r":0, "t":0, "b":0},
        }
    }

# -------------------------
# Dummy data
# -------------------------
rng = np.random.default_rng(42)
X = rng.normal(size=(40, 25)) * 2 + rng.normal(size=(40, 1)) * 3
rows = [f"Gene_{i}" for i in range(X.shape[0])]
cols = [f"Cell_{j}" for j in range(X.shape[1])]
point_labels = [f"C{i}" for i in range(X.shape[0])]

# -------------------------
# Dash app
# -------------------------
app = Dash(__name__)
app.layout = html.Div([
    html.Div([
        html.Label("Heatmap backend"),
        dcc.Dropdown(["matplotlib", "plotly"], "plotly", id="heat-backend"),
        html.Label("UMAP backend"),
        dcc.Dropdown(["matplotlib", "plotly_express"], "plotly_express", id="umap-backend"),
        html.Label("Normalizer"),
        dcc.Dropdown(list(REG.normalizers.keys()), "zscore_rows", id="normalizer"),
        html.Label("Clustering"),
        dcc.Dropdown(list(REG.clusterers.keys()), "simple", id="clusterer"),
        html.Label("Colormap scale"),
        dcc.Dropdown(list(REG.cmap_scales.keys()), "center_zero", id="cmap_scale"),
        html.Label("Pre-UMAP PCA components"),
        dcc.Slider(min=0, max=25, step=1, value=10,
                   marks={i:str(i) for i in range(0,26)}, id="pca-components")
    ], style={"width": "30%", "display": "inline-block", "verticalAlign": "top", "padding": "10px"}),

    html.Div([
        html.H4("Heatmap"),
        dcc.Graph(id="heatmap-graph"),
        html.H4("UMAP"),
        dcc.Graph(id="umap-graph"),
    ], style={"width": "65%", "display": "inline-block", "paddingLeft": "20px"})
])

# -------------------------
# Robust callback
# -------------------------
@app.callback(
    Output("heatmap-graph", "figure"),
    Output("umap-graph", "figure"),
    Input("heat-backend", "value"),
    Input("umap-backend", "value"),
    Input("normalizer", "value"),
    Input("clusterer", "value"),
    Input("cmap_scale", "value"),
    Input("pca-components", "value")
)
def update_plots(heat_backend, umap_backend, normalizer, clusterer, cmap_scale, pca_components):
    try:
        # Always use plotly_express internally for UMAP in CombinedPolicy
        combo = CombinedPolicy(
            normalizer=normalizer,
            heat_clusterer=clusterer,
            heat_cmap_scale=cmap_scale,
            heat_plotter=heat_backend,
            heat_plot_kwargs=dict(show=False),
            pre_umap_pca=pca_components if pca_components > 0 else None,
            umap_scatter="plotly_express",
            umap_scatter_kwargs=dict(show=False),
        )

        (heat_fig, umap_fig), meta = combo.apply(X, row_labels=rows, col_labels=cols, point_labels=point_labels)
        
        # --- Heatmap backend ---
        if heat_backend == "matplotlib":
            heat_fig = mpl_to_dash_figure(plot_heat_mpl_dash(meta["X_clustered"], rows, cols))
        
        # --- UMAP backend ---
        if umap_backend == "matplotlib":
            emb = meta.get("embedding")
            if emb is None or emb.shape[1] < 2:
                print("UMAP embedding invalid:", emb)
                umap_fig = {"data": [], "layout": {"title":"UMAP embedding not available"}}
            else:
                umap_fig = mpl_to_dash_figure(plot_umap_mpl_dash(emb, point_labels))
        # else keep Plotly figure as is
        
        return heat_fig, umap_fig

    except Exception as e:
        print("Callback error:", e)
        fallback_fig = {"data": [], "layout": {"title": "Error in figure"}}
        return fallback_fig, fallback_fig

# -------------------------
if __name__ == "__main__":
    import webbrowser
    webbrowser.open("http://127.0.0.1:8050/")
    app.run(debug=True)
