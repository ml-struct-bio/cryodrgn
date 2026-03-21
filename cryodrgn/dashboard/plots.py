"""Plotly figures for the dashboard."""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from cryodrgn.dashboard.data import DashboardExperiment

# Pastel blue for upper-triangle scatters and hex density.
_UPPER_SCATTER_BLUE = "rgba(154, 196, 232, 0.42)"
_UPPER_HEX_COLORSCALE: list[list] = [
    [0.0, "rgb(248, 252, 255)"],
    [0.3, "rgb(220, 236, 250)"],
    [0.55, "rgb(190, 220, 245)"],
    [0.8, "rgb(150, 200, 235)"],
    [1.0, "rgb(115, 178, 228)"],
]


def scatter_json(
    exp: DashboardExperiment,
    xcol: str,
    ycol: str,
    color_col: Optional[str],
    max_points: int = 120_000,
) -> str:
    df = exp.plot_df
    n = len(df)
    rng = np.random.default_rng(0)
    if n > max_points:
        pick = rng.choice(n, size=max_points, replace=False)
        sub = df.iloc[pick].reset_index(drop=True)
        row_map = pick.tolist()
    else:
        sub = df
        row_map = list(range(n))

    if color_col and color_col != "none" and color_col in sub.columns:
        marker = dict(size=4, opacity=0.35, color=sub[color_col], colorscale="Viridis")
    else:
        marker = dict(size=4, opacity=0.35, color="#4a5568")

    text = [
        f"row {r} · index {int(sub.iloc[i]['index'])}" for i, r in enumerate(row_map)
    ]

    sc = go.Scattergl(
        x=sub[xcol],
        y=sub[ycol],
        mode="markers",
        text=text,
        customdata=np.column_stack([sub["index"].values, row_map]),
        hovertemplate="%{text}<extra></extra>",
        marker=marker,
    )
    fig = go.Figure(sc)
    fig.update_layout(
        template="plotly_white",
        margin=dict(l=50, r=20, t=40, b=50),
        title=f"{ycol} vs {xcol}",
        xaxis_title=xcol,
        yaxis_title=ycol,
        dragmode="lasso",
        uirevision="scatter",
    )
    return fig.to_json()


def pair_grid_json(
    exp: DashboardExperiment,
    lower_color_col: str,
    diagonal_emb: str,
    upper_style: str,
) -> str:
    """z_dim × z_dim grid: axes are z_j (x) vs z_i (y); lower triangle shares one colour covariate."""
    df = exp.plot_df
    zdim = int(exp.z.shape[1])
    zcols = [f"z{i}" for i in range(zdim)]
    for c in zcols:
        if c not in df.columns:
            raise ValueError(f"Missing latent column {c} in analysis table.")

    if lower_color_col not in df.columns:
        raise ValueError(f"Unknown colour column: {lower_color_col}")
    color_s = df[lower_color_col]
    if not pd.api.types.is_numeric_dtype(color_s):
        color_s = pd.to_numeric(color_s, errors="coerce")
    if color_s.isna().all():
        raise ValueError(
            f"Lower-triangle colour column `{lower_color_col}` has no numeric values."
        )
    cvals = np.asarray(color_s, dtype=np.float64)
    cfinite = cvals[np.isfinite(cvals)]
    if cfinite.size == 0:
        cmin, cmax = 0.0, 1.0
    elif np.isclose(cfinite.min(), cfinite.max()):
        cmin = float(cfinite.min()) - 0.5
        cmax = float(cfinite.max()) + 0.5
    else:
        cmin = float(np.min(cfinite))
        cmax = float(np.max(cfinite))

    c_mid = 0.5 * (cmin + cmax)
    cvals_plot = np.where(np.isfinite(cvals), cvals, c_mid)

    emb = (diagonal_emb or "pc").lower()
    if emb == "umap":
        if exp.umap is None or "UMAP1" not in df.columns or "UMAP2" not in df.columns:
            raise ValueError("UMAP embedding is not available for this run.")
        emb_x, emb_y = df["UMAP1"], df["UMAP2"]
    elif emb == "pc":
        if "PC1" not in df.columns or "PC2" not in df.columns:
            raise ValueError("PCA components PC1/PC2 are not available.")
        emb_x, emb_y = df["PC1"], df["PC2"]
    else:
        raise ValueError("diagonal_emb must be 'pc' or 'umap'.")

    upper = (upper_style or "scatter").lower()
    if upper not in ("scatter", "hex"):
        raise ValueError("upper_style must be 'scatter' or 'hex'.")

    n_pts = len(df)
    hex_bins = int(np.clip(25 + np.sqrt(max(n_pts, 1)), 35, 80))

    fig = make_subplots(
        rows=zdim,
        cols=zdim,
        horizontal_spacing=0.04,
        vertical_spacing=0.04,
    )

    lower_colorbar_done = False
    for i in range(zdim):
        for j in range(zdim):
            row_, col_ = i + 1, j + 1
            xi = df[zcols[j]]
            yi = df[zcols[i]]
            if i == j:
                fig.add_trace(
                    go.Scattergl(
                        x=emb_x,
                        y=emb_y,
                        mode="markers",
                        marker=dict(
                            size=2,
                            color=df[zcols[i]],
                            colorscale="Viridis",
                            opacity=0.5,
                            showscale=False,
                        ),
                        showlegend=False,
                    ),
                    row=row_,
                    col=col_,
                )
            elif i < j:
                if upper == "hex":
                    fig.add_trace(
                        go.Histogram2d(
                            x=xi,
                            y=yi,
                            colorscale=_UPPER_HEX_COLORSCALE,
                            showscale=False,
                            nbinsx=hex_bins,
                            nbinsy=hex_bins,
                        ),
                        row=row_,
                        col=col_,
                    )
                else:
                    fig.add_trace(
                        go.Scattergl(
                            x=xi,
                            y=yi,
                            mode="markers",
                            marker=dict(size=2, color=_UPPER_SCATTER_BLUE),
                            showlegend=False,
                        ),
                        row=row_,
                        col=col_,
                    )
            else:
                show_cb = not lower_colorbar_done
                lower_colorbar_done = True
                mk_lower: dict = dict(
                    size=3,
                    color=cvals_plot,
                    cmin=cmin,
                    cmax=cmax,
                    colorscale="Viridis",
                    opacity=0.55,
                    showscale=show_cb,
                )
                if show_cb:
                    mk_lower["colorbar"] = dict(
                        orientation="h",
                        x=0.5,
                        xanchor="center",
                        y=-0.04,
                        yanchor="top",
                        len=0.48,
                        thickness=14,
                        outlinewidth=0,
                        title=dict(
                            text=f"Lower ▽ colour · {lower_color_col}",
                            side="bottom",
                            font=dict(size=11, color="#243b53"),
                        ),
                        tickfont=dict(size=10, color="#243b53"),
                    )
                fig.add_trace(
                    go.Scattergl(
                        x=xi,
                        y=yi,
                        mode="markers",
                        marker=mk_lower,
                        showlegend=False,
                    ),
                    row=row_,
                    col=col_,
                )

    fig.update_xaxes(
        showticklabels=False,
        showline=True,
        linewidth=1.35,
        linecolor="rgba(36, 59, 83, 0.55)",
        mirror=True,
        zeroline=False,
        showgrid=False,
    )
    fig.update_yaxes(
        showticklabels=False,
        showline=True,
        linewidth=1.35,
        linecolor="rgba(36, 59, 83, 0.55)",
        mirror=True,
        zeroline=False,
        showgrid=False,
    )

    diag_lbl = "UMAP" if emb == "umap" else "PC1–PC2"
    upper_lbl = "hex density" if upper == "hex" else "scatter"
    legend_html = (
        f"<b>Regions</b><br>"
        f"• Upper △ — pastel blue ({upper_lbl})<br>"
        f"• Diagonal — {diag_lbl}, marker colour = zᵢ<br>"
        f"• Lower ▽ — axes zⱼ vs zᵢ, marker colour = {lower_color_col}"
    )
    fig.update_layout(
        template="plotly_white",
        margin=dict(l=36, r=28, t=52, b=132),
        title=(
            f"z_dim × z_dim · diagonal {diag_lbl} · upper pastel {upper_lbl} · "
            f"lower by {lower_color_col}"
        ),
        autosize=True,
        uirevision="pair",
        annotations=[
            dict(
                text=legend_html,
                xref="paper",
                yref="paper",
                x=0.5,
                y=-0.19,
                xanchor="center",
                yanchor="top",
                showarrow=False,
                align="left",
                font=dict(size=11, color="#243b53", family="Arial, sans-serif"),
                bgcolor="rgba(250, 248, 244, 0.96)",
                bordercolor="rgba(36, 59, 83, 0.28)",
                borderwidth=1,
                borderpad=8,
            )
        ],
    )

    return fig.to_json()
