"""Plotly figures for the dashboard; pair grid uses Matplotlib + Seaborn (PNG)."""

from __future__ import annotations

import io
from collections.abc import Collection
from typing import Any, cast

import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import seaborn as sns
from matplotlib.offsetbox import (
    AnchoredOffsetbox,
    DrawingArea,
    HPacker,
    TextArea,
)
from matplotlib.patches import Rectangle
from seaborn.distributions import (
    _freedman_diaconis_bins,  # type: ignore[attr-defined]
)
from seaborn.palettes import blend_palette
from seaborn.utils import set_hls_values

from cryodrgn import analysis
from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.mpl_style import ezlab_matplotlib_rc


def _lower_color_series_is_discrete(s: pd.Series) -> bool:
    """Treat integer / categorical / few whole-valued floats as discrete colouring."""
    s = s.dropna()
    if len(s) == 0:
        return False
    if pd.api.types.is_bool_dtype(s):
        return True
    if isinstance(s.dtype, pd.CategoricalDtype):
        return True
    if pd.api.types.is_integer_dtype(s) and not pd.api.types.is_float_dtype(s):
        return True
    if s.dtype == object:
        return True
    if pd.api.types.is_float_dtype(s):
        u = np.unique(s.to_numpy())
        if len(u) > 40:
            return False
        if len(u) == 0:
            return False
        return bool(
            np.all(np.isfinite(u)) and np.allclose(u, np.round(u), rtol=0, atol=1e-8)
        )
    return False


# Pastel blue for upper-triangle scatter (matplotlib).
_UPPER_SCATTER_BLUE = "#9ec5e8"
# Matches ``base.html`` ``--cream`` (page background).
_DASHBOARD_CREAM = "#faf8f4"

# Pair grid: gutter between cells (fraction of average subplot width/height).
_PAIR_CELL_WSPACE = 0.06
_PAIR_CELL_HSPACE = 0.06

# Figure margins for the pair subplot block (before ``_refine_pair_grid_right_margin``).
_PAIR_GRID_LEFT = 0.035
_PAIR_GRID_TOP = 0.98
_PAIR_GRID_BOTTOM = 0.06
_PAIR_GRID_RIGHT_INITIAL = 0.68
_PAIR_GRID_RIGHT_PLACEHOLDER = 0.712

# Shared Plotly font (matches base.html --font-sans).
_PLOTLY_FONT = dict(
    family="Barlow, Roboto, -apple-system, BlinkMacSystemFont, sans-serif",
    size=13,
)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _subsample(
    df: pd.DataFrame,
    max_points: int | None,
    *,
    seed: int = 0,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Random subsample of *df*; returns ``(sub_df, original_row_indices)``."""
    n = len(df)
    if max_points is not None and n > max_points:
        pick = np.random.default_rng(seed).choice(
            n,
            size=max_points,
            replace=False,
        )
        return (
            df.iloc[pick].reset_index(drop=True),
            pick.astype(np.int64, copy=False),
        )
    return df, np.arange(n, dtype=np.int64)


def _axes_cell_bboxes(axes: np.ndarray) -> list[dict[str, float]]:
    """Per-cell axis bounding boxes in figure coordinates (row-major)."""
    zdim = int(axes.shape[0])
    out: list[dict[str, float]] = []
    for i in range(zdim):
        for j in range(zdim):
            bb = axes[i, j].get_position()
            out.append(
                {
                    "x0": float(bb.x0),
                    "y0": float(bb.y0),
                    "x1": float(bb.x1),
                    "y1": float(bb.y1),
                }
            )
    return out


def _plotly_to_json(fig: go.Figure) -> str:
    sj = fig.to_json()
    if sj is None:
        raise RuntimeError("Plotly failed to serialize figure.")
    return sj


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def pair_grid_skeleton_placeholder_layout(
    zdim: int,
) -> list[dict[str, float]]:
    """Cell bboxes for the HTML skeleton before the first PNG loads."""
    if zdim < 1:
        return []
    with ezlab_matplotlib_rc():
        fig, axes = plt.subplots(
            zdim,
            zdim,
            figsize=(4, 4),
            squeeze=False,
        )
        fig.subplots_adjust(
            left=_PAIR_GRID_LEFT,
            right=_PAIR_GRID_RIGHT_PLACEHOLDER,
            top=_PAIR_GRID_TOP,
            bottom=_PAIR_GRID_BOTTOM,
            wspace=_PAIR_CELL_WSPACE,
            hspace=_PAIR_CELL_HSPACE,
        )
        out = _axes_cell_bboxes(np.asarray(axes))
        plt.close(fig)
        return out


def _pair_jointplot_hex_cmap(color: str | None = None):
    """Same sequential colormap as ``sns.jointplot(..., kind=\"hex\")`` (default ``color``)."""
    if color is None:
        color = "C0"
    color_rgb = mcolors.to_rgb(color)
    ramp = [set_hls_values(color_rgb, l=val) for val in np.linspace(1, 0, 12)]
    return blend_palette(ramp, as_cmap=True)


def _pair_jointplot_hex_gridsize(x: np.ndarray, y: np.ndarray) -> int:
    """Match ``sns.jointplot`` hex ``gridsize`` (Freedman–Diaconis, cap 50)."""
    xb = min(_freedman_diaconis_bins(x), 50)
    yb = min(_freedman_diaconis_bins(y), 50)
    return max(int(np.mean([xb, yb])), 1)


def _refine_pair_grid_right_margin(
    fig: plt.Figure,
    *,
    grid_side_in: float,
    left_m: float,
    top_m: float,
    bottom_m: float,
    right_artist,
    max_gap_inches: float,
) -> None:
    """Set subplot ``right`` so gap to legend / colour bar is at most ``max_gap_inches`` wide."""
    h_frac = top_m - bottom_m
    fig_h = grid_side_in / h_frac
    for _ in range(2):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        if hasattr(right_artist, "get_tightbbox"):
            bb = right_artist.get_tightbbox(renderer)
        else:
            bb = right_artist.get_window_extent(renderer)
        bb_fig = bb.transformed(fig.transFigure.inverted())
        x0 = float(bb_fig.x0)
        fig_w_in = float(fig.get_figwidth())
        max_gap_frac = max_gap_inches / fig_w_in
        eps = 2.0 / (fig.dpi * fig_w_in)
        right_m2 = x0 - max(max_gap_frac, eps)
        right_m2 = max(right_m2, left_m + 0.12)
        right_m2 = min(right_m2, x0 - eps)
        w_sub = right_m2 - left_m
        if w_sub <= 0.02:
            break
        fig_w_new = grid_side_in / w_sub
        fig.set_size_inches(fig_w_new, fig_h)
        fig.subplots_adjust(
            left=left_m,
            right=right_m2,
            top=top_m,
            bottom=bottom_m,
            wspace=_PAIR_CELL_WSPACE,
            hspace=_PAIR_CELL_HSPACE,
        )


def _pair_lower_triangle_pictogram_da(px: float = 28.0) -> DrawingArea:
    """4×4 lower-triangle mini-grid (matches dashboard ``picto_lower`` SVG)."""
    gap = 0.75
    inner = max(px - 3 * gap, 1.0)
    cs = inner / 4.0
    da = DrawingArea(px, px, clip=False)
    lo = "#e8ecf1"
    hi = "#3d5a80"
    # Rows bottom → top in figure space (matches SVG rows bottom → top).
    fills_from_bottom = (
        (hi, hi, hi, lo),
        (hi, hi, lo, lo),
        (hi, lo, lo, lo),
        (lo, lo, lo, lo),
    )
    for row in range(4):
        for col in range(4):
            x = col * (cs + gap)
            y = row * (cs + gap)
            da.add_artist(
                Rectangle(
                    (x, y),
                    cs,
                    cs,
                    facecolor=fills_from_bottom[row][col],
                    edgecolor="#8899aa",
                    linewidth=0.35,
                )
            )
    return da


def _lower_legend_covariate_title(
    lower_color_col: str, exp: DashboardExperiment
) -> str:
    if lower_color_col == "labels":
        return f"kmeans{exp.kmeans_folder_id} labels"
    return lower_color_col


def _lower_legend_entry_label(lower_color_col: str, u: Any) -> str:
    """1-based cluster indices for the ``labels`` (k-means) column."""
    if pd.isna(u):
        return str(u)
    if lower_color_col == "labels":
        if isinstance(u, (bool, np.bool_)):
            return str(u)
        if isinstance(u, (int, np.integer)):
            return str(int(u) + 1)
        if isinstance(u, float) and np.isfinite(u) and u == int(u):
            return str(int(u) + 1)
        try:
            fu = float(u)
            if np.isfinite(fu) and fu == int(fu):
                return str(int(fu) + 1)
        except (TypeError, ValueError):
            pass
    return str(u)


def _attach_pair_lower_legend_caption(
    fig: plt.Figure,
    leg: Any,
    covariate_title: str,
    *,
    title_fontsize: float = 24.0,
) -> None:
    """Legend title row: lower-triangle pictogram + covariate name (no ▽ glyph)."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bb_disp = leg.get_window_extent(renderer)
    bb_fig = bb_disp.transformed(fig.transFigure.inverted())
    pict = _pair_lower_triangle_pictogram_da(26.0)
    text_area = TextArea(
        f"  {covariate_title}",
        textprops=dict(
            fontsize=title_fontsize,
            fontweight="bold",
            color="#243b53",
        ),
    )
    pack = HPacker(children=[pict, text_area], align="center", pad=0, sep=2)
    y_top = min(0.998, float(bb_fig.y1) + 0.014)
    ab = AnchoredOffsetbox(
        "lower left",
        child=pack,
        frameon=False,
        bbox_to_anchor=(float(bb_fig.x0), y_top),
        bbox_transform=fig.transFigure,
        borderpad=0,
    )
    fig.add_artist(ab)


def scatter_json(
    exp: DashboardExperiment,
    xcol: str,
    ycol: str,
    color_col: str | None,
    max_points: int | None = 120_000,
    *,
    preselect_plot_df_rows: Collection[int] | None = None,
    use_webgl: bool = True,
) -> str:
    sub, row_indices = _subsample(exp.plot_df, max_points, seed=0)

    if color_col and color_col != "none" and color_col in sub.columns:
        marker = dict(
            size=4,
            opacity=0.35,
            color=sub[color_col],
            colorscale="Viridis",
        )
    else:
        marker = dict(size=4, opacity=0.35, color="#4a5568")

    customdata = np.column_stack(
        [sub["index"].to_numpy(), row_indices],
    )
    # Scattergl can leave Plotly.react() pending on some browsers/GPUs.
    trace_cls = go.Scattergl if use_webgl else go.Scatter
    sc = trace_cls(
        x=sub[xcol],
        y=sub[ycol],
        mode="markers",
        customdata=customdata,
        hovertemplate=(
            "row %{customdata[1]} · index %{customdata[0]}" "<extra></extra>"
        ),
        marker=marker,
    )

    layout_kw: dict[str, Any] = dict(
        template="plotly_white",
        paper_bgcolor=_DASHBOARD_CREAM,
        margin=dict(l=50, r=20, t=40, b=50),
        title=f"{ycol} vs {xcol}",
        xaxis_title=xcol,
        yaxis_title=ycol,
        dragmode="lasso",
        uirevision="scatter",
        font=_PLOTLY_FONT,
    )
    if preselect_plot_df_rows is not None:
        want = frozenset(int(x) for x in preselect_plot_df_rows)
        layout_kw["meta"] = dict(
            cdrgn_preselected=[
                int(i) for i in range(len(row_indices)) if int(row_indices[i]) in want
            ]
        )

    fig = go.Figure(sc)
    fig.update_layout(**layout_kw)
    return _plotly_to_json(fig)


def scatter3d_z_json(
    exp: DashboardExperiment,
    xcol: str,
    ycol: str,
    zcol: str,
    color_col: str | None,
    max_points: int = 120_000,
) -> str:
    """Interactive 3D scatter of three latent ``z*`` columns."""
    df = exp.plot_df
    z_allowed = {f"z{i}" for i in range(int(exp.z.shape[1]))}
    for c in (xcol, ycol, zcol):
        if c not in z_allowed:
            raise ValueError(f"Axis {c!r} is not a latent dimension for this run.")
        if c not in df.columns:
            raise ValueError(f"Missing column {c!r} in analysis table.")
    if len({xcol, ycol, zcol}) < 3:
        raise ValueError("Choose three distinct latent axes.")

    sub, row_indices = _subsample(df, max_points, seed=1)

    if color_col and color_col != "none" and color_col in sub.columns:
        marker = dict(
            size=0.75,
            opacity=0.4,
            color=sub[color_col],
            colorscale="Viridis",
        )
    else:
        marker = dict(size=0.75, opacity=0.4, color="#4a5568")

    customdata = np.column_stack(
        [sub["index"].to_numpy(), row_indices],
    )
    sc = go.Scatter3d(
        x=sub[xcol],
        y=sub[ycol],
        z=sub[zcol],
        mode="markers",
        customdata=customdata,
        hovertemplate=(
            "row %{customdata[1]} · index %{customdata[0]}" "<extra></extra>"
        ),
        marker=marker,
    )
    transparent = "rgba(0,0,0,0)"
    fig = go.Figure(sc)
    fig.update_layout(
        template="plotly_white",
        autosize=True,
        paper_bgcolor=_DASHBOARD_CREAM,
        margin=dict(l=0, r=0, t=44, b=0),
        title=f"{zcol} · {ycol} · {xcol}",
        scene=dict(
            xaxis_title=xcol,
            yaxis_title=ycol,
            zaxis_title=zcol,
            aspectmode="data",
            bgcolor="rgba(250,248,244,0.95)",
            xaxis=dict(backgroundcolor=transparent),
            yaxis=dict(backgroundcolor=transparent),
            zaxis=dict(backgroundcolor=transparent),
        ),
        uirevision="scatter3d_z",
        font=_PLOTLY_FONT,
    )
    return _plotly_to_json(fig)


def pair_grid_png(
    exp: DashboardExperiment,
    lower_color_col: str,
    diagonal_emb: str,
    upper_style: str,
    dpi: int = 120,
) -> tuple[bytes, list[dict[str, float]]]:
    """z_dim × z_dim Matplotlib grid (square cells); upper hex matches ``sns.jointplot(..., kind=\"hex\")``.

    Returns PNG bytes and per-cell axis bboxes in figure coordinates (for HTML overlay alignment).
    """
    df = exp.plot_df
    zdim = int(exp.z.shape[1])
    zcols = [f"z{i}" for i in range(zdim)]
    for c in zcols:
        if c not in df.columns:
            raise ValueError(f"Missing latent column {c} in analysis table.")

    if lower_color_col not in df.columns:
        raise ValueError(f"Unknown colour column: {lower_color_col}")
    raw_color = df[lower_color_col]
    discrete = _lower_color_series_is_discrete(raw_color)

    pal: list[str] = []
    uniques: np.ndarray | pd.Index = np.array([])
    codes: np.ndarray | None = None
    point_colors: list[str] | None = None
    cvals_plot: np.ndarray | None = None
    cmin: float | None = None
    cmax: float | None = None

    if discrete:
        codes_arr, uniques = pd.factorize(raw_color, sort=True)
        codes = np.asarray(codes_arr, dtype=np.int64)
        if len(uniques) == 0:
            raise ValueError(
                f"Lower-triangle colour column `{lower_color_col}` has no values."
            )
        n_u = len(uniques)
        pal = list(analysis._get_chimerax_colors(max(n_u, 1)))
        point_colors = []
        for k in range(len(df)):
            code = int(codes[k])
            if code < 0:
                point_colors.append("#aab4bf")
            else:
                point_colors.append(pal[code % len(pal)])
    else:
        color_num = cast(pd.Series, pd.to_numeric(raw_color, errors="coerce"))
        if bool(color_num.isna().all()):
            raise ValueError(
                f"Lower-triangle colour column `{lower_color_col}` has no numeric values."
            )
        cvals = np.asarray(color_num, dtype=np.float64)
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
        emb_x = df["UMAP1"].to_numpy(dtype=np.float64)
        emb_y = df["UMAP2"].to_numpy(dtype=np.float64)
    elif emb == "pc":
        if "PC1" not in df.columns or "PC2" not in df.columns:
            raise ValueError("PCA components PC1/PC2 are not available.")
        emb_x = df["PC1"].to_numpy(dtype=np.float64)
        emb_y = df["PC2"].to_numpy(dtype=np.float64)
    else:
        raise ValueError("diagonal_emb must be 'pc' or 'umap'.")

    upper = (upper_style or "scatter").lower()
    if upper not in ("scatter", "hex"):
        raise ValueError("upper_style must be 'scatter' or 'hex'.")

    n_pts = len(df)

    with ezlab_matplotlib_rc():
        hex_cmap = _pair_jointplot_hex_cmap()

        inch_per = max(2.35, min(2.85, 11.0 / max(zdim, 1)))
        left_m, top_m = _PAIR_GRID_LEFT, _PAIR_GRID_TOP
        bottom_m = _PAIR_GRID_BOTTOM
        # Subplot grid right edge (legend / vertical colour bar sit to the right).
        right_axes = _PAIR_GRID_RIGHT_INITIAL
        legend_handles: list[mlines.Line2D] | None = None
        if discrete:
            legend_handles = []
            for idx, u in enumerate(uniques):
                if pd.isna(u):
                    continue
                legend_handles.append(
                    mlines.Line2D(
                        [0],
                        [0],
                        marker="o",
                        color="w",
                        label=_lower_legend_entry_label(lower_color_col, u),
                        markerfacecolor=pal[idx % len(pal)],
                        markersize=16.5,
                        linestyle="None",
                    )
                )
            if bool(raw_color.isna().any()):
                legend_handles.append(
                    mlines.Line2D(
                        [0],
                        [0],
                        marker="o",
                        color="w",
                        label="(missing)",
                        markerfacecolor="#aab4bf",
                        markersize=16.5,
                        linestyle="None",
                    )
                )

        # Square axes region; reserve a strip on the right for legend or colour bar.
        w_frac = right_axes - left_m
        h_frac = top_m - bottom_m
        grid_side_in = inch_per * zdim
        fig_w = grid_side_in / w_frac
        fig_h = grid_side_in / h_frac

        with sns.axes_style(
            "white",
            rc={
                "axes.facecolor": _DASHBOARD_CREAM,
                "figure.facecolor": _DASHBOARD_CREAM,
            },
        ):
            fig, axes = plt.subplots(
                zdim,
                zdim,
                figsize=(fig_w, fig_h),
                squeeze=False,
                constrained_layout=False,
            )
        fig.patch.set_facecolor(_DASHBOARD_CREAM)

        lower_cbar_mappable = None

        for i in range(zdim):
            for j in range(zdim):
                ax = axes[i, j]
                xi = df[zcols[j]].to_numpy(dtype=np.float64)
                yi = df[zcols[i]].to_numpy(dtype=np.float64)
                if i == j:
                    zi = df[zcols[i]].to_numpy(dtype=np.float64)
                    ax.scatter(
                        emb_x,
                        emb_y,
                        c=zi,
                        cmap="viridis",
                        s=2,
                        alpha=0.55,
                        linewidths=0,
                        rasterized=True,
                    )
                elif i < j:
                    if upper == "hex":
                        mxy = np.isfinite(xi) & np.isfinite(yi)
                        xh, yh = xi[mxy], yi[mxy]
                        if xh.size == 0:
                            pass
                        else:
                            try:
                                gs = _pair_jointplot_hex_gridsize(xh, yh)
                                ax.hexbin(xh, yh, gridsize=gs, cmap=hex_cmap)
                            except ZeroDivisionError:
                                ax.hexbin(
                                    xh,
                                    yh,
                                    gridsize=int(
                                        np.clip(22 + np.sqrt(max(n_pts, 1)), 30, 75)
                                    ),
                                    cmap=hex_cmap,
                                    mincnt=1,
                                )
                    else:
                        ax.scatter(
                            xi,
                            yi,
                            c=_UPPER_SCATTER_BLUE,
                            s=1.2,
                            alpha=0.05,
                            linewidths=0,
                            rasterized=True,
                        )
                else:
                    if discrete:
                        assert point_colors is not None
                        ax.scatter(
                            xi,
                            yi,
                            c=point_colors,
                            s=2.5,
                            alpha=0.2,
                            linewidths=0,
                            rasterized=True,
                        )
                    else:
                        assert (
                            cvals_plot is not None
                            and cmin is not None
                            and cmax is not None
                        )
                        sc = ax.scatter(
                            xi,
                            yi,
                            c=cvals_plot,
                            cmap="viridis",
                            vmin=cmin,
                            vmax=cmax,
                            s=2.5,
                            alpha=0.2,
                            linewidths=0,
                            rasterized=True,
                        )
                        if lower_cbar_mappable is None:
                            lower_cbar_mappable = sc

                ax.set_xticks([])
                ax.set_yticks([])
                # Equal data aspect without shrinking panels unevenly (``box`` gives ragged cell sizes).
                ax.set_aspect("equal", adjustable="datalim")
                ax.set_box_aspect(1)
                ax.set_facecolor(_DASHBOARD_CREAM)
                for spine in ax.spines.values():
                    spine.set_linewidth(1.05)
                    spine.set_color("#5a6b7a")

        fig.subplots_adjust(
            left=left_m,
            right=right_axes,
            top=top_m,
            bottom=bottom_m,
            wspace=_PAIR_CELL_WSPACE,
            hspace=_PAIR_CELL_HSPACE,
        )

        max_side_gap_in = inch_per / 6.0
        leg = None
        _plot_center_y_fig = (bottom_m + top_m) / 2.0
        _legend_x_fig = right_axes + 0.018
        if discrete and legend_handles is not None:
            leg = fig.legend(
                handles=legend_handles,
                loc="center left",
                bbox_to_anchor=(_legend_x_fig, _plot_center_y_fig),
                bbox_transform=fig.transFigure,
                ncol=2,
                frameon=True,
                fancybox=False,
                fontsize=22.5,
                borderaxespad=0.2,
                columnspacing=1.5,
                handletextpad=0.65,
            )
            leg.get_frame().set_linewidth(0.8)
            leg.get_frame().set_edgecolor("#8899aa")
            leg.get_frame().set_facecolor(_DASHBOARD_CREAM)
            _refine_pair_grid_right_margin(
                fig,
                grid_side_in=grid_side_in,
                left_m=left_m,
                top_m=top_m,
                bottom_m=bottom_m,
                right_artist=leg,
                max_gap_inches=max_side_gap_in,
            )
            _attach_pair_lower_legend_caption(
                fig,
                leg,
                _lower_legend_covariate_title(lower_color_col, exp),
                title_fontsize=24.0,
            )

        if not discrete and lower_cbar_mappable is not None:
            flat_ax = list(axes.ravel())
            cb = fig.colorbar(
                lower_cbar_mappable,
                ax=flat_ax,
                orientation="vertical",
                fraction=0.04,
                pad=0.03,
                aspect=22,
                shrink=1.0,
                anchor=(0.5, 0.5),
            )
            cb.ax.tick_params(labelsize=22.5, length=4, pad=2)
            cb.ax.set_facecolor(_DASHBOARD_CREAM)
            cb.set_label(
                _lower_legend_covariate_title(lower_color_col, exp),
                fontsize=24,
                rotation=270,
                labelpad=28,
            )
            _refine_pair_grid_right_margin(
                fig,
                grid_side_in=grid_side_in,
                left_m=left_m,
                top_m=top_m,
                bottom_m=bottom_m,
                right_artist=cb.ax,
                max_gap_inches=max_side_gap_in,
            )

        _zk_h = max(11.0, min(16.0, 4.9 * inch_per / 3.05))
        _zk_v = _zk_h
        for k in range(zdim):
            ax = axes[k, k]
            _zlab = f"zdim {k + 1}"
            # Horizontal title inside (top-right).
            ax.text(
                0.97,
                0.96,
                _zlab,
                transform=ax.transAxes,
                fontsize=_zk_h,
                color="#243b53",
                fontweight="bold",
                ha="right",
                va="top",
                clip_on=True,
            )
            # Vertical title along the top of the left edge (row-style).
            ax.text(
                0.05,
                0.96,
                _zlab,
                transform=ax.transAxes,
                fontsize=_zk_v,
                color="#243b53",
                fontweight="bold",
                ha="right",
                va="top",
                rotation=90,
                clip_on=True,
            )

        cell_layout = _axes_cell_bboxes(np.asarray(axes))

        buf = io.BytesIO()
        # Full figure bbox (no tight crop) so figure fractions match PNG pixels for the overlay.
        fig.savefig(
            buf,
            format="png",
            dpi=dpi,
            facecolor=_DASHBOARD_CREAM,
            edgecolor="none",
            pad_inches=0,
        )
        plt.close(fig)
        buf.seek(0)
        return buf.getvalue(), cell_layout
