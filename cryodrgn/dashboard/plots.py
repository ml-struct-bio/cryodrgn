"""Plotly figures for the dashboard; pair grid uses Matplotlib + Seaborn (PNG)."""

from __future__ import annotations

import io
import math
import re
from collections.abc import Collection
from typing import Any, cast

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import seaborn as sns
from seaborn.distributions import (
    _freedman_diaconis_bins,  # type: ignore[attr-defined]
)
from seaborn.palettes import blend_palette
from seaborn.utils import set_hls_values

from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.mpl_style import ezlab_matplotlib_rc


def _lower_color_series_is_discrete(s: pd.Series) -> bool:
    """Treat integer / categorical / few whole-valued floats as discrete coloring."""
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
# Default sequential scale for continuous covariates (Plotly name + Matplotlib cmap).
DEFAULT_DASHBOARD_CONTINUOUS_PALETTE = "Viridis"
_DASHBOARD_PALETTE_TO_MPL: dict[str, str] = {
    "Viridis": "viridis",
    "Plasma": "plasma",
    "Inferno": "inferno",
    "Magma": "magma",
    "Cividis": "cividis",
    "Turbo": "turbo",
    "Blues": "Blues",
    "Greens": "Greens",
    "Greys": "Greys",
    "Oranges": "Oranges",
    "Purples": "Purples",
    "Reds": "Reds",
    "YlGnBu": "YlGnBu",
    "YlOrRd": "YlOrRd",
    "RdBu": "RdBu",
    "Portland": "turbo",
    "Jet": "jet",
    "Hot": "hot",
    "Blackbody": "afmhot",
    "Electric": "magma",
    "Rainbow": "rainbow",
    "Earth": "gist_earth",
}


def normalize_continuous_palette(raw: str | None) -> str:
    """Map a user/API string to a supported Plotly colorscale name."""
    if not raw or not isinstance(raw, str):
        return DEFAULT_DASHBOARD_CONTINUOUS_PALETTE
    s = raw.strip()
    for pk in _DASHBOARD_PALETTE_TO_MPL:
        if pk.lower() == s.lower():
            return pk
    return DEFAULT_DASHBOARD_CONTINUOUS_PALETTE


def mpl_cmap_for_palette(plotly_palette: str) -> str:
    """Matplotlib colormap name for a Plotly palette (falls back to Viridis)."""
    return _DASHBOARD_PALETTE_TO_MPL.get(
        plotly_palette,
        _DASHBOARD_PALETTE_TO_MPL[DEFAULT_DASHBOARD_CONTINUOUS_PALETTE],
    )


# Pair grid: gutter between cells (fraction of average subplot width/height).
_PAIR_CELL_WSPACE = 0.06
_PAIR_CELL_HSPACE = 0.06

# Figure margins for the pair subplot block (no matplotlib covariate legend — UI).
_PAIR_GRID_LEFT = 0.035
_PAIR_GRID_TOP = 0.98
_PAIR_GRID_BOTTOM = 0.06
_PAIR_GRID_RIGHT = 0.965
_PAIR_GRID_RIGHT_PLACEHOLDER = _PAIR_GRID_RIGHT

_UMAP_RE = re.compile(r"(?i)^umap")

_PLOTLY_FONT = dict(
    family="Barlow, Roboto, -apple-system, BlinkMacSystemFont, sans-serif",
    size=13,
)

# Dashboard-only discrete palette for labels/k-means (kept out of non-dashboard tools).
_DASHBOARD_CHIMERAX_COLORS: tuple[str, ...] = (
    "#949494",
    "#ffd84d",
    "#81bfd1",
    "#929ad1",
    "#c39bc3",
    "#d1a0a0",
    "#92c29a",
    "#b5a086",
    "#829db1",
    "#b5b547",
)


def _dashboard_chimerax_colors(k: int) -> list[str]:
    if k <= 0:
        return []
    return [
        _DASHBOARD_CHIMERAX_COLORS[i % len(_DASHBOARD_CHIMERAX_COLORS)]
        for i in range(k)
    ]


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


def _scatter3d_marker_size_opacity(
    n_points: int, *, point_cap: int
) -> tuple[float, float]:
    """Scale 3D point glyphs for readability — larger / more opaque when fewer points."""
    n = max(1, int(n_points))
    hi_n = float(max(2, int(point_cap)))
    if n >= hi_n:
        return 0.75, 0.4
    t = math.log10(float(n)) / math.log10(hi_n)
    t = max(0.0, min(1.0, t))
    size = 6.0 + t * (0.75 - 6.0)
    opacity = 0.95 + t * (0.4 - 0.95)
    return float(size), float(opacity)


def _mpl_pair_grid_marker_s_alpha(
    n_points: int,
    *,
    point_cap: int,
    dense_s: float,
    dense_a: float,
    sparse_s: float,
    sparse_a: float,
) -> tuple[float, float]:
    """Matplotlib scatter ``s`` / ``alpha`` vs row count (pair-grid cells)."""
    n = max(1, int(n_points))
    hi_n = float(max(2, int(point_cap)))
    if n >= hi_n:
        return dense_s, dense_a
    t = math.log10(float(n)) / math.log10(hi_n)
    t = max(0.0, min(1.0, t))
    s = sparse_s + t * (dense_s - sparse_s)
    a = sparse_a + t * (dense_a - sparse_a)
    return float(s), float(a)


def plot_df_row_indices_for_explorer_scatter(
    plot_df: pd.DataFrame,
    max_points: int = 200_000,
    *,
    seed: int = 0,
) -> np.ndarray:
    """``plot_df`` row positions shown in the particle explorer scatter (``api_scatter`` defaults).

    Matches :func:`scatter_json` / ``api_scatter`` when the request does not set
    ``filter_ui`` or ``full`` — same ``max_points`` cap and fixed RNG ``seed``.

    Returns
    -------
    np.ndarray
        Integer indices into ``plot_df`` (same as ``customdata[:, 1]`` for those points).
    """
    _, row_indices = _subsample(plot_df, max_points, seed=seed)
    return row_indices


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


def _labels_colors_and_legend_items(
    values: pd.Series,
) -> tuple[list[str], list[dict[str, str]]]:
    """Point colors and legend items for k-means-style integer labels."""
    codes_arr, uniques = pd.factorize(values, sort=True)
    codes = np.asarray(codes_arr, dtype=np.int64)
    pal = _dashboard_chimerax_colors(max(len(uniques), 1))
    colors: list[str] = []
    for code in codes:
        if int(code) < 0:
            colors.append("#aab4bf")
        else:
            colors.append(pal[int(code) % len(pal)])
    items: list[dict[str, str]] = []
    for idx, u in enumerate(uniques):
        if pd.isna(u):
            continue
        items.append(
            {
                "label": _lower_legend_entry_label("labels", u),
                "color": pal[idx % len(pal)],
            }
        )
    if bool(values.isna().any()):
        items.append({"label": "(missing)", "color": "#aab4bf"})
    return colors, items


def _continuous_series_stats(values: pd.Series) -> tuple[np.ndarray, float, float]:
    """Return numeric array plus robust (min,max) bounds for continuous coloring."""
    color_num = cast(pd.Series, pd.to_numeric(values, errors="coerce"))
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
    return cvals, cmin, cmax


def _pair_grid_finite_axis_lim(values: pd.Series) -> tuple[float, float]:
    """Min/max of finite entries (for stable pair-grid axes across colour filters)."""
    v = values.to_numpy(dtype=np.float64)
    m = np.isfinite(v)
    if not np.any(m):
        return -1.0, 1.0
    lo = float(np.min(v[m]))
    hi = float(np.max(v[m]))
    if hi <= lo:
        hi = lo + 1e-9
    return lo, hi


# Fraction of the axis span added as padding on each side (pair grid + pinned 3-D axes).
_PAIR_GRID_AXIS_PAD_FRAC = 0.045


def _pair_grid_axis_lim_padded(
    values: pd.Series,
    *,
    pad_frac: float | None = None,
) -> tuple[float, float]:
    """Finite axis bounds widened symmetrically so points clear the panel edges."""
    pf = _PAIR_GRID_AXIS_PAD_FRAC if pad_frac is None else pad_frac
    pf = max(0.0, float(pf))
    lo, hi = _pair_grid_finite_axis_lim(values)
    span = hi - lo
    pad = span * pf if span > 0 else max(abs(lo), abs(hi), 1.0) * pf
    return lo - pad, hi + pad


def _pair_grid_square_xy_limits(
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Align x/y data ranges so span matches on square ``equal``-aspect axes.

    Matplotlib ≥3.10 logs when ``set_xlim`` / ``set_ylim`` are later overridden by
    ``set_aspect(..., adjustable=\"datalim\")``. Matching spans avoids that adjustment.
    """
    lo_x, hi_x = float(xlim[0]), float(xlim[1])
    lo_y, hi_y = float(ylim[0]), float(ylim[1])
    span_x = hi_x - lo_x
    span_y = hi_y - lo_y
    if not math.isfinite(span_x) or span_x <= 0:
        span_x = max(abs(lo_x), abs(hi_x), 1.0) * 1e-9 + 1e-15
        hi_x = lo_x + span_x
    if not math.isfinite(span_y) or span_y <= 0:
        span_y = max(abs(lo_y), abs(hi_y), 1.0) * 1e-9 + 1e-15
        hi_y = lo_y + span_y
    cx = 0.5 * (lo_x + hi_x)
    cy = 0.5 * (lo_y + hi_y)
    if span_x >= span_y:
        half = 0.5 * span_x
        return (lo_x, hi_x), (cy - half, cy + half)
    half = 0.5 * span_y
    return (cx - half, cx + half), (lo_y, hi_y)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def pair_grid_skeleton_placeholder_layout(
    zdim: int,
) -> list[dict[str, float]]:
    """Cell bboxes for the HTML skeleton before the first PNG loads.

    Mirrors the Matplotlib grid setup in :func:`pair_grid_png` (figure size, per-axis
    ``set_aspect`` / ``set_box_aspect``, spine widths, ``subplots_adjust``) so
    normalized axes positions match the delivered PNG. Uses symmetric axis limits
    instead of real data so ``adjustable=\"datalim\"`` matches typical behaviour.
    """
    if zdim < 1:
        return []
    inch_per = max(2.45, min(3.35, 13.0 / max(zdim, 1)))
    left_m = _PAIR_GRID_LEFT
    top_m = min(_PAIR_GRID_TOP, 0.955)
    bottom_m = _PAIR_GRID_BOTTOM
    right_axes = _PAIR_GRID_RIGHT
    w_frac = right_axes - left_m
    h_frac = top_m - bottom_m
    grid_side_in = inch_per * zdim
    fig_w = grid_side_in / w_frac
    fig_h = grid_side_in / h_frac
    with ezlab_matplotlib_rc():
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
        for i in range(zdim):
            for j in range(zdim):
                ax = axes[i, j]
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xlim(-6.0, 6.0)
                ax.set_ylim(-6.0, 6.0)
                ax.set_aspect("equal", adjustable="datalim")
                ax.set_box_aspect(1)
                ax.set_facecolor(_DASHBOARD_CREAM)
                is_diag = i == j
                for spine in ax.spines.values():
                    spine.set_linewidth(2.1 if is_diag else 1.05)
                    spine.set_color("#334e68" if is_diag else "#5a6b7a")
        fig.subplots_adjust(
            left=left_m,
            right=right_axes,
            top=top_m,
            bottom=bottom_m,
            wspace=_PAIR_CELL_WSPACE,
            hspace=_PAIR_CELL_HSPACE,
        )
        cells = _axes_cell_bboxes(np.asarray(axes))
        plt.close(fig)
    return cells


def pair_grid_margin_fractions_for_js() -> dict[str, float]:
    """Figure-edge parameters for pair-plot skeleton labels (aligned with ``_draw_pair_grid_edge_labels``)."""
    edge_gap = 0.008
    return {
        "left_m": float(_PAIR_GRID_LEFT),
        "top_m": float(min(_PAIR_GRID_TOP, 0.955)),
        "bottom_m": float(_PAIR_GRID_BOTTOM),
        "right_axes": float(_PAIR_GRID_RIGHT),
        "edge_gap": float(edge_gap),
        "right_gap": float(edge_gap * 0.6),
    }


def pair_grid_figure_aspect_ratio(zdim: int) -> float:
    """Figure width ÷ height for :func:`pair_grid_png` (stable layout before the PNG loads)."""
    if zdim < 1:
        return 1.0
    inch_per = max(2.45, min(3.35, 13.0 / max(zdim, 1)))
    left_m = _PAIR_GRID_LEFT
    top_m = min(_PAIR_GRID_TOP, 0.955)
    bottom_m = _PAIR_GRID_BOTTOM
    right_axes = _PAIR_GRID_RIGHT
    w_frac = right_axes - left_m
    h_frac = top_m - bottom_m
    grid_side_in = inch_per * zdim
    fig_w = grid_side_in / w_frac
    fig_h = grid_side_in / h_frac
    return float(fig_w / fig_h) if fig_h > 0 else 1.0


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


def covariate_row_filter_key(color_col: str, u: Any) -> str:
    """Stable string key for a covariate value (matches dashboard scatter customdata)."""
    if pd.isna(u):
        return "__na__"
    if color_col == "labels":
        return str(_lower_legend_entry_label(color_col, u))
    num = pd.to_numeric(pd.Series([u], dtype=object), errors="coerce").iloc[0]
    if pd.notna(num):
        fv = float(num)
        if np.isfinite(fv) and abs(fv - round(fv)) < 1e-8:
            return str(int(round(fv)))
        return str(fv)
    return str(u)


def plot_df_color_filter_mask(
    df: pd.DataFrame,
    color_col: str,
    color_filter: dict[str, Any] | None,
) -> np.ndarray:
    """Boolean mask aligned with ``df`` rows (True = include)."""
    if not color_filter:
        return np.ones(len(df), dtype=bool)
    kind = color_filter.get("kind")
    if kind == "threshold":
        level = float(color_filter["level"])
        use_max = bool(color_filter.get("use_max"))
        series = pd.to_numeric(df[color_col], errors="coerce")
        if use_max:
            mask = (series <= level) & series.notna()
        else:
            mask = (series >= level) & series.notna()
        return mask.to_numpy(dtype=bool)
    if kind == "range":
        r0 = float(color_filter["range_min"])
        r1 = float(color_filter["range_max"])
        if r0 > r1:
            r0, r1 = r1, r0
        series = pd.to_numeric(df[color_col], errors="coerce")
        inside = (series >= r0) & (series <= r1) & series.notna()
        if bool(color_filter.get("invert_range")):
            mask = (~inside) & series.notna()
        else:
            mask = inside
        return mask.to_numpy(dtype=bool)
    if kind == "discrete":
        keys_inc = frozenset(str(k) for k in color_filter.get("keys") or ())
        ser = df[color_col]
        keys_series = ser.map(lambda u, cc=color_col: covariate_row_filter_key(cc, u))
        return keys_series.isin(keys_inc).to_numpy(dtype=bool)
    return np.ones(len(df), dtype=bool)


def _discrete_legend_sort_tuple(key: str) -> tuple:
    if key == "__na__":
        return (2, "")
    try:
        return (0, float(key))
    except ValueError:
        return (1, str(key))


def covariate_legend_context_payload(
    exp: DashboardExperiment,
    column: str,
    *,
    violin_sample_max: int = 12_000,
) -> dict[str, Any]:
    """JSONable payload for building the colour histogram / discrete toggles."""
    if column not in exp.plot_df.columns:
        raise ValueError("Unknown column.")
    s = exp.plot_df[column]
    if _lower_color_series_is_discrete(s):
        codes_arr, uniques = pd.factorize(s, sort=True)
        codes = np.asarray(codes_arr, dtype=np.int64)
        pal = _dashboard_chimerax_colors(max(len(uniques), 1))
        cats: list[dict[str, Any]] = []
        for idx, u in enumerate(uniques):
            if pd.isna(u):
                continue
            cnt = int(np.sum(codes == idx))
            if cnt == 0:
                continue
            cats.append(
                {
                    "key": covariate_row_filter_key(column, u),
                    "label": _lower_legend_entry_label(column, u),
                    "count": cnt,
                    "color": pal[idx % len(pal)],
                }
            )
        if bool(s.isna().any()):
            cats.append(
                {
                    "key": "__na__",
                    "label": "(missing)",
                    "count": int(s.isna().sum()),
                    "color": "#aab4bf",
                }
            )
        cats.sort(key=lambda c: _discrete_legend_sort_tuple(str(c["key"])))
        n_pts = len(exp.plot_df)
        point_cap = max(n_pts, 1)
        _, chip_blend_alpha = _mpl_pair_grid_marker_s_alpha(
            n_pts,
            point_cap=point_cap,
            dense_s=2.5,
            dense_a=0.2,
            sparse_s=10.0,
            sparse_a=0.78,
        )
        return {
            "mode": "discrete",
            "categories": cats,
            "chip_blend_alpha": float(chip_blend_alpha),
            "chip_blend_bg": _DASHBOARD_CREAM,
        }
    vals = pd.to_numeric(s, errors="coerce").dropna().to_numpy(dtype=np.float64)
    if vals.size == 0:
        return {"mode": "continuous", "values": []}
    rng = np.random.default_rng(0)
    if vals.size > violin_sample_max:
        vals = rng.choice(vals, size=violin_sample_max, replace=False)
    return {"mode": "continuous", "values": vals.tolist()}


_EDGE_LABEL_COLOR = "#243b53"


def _draw_pair_grid_edge_labels(
    fig: plt.Figure,
    axes: np.ndarray,
    *,
    zdim: int,
    left_m: float,
    top_m: float,
    bottom_m: float,
    right_axes: float,
    fontsize: float,
    edge_gap: float = 0.008,
) -> None:
    """Draw ``zdim N`` labels on all four edges of the pair grid."""
    right_gap = edge_gap * 0.6
    common = dict(fontsize=fontsize, color=_EDGE_LABEL_COLOR, fontweight="bold")
    # Top edge (skip col 0 to avoid the diagonal).
    for j in range(1, zdim):
        bb = axes[0, j].get_position()
        fig.text(
            bb.x0 + bb.width / 2.0,
            top_m + edge_gap,
            f"zdim {j + 1}",
            ha="center",
            va="bottom",
            **common,
        )
    # Bottom edge (skip the last column to avoid the diagonal).
    for j in range(zdim - 1):
        bb = axes[zdim - 1, j].get_position()
        fig.text(
            bb.x0 + bb.width / 2.0,
            bottom_m - edge_gap,
            f"zdim {j + 1}",
            ha="center",
            va="top",
            **common,
        )
    # Right edge (skip the last row to avoid the diagonal).
    for i in range(zdim - 1):
        bb = axes[i, zdim - 1].get_position()
        fig.text(
            right_axes + right_gap,
            bb.y0 + bb.height / 2.0,
            f"zdim {i + 1}",
            ha="left",
            va="center",
            rotation=270,
            **common,
        )
    # Left edge (skip the first row to avoid the diagonal).
    for i in range(1, zdim):
        bb = axes[i, 0].get_position()
        fig.text(
            left_m - edge_gap,
            bb.y0 + bb.height / 2.0,
            f"zdim {i + 1}",
            ha="right",
            va="center",
            rotation=90,
            **common,
        )


def _draw_pair_grid_diagonal_legends(
    axes: np.ndarray,
    *,
    zdim: int,
    mpl_cmap: str,
    diagonal_color_ranges: list[tuple[float, float] | None],
    label_fontsize: float,
) -> None:
    """Inset color-bar + label in the bottom-left of each diagonal panel."""
    grad_left, grad_bottom = 0.115, 0.028
    grad_width, grad_height = 0.24, 0.044
    side_gap = grad_bottom * 0.5
    bbox = {
        "boxstyle": "round,pad=0.16",
        "facecolor": (250 / 255, 248 / 255, 244 / 255, 0.9),
        "edgecolor": "none",
    }
    for k in range(zdim):
        ax = axes[k, k]
        ax.text(
            grad_left + grad_width / 2.0,
            0.10,
            f"zdim {k + 1}",
            transform=ax.transAxes,
            fontsize=max(9.6, label_fontsize * 0.864),
            color=_EDGE_LABEL_COLOR,
            fontweight="bold",
            ha="center",
            va="bottom",
            clip_on=True,
            bbox=bbox,
        )
        zrange = diagonal_color_ranges[k]
        if zrange is None:
            continue
        zmin, zmax = zrange
        grad_ax = ax.inset_axes(
            [grad_left, grad_bottom, grad_width, grad_height],
            transform=ax.transAxes,
        )
        grad = np.linspace(0.0, 1.0, 128, dtype=np.float64).reshape(1, -1)
        grad_ax.imshow(grad, cmap=mpl_cmap, aspect="auto", origin="lower")
        grad_ax.set_xticks([])
        grad_ax.set_yticks([])
        for spine in grad_ax.spines.values():
            spine.set_linewidth(0.6)
            spine.set_color("#6b7c8d")
        y_mid = grad_bottom + grad_height / 2.0
        tick_kwargs = dict(
            transform=ax.transAxes,
            fontsize=max(8.4, label_fontsize * 0.672),
            color=_EDGE_LABEL_COLOR,
            va="center",
            clip_on=True,
        )
        ax.text(grad_left - side_gap, y_mid, f"{zmin:.2g}", ha="right", **tick_kwargs)
        ax.text(
            grad_left + grad_width + side_gap,
            y_mid,
            f"{zmax:.2g}",
            ha="left",
            **tick_kwargs,
        )


def scatter_json(
    exp: DashboardExperiment,
    xcol: str,
    ycol: str,
    color_col: str | None,
    max_points: int | None = 120_000,
    *,
    preselect_plot_df_rows: Collection[int] | None = None,
    use_webgl: bool = True,
    marker_size: float = 4,
    continuous_palette: str | None = None,
) -> str:
    plotly_cs = normalize_continuous_palette(continuous_palette)
    sub, row_indices = _subsample(exp.plot_df, max_points, seed=0)

    if color_col and color_col != "none" and color_col in sub.columns:
        if color_col == "labels":
            colors, _legend_items = _labels_colors_and_legend_items(sub[color_col])
            marker = dict(
                size=marker_size,
                opacity=0.35,
                color=colors,
            )
        else:
            marker = dict(
                size=marker_size,
                opacity=0.35,
                color=sub[color_col],
                colorscale=plotly_cs,
            )
    else:
        marker = dict(size=marker_size, opacity=0.35, color="#4a5568")

    idx_arr = sub["index"].to_numpy()
    row_arr = np.asarray(row_indices, dtype=np.int64)
    n_pts = len(sub)
    if color_col and color_col != "none" and color_col in sub.columns:
        disp = np.empty(n_pts, dtype=object)
        if color_col == "labels":
            cv = sub[color_col]
            for i in range(n_pts):
                disp[i] = _lower_legend_entry_label("labels", cv.iloc[i])
        else:
            color_num = cast(pd.Series, pd.to_numeric(sub[color_col], errors="coerce"))
            for i in range(n_pts):
                v = color_num.iloc[i]
                if pd.isna(v):
                    disp[i] = None
                else:
                    fv = float(v)
                    disp[i] = fv if np.isfinite(fv) else None
        customdata = np.column_stack([idx_arr, row_arr, disp])
    else:
        customdata = np.column_stack([idx_arr, row_arr])
    # Scattergl can leave Plotly.react() pending on some browsers/GPUs.
    trace_cls = go.Scattergl if use_webgl else go.Scatter
    sc = trace_cls(
        x=sub[xcol],
        y=sub[ycol],
        mode="markers",
        customdata=customdata,
        hovertemplate="row %{customdata[1]} · index %{customdata[0]}<extra></extra>",
        marker=marker,
    )

    xaxis_kw: dict[str, Any] = dict(title=xcol)
    yaxis_kw: dict[str, Any] = dict(title=ycol)
    if _UMAP_RE.match(xcol):
        xaxis_kw["showticklabels"] = False
        xaxis_kw["showgrid"] = False
    if _UMAP_RE.match(ycol):
        yaxis_kw["showticklabels"] = False
        yaxis_kw["showgrid"] = False

    layout_kw: dict[str, Any] = dict(
        template="plotly_white",
        paper_bgcolor=_DASHBOARD_CREAM,
        margin=dict(l=50, r=20, t=16, b=50),
        title=None,
        xaxis=xaxis_kw,
        yaxis=yaxis_kw,
        dragmode="lasso",
        uirevision="scatter",
        font=_PLOTLY_FONT,
        showlegend=False,
    )
    layout_meta: dict[str, Any] = {}
    if preselect_plot_df_rows is not None:
        want = frozenset(int(x) for x in preselect_plot_df_rows)
        layout_meta["cdrgn_preselected"] = [
            int(i) for i in range(len(row_indices)) if int(row_indices[i]) in want
        ]
    if color_col and color_col != "none" and color_col in sub.columns:
        layout_meta["cdrgn_color_mode"] = (
            "discrete" if color_col == "labels" else "continuous"
        )
    if layout_meta:
        layout_kw["meta"] = layout_meta

    fig = go.Figure(sc)
    fig.update_layout(**layout_kw)
    return _plotly_to_json(fig)


def _validate_three_latent_axes(
    exp: DashboardExperiment,
    df: pd.DataFrame,
    axes: tuple[str, str, str],
) -> None:
    """Validate three distinct latent ``z*`` axes against the experiment table."""
    z_allowed = {f"z{i}" for i in range(int(exp.z.shape[1]))}
    for col in axes:
        if col not in z_allowed:
            raise ValueError(f"Axis {col!r} is not a latent dimension for this run.")
        if col not in df.columns:
            raise ValueError(f"Missing column {col!r} in analysis table.")
    if len(set(axes)) < 3:
        raise ValueError("Choose three distinct latent axes.")


def scatter3d_z_json(
    exp: DashboardExperiment,
    xcol: str,
    ycol: str,
    zcol: str,
    color_col: str | None,
    max_points: int = 120_000,
    *,
    continuous_palette: str | None = None,
    color_filter: dict[str, Any] | None = None,
) -> str:
    """Interactive 3D scatter of three latent ``z*`` columns."""
    plotly_cs = normalize_continuous_palette(continuous_palette)
    df_all = exp.plot_df
    df = df_all
    _validate_three_latent_axes(exp, df, (xcol, ycol, zcol))

    pinned_color_range: tuple[float, float] | None = None
    if color_filter and color_col and color_col != "none":
        if color_col != "labels" and color_col in df.columns:
            _, pcmin, pcmax = _continuous_series_stats(df_all[color_col])
            pinned_color_range = (pcmin, pcmax)

    if color_filter and color_col and color_col != "none":
        mask = plot_df_color_filter_mask(df_all, color_col, color_filter)
        if not bool(np.any(mask)):
            raise ValueError("No particles match the current colour covariate filter.")
        df = df_all.loc[mask].reset_index(drop=True)

    sub, row_indices = _subsample(df, max_points, seed=1)
    cap = max_points if max_points is not None else max(len(sub), 1)
    msize, mopacity = _scatter3d_marker_size_opacity(len(sub), point_cap=cap)

    legend_meta: dict[str, Any] | None = None
    if color_col and color_col != "none" and color_col in sub.columns:
        if color_col == "labels":
            colors, items = _labels_colors_and_legend_items(sub[color_col])
            marker = dict(
                size=msize,
                opacity=mopacity,
                color=colors,
            )
            legend_meta = {
                "type": "discrete",
                "title": "k-means labels",
                "items": items,
            }
        else:
            _cvals, cmin, cmax = _continuous_series_stats(sub[color_col])
            marker = dict(
                size=msize,
                opacity=mopacity,
                color=sub[color_col],
                colorscale=plotly_cs,
            )
            if pinned_color_range is not None:
                cmin, cmax = pinned_color_range
                marker["cmin"] = cmin
                marker["cmax"] = cmax
            legend_meta = {
                "type": "continuous",
                "title": color_col,
                "min": cmin,
                "max": cmax,
            }
    else:
        marker = dict(size=msize, opacity=mopacity, color="#4a5568")

    idx_cd = sub["index"].to_numpy()
    row_cd = np.asarray(row_indices, dtype=np.int64)
    if (
        color_col
        and color_col != "none"
        and color_col in sub.columns
        and color_col == "labels"
    ):
        cv3 = sub[color_col]
        disp_ca = np.empty(len(sub), dtype=object)
        for ii in range(len(sub)):
            disp_ca[ii] = _lower_legend_entry_label("labels", cv3.iloc[ii])
        customdata = np.column_stack([idx_cd, row_cd, disp_ca])
    else:
        customdata = np.column_stack([idx_cd, row_cd])
    sc = go.Scatter3d(
        x=sub[xcol],
        y=sub[ycol],
        z=sub[zcol],
        mode="markers",
        customdata=customdata,
        hovertemplate="row %{customdata[1]} · index %{customdata[0]}<extra></extra>",
        marker=marker,
    )
    transparent = "rgba(0,0,0,0)"

    def _scene_axis(extra: dict[str, Any] | None = None) -> dict[str, Any]:
        d: dict[str, Any] = dict(backgroundcolor=transparent)
        if extra:
            d.update(extra)
        return d

    scene_axes_extra: dict[str, dict[str, Any]] = {}
    if color_filter:
        scene_axes_extra["xaxis"] = _scene_axis(
            {"range": list(_pair_grid_axis_lim_padded(df_all[xcol]))}
        )
        scene_axes_extra["yaxis"] = _scene_axis(
            {"range": list(_pair_grid_axis_lim_padded(df_all[ycol]))}
        )
        scene_axes_extra["zaxis"] = _scene_axis(
            {"range": list(_pair_grid_axis_lim_padded(df_all[zcol]))}
        )
    else:
        scene_axes_extra["xaxis"] = _scene_axis()
        scene_axes_extra["yaxis"] = _scene_axis()
        scene_axes_extra["zaxis"] = _scene_axis()

    fig = go.Figure(sc)
    fig.update_layout(
        template="plotly_white",
        autosize=True,
        paper_bgcolor=_DASHBOARD_CREAM,
        margin=dict(l=0, r=0, t=0, b=0),
        scene=dict(
            xaxis_title=xcol,
            yaxis_title=ycol,
            zaxis_title=zcol,
            aspectmode="data",
            bgcolor="rgba(250,248,244,0.95)",
            **scene_axes_extra,
        ),
        uirevision="scatter3d_z",
        font=_PLOTLY_FONT,
        meta={"cdrgn_color_legend": legend_meta},
    )
    return _plotly_to_json(fig)


def scatter3d_z_preview_png(
    exp: DashboardExperiment,
    xcol: str,
    ycol: str,
    zcol: str,
    color_col: str | None,
    max_points: int = 120_000,
    *,
    continuous_palette: str | None = None,
    dpi: int = 100,
    elev: float = 22.0,
    azim: float = -65.0,
) -> bytes:
    """Matplotlib 3D scatter PNG using the same subsample/coloring rules as ``scatter3d_z_json``.

    Used for dashboard GIF capture where headless browsers often fail to composite WebGL.
    """
    plotly_cs = normalize_continuous_palette(continuous_palette)
    mpl_cmap_name = mpl_cmap_for_palette(plotly_cs)
    df = exp.plot_df
    _validate_three_latent_axes(exp, df, (xcol, ycol, zcol))

    sub, _row_indices = _subsample(df, max_points, seed=1)
    xs = sub[xcol].to_numpy(dtype=np.float64)
    ys = sub[ycol].to_numpy(dtype=np.float64)
    zs = sub[zcol].to_numpy(dtype=np.float64)

    with ezlab_matplotlib_rc():
        fig = plt.figure(figsize=(5.2, 4.6), facecolor=_DASHBOARD_CREAM)
        ax = fig.add_subplot(111, projection="3d", facecolor=_DASHBOARD_CREAM)

        if color_col and color_col != "none" and color_col in sub.columns:
            if color_col == "labels":
                colors, _legend_items = _labels_colors_and_legend_items(sub[color_col])
                ax.scatter(
                    xs,
                    ys,
                    zs,
                    c=colors,
                    s=0.8,
                    alpha=0.1,
                    linewidths=0,
                    depthshade=True,
                )
            else:
                cvals, cmin, cmax = _continuous_series_stats(sub[color_col])
                norm = mcolors.Normalize(vmin=cmin, vmax=cmax)
                cmap = plt.get_cmap(mpl_cmap_name)
                cplot = np.where(np.isfinite(cvals), cvals, 0.5 * (cmin + cmax))
                m = ax.scatter(
                    xs,
                    ys,
                    zs,
                    c=cplot,
                    s=0.8,
                    cmap=cmap,
                    norm=norm,
                    alpha=0.1,
                    linewidths=0,
                    depthshade=True,
                )
                fig.colorbar(m, ax=ax, shrink=0.5, pad=0.12, label=color_col)
        else:
            ax.scatter(
                xs,
                ys,
                zs,
                c="#4a5568",
                s=0.8,
                alpha=0.1,
                linewidths=0,
                depthshade=True,
            )

        ax.set_xlabel(xcol, fontsize=10)
        ax.set_ylabel(ycol, fontsize=10)
        ax.set_zlabel(zcol, fontsize=10)
        ax.tick_params(axis="both", labelsize=7)
        ax.view_init(elev=float(elev), azim=float(azim))
        fig.tight_layout(pad=0.6)
        buf = io.BytesIO()
        fig.savefig(
            buf,
            format="png",
            dpi=dpi,
            facecolor=_DASHBOARD_CREAM,
            bbox_inches="tight",
        )
        plt.close(fig)

    return buf.getvalue()


def pair_grid_png(
    exp: DashboardExperiment,
    lower_color_col: str,
    diagonal_emb: str,
    upper_style: str,
    dpi: int = 120,
    *,
    continuous_palette: str | None = None,
    color_filter: dict[str, Any] | None = None,
    discrete_label_colors: dict[str, str] | None = None,
) -> tuple[bytes, list[dict[str, float]]]:
    """z_dim × z_dim Matplotlib grid (square cells); upper hex matches ``sns.jointplot(..., kind=\"hex\")``.

    Returns PNG bytes and per-cell axis bboxes in figure coordinates (for HTML overlay alignment).
    """
    mpl_cmap = mpl_cmap_for_palette(normalize_continuous_palette(continuous_palette))
    full_df = exp.plot_df
    zdim = int(exp.z.shape[1])
    zcols = [f"z{i}" for i in range(zdim)]
    for c in zcols:
        if c not in full_df.columns:
            raise ValueError(f"Missing latent column {c} in analysis table.")

    z_axis_lim = [_pair_grid_axis_lim_padded(full_df[zc]) for zc in zcols]

    if lower_color_col not in full_df.columns:
        raise ValueError(f"Unknown color column: {lower_color_col}")
    if color_filter:
        mask = plot_df_color_filter_mask(full_df, lower_color_col, color_filter)
        if not bool(np.any(mask)):
            raise ValueError("No particles match the current colour covariate filter.")
        df = full_df.loc[mask].reset_index(drop=True)
    else:
        df = full_df
    n_pts = len(df)
    point_cap = max(len(full_df), 1)
    lo_s, lo_a = _mpl_pair_grid_marker_s_alpha(
        n_pts,
        point_cap=point_cap,
        dense_s=2.5,
        dense_a=0.2,
        sparse_s=10.0,
        sparse_a=0.78,
    )
    diag_s, diag_a = _mpl_pair_grid_marker_s_alpha(
        n_pts,
        point_cap=point_cap,
        dense_s=2.0,
        dense_a=0.55,
        sparse_s=8.0,
        sparse_a=0.9,
    )
    # Upper-triangle scatter (fixed pastel blue): scale like lower cells so small
    # selections stay readable; keep dense plots lighter than the coloured lower panels.
    up_scatter_s, up_scatter_a = _mpl_pair_grid_marker_s_alpha(
        n_pts,
        point_cap=point_cap,
        dense_s=1.85,
        dense_a=0.085,
        sparse_s=14.5,
        sparse_a=0.88,
    )
    raw_color = df[lower_color_col]
    discrete = _lower_color_series_is_discrete(raw_color)

    point_colors: list[str] | None = None
    cvals_plot: np.ndarray | None = None
    cmin: float | None = None
    cmax: float | None = None

    if discrete:
        raw_full_color = full_df[lower_color_col]
        _, uniques_f = pd.factorize(raw_full_color, sort=True)
        if len(uniques_f) == 0:
            raise ValueError(
                f"Lower-triangle color column `{lower_color_col}` has no values."
            )
        pal = _dashboard_chimerax_colors(max(len(uniques_f), 1))

        def _discrete_hex_for_key(fk: str, palette_idx: int) -> str:
            override = discrete_label_colors.get(fk) if discrete_label_colors else None
            if isinstance(override, str) and re.match(
                r"^#[0-9a-fA-F]{6}$", override.strip()
            ):
                return override.strip().lower()
            return pal[palette_idx % len(pal)]

        # Stable colours keyed like the legend API — never re-factorize the filtered
        # subset (that would remap the sole remaining category to palette slot 0).
        stable_hex_by_key: dict[str, str] = {}
        for idx, u in enumerate(uniques_f):
            if pd.isna(u):
                continue
            fk = covariate_row_filter_key(lower_color_col, u)
            stable_hex_by_key[fk] = _discrete_hex_for_key(fk, idx)

        if bool(raw_full_color.isna().any()):
            fk_na = "__na__"
            ov_na = discrete_label_colors.get(fk_na) if discrete_label_colors else None
            if isinstance(ov_na, str) and re.match(r"^#[0-9a-fA-F]{6}$", ov_na.strip()):
                stable_hex_by_key[fk_na] = ov_na.strip().lower()
            else:
                stable_hex_by_key[fk_na] = "#aab4bf"

        point_colors = []
        for u in raw_color.to_numpy(dtype=object):
            if pd.isna(u):
                point_colors.append(stable_hex_by_key.get("__na__", "#aab4bf"))
            else:
                fk = covariate_row_filter_key(lower_color_col, u)
                point_colors.append(stable_hex_by_key.get(fk, "#aab4bf"))
    else:
        if bool(pd.to_numeric(raw_color, errors="coerce").isna().all()):
            raise ValueError(
                f"Lower-triangle color column `{lower_color_col}` has no numeric values."
            )
        _, cmin, cmax = _continuous_series_stats(full_df[lower_color_col])
        color_num = cast(pd.Series, pd.to_numeric(raw_color, errors="coerce"))
        cvals = np.asarray(color_num, dtype=np.float64)
        cvals_plot = np.where(np.isfinite(cvals), cvals, 0.5 * (cmin + cmax))

    emb = (diagonal_emb or "pc").lower()
    if emb == "umap":
        if exp.umap is None or "UMAP1" not in df.columns or "UMAP2" not in df.columns:
            raise ValueError("UMAP embedding is not available for this run.")
        emb_x = df["UMAP1"].to_numpy(dtype=np.float64)
        emb_y = df["UMAP2"].to_numpy(dtype=np.float64)
        emb_x_lim = _pair_grid_axis_lim_padded(full_df["UMAP1"])
        emb_y_lim = _pair_grid_axis_lim_padded(full_df["UMAP2"])
    elif emb == "pc":
        if "PC1" not in df.columns or "PC2" not in df.columns:
            raise ValueError("PCA components PC1/PC2 are not available.")
        emb_x = df["PC1"].to_numpy(dtype=np.float64)
        emb_y = df["PC2"].to_numpy(dtype=np.float64)
        emb_x_lim = _pair_grid_axis_lim_padded(full_df["PC1"])
        emb_y_lim = _pair_grid_axis_lim_padded(full_df["PC2"])
    else:
        raise ValueError("diagonal_emb must be 'pc' or 'umap'.")

    upper = (upper_style or "scatter").lower()
    if upper not in ("scatter", "hex"):
        raise ValueError("upper_style must be 'scatter' or 'hex'.")

    with ezlab_matplotlib_rc():
        hex_cmap = _pair_jointplot_hex_cmap()

        inch_per = max(2.45, min(3.35, 13.0 / max(zdim, 1)))
        left_m = _PAIR_GRID_LEFT
        # Reserve extra headroom for enlarged top edge labels.
        top_m = min(_PAIR_GRID_TOP, 0.955)
        bottom_m = _PAIR_GRID_BOTTOM
        right_axes = _PAIR_GRID_RIGHT

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

        diagonal_color_ranges: list[tuple[float, float] | None] = [None] * zdim

        for i in range(zdim):
            for j in range(zdim):
                ax = axes[i, j]
                xi = df[zcols[j]].to_numpy(dtype=np.float64)
                yi = df[zcols[i]].to_numpy(dtype=np.float64)
                if i == j:
                    zi = df[zcols[i]].to_numpy(dtype=np.float64)
                    zi_full = full_df[zcols[i]].to_numpy(dtype=np.float64)
                    zm_full = np.isfinite(zi_full)
                    z_scatter_kw: dict[str, Any] = {}
                    if np.any(zm_full):
                        zmin_g = float(np.min(zi_full[zm_full]))
                        zmax_g = float(np.max(zi_full[zm_full]))
                        if zmax_g <= zmin_g:
                            zmax_g = zmin_g + 1e-9
                        diagonal_color_ranges[i] = (zmin_g, zmax_g)
                        z_scatter_kw["vmin"] = zmin_g
                        z_scatter_kw["vmax"] = zmax_g
                    ax.scatter(
                        emb_x,
                        emb_y,
                        c=zi,
                        cmap=mpl_cmap,
                        s=diag_s,
                        alpha=diag_a,
                        linewidths=0,
                        rasterized=True,
                        **z_scatter_kw,
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
                            s=up_scatter_s,
                            alpha=up_scatter_a,
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
                            s=lo_s,
                            alpha=lo_a,
                            linewidths=0,
                            rasterized=True,
                        )
                    else:
                        assert (
                            cvals_plot is not None
                            and cmin is not None
                            and cmax is not None
                        )
                        ax.scatter(
                            xi,
                            yi,
                            c=cvals_plot,
                            cmap=mpl_cmap,
                            vmin=cmin,
                            vmax=cmax,
                            s=lo_s,
                            alpha=lo_a,
                            linewidths=0,
                            rasterized=True,
                        )

                ax.set_xticks([])
                ax.set_yticks([])
                # Freeze axes to the unfiltered table so legend selections only subset points.
                if i == j:
                    xl, yl = _pair_grid_square_xy_limits(emb_x_lim, emb_y_lim)
                    ax.set_xlim(xl[0], xl[1])
                    ax.set_ylim(yl[0], yl[1])
                else:
                    xl, yl = _pair_grid_square_xy_limits(z_axis_lim[j], z_axis_lim[i])
                    ax.set_xlim(xl[0], xl[1])
                    ax.set_ylim(yl[0], yl[1])
                # Equal data aspect without shrinking panels unevenly (``box`` gives ragged cell sizes).
                ax.set_aspect("equal", adjustable="datalim")
                ax.set_box_aspect(1)
                ax.set_facecolor(_DASHBOARD_CREAM)
                is_diag = i == j
                for spine in ax.spines.values():
                    spine.set_linewidth(2.1 if is_diag else 1.05)
                    spine.set_color("#334e68" if is_diag else "#5a6b7a")

        fig.subplots_adjust(
            left=left_m,
            right=right_axes,
            top=top_m,
            bottom=bottom_m,
            wspace=_PAIR_CELL_WSPACE,
            hspace=_PAIR_CELL_HSPACE,
        )

        label_fs = max(11.0, min(16.0, 4.9 * inch_per / 3.05))
        _draw_pair_grid_edge_labels(
            fig,
            axes,
            zdim=zdim,
            left_m=left_m,
            top_m=top_m,
            bottom_m=bottom_m,
            right_axes=right_axes,
            fontsize=label_fs * 1.5,
        )
        _draw_pair_grid_diagonal_legends(
            axes,
            zdim=zdim,
            mpl_cmap=mpl_cmap,
            diagonal_color_ranges=diagonal_color_ranges,
            label_fontsize=label_fs,
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
