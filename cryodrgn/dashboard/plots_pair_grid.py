"""Latent pair-grid Matplotlib figures (PNG) and layout helpers."""

from __future__ import annotations

import functools
import io
import json
import threading
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, cast

import matplotlib
import matplotlib.colors as mcolors
import matplotlib.image as mimage
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from PIL import Image
from seaborn.distributions import (
    _freedman_diaconis_bins,  # type: ignore[attr-defined]
)
from seaborn.palettes import blend_palette
from seaborn.utils import set_hls_values

from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.palette_config import (
    mpl_cmap_for_palette,
    normalize_continuous_palette,
)
from cryodrgn.dashboard.plots_color_covariate import (
    _continuous_series_stats,
    _lower_color_series_is_discrete,
    _scatter_discrete_marker_arrays,
    plot_df_color_filter_mask,
)
from cryodrgn.dashboard.plots_figure_utils import (
    _UPPER_SCATTER_BLUE,
    _axes_cell_bboxes,
    _DASHBOARD_CREAM,
    _pair_grid_axis_lim_padded,
    _pair_grid_square_xy_limits,
    _mpl_pair_grid_marker_s_alpha,
    ezlab_matplotlib_rc,
)

# Pair grid: gutter between cells (fraction of average subplot width/height).
_PAIR_CELL_WSPACE = 0.06
_PAIR_CELL_HSPACE = 0.06

# Figure margins for the pair subplot block (no matplotlib covariate legend — UI).
_PAIR_GRID_LEFT = 0.035
_PAIR_GRID_TOP = 0.98
_PAIR_GRID_BOTTOM = 0.06
_PAIR_GRID_RIGHT = 0.965


def _pair_jointplot_hex_cmap(color: str | None = None):
    """Same sequential colormap as ``sns.jointplot(..., kind=\"hex\")``.

    Uses default ``color`` when omitted.
    """
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


def _pair_grid_geometry(zdim: int) -> dict[str, float]:
    """Figure size (inches) and subplot margins (figure fractions) of the pair grid."""
    inch_per = max(2.45, min(3.35, 13.0 / max(zdim, 1)))
    left, right = _PAIR_GRID_LEFT, _PAIR_GRID_RIGHT
    # Extra headroom for the enlarged top edge labels.
    top, bottom = min(_PAIR_GRID_TOP, 0.955), _PAIR_GRID_BOTTOM
    grid_side_in = inch_per * zdim
    return {
        "inch_per": inch_per,
        "left": left,
        "right": right,
        "top": top,
        "bottom": bottom,
        "fig_w": grid_side_in / (right - left),
        "fig_h": grid_side_in / (top - bottom),
    }


def _new_pair_grid_figure(
    zdim: int,
    cell_limits: Callable[[int, int], tuple[tuple[float, float], tuple[float, float]]],
) -> tuple[Any, np.ndarray]:
    """Figure and ``zdim`` × ``zdim`` empty axes, laid out as in the delivered PNG.

    ``cell_limits(i, j)`` gives the x and y limits of the cell in row ``i``, column
    ``j``. Call inside :func:`ezlab_matplotlib_rc`.
    """
    geom = _pair_grid_geometry(zdim)
    with sns.axes_style(
        "white",
        rc={"axes.facecolor": _DASHBOARD_CREAM, "figure.facecolor": _DASHBOARD_CREAM},
    ):
        fig, axes = plt.subplots(
            zdim,
            zdim,
            figsize=(geom["fig_w"], geom["fig_h"]),
            squeeze=False,
            constrained_layout=False,
        )
    fig.patch.set_facecolor(_DASHBOARD_CREAM)
    for i in range(zdim):
        for j in range(zdim):
            ax = axes[i, j]
            ax.set_xticks([])
            ax.set_yticks([])
            xl, yl = cell_limits(i, j)
            ax.set_xlim(xl[0], xl[1])
            ax.set_ylim(yl[0], yl[1])
            # Equal data aspect without shrinking panels unevenly
            # (``box`` gives ragged cell sizes).
            ax.set_aspect("equal", adjustable="datalim")
            ax.set_box_aspect(1)
            ax.set_facecolor(_DASHBOARD_CREAM)
            is_diag = i == j
            for spine in ax.spines.values():
                spine.set_linewidth(2.1 if is_diag else 1.05)
                spine.set_color("#334e68" if is_diag else "#5a6b7a")
    fig.subplots_adjust(
        left=geom["left"],
        right=geom["right"],
        top=geom["top"],
        bottom=geom["bottom"],
        wspace=_PAIR_CELL_WSPACE,
        hspace=_PAIR_CELL_HSPACE,
    )
    return fig, axes


def pair_grid_skeleton_placeholder_layout(
    zdim: int,
) -> list[dict[str, float]]:
    """Cell bboxes for the HTML skeleton before the first PNG loads.

    Uses the figure of :func:`pair_grid_png` with symmetric axis limits instead of
    real data, so ``adjustable="datalim"`` behaves as it typically does there.
    """
    if zdim < 1:
        return []
    return [dict(c) for c in _pair_grid_skeleton_cells(int(zdim))]


@functools.lru_cache(maxsize=8)
def _pair_grid_skeleton_cells(zdim: int) -> tuple[dict[str, float], ...]:
    # Building the figure takes about 0.3 s; the page asks on every visit.
    with ezlab_matplotlib_rc():
        fig, axes = _new_pair_grid_figure(zdim, lambda i, j: ((-6.0, 6.0), (-6.0, 6.0)))
        cells = _axes_cell_bboxes(np.asarray(axes))
        plt.close(fig)
    return tuple(cells)


def pair_grid_margin_fractions_for_js() -> dict[str, float]:
    """Figure-edge parameters for pair-plot skeleton labels.

    Aligned with ``_draw_pair_grid_edge_labels``.
    """
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
    """Figure width ÷ height for :func:`pair_grid_png`.

    Stable layout before the PNG loads.
    """
    if zdim < 1:
        return 1.0
    geom = _pair_grid_geometry(zdim)
    return float(geom["fig_w"] / geom["fig_h"]) if geom["fig_h"] > 0 else 1.0


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


@functools.lru_cache(maxsize=64)
def _marker_footprint(diameter_px: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Pixel offsets and covered fractions of a round marker centred on a pixel.

    The fractions are what an anti-aliased renderer gives a disc of this diameter
    (8 × 8 samples per pixel); pixels covered by less than 1% are left out.
    """
    radius = diameter_px / 2.0
    k = int(np.ceil(radius + 0.5))
    sub = (np.arange(8) + 0.5) / 8.0 - 0.5
    offs = np.arange(-k, k + 1)
    sx = (offs[:, None] + sub[None, :]).ravel()  # sample positions along one axis
    inside = (sx[:, None] ** 2 + sx[None, :] ** 2) <= radius * radius
    n = offs.size
    frac = inside.reshape(n, 8, n, 8).mean(axis=(1, 3))  # [row offset, column offset]
    dy, dx = np.meshgrid(offs, offs, indexing="ij")
    keep = frac > 0.01
    return dx[keep].astype(np.intp), dy[keep].astype(np.intp), frac[keep]


def _panel_pixels(ax: Any) -> dict[str, float] | None:
    """Data limits and output-pixel box of ``ax``, as the PNG will lay them out."""
    ax.apply_aspect()
    bb = ax.get_window_extent()
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    c0, c1 = int(np.floor(bb.x0)), int(np.ceil(bb.x1))
    r0, r1 = int(np.floor(bb.y0)), int(np.ceil(bb.y1))
    if c1 <= c0 or r1 <= r0 or x1 == x0 or y1 == y0:
        return None
    return {
        "x0": x0,
        "x1": x1,
        "y0": y0,
        "y1": y1,
        "bx0": bb.x0,
        "by0": bb.y0,
        "bw": bb.width,
        "bh": bb.height,
        "c0": c0,
        "c1": c1,
        "r0": r0,
        "r1": r1,
    }


def _composite_points(
    px: dict[str, float],
    x: np.ndarray,
    y: np.ndarray,
    colour: tuple,
    alpha: float,
    size_pt2: float,
    dpi: float,
) -> np.ndarray:
    """RGBA image (straight alpha) of a scatter, for the pixel box ``px`` of a panel.

    Replaces Matplotlib drawing each of tens of thousands of markers as a path in each
    of up to 64 panels, which took about 9 of the 10 s of a pair grid. Each point covers
    the pixels of its marker (``s`` is the marker area in points², so the diameter is
    ``sqrt(s)`` points), with the anti-aliased fraction at the rim, and each pixel
    composites the points over it in drawing order, the last on top, as Matplotlib's
    renderer does: a point adds its colour with opacity ``a = alpha * covered
    fraction``, times the transparency ``1 - a`` of each point drawn after it there.
    Points are snapped to pixel centres.

    ``colour`` is ``("one", rgb)``, ``("lut", codes, lut)`` (a palette index per point)
    or ``("rgb", rgb_per_point)``. Points with a non-finite position or colour are not
    drawn (Matplotlib's "bad" colour). Pure NumPy, so panels can run in threads.
    """
    width, height = int(px["c1"] - px["c0"]), int(px["r1"] - px["r0"])
    n_pix = width * height
    kind = colour[0]
    ok = np.isfinite(x) & np.isfinite(y)
    if kind == "rgb":
        ok &= np.all(np.isfinite(colour[1]), axis=1)
    img = np.zeros((n_pix, 4), dtype=np.float32)
    n = int(np.count_nonzero(ok))
    if n == 0:
        return img.reshape(height, width, 4)
    sx = px["bw"] / (px["x1"] - px["x0"])
    sy = px["bh"] / (px["y1"] - px["y0"])
    col = np.floor(px["bx0"] + (x[ok] - px["x0"]) * sx).astype(np.intp) - int(px["c0"])
    row = (int(px["r1"]) - 1) - np.floor(px["by0"] + (y[ok] - px["y0"]) * sy).astype(
        np.intp
    )
    dx, dy, frac = _marker_footprint(round(float(np.sqrt(size_pt2)) * dpi / 72.0, 3))
    n_off = dx.size
    cols = (col[:, None] + dx[None, :]).ravel()
    rows = (row[:, None] - dy[None, :]).ravel()
    inside = (cols >= 0) & (cols < width) & (rows >= 0) & (rows < height)
    pix = (rows * width + cols)[inside]
    # Entries go point after point, i.e. in drawing order; each is one marker pixel.
    offset = np.tile(np.arange(n_off, dtype=np.int16), n)[inside]
    opacity_of = np.minimum(alpha * frac, 1.0 - 1e-12)
    log_keep_of = np.log1p(-opacity_of)
    cover = 1.0 - np.exp(np.bincount(pix, weights=log_keep_of[offset], minlength=n_pix))
    img[:, 3] = cover
    if kind == "one":
        img[:, :3] = colour[1]
        return img.reshape(height, width, 4)
    point = np.repeat(np.arange(n, dtype=np.int32), n_off)[inside]
    # Sorting pixel * m + entry groups entries by pixel and keeps their drawing order
    # (several times faster than a stable argsort; the key carries both parts back).
    m = pix.size
    key = pix.astype(np.int64) * m + np.arange(m, dtype=np.int64)
    key.sort()
    pix_sorted, entry = np.divmod(key, m)
    off_sorted = offset[entry]
    log_keep = np.cumsum(log_keep_of[off_sorted])
    group_last = np.cumsum(np.bincount(pix_sorted, minlength=n_pix)) - 1
    # Transparency of everything drawn after each entry on its pixel.
    weight = opacity_of[off_sorted] * np.exp(
        log_keep[group_last[pix_sorted]] - log_keep
    )
    point_sorted = point[entry]
    if kind == "lut":
        codes, lut = colour[1][ok], colour[2]
        k = lut.shape[0]
        per_code = np.bincount(
            pix_sorted * k + codes[point_sorted], weights=weight, minlength=n_pix * k
        )
        acc = per_code.reshape(n_pix, k) @ lut
    else:
        rgb = colour[1][ok]
        acc = np.stack(
            [
                np.bincount(
                    pix_sorted,
                    weights=weight * np.ascontiguousarray(rgb[:, c])[point_sorted],
                    minlength=n_pix,
                )
                for c in range(3)
            ],
            axis=1,
        )
    covered = cover > 0
    img[covered, :3] = acc[covered] / cover[covered, None]
    return img.reshape(height, width, 4)


def _add_points_image(ax: Any, px: dict[str, float], img: np.ndarray) -> None:
    """Place ``img`` (from :func:`_composite_points`) on its panel, pixel for pixel."""
    image = mimage.AxesImage(ax, interpolation="nearest", origin="upper")
    image.set_data(img)
    dx = (px["x1"] - px["x0"]) / px["bw"]
    dy = (px["y1"] - px["y0"]) / px["bh"]
    image.set_extent(
        (
            px["x0"] + (px["c0"] - px["bx0"]) * dx,
            px["x0"] + (px["c1"] - px["bx0"]) * dx,
            px["y0"] + (px["r0"] - px["by0"]) * dy,
            px["y0"] + (px["r1"] - px["by0"]) * dy,
        )
    )
    image.set_zorder(1)  # where the scatter collection used to be
    ax.add_image(image)
    ax.set_xlim(px["x0"], px["x1"])
    ax.set_ylim(px["y0"], px["y1"])


def _colormap_rgb(cmap_name: Any, values: np.ndarray, vmin: float, vmax: float):
    """Per-value RGB ``(n, 3)`` as ``scatter(c=values, cmap, vmin, vmax)`` would draw.

    Non-finite values come out as NaN, which :func:`_composite_points` skips.
    """
    cmap = matplotlib.colormaps[cmap_name] if isinstance(cmap_name, str) else cmap_name
    rgb = cmap(mcolors.Normalize(vmin=vmin, vmax=vmax)(values))[:, :3]
    rgb[~np.isfinite(values)] = np.nan
    return rgb


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
    workers: int = 1,
) -> tuple[bytes, list[dict[str, float]]]:
    """z_dim × z_dim Matplotlib grid (square cells).

    Lower panels: latent pairs coloured by ``lower_color_col``; diagonal: the 2-D
    embedding coloured by that latent dimension; upper: latent pairs in one colour, or
    hexbin densities that match ``sns.jointplot(..., kind="hex")``. Points are drawn by
    :func:`_composite_points`, one panel per thread on up to ``workers`` threads
    (the result does not depend on ``workers``).

    Returns PNG bytes and per-cell axis bboxes in figure coordinates
    (for HTML overlay alignment).
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
    if _lower_color_series_is_discrete(raw_color):
        if not bool(full_df[lower_color_col].notna().any()):
            raise ValueError(
                f"Lower-triangle color column `{lower_color_col}` has no values."
            )
        # Stable colours keyed like the legend API, from the unfiltered table: never
        # re-factorize the filtered subset (that would remap its categories).
        point_hex, _ = _scatter_discrete_marker_arrays(
            full_df, df, lower_color_col, discrete_label_colors
        )
        codes, hexes = pd.factorize(pd.Series(point_hex, dtype=object))
        lut = np.array([mcolors.to_rgb(h) for h in hexes], dtype=np.float64)
        lower_colour: tuple = ("lut", codes.astype(np.int64), lut)
    else:
        if bool(pd.to_numeric(raw_color, errors="coerce").isna().all()):
            raise ValueError(
                f"Lower-triangle color column `{lower_color_col}` "
                f"has no numeric values."
            )
        _, cmin, cmax = _continuous_series_stats(full_df[lower_color_col])
        color_num = cast(pd.Series, pd.to_numeric(raw_color, errors="coerce"))
        cvals = np.asarray(color_num, dtype=np.float64)
        cvals_plot = np.where(np.isfinite(cvals), cvals, 0.5 * (cmin + cmax))
        lower_colour = ("rgb", _colormap_rgb(mpl_cmap, cvals_plot, cmin, cmax))

    emb = (diagonal_emb or "pc").lower()
    if emb == "umap":
        if exp.umap is None or "UMAP1" not in df.columns or "UMAP2" not in df.columns:
            raise ValueError("UMAP embedding is not available for this run.")
        emb_x_col, emb_y_col = "UMAP1", "UMAP2"
    elif emb == "pc":
        if "PC1" not in df.columns or "PC2" not in df.columns:
            raise ValueError("PCA components PC1/PC2 are not available.")
        emb_x_col, emb_y_col = "PC1", "PC2"
    else:
        raise ValueError("diagonal_emb must be 'pc' or 'umap'.")
    emb_x = df[emb_x_col].to_numpy(dtype=np.float64)
    emb_y = df[emb_y_col].to_numpy(dtype=np.float64)
    emb_x_lim = _pair_grid_axis_lim_padded(full_df[emb_x_col])
    emb_y_lim = _pair_grid_axis_lim_padded(full_df[emb_y_col])

    upper = (upper_style or "scatter").lower()
    if upper not in ("scatter", "hex"):
        raise ValueError("upper_style must be 'scatter' or 'hex'.")

    def cell_limits(i: int, j: int):
        # Frozen to the unfiltered table, so legend selections only subset points.
        if i == j:
            return _pair_grid_square_xy_limits(emb_x_lim, emb_y_lim)
        return _pair_grid_square_xy_limits(z_axis_lim[j], z_axis_lim[i])

    z_values = [df[zc].to_numpy(dtype=np.float64) for zc in zcols]
    upper_rgb = np.asarray(mcolors.to_rgb(_UPPER_SCATTER_BLUE), dtype=np.float64)

    with ezlab_matplotlib_rc():
        hex_cmap = _pair_jointplot_hex_cmap()
        geom = _pair_grid_geometry(zdim)
        fig, axes = _new_pair_grid_figure(zdim, cell_limits)
        # Points are drawn at output resolution: lay the figure out at that dpi.
        fig.set_dpi(dpi)

        diagonal_color_ranges: list[tuple[float, float] | None] = [None] * zdim
        layers: list[tuple[Any, dict[str, float], tuple]] = []
        for i in range(zdim):
            for j in range(zdim):
                ax = axes[i, j]
                xi, yi = z_values[j], z_values[i]
                if i == j:
                    zi_full = full_df[zcols[i]].to_numpy(dtype=np.float64)
                    zm_full = np.isfinite(zi_full)
                    if not np.any(zm_full):
                        continue
                    zmin_g = float(np.min(zi_full[zm_full]))
                    zmax_g = float(np.max(zi_full[zm_full]))
                    if zmax_g <= zmin_g:
                        zmax_g = zmin_g + 1e-9
                    diagonal_color_ranges[i] = (zmin_g, zmax_g)
                    diag_rgb = _colormap_rgb(mpl_cmap, z_values[i], zmin_g, zmax_g)
                    points = (emb_x, emb_y, ("rgb", diag_rgb), diag_a, diag_s)
                elif i < j and upper == "hex":
                    mxy = np.isfinite(xi) & np.isfinite(yi)
                    xh, yh = xi[mxy], yi[mxy]
                    if xh.size:
                        try:
                            gs = _pair_jointplot_hex_gridsize(xh, yh)
                            ax.hexbin(xh, yh, gridsize=gs, cmap=hex_cmap)
                        except ZeroDivisionError:
                            gs = int(np.clip(22 + np.sqrt(max(n_pts, 1)), 30, 75))
                            ax.hexbin(xh, yh, gridsize=gs, cmap=hex_cmap, mincnt=1)
                        xl, yl = cell_limits(i, j)
                        ax.set_xlim(xl[0], xl[1])
                        ax.set_ylim(yl[0], yl[1])
                    continue
                elif i < j:
                    points = (xi, yi, ("one", upper_rgb), up_scatter_a, up_scatter_s)
                else:
                    points = (xi, yi, lower_colour, lo_a, lo_s)
                px = _panel_pixels(ax)
                if px is not None:
                    layers.append((ax, px, points))

        def composite(layer: tuple) -> np.ndarray:
            _, px, (x, y, colour, alpha, size) = layer
            return _composite_points(px, x, y, colour, alpha, size, dpi)

        if workers > 1 and len(layers) > 1:
            with ThreadPoolExecutor(max_workers=min(workers, len(layers))) as pool:
                images = list(pool.map(composite, layers))
        else:
            images = [composite(layer) for layer in layers]
        for (ax, px, _), img in zip(layers, images):
            _add_points_image(ax, px, img)

        label_fs = max(11.0, min(16.0, 4.9 * geom["inch_per"] / 3.05))
        _draw_pair_grid_edge_labels(
            fig,
            axes,
            zdim=zdim,
            left_m=geom["left"],
            top_m=geom["top"],
            bottom_m=geom["bottom"],
            right_axes=geom["right"],
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
        png = _opaque_figure_png(fig)
        plt.close(fig)
        return png, cell_layout


def _opaque_figure_png(fig: Any, compress_level: int = 3) -> bytes:
    """The whole figure (no tight crop, so figure fractions match PNG pixels for the
    overlay) as an RGB PNG at the figure's dpi.

    Same pixels as ``savefig(format="png")`` with the figure's own facecolor, without
    its always-opaque alpha channel; ``compress_level`` 3 instead of Pillow's 6 saves
    0.2 s per grid for 10% more bytes.
    """
    fig.canvas.draw()
    rgb = np.asarray(fig.canvas.buffer_rgba())[..., :3]
    buf = io.BytesIO()
    Image.fromarray(np.ascontiguousarray(rgb)).save(
        buf, format="png", compress_level=compress_level
    )
    return buf.getvalue()


# Recent pair grids by request: switching back and "Save PNG" reuse them instead of
# rendering again. Columns of ``plot_df`` never change once loaded (covariates loaded
# later get new names), so the experiment object and the request identify the image.
_PAIR_GRID_CACHE: "OrderedDict[tuple, tuple[bytes, list[dict[str, float]]]]" = (
    OrderedDict()
)
_PAIR_GRID_CACHE_LOCK = threading.Lock()
_PAIR_GRID_CACHE_SIZE = 8


def pair_grid_png_cached(
    exp: DashboardExperiment, **kwargs: Any
) -> tuple[bytes, list[dict[str, float]]]:
    """:func:`pair_grid_png`, remembering the last few results (about 2-4 MB each)."""
    # ``workers`` changes how fast, not what: leave it out of the key.
    request = {k: v for k, v in kwargs.items() if k != "workers"}
    key = (
        id(exp),
        exp.workdir,
        int(exp.epoch),
        json.dumps(request, sort_keys=True, default=str),
    )
    with _PAIR_GRID_CACHE_LOCK:
        hit = _PAIR_GRID_CACHE.get(key)
        if hit is not None:
            _PAIR_GRID_CACHE.move_to_end(key)
            return hit[0], [dict(c) for c in hit[1]]
    png, cells = pair_grid_png(exp, **kwargs)
    with _PAIR_GRID_CACHE_LOCK:
        _PAIR_GRID_CACHE[key] = (png, [dict(c) for c in cells])
        while len(_PAIR_GRID_CACHE) > _PAIR_GRID_CACHE_SIZE:
            _PAIR_GRID_CACHE.popitem(last=False)
    return png, cells
