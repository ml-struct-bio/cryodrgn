"""Latent pair-grid Matplotlib figures (PNG) and layout helpers."""

from __future__ import annotations

import io
import re
from typing import Any, cast

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
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
    covariate_row_filter_key,
    plot_df_color_filter_mask,
)
from cryodrgn.dashboard.plots_figure_utils import (
    _UPPER_SCATTER_BLUE,
    _axes_cell_bboxes,
    _dashboard_chimerax_colors,
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
_PAIR_GRID_RIGHT_PLACEHOLDER = _PAIR_GRID_RIGHT


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
    """z_dim × z_dim Matplotlib grid (square cells).

    Upper hex matches ``sns.jointplot(..., kind=\"hex\")``.

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
                f"Lower-triangle color column `{lower_color_col}` "
                f"has no numeric values."
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
                # Freeze axes to the unfiltered table so legend selections
                # only subset points.
                if i == j:
                    xl, yl = _pair_grid_square_xy_limits(emb_x_lim, emb_y_lim)
                    ax.set_xlim(xl[0], xl[1])
                    ax.set_ylim(yl[0], yl[1])
                else:
                    xl, yl = _pair_grid_square_xy_limits(z_axis_lim[j], z_axis_lim[i])
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
        # Full figure bbox (no tight crop) so figure fractions match
        # PNG pixels for the overlay.
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
