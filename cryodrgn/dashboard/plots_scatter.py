"""Plotly scatter figures (2-D explorer and 3-D latent / volume-PCA views)."""

from __future__ import annotations

from collections.abc import Collection
from typing import Any, cast

import io
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from cryodrgn.dashboard.column_names import VOL_LANDSCAPE_3D_PLOT_DF_ROW
from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.palette_config import (
    mpl_cmap_for_palette,
    normalize_continuous_palette,
)
from cryodrgn.dashboard.plots_color_covariate import (
    _continuous_series_stats,
    _discrete_legend_sort_tuple,
    _labels_colors_and_legend_items,
    _lower_color_series_is_discrete,
    _scatter_discrete_marker_arrays,
    _stable_discrete_covariate_hex_map,
    discrete_category_counts_by_filter_key,
    plot_df_color_filter_mask,
)
from cryodrgn.dashboard.plots_figure_utils import (
    _DASHBOARD_CREAM,
    _PLOTLY_FONT,
    _pair_grid_axis_lim_padded,
    _scatter3d_marker_size_opacity,
    _subsample,
    _UMAP_RE,
    ezlab_matplotlib_rc,
    _plotly_to_json,
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

    idx_arr = sub["index"].to_numpy()
    row_arr = np.asarray(row_indices, dtype=np.int64)
    n_pts = len(sub)

    discrete_trace = (
        color_col
        and color_col != "none"
        and color_col in exp.plot_df.columns
        and color_col in sub.columns
        and _lower_color_series_is_discrete(exp.plot_df[color_col])
    )

    if discrete_trace:
        hex_colors, fk_list = _scatter_discrete_marker_arrays(
            exp.plot_df,
            sub,
            color_col,
            None,
        )
        marker = dict(
            size=marker_size,
            opacity=0.35,
            color=hex_colors,
        )
        fk_arr = np.asarray(fk_list, dtype=object)
        customdata = np.column_stack([idx_arr, row_arr, fk_arr])
    elif color_col and color_col != "none" and color_col in sub.columns:
        marker = dict(
            size=marker_size,
            opacity=0.35,
            color=sub[color_col],
            colorscale=plotly_cs,
        )
        disp = np.empty(n_pts, dtype=object)
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
        marker = dict(size=marker_size, opacity=0.35, color="#4a5568")
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
        layout_meta["cdrgn_color_mode"] = "discrete" if discrete_trace else "continuous"
        if discrete_trace and color_col:
            layout_meta[
                "cdrgn_discrete_category_counts"
            ] = discrete_category_counts_by_filter_key(exp.plot_df, color_col)
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


def _validate_three_xyz_allowed(
    df: pd.DataFrame,
    axes: tuple[str, str, str],
    allowed: frozenset[str],
) -> None:
    """Validate three distinct axes, each in ``allowed`` and present in ``df``."""
    for col in axes:
        if col not in allowed:
            raise ValueError(f"Axis {col!r} is not allowed for this plot.")
        if col not in df.columns:
            raise ValueError(f"Missing column {col!r} in analysis table.")
    if len(set(axes)) < 3:
        raise ValueError("Choose three distinct axes.")


def scatter3d_z_json(
    exp: DashboardExperiment,
    xcol: str,
    ycol: str,
    zcol: str,
    color_col: str | None,
    max_points: int = 120_000,
    *,
    plot_df: pd.DataFrame | None = None,
    continuous_palette: str | None = None,
    color_filter: dict[str, Any] | None = None,
    discrete_label_colors: dict[str, str] | None = None,
    uirevision: str = "scatter3d_z",
    xyz_axes_allowed: frozenset[str] | None = None,
) -> str:
    """Interactive 3D scatter of three latent ``z*`` columns.

    Pass ``plot_df`` to use an alternate table (e.g. ``analyze_landscape_full`` sampled rows)
    while still validating axes against ``exp.z`` shape.

    When ``xyz_axes_allowed`` is set (e.g. volume PCA column names), ``xcol``/``ycol``/``zcol``
    must be three distinct members of that set instead of latent ``z*`` dimensions.
    """
    plotly_cs = normalize_continuous_palette(continuous_palette)
    df_all = exp.plot_df if plot_df is None else plot_df
    df = df_all
    if xyz_axes_allowed is not None:
        _validate_three_xyz_allowed(df_all, (xcol, ycol, zcol), xyz_axes_allowed)
    else:
        _validate_three_latent_axes(exp, df_all, (xcol, ycol, zcol))

    pinned_color_range: tuple[float, float] | None = None
    if (
        color_filter
        and color_col
        and color_col != "none"
        and color_col in df.columns
        and not _lower_color_series_is_discrete(df_all[color_col])
    ):
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

    discrete_trace = (
        color_col
        and color_col != "none"
        and color_col in df_all.columns
        and color_col in sub.columns
        and _lower_color_series_is_discrete(df_all[color_col])
    )

    legend_meta: dict[str, Any] | None = None
    if color_col and color_col != "none" and color_col in sub.columns:
        if discrete_trace:
            hex_colors, fk_list = _scatter_discrete_marker_arrays(
                df_all,
                sub,
                color_col,
                discrete_label_colors,
            )
            marker = dict(
                size=msize,
                opacity=mopacity,
                color=hex_colors,
            )
            lookup_preview = _stable_discrete_covariate_hex_map(
                df_all, color_col, discrete_label_colors
            )
            counts_map = discrete_category_counts_by_filter_key(df_all, color_col)
            sort_keys = sorted(lookup_preview.keys(), key=_discrete_legend_sort_tuple)
            legend_meta = {
                "type": "discrete",
                "title": color_col,
                "items": [
                    {
                        "label": k,
                        "color": lookup_preview[k],
                        "count": counts_map.get(k, 0),
                    }
                    for k in sort_keys
                ],
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
    if VOL_LANDSCAPE_3D_PLOT_DF_ROW in sub.columns:
        row_cd = np.asarray(sub[VOL_LANDSCAPE_3D_PLOT_DF_ROW], dtype=np.int64)
    else:
        row_cd = np.asarray(row_indices, dtype=np.int64)
    if discrete_trace:
        fk_arr = np.asarray(fk_list, dtype=object)
        customdata = np.column_stack([idx_cd, row_cd, fk_arr])
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

    plot_meta: dict[str, Any] = {}
    if legend_meta is not None:
        plot_meta["cdrgn_color_legend"] = legend_meta
    if discrete_trace and color_col and color_col != "none":
        plot_meta[
            "cdrgn_discrete_category_counts"
        ] = discrete_category_counts_by_filter_key(df_all, color_col)
    if color_col and color_col != "none" and color_col in sub.columns:
        plot_meta["cdrgn_color_mode"] = "discrete" if discrete_trace else "continuous"

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
        uirevision=uirevision,
        font=_PLOTLY_FONT,
        meta=plot_meta,
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
    plot_df: pd.DataFrame | None = None,
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
    df = exp.plot_df if plot_df is None else plot_df
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
