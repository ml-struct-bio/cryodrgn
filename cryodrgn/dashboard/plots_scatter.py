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

from cryodrgn.dashboard.column_names import (
    VOL_LANDSCAPE_3D_PLOT_DF_ROW,
    VOL_LANDSCAPE_IS_SKETCH_CENTROID,
    VOL_LANDSCAPE_NEAREST_SKETCH_VOL,
)
from cryodrgn.dashboard.covariate_labels import (
    covariate_display_name,
    landscape_vol_pc_column_pretty_label,
)
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
    _lower_legend_entry_label,
    _scatter_discrete_marker_arrays,
    _stable_discrete_covariate_hex_map,
    discrete_category_counts_by_filter_key,
    numeric_array_to_plotly_hex,
    plot_df_color_filter_mask,
)
from cryodrgn.dashboard.plots_figure_utils import (
    _DASHBOARD_CREAM,
    _PLOTLY_FONT,
    DASHBOARD_SCATTER_HOVERLABEL_FONT_SIZE,
    _pair_grid_axis_lim_padded,
    _scatter3d_marker_size_opacity,
    _subsample,
    _UMAP_RE,
    ezlab_matplotlib_rc,
    _plotly_to_json,
)


def _dashboard_scatter3d_glyph_visual_scale(
    size: float, opacity: float
) -> tuple[float, float]:
    """In-dashboard Plotly scatter3d: ~30% larger glyphs, ~20% lower alpha (more transparent)."""
    return (
        float(size * 1.3),
        float(max(0.06, min(0.98, opacity * 0.8))),
    )


def _scatter_color_hover_texts(
    sub: pd.DataFrame,
    color_col: str,
    *,
    discrete_trace: bool,
    vol_pc_explained_variance_ratio: np.ndarray | None = None,
) -> list[str]:
    """Per-point second line for scatter hover (covariate text; continuous uses ``name: value``)."""
    ser = sub[color_col]
    n = len(sub)
    out: list[str] = []
    if discrete_trace:
        for i in range(n):
            v = ser.iloc[i]
            if pd.isna(v):
                out.append("(missing)")
            elif color_col == "labels":
                out.append(f"kmeans={_lower_legend_entry_label('labels', v)}")
            else:
                out.append(str(v))
    else:
        cov_label = landscape_vol_pc_column_pretty_label(
            color_col,
            vol_pc_explained_variance_ratio,
        )
        if cov_label == color_col:
            cov_label = covariate_display_name(color_col)
        color_num = cast(pd.Series, pd.to_numeric(ser, errors="coerce"))
        for i in range(n):
            v = color_num.iloc[i]
            if pd.isna(v):
                out.append(f"{cov_label}: —")
            else:
                fv = float(v)
                if np.isfinite(fv):
                    out.append(f"{cov_label}: {fv:.5g}")
                else:
                    out.append(f"{cov_label}: —")
    return out


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

    has_cov_color = bool(color_col and color_col != "none" and color_col in sub.columns)
    # Scattergl does not reliably show ``hovertext`` in ``hovertemplate``;
    # put the covariate line in ``customdata`` so the second line renders in WebGL too.
    if has_cov_color:
        cc = cast(str, color_col)
        hov2 = np.asarray(
            _scatter_color_hover_texts(sub, cc, discrete_trace=bool(discrete_trace)),
            dtype=object,
        ).reshape(-1, 1)
        customdata = np.column_stack([customdata, hov2])
        hover_kwargs: dict[str, Any] = dict(
            hovertemplate="particle: %{customdata[0]}<br>%{customdata[3]}<extra></extra>",
        )
    else:
        hover_kwargs = dict(hovertemplate="particle: %{customdata[0]}<extra></extra>")

    # Scattergl can leave Plotly.react() pending on some GPUs; SVG Scatter is unusably slow
    # at dashboard point caps (~1e5+). The explorer client uses ``plotly_afterplot`` plus a timeout.
    trace_cls = go.Scattergl if use_webgl else go.Scatter
    sc = trace_cls(
        x=sub[xcol],
        y=sub[ycol],
        mode="markers",
        customdata=customdata,
        marker=marker,
        **hover_kwargs,
    )

    xaxis_kw: dict[str, Any] = dict(title=xcol)
    yaxis_kw: dict[str, Any] = dict(title=ycol)
    if _UMAP_RE.match(xcol):
        xaxis_kw["showticklabels"] = False
        xaxis_kw["showgrid"] = False
    if _UMAP_RE.match(ycol):
        yaxis_kw["showticklabels"] = False
        yaxis_kw["showgrid"] = False

    hoverlabel_font = dict(_PLOTLY_FONT)
    hoverlabel_font["size"] = DASHBOARD_SCATTER_HOVERLABEL_FONT_SIZE
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
        hoverlabel=dict(font=hoverlabel_font, align="left"),
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
    no_subsample: bool = False,
    plot_df: pd.DataFrame | None = None,
    continuous_palette: str | None = None,
    color_filter: dict[str, Any] | None = None,
    discrete_label_colors: dict[str, str] | None = None,
    uirevision: str = "scatter3d_z",
    xyz_axes_allowed: frozenset[str] | None = None,
    scene_axis_titles: tuple[str, str, str] | None = None,
    vol_pc_explained_variance_ratio: np.ndarray | None = None,
    volume_landscape_3d_style: bool = False,
) -> str:
    """Interactive 3D scatter of three latent ``z*`` columns.

    Pass ``plot_df`` to use an alternate table (e.g. ``analyze_landscape_full`` sampled rows)
    while still validating axes against ``exp.z`` shape.

    When ``xyz_axes_allowed`` is set (e.g. volume PCA column names), ``xcol``/``ycol``/``zcol``
    must be three distinct members of that set instead of latent ``z*`` dimensions.

    ``scene_axis_titles`` overrides Plotly scene axis titles (defaults to ``xcol``/``ycol``/``zcol``).

    ``vol_pc_explained_variance_ratio`` improves hover and colour-legend titles for
    ``landscape_vol_PC*`` covariates when present.

    ``volume_landscape_3d_style`` tightens marker size (vs the usual n-points curve, then two
    further ~13% steps for this 3D UI), emboldens axis titles, and bumps scene tick label size for the
    full-volume-landscape 3D view.

    When ``no_subsample`` is true, every row in the (possibly colour-filtered) table is plotted
    instead of capping at ``max_points`` — used for discrete-level GIF frames so rare categories
    are not randomly thinned away.
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

    if no_subsample:
        sub, row_indices = _subsample(df, None, seed=1)
    elif VOL_LANDSCAPE_IS_SKETCH_CENTROID in df.columns:
        from cryodrgn.dashboard.plots_figure_utils import (
            _subsample_preserving_sketch_centroids,
        )

        sub, row_indices = _subsample_preserving_sketch_centroids(
            df, max_points, seed=1, centroid_col=VOL_LANDSCAPE_IS_SKETCH_CENTROID
        )
    else:
        sub, row_indices = _subsample(df, max_points, seed=1)
    cap = max(len(sub), 1) if no_subsample else max_points
    msize, mopacity = _scatter3d_marker_size_opacity(len(sub), point_cap=cap)
    if volume_landscape_3d_style:
        # ~31% smaller vs the usual n-points curve, then two ~13% reductions for this 3D UI,
        # then a further 20% shrink for the volume-landscape glyph size.
        msize *= (1.0 - 0.31) * (1.0 - 0.13) * (1.0 - 0.13) * 0.8
    msize, mopacity = _dashboard_scatter3d_glyph_visual_scale(msize, mopacity)
    if volume_landscape_3d_style:
        # ~13% larger points, then +73% radius vs pre-landscape baseline (three +20% steps); ~19% more transparent.
        msize *= 1.13 * 1.728 * 1.11
        mopacity = float(max(0.0, min(1.0, mopacity * 0.81)))

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
            color_legend_title = landscape_vol_pc_column_pretty_label(
                color_col,
                vol_pc_explained_variance_ratio,
            )
            if color_legend_title == color_col:
                color_legend_title = covariate_display_name(color_col)
            legend_meta = {
                "type": "discrete",
                "title": color_legend_title,
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
            if pinned_color_range is not None:
                cmin, cmax = pinned_color_range
            hex_colors = numeric_array_to_plotly_hex(
                _cvals, plotly_cs, vmin=cmin, vmax=cmax
            )
            # Per-point hex matches Plotly's colorscale exactly (WebGL + letter overlays).
            marker = dict(
                size=msize,
                opacity=mopacity,
                color=hex_colors,
            )
            color_legend_title = landscape_vol_pc_column_pretty_label(
                color_col,
                vol_pc_explained_variance_ratio,
            )
            if color_legend_title == color_col:
                color_legend_title = covariate_display_name(color_col)
            legend_meta = {
                "type": "continuous",
                "title": color_legend_title,
                "min": cmin,
                "max": cmax,
            }
    else:
        marker = dict(size=msize, opacity=mopacity, color="#4a5568")

    # Plotly Scatter3d otherwise picks up a default marker outline in WebGL.
    marker["line"] = dict(width=0)

    idx_cd = sub["index"].to_numpy()
    if VOL_LANDSCAPE_3D_PLOT_DF_ROW in sub.columns:
        row_cd = np.asarray(sub[VOL_LANDSCAPE_3D_PLOT_DF_ROW], dtype=np.int64)
    else:
        row_cd = np.asarray(row_indices, dtype=np.int64)

    has_cov_color = bool(color_col and color_col != "none" and color_col in sub.columns)
    if discrete_trace:
        fk_arr = np.asarray(fk_list, dtype=object)
        customdata = np.column_stack([idx_cd, row_cd, fk_arr])
    elif has_cov_color:
        n_sub = len(sub)
        disp_cd = np.empty(n_sub, dtype=object)
        color_num = cast(pd.Series, pd.to_numeric(sub[color_col], errors="coerce"))
        for i in range(n_sub):
            v = color_num.iloc[i]
            if pd.isna(v):
                disp_cd[i] = None
            else:
                fv = float(v)
                disp_cd[i] = fv if np.isfinite(fv) else None
        customdata = np.column_stack([idx_cd, row_cd, disp_cd])
    else:
        customdata = np.column_stack([idx_cd, row_cd])

    if has_cov_color:
        cc = cast(str, color_col)
        hov2_3d = np.asarray(
            _scatter_color_hover_texts(
                sub,
                cc,
                discrete_trace=bool(discrete_trace),
                vol_pc_explained_variance_ratio=vol_pc_explained_variance_ratio,
            ),
            dtype=object,
        ).reshape(-1, 1)
        customdata = np.column_stack([customdata, hov2_3d])
        hover_kwargs_3d: dict[str, Any] = dict(
            hovertemplate="particle: %{customdata[0]}<br>%{customdata[3]}<extra></extra>",
        )
    else:
        hover_kwargs_3d = dict(
            hovertemplate="particle: %{customdata[0]}<extra></extra>",
        )

    if VOL_LANDSCAPE_NEAREST_SKETCH_VOL in sub.columns:
        nv = np.asarray(sub[VOL_LANDSCAPE_NEAREST_SKETCH_VOL], dtype=np.int64)
        if VOL_LANDSCAPE_IS_SKETCH_CENTROID in sub.columns:
            cd_cent = np.asarray(sub[VOL_LANDSCAPE_IS_SKETCH_CENTROID], dtype=np.int64)
            customdata = np.column_stack([customdata, cd_cent])
        customdata = np.column_stack([customdata, nv])

        # Volume-landscape hover: show which k-means volume each particle belongs to.
        # We append `nv` as the last customdata column, so use `customdata[-1]`.
        if volume_landscape_3d_style:
            nv_cd_idx = int(customdata.shape[1] - 1)
            if has_cov_color:
                hover_kwargs_3d["hovertemplate"] = (
                    "particle: %{customdata[0]}"
                    + "<br>%{customdata[3]}"
                    + "<br>k_means_vol=%{customdata["
                    + str(nv_cd_idx)
                    + "]}"
                    + "<extra></extra>"
                )
            else:
                hover_kwargs_3d["hovertemplate"] = (
                    "particle: %{customdata[0]}"
                    + "<br>k_means_vol=%{customdata["
                    + str(nv_cd_idx)
                    + "]}"
                    + "<extra></extra>"
                )

    # Match the volume-selection overlay (``latent3d_landscape_vol_animations.js``): keep
    # ``markers+text`` from the first draw so the first selection only restyles marker/text
    # fields — not ``mode`` — which avoids Plotly resetting orbit, zoom, or axis limits.
    vol_landscape_sel_overlay = VOL_LANDSCAPE_NEAREST_SKETCH_VOL in sub.columns
    n_pts = len(sub)
    scatter3d_mode = "markers+text" if vol_landscape_sel_overlay else "markers"
    scatter3d_extras: dict[str, Any] = {}
    if vol_landscape_sel_overlay:
        scatter3d_extras["text"] = [""] * n_pts
        # Bring montage letters closer to the points (previously "top center").
        scatter3d_extras["textposition"] = "middle center"
        scatter3d_extras["textfont"] = dict(
            size=36 * 0.8 * 0.75 * 1.5 * 0.8 * 0.8 * 0.8 * 0.8,
            color="#1a1a1a",
            family="system-ui, Segoe UI, sans-serif",
        )

    sc = go.Scatter3d(
        x=sub[xcol],
        y=sub[ycol],
        z=sub[zcol],
        mode=scatter3d_mode,
        customdata=customdata,
        marker=marker,
        **scatter3d_extras,
        **hover_kwargs_3d,
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
        if not discrete_trace:
            plot_meta["cdrgn_continuous_palette"] = plotly_cs
    if VOL_LANDSCAPE_NEAREST_SKETCH_VOL in sub.columns:
        plot_meta["cdrgn_landscape_vol_animation"] = True
        if VOL_LANDSCAPE_IS_SKETCH_CENTROID in sub.columns:
            plot_meta["cdrgn_landscape_sketch_centroid_cd"] = True

    hoverlabel_font_3d = dict(_PLOTLY_FONT)
    hoverlabel_font_3d["size"] = DASHBOARD_SCATTER_HOVERLABEL_FONT_SIZE
    x_title, y_title, z_title = (
        scene_axis_titles if scene_axis_titles is not None else (xcol, ycol, zcol)
    )
    fig = go.Figure(sc)
    scene_kw: dict[str, Any] = dict(
        aspectmode="data",
        bgcolor="rgba(250,248,244,0.95)",
        dragmode="orbit",
    )
    if volume_landscape_3d_style:
        # Reference tick size 12px → +33%; layout font → +20% for titles (bold).
        tick_px = 12.0 * (1.0 + 0.33)
        title_px = float(_PLOTLY_FONT["size"]) * (1.0 + 0.20)
        title_font = dict(
            family=_PLOTLY_FONT["family"],
            size=title_px,
            color="#1a202c",
            weight="bold",
        )
        tick_font = dict(
            family=_PLOTLY_FONT["family"],
            size=tick_px,
            color="#4a5568",
        )

        def _scene_axis_vol_land(
            axis_extra: dict[str, Any], title_text: str
        ) -> dict[str, Any]:
            d = dict(axis_extra)
            d["tickfont"] = tick_font
            d["title"] = dict(text=title_text, font=title_font)
            return d

        scene_kw["xaxis"] = _scene_axis_vol_land(scene_axes_extra["xaxis"], x_title)
        scene_kw["yaxis"] = _scene_axis_vol_land(scene_axes_extra["yaxis"], y_title)
        scene_kw["zaxis"] = _scene_axis_vol_land(scene_axes_extra["zaxis"], z_title)
    else:
        scene_kw["xaxis_title"] = x_title
        scene_kw["yaxis_title"] = y_title
        scene_kw["zaxis_title"] = z_title
        scene_kw.update(scene_axes_extra)
    fig.update_layout(
        template="plotly_white",
        autosize=True,
        paper_bgcolor=_DASHBOARD_CREAM,
        margin=dict(l=0, r=0, t=0, b=0),
        scene=scene_kw,
        uirevision=uirevision,
        font=_PLOTLY_FONT,
        meta=plot_meta,
        hoverlabel=dict(font=hoverlabel_font_3d, align="left"),
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


def scatter3d_discrete_level_png_bytes(
    plot_df: pd.DataFrame,
    xcol: str,
    ycol: str,
    zcol: str,
    color_col: str,
    filter_key: str,
    *,
    discrete_label_colors: dict[str, str] | None,
    scene_axis_titles: tuple[str, str, str],
    exp: DashboardExperiment | None,
    xyz_axes_allowed: frozenset[str] | None,
    volume_landscape_marker_shrink: bool = False,
    elev: float = 22.0,
    azim: float = -65.0,
    dpi: int = 100,
) -> bytes:
    """One Matplotlib 3D frame: every row in a single discrete category (no subsampling).

    Used for server-side discrete-level GIFs. Axis limits use padded min/max on the full
    ``plot_df`` (same idea as the colour-filtered Plotly scene).

    Pass ``xyz_axes_allowed`` for volume-PCA landscape tables (with ``exp=None``), or
    ``exp`` with ``xyz_axes_allowed=None`` for standard ``z*`` latent axes.
    """
    if xyz_axes_allowed is not None:
        _validate_three_xyz_allowed(plot_df, (xcol, ycol, zcol), xyz_axes_allowed)
    else:
        if exp is None:
            raise ValueError(
                "exp is required for latent-axis discrete PNG export when xyz_axes_allowed is unset."
            )
        _validate_three_latent_axes(exp, plot_df, (xcol, ycol, zcol))
    if not color_col or color_col == "none" or color_col not in plot_df.columns:
        raise ValueError("A discrete colour column is required for this export.")
    df_all = plot_df
    if not _lower_color_series_is_discrete(df_all[color_col]):
        raise ValueError(
            "Colour column must be discrete for discrete-level GIF export."
        )

    color_filter: dict[str, Any] = {"kind": "discrete", "keys": [filter_key]}
    mask = plot_df_color_filter_mask(df_all, color_col, color_filter)
    if not bool(np.any(mask)):
        raise ValueError(f"No particles match discrete level {filter_key!r}.")

    df = df_all.loc[mask].reset_index(drop=True)
    sub, _ = _subsample(df, None, seed=1)

    lookup = _stable_discrete_covariate_hex_map(
        df_all, color_col, discrete_label_colors
    )
    if filter_key not in lookup:
        raise ValueError(f"Unknown discrete level {filter_key!r}.")
    hex_c = lookup[filter_key]

    cap = max(len(df), 1)
    msize, mopacity = _scatter3d_marker_size_opacity(len(sub), point_cap=cap)
    if volume_landscape_marker_shrink:
        msize *= 1.0 - 0.31
    msize, mopacity = _dashboard_scatter3d_glyph_visual_scale(msize, mopacity)
    mpl_marker_area = float(max(6.0, (msize * 2.2) ** 2))

    xs = sub[xcol].to_numpy(dtype=np.float64)
    ys = sub[ycol].to_numpy(dtype=np.float64)
    zs = sub[zcol].to_numpy(dtype=np.float64)

    xr = list(_pair_grid_axis_lim_padded(df_all[xcol]))
    yr = list(_pair_grid_axis_lim_padded(df_all[ycol]))
    zr = list(_pair_grid_axis_lim_padded(df_all[zcol]))

    x_title, y_title, z_title = scene_axis_titles

    with ezlab_matplotlib_rc():
        fig = plt.figure(figsize=(5.2, 4.6), facecolor=_DASHBOARD_CREAM)
        ax = fig.add_subplot(111, projection="3d", facecolor=_DASHBOARD_CREAM)
        ax.scatter(
            xs,
            ys,
            zs,
            c=hex_c,
            s=mpl_marker_area,
            alpha=float(mopacity),
            linewidths=0,
            depthshade=True,
        )
        ax.set_xlim(xr)
        ax.set_ylim(yr)
        ax.set_zlim(zr)
        ax.set_xlabel(x_title, fontsize=10)
        ax.set_ylabel(y_title, fontsize=10)
        ax.set_zlabel(z_title, fontsize=10)
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


def scatter3d_landscape_full_discrete_level_png_bytes(
    plot_df: pd.DataFrame,
    xcol: str,
    ycol: str,
    zcol: str,
    color_col: str,
    filter_key: str,
    *,
    discrete_label_colors: dict[str, str] | None,
    xyz_axes_allowed: frozenset[str],
    scene_axis_titles: tuple[str, str, str],
    elev: float = 22.0,
    azim: float = -65.0,
    dpi: int = 100,
) -> bytes:
    """Volume-landscape discrete frame (smaller markers); prefer :func:`scatter3d_discrete_level_png_bytes`."""
    return scatter3d_discrete_level_png_bytes(
        plot_df,
        xcol,
        ycol,
        zcol,
        color_col,
        filter_key,
        discrete_label_colors=discrete_label_colors,
        scene_axis_titles=scene_axis_titles,
        exp=None,
        xyz_axes_allowed=xyz_axes_allowed,
        volume_landscape_marker_shrink=True,
        elev=elev,
        azim=azim,
        dpi=dpi,
    )
