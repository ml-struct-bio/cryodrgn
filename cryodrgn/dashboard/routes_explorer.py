"""HTTP routes: landing, explorer, latent 3-D, pair-grid-adjacent APIs through preload."""

from __future__ import annotations

import base64
import logging
import os
import pickle
import time
from concurrent.futures import ThreadPoolExecutor
from typing import cast

import numpy as np
import pandas as pd
from flask import (
    Response,
    current_app,
    g,
    jsonify,
    render_template,
    request,
    url_for,
)

from cryodrgn.dashboard.context import (
    PRELOAD_CACHE,
    active_workdir,
    clear_preload_cache_for_experiment,
    command_builder_template_kwargs,
    _request_json_dict,
)
from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.covariate_labels import landscape_vol_pc_column_pretty_label
from cryodrgn.dashboard.landscape_full_3d import (
    LANDSCAPE_FULL_3D_LEAD_HTML,
    LANDSCAPE_FULL_3D_LEGEND_CONTEXT_EXTRA,
    landscape_full_3d_ready,
    landscape_full_ready,
    landscape_full_sampled_numeric_covariates,
    landscape_full_sampled_plot_df,
    landscape_full_vol_pca_axis_columns,
)
from cryodrgn.dashboard.landscape_volpca import (
    landscape_analysis_ready,
    landscape_dir_for_epoch,
    load_pca_explained_variance,
)
from cryodrgn.dashboard.plot_gif_utils import png_base64_frames_to_gif_bytes
from cryodrgn.dashboard.particle_explorer import (
    DEFAULT_CHIMERAX_PARALLEL,
    DEFAULT_GIF_FRAMES,
    explorer_volumes_eligible,
    generate_montage_volume_pngs,
    volume_cell_gif_from_cache,
)
from cryodrgn.dashboard.plots import (
    covariate_legend_context_payload,
    plot_df_color_filter_mask,
    plot_df_row_indices_for_explorer_scatter,
    scatter3d_landscape_full_discrete_level_png_bytes,
    scatter3d_z_json,
    scatter3d_z_preview_png,
    scatter_json,
)
from cryodrgn.dashboard.plots_color_covariate import _lower_color_series_is_discrete
from cryodrgn.dashboard.preload import (
    encode_particle_batch,
    explorer_cache_size_power10_step,
    explorer_initial_preload_image_limit,
    format_preload_cache_time_hint,
    load_plot_df_rows_from_plot_inds_file,
    montage_bytes,
    sample_plot_df_rows_for_preload,
)
from cryodrgn.dashboard.route_helpers import (
    _EXPLORER_VOLUMES_INELIGIBLE_MSG,
    _covariate_display_map,
    _default_xy_cols,
    _filter_ui_scatter_max_points,
    _parse_color_filter_for_column,
    _parse_optional_discrete_label_colors,
    _parse_preload_image_limit,
    _parse_preselect_rows_param,
    _particle_explorer_scatter_cap_from_env,
    _particle_explorer_scatter_max_points,
    _redirect,
)

logger = logging.getLogger(__name__)


def _scatter3d_no_subsample_for_discrete_gif_frame(
    payload: dict, color_filter: dict | None
) -> bool:
    """Return True when the client is capturing one discrete level at a time for a GIF."""
    if not payload.get("discrete_level_gif_frame"):
        return False
    if not color_filter or color_filter.get("kind") != "discrete":
        return False
    keys = color_filter.get("keys") or ()
    return len(tuple(keys)) == 1


def _landscape_full_vol_pc_explained_variance(
    workdir: str, epoch: int
) -> np.ndarray | None:
    """Explained variance ratios from ``landscape.{epoch}/vol_pca_obj.pkl``, if present."""
    return load_pca_explained_variance(landscape_dir_for_epoch(workdir, int(epoch)))


_LEAD_LATENT_3D = (
    '<p class="cryo-dash-lead">Each point is a particle in the space of three latent coordinates '
    "you choose (drag to rotate, scroll to zoom). Optional colour encodes another quantity from "
    "the analysis table. Use the histogram or discrete toggles in the legend beside the plot to "
    "filter particles; the continuous palette lives in the drop-down under the histogram.</p>"
)

_LATENT_3D_AXES_NOTE = (
    '<p class="cryo-dash-legend-note" style="margin:0 0 0.35rem;">Pick three <strong>different</strong> coordinates '
    "(z0, z1, …).</p>"
)

# ---------------------------------------------------------------------------
# Top-level pages
# ---------------------------------------------------------------------------


def index():
    """Landing page: feature flags and links depend on workdir / GPU / landscape outputs."""
    if not active_workdir(current_app):
        return render_template(
            "index.html",
            can_images=False,
            zdim=0,
            show_trajectory_creator=False,
            landscape_volpca_active=False,
            landscape_full_3d_active=False,
            exp_epoch=0,
            command_builder_only=True,
        )
    e: DashboardExperiment = g.dashboard_exp
    zdim = int(e.z.shape[1])
    landscape_full_3d_active = landscape_full_3d_ready(e.workdir, e.epoch)
    return render_template(
        "index.html",
        can_images=e.can_preview_particles,
        zdim=zdim,
        show_trajectory_creator=explorer_volumes_eligible(e),
        landscape_volpca_active=landscape_analysis_ready(e.workdir, e.epoch),
        landscape_full_3d_active=landscape_full_3d_active,
        exp_epoch=int(e.epoch),
        command_builder_only=False,
    )


def command_builder_page():
    """CryoDRGN CLI command builder (optional experiment defaults when a workdir is selected)."""
    e: DashboardExperiment | None = None
    if active_workdir(current_app):
        e = g.dashboard_exp
    return render_template(
        "command_builder.html",
        **command_builder_template_kwargs(e),
    )


def abinit_builder_redirect():
    """Legacy URL → command builder."""
    return _redirect("command_builder_page")


def filter_page_redirect():
    """Legacy ``/filter`` URL → scatter explorer."""
    return _redirect("explorer")


def api_save_selection():
    """Persist selected ``plot_df`` rows as a dataset-index pickle (and optional inverse)."""
    e: DashboardExperiment = g.dashboard_exp
    data = _request_json_dict()
    rows_raw = data.get("rows")
    if not isinstance(rows_raw, list) or len(rows_raw) == 0:
        return jsonify(error="No particles selected."), 400
    try:
        rows = sorted({int(r) for r in rows_raw})
    except (TypeError, ValueError):
        return jsonify(error="Invalid row list."), 400
    n_df = len(e.plot_df)
    if any(r < 0 or r >= n_df for r in rows):
        return jsonify(error="Row index out of range."), 400
    selected_ds = np.asarray(e.all_indices[rows], dtype=int)
    force = bool(data.get("force"))
    save_inverse = bool(data.get("save_inverse"))
    basename = (data.get("basename") or "indices").strip() or "indices"
    basename = os.path.basename(basename)
    if force:
        out_base = os.path.join(e.workdir, "indices")
    else:
        sel_dir = (data.get("sel_dir") or "").strip()
        dir_abs = os.path.abspath(sel_dir) if sel_dir else os.path.abspath(e.workdir)
        if not os.path.isdir(dir_abs):
            return jsonify(error=f"Directory does not exist: {dir_abs}"), 400
        out_base = os.path.join(dir_abs, basename)
    path_main = out_base + ".pkl"
    path_inv = out_base + "_inverse.pkl"
    parent = os.path.dirname(path_main)
    if parent:
        os.makedirs(parent, exist_ok=True)
    inverse_ds = np.setdiff1d(np.asarray(e.all_indices, dtype=int), selected_ds)
    try:
        with open(path_main, "wb") as fh:
            pickle.dump(selected_ds, fh)
        if save_inverse:
            with open(path_inv, "wb") as fh:
                pickle.dump(inverse_ds, fh)
    except OSError as err:
        logger.exception("save selection failed")
        return jsonify(error=str(err)), 500
    payload = {"ok": True, "path": path_main, "n_selected": int(selected_ds.size)}
    if save_inverse:
        payload["inverse_path"] = path_inv
    return jsonify(payload)


def api_covariate_threshold_rows():
    """All ``plot_df`` row indices where a numeric covariate passes a threshold.

    Matches the color-by histogram semantics (≥ or ≤ level), or an inclusive
    ``range_min`` … ``range_max`` interval (optionally inverted). Used so histogram
    selections cover the full dataset like saved lasso indices, not only points
    on the scatter subsample.
    """
    e: DashboardExperiment = g.dashboard_exp
    data = _request_json_dict()
    col = data.get("column")

    if not col or not isinstance(col, str) or col not in e.plot_df.columns:
        return jsonify(error="Invalid column.", rows=[]), 400
    if col == "index":
        return jsonify(error="Cannot threshold index column.", rows=[]), 400

    series = pd.Series(pd.to_numeric(e.plot_df[col], errors="coerce"))
    if not series.notna().any():
        return (
            jsonify(error="Covariate has no numeric values for thresholding.", rows=[]),
            400,
        )

    color_filter: dict[str, object]
    raw_r0 = data.get("range_min")
    raw_r1 = data.get("range_max")
    if raw_r0 is not None and raw_r1 is not None:
        try:
            r0 = float(raw_r0)
            r1 = float(raw_r1)
        except (TypeError, ValueError):
            return jsonify(error="Invalid range bounds.", rows=[]), 400
        color_filter = {
            "kind": "range",
            "range_min": min(r0, r1),
            "range_max": max(r0, r1),
            "invert_range": bool(data.get("invert_range")),
        }
        mask = plot_df_color_filter_mask(e.plot_df, col, color_filter)
        rows_arr = np.flatnonzero(mask).astype(np.int64, copy=False)
        return jsonify(rows=rows_arr.tolist(), n=int(rows_arr.size))

    raw_level = data.get("level")
    try:
        if raw_level is None:
            raise ValueError
        level = float(raw_level)
    except (TypeError, ValueError):
        return jsonify(error="Invalid threshold level.", rows=[]), 400

    color_filter = {
        "kind": "threshold",
        "level": level,
        "use_max": bool(data.get("use_max")),
    }
    mask = plot_df_color_filter_mask(e.plot_df, col, color_filter)
    rows_arr = np.flatnonzero(mask).astype(np.int64, copy=False)

    return jsonify(rows=rows_arr.tolist(), n=int(rows_arr.size))


def api_covariate_legend_context():
    """Histogram sample or discrete category list for covariate colour UI (pair-grid, 3D)."""
    e: DashboardExperiment = g.dashboard_exp
    data = _request_json_dict()
    col = data.get("column")
    scope = data.get("scope")

    if not col or not isinstance(col, str):
        return jsonify(error="Invalid column."), 400

    plot_df_override: pd.DataFrame | None = None
    if scope == "landscape_full_sampled":
        if not landscape_full_ready(e.workdir, e.epoch):
            return (
                jsonify(
                    error="Landscape full outputs not found for this epoch "
                    "(run cryodrgn analyze_landscape_full)."
                ),
                400,
            )
        try:
            plot_df_override = landscape_full_sampled_plot_df(e)
        except (FileNotFoundError, ValueError, OSError) as err:
            return jsonify(error=str(err)), 400
        if col not in plot_df_override.columns:
            return jsonify(error="Invalid column."), 400
    elif col not in e.plot_df.columns:
        return jsonify(error="Invalid column."), 400

    try:
        payload = covariate_legend_context_payload(e, col, plot_df=plot_df_override)
    except ValueError as err:
        return jsonify(error=str(err)), 400

    return jsonify(payload)


# ---------------------------------------------------------------------------
# Explorer + scatter routes
# ---------------------------------------------------------------------------


def explorer():
    """Interactive 2-D scatter with optional particle montage and volume explorer."""
    e: DashboardExperiment = g.dashboard_exp
    if not e.can_preview_particles:
        return (
            render_template(
                "no_images.html",
                reason=(
                    "Tilt-series data does not include single-particle thumbnails "
                    "in this view."
                ),
            ),
            200,
        )
    cols = e.numeric_columns
    dx, dy = _default_xy_cols(cols)
    initial_rows = load_plot_df_rows_from_plot_inds_file(
        e, current_app.config.get("FILTER_PLOT_INDS")
    )
    pc = int(current_app.config["PRELOAD_CPUS"])
    scatter_cap = _particle_explorer_scatter_max_points()
    scatter_plotted_n = min(int(len(e.plot_df)), scatter_cap)
    preload_limit = explorer_initial_preload_image_limit(scatter_plotted_n)
    preload_cache_step = explorer_cache_size_power10_step(scatter_plotted_n)
    return render_template(
        "particle_explorer.html",
        numeric_cols=cols,
        covariate_display_map=_covariate_display_map(cols),
        default_x=dx,
        default_y=dy,
        initial_rows=initial_rows,
        total_particles=int(len(e.all_indices)),
        workdir=e.workdir,
        preload_image_limit=preload_limit,
        preload_cache_step=preload_cache_step,
        preload_cpus=pc,
        preload_cache_time_hint=format_preload_cache_time_hint(pc),
        show_volume_explorer=explorer_volumes_eligible(e),
        explorer_scatter_max_points=_particle_explorer_scatter_max_points(),
        explorer_scatter_cap_from_env=_particle_explorer_scatter_cap_from_env(),
        scatter_plotted_n=scatter_plotted_n,
    )


def api_explorer_volume_media():
    """Decode montage cells to PNG or per-cell rotating GIF (ChimeraX); returns cache id for GIFs."""
    e: DashboardExperiment = g.dashboard_exp
    if not explorer_volumes_eligible(e):
        return jsonify(error=_EXPLORER_VOLUMES_INELIGIBLE_MSG), 400
    data = _request_json_dict()
    rows_raw = data.get("rows")
    if not isinstance(rows_raw, list) or len(rows_raw) == 0:
        return jsonify(error="No montage rows supplied."), 400
    try:
        rows = [int(r) for r in rows_raw]
    except (TypeError, ValueError):
        return jsonify(error="Invalid row indices."), 400
    if len(rows) > 64:
        return jsonify(error="At most 64 montage cells are supported."), 400
    mode = str(data.get("mode") or "static").strip().lower()
    try:
        if mode == "static":
            cc = int(data.get("chimerax_cpus", DEFAULT_CHIMERAX_PARALLEL))
            blobs, cache_token = generate_montage_volume_pngs(e, rows, chimerax_cpus=cc)
            b64s = [base64.standard_b64encode(b).decode("ascii") for b in blobs]
            return jsonify(
                ok=True,
                format="png",
                images=b64s,
                rows=rows,
                volume_cache_id=cache_token,
            )
        if mode in ("animate", "gif"):
            cache_id = data.get("volume_cache_id")
            if not cache_id or not isinstance(cache_id, str):
                return (
                    jsonify(
                        error="Animate requires volume_cache_id from Generate volumes.",
                    ),
                    400,
                )
            raw_ci = data.get("cell_index")
            if raw_ci is None:
                return jsonify(error="cell_index is required for animate."), 400
            try:
                cell_index = int(raw_ci)
            except (TypeError, ValueError):
                return jsonify(error="cell_index must be an integer."), 400
            gf = int(data.get("gif_frames", DEFAULT_GIF_FRAMES))
            cc = int(data.get("chimerax_cpus", DEFAULT_CHIMERAX_PARALLEL))
            gif_bytes = volume_cell_gif_from_cache(
                cache_id,
                cell_index,
                rows_expected=tuple(rows),
                gif_frames=gf,
                chimerax_cpus=cc,
            )
            b64_one = base64.standard_b64encode(gif_bytes).decode("ascii")
            return jsonify(
                ok=True,
                format="gif",
                image=b64_one,
                cell_index=cell_index,
                rows=rows,
            )
        return jsonify(error='mode must be "static" or "animate".'), 400
    except EnvironmentError as err:
        return jsonify(error=str(err), need_chimerax=True), 503
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("explorer volume media failed")
        return jsonify(error=str(err)), 500


def api_scatter():
    """Plotly JSON for the explorer / filter scatter (subsample caps via query flags)."""
    e: DashboardExperiment = g.dashboard_exp
    xcol = request.args.get("x", e.numeric_columns[0])
    ycol = request.args.get("y", e.numeric_columns[0])
    ccol = request.args.get("color") or "none"
    filter_ui = request.args.get("filter_ui") == "1"
    explorer_scatter = request.args.get("explorer_scatter") == "1"
    full = request.args.get("full") == "1"
    if xcol not in e.plot_df.columns or ycol not in e.plot_df.columns:
        return jsonify(error="bad axis column"), 400
    if ccol != "none" and ccol not in e.plot_df.columns:
        return jsonify(error="bad color column"), 400
    if filter_ui:
        max_pts = _filter_ui_scatter_max_points()
    elif full:
        max_pts = None
    elif explorer_scatter:
        max_pts = _particle_explorer_scatter_max_points()
    else:
        max_pts = 200_000
        raw_max_pts = request.args.get("max_points")
        if raw_max_pts:
            try:
                max_pts = max(1, min(int(raw_max_pts), 200_000))
            except ValueError:
                pass
    preselect_rows, pre_err = _parse_preselect_rows_param(
        request.args.get("preselect_rows"),
    )
    if pre_err:
        return jsonify(error=pre_err), 400
    marker_size = 4.0
    raw_ms = request.args.get("marker_size")
    if raw_ms:
        try:
            marker_size = max(0.5, min(float(raw_ms), 20))
        except ValueError:
            pass
    try:
        js = scatter_json(
            e,
            xcol,
            ycol,
            None if ccol == "none" else ccol,
            max_points=max_pts,
            preselect_plot_df_rows=preselect_rows,
            use_webgl=False,
            marker_size=marker_size,
            continuous_palette=request.args.get("palette"),
        )
    except Exception as err:
        logger.exception("scatter plot failed")
        return jsonify(error=str(err)), 500
    return Response(js, mimetype="application/json")


def latent_3d_page():
    """3-D latent visualizer shell (requires ``zdim >= 3``)."""
    e: DashboardExperiment = g.dashboard_exp
    zdim = int(e.z.shape[1])
    if zdim < 3:
        return (
            render_template("pair_grid_need_more_cols.html", kind="z3", n=zdim),
            200,
        )
    axis_cols = [f"z{i}" for i in range(zdim)]
    cols = e.numeric_columns
    return render_template(
        "latent_3d.html",
        page_title="3D latent space visualizer · cryoDRGN",
        nav_bar_title="3D visualizer",
        lead_html=_LEAD_LATENT_3D,
        axis_cols=axis_cols,
        axes_fieldset_legend="Latent axes",
        axes_fieldset_note_html=_LATENT_3D_AXES_NOTE,
        numeric_cols=cols,
        covariate_display_map=_covariate_display_map(cols),
        default_x="z0",
        default_y="z1",
        default_z="z2",
        scatter3d_url=url_for("api_scatter3d_z"),
        legend_context_body_extra=None,
        show_landscape_vol_animations=False,
        show_vol_landscape_quick_actions=False,
    )


def landscape_full_3d_page():
    """3-D scatter of sampled particles in ``vol_pca_all.pkl`` space (optional ChimeraX GIF previews)."""
    e: DashboardExperiment = g.dashboard_exp
    if not landscape_full_3d_ready(e.workdir, e.epoch):
        return render_template(
            "volume_landscape_3d_need_outputs.html",
            exp_epoch=int(e.epoch),
            error_message="",
        )
    try:
        sampled = landscape_full_sampled_plot_df(e)
    except (FileNotFoundError, ValueError, OSError) as err:
        logger.exception("landscape full 3D: failed to build sampled table")
        return render_template(
            "volume_landscape_3d_need_outputs.html",
            exp_epoch=int(e.epoch),
            error_message=str(err),
        )
    vol_axes = landscape_full_vol_pca_axis_columns(sampled)
    if len(vol_axes) < 3:
        return render_template(
            "volume_landscape_3d_need_outputs.html",
            exp_epoch=int(e.epoch),
            error_message=(
                "Expected at least three columns in vol_pca_all.pkl (landscape_vol_PC1 …) "
                "for this epoch."
            ),
        )
    dx, dy, dz = vol_axes[0], vol_axes[1], vol_axes[2]
    cols = landscape_full_sampled_numeric_covariates(sampled)
    vol_evr = _landscape_full_vol_pc_explained_variance(e.workdir, e.epoch)
    vol_anim = landscape_analysis_ready(e.workdir, e.epoch)
    return render_template(
        "latent_3d.html",
        page_title="3D volume landscapes · cryoDRGN",
        nav_bar_title="3D volume landscapes",
        lead_html=LANDSCAPE_FULL_3D_LEAD_HTML,
        axis_cols=vol_axes,
        axes_fieldset_legend="Volume PCA axes",
        axes_fieldset_note_html="",
        numeric_cols=cols,
        covariate_display_map=_covariate_display_map(
            cols,
            vol_pc_explained_variance_ratio=vol_evr,
        ),
        default_x=dx,
        default_y=dy,
        default_z=dz,
        scatter3d_url=url_for("api_scatter3d_z_landscape_full"),
        legend_context_body_extra=LANDSCAPE_FULL_3D_LEGEND_CONTEXT_EXTRA,
        show_landscape_vol_animations=vol_anim,
        show_vol_landscape_quick_actions=True,
    )


def api_scatter3d_z():
    """Plotly JSON for the 3-D latent scatter."""
    e: DashboardExperiment = g.dashboard_exp
    color_filter = None
    if request.method == "POST":
        payload = _request_json_dict()
        xcol = str(payload.get("x") or "z0")
        ycol = str(payload.get("y") or "z1")
        zcol = str(payload.get("z") or "z2")
        ccol = str(payload.get("color") or "none")
        raw_palette = payload.get("palette")
        if ccol != "none" and ccol not in e.plot_df.columns:
            return jsonify(error="bad color column"), 400
        try:
            if ccol != "none":
                color_filter = _parse_color_filter_for_column(ccol, payload)
            discrete_label_colors = _parse_optional_discrete_label_colors(
                payload.get("discrete_label_colors")
            )
            no_sub = _scatter3d_no_subsample_for_discrete_gif_frame(
                payload, color_filter
            )
            js = scatter3d_z_json(
                e,
                xcol,
                ycol,
                zcol,
                None if ccol == "none" else ccol,
                continuous_palette=raw_palette,
                color_filter=color_filter,
                discrete_label_colors=discrete_label_colors,
                no_subsample=no_sub,
            )
        except ValueError as err:
            return jsonify(error=str(err)), 400
        except Exception as err:
            logger.exception("3D latent space visualizer failed")
            return jsonify(error=str(err)), 500
        return Response(js, mimetype="application/json")

    xcol = request.args.get("x", "z0")
    ycol = request.args.get("y", "z1")
    zcol = request.args.get("z", "z2")
    ccol = request.args.get("color") or "none"
    if ccol != "none" and ccol not in e.plot_df.columns:
        return jsonify(error="bad color column"), 400
    try:
        js = scatter3d_z_json(
            e,
            xcol,
            ycol,
            zcol,
            None if ccol == "none" else ccol,
            continuous_palette=request.args.get("palette"),
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("3D latent space visualizer failed")
        return jsonify(error=str(err)), 500
    return Response(js, mimetype="application/json")


def api_scatter3d_z_landscape_full():
    """Plotly JSON for 3D scatter: axes are ``vol_pca_all.pkl`` columns (``landscape_vol_PC*``)."""
    e: DashboardExperiment = g.dashboard_exp
    if not landscape_full_3d_ready(e.workdir, e.epoch):
        return (
            jsonify(
                error=(
                    f"Need landscape.{int(e.epoch)}/landscape_full/ with z.sampled.pkl, ind.sampled.pkl, "
                    "and vol_pca_all.pkl (at least three PCA columns)."
                )
            ),
            400,
        )
    try:
        sampled = landscape_full_sampled_plot_df(e)
    except (FileNotFoundError, ValueError, OSError) as err:
        return jsonify(error=str(err)), 400

    vol_axes = landscape_full_vol_pca_axis_columns(sampled)
    if len(vol_axes) < 3:
        return (
            jsonify(
                error="vol_pca_all.pkl must yield at least three landscape_vol_PC* columns in the sampled table."
            ),
            400,
        )
    ax_allow = frozenset(vol_axes)
    d0, d1, d2 = vol_axes[0], vol_axes[1], vol_axes[2]
    evr = _landscape_full_vol_pc_explained_variance(e.workdir, e.epoch)

    def _vol_landscape_scene_titles(xc: str, yc: str, zc: str) -> tuple[str, str, str]:
        return (
            landscape_vol_pc_column_pretty_label(xc, evr),
            landscape_vol_pc_column_pretty_label(yc, evr),
            landscape_vol_pc_column_pretty_label(zc, evr),
        )

    color_filter = None
    if request.method == "POST":
        payload = _request_json_dict()
        xcol = str(payload.get("x") or d0)
        ycol = str(payload.get("y") or d1)
        zcol = str(payload.get("z") or d2)
        ccol = str(payload.get("color") or "none")
        raw_palette = payload.get("palette")
        if ccol != "none" and ccol not in sampled.columns:
            return jsonify(error="bad color column"), 400
        try:
            if ccol != "none":
                color_filter = _parse_color_filter_for_column(ccol, payload)
            discrete_label_colors = _parse_optional_discrete_label_colors(
                payload.get("discrete_label_colors")
            )
            no_sub = _scatter3d_no_subsample_for_discrete_gif_frame(
                payload, color_filter
            )
            js = scatter3d_z_json(
                e,
                xcol,
                ycol,
                zcol,
                None if ccol == "none" else ccol,
                continuous_palette=raw_palette,
                color_filter=color_filter,
                discrete_label_colors=discrete_label_colors,
                plot_df=sampled,
                uirevision="scatter3d_z_landscape_full",
                xyz_axes_allowed=ax_allow,
                scene_axis_titles=_vol_landscape_scene_titles(xcol, ycol, zcol),
                vol_pc_explained_variance_ratio=evr,
                volume_landscape_3d_style=True,
                no_subsample=no_sub,
            )
        except ValueError as err:
            return jsonify(error=str(err)), 400
        except Exception as err:
            logger.exception("3D volume landscapes scatter failed")
            return jsonify(error=str(err)), 500
        return Response(js, mimetype="application/json")

    xcol = request.args.get("x", d0)
    ycol = request.args.get("y", d1)
    zcol = request.args.get("z", d2)
    ccol = request.args.get("color") or "none"
    if ccol != "none" and ccol not in sampled.columns:
        return jsonify(error="bad color column"), 400
    try:
        js = scatter3d_z_json(
            e,
            xcol,
            ycol,
            zcol,
            None if ccol == "none" else ccol,
            continuous_palette=request.args.get("palette"),
            plot_df=sampled,
            uirevision="scatter3d_z_landscape_full",
            xyz_axes_allowed=ax_allow,
            scene_axis_titles=_vol_landscape_scene_titles(xcol, ycol, zcol),
            vol_pc_explained_variance_ratio=evr,
            volume_landscape_3d_style=True,
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("3D volume landscapes scatter failed")
        return jsonify(error=str(err)), 500
    return Response(js, mimetype="application/json")


def api_latent3d_preview_png():
    """PNG snapshot of the 3D latent scatter (matplotlib)."""
    e: DashboardExperiment = g.dashboard_exp
    xcol = request.args.get("x", "z0")
    ycol = request.args.get("y", "z1")
    zcol = request.args.get("z", "z2")
    ccol = request.args.get("color") or "none"
    try:
        elev = float(request.args.get("elev", "22"))
        azim = float(request.args.get("azim", "-65"))
    except ValueError:
        return jsonify(error="elev/azim must be numeric"), 400
    if ccol != "none" and ccol not in e.plot_df.columns:
        return jsonify(error="bad color column"), 400
    try:
        png = scatter3d_z_preview_png(
            e,
            xcol,
            ycol,
            zcol,
            None if ccol == "none" else ccol,
            continuous_palette=request.args.get("palette"),
            elev=elev,
            azim=azim,
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("3D latent preview PNG failed")
        return jsonify(error=str(err)), 500
    return Response(png, mimetype="image/png")


def api_latent3d_plot_gif_from_png_frames():
    """Assemble base64-encoded PNG frames into an animated GIF (browser-captured Plotly views)."""
    _ = g.dashboard_exp
    data = _request_json_dict()
    frames = data.get("frames")
    if not isinstance(frames, list):
        return jsonify(error="frames must be a list of base64 PNG strings."), 400
    str_frames = [str(x) for x in frames if isinstance(x, (str, bytes))]
    if len(str_frames) != len(frames):
        return jsonify(error="each frame must be a string."), 400
    dur_raw = data.get("durations_ms", data.get("duration_ms", 100))
    try:
        gif_bytes = png_base64_frames_to_gif_bytes(str_frames, durations_ms=dur_raw)
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except RuntimeError as err:
        return jsonify(error=str(err)), 500
    except Exception as err:
        logger.exception("latent3d GIF assembly failed")
        return jsonify(error=str(err)), 500
    return jsonify(
        gif_b64=base64.standard_b64encode(gif_bytes).decode("ascii"),
    )


def api_latent3d_landscape_full_discrete_gif():
    """Discrete-level GIF via Matplotlib frames (parallel); avoids browser WebGL capture."""
    e: DashboardExperiment = g.dashboard_exp
    if not landscape_full_3d_ready(e.workdir, e.epoch):
        return (
            jsonify(
                error=(
                    f"Need landscape.{int(e.epoch)}/landscape_full/ with z.sampled.pkl, ind.sampled.pkl, "
                    "and vol_pca_all.pkl (at least three PCA columns)."
                )
            ),
            400,
        )
    try:
        sampled = landscape_full_sampled_plot_df(e)
    except (FileNotFoundError, ValueError, OSError) as err:
        return jsonify(error=str(err)), 400

    vol_axes = landscape_full_vol_pca_axis_columns(sampled)
    if len(vol_axes) < 3:
        return (
            jsonify(
                error="vol_pca_all.pkl must yield at least three landscape_vol_PC* columns in the sampled table."
            ),
            400,
        )
    ax_allow = frozenset(vol_axes)
    d0, d1, d2 = vol_axes[0], vol_axes[1], vol_axes[2]
    evr = _landscape_full_vol_pc_explained_variance(e.workdir, e.epoch)

    def _scene_titles(xc: str, yc: str, zc: str) -> tuple[str, str, str]:
        return (
            landscape_vol_pc_column_pretty_label(xc, evr),
            landscape_vol_pc_column_pretty_label(yc, evr),
            landscape_vol_pc_column_pretty_label(zc, evr),
        )

    data = _request_json_dict()
    keys_raw = data.get("discrete_keys")
    if not isinstance(keys_raw, list) or len(keys_raw) == 0:
        return jsonify(error="discrete_keys must be a non-empty list."), 400
    keys = [str(x) for x in keys_raw]
    if len(keys) > 120:
        return jsonify(error="Too many discrete levels (max 120)."), 400

    xcol = str(data.get("x") or d0)
    ycol = str(data.get("y") or d1)
    zcol = str(data.get("z") or d2)
    ccol = str(data.get("color") or "none")
    if ccol == "none" or ccol not in sampled.columns:
        return jsonify(error="A discrete colour covariate is required."), 400
    if not _lower_color_series_is_discrete(cast(pd.Series, sampled[ccol])):
        return jsonify(error="Colour column must be discrete for this export."), 400

    discrete_label_colors = _parse_optional_discrete_label_colors(
        data.get("discrete_label_colors")
    )

    dur_raw = data.get("frame_duration_ms", data.get("durations_ms", 960))

    titles = _scene_titles(xcol, ycol, zcol)

    def _png_b64_for_key(k: str) -> str:
        png = scatter3d_landscape_full_discrete_level_png_bytes(
            sampled,
            xcol,
            ycol,
            zcol,
            ccol,
            k,
            discrete_label_colors=discrete_label_colors,
            xyz_axes_allowed=ax_allow,
            scene_axis_titles=titles,
        )
        return base64.standard_b64encode(png).decode("ascii")

    n_workers = min(max(1, os.cpu_count() or 4), len(keys), 16)
    try:
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            png_b64_list = list(pool.map(_png_b64_for_key, keys))
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("discrete-level landscape GIF frame render failed")
        return jsonify(error=str(err)), 500

    if len(png_b64_list) < 2:
        png_b64_list = png_b64_list + [png_b64_list[-1]]

    try:
        gif_bytes = png_base64_frames_to_gif_bytes(png_b64_list, durations_ms=dur_raw)
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except RuntimeError as err:
        return jsonify(error=str(err)), 500
    except Exception as err:
        logger.exception("discrete-level landscape GIF assembly failed")
        return jsonify(error=str(err)), 500

    return jsonify(
        gif_b64=base64.standard_b64encode(gif_bytes).decode("ascii"),
    )


def api_preview_montage():
    e: DashboardExperiment = g.dashboard_exp
    raw = request.args.get("rows", "")
    if not raw.strip():
        return Response(montage_bytes(e, []), mimetype="image/png")
    parts = [p.strip() for p in raw.split(",") if p.strip()]
    try:
        idxs = [int(p) for p in parts]
    except ValueError:
        return jsonify(error="rows must be integers"), 400
    return Response(montage_bytes(e, idxs), mimetype="image/png")


def api_preload_images():
    """Return base64 JPEG thumbnails for a preload cache.

    Selection (see :func:`sample_plot_df_rows_for_preload`): half of the picks are
    uniform random; half favor large k-th nearest-neighbor distance in the
    scatter plane, with a coarse XY grid so outliers do not pile into one bin.

    Use **POST** with a JSON body when ``selected_rows`` is large (lasso
    selections); query strings hit proxy / server URI length limits (414).

    With ``restrict_to_scatter_plot`` (particle explorer), thumbnails are drawn only
    from rows visible in the scatter subsample (same rule as ``api_scatter``, default
    200k points, seed 0).

    Delta responses (``response_mode: "delta"``) include ``batch_elapsed``: wall seconds
    for this request only (sampling + encode). ``elapsed`` remains cumulative for the
    server-side cache entry.

    POST with ``{"invalidate_cache": true}`` drops server-side thumbnails for this
    experiment epoch (used when the browser clears its local preload state).
    """
    e: DashboardExperiment = g.dashboard_exp
    if request.method == "POST":
        _inv_payload = _request_json_dict()
        if _inv_payload.get("invalidate_cache"):
            n_cleared = clear_preload_cache_for_experiment(e)
            return jsonify(ok=True, cleared=n_cleared)

    if not e.can_preview_particles:
        return jsonify(rows=[], images=[], elapsed=0)

    cols = e.plot_df.columns
    restrict_list: list[int] | None = None
    initial_list: list[int] | None = None
    restrict_to_scatter = False
    scatter_max_points = 200_000

    if request.method == "POST":
        data = _request_json_dict()
        xcol = str(data.get("x") or "")
        ycol = str(data.get("y") or "")
        response_mode = str(data.get("response_mode") or "full")
        restrict_to_scatter = bool(data.get("restrict_to_scatter_plot"))
        smp = data.get("scatter_max_points")
        if smp is not None and smp != "":
            try:
                scatter_max_points = max(1, min(int(smp), 2_000_000))
            except (TypeError, ValueError):
                scatter_max_points = 200_000
        try:
            max_images = _parse_preload_image_limit(data.get("cache_size"))
        except ValueError as err:
            return jsonify(error=str(err)), 400
        raw_list = data.get("selected_rows")
        if raw_list is None:
            restrict_list = None
        elif not isinstance(raw_list, list):
            return jsonify(error="selected_rows must be a JSON array of integers."), 400
        else:
            try:
                restrict_list = [int(r) for r in raw_list]
            except (TypeError, ValueError):
                return (
                    jsonify(error="selected_rows must be a JSON array of integers."),
                    400,
                )
            if not restrict_list:
                restrict_list = None
        raw_initial_list = data.get("initial_rows")
        if raw_initial_list is None:
            initial_list = None
        elif not isinstance(raw_initial_list, list):
            return jsonify(error="initial_rows must be a JSON array of integers."), 400
        else:
            try:
                initial_list = [int(r) for r in raw_initial_list]
            except (TypeError, ValueError):
                return (
                    jsonify(error="initial_rows must be a JSON array of integers."),
                    400,
                )
            if not initial_list:
                initial_list = None
    else:
        xcol = request.args.get("x", "")
        ycol = request.args.get("y", "")
        response_mode = str(request.args.get("response_mode") or "full")
        restrict_to_scatter = request.args.get("restrict_to_scatter_plot") == "1"
        sms = (request.args.get("scatter_max_points") or "").strip()
        if sms:
            try:
                scatter_max_points = max(1, min(int(sms), 2_000_000))
            except ValueError:
                scatter_max_points = 200_000
        try:
            max_images = _parse_preload_image_limit(request.args.get("cache_size"))
        except ValueError as err:
            return jsonify(error=str(err)), 400
        raw_sel = (request.args.get("selected_rows") or "").strip()
        if raw_sel:
            try:
                restrict_list = [
                    int(x.strip()) for x in raw_sel.split(",") if x.strip()
                ]
            except ValueError:
                return (
                    jsonify(error="selected_rows must be comma-separated integers."),
                    400,
                )
            if not restrict_list:
                restrict_list = None

    if xcol not in cols:
        xcol = str(cols[0])
    if ycol not in cols:
        ycol = str(cols[min(1, len(cols) - 1)])
    delta_response = response_mode == "delta"

    scatter_pool: list[int] | None = None
    if restrict_to_scatter:
        scatter_pool = [
            int(x)
            for x in plot_df_row_indices_for_explorer_scatter(
                e.plot_df, scatter_max_points
            )
        ]
    sp_set: set[int] | None = set(scatter_pool) if scatter_pool is not None else None

    if restrict_list is not None:
        if sp_set is not None:
            preload_restrict_list = [r for r in restrict_list if r in sp_set]
        else:
            preload_restrict_list = list(restrict_list)
    elif scatter_pool is not None:
        preload_restrict_list = scatter_pool
    else:
        preload_restrict_list = initial_list

    def _preload_restriction_key() -> tuple:
        if restrict_list is not None and sp_set is not None:
            return (
                "scatter_sel",
                scatter_max_points,
                tuple(sorted({int(r) for r in restrict_list})),
            )
        if restrict_list is not None:
            return ("sel", tuple(sorted({int(r) for r in restrict_list})))
        if restrict_to_scatter:
            return ("scatter", scatter_max_points)
        if initial_list is not None:
            return ("initial", tuple(sorted({int(r) for r in initial_list})))
        return ("full",)

    key = (e.epoch, e.kmeans_folder_id, xcol, ycol, _preload_restriction_key())
    cpus = int(current_app.config.get("PRELOAD_CPUS") or 4)

    def _encode_indices(global_indices: list[int]) -> list[str]:
        parallel_threshold = max(128, cpus * 32)
        if cpus > 1 and len(global_indices) >= parallel_threshold:
            from concurrent.futures import ProcessPoolExecutor

            chunk_sz = -(-len(global_indices) // cpus)
            chunks = [
                global_indices[i : i + chunk_sz]
                for i in range(0, len(global_indices), chunk_sz)
            ]
            with ProcessPoolExecutor(max_workers=len(chunks)) as pool:
                futures = [
                    pool.submit(
                        encode_particle_batch, e.particles_path, e.datadir, ch, 96
                    )
                    for ch in chunks
                ]
                imgs: list[str] = []
                for f in futures:
                    imgs.extend(f.result())
                return imgs
        return encode_particle_batch(e.particles_path, e.datadir, global_indices, 96)

    cached = PRELOAD_CACHE.get(key)
    if cached:
        cached_rows, cached_imgs, cached_elapsed = cached
        if len(cached_rows) >= max_images:
            if delta_response:
                return jsonify(
                    rows=[],
                    images=[],
                    elapsed=cached_elapsed,
                    total_cached=len(cached_rows),
                    batch_elapsed=0.0,
                )
            return jsonify(
                rows=cached_rows[:max_images],
                images=cached_imgs[:max_images],
                elapsed=cached_elapsed,
                total_cached=min(len(cached_rows), max_images),
            )
        t0 = time.monotonic()
        add_rows, add_global_indices = sample_plot_df_rows_for_preload(
            e,
            xcol,
            ycol,
            restrict_to_rows=preload_restrict_list,
            exclude_rows=cached_rows,
            max_images=max_images - len(cached_rows),
        )
        if not add_global_indices:
            if delta_response:
                return jsonify(
                    rows=[],
                    images=[],
                    elapsed=cached_elapsed,
                    total_cached=len(cached_rows),
                    batch_elapsed=0.0,
                )
            return jsonify(rows=cached_rows, images=cached_imgs, elapsed=cached_elapsed)
        add_imgs = _encode_indices(add_global_indices)
        rows = cached_rows + add_rows
        imgs = cached_imgs + add_imgs
        batch_elapsed = round(time.monotonic() - t0, 1)
        elapsed = round(cached_elapsed + batch_elapsed, 1)
        PRELOAD_CACHE[key] = (rows, imgs, elapsed)
        if delta_response:
            return jsonify(
                rows=add_rows,
                images=add_imgs,
                elapsed=elapsed,
                total_cached=len(rows),
                batch_elapsed=batch_elapsed,
            )
        return jsonify(rows=rows, images=imgs, elapsed=elapsed)

    t0 = time.monotonic()
    rows, global_indices = sample_plot_df_rows_for_preload(
        e, xcol, ycol, restrict_to_rows=preload_restrict_list, max_images=max_images
    )
    if not global_indices:
        return jsonify(rows=[], images=[], elapsed=0.0)

    imgs = _encode_indices(global_indices)

    elapsed = round(time.monotonic() - t0, 1)
    PRELOAD_CACHE[key] = (rows, imgs, elapsed)
    return jsonify(rows=rows, images=imgs, elapsed=elapsed, total_cached=len(rows))
