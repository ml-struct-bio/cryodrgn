"""Flask app for the cryoDRGN analysis dashboard.

Responsibilities are split across small sibling modules:

* :mod:`cryodrgn.dashboard.context` — workdir/epoch resolution + template
  context injectors + long-lived caches.
* :mod:`cryodrgn.dashboard.route_helpers` — small request/parsing helpers for
  route handlers.
* :mod:`cryodrgn.dashboard.trajectory` — pure-logic helpers for the trajectory
  creator (axis validation, anchor parsing, graph traversal, JSON payloads).
* :mod:`cryodrgn.dashboard.preload` — thumbnail sampling / encoding for the
  explorer preview montage and hover pre-load cache.

Shared CSS and JavaScript for templates live under ``cryodrgn/dashboard/static/``
(and are linked with :func:`flask.url_for` ``'static'``).

This file wires HTTP routes onto those helpers and keeps the request-handling
layer thin.
"""

from __future__ import annotations

import base64
import logging
import os
import pickle
import time
import uuid

import numpy as np
import pandas as pd
from flask import (
    Flask,
    Response,
    current_app,
    g,
    jsonify,
    render_template,
    request,
)

from cryodrgn.dashboard.context import (
    PRELOAD_CACHE,
    abbrev_middle,
    active_workdir,
    api_set_epoch,
    clear_preload_cache_for_experiment,
    api_set_workdir,
    bind_dashboard_exp,
    command_builder_template_kwargs,
    discover_cryodrgn_workdirs,
    inject_meta,
    inject_meta_command_builder_only,
    _request_json_dict,
)
from cryodrgn.dashboard.data import DashboardExperiment, list_z_epochs
from cryodrgn.dashboard.particle_explorer import (
    DEFAULT_CHIMERAX_PARALLEL,
    DEFAULT_GIF_FRAMES,
    explorer_volumes_eligible,
    generate_montage_volume_pngs,
    generate_trajectory_volume_pngs,
    save_cached_volumes_to_dir,
    volume_cell_gif_from_cache,
)
from cryodrgn.dashboard.landscape_volpca import (
    animation_payload_b64,
    generate_landscape_volume_animations,
    landscape_analysis_ready,
    landscape_dir_for_epoch,
    landscape_volpca_scatter_json,
    list_landscape_epochs,
    load_sketch_state_labels,
    meta_for_api,
    save_landscape_animations,
)

from cryodrgn.dashboard.plots import (
    covariate_legend_context_payload,
    pair_grid_figure_aspect_ratio,
    pair_grid_margin_fractions_for_js,
    pair_grid_png,
    pair_grid_skeleton_placeholder_layout,
    plot_df_color_filter_mask,
    plot_df_row_indices_for_explorer_scatter,
    scatter3d_z_json,
    scatter3d_z_preview_png,
    scatter_json,
)
from cryodrgn.dashboard.preload import (
    encode_particle_batch,
    explorer_cache_size_power10_step,
    explorer_initial_preload_image_limit,
    format_preload_cache_time_hint,
    load_plot_df_rows_from_plot_inds_file,
    montage_bytes,
    particle_thumbnail_b64_from_row,
    sample_plot_df_rows_for_preload,
)
from cryodrgn.dashboard.route_helpers import (
    _EXPLORER_VOLUMES_INELIGIBLE_MSG,
    _add_direct_anchor_pidx,
    _covariate_display_map,
    _default_xy_cols,
    _filter_ui_scatter_max_points,
    _parse_color_filter_for_column,
    _parse_optional_discrete_label_colors,
    _parse_pairplot_request,
    _parse_preload_image_limit,
    _parse_preselect_rows_param,
    _particle_explorer_scatter_cap_from_env,
    _particle_explorer_scatter_max_points,
    _redirect,
    _trajectory_eligibility_error,
)
from cryodrgn.dashboard.trajectory import (
    compute_trajectory_latent_path,
    default_trajectory_endpoints_xy,
    has_pc_columns,
    has_umap_columns,
    load_kmeans_center_indices,
    parse_anchor_indices_txt,
    parse_trajectory_request_body,
    random_dataset_indices,
    trajectory_anchor_mode_params,
    trajectory_anchor_payload_from_indices,
    trajectory_axes_from_payload,
    trajectory_default_xy_cols,
    trajectory_plot_axis_columns,
    trajectory_shared_json_payload,
    validate_trajectory_plot_axes,
)

logger = logging.getLogger(__name__)

# Default rotation-frame count for the volume sketched landscape explorer GIFs only.
LANDSCAPE_SKETCH_GIF_FRAMES_DEFAULT = 20

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_TEMPLATE_DIR = os.path.join(_THIS_DIR, "templates")
_STATIC_DIR = os.path.join(_THIS_DIR, "static")


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
            exp_epoch=0,
            command_builder_only=True,
        )
    e: DashboardExperiment = g.dashboard_exp
    return render_template(
        "index.html",
        can_images=e.can_preview_particles,
        zdim=int(e.z.shape[1]),
        show_trajectory_creator=explorer_volumes_eligible(e),
        landscape_volpca_active=landscape_analysis_ready(e.workdir, e.epoch),
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

    if not col or not isinstance(col, str) or col not in e.plot_df.columns:
        return jsonify(error="Invalid column."), 400
    try:
        payload = covariate_legend_context_payload(e, col)
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
    cols = e.numeric_columns
    return render_template(
        "latent_3d.html",
        z_cols=[f"z{i}" for i in range(zdim)],
        numeric_cols=cols,
        covariate_display_map=_covariate_display_map(cols),
        default_x="z0",
        default_y="z1",
        default_z="z2",
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
            js = scatter3d_z_json(
                e,
                xcol,
                ycol,
                zcol,
                None if ccol == "none" else ccol,
                continuous_palette=raw_palette,
                color_filter=color_filter,
            )
        except ValueError as err:
            return jsonify(error=str(err)), 400
        except Exception as err:
            logger.exception("3-D Latent Space Visualizer failed")
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
        logger.exception("3-D Latent Space Visualizer failed")
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
        logger.exception("3-D latent preview PNG failed")
        return jsonify(error=str(err)), 500
    return Response(png, mimetype="image/png")


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


# ---------------------------------------------------------------------------
# Pair-grid routes
# ---------------------------------------------------------------------------


def pairplot_page():
    """Pair-grid UI: latent vs latent with chosen color covariate (needs PCA + numeric color column)."""
    e: DashboardExperiment = g.dashboard_exp
    zdim = int(e.z.shape[1])
    if zdim < 2:
        return (
            render_template("pair_grid_need_more_cols.html", kind="zdim", n=zdim),
            200,
        )
    if not has_pc_columns(e):
        return (
            render_template("pair_grid_need_more_cols.html", kind="pca", n=zdim),
            200,
        )
    z_names = frozenset(f"z{i}" for i in range(zdim))
    color_choices = [c for c in e.numeric_columns if c not in z_names]
    if not color_choices:
        return (
            render_template("pair_grid_need_more_cols.html", kind="numeric", n=0),
            200,
        )
    default_color = next(
        (c for c in ("labels", "znorm", "UMAP1", "PC1") if c in color_choices),
        color_choices[0],
    )
    return render_template(
        "pair_grid.html",
        color_choices=color_choices,
        covariate_display_map=_covariate_display_map(color_choices),
        default_color=default_color,
        has_umap=has_umap_columns(e),
        zdim=zdim,
        skeleton_placeholder_cells=pair_grid_skeleton_placeholder_layout(zdim),
        pairplot_margin_fractions=pair_grid_margin_fractions_for_js(),
        pairplot_fig_aspect=pair_grid_figure_aspect_ratio(zdim),
        pairplot_save_default_name="zdim_pairplot.png",
        pairplot_save_default_dir=os.path.join(e.workdir, f"analyze.{e.epoch}"),
    )


def api_pairplot():
    """Render the z-dimensional pair grid to an in-memory PNG; returns base64 + cell layout JSON."""
    e: DashboardExperiment = g.dashboard_exp
    payload = _request_json_dict()
    try:
        color_col, diagonal_emb, upper_style, palette = _parse_pairplot_request(
            e, payload
        )
        color_filter = _parse_color_filter_for_column(color_col, payload)
        discrete_label_colors = _parse_optional_discrete_label_colors(
            payload.get("discrete_label_colors")
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400
    try:
        png, cells = pair_grid_png(
            e,
            lower_color_col=color_col,
            diagonal_emb=diagonal_emb,
            upper_style=upper_style,
            continuous_palette=palette,
            color_filter=color_filter,
            discrete_label_colors=discrete_label_colors,
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("pair grid failed")
        return jsonify(error=str(err)), 500

    return jsonify(
        png_b64=base64.standard_b64encode(png).decode("ascii"),
        cells=cells,
    )


def api_save_pairplot_png():
    """Write the same pair grid as ``api_pairplot`` to ``analyze.{epoch}/`` on disk."""
    e: DashboardExperiment = g.dashboard_exp
    payload = _request_json_dict()
    try:
        color_col, diagonal_emb, upper_style, palette = _parse_pairplot_request(
            e, payload
        )
        color_filter = _parse_color_filter_for_column(color_col, payload)
        discrete_label_colors = _parse_optional_discrete_label_colors(
            payload.get("discrete_label_colors")
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400

    raw_name = str(payload.get("filename") or "zdim_pairplot.png").strip()
    filename = os.path.basename(raw_name) or "zdim_pairplot.png"
    if not filename.lower().endswith(".png"):
        filename += ".png"
    out_dir = os.path.join(e.workdir, f"analyze.{e.epoch}")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, filename)
    try:
        png, _ = pair_grid_png(
            e,
            lower_color_col=color_col,
            diagonal_emb=diagonal_emb,
            upper_style=upper_style,
            continuous_palette=palette,
            color_filter=color_filter,
            discrete_label_colors=discrete_label_colors,
        )
        with open(out_path, "wb") as fh:
            fh.write(png)
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("save pairplot failed")
        return jsonify(error=str(err)), 500

    return jsonify(ok=True, path=out_path, filename=filename)


# ---------------------------------------------------------------------------
# Trajectory routes
# ---------------------------------------------------------------------------


def api_default_trajectory_endpoints():
    """Start/end in plot space along the long axis of the point cloud (see helper)."""
    e: DashboardExperiment = g.dashboard_exp
    err = _trajectory_eligibility_error(e)
    if err is not None:
        return err
    xcol = request.args.get("x", "")
    ycol = request.args.get("y", "")
    if xcol not in e.plot_df.columns or ycol not in e.plot_df.columns:
        return jsonify(error="bad axis column"), 400
    try:
        validate_trajectory_plot_axes(e, xcol, ycol)
        start, end = default_trajectory_endpoints_xy(e, xcol, ycol)
    except ValueError as err:
        return jsonify(error=str(err)), 400

    return jsonify(ok=True, start=start, end=end)


def trajectory_creator_page():
    """Trajectory creator UI (requires CUDA + weights like the volume explorer)."""
    e: DashboardExperiment = g.dashboard_exp
    if not explorer_volumes_eligible(e):
        return (
            render_template(
                "no_images.html",
                reason=(
                    "Trajectory creator needs a CUDA GPU and model weights for the "
                    "current epoch."
                ),
            ),
            200,
        )
    zdim = int(e.z.shape[1])
    traj_cols = trajectory_plot_axis_columns(e)
    color_cols = e.numeric_columns
    dx, dy = trajectory_default_xy_cols(traj_cols, zdim)
    cov_keys = list(dict.fromkeys(traj_cols + color_cols))

    return render_template(
        "trajectory_creator.html",
        traj_axis_cols=traj_cols,
        numeric_cols=color_cols,
        covariate_display_map=_covariate_display_map(cov_keys),
        default_x=dx,
        default_y=dy,
        zdim=zdim,
        exp_workdir=e.workdir,
        chimerax_cpus_default=DEFAULT_CHIMERAX_PARALLEL,
    )


def api_trajectory_save_zpath():
    """Write z-path text to a server-side path (defaults to ``workdir/z-path.txt``)."""
    e: DashboardExperiment = g.dashboard_exp
    err = _trajectory_eligibility_error(e)
    if err is not None:
        return err
    data = _request_json_dict()
    txt = data.get("z_path_txt")
    if not isinstance(txt, str):
        return jsonify(error="z_path_txt must be a string"), 400
    raw_out_path = data.get("out_path")
    if raw_out_path is None or str(raw_out_path).strip() == "":
        out_path = os.path.join(os.path.abspath(e.workdir), "z-path.txt")
    else:
        if not isinstance(raw_out_path, str):
            return jsonify(error="out_path must be a string path"), 400
        req_path = raw_out_path.strip()
        if not req_path:
            return jsonify(error="out_path must not be empty"), 400
        if not req_path.lower().endswith(".txt"):
            req_path = req_path + ".txt"
        if os.path.isabs(req_path):
            out_path = os.path.abspath(req_path)
        else:
            out_path = os.path.abspath(os.path.join(e.workdir, req_path))
    out_dir = os.path.dirname(out_path) or os.path.abspath(e.workdir)
    try:
        os.makedirs(out_dir, exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(txt)
    except OSError as err:
        return jsonify(error=str(err)), 500

    return jsonify(ok=True, path=out_path)


def api_trajectory_save_volumes():
    """Save cached trajectory ``.mrc`` files into a chosen server-side folder."""
    e: DashboardExperiment = g.dashboard_exp
    err = _trajectory_eligibility_error(e)
    if err is not None:
        return err
    data = _request_json_dict()
    token = str(data.get("volume_cache_id", "") or "").strip()
    out_dir = str(data.get("out_dir", "") or "").strip()

    if not token:
        return jsonify(error="Missing volume_cache_id. Generate volumes first."), 400
    if not out_dir:
        return jsonify(error="Choose an output folder."), 400

    try:
        saved = save_cached_volumes_to_dir(
            token, out_dir, filename_prefix="trajectory_volume"
        )
        return jsonify(
            ok=True, out_dir=os.path.abspath(out_dir), files=saved, n_saved=len(saved)
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except OSError as err:
        return jsonify(error=str(err)), 500


def _trajectory_anchor_driven_json(
    e: DashboardExperiment, anchor_indices: list[int], data: dict
):
    """Shared tail for k-means/random/import anchor endpoints."""
    xcol, ycol = trajectory_axes_from_payload(e, data)
    mode, n_points, max_neighbors, avg_neighbors = trajectory_anchor_mode_params(data)
    return jsonify(
        trajectory_anchor_payload_from_indices(
            e,
            anchor_indices,
            xcol,
            ycol,
            mode=mode,
            n_points=n_points,
            max_neighbors=max_neighbors,
            avg_neighbors=avg_neighbors,
        )
    )


def api_trajectory_import_anchors():
    """Import anchor indices from a server-side ``.txt`` path."""
    e: DashboardExperiment = g.dashboard_exp

    err = _trajectory_eligibility_error(e)
    if err is not None:
        return err

    data = _request_json_dict()
    server_path = str(data.get("server_path", "") or "").strip()
    if not server_path:
        return jsonify(error="no file path provided"), 400
    if not server_path.lower().endswith(".txt"):
        return jsonify(error="select a .txt file"), 400
    abs_path = os.path.abspath(server_path)
    if not os.path.isfile(abs_path):
        return jsonify(error="file not found on server"), 400

    try:
        from pathlib import Path as _Path

        anchor_indices = parse_anchor_indices_txt(_Path(abs_path).read_bytes())
        return _trajectory_anchor_driven_json(e, anchor_indices, data)
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("trajectory anchor import failed")
        return jsonify(error=str(err)), 500


def api_list_server_files():
    """List directories and ``.txt`` files for the server-side file browser."""
    e: DashboardExperiment = g.dashboard_exp
    req_dir = request.args.get("dir", "").strip()
    root = os.path.abspath(e.workdir)
    if req_dir:
        browse = os.path.abspath(req_dir)
    else:
        analyze_dir = os.path.join(root, f"analyze.{e.epoch}")
        browse = analyze_dir if os.path.isdir(analyze_dir) else root
    if not os.path.isdir(browse):
        return jsonify(error="directory not found"), 400
    entries: list[dict] = []
    try:
        with os.scandir(browse) as it:
            for entry in it:
                name = entry.name
                if entry.is_dir():
                    entries.append({"name": name, "type": "dir"})
                elif name.lower().endswith(".txt"):
                    entries.append({"name": name, "type": "file"})
        entries.sort(key=lambda row: row["name"])
    except PermissionError:
        return jsonify(error="permission denied"), 403
    parent = os.path.dirname(browse) if browse != "/" else None

    return jsonify(ok=True, dir=browse, parent=parent, entries=entries)


def api_trajectory_kmeans_centers():
    """Load ``kmeansK/centers_ind.txt`` for the current ``analyze.N``."""
    e: DashboardExperiment = g.dashboard_exp
    err = _trajectory_eligibility_error(e)
    if err is not None:
        return err

    try:
        return _trajectory_anchor_driven_json(
            e, load_kmeans_center_indices(e), _request_json_dict()
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("trajectory k-means centers failed")
        return jsonify(error=str(err)), 500


def api_trajectory_random_indices():
    """Choose up to 10 random dataset indices (fewer if the stack is smaller)."""
    e: DashboardExperiment = g.dashboard_exp
    err = _trajectory_eligibility_error(e)
    if err is not None:
        return err

    try:
        return _trajectory_anchor_driven_json(
            e, random_dataset_indices(e, 10), _request_json_dict()
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("trajectory random indices failed")
        return jsonify(error=str(err)), 500


def api_trajectory_coords():
    """Latent z along the trajectory line (no ChimeraX / volume decode)."""
    e: DashboardExperiment = g.dashboard_exp
    err = _trajectory_eligibility_error(e)
    if err is not None:
        return err

    try:
        p = parse_trajectory_request_body(e, _request_json_dict())
    except ValueError as err:
        return jsonify(error=str(err)), 400
    try:
        z_traj, traj_rows, traj_xy = compute_trajectory_latent_path(e, p)
        payload = trajectory_shared_json_payload(
            e,
            z_traj,
            traj_rows,
            traj_xy,
            mode=p["mode"],
            n_points=p["n_points"],
            xcol=p["xcol"],
            ycol=p["ycol"],
        )
        _add_direct_anchor_pidx(payload, p, z_traj)
        return jsonify(payload)
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("trajectory coords failed")
        return jsonify(error=str(err)), 500


def api_trajectory_volumes():
    e: DashboardExperiment = g.dashboard_exp
    err = _trajectory_eligibility_error(e)
    if err is not None:
        return err

    try:
        data = _request_json_dict()
        p = parse_trajectory_request_body(e, data)
    except ValueError as err:
        return jsonify(error=str(err)), 400
    try:
        z_traj, traj_rows, traj_xy = compute_trajectory_latent_path(e, p)
        cc = int(data.get("chimerax_cpus", DEFAULT_CHIMERAX_PARALLEL))
        blobs, cache_token = generate_trajectory_volume_pngs(
            e, z_traj, chimerax_cpus=cc
        )
        payload = trajectory_shared_json_payload(
            e,
            z_traj,
            traj_rows,
            traj_xy,
            mode=p["mode"],
            n_points=p["n_points"],
            xcol=p["xcol"],
            ycol=p["ycol"],
        )
        _add_direct_anchor_pidx(payload, p, z_traj)
        payload["images"] = [
            base64.standard_b64encode(b).decode("ascii") for b in blobs
        ]
        payload["volume_cache_id"] = cache_token
        if traj_rows is not None and p["mode"] in ("nearest", "graph"):
            payload["particle_thumbs"] = [
                particle_thumbnail_b64_from_row(e, int(r)) for r in traj_rows
            ]
        return jsonify(payload)
    except EnvironmentError as err:
        return jsonify(error=str(err), need_chimerax=True), 503
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("trajectory volume generation failed")
        return jsonify(error=str(err)), 500


# ---------------------------------------------------------------------------
# Landscape volume PCA (analyze_landscape outputs)
# ---------------------------------------------------------------------------


def landscape_volpca_page():
    """Volume PCA / sketch landscape view (``analyze_landscape`` outputs)."""
    e: DashboardExperiment = g.dashboard_exp
    if not list_landscape_epochs(e.workdir):
        return render_template(
            "landscape_volpca.html",
            no_landscape=True,
            reason=(
                "No landscape.N folders in this output directory — run "
                "<code>cryodrgn analyze_landscape</code> first."
            ),
            exp_epoch=int(e.epoch),
        )
    if not landscape_analysis_ready(e.workdir, e.epoch):
        return render_template(
            "landscape_volpca.html",
            no_landscape=True,
            reason=(
                f"No complete landscape output for <strong>epoch {e.epoch}</strong>. "
                "Switch epoch in the nav or run "
                "<code>cryodrgn analyze_landscape</code> for this epoch."
            ),
            exp_epoch=int(e.epoch),
        )

    return render_template(
        "landscape_volpca.html",
        no_landscape=False,
        reason="",
        exp_epoch=int(e.epoch),
    )


def api_landscape_volpca_meta():
    e: DashboardExperiment = g.dashboard_exp

    try:
        return jsonify(meta_for_api(e))
    except FileNotFoundError as err:
        return jsonify(ok=False, error=str(err)), 400
    except Exception as err:
        logger.exception("landscape volpca meta failed")
        return jsonify(ok=False, error=str(err)), 500


def api_landscape_volpca_scatter():
    e: DashboardExperiment = g.dashboard_exp
    epochs = list_landscape_epochs(e.workdir)
    if not epochs:
        return jsonify(error="No landscape outputs for this workdir."), 400
    le = int(e.epoch)
    if le not in epochs:
        return (
            jsonify(
                error=(
                    f"No landscape.{le} for the current dashboard epoch. "
                    "Switch epoch in the nav or run analyze_landscape."
                ),
            ),
            400,
        )
    ax_x_q = (request.args.get("axis_x") or "").strip()
    ax_y_q = (request.args.get("axis_y") or "").strip()
    if ax_x_q and ax_y_q:
        axis_x = ax_x_q
        axis_y = ax_y_q
    elif ax_x_q or ax_y_q:
        return (
            jsonify(
                error="Provide both axis_x and axis_y (e.g. pc:0 and umap:1), "
                "or use legacy pc_x and pc_y.",
            ),
            400,
        )
    else:
        axis_x = f"pc:{int(request.args.get('pc_x', '0'))}"
        axis_y = f"pc:{int(request.args.get('pc_y', '1'))}"
    color = (request.args.get("color") or "none").strip().lower()
    landscape_dir = landscape_dir_for_epoch(e.workdir, le)
    if color not in ("none", "state") and color not in e.numeric_columns:
        return jsonify(error="bad color column"), 400
    if color == "state":
        if load_sketch_state_labels(landscape_dir) is None:
            return (
                jsonify(
                    error=(
                        "Agglomerative state coloring is not available for this landscape."
                    ),
                ),
                400,
            )
    try:
        js = landscape_volpca_scatter_json(
            landscape_dir,
            e,
            axis_x=axis_x,
            axis_y=axis_y,
            color_mode=color,
            continuous_palette=request.args.get("palette"),
        )
    except FileNotFoundError as err:
        return jsonify(error=str(err)), 400
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("landscape volpca scatter failed")
        return jsonify(error=str(err)), 500

    return Response(js, mimetype="application/json")


def api_landscape_volpca_generate_animations():
    e: DashboardExperiment = g.dashboard_exp
    epochs = list_landscape_epochs(e.workdir)
    if not epochs:
        return jsonify(error="No landscape outputs for this workdir."), 400
    data = _request_json_dict()
    le = int(e.epoch)
    if le not in epochs:
        return jsonify(error=f"No folder landscape.{le} for the current epoch."), 400
    vol_raw = data.get("vol_indices") or data.get("volumes")
    if not isinstance(vol_raw, list) or not vol_raw:
        return jsonify(error="Provide vol_indices (integers)."), 400
    try:
        vol_indices = sorted({int(v) for v in vol_raw})
    except (TypeError, ValueError):
        return jsonify(error="vol_indices must be integers."), 400
    mode = str(data.get("mode") or "cycle").strip().lower()
    gf = int(data.get("gif_frames", LANDSCAPE_SKETCH_GIF_FRAMES_DEFAULT))
    cc = int(data.get("chimerax_cpus", DEFAULT_CHIMERAX_PARALLEL))
    cf = int(data.get("cycle_frames_per_vol", 8))
    plot_color_raw = (
        str(data.get("color_mode") or data.get("color") or "none").strip().lower()
    )
    cm = "state" if plot_color_raw == "state" else "none"
    palette = data.get("palette")
    if palette is not None:
        palette = str(palette).strip() or None
    landscape_dir = landscape_dir_for_epoch(e.workdir, le)
    if cm == "state" and load_sketch_state_labels(landscape_dir) is None:
        cm = "none"
    pcm = plot_color_raw
    if pcm == "state" and load_sketch_state_labels(landscape_dir) is None:
        pcm = "none"
    reuse_tok = data.get("reuse_rotate_keyframes_token")
    if reuse_tok is not None:
        reuse_tok = str(reuse_tok).strip() or None
    view_rotations = data.get("view_rotations")
    if view_rotations is not None and not isinstance(view_rotations, dict):
        return (
            jsonify(error="view_rotations must be an object with x, y, and z degrees."),
            400,
        )

    try:
        t0 = time.perf_counter()
        token, _files, rendered_vol_indices = generate_landscape_volume_animations(
            landscape_dir,
            vol_indices,
            exp=e,
            mode=mode,
            gif_frames=gf,
            chimerax_cpus=cc,
            cycle_frames_per_vol=cf,
            color_mode=cm,
            plot_color_mode=pcm,
            continuous_palette=palette,
            reuse_rotate_keyframes_token=reuse_tok,
            view_rotations=view_rotations,
        )
        items = animation_payload_b64(token)
        view_matrix = next(
            (it.get("view_matrix") for it in items if it.get("view_matrix")),
            None,
        )
        elapsed_s = time.perf_counter() - t0
        return jsonify(
            ok=True,
            token=token,
            items=items,
            landscape_epoch=le,
            duration_s=int(round(elapsed_s)),
            rendered_vol_indices=rendered_vol_indices,
            batch_mode=mode,
            view_matrix=view_matrix,
        )
    except EnvironmentError as err:
        return jsonify(error=str(err), need_chimerax=True), 503
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("landscape volpca animation generate failed")
        return jsonify(error=str(err)), 500


def api_landscape_volpca_save_animations():
    """Persist generated GIF bytes from ``generate_animations`` to disk (JSON: token, out_dir)."""
    e: DashboardExperiment = g.dashboard_exp
    data = _request_json_dict()
    token = data.get("token")
    if not token or not isinstance(token, str):
        return jsonify(error="Missing token."), 400
    le = int(e.epoch)
    out_dir = data.get("out_dir")
    if out_dir is not None and not isinstance(out_dir, str):
        return jsonify(error="out_dir must be a string or omitted."), 400

    try:
        paths = save_landscape_animations(
            token,
            out_dir if isinstance(out_dir, str) else None,
            exp=e,
            landscape_epoch=le,
        )
        return jsonify(ok=True, paths=paths)
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("landscape volpca save animations failed")
        return jsonify(error=str(err)), 500


# ---------------------------------------------------------------------------
# App factory + route table
# ---------------------------------------------------------------------------


_ROUTES = (
    ("/api/set_epoch", api_set_epoch, ("POST",)),
    ("/api/set_workdir", api_set_workdir, ("POST",)),
    ("/", index, ("GET",)),
    ("/command-builder", command_builder_page, ("GET",)),
    ("/abinit-builder", abinit_builder_redirect, ("GET",)),
    ("/filter", filter_page_redirect, ("GET",)),
    ("/api/save_selection", api_save_selection, ("POST",)),
    ("/api/covariate_threshold_rows", api_covariate_threshold_rows, ("POST",)),
    ("/api/covariate_legend_context", api_covariate_legend_context, ("POST",)),
    ("/explorer", explorer, ("GET",)),
    ("/api/explorer_volume_media", api_explorer_volume_media, ("POST",)),
    ("/api/scatter", api_scatter, ("GET",)),
    ("/latent-3d", latent_3d_page, ("GET",)),
    ("/api/scatter3d_z", api_scatter3d_z, ("GET", "POST")),
    ("/api/latent3d_preview.png", api_latent3d_preview_png, ("GET",)),
    ("/api/preview_montage", api_preview_montage, ("GET",)),
    ("/api/preload_images", api_preload_images, ("GET", "POST")),
    ("/pairplot", pairplot_page, ("GET",)),
    ("/api/pairplot", api_pairplot, ("POST",)),
    ("/api/save_pairplot_png", api_save_pairplot_png, ("POST",)),
    ("/trajectory", trajectory_creator_page, ("GET",)),
    ("/api/trajectory_volumes", api_trajectory_volumes, ("POST",)),
    ("/api/trajectory_coords", api_trajectory_coords, ("POST",)),
    ("/api/trajectory_save_zpath", api_trajectory_save_zpath, ("POST",)),
    ("/api/trajectory_save_volumes", api_trajectory_save_volumes, ("POST",)),
    ("/api/trajectory_import_anchors", api_trajectory_import_anchors, ("POST",)),
    ("/api/list_server_files", api_list_server_files, ("GET",)),
    ("/api/trajectory_kmeans_centers", api_trajectory_kmeans_centers, ("POST",)),
    ("/api/trajectory_random_indices", api_trajectory_random_indices, ("POST",)),
    ("/api/default_trajectory_endpoints", api_default_trajectory_endpoints, ("GET",)),
    ("/landscape-volpca", landscape_volpca_page, ("GET",)),
    ("/api/landscape_volpca/meta", api_landscape_volpca_meta, ("GET",)),
    ("/api/landscape_volpca/scatter", api_landscape_volpca_scatter, ("GET",)),
    (
        "/api/landscape_volpca/generate_animations",
        api_landscape_volpca_generate_animations,
        ("POST",),
    ),
    (
        "/api/landscape_volpca/save_animations",
        api_landscape_volpca_save_animations,
        ("POST",),
    ),
)


def create_app(
    workdir: str | None,
    epoch: int = -1,
    kmeans: int = -1,
    filter_plot_inds: str | None = None,
    cpus: int = 4,
) -> Flask:
    """Construct the Flask app: config, discovery mode, ``before_request``, and all routes."""
    app = Flask(
        __name__,
        template_folder=_TEMPLATE_DIR,
        static_folder=_STATIC_DIR,
        static_url_path="/static",
    )
    app.jinja_env.filters["abbrev_middle"] = abbrev_middle
    app.config["PRELOAD_CPUS"] = max(1, cpus)
    app.secret_key = os.environ.get(
        "CRYODRGN_DASHBOARD_SECRET", "cryodrgn-dashboard-dev-key"
    )
    # Invalidates Flask session epoch when the dashboard process restarts so CLI
    # ``--epoch`` matches the first load (see :func:`resolve_epoch`).
    app.config["DASHBOARD_SESSION_BOOT_ID"] = str(uuid.uuid4())
    command_builder_only = workdir is None
    app.config["COMMAND_BUILDER_ONLY"] = command_builder_only
    app.config["DASHBOARD_DISCOVERY_CWD"] = os.getcwd()
    discovered = discover_cryodrgn_workdirs(os.getcwd()) if command_builder_only else []
    app.config["DASHBOARD_DISCOVERED_WORKDIRS"] = discovered
    app.config["DASHBOARD_DISCOVERY_BOOT_ID"] = (
        str(uuid.uuid4()) if command_builder_only and discovered else None
    )
    app.before_request(bind_dashboard_exp)

    if command_builder_only:
        app.config["DASHBOARD_WORKDIR"] = None
        app.config["DASHBOARD_KMEANS"] = -1
        app.config["DASHBOARD_EPOCHS"] = [0]
        app.config["DASHBOARD_START_EPOCH"] = 0
        app.config["FILTER_PLOT_INDS"] = None
        app.context_processor(inject_meta_command_builder_only)
    else:
        workdir = os.path.abspath(workdir)
        epochs = list_z_epochs(workdir)
        if not epochs:
            raise ValueError(
                f"No analyzed epochs under {workdir!r} — need z.N.pkl and analyze.N/ "
                "(run `cryodrgn analyze` first)."
            )
        start_epoch = epoch if epoch != -1 else max(epochs)
        if start_epoch not in epochs:
            start_epoch = max(epochs)
        app.config["DASHBOARD_WORKDIR"] = workdir
        app.config["DASHBOARD_KMEANS"] = kmeans
        app.config["DASHBOARD_EPOCHS"] = epochs
        app.config["DASHBOARD_START_EPOCH"] = start_epoch
        app.config["FILTER_PLOT_INDS"] = filter_plot_inds
        app.context_processor(inject_meta)

    for rule, view_func, methods in _ROUTES:
        app.add_url_rule(rule, view_func=view_func, methods=list(methods))

    return app


def run_server(
    workdir: str | None,
    epoch: int = -1,
    kmeans: int = -1,
    plot_inds: str | None = None,
    host: str = "127.0.0.1",
    port: int = 5050,
    debug: bool = False,
    cpus: int = 4,
) -> None:
    """Entry point for ``cryodrgn dashboard``: threaded Werkzeug development server."""
    app = create_app(
        workdir=workdir,
        epoch=epoch,
        kmeans=kmeans,
        filter_plot_inds=plot_inds,
        cpus=cpus,
    )

    app.run(host=host, port=port, debug=debug, threaded=True)
