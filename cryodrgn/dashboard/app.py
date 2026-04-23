"""Flask app for the cryoDRGN analysis dashboard.

Responsibilities are split across small sibling modules:

* :mod:`cryodrgn.dashboard.context` — workdir/epoch resolution + template
  context injectors + long-lived caches.
* :mod:`cryodrgn.dashboard.trajectory` — pure-logic helpers for the trajectory
  creator (axis validation, anchor parsing, graph traversal, JSON payloads).
* :mod:`cryodrgn.dashboard.preload` — thumbnail sampling / encoding for the
  explorer preview montage and hover pre-load cache.

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
from flask import (
    Flask,
    Response,
    current_app,
    g,
    jsonify,
    redirect,
    render_template,
    request,
    url_for,
)

from cryodrgn.dashboard.context import (
    PRELOAD_CACHE,
    abbrev_middle,
    active_workdir,
    api_set_epoch,
    api_set_workdir,
    bind_dashboard_exp,
    command_builder_template_kwargs,
    discover_cryodrgn_workdirs,
    inject_meta,
    inject_meta_command_builder_only,
)
from cryodrgn.dashboard.data import DashboardExperiment, list_z_epochs
from cryodrgn.dashboard.explorer_volumes import (
    DEFAULT_CHIMERAX_PARALLEL,
    DEFAULT_GIF_FRAMES,
    explorer_volumes_eligible,
    generate_montage_volume_pngs,
    generate_trajectory_volume_pngs,
    save_cached_volumes_to_dir,
    volume_cell_gif_from_cache,
)
from cryodrgn.dashboard.plots import (
    normalize_continuous_palette,
    pair_grid_png,
    pair_grid_skeleton_placeholder_layout,
    scatter3d_z_json,
    scatter3d_z_preview_png,
    scatter_json,
)
from cryodrgn.dashboard.preload import (
    encode_particle_batch,
    format_preload_cache_time_hint,
    load_plot_df_rows_from_plot_inds_file,
    montage_bytes,
    particle_thumbnail_b64_from_row,
    sample_plot_df_rows_for_preload,
)
from cryodrgn.dashboard.trajectory import (
    compute_trajectory_latent_path,
    default_trajectory_endpoints_xy,
    direct_anchor_particle_indices_payload,
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

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_TEMPLATE_DIR = os.path.join(_THIS_DIR, "templates")
_STATIC_DIR = os.path.join(_THIS_DIR, "static")

_TRAJECTORY_INELIGIBLE_MSG = (
    "Trajectory creator needs a CUDA GPU, single-particle data, "
    "and weights for the current epoch."
)
_EXPLORER_VOLUMES_INELIGIBLE_MSG = (
    "Volume explorer needs a CUDA GPU, single-particle data, "
    "and weights for the current epoch."
)


# ---------------------------------------------------------------------------
# Small request / axis / covariate helpers used only by routes in this file
# ---------------------------------------------------------------------------


def _request_json_dict() -> dict:
    """Return request JSON payload as a mapping, defaulting to ``{}``."""
    data = request.get_json(force=True, silent=True)
    return data if isinstance(data, dict) else {}


def _filter_ui_scatter_max_points() -> int:
    """Cap for ``/api/scatter`` when ``filter_ui=1`` (env override available)."""
    try:
        v = int(os.environ.get("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", "500000"))
    except ValueError:
        v = 500_000
    return max(50_000, min(v, 2_000_000))


def _default_xy_cols(cols: list[str]) -> tuple[str, str]:
    """Pick sensible default X/Y axes (UMAP if available, else first two)."""
    x = "UMAP1" if "UMAP1" in cols else cols[0]
    y = "UMAP2" if "UMAP2" in cols else cols[min(1, len(cols) - 1)]
    return x, y


def _covariate_display_name(name: str) -> str:
    """Human-friendly covariate names in dashboard selectors."""
    if name == "labels":
        return "k-means labels"
    return name


def _parse_preselect_rows_param(raw: str | None) -> tuple[list[int] | None, str | None]:
    """Parse ``preselect_rows`` query (comma-separated ints). Returns ``(rows, err)``."""
    s = (raw or "").strip()
    if not s:
        return None, None
    try:
        return [int(p) for p in s.split(",") if p.strip()][:5000], None
    except ValueError as exc:
        return None, f"invalid preselect_rows: {exc}"


def _redirect(endpoint: str):
    return redirect(url_for(endpoint), code=302)


def _trajectory_eligibility_error(e: DashboardExperiment):
    """JSON 400 when the trajectory creator is not available, else ``None``."""
    if explorer_volumes_eligible(e):
        return None
    return jsonify(error=_TRAJECTORY_INELIGIBLE_MSG), 400


def _parse_pairplot_request(
    e: DashboardExperiment, payload: dict
) -> tuple[str, str, str, str]:
    color_col = payload.get("color_col") or payload.get("lower_color_col")
    if not color_col or not isinstance(color_col, str):
        raise ValueError("Choose a color covariate.")
    if color_col not in e.plot_df.columns:
        raise ValueError("Invalid color column.")
    if color_col in {f"z{i}" for i in range(int(e.z.shape[1]))}:
        raise ValueError("Latent z columns cannot be used as the color covariate.")
    raw_diag = payload.get("diagonal_emb")
    if raw_diag is None or (isinstance(raw_diag, str) and raw_diag.strip() == ""):
        diagonal_emb = "umap" if has_umap_columns(e) else "pc"
    else:
        diagonal_emb = str(raw_diag).lower()
    upper_style = (payload.get("upper_style") or "scatter").lower()
    if diagonal_emb not in ("pc", "umap"):
        raise ValueError("diagonal_emb must be pc or umap.")
    if upper_style not in ("scatter", "hex"):
        raise ValueError("upper_style must be scatter or hex.")
    if diagonal_emb == "umap" and not has_umap_columns(e):
        raise ValueError("UMAP is not available for this run.")
    if diagonal_emb == "pc" and not has_pc_columns(e):
        raise ValueError("PCA components are not available.")
    raw_palette = payload.get("palette")
    pair_palette = normalize_continuous_palette(
        str(raw_palette) if raw_palette is not None else None,
    )
    return color_col, diagonal_emb, upper_style, pair_palette


def _add_direct_anchor_pidx(payload: dict, p: dict, z_traj: np.ndarray) -> None:
    """Merge direct-anchor particle IDs into ``payload`` when applicable."""
    if not (p.get("use_anchors") and p["mode"] == "direct"):
        return
    pidx = direct_anchor_particle_indices_payload(
        anchor_indices=p["anchor_indices"],
        interpolation_points=p["n_points"],
        n_total=int(np.asarray(z_traj).shape[0]),
    )
    if pidx is not None:
        payload["traj_particle_indices"] = pidx


# ---------------------------------------------------------------------------
# Top-level pages
# ---------------------------------------------------------------------------


def index():
    if not active_workdir(current_app):
        return render_template(
            "index.html",
            can_images=False,
            zdim=0,
            show_trajectory_creator=False,
            command_builder_only=True,
        )
    e: DashboardExperiment = g.dashboard_exp
    return render_template(
        "index.html",
        can_images=e.can_preview_particles,
        zdim=int(e.z.shape[1]),
        show_trajectory_creator=explorer_volumes_eligible(e),
        command_builder_only=False,
    )


def command_builder_page():
    e: DashboardExperiment | None = None
    if active_workdir(current_app):
        e = g.dashboard_exp
    return render_template(
        "command_builder.html",
        **command_builder_template_kwargs(e),
    )


def abinit_builder_redirect():
    return _redirect("command_builder_page")


def filter_page_redirect():
    return _redirect("explorer")


def api_save_selection():
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


# ---------------------------------------------------------------------------
# Explorer + scatter routes
# ---------------------------------------------------------------------------


def explorer():
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
    return render_template(
        "scatter_explorer.html",
        numeric_cols=cols,
        covariate_display_map={c: _covariate_display_name(c) for c in cols},
        default_x=dx,
        default_y=dy,
        initial_rows=initial_rows,
        total_particles=int(len(e.all_indices)),
        workdir=e.workdir,
        preload_cache_time_hint=format_preload_cache_time_hint(pc),
        show_volume_explorer=explorer_volumes_eligible(e),
    )


def api_explorer_volume_media():
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
            blobs, cache_token = generate_montage_volume_pngs(e, rows)
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
    e: DashboardExperiment = g.dashboard_exp
    xcol = request.args.get("x", e.numeric_columns[0])
    ycol = request.args.get("y", e.numeric_columns[0])
    ccol = request.args.get("color") or "none"
    filter_ui = request.args.get("filter_ui") == "1"
    full = request.args.get("full") == "1"
    if xcol not in e.plot_df.columns or ycol not in e.plot_df.columns:
        return jsonify(error="bad axis column"), 400
    if ccol != "none" and ccol not in e.plot_df.columns:
        return jsonify(error="bad color column"), 400
    if filter_ui:
        max_pts = _filter_ui_scatter_max_points()
    elif full:
        max_pts = None
    else:
        max_pts = 200_000
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
        covariate_display_map={c: _covariate_display_name(c) for c in cols},
        default_x="z0",
        default_y="z1",
        default_z="z2",
    )


def api_scatter3d_z():
    e: DashboardExperiment = g.dashboard_exp
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
    """Return a stratified subsample of particle images as base64 JPEGs.

    Selection: ~2/3 random, ~1/6 high mean distance to reference points, ~1/6
    high nearest-neighbor distance. Use **POST** with a JSON body when
    ``selected_rows`` is large (lasso selections); query strings hit proxy /
    server URI length limits (414).
    """
    e: DashboardExperiment = g.dashboard_exp
    if not e.can_preview_particles:
        return jsonify(rows=[], images=[], elapsed=0)

    cols = e.plot_df.columns
    restrict_list: list[int] | None = None

    if request.method == "POST":
        data = _request_json_dict()
        xcol = str(data.get("x") or "")
        ycol = str(data.get("y") or "")
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
    else:
        xcol = request.args.get("x", "")
        ycol = request.args.get("y", "")
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
    sel_key: tuple[int, ...] | None = (
        tuple(sorted(set(restrict_list))) if restrict_list else None
    )
    key = (e.epoch, e.kmeans_folder_id, xcol, ycol, sel_key)
    if key in PRELOAD_CACHE:
        rows, imgs, elapsed = PRELOAD_CACHE[key]
        return jsonify(rows=rows, images=imgs, elapsed=elapsed)

    t0 = time.monotonic()
    cpus = int(current_app.config.get("PRELOAD_CPUS") or 4)
    rows, global_indices = sample_plot_df_rows_for_preload(
        e, xcol, ycol, restrict_to_rows=restrict_list
    )
    if not global_indices:
        return jsonify(rows=[], images=[], elapsed=0.0)

    if cpus > 1 and len(global_indices) > cpus:
        from concurrent.futures import ProcessPoolExecutor

        chunk_sz = -(-len(global_indices) // cpus)
        chunks = [
            global_indices[i : i + chunk_sz]
            for i in range(0, len(global_indices), chunk_sz)
        ]
        with ProcessPoolExecutor(max_workers=len(chunks)) as pool:
            futures = [
                pool.submit(encode_particle_batch, e.particles_path, e.datadir, ch, 96)
                for ch in chunks
            ]
            imgs: list[str] = []
            for f in futures:
                imgs.extend(f.result())
    else:
        imgs = encode_particle_batch(e.particles_path, e.datadir, global_indices, 96)

    elapsed = round(time.monotonic() - t0, 1)
    PRELOAD_CACHE[key] = (rows, imgs, elapsed)
    return jsonify(rows=rows, images=imgs, elapsed=elapsed)


# ---------------------------------------------------------------------------
# Pair-grid routes
# ---------------------------------------------------------------------------


def pairplot_page():
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
    z_names = {f"z{i}" for i in range(zdim)}
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
        covariate_display_map={c: _covariate_display_name(c) for c in color_choices},
        default_color=default_color,
        has_umap=has_umap_columns(e),
        zdim=zdim,
        skeleton_placeholder_cells=pair_grid_skeleton_placeholder_layout(zdim),
        pairplot_save_default_name="zdim_pairplot.png",
        pairplot_save_default_dir=os.path.join(e.workdir, f"analyze.{e.epoch}"),
    )


def api_pairplot():
    e: DashboardExperiment = g.dashboard_exp
    try:
        color_col, diagonal_emb, upper_style, palette = _parse_pairplot_request(
            e, _request_json_dict()
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
    e: DashboardExperiment = g.dashboard_exp
    payload = _request_json_dict()
    try:
        color_col, diagonal_emb, upper_style, palette = _parse_pairplot_request(
            e, payload
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
        covariate_display_map={c: _covariate_display_name(c) for c in cov_keys},
        default_x=dx,
        default_y=dy,
        zdim=zdim,
        exp_workdir=e.workdir,
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
        for name in sorted(os.listdir(browse)):
            full = os.path.join(browse, name)
            if os.path.isdir(full):
                entries.append({"name": name, "type": "dir"})
            elif name.lower().endswith(".txt"):
                entries.append({"name": name, "type": "file"})
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
        p = parse_trajectory_request_body(e, _request_json_dict())
    except ValueError as err:
        return jsonify(error=str(err)), 400
    try:
        z_traj, traj_rows, traj_xy = compute_trajectory_latent_path(e, p)
        blobs, cache_token = generate_trajectory_volume_pngs(e, z_traj)
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
    ("/explorer", explorer, ("GET",)),
    ("/api/explorer_volume_media", api_explorer_volume_media, ("POST",)),
    ("/api/scatter", api_scatter, ("GET",)),
    ("/latent-3d", latent_3d_page, ("GET",)),
    ("/api/scatter3d_z", api_scatter3d_z, ("GET",)),
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
)


def create_app(
    workdir: str | None,
    epoch: int = -1,
    kmeans: int = -1,
    filter_plot_inds: str | None = None,
    cpus: int = 4,
) -> Flask:
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
    app = create_app(
        workdir=workdir,
        epoch=epoch,
        kmeans=kmeans,
        filter_plot_inds=plot_inds,
        cpus=cpus,
    )
    app.run(host=host, port=port, debug=debug, threaded=True)
