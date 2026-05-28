"""HTTP routes: pair plot, trajectory creator, landscape vol-PCA APIs."""

from __future__ import annotations

import base64
import logging
import os
import time

from flask import (
    Response,
    g,
    jsonify,
    render_template,
    request,
)

from cryodrgn.dashboard.context import _request_json_dict
from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.landscape_volpca import (
    LANDSCAPE_SKETCH_GIF_FRAMES_DEFAULT,
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
from cryodrgn.dashboard.particle_explorer import (
    DEFAULT_CHIMERAX_PARALLEL,
    explorer_volumes_eligible,
    generate_trajectory_volume_pngs,
    save_cached_volumes_to_dir,
)
from cryodrgn.dashboard.preload import particle_thumbnail_b64_from_row
from cryodrgn.dashboard.plots import (
    pair_grid_figure_aspect_ratio,
    pair_grid_margin_fractions_for_js,
    pair_grid_png,
    pair_grid_skeleton_placeholder_layout,
)
from cryodrgn.dashboard.route_helpers import (
    _add_direct_anchor_pidx,
    _covariate_display_map,
    _parse_color_filter_for_column,
    _parse_optional_discrete_label_colors,
    _parse_pairplot_request,
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

# ---------------------------------------------------------------------------
# Pair-grid routes
# ---------------------------------------------------------------------------


def pairplot_page():
    """Pair-grid UI: latent vs latent with chosen color covariate.

    Needs PCA + numeric color column.
    """
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
    """Render the z-dimensional pair grid to an in-memory PNG.

    Returns base64 + cell layout JSON.
    """
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
        data = _request_json_dict()
        p = parse_trajectory_request_body(e, data)
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
            color_col=str(data.get("color") or "none"),
            continuous_palette=data.get("palette"),
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
            color_col=str(data.get("color") or "none"),
            continuous_palette=data.get("palette"),
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
                        "Agglomerative state coloring is not available "
                        "for this landscape."
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
    plot_color_raw = str(data.get("color_mode") or data.get("color") or "none").strip()
    pcm_lower = plot_color_raw.lower()
    if pcm_lower == "none":
        pcm = "none"
    elif pcm_lower == "state":
        pcm = "state"
    else:
        pcm = plot_color_raw
    cm = "state" if pcm_lower == "state" else "none"
    palette = data.get("palette")
    if palette is not None:
        palette = str(palette).strip() or None
    landscape_dir = landscape_dir_for_epoch(e.workdir, le)
    if cm == "state" and load_sketch_state_labels(landscape_dir) is None:
        cm = "none"
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
    view_matrix = data.get("view_matrix")
    if view_matrix is not None and not isinstance(view_matrix, str):
        return jsonify(error="view_matrix must be a string."), 400
    if view_matrix is not None:
        view_matrix = str(view_matrix).strip() or None

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
            view_matrix=view_matrix,
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
    """Persist generated GIF bytes from ``generate_animations`` to disk.

    JSON body: token, out_dir.
    """
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
