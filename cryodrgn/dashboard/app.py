"""Flask app for the cryoDRGN analysis dashboard."""

from __future__ import annotations

import base64
import io
import logging
import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
from flask import (
    Flask,
    Response,
    current_app,
    g,
    jsonify,
    render_template,
    request,
    session,
)

from cryodrgn import utils
from cryodrgn.dashboard.data import (
    DashboardExperiment,
    list_z_epochs,
    load_experiment,
    particle_image_array,
)
from cryodrgn.dashboard.mpl_style import ezlab_matplotlib_rc
from cryodrgn.dashboard.plots import (
    pair_grid_png,
    pair_grid_skeleton_placeholder_layout,
    scatter3d_z_json,
    scatter_json,
)

logger = logging.getLogger(__name__)

_TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")

# (epoch, kmeans) -> loaded experiment (invalidated when epoch changes in session).
_EXP_CACHE: dict[tuple[int, int], DashboardExperiment] = {}


def _filter_ui_scatter_max_points() -> int:
    """Max points for the filter scatter (random subsample, re-read per request)."""
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


def _montage_bytes(exp: DashboardExperiment, row_indices: list[int]) -> bytes:
    rows = [int(r) for r in row_indices[:25]]
    with ezlab_matplotlib_rc():
        if not rows:
            fig, ax = plt.subplots(figsize=(4, 4))
            ax.text(
                0.5, 0.5, "Hover points to preview particles", ha="center", va="center"
            )
            ax.axis("off")
            buf = io.BytesIO()
            fig.savefig(buf, format="png", bbox_inches="tight")
            plt.close(fig)
            buf.seek(0)
            return buf.getvalue()

        imgs = [particle_image_array(exp, r) for r in rows]
        n = len(imgs)
        ncol = int(np.ceil(n**0.5))
        nrow = int(np.ceil(n / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(1.8 * ncol, 1.8 * nrow))
        axes_flat = np.atleast_1d(axes).ravel()
        for ax in axes_flat[n:]:
            ax.axis("off")
        for i, img in enumerate(imgs):
            ax = axes_flat[i]
            lo, hi = np.percentile(img, (2, 98))
            disp = np.clip((img - lo) / (hi - lo + 1e-9), 0, 1)
            ax.imshow(disp, cmap="gray")
            ax.set_title(f"idx {int(exp.plot_df.iloc[rows[i]]['index'])}", fontsize=8)
            ax.axis("off")
        plt.tight_layout()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight", dpi=120)
        plt.close(fig)
        buf.seek(0)
        return buf.getvalue()


def _plot_df_rows_for_dataset_indices(
    exp: DashboardExperiment, dataset_indices: np.ndarray
) -> list[int]:
    """Map saved dataset indices (values in ``indices.pkl``) to ``plot_df`` row positions."""
    want = {int(x) for x in np.asarray(dataset_indices).ravel()}
    ai = np.asarray(exp.all_indices)
    return [int(i) for i in range(len(ai)) if int(ai[i]) in want]


def _resolve_epoch(app: Flask) -> int:
    epochs: list[int] = app.config["DASHBOARD_EPOCHS"]
    if not epochs:
        raise RuntimeError("No z.N.pkl epochs in workdir.")
    sess = session.get("dashboard_epoch")
    if sess is None:
        return int(app.config["DASHBOARD_START_EPOCH"])
    try:
        ep = int(sess)
    except (TypeError, ValueError):
        return int(app.config["DASHBOARD_START_EPOCH"])
    if ep not in epochs:
        return max(epochs)
    return ep


def _get_dashboard_exp(app: Flask) -> DashboardExperiment:
    ep = _resolve_epoch(app)
    km = int(app.config["DASHBOARD_KMEANS"])
    key = (ep, km)
    if key not in _EXP_CACHE:
        _EXP_CACHE[key] = load_experiment(
            app.config["DASHBOARD_WORKDIR"],
            epoch=ep,
            kmeans=km,
        )
    return _EXP_CACHE[key]


def create_app(
    workdir: str,
    epoch: int = -1,
    kmeans: int = -1,
    filter_plot_inds: str | None = None,
) -> Flask:
    app = Flask(__name__, template_folder=_TEMPLATE_DIR)
    app.secret_key = os.environ.get(
        "CRYODRGN_DASHBOARD_SECRET", "cryodrgn-dashboard-dev-key"
    )
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

    @app.before_request
    def _bind_dashboard_exp() -> None:
        g.dashboard_exp = _get_dashboard_exp(current_app)

    @app.context_processor
    def inject_meta():
        e: DashboardExperiment = g.dashboard_exp
        return {
            "exp_workdir": e.workdir,
            "exp_epoch": e.epoch,
            "exp_kmeans": e.kmeans_folder_id,
            "filter_plot_inds_default": app.config["FILTER_PLOT_INDS"] or "",
            "dashboard_epochs": app.config["DASHBOARD_EPOCHS"],
        }

    @app.post("/api/set_epoch")
    def api_set_epoch():
        data = request.get_json(force=True, silent=True) or {}
        try:
            ep = int(data.get("epoch"))
        except (TypeError, ValueError):
            return jsonify(error="Invalid epoch."), 400
        if ep not in current_app.config["DASHBOARD_EPOCHS"]:
            return jsonify(error="Epoch not available for this output folder."), 400
        session["dashboard_epoch"] = ep
        _EXP_CACHE.clear()
        return jsonify(ok=True, epoch=ep)

    @app.get("/")
    def index():
        e: DashboardExperiment = g.dashboard_exp
        return render_template(
            "index.html",
            can_images=e.can_preview_particles,
            zdim=int(e.z.shape[1]),
        )

    @app.get("/filter")
    def filter_page():
        e: DashboardExperiment = g.dashboard_exp
        initial_rows: list[int] = []
        plot_inds_path = app.config.get("FILTER_PLOT_INDS") or ""
        if plot_inds_path.strip() and os.path.isfile(plot_inds_path):
            try:
                pre = utils.load_pkl(plot_inds_path)
                initial_rows = _plot_df_rows_for_dataset_indices(e, np.asarray(pre))
            except Exception as err:
                logger.warning("Could not load --plot-inds for filter page: %s", err)
        cols = e.numeric_columns
        dx, dy = _default_xy_cols(cols)
        return render_template(
            "filter.html",
            numeric_cols=cols,
            default_x=dx,
            default_y=dy,
            initial_rows=initial_rows,
            workdir=e.workdir,
        )

    @app.post("/api/save_selection")
    def api_save_selection():
        e: DashboardExperiment = g.dashboard_exp
        data = request.get_json(force=True, silent=True) or {}
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
        basename = (data.get("basename") or "indices").strip() or "indices"
        basename = os.path.basename(basename)
        if force:
            out_base = os.path.join(e.workdir, "indices")
        else:
            sel_dir = (data.get("sel_dir") or "").strip()
            dir_abs = (
                os.path.abspath(sel_dir) if sel_dir else os.path.abspath(e.workdir)
            )
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
            with open(path_inv, "wb") as fh:
                pickle.dump(inverse_ds, fh)
        except OSError as err:
            logger.exception("save selection failed")
            return jsonify(error=str(err)), 500
        return jsonify(
            ok=True,
            path=path_main,
            inverse_path=path_inv,
            n_selected=int(selected_ds.size),
        )

    @app.get("/explorer")
    def explorer():
        e: DashboardExperiment = g.dashboard_exp
        if not e.can_preview_particles:
            return (
                render_template(
                    "no_images.html",
                    reason="Tilt-series data does not include single-particle thumbnails in this view.",
                ),
                200,
            )
        cols = e.numeric_columns
        dx, dy = _default_xy_cols(cols)
        return render_template(
            "scatter_explorer.html",
            numeric_cols=cols,
            default_x=dx,
            default_y=dy,
        )

    @app.get("/api/scatter")
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
            return jsonify(error="bad colour column"), 400
        if filter_ui:
            max_pts = _filter_ui_scatter_max_points()
        elif full:
            max_pts = None
        else:
            max_pts = 120_000
        preselect_rows: list[int] | None = None
        if filter_ui:
            raw_pre = (request.args.get("preselect_rows") or "").strip()
            if raw_pre:
                try:
                    preselect_rows = [int(p) for p in raw_pre.split(",") if p.strip()][
                        :5000
                    ]
                except ValueError as exc:
                    return jsonify(error=f"invalid preselect_rows: {exc}"), 400
        try:
            js = scatter_json(
                e,
                xcol,
                ycol,
                None if ccol == "none" else ccol,
                max_points=max_pts,
                preselect_plot_df_rows=preselect_rows,
                # Scattergl + Plotly.react often leaves the returned promise pending; canvas Scatter
                # matches the filter UI and clears reliably (see plotly_afterplot in templates).
                use_webgl=False,
            )
        except Exception as err:
            logger.exception("scatter plot failed")
            return jsonify(error=str(err)), 500
        return Response(js, mimetype="application/json")

    @app.get("/latent-3d")
    def latent_3d_page():
        e: DashboardExperiment = g.dashboard_exp
        zdim = int(e.z.shape[1])
        if zdim < 3:
            return (
                render_template(
                    "pair_grid_need_more_cols.html",
                    kind="z3",
                    n=zdim,
                ),
                200,
            )
        z_cols = [f"z{i}" for i in range(zdim)]
        cols = e.numeric_columns
        return render_template(
            "latent_3d.html",
            z_cols=z_cols,
            numeric_cols=cols,
            default_x="z0",
            default_y="z1",
            default_z="z2",
        )

    @app.get("/api/scatter3d_z")
    def api_scatter3d_z():
        e: DashboardExperiment = g.dashboard_exp
        xcol = request.args.get("x", "z0")
        ycol = request.args.get("y", "z1")
        zcol = request.args.get("z", "z2")
        ccol = request.args.get("color") or "none"
        if ccol != "none" and ccol not in e.plot_df.columns:
            return jsonify(error="bad colour column"), 400
        try:
            js = scatter3d_z_json(
                e,
                xcol,
                ycol,
                zcol,
                None if ccol == "none" else ccol,
            )
        except ValueError as err:
            return jsonify(error=str(err)), 400
        except Exception as err:
            logger.exception("3D latent scatter failed")
            return jsonify(error=str(err)), 500
        return Response(js, mimetype="application/json")

    @app.get("/api/preview_montage")
    def api_preview_montage():
        e: DashboardExperiment = g.dashboard_exp
        raw = request.args.get("rows", "")
        if not raw.strip():
            data = _montage_bytes(e, [])
        else:
            parts = [p.strip() for p in raw.split(",") if p.strip()]
            try:
                idxs = [int(p) for p in parts]
            except ValueError:
                return jsonify(error="rows must be integers"), 400
            data = _montage_bytes(e, idxs)
        return Response(data, mimetype="image/png")

    @app.get("/pairplot")
    def pairplot_page():
        e: DashboardExperiment = g.dashboard_exp
        zdim = int(e.z.shape[1])
        if zdim < 2:
            return (
                render_template(
                    "pair_grid_need_more_cols.html",
                    kind="zdim",
                    n=zdim,
                ),
                200,
            )
        has_pc = "PC1" in e.plot_df.columns and "PC2" in e.plot_df.columns
        has_umap = (
            e.umap is not None
            and "UMAP1" in e.plot_df.columns
            and "UMAP2" in e.plot_df.columns
        )
        if not has_pc:
            return (
                render_template(
                    "pair_grid_need_more_cols.html",
                    kind="pca",
                    n=zdim,
                ),
                200,
            )
        z_names = {f"z{i}" for i in range(zdim)}
        color_choices = [c for c in e.numeric_columns if c not in z_names]
        if not color_choices:
            return (
                render_template(
                    "pair_grid_need_more_cols.html",
                    kind="numeric",
                    n=0,
                ),
                200,
            )
        default_color = next(
            (c for c in ("labels", "znorm", "UMAP1", "PC1") if c in color_choices),
            color_choices[0],
        )
        skeleton_cells = pair_grid_skeleton_placeholder_layout(zdim)
        return render_template(
            "pair_grid.html",
            color_choices=color_choices,
            default_color=default_color,
            has_umap=has_umap,
            zdim=zdim,
            skeleton_placeholder_cells=skeleton_cells,
        )

    @app.post("/api/pairplot")
    def api_pairplot():
        e: DashboardExperiment = g.dashboard_exp
        payload = request.get_json(force=True, silent=True) or {}
        color_col = payload.get("color_col") or payload.get("lower_color_col")
        if not color_col or not isinstance(color_col, str):
            return jsonify(error="Choose a colour covariate."), 400
        if color_col not in e.plot_df.columns:
            return jsonify(error="Invalid colour column."), 400
        if color_col in {f"z{i}" for i in range(int(e.z.shape[1]))}:
            return (
                jsonify(
                    error="Latent z columns cannot be used as the colour covariate."
                ),
                400,
            )
        raw_diag = payload.get("diagonal_emb")
        if raw_diag is None or (isinstance(raw_diag, str) and raw_diag.strip() == ""):
            diagonal_emb = (
                "umap"
                if e.umap is not None
                and "UMAP1" in e.plot_df.columns
                and "UMAP2" in e.plot_df.columns
                else "pc"
            )
        else:
            diagonal_emb = str(raw_diag).lower()
        upper_style = (payload.get("upper_style") or "scatter").lower()
        if diagonal_emb not in ("pc", "umap"):
            return jsonify(error="diagonal_emb must be pc or umap."), 400
        if upper_style not in ("scatter", "hex"):
            return jsonify(error="upper_style must be scatter or hex."), 400
        if diagonal_emb == "umap" and (
            e.umap is None
            or "UMAP1" not in e.plot_df.columns
            or "UMAP2" not in e.plot_df.columns
        ):
            return jsonify(error="UMAP is not available for this run."), 400
        if diagonal_emb == "pc" and (
            "PC1" not in e.plot_df.columns or "PC2" not in e.plot_df.columns
        ):
            return jsonify(error="PCA components are not available."), 400
        try:
            png, cells = pair_grid_png(
                e,
                lower_color_col=color_col,
                diagonal_emb=diagonal_emb,
                upper_style=upper_style,
            )
        except ValueError as err:
            return jsonify(error=str(err)), 400
        except Exception as err:
            logger.exception("pair grid failed")
            return jsonify(error=str(err)), 500
        b64 = base64.standard_b64encode(png).decode("ascii")
        return jsonify(png_b64=b64, cells=cells)

    return app


def run_server(
    workdir: str,
    epoch: int = -1,
    kmeans: int = -1,
    plot_inds: str | None = None,
    host: str = "127.0.0.1",
    port: int = 5050,
    debug: bool = False,
) -> None:
    app = create_app(
        workdir=workdir,
        epoch=epoch,
        kmeans=kmeans,
        filter_plot_inds=plot_inds,
    )
    app.run(host=host, port=port, debug=debug, threaded=True)
