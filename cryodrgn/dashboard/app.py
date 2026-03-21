"""Flask app for the cryoDRGN analysis dashboard."""

from __future__ import annotations

import io
import logging
import os
import subprocess
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from flask import Flask, Response, jsonify, render_template, request

from cryodrgn.dashboard.data import (
    DashboardExperiment,
    load_experiment,
    particle_image_array,
)
from cryodrgn.dashboard.plots import pair_grid_json, scatter_json

matplotlib.use("Agg")
logger = logging.getLogger(__name__)

_TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")


def _montage_bytes(exp: DashboardExperiment, row_indices: list[int]) -> bytes:
    rows = [int(r) for r in row_indices[:25]]
    if not rows:
        fig, ax = plt.subplots(figsize=(4, 4))
        ax.text(0.5, 0.5, "Hover points to preview particles", ha="center", va="center")
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


def create_app(
    exp: DashboardExperiment,
    filter_plot_inds: str | None = None,
) -> Flask:
    app = Flask(__name__, template_folder=_TEMPLATE_DIR)
    app.config["DASHBOARD_EXP"] = exp
    app.config["FILTER_PLOT_INDS"] = filter_plot_inds

    @app.context_processor
    def inject_meta():
        e: DashboardExperiment = app.config["DASHBOARD_EXP"]
        return {
            "exp_workdir": e.workdir,
            "exp_epoch": e.epoch,
            "exp_kmeans": e.kmeans_folder_id,
            "filter_plot_inds_default": app.config["FILTER_PLOT_INDS"] or "",
        }

    @app.get("/")
    def index():
        e: DashboardExperiment = app.config["DASHBOARD_EXP"]
        return render_template(
            "index.html",
            can_images=e.can_preview_particles,
        )

    @app.get("/filter")
    def filter_page():
        return render_template("filter_launch.html")

    @app.post("/api/launch-filter")
    def api_launch_filter():
        e: DashboardExperiment = app.config["DASHBOARD_EXP"]
        force = request.form.get("force") == "on"
        sel_dir = request.form.get("sel_dir") or None
        plot_inds = request.form.get("plot_inds") or None
        if sel_dir == "":
            sel_dir = None
        if plot_inds == "":
            plot_inds = None

        script = f"""from argparse import Namespace
from cryodrgn.commands.filter import main
args = Namespace(
    outdir={e.workdir!r},
    epoch={e.epoch},
    kmeans={e.kmeans_folder_id},
    force={force!r},
    plot_inds={plot_inds!r},
    sel_dir={sel_dir!r},
)
main(args)
"""
        subprocess.Popen([sys.executable, "-c", script], close_fds=True)
        return jsonify(ok=True, message="Desktop filter window should open shortly.")

    @app.get("/explorer")
    def explorer():
        e: DashboardExperiment = app.config["DASHBOARD_EXP"]
        if not e.can_preview_particles:
            return (
                render_template(
                    "no_images.html",
                    reason="Tilt-series experiments do not support particle thumbnails here.",
                ),
                200,
            )
        cols = e.numeric_columns
        default_x = "UMAP1" if "UMAP1" in cols else cols[0]
        default_y = "UMAP2" if "UMAP2" in cols else cols[min(1, len(cols) - 1)]
        return render_template(
            "scatter_explorer.html",
            numeric_cols=cols,
            default_x=default_x,
            default_y=default_y,
        )

    @app.get("/api/scatter")
    def api_scatter():
        e: DashboardExperiment = app.config["DASHBOARD_EXP"]
        xcol = request.args.get("x", e.numeric_columns[0])
        ycol = request.args.get("y", e.numeric_columns[0])
        ccol = request.args.get("color") or "none"
        if xcol not in e.plot_df.columns or ycol not in e.plot_df.columns:
            return jsonify(error="bad axis column"), 400
        if ccol != "none" and ccol not in e.plot_df.columns:
            return jsonify(error="bad colour column"), 400
        try:
            js = scatter_json(e, xcol, ycol, None if ccol == "none" else ccol)
        except Exception as err:
            logger.exception("scatter plot failed")
            return jsonify(error=str(err)), 500
        return Response(js, mimetype="application/json")

    @app.get("/api/preview_montage")
    def api_preview_montage():
        e: DashboardExperiment = app.config["DASHBOARD_EXP"]
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
        e: DashboardExperiment = app.config["DASHBOARD_EXP"]
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
        color_choices = e.numeric_columns
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
        return render_template(
            "pair_grid.html",
            color_choices=color_choices,
            default_color=default_color,
            has_umap=has_umap,
            zdim=zdim,
        )

    @app.post("/api/pairplot")
    def api_pairplot():
        e: DashboardExperiment = app.config["DASHBOARD_EXP"]
        payload = request.get_json(force=True, silent=True) or {}
        color_col = payload.get("color_col") or payload.get("lower_color_col")
        if not color_col or not isinstance(color_col, str):
            return jsonify(error="Choose a colour covariate."), 400
        if color_col not in e.plot_df.columns:
            return jsonify(error="Invalid colour column."), 400
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
            js = pair_grid_json(
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
        return Response(js, mimetype="application/json")

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
    exp = load_experiment(workdir, epoch=epoch, kmeans=kmeans)
    app = create_app(exp, filter_plot_inds=plot_inds)
    app.run(host=host, port=port, debug=debug, threaded=True)
