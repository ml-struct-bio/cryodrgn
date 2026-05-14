"""Flask app for the cryoDRGN analysis dashboard.

Route handlers live in :mod:`cryodrgn.dashboard.routes_explorer` and
:mod:`cryodrgn.dashboard.routes_analysis`; this module defines the app factory,
static/template paths, and URL wiring.
"""

from __future__ import annotations

import logging
import os
import uuid

from flask import Flask

from cryodrgn.dashboard.context import (
    abbrev_middle,
    api_set_epoch,
    api_set_workdir,
    bind_dashboard_exp,
    discover_cryodrgn_workdirs,
    inject_meta,
    inject_meta_command_builder_only,
)
from cryodrgn.dashboard.data import discover_analyzed_workdirs, list_z_epochs

# Scatter caps: re-exported for tests that use ``from cryodrgn.dashboard import app`` (logic in route_helpers).
from cryodrgn.dashboard.route_helpers import (
    _filter_ui_scatter_max_points,
    _particle_explorer_scatter_cap_from_env,
    _particle_explorer_scatter_max_points,
)

__all__ = [
    "create_app",
    "run_server",
    "_ROUTES",
    "_filter_ui_scatter_max_points",
    "_particle_explorer_scatter_cap_from_env",
    "_particle_explorer_scatter_max_points",
]
from cryodrgn.dashboard import routes_analysis as ra
from cryodrgn.dashboard import routes_explorer as re

logger = logging.getLogger(__name__)

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_TEMPLATE_DIR = os.path.join(_THIS_DIR, "templates")
_STATIC_DIR = os.path.join(_THIS_DIR, "static")

# App factory + route table
# ---------------------------------------------------------------------------


_ROUTES = (
    ("/api/set_epoch", api_set_epoch, ("POST",)),
    ("/api/set_workdir", api_set_workdir, ("POST",)),
    ("/", re.index, ("GET",)),
    ("/command-builder", re.command_builder_page, ("GET",)),
    ("/abinit-builder", re.abinit_builder_redirect, ("GET",)),
    ("/filter", re.filter_page_redirect, ("GET",)),
    ("/api/save_selection", re.api_save_selection, ("POST",)),
    ("/api/covariate_threshold_rows", re.api_covariate_threshold_rows, ("POST",)),
    ("/api/covariate_legend_context", re.api_covariate_legend_context, ("POST",)),
    ("/explorer", re.explorer, ("GET",)),
    ("/api/explorer_volume_media", re.api_explorer_volume_media, ("POST",)),
    ("/api/scatter", re.api_scatter, ("GET",)),
    ("/latent-3d", re.latent_3d_page, ("GET",)),
    ("/landscape-full-3d", re.landscape_full_3d_page, ("GET",)),
    ("/api/scatter3d_z", re.api_scatter3d_z, ("GET", "POST")),
    (
        "/api/scatter3d_z_landscape_full",
        re.api_scatter3d_z_landscape_full,
        ("GET", "POST"),
    ),
    ("/api/latent3d_preview.png", re.api_latent3d_preview_png, ("GET",)),
    ("/api/preview_montage", re.api_preview_montage, ("GET",)),
    ("/api/preload_images", re.api_preload_images, ("GET", "POST")),
    ("/pairplot", ra.pairplot_page, ("GET",)),
    ("/api/pairplot", ra.api_pairplot, ("POST",)),
    ("/api/save_pairplot_png", ra.api_save_pairplot_png, ("POST",)),
    ("/trajectory", ra.trajectory_creator_page, ("GET",)),
    ("/api/trajectory_volumes", ra.api_trajectory_volumes, ("POST",)),
    ("/api/trajectory_coords", ra.api_trajectory_coords, ("POST",)),
    ("/api/trajectory_save_zpath", ra.api_trajectory_save_zpath, ("POST",)),
    ("/api/trajectory_save_volumes", ra.api_trajectory_save_volumes, ("POST",)),
    ("/api/trajectory_import_anchors", ra.api_trajectory_import_anchors, ("POST",)),
    ("/api/list_server_files", ra.api_list_server_files, ("GET",)),
    ("/api/trajectory_kmeans_centers", ra.api_trajectory_kmeans_centers, ("POST",)),
    ("/api/trajectory_random_indices", ra.api_trajectory_random_indices, ("POST",)),
    (
        "/api/default_trajectory_endpoints",
        ra.api_default_trajectory_endpoints,
        ("GET",),
    ),
    ("/landscape-volpca", ra.landscape_volpca_page, ("GET",)),
    ("/api/landscape_volpca/meta", ra.api_landscape_volpca_meta, ("GET",)),
    ("/api/landscape_volpca/scatter", ra.api_landscape_volpca_scatter, ("GET",)),
    (
        "/api/landscape_volpca/generate_animations",
        ra.api_landscape_volpca_generate_animations,
        ("POST",),
    ),
    (
        "/api/landscape_volpca/save_animations",
        ra.api_landscape_volpca_save_animations,
        ("POST",),
    ),
)


def create_app(
    workdir: str | None,
    epoch: int = -1,
    kmeans: int = -1,
    filter_plot_inds: str | None = None,
    cpus: int = 4,
    discovery_root: str | None = None,
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

    if command_builder_only:
        base = (
            os.path.abspath(discovery_root)
            if discovery_root is not None
            else os.getcwd()
        )
        app.config["DASHBOARD_DISCOVERY_CWD"] = base
        discovered = (
            discover_analyzed_workdirs(base)
            if discovery_root is not None
            else discover_cryodrgn_workdirs(base)
        )
        if discovery_root is not None and not discovered:
            raise ValueError(
                f"No cryoDRGN analyzed outputs found under {base!r}. Expected a run "
                "folder with z.N.pkl and analyze.N/, or a directory whose immediate "
                "subfolders are such outputs (after `cryodrgn analyze`)."
            )
    else:
        app.config["DASHBOARD_DISCOVERY_CWD"] = os.getcwd()
        discovered = []
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
    discovery_root: str | None = None,
) -> None:
    """Entry point for ``cryodrgn dashboard``: threaded Werkzeug development server."""
    app = create_app(
        workdir=workdir,
        epoch=epoch,
        kmeans=kmeans,
        filter_plot_inds=plot_inds,
        cpus=cpus,
        discovery_root=discovery_root,
    )

    app.run(host=host, port=port, debug=debug, threaded=True)
