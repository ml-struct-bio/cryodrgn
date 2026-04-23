"""Launch a local web dashboard for cryoDRGN interactive analyses.

The dashboard opens in your browser with a particle explorer (scatter, selection,
optional thumbnails), pairwise latent panels, the 3-D Latent Space Visualizer, and a command builder.

Example
-------
$ cryodrgn dashboard
$ cryodrgn dashboard 00_trainvae
$ cryodrgn dashboard 00_trainvae --epoch 30 --port 8080
$ cryodrgn dashboard 00_trainvae --image-viewer
$ cryodrgn dashboard 00_trainvae --particle --filter-max 10000
"""

from __future__ import annotations

import argparse
import logging
import webbrowser
from threading import Timer

logger = logging.getLogger(__name__)


def _configure_dashboard_logging(verbosity: int) -> None:
    """Configure dashboard-related logger levels from ``-v`` count.

    Levels:
    - 0 (default): quiet; suppress werkzeug internal/request chatter.
    - 1 (``-v``): dashboard info logs.
    - 2 (``-vv``): dashboard debug logs.
    - 3+ (``-vvv``): include werkzeug debug logs.
    """

    verbosity = max(0, int(verbosity or 0))
    if verbosity <= 0:
        dashboard_level = logging.WARNING
        web_level = logging.WARNING
    elif verbosity == 1:
        dashboard_level = logging.INFO
        web_level = logging.WARNING
    elif verbosity == 2:
        dashboard_level = logging.DEBUG
        web_level = logging.INFO
    else:
        dashboard_level = logging.DEBUG
        web_level = logging.DEBUG

    root = logging.getLogger()
    if not root.handlers:
        logging.basicConfig(
            level=min(dashboard_level, web_level),
            format="%(levelname)s:%(name)s:%(message)s",
        )

    logging.getLogger("cryodrgn.dashboard").setLevel(dashboard_level)
    logging.getLogger(__name__).setLevel(dashboard_level)
    # Default suppression target requested by user.
    logging.getLogger("werkzeug").setLevel(web_level)
    logging.getLogger("werkzeug._internal").setLevel(web_level)


def add_args(parser: argparse.ArgumentParser) -> None:
    import os

    parser.add_argument(
        "outdir",
        nargs="?",
        default=None,
        type=os.path.abspath,
        help="experiment output directory (optional; omit for command-builder-only launch mode)",
    )
    parser.add_argument(
        "--epoch",
        "-e",
        type=int,
        default=-1,
        help="train epoch to load (default: latest epoch with analyze.N/ output)",
    )
    parser.add_argument(
        "--kmeans",
        "-k",
        type=int,
        default=-1,
        help="k-means folder kmeansK to use (default: first found under analyze.*)",
    )
    parser.add_argument(
        "--plot-inds",
        type=str,
        default=None,
        help="optional path to indices.pkl to pre-highlight on the particle explorer",
    )
    parser.add_argument(
        "--filter-max-points",
        "--filter-max",
        type=int,
        default=None,
        metavar="N",
        dest="filter_max_points",
        help=(
            "max particles drawn in the web filter scatter (default 500000; "
            "clamped 50k–2M; same as env CRYODRGN_DASHBOARD_FILTER_MAX_POINTS)"
        ),
    )
    parser.add_argument(
        "--host",
        type=str,
        default="127.0.0.1",
        help="bind address (default: %(default)s)",
    )
    parser.add_argument(
        "--port",
        "-p",
        type=int,
        default=5050,
        help="port for the dashboard server (default: %(default)s)",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="do not open a browser tab automatically",
    )
    view = parser.add_mutually_exclusive_group()
    view.add_argument(
        "--particle-selection",
        "--filter",
        "--image-viewer",
        action="store_true",
        dest="particle_selection",
        help="open the particle explorer instead of the launch menu",
    )
    view.add_argument(
        "--pair-grid",
        action="store_true",
        help="open the pairwise latent grid instead of the launch menu",
    )
    view.add_argument(
        "--three-dimensional",
        action="store_true",
        help="open the 3-D Latent Space Visualizer instead of the launch menu",
    )
    view.add_argument(
        "--command-builder",
        action="store_true",
        help="open the command builder instead of the launch menu",
    )
    view.add_argument(
        "--trajectory-creator",
        action="store_true",
        help="open the trajectory creator instead of the launch menu",
    )
    parser.add_argument(
        "--cpus",
        "-c",
        type=int,
        default=4,
        help="CPU cores for image cache generation (default: %(default)s)",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Flask debug mode (reloads; not for production)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="count",
        default=0,
        help=(
            "increase dashboard log verbosity (repeat: -v info, -vv debug, "
            "-vvv includes werkzeug internals)"
        ),
    )


def main(args: argparse.Namespace) -> None:
    import os

    _configure_dashboard_logging(args.verbose)

    if args.filter_max_points is not None:
        os.environ["CRYODRGN_DASHBOARD_FILTER_MAX_POINTS"] = str(args.filter_max_points)

    outdir = args.outdir
    command_builder_only = outdir is None

    view_paths = {
        "particle_selection": "/explorer",
        "pair_grid": "/pairplot",
        "three_dimensional": "/latent-3d",
        "command_builder": "/command-builder",
        "trajectory_creator": "/trajectory",
    }
    experiment_views = set(view_paths) - {"command_builder"}
    if command_builder_only and any(getattr(args, f, False) for f in experiment_views):
        raise ValueError(
            "Dashboard views that require experiment data need an output directory. "
            "Use `cryodrgn dashboard <outdir>` for explorer/pair-grid/3D/trajectory."
        )
    initial_path = next(
        (p for flag, p in view_paths.items() if getattr(args, flag, False)),
        "/",
    )

    if not args.no_browser:

        def _open() -> None:
            webbrowser.open(f"http://{args.host}:{args.port}{initial_path}")

        Timer(1.0, _open).start()

    from cryodrgn.dashboard.app import run_server

    if command_builder_only:
        logger.info(
            "Starting dashboard in command-builder-only mode at http://%s:%s%s",
            args.host,
            args.port,
            initial_path,
        )
    else:
        logger.info(
            "Starting dashboard for %s (epoch=%s) at http://%s:%s%s",
            outdir,
            args.epoch,
            args.host,
            args.port,
            initial_path,
        )
    run_server(
        workdir=outdir,
        epoch=args.epoch,
        kmeans=args.kmeans,
        plot_inds=args.plot_inds,
        host=args.host,
        port=args.port,
        debug=args.debug,
        cpus=args.cpus,
    )
