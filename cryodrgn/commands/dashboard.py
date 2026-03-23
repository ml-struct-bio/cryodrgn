"""Launch a local web dashboard for cryoDRGN interactive analyses.

The dashboard opens in your browser with a particle explorer (scatter, selection,
optional thumbnails), pairwise latent panels, a 3D latent scatter, and a command builder.

Example
-------
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


def add_args(parser: argparse.ArgumentParser) -> None:
    import os

    parser.add_argument(
        "outdir",
        type=os.path.abspath,
        help="experiment output directory (same as ``cryodrgn filter``)",
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
        "--particle",
        "--filter",
        action="store_true",
        dest="particle_selection",
        help="open the particle explorer instead of the launch menu",
    )
    view.add_argument(
        "--image-viewer",
        "--scatter",
        action="store_true",
        dest="image_viewer",
        help="open the image viewer instead of the launch menu",
    )
    view.add_argument(
        "--pair-grid",
        action="store_true",
        help="open the pairwise latent grid instead of the launch menu",
    )
    view.add_argument(
        "--three-dimensional",
        action="store_true",
        help="open the 3D latent scatter instead of the launch menu",
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


def main(args: argparse.Namespace) -> None:
    import os

    if args.filter_max_points is not None:
        os.environ["CRYODRGN_DASHBOARD_FILTER_MAX_POINTS"] = str(args.filter_max_points)

    outdir = args.outdir
    _VIEW_PATHS = {
        "particle_selection": "/explorer",
        "image_viewer": "/explorer",
        "pair_grid": "/pairplot",
        "three_dimensional": "/latent-3d",
    }
    initial_path = next(
        (p for flag, p in _VIEW_PATHS.items() if getattr(args, flag, False)),
        "/",
    )

    if not args.no_browser:

        def _open() -> None:
            webbrowser.open(f"http://{args.host}:{args.port}{initial_path}")

        Timer(1.0, _open).start()

    from cryodrgn.dashboard.app import run_server

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
