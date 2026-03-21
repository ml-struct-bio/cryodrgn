"""Launch a local web dashboard for cryoDRGN interactive analyses.

The dashboard opens in your browser with three entry points:

1. **Desktop filter** — spawns the classic ``cryodrgn filter`` Tk/matplotlib UI.
2. **Scatter explorer** — Plotly scatter with particle thumbnails on hover (SPA only).
3. **Pair grid** — Interactive latent-style grid (as on the ``new_plots`` branch) with
   optional custom X/Y for all lower-triangle panels.

Example
-------
$ cryodrgn dashboard 00_trainvae
$ cryodrgn dashboard 00_trainvae --epoch 30 --port 8080
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
        help="train epoch to load (default: latest z.N.pkl)",
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
        help="optional path to indices.pkl (passed through to desktop filter only)",
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
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Flask debug mode (reloads; not for production)",
    )


def main(args: argparse.Namespace) -> None:
    outdir = args.outdir
    if not args.no_browser:

        def _open() -> None:
            webbrowser.open(f"http://{args.host}:{args.port}/")

        Timer(1.0, _open).start()

    from cryodrgn.dashboard.app import run_server

    logger.info(
        "Starting dashboard for %s (epoch=%s) at http://%s:%s/",
        outdir,
        args.epoch,
        args.host,
        args.port,
    )
    run_server(
        workdir=outdir,
        epoch=args.epoch,
        kmeans=args.kmeans,
        plot_inds=args.plot_inds,
        host=args.host,
        port=args.port,
        debug=args.debug,
    )
