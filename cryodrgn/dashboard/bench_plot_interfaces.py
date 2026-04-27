#!/usr/bin/env python3
"""Time initial plot payloads for each dashboard analysis view (same paths as the Flask app).

Usage::
    python -m cryodrgn.dashboard.bench_plot_interfaces /path/to/train_outdir

Uses CRYODRGN_DASHBOARD_FILTER_MAX_POINTS like ``/api/scatter?filter_ui=1``.
"""
from __future__ import annotations

import argparse
import cProfile
import io
import os
import pstats
import sys
import time


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("workdir", help="training output directory (with analyze.N/)")
    p.add_argument(
        "--profile-pairplot",
        action="store_true",
        help="print cProfile top functions for pair_grid_png only",
    )
    args = p.parse_args()

    workdir = os.path.abspath(args.workdir)
    os.environ.setdefault("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", "500000")

    from cryodrgn.dashboard.data import list_z_epochs, load_experiment
    from cryodrgn.dashboard.plots import (
        pair_grid_png,
        scatter3d_z_json,
        scatter_json,
    )

    def pairplot_color_and_diagonal(exp):
        """Default pair-grid colour column and diagonal (UMAP vs PC) for ``exp``."""
        cols = exp.numeric_columns
        color = next(
            (
                c
                for c in ("labels", "znorm", "UMAP1", "PC1")
                if c in exp.plot_df.columns
            ),
            cols[0],
        )
        diag = (
            "umap" if exp.umap is not None and "UMAP1" in exp.plot_df.columns else "pc"
        )
        return color, diag

    epochs = list_z_epochs(workdir)
    if not epochs:
        print("No epochs with analyze.N/ found.", file=sys.stderr)
        return 1
    ep = max(epochs)

    def bench(name: str, fn) -> None:
        t0 = time.perf_counter()
        out = fn()
        dt = time.perf_counter() - t0
        extra = ""
        if isinstance(out, str):
            extra = f"  JSON {len(out) / 1e6:.2f} MB"
        elif isinstance(out, tuple) and out and isinstance(out[0], bytes):
            extra = f"  PNG {len(out[0]) / 1e6:.2f} MB"
        print(f"  {name:42s} {dt:7.3f}s{extra}")

    def run_round(label: str) -> None:
        print(f"\n=== {label} ===")
        t0 = time.perf_counter()
        exp = load_experiment(workdir, ep, kmeans=-1)
        print(
            f"load_experiment(epoch={ep}): "
            f"{time.perf_counter() - t0:.3f}s  "
            f"(rows={len(exp.plot_df):,}, zdim={int(exp.z.shape[1])})"
        )

        cols = exp.numeric_columns
        x = "UMAP1" if "UMAP1" in cols else cols[0]
        y = "UMAP2" if "UMAP2" in cols else cols[min(1, len(cols) - 1)]

        bench(
            "Scatter (explorer, 120k cap, canvas)",
            lambda: scatter_json(exp, x, y, None, 120_000, use_webgl=False),
        )
        max_f = int(
            os.environ.get(
                "CRYODRGN_DASHBOARD_FILTER_MAX_POINTS",
                "500000",
            )
        )
        bench(
            f"Scatter (filter_ui=1, {max_f:,} cap, canvas)",
            lambda: scatter_json(exp, x, y, None, max_f, use_webgl=False),
        )
        if int(exp.z.shape[1]) >= 3:
            bench(
                "3D latent (z0,z1,z2, 120k cap)",
                lambda: scatter3d_z_json(
                    exp,
                    "z0",
                    "z1",
                    "z2",
                    None,
                    120_000,
                ),
            )
        color, diag = pairplot_color_and_diagonal(exp)
        bench(
            "Pair grid PNG (default-ish color)",
            lambda: pair_grid_png(
                exp,
                lower_color_col=color,
                diagonal_emb=diag,
                upper_style="scatter",
            ),
        )

    run_round("Cold (fresh load_experiment per round)")
    run_round("Warm (second round, new load_experiment — filesystem cache)")

    if args.profile_pairplot:
        exp = load_experiment(workdir, ep, kmeans=-1)
        color, diag = pairplot_color_and_diagonal(exp)
        pr = cProfile.Profile()
        pr.enable()
        pair_grid_png(
            exp,
            lower_color_col=color,
            diagonal_emb=diag,
            upper_style="scatter",
        )
        pr.disable()
        buf = io.StringIO()
        pstats.Stats(pr, stream=buf).sort_stats(
            pstats.SortKey.CUMULATIVE,
        ).print_stats(25)
        print("\n=== cProfile: pair_grid_png (top 25 cumulative) ===\n")
        print(buf.getvalue())

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
