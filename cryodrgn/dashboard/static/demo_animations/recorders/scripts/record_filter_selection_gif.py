#!/usr/bin/env python3
"""Headless GIF: ``cryodrgn filter``-style matplotlib UI — lasso selection + saving indices (~5 s default).

Uses the matplotlib ``Agg`` backend (no Tk window). Loads the same workspace as ``cryodrgn filter``
(``analyze.{epoch}/umap.pkl``, k-means labels, …) via :func:`cryodrgn.commands.filter.prepare_filter_workspace`.

Bootstrap (``CRYODRGN_FILTER_FORCE_TKAGG`` + ``Agg``) lives in
:mod:`cryodrgn.dashboard.static.demo_animations.recorders.scripts._record_filter_selection_gif_bootstrap`,
imported before :mod:`cryodrgn.commands.filter`.

Example::

    PYTHONPATH=/path/to/cryodrgn_beta \\
        python -m cryodrgn.dashboard.static.demo_animations.recorders.scripts.record_filter_selection_gif \\
        /scratch/.../004_train-vae_1gpu_dim.1024/ \\
        -o cryodrgn/dashboard/static/demo_animations/recorders/filter_selection_demo.gif
"""

from __future__ import annotations

import cryodrgn.dashboard.static.demo_animations.recorders.scripts._record_filter_selection_gif_bootstrap  # noqa: F401

import argparse
import io
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from cryodrgn.commands.filter import SelectFromScatter, prepare_filter_workspace
from cryodrgn.dashboard.static.demo_animations.recorders.scripts.record_dashboard_interactions_gif import (
    DEMO_ANIMATIONS_RECORDERS_DIR,
)
from cryodrgn.dashboard.static.demo_animations.recorders.scripts.record_particle_explorer_actions_gif import (
    _composite_fixed_footer_below,
    _max_footer_height_for_captions,
    normalize_frames_fixed_fps,
    save_pngs_duration_gif,
)

DEFAULT_OUTPUT = DEMO_ANIMATIONS_RECORDERS_DIR / "filter_selection_demo.gif"

# Match explorer caption band hue so plain padding above the footer isn’t visibly different.
_FOOTER_TOP_RGB = (250, 248, 244)
_SAVE_HINT_FONT_PX = 19
# Hold the final "Save Selection" footer on screen past the trimmed body clip.
_SAVE_HINT_TAIL_LINGER_SECONDS = 2.8


def _append_tail_hold_on_last_frame(
    pngs: list[bytes],
    durs_ms: list[int],
    *,
    linger_seconds: float,
    fps: float,
) -> None:
    """Append copies of ``pngs[-1]`` with durations summing to about ``linger_seconds``."""

    if linger_seconds <= 0 or len(pngs) != len(durs_ms) or not pngs:
        return
    target_ms_total = linger_seconds * 1000.0
    tick = max(10, int(round(1000.0 / max(float(fps), 1e-9))))
    tail = pngs[-1]
    added_total = 0.0
    while added_total + 1e-6 < target_ms_total:
        left_ms = target_ms_total - added_total
        chunk = tick if tick <= left_ms + 1e-6 else int(max(10, round(left_ms)))
        pngs.append(tail)
        durs_ms.append(int(chunk))
        added_total += int(chunk)


def _pad_plain_footer_band(im_rgb, footer_h: int):
    """Add ``footer_h`` px of ``_FOOTER_TOP_RGB`` below ``im_rgb`` (no inner border line)."""

    from PIL import Image

    if footer_h <= 0:
        return im_rgb
    w, h = im_rgb.size
    canvas = Image.new("RGB", (w, h + footer_h), _FOOTER_TOP_RGB)
    canvas.paste(im_rgb, (0, 0))
    return canvas


def _frame_bytes_with_optional_footer(
    fig: Figure, *, footer_h: int, caption: str | None
) -> bytes:
    """Plain pad (`caption=None`) vs explorer-style caption band (when ``caption`` is set); same height."""

    from PIL import Image

    raw = _fig_to_png(fig)
    im = Image.open(io.BytesIO(raw)).convert("RGB")
    composed = None
    try:
        if not caption:
            composed = _pad_plain_footer_band(im, footer_h)
        else:
            composed = _composite_fixed_footer_below(
                im,
                caption=caption,
                footer_h=footer_h,
                font_px=_SAVE_HINT_FONT_PX,
            )
        bio = io.BytesIO()
        composed.save(bio, format="PNG")
        return bio.getvalue()
    finally:
        im.close()
        if composed is not None and composed is not im:
            composed.close()


def _densify_closed(poly: np.ndarray, subdivisions_between: int = 10) -> np.ndarray:
    if len(poly) < 2:
        return np.asarray(poly, dtype=float)
    out: list[np.ndarray] = []
    n = len(poly)
    for i in range(n):
        p0 = poly[i]
        p1 = poly[(i + 1) % n]
        out.append(np.asarray(p0, dtype=float))
        for s in range(1, subdivisions_between + 1):
            t = s / float(subdivisions_between)
            out.append((1.0 - t) * p0 + t * p1)
    return np.stack(out, axis=0)


def _irregular_lasso_demo_polygon(
    x: np.ndarray,
    y: np.ndarray,
    *,
    sx_scale: float = 1.12,
    sy_scale: float = 1.12,
) -> np.ndarray:
    """Closed polygon in plot space — uneven angles and radii so the stroke looks hand-lassoed."""

    xc = float(np.median(x))
    yc = float(np.median(y))
    sx = max(float(np.std(x)) * sx_scale, 1e-6)
    sy = max(float(np.std(y)) * sy_scale, 1e-6)
    # Irregular angle steps (sum of gaps = one full turn).
    dtheta = np.array(
        [
            0.42,
            0.71,
            0.39,
            0.88,
            0.52,
            0.61,
            0.46,
            0.93,
            0.54,
            0.67,
            0.51,
            0.36,
            0.74,
            0.38,
            0.60,
            0.40,
        ],
        dtype=float,
    )
    dtheta *= float(2.0 * np.pi / dtheta.sum())
    theta0 = 0.19
    thetas = theta0 + np.concatenate([[0.0], np.cumsum(dtheta[:-1])])
    rad = np.array(
        [
            1.03,
            0.74,
            1.06,
            0.78,
            1.09,
            0.71,
            0.93,
            1.07,
            0.76,
            1.06,
            0.86,
            1.07,
            0.80,
            1.06,
            0.86,
            0.93,
        ],
        dtype=float,
    )
    return np.column_stack(
        [
            xc + sx * rad * np.cos(thetas),
            yc + sy * rad * np.sin(thetas),
        ]
    )


def _fig_to_png(fig: Figure) -> bytes:
    bio = io.BytesIO()
    fig.savefig(bio, format="png", facecolor=fig.get_facecolor())
    return bio.getvalue()


def render_filter_selection_gif_frames(
    plot_df,
    *,
    fps: float,
    duration_s: float,
) -> list[bytes]:
    """Build PNG sequence; ``normalize_frames_fixed_fps`` trims/pads to duration×fps."""

    base = max(60.0, float(duration_s) * float(fps))
    n_intro = max(15, int(0.15 * base))
    n_lasso = max(28, int(0.48 * base))
    n_hold = max(18, int(0.22 * base))
    n_save = max(14, int(0.22 * base))

    selector = SelectFromScatter(
        plot_df,
        pre_indices=None,
        interactive=False,
        figure_kwargs={"figsize": (11.9, 8.8), "dpi": 160},
    )
    pngs: list[bytes] = []

    sel = selector
    _save_hint = (
        "When you hit Save Selection, `cryodrgn filter` writes indices.pkl and *_inverse.pkl in the terminal.\n"
        "The first pickle is your selected particle IDs — the inverse file is the complementary set."
    )
    cw, _ch = sel.fig.canvas.get_width_height()
    footer_h = _max_footer_height_for_captions(
        [_save_hint],
        content_width_px=max(240, int(cw)),
        font_px=_SAVE_HINT_FONT_PX,
        cap_max_px=168,
    )

    def _snap_to_png_pad_or_caption(include_save_hint: bool) -> None:
        pngs.append(
            _frame_bytes_with_optional_footer(
                sel.fig,
                footer_h=footer_h,
                caption=_save_hint if include_save_hint else None,
            )
        )

    xcol, ycol = sel.xcol, sel.ycol
    xv = sel.data_table[xcol].to_numpy(dtype=float, copy=False)
    yv = sel.data_table[ycol].to_numpy(dtype=float, copy=False)

    lasso_verts = _irregular_lasso_demo_polygon(xv, yv)
    path_pts = _densify_closed(lasso_verts, subdivisions_between=12)
    lamax = max(4, len(path_pts) - 1)

    for _ in range(n_intro):
        _snap_to_png_pad_or_caption(False)

    ln = None
    for j in range(n_lasso):
        if ln is not None:
            ln.remove()
        upto = 3 + int((lamax - 3) * (j + 1) / max(1, n_lasso))
        seg = path_pts[:upto]
        (ln,) = sel.main_ax.plot(
            seg[:, 0],
            seg[:, 1],
            "-",
            color="tab:purple",
            lw=2.2,
            zorder=30,
            alpha=0.92,
        )
        sel.fig.canvas.draw()
        _snap_to_png_pad_or_caption(False)

    if ln is not None:
        ln.remove()
        sel.fig.canvas.draw()

    sel.choose_points(np.asarray(path_pts, dtype=float))
    _snap_to_png_pad_or_caption(False)

    for _ in range(n_hold):
        _snap_to_png_pad_or_caption(False)

    for _ in range(n_save):
        _snap_to_png_pad_or_caption(True)

    plt.close(sel.fig)
    return pngs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("outdir", type=Path, help="same as cryodrgn filter workdir")
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
    )
    ap.add_argument("--fps", type=float, default=25.0)
    ap.add_argument("--duration", type=float, default=5.0)
    ap.add_argument("--epoch", "-e", type=int, default=-1)
    ap.add_argument("--kmeans", "-k", type=int, default=-1)
    ap.add_argument("--max-width", type=int, default=720)
    args = ap.parse_args()

    outdir = args.outdir.expanduser().resolve()
    if not outdir.is_dir():
        print(f"error: not a directory: {outdir}", file=sys.stderr)
        return 1

    try:
        import PIL.Image  # noqa: F401
    except ImportError:
        print("error: Pillow required", file=sys.stderr)
        return 1

    try:
        plot_df, _, _ = prepare_filter_workspace(
            str(outdir),
            epoch=args.epoch,
            kmeans=args.kmeans,
            plot_inds=None,
        )
        pngs_raw = render_filter_selection_gif_frames(
            plot_df,
            fps=args.fps,
            duration_s=args.duration,
        )
        pngs_norm, durs_norm, _ = normalize_frames_fixed_fps(
            pngs_raw,
            fps=float(args.fps),
            wall_seconds=float(args.duration),
            captions=None,
        )
        _append_tail_hold_on_last_frame(
            pngs_norm,
            durs_norm,
            linger_seconds=_SAVE_HINT_TAIL_LINGER_SECONDS,
            fps=float(args.fps),
        )
        out_path = Path(args.output).expanduser().resolve()
        save_pngs_duration_gif(
            pngs_norm,
            durs_norm,
            out_path,
            max_width=max(0, int(args.max_width)),
            log=lambda _msg: None,
        )
        play_s = sum(durs_norm) / 1000.0
        nf = len(pngs_norm)
        print(
            f"Wrote {out_path} ({nf} frames, ~{float(args.fps):g} fps avg, {play_s:.2f} s)"
        )
        return 0
    except Exception as e:
        print(f"error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
