#!/usr/bin/env python3
"""Headless GIF: pairwise latent grid (k-means + PCA diagonal + hex upstairs + discrete toggles).

Opens ``/pairplot``: switches to **k-means** colouring first, sets the diagonal to **PCA**, **hex**
densities above diagonal, then **toggles cluster visibility** via the discrete legend.

Same encoding style as ``record_particle_explorer_actions_gif`` — optional ``--annotate-caption`` band,
``--annotate-min-segment-seconds``, and ``--high-resolution-output``. Defaults **~22 s × 40 fps** so
captions breathe before caption-extension.

Depends on Pillow, Playwright/Chromium — same assumptions as ``record_dashboard_interactions_gif``.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from cryodrgn.dashboard.record_dashboard_interactions_gif import (
    DEMO_ANIMATIONS_DIR,
    PLOTLY_CDN_ROUTE_GLOB,
    PAIRPLOT_SNAP_MAX_STEPS,
    PAIRPLOT_WAIT_MS,
    REPO_ROOT,
    _chromium_launch_options,
    _ensure_plotly_cached,
    _free_port,
    _snap_until_pairplot_ready,
    _timed_step,
    _wait_http,
    _wait_pairplot_ready,
)
from cryodrgn.dashboard.record_particle_explorer_actions_gif import (
    CaptionedExplorerBuffer,
    extend_caption_segments_min_duration,
    normalize_frames_fixed_fps,
    save_pngs_duration_gif,
    _stamp_step,
)

DEFAULT_CONDA_PREFIX = Path("/projects/CRYOEM/zhonglab/mg2332/conda_envs/cdrgn_beta")
DEFAULT_OUTPUT = DEMO_ANIMATIONS_DIR / "pair_grid_actions_demo.gif"
_LOG_PREFIX = "[record-pairplot-actions-gif]"


def _make_logger(quiet: bool):
    def log(msg: str) -> None:
        if not quiet:
            print(f"{_LOG_PREFIX} {msg}", flush=True)

    return log


def _pair_click_named_radio(page, name: str, value: str) -> None:
    """Click ``input[name][value]`` in pair controls via JS — safe for awkward ``value`` strings."""
    page.evaluate(
        """([n, v]) => {
          const inputs = document.querySelectorAll(
            '#pair-controls-card input[type="radio"][name="' + n + '"]');
          for (let i = 0; i < inputs.length; i++) {
            if ((inputs[i].value || '') === v) {
              inputs[i].click();
              return;
            }
          }
        }""",
        [name, value],
    )


def _pair_refresh(buf: CaptionedExplorerBuffer, page) -> None:
    rms = _snap_until_pairplot_ready(buf, page, max_steps=PAIRPLOT_SNAP_MAX_STEPS)
    _wait_pairplot_ready(page, PAIRPLOT_WAIT_MS)
    buf.sleep_snap(page, max(0.001, rms / 1000.0))


def _pair_has_umap(page) -> bool:
    return bool(
        page.evaluate(
            """() => {
          const u = document.querySelector('input[name="diag_emb"][value="umap"]');
          return !!(u && !u.disabled);
        }"""
        )
    )


def _pair_kmeans_radio_value(page) -> str | None:
    """Return ``value`` of the colour covariate radio matching k-means in its label."""
    return page.evaluate(
        """() => {
          const inputs = document.querySelectorAll(
              '#pair-controls-card input[type="radio"][name="color_cov"]');
          for (let i = 0; i < inputs.length; i++) {
            const inp = inputs[i];
            const lab = inp.closest('label');
            const t = ('' + (lab && lab.textContent || '') + ' ' + (inp.value || '')).trim();
            if (/k[-\\s_]*means|kmeans/i.test(t)) return inp.value || null;
          }
          return null;
        }"""
    )


def record_pairplot_sequence(
    page,
    base: str,
    buf: CaptionedExplorerBuffer,
    *,
    frame_ms: int,
    log,
) -> None:
    """K-means colour → PCA diagonal → hex upper triangle → discrete cluster toggles."""

    page.goto(base + "/pairplot", wait_until="domcontentloaded", timeout=240_000)
    time.sleep(0.35)

    _stamp_step(
        buf,
        "Pair-grid generator — latent z_i vs z_j panels; diagonal shows a 2D embedding.",
    )
    for _ in range(max(1, int(round(0.4 * 1000 / frame_ms)))):
        buf.snap(page)
        time.sleep(frame_ms / 1000.0)

    _stamp_step(
        buf,
        "Initial PNG grid — skeleton cells resolve when the matplotlib render finishes.",
    )
    _pair_refresh(buf, page)

    km_val = _pair_kmeans_radio_value(page)

    has_umap = _pair_has_umap(page)

    if km_val:
        _stamp_step(
            buf,
            "Discrete k-means colour — below-diagonal panels use clustering labels.",
        )
        try:
            _pair_click_named_radio(page, "color_cov", km_val)
            _pair_refresh(buf, page)
            page.wait_for_function(
                """() => {
                  const w = document.getElementById('pair-color-discrete-wrap');
                  const sw = document.getElementById('pair-color-discrete-switches');
                  return !!(w && w.classList.contains('cryo-cc-discrete-wrap--show')
                    && sw && sw.querySelector('label.cryo-cc-discrete-switch'));
                }""",
                timeout=120_000,
            )
        except Exception as exc:
            log(f"k-means / discrete legend readiness: {exc}")
    else:
        _stamp_step(
            buf,
            "No k-means colour column detected — PCA / hex segments still shown; skips cluster toggles.",
        )
        for _ in range(max(1, int(round(0.45 * 1000 / frame_ms)))):
            buf.snap(page)
            time.sleep(frame_ms / 1000.0)

    _stamp_step(
        buf,
        "Diagonal PCA — PCA 2-D embedding replaces UMAP tiles when analyse includes UMAP.",
    )
    if has_umap:
        _pair_click_named_radio(page, "diag_emb", "pc")
        _pair_refresh(buf, page)
    else:
        for _ in range(max(1, int(round(0.45 * 1000 / frame_ms)))):
            buf.snap(page)
            time.sleep(frame_ms / 1000.0)

    _stamp_step(
        buf,
        "Above diagonal — hex densities emphasise latent pair occupancy instead of scatter points.",
    )
    _pair_click_named_radio(page, "upper_style", "hex")
    _pair_refresh(buf, page)

    if km_val:
        sw = page.locator("#pair-color-discrete-switches label.cryo-cc-discrete-switch")
        if sw.count():
            _stamp_step(
                buf,
                "Cluster toggles alone — mute or revive individual classes across the grid.",
            )
            try:
                sw.first.scroll_into_view_if_needed()
                sw.first.click(force=True)
                _pair_refresh(buf, page)
                if sw.count() >= 2:
                    sw.nth(1).scroll_into_view_if_needed()
                    sw.nth(1).click(force=True)
                    _pair_refresh(buf, page)
            except Exception as e:
                log(f"Discrete toggle interaction: {e}")
            for _ in range(max(1, int(round(0.75 * 1000 / frame_ms)))):
                buf.snap(page)
                time.sleep(frame_ms / 1000.0)
        else:
            log("Discrete switches missing after selecting k-means; skip toggling.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "outdir", type=Path, help="experiment output dir for cryodrgn dashboard"
    )
    ap.add_argument("-o", "--output", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--port", type=int, default=0)
    ap.add_argument("--conda-prefix", type=Path, default=DEFAULT_CONDA_PREFIX)
    ap.add_argument("--fps", type=float, default=40.0)
    ap.add_argument(
        "--duration",
        type=float,
        default=22.0,
        help="GIF playback duration (nominal wall before annotate extension; default 22)",
    )
    ap.add_argument("--cpus", type=int, default=8)
    ap.add_argument("--epoch", "-e", type=int, default=-1)
    ap.add_argument("--kmeans", "-k", type=int, default=-1)
    ap.add_argument("--max-width", type=int, default=512)
    ap.add_argument("--annotate-caption", action="store_true")
    ap.add_argument("--high-resolution-output", type=Path, default=None, metavar="PATH")
    ap.add_argument("--high-resolution-max-width", type=int, default=0, metavar="PX")
    ap.add_argument("--caption-font-size", type=int, default=21, metavar="PX")
    ap.add_argument(
        "--annotate-min-segment-seconds",
        type=float,
        default=3.5,
        metavar="SEC",
    )
    ap.add_argument("--headed", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()

    fps = float(args.fps)
    wall_s = float(args.duration)
    if fps <= 0 or wall_s <= 0:
        print("error: fps and duration must be positive", file=sys.stderr)
        return 1
    frame_ms = max(1, int(round(1000.0 / fps)))
    log = _make_logger(args.quiet)
    outdir = args.outdir.expanduser().resolve()
    if not outdir.is_dir():
        print(f"error: not a directory: {outdir}", file=sys.stderr)
        return 1
    conda_exe = shutil.which("conda")
    if conda_exe is None:
        print("error: conda not found on PATH", file=sys.stderr)
        return 1
    try:
        import PIL.Image  # noqa: F401
        from playwright.sync_api import sync_playwright
    except ImportError:
        print("error: pillow and playwright required", file=sys.stderr)
        return 1

    port = args.port or _free_port()
    base = f"http://127.0.0.1:{port}"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO_ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    dash = [
        conda_exe,
        "run",
        "-p",
        str(Path(args.conda_prefix).expanduser().resolve()),
        "--no-capture-output",
        "cryodrgn",
        "dashboard",
        str(outdir),
        "--no-browser",
        "--host",
        "127.0.0.1",
        "--port",
        str(port),
        "--cpus",
        str(int(args.cpus)),
    ]
    if args.epoch >= 0:
        dash.extend(["--epoch", str(int(args.epoch))])
    if args.kmeans >= 0:
        dash.extend(["--kmeans", str(int(args.kmeans))])

    proc = subprocess.Popen(
        dash,
        env=env,
        cwd=str(REPO_ROOT),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        text=True,
    )

    buf = CaptionedExplorerBuffer(frame_ms)

    try:
        log(f"Dashboard subprocess conda -p {args.conda_prefix} …")
        _wait_http(base + "/", timeout_s=420.0, log=log)
        plotly_cached = _ensure_plotly_cached(log)

        ph, launch_args = _chromium_launch_options(headed=args.headed)
        with sync_playwright() as p:
            browser = p.chromium.launch(headless=ph, args=launch_args)
            context = browser.new_context(
                viewport={"width": 1440, "height": 900},
                device_scale_factor=1,
            )
            if plotly_cached is not None:
                pb = plotly_cached.read_bytes()

                def _fulfill(route):
                    route.fulfill(
                        status=200,
                        headers={
                            "content-type": "application/javascript; charset=utf-8",
                            "cache-control": "public, max-age=86400",
                        },
                        body=pb,
                    )

                context.route(PLOTLY_CDN_ROUTE_GLOB, _fulfill)

            page = context.new_page()
            with _timed_step(log, "Record pair-grid clip", slow_note=True):
                record_pairplot_sequence(page, base, buf, frame_ms=frame_ms, log=log)
            browser.close()

        pngs_norm, durs_norm, caps_norm = normalize_frames_fixed_fps(
            buf.pngs,
            fps=fps,
            wall_seconds=wall_s,
            captions=list(buf.captions_per_frame) if args.annotate_caption else None,
        )
        annotate = bool(args.annotate_caption)
        cap_lines = caps_norm if annotate else None
        caption_font_px = max(8, int(args.caption_font_size))
        min_seg = float(args.annotate_min_segment_seconds)
        if annotate and cap_lines is not None and min_seg > 0:
            (
                pngs_norm,
                durs_norm,
                cap_lines,
                n_seg_ext,
            ) = extend_caption_segments_min_duration(
                pngs_norm,
                durs_norm,
                cap_lines,
                min_segment_s=min_seg,
                fps=fps,
            )
            if n_seg_ext:
                new_play_s = sum(durs_norm) / 1000.0
                extra_s = new_play_s - float(wall_s)
                tail = (
                    f", +{extra_s:.1f}s vs --duration target)"
                    if extra_s > 0.05
                    else ")"
                )
                log(
                    f"Annotated captions: ensured ≥ {min_seg:g}s per caption span "
                    f"({n_seg_ext} span(s) extended; playback now ~{new_play_s:.1f}s{tail}"
                )

        std_path = Path(args.output).resolve()
        save_pngs_duration_gif(
            pngs_norm,
            durs_norm,
            std_path,
            max_width=max(0, int(args.max_width)),
            log=log,
            annotate=annotate,
            step_captions=cap_lines,
            caption_font_px=caption_font_px,
        )
        nf = len(pngs_norm)
        play_s = sum(durs_norm) / 1000.0
        print(
            f"Wrote {std_path} ({nf} frames, variable delay ~{fps:g} fps avg, "
            f"{play_s:.2f} s playback)"
        )
        if args.high_resolution_output is not None:
            hr_path = Path(args.high_resolution_output).expanduser().resolve()
            hr_mw = max(0, int(args.high_resolution_max_width))
            save_pngs_duration_gif(
                pngs_norm,
                durs_norm,
                hr_path,
                max_width=hr_mw,
                log=log,
                annotate=annotate,
                step_captions=cap_lines,
                caption_font_px=caption_font_px,
            )
            tag = "full native width" if hr_mw == 0 else f"max width {hr_mw}px"
            print(f"Wrote {hr_path} (high-resolution copy, {tag})")
        return 0
    except Exception as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=20)
        except subprocess.TimeoutExpired:
            proc.kill()


if __name__ == "__main__":
    raise SystemExit(main())
