#!/usr/bin/env python3
"""Headless GIF: 3-D latent scatter — camera sweeps and Colour-by covariates.

The clip composites matplotlib preview PNGs over the Plotly viewport (WebGL blanks in headless Chromium).
Aside **covariate legend** (histogram / discrete toggles) is hidden — only ``#sc`` and camera motion are
shown so the GIF never centres cluster-toggle UI.

Depends on Pillow, Playwright/Chromium — same assumptions as ``record_dashboard_interactions_gif``.
"""

from __future__ import annotations

import argparse
import io
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from cryodrgn.dashboard.static.demo_animations.recorders.scripts.record_dashboard_interactions_gif import (
    DEMO_ANIMATIONS_RECORDERS_DIR,
    PLOTLY_CDN_ROUTE_GLOB,
    PAIRPLOT_RENDER_POLL_SNAPSHOT_STRIDE,
    PAIRPLOT_SNAP_MAX_STEPS,
    PAIRPLOT_WAIT_MS,
    REPO_ROOT,
    RENDER_POLL_FRAME_DURATION_MULT,
    _append_latent3d_composite,
    _chromium_launch_options,
    _ensure_plotly_cached,
    _free_port,
    _hide_latent3d_plotly_colorbar,
    _hide_latent3d_color_legend,
    _timed_step,
    _wait_http,
    _wait_latent3d_ready,
    sleep_snap_composite_3d,
)
from cryodrgn.dashboard.static.demo_animations.recorders.scripts.record_particle_explorer_actions_gif import (
    CaptionedExplorerBuffer,
    extend_caption_segments_min_duration,
    normalize_frames_fixed_fps,
    save_pngs_duration_gif,
    _pil_ui_font,
    _stamp_step,
)

DEFAULT_CONDA_PREFIX = Path("/scratch/gpfs/ZHONGE/mg2332/conda-envs/cdrgn_gifs")
DEFAULT_OUTPUT = DEMO_ANIMATIONS_RECORDERS_DIR / "latent3d_actions_demo.gif"
_LOG_PREFIX = "[record-latent3d-actions-gif]"
# Matches ``_DASHBOARD_CREAM`` and dashboard overlay heading colour (~``#243b53``).
_LATENT3D_PANEL_BG_RGB = (250, 248, 244)
_LATENT3D_RENDERING_LABEL_RGB = (36, 59, 83)


def _make_logger(quiet: bool):
    def log(msg: str) -> None:
        if not quiet:
            print(f"{_LOG_PREFIX} {msg}", flush=True)

    return log


def _append_latent3d_rendering_placeholder(
    buf: CaptionedExplorerBuffer,
    page,
    *,
    duration_ms: int | None = None,
) -> None:
    """Replace ``#latent3d`` in a viewport screenshot with a cream **Rendering…** panel."""

    from PIL import Image, ImageDraw

    shot = page.screenshot(type="png", full_page=False)
    canvas_im = Image.open(io.BytesIO(shot)).convert("RGBA")
    host = page.locator("#latent3d")
    if not host.count():
        buf.add_png(shot, duration_ms=duration_ms)
        return
    box = host.bounding_box()
    if box is None or box.get("width", 0) < 16 or box.get("height", 0) < 16:
        buf.add_png(shot, duration_ms=duration_ms)
        return
    bx, by = int(box["x"]), int(box["y"])
    bw, bh = int(box["width"]), int(box["height"])
    panel = Image.new("RGB", (bw, bh), _LATENT3D_PANEL_BG_RGB)
    draw = ImageDraw.Draw(panel)
    font_px = min(26, max(14, bw // 45))
    font = _pil_ui_font(font_px)
    msg = "Rendering…"
    bbox_txt = draw.textbbox((0, 0), msg, font=font)
    tw = bbox_txt[2] - bbox_txt[0]
    th = bbox_txt[3] - bbox_txt[1]
    ty = max(12, (bh - th) // 2 - bh // 22)
    tx = max(8, (bw - tw) // 2)
    draw.text(
        (tx, ty),
        msg,
        fill=_LATENT3D_RENDERING_LABEL_RGB,
        font=font,
    )
    track_margin = max(40, bw // 10)
    track_w = max(80, bw - 2 * track_margin)
    track_h = max(5, bh // 72)
    track_x = (bw - track_w) // 2
    track_y = min(bh - 52, ty + th + bh // 16)
    track_bg = (220, 214, 198)
    track_fg = (72, 108, 142)
    draw.rounded_rectangle(
        [track_x, track_y, track_x + track_w, track_y + track_h],
        radius=max(2, track_h // 2),
        fill=track_bg,
    )
    seg = max(int(track_w * 0.36), track_h + 8)
    draw.rounded_rectangle(
        [track_x, track_y, track_x + seg, track_y + track_h],
        radius=max(2, track_h // 2),
        fill=track_fg,
    )
    canvas_im.paste(panel, (bx, by))
    out = io.BytesIO()
    canvas_im.convert("RGB").save(out, format="PNG")
    buf.add_png(out.getvalue(), duration_ms=duration_ms)


def _latent3d_ready_js() -> str:
    return """() => {
      const o = document.getElementById('latent3d-rendering-overlay');
      return o && !o.classList.contains('cryo-plot-rendering-overlay--show');
    }"""


def _snap_until_latent3d_ready(
    buf: CaptionedExplorerBuffer, page, max_steps: int = PAIRPLOT_SNAP_MAX_STEPS
) -> int:
    """Poll until the 3-D overlay clears; return total ms of render-poll frames appended."""

    render_ms = max(1, int(buf.frame_ms * RENDER_POLL_FRAME_DURATION_MULT))
    rendering_ms = 0
    for step in range(max_steps):
        if step % PAIRPLOT_RENDER_POLL_SNAPSHOT_STRIDE == 0:
            _append_latent3d_rendering_placeholder(buf, page, duration_ms=render_ms)
            rendering_ms += render_ms
        if page.evaluate(_latent3d_ready_js()):
            return rendering_ms
        time.sleep(buf.frame_ms / 1000.0)
    _wait_latent3d_ready(page, PAIRPLOT_WAIT_MS)
    return rendering_ms


def _latent3d_sleep_snap_composite(
    buf: CaptionedExplorerBuffer, page, base: str, seconds: float
) -> None:
    """Screenshot chrome and paste ``/api/latent3d_preview.png`` — works when WebGL is blank."""

    n = max(1, int(round(seconds * 1000 / buf.frame_ms)))
    for _ in range(n):
        _append_latent3d_composite(buf, page, base)
        time.sleep(buf.frame_ms / 1000.0)


def _latent3d_refresh(buf: CaptionedExplorerBuffer, page, base: str) -> None:
    rms = _snap_until_latent3d_ready(buf, page)
    _wait_latent3d_ready(page, PAIRPLOT_WAIT_MS)
    _latent3d_sleep_snap_composite(buf, page, base, max(0.001, rms / 1000.0))


def _latent3d_covariate_demo_values(page, *, max_pick: int = 6) -> list[str]:
    """``#sc`` values for the GIF: skip ``none``, ``labels``, and k-means labelling."""

    raw = page.evaluate(
        """([maxPick]) => {
          const s = document.getElementById('sc');
          if (!s) return [];
          const cand = [];
          for (let i = 0; i < s.options.length; i++) {
            const o = s.options[i];
            const v = (o.value || '').trim();
            if (!v || v === 'none') continue;
            const t = ('' + (o.textContent || '') + ' ' + v).trim();
            if (/k[-\\s_]*means|\\bkmeans\\b/i.test(t)) continue;
            if (/^labels$/i.test(v)) continue;
            cand.push(v);
          }
          const seen = Object.create(null);
          const out = [];
          for (let j = 0; j < cand.length; j++) {
            const k = cand[j];
            if (seen[k]) continue;
            seen[k] = 1;
            out.push(k);
            if (out.length >= maxPick) break;
          }
          return out;
        }""",
        [int(max_pick)],
    )
    return list(raw or [])


def _latent3d_palette_click(page, name: str) -> None:
    page.evaluate(
        """(nm) => {
          const inp = document.querySelector(
            '#latent3d-palette-options input[name="latent3d_palette"][value="' + nm + '"]'
          );
          if (inp) inp.click();
        }""",
        name,
    )


def record_latent3d_sequence(
    page,
    base: str,
    buf: CaptionedExplorerBuffer,
    *,
    frame_ms: int,
    log,
) -> None:
    """Camera sweeps + ``#sc`` covariate cycling; colour-legend aside stays display:none."""

    page.goto(base + "/latent-3d", wait_until="domcontentloaded", timeout=240_000)
    time.sleep(0.35)

    _stamp_step(
        buf,
        "3-D latent scatter — drag rotates the Plotly scene; scroll zooms.",
    )
    for _ in range(max(1, int(round(0.4 * 1000 / frame_ms)))):
        _append_latent3d_rendering_placeholder(buf, page, duration_ms=frame_ms)
        time.sleep(frame_ms / 1000.0)

    _stamp_step(
        buf,
        "Rendering overlay hides once the latent scatter response is applied in the viewer.",
    )
    _latent3d_refresh(buf, page, base)
    _hide_latent3d_plotly_colorbar(page)
    _hide_latent3d_color_legend(page)

    _stamp_step(
        buf,
        "Opening orbit — tilt down while sweeping azimuth opposite the finale move.",
    )
    sleep_snap_composite_3d(
        buf,
        page,
        base,
        3.1,
        pan_start=(0.1, -0.06),
        pan_end=(-0.09, 0.08),
        zoom_start=1.06,
        zoom_end=1.0,
        elev_start=34.0,
        elev_end=14.5,
        azim_start=95.0,
        azim_end=-78.0,
    )

    _stamp_step(
        buf,
        "Second sweep — climb and wrap azimuth from the far side of the cloud.",
    )
    sleep_snap_composite_3d(
        buf,
        page,
        base,
        2.75,
        pan_start=(-0.06, 0.04),
        pan_end=(0.07, -0.03),
        zoom_start=1.0,
        zoom_end=1.04,
        elev_start=10.0,
        elev_end=38.0,
        azim_start=-118.0,
        azim_end=44.0,
    )

    demo_cols = _latent3d_covariate_demo_values(page, max_pick=5)
    _stamp_step(
        buf,
        "Color by — pick different table columns in the controls; no covariate legend strip here.",
    )
    if demo_cols:
        try:
            for i, col in enumerate(demo_cols):
                page.select_option("#sc", col)
                _latent3d_refresh(buf, page, base)
                _hide_latent3d_plotly_colorbar(page)
                _hide_latent3d_color_legend(page)
                _latent3d_sleep_snap_composite(buf, page, base, 0.42)
                if i == 1:
                    sleep_snap_composite_3d(
                        buf,
                        page,
                        base,
                        2.35,
                        elev_start=28.0,
                        elev_end=17.0,
                        azim_start=48.0,
                        azim_end=-98.0,
                        zoom_start=1.02,
                        zoom_end=0.98,
                    )
                if i == 2:
                    _latent3d_palette_click(page, "Plasma")
                    _latent3d_refresh(buf, page, base)
                    _hide_latent3d_plotly_colorbar(page)
                    _hide_latent3d_color_legend(page)
                    _latent3d_sleep_snap_composite(buf, page, base, 0.4)
        except Exception as exc:
            log(f"Color by sequence: {exc}")
    else:
        _stamp_step(
            buf,
            "No suitable colour columns in the menu — stay on extra camera motion.",
        )
        sleep_snap_composite_3d(
            buf,
            page,
            base,
            2.85,
            pan_start=(0.0, 0.0),
            pan_end=(0.05, 0.06),
            zoom_start=1.04,
            zoom_end=1.0,
            elev_start=22.0,
            elev_end=30.0,
            azim_start=-40.0,
            azim_end=130.0,
        )

    _stamp_step(
        buf,
        "Third pan — zoom in slightly while pulling the camera under the cloud.",
    )
    sleep_snap_composite_3d(
        buf,
        page,
        base,
        2.45,
        pan_start=(0.04, -0.02),
        pan_end=(-0.05, 0.05),
        zoom_start=0.98,
        zoom_end=1.09,
        elev_start=18.0,
        elev_end=32.0,
        azim_start=75.0,
        azim_end=-105.0,
    )

    _stamp_step(
        buf,
        "Closing orbit pair — eased azimuth and elevation while chrome stays fixed.",
    )
    sleep_snap_composite_3d(
        buf,
        page,
        base,
        3.4,
        elev_start=20.0,
        elev_end=34.0,
        azim_start=-72.0,
        azim_end=58.0,
    )
    sleep_snap_composite_3d(
        buf,
        page,
        base,
        2.6,
        zoom_start=1.0,
        zoom_end=1.08,
        elev_start=34.0,
        elev_end=24.0,
        azim_start=58.0,
        azim_end=-118.0,
    )


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
    ap.add_argument(
        "--kmeans",
        "-k",
        type=int,
        default=-1,
        help="forwarded to ``cryodrgn dashboard --kmeans`` when ≥ 0 (default -1: workspace default)",
    )
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
    if int(args.kmeans) < -1:
        print(
            "error: --kmeans / -k must be ≥ -1 (-1 skips passing --kmeans to the dashboard)",
            file=sys.stderr,
        )
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
    if int(args.kmeans) >= 0:
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
            with _timed_step(log, "Record 3-D visualizer clip", slow_note=True):
                record_latent3d_sequence(page, base, buf, frame_ms=frame_ms, log=log)
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
