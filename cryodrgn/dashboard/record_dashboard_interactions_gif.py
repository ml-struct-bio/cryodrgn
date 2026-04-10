#!/usr/bin/env python3
"""Headless capture: particle explorer (lasso + image cache), 3D latent view, pair-plot covariates.

Uses a uniform ~33 ms frame interval (30 fps). Total GIF duration is capped at **16 s**;
long segments (image preload, pair-plot re-renders) are trimmed with frame caps and shorter holds.

**3D plot:** WebGL screenshots are often blank in automated Chromium. Frames for the latent-3D page
composite a server-rendered matplotlib PNG (``GET /api/latent3d_preview.png``) over the ``#latent3d``
panel, so the GIF shows the same data as the live Plotly view without needing ``--headed`` / xvfb.

Depends on the same setup as ``record_dashboard_demo_gif.py`` (playwright, pillow, chromium).

Optional ``--pairplot-companion`` also writes ``dashboard_interactions_pairplot_after_selection.gif``
(long holds after covariate changes) as a separate file.

Example::

    PYTHONPATH=/path/to/cryodrgn_beta conda run -n cdrgn-4.2.0_py-3.13_beta \\
        python scripts/record_dashboard_interactions_gif.py \\
        /scratch/.../004_train-vae_1gpu_dim.1024/

By default the script logs progress and elapsed time for slow steps. Pass ``-q`` / ``--quiet`` for
minimal output (final ``Wrote …`` lines only, plus errors on stderr).
"""

from __future__ import annotations

import argparse
import contextlib
import io
import math
import os
import shutil
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
FRAME_MS_30FPS = 33
MAX_GIF_WALL_MS = 16_000
# Companion GIF can run longer — mostly static holds on the finished pairplot PNG.
COMPANION_MAX_WALL_MS = 20_000
DEFAULT_PAIRPLOT_COMPANION = (
    REPO_ROOT
    / "cryodrgn/dashboard/static/dashboard_interactions_pairplot_after_selection.gif"
)

_LOG_PREFIX = "[record-dashboard-interactions-gif]"


def _make_logger(quiet: bool):
    """Return a log function; no-op when ``quiet`` is True."""

    def log(msg: str) -> None:
        if not quiet:
            print(f"{_LOG_PREFIX} {msg}", flush=True)

    return log


@contextlib.contextmanager
def _timed_step(log, intro: str, *, slow_note: bool = True):
    """Log start, run block, then log elapsed time. ``intro`` is a short phrase (no trailing ellipsis)."""
    suffix = " — this can take a while" if slow_note else ""
    log(f"{intro} …{suffix}")
    t0 = time.monotonic()
    try:
        yield
    finally:
        dt = time.monotonic() - t0
        log(f"{intro} done in {dt:.2f}s")


def _chromium_args_headless_webgl() -> list[str]:
    """Best-effort flags for WebGL in **headless** Chromium (often still blank for Plotly 3D)."""
    return [
        "--enable-unsafe-swiftshader",
        "--use-gl=angle",
        "--use-angle=swiftshader-webgl",
        "--ignore-gpu-blocklist",
        "--headless=new",
        "--no-sandbox",
        "--disable-setuid-sandbox",
        "--disable-dev-shm-usage",
    ]


def _chromium_launch_options(*, headed: bool) -> tuple[bool, list[str]]:
    """(headless, args). Use ``headed=True`` with ``xvfb-run`` or a real DISPLAY for working 3D."""
    common = [
        "--no-sandbox",
        "--disable-setuid-sandbox",
        "--disable-dev-shm-usage",
    ]
    if headed:
        return False, common
    return True, _chromium_args_headless_webgl()


def _free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return int(s.getsockname()[1])


def _wait_http(url: str, timeout_s: float = 180.0, *, log=None) -> None:
    if log:
        with _timed_step(log, f"Waiting for dashboard at {url}", slow_note=True):
            _wait_http(url, timeout_s=timeout_s, log=None)
        return
    deadline = time.monotonic() + timeout_s
    last_err: Exception | None = None
    while time.monotonic() < deadline:
        try:
            with urllib.request.urlopen(url, timeout=5) as r:
                if r.status == 200:
                    return
        except (urllib.error.URLError, OSError) as e:
            last_err = e
        time.sleep(0.4)
    raise RuntimeError(f"Server did not respond at {url!r}: {last_err}")


class FrameBuffer:
    def __init__(self, frame_ms: int) -> None:
        self.frame_ms = frame_ms
        self.pngs: list[bytes] = []

    def snap(self, page) -> None:
        self.pngs.append(page.screenshot(type="png", full_page=False))

    def sleep_snap(self, page, seconds: float) -> None:
        n = max(1, int(round(seconds * 1000 / self.frame_ms)))
        for _ in range(n):
            self.snap(page)
            time.sleep(self.frame_ms / 1000.0)


def _scatter_plot_box(page):
    host = page.locator("#scatter")
    if not host.count():
        return None
    return host.bounding_box()


def _wait_scatter_ready(page, timeout_ms: float) -> None:
    page.wait_for_function(
        """() => {
          const el = document.getElementById('scatter-rendering-overlay');
          return el && !el.classList.contains('cryo-plot-rendering-overlay--show');
        }""",
        timeout=timeout_ms,
    )


def _wait_latent3d_ready(page, timeout_ms: float) -> None:
    page.wait_for_function(
        """() => {
          const el = document.getElementById('latent3d-rendering-overlay');
          return el && !el.classList.contains('cryo-plot-rendering-overlay--show');
        }""",
        timeout=timeout_ms,
    )


def _pairplot_ready_js() -> str:
    return """() => {
      const o = document.getElementById('pairplot-rendering-overlay');
      const img = document.getElementById('pairplot');
      const ready = o && !o.classList.contains('cryo-plot-rendering-overlay--show');
      const drawn = img && img.naturalWidth > 0;
      return !!(ready && drawn);
    }"""


def _wait_pairplot_ready(page, timeout_ms: float) -> None:
    page.wait_for_function(_pairplot_ready_js(), timeout=timeout_ms)


def _snap_until_pairplot_ready(buf: FrameBuffer, page, max_steps: int = 120) -> None:
    for _ in range(max_steps):
        buf.snap(page)
        if page.evaluate(_pairplot_ready_js()):
            return
        time.sleep(buf.frame_ms / 1000.0)
    _wait_pairplot_ready(page, 120_000)


def _latent3d_query_from_dom(page) -> dict[str, str] | None:
    return page.evaluate(
        """() => {
      const sx = document.getElementById('sx');
      const sy = document.getElementById('sy');
      const sz = document.getElementById('sz');
      const sc = document.getElementById('sc');
      const p = document.querySelector('input[name="latent3d_palette"]:checked');
      if (!sx || !sy || !sz || !sc) return null;
      return {
        x: sx.value,
        y: sy.value,
        z: sz.value,
        color: sc.value,
        palette: (p && p.value) ? p.value : 'Viridis',
      };
    }"""
    )


def _append_latent3d_composite(
    buf: FrameBuffer,
    page,
    base: str,
    *,
    zoom: float = 1.0,
    pan_x: float = 0.0,
    pan_y: float = 0.0,
    elev: float = 22.0,
    azim: float = -65.0,
    plot_only: bool = False,
) -> None:
    from PIL import Image

    params = _latent3d_query_from_dom(page)
    if not params:
        buf.snap(page)
        return
    req_params = dict(params)
    req_params["elev"] = f"{float(elev):.4f}"
    req_params["azim"] = f"{float(azim):.4f}"
    q = urllib.parse.urlencode(req_params)
    try:
        with urllib.request.urlopen(
            f"{base}/api/latent3d_preview.png?{q}", timeout=120
        ) as r:
            prev_bytes = r.read()
    except (urllib.error.URLError, OSError):
        buf.snap(page)
        return

    preview = Image.open(io.BytesIO(prev_bytes)).convert("RGBA")
    if plot_only:
        vp = page.viewport_size or {"width": 1440, "height": 900}
        bw, bh = int(vp["width"]), int(vp["height"])
        bx, by = 0, 0
        canvas = Image.new("RGBA", (bw, bh), (250, 248, 244, 255))
    else:
        shot = page.screenshot(type="png", full_page=False)
        canvas = Image.open(io.BytesIO(shot)).convert("RGBA")
        host = page.locator("#latent3d")
        if not host.count():
            buf.pngs.append(shot)
            return
        box = host.bounding_box()
        if not box:
            buf.pngs.append(shot)
            return
        bx, by = int(box["x"]), int(box["y"])
        bw, bh = int(box["width"]), int(box["height"])
    pw, ph = preview.size
    if bw < 8 or bh < 8:
        buf.pngs.append(page.screenshot(type="png", full_page=False))
        return
    panel = Image.new("RGBA", (bw, bh), (250, 248, 244, 255))
    scale = min(bw / pw, bh / ph) * 0.97 * max(0.5, zoom)
    nw, nh = max(1, int(pw * scale)), max(1, int(ph * scale))
    pr = preview.resize((nw, nh), Image.Resampling.LANCZOS)
    ox = (bw - nw) // 2 + int(pan_x * bw * 0.14)
    oy = (bh - nh) // 2 + int(pan_y * bh * 0.14)
    panel.paste(pr, (ox, oy), pr)
    canvas.paste(panel, (bx, by), panel)
    out = io.BytesIO()
    canvas.convert("RGB").save(out, format="PNG")
    buf.pngs.append(out.getvalue())


def sleep_snap_composite_3d(
    buf: FrameBuffer,
    page,
    base: str,
    seconds: float,
    *,
    zoom_start: float = 1.0,
    zoom_end: float = 1.0,
    pan_start: tuple[float, float] = (0.0, 0.0),
    pan_end: tuple[float, float] = (0.0, 0.0),
    elev_start: float = 22.0,
    elev_end: float = 22.0,
    azim_start: float = -65.0,
    azim_end: float = -65.0,
    plot_only: bool = False,
) -> None:
    n = max(1, int(round(seconds * 1000 / buf.frame_ms)))
    for i in range(n):
        t = i / max(1, n - 1)
        z = zoom_start + (zoom_end - zoom_start) * t
        px = pan_start[0] + (pan_end[0] - pan_start[0]) * t
        py = pan_start[1] + (pan_end[1] - pan_start[1]) * t
        # Smooth ease for camera sweep through both y/elevation and z/azimuth angles.
        tw = 0.5 - 0.5 * math.cos(math.pi * t)
        elev = elev_start + (elev_end - elev_start) * tw
        azim = azim_start + (azim_end - azim_start) * tw
        _append_latent3d_composite(
            buf,
            page,
            base,
            zoom=z,
            pan_x=px,
            pan_y=py,
            elev=elev,
            azim=azim,
            plot_only=plot_only,
        )
        time.sleep(buf.frame_ms / 1000.0)


def _montage_cache_ready_js() -> str:
    return """() => {
      const o = document.getElementById('montage-preload-overlay');
      const busy = o && o.classList.contains('cryo-plot-rendering-overlay--show');
      if (busy) return false;
      const imgs = document.querySelectorAll('#montage-grid img');
      if (!imgs.length) return false;
      for (let i = 0; i < imgs.length; i++) {
        const im = imgs[i];
        if (im && im.src && im.src.indexOf('data:image') === 0 && im.naturalWidth > 0) return true;
      }
      return false;
    }"""


def record_sequence(page, base: str, buf: FrameBuffer, *, log) -> None:
    # --- Particle explorer ---
    with _timed_step(log, "Particle explorer: load page and wait for scatter plot"):
        page.goto(base + "/explorer", wait_until="domcontentloaded", timeout=120_000)
        _wait_scatter_ready(page, 180_000)
    buf.sleep_snap(page, 0.55)

    box = _scatter_plot_box(page)
    if not box:
        buf.snap(page)
    else:
        w, h = box["width"], box["height"]
        page.evaluate(
            """() => {
          const gd = document.getElementById('scatter');
          if (!gd || !window.Plotly) return false;
          Plotly.relayout(gd, { dragmode: 'lasso' });
          return true;
        }"""
        )
        # Freehand lasso path (Plotly closes the region on mouseup).
        rel = [
            (0.32, 0.48),
            (0.44, 0.38),
            (0.58, 0.42),
            (0.64, 0.55),
            (0.52, 0.64),
            (0.38, 0.58),
            (0.30, 0.50),
        ]
        page.mouse.move(box["x"] + w * rel[0][0], box["y"] + h * rel[0][1])
        page.mouse.down()
        for i in range(1, len(rel)):
            page.mouse.move(box["x"] + w * rel[i][0], box["y"] + h * rel[i][1])
            time.sleep(buf.frame_ms / 1000.0)
            buf.snap(page)
        page.mouse.up()
        buf.sleep_snap(page, 0.55)

    with _timed_step(
        log,
        "Particle explorer: preload image montage after View images (server decode + grid)",
        slow_note=True,
    ):
        page.locator("#btn-view-images").click()
        max_preload_frames = 40
        for _ in range(max_preload_frames):
            if page.evaluate(_montage_cache_ready_js()):
                break
            buf.snap(page)
            time.sleep(buf.frame_ms / 1000.0)
        try:
            page.wait_for_function(_montage_cache_ready_js(), timeout=120_000)
        except Exception:
            pass

    with _timed_step(
        log,
        "Particle explorer: recolour scatter + palette (may re-fetch plot)",
        slow_note=True,
    ):
        sc = page.locator("#sc")
        if sc.count():
            opts = page.locator("#sc option")
            n = opts.count()
            if n > 2:
                val = opts.nth(2).get_attribute("value")
                if val:
                    sc.select_option(value=val)
                    _wait_scatter_ready(page, 180_000)
                    turbo = page.locator('input[name="scatter_palette"][value="Turbo"]')
                    if turbo.count():
                        turbo.scroll_into_view_if_needed()
                        turbo.click(force=True)
                        _wait_scatter_ready(page, 180_000)
    buf.sleep_snap(page, 1.35)

    # --- 3D visualizer (matplotlib composite over #latent3d; WebGL capture is unreliable) ---
    with _timed_step(
        log, "3D latent page: load and wait for Plotly overlay", slow_note=True
    ):
        page.goto(base + "/latent-3d", wait_until="domcontentloaded", timeout=120_000)
        _wait_latent3d_ready(page, 180_000)
    # Keep this section concise since we now only retain the latter half of the 3D sweep.
    sleep_snap_composite_3d(
        buf,
        page,
        base,
        0.35,
        elev_start=18.0,
        elev_end=22.0,
        azim_start=-80.0,
        azim_end=-65.0,
        plot_only=True,
    )

    with _timed_step(
        log, "3D: apply colour covariate and wait for redraw", slow_note=True
    ):
        sc = page.locator("#sc")
        if sc.count():
            opts = page.locator("#sc option")
            n = opts.count()
            if n > 2:
                val = opts.nth(2).get_attribute("value")
                if val:
                    sc.select_option(value=val)
        _wait_latent3d_ready(page, 180_000)
    sleep_snap_composite_3d(
        buf,
        page,
        base,
        0.45,
        elev_start=22.0,
        elev_end=30.0,
        azim_start=-65.0,
        azim_end=-38.0,
        plot_only=True,
    )

    with _timed_step(log, "3D: Plasma palette and wait for redraw", slow_note=True):
        plasma = page.locator('input[name="latent3d_palette"][value="Plasma"]')
        if plasma.count():
            plasma.scroll_into_view_if_needed()
            plasma.click(force=True)
        _wait_latent3d_ready(page, 180_000)
    sleep_snap_composite_3d(
        buf,
        page,
        base,
        0.35,
        elev_start=30.0,
        elev_end=28.0,
        azim_start=-38.0,
        azim_end=-22.0,
        plot_only=True,
    )

    with _timed_step(
        log,
        "3D: capture plot-only frames (many server PNG composites — slow)",
        slow_note=True,
    ):
        sleep_snap_composite_3d(
            buf,
            page,
            base,
            2.4,
            zoom_start=1.2,
            zoom_end=1.42,
            pan_start=(0.28, 0.16),
            pan_end=(-0.12, 0.26),
            elev_start=14.0,
            elev_end=34.0,
            azim_start=42.0,
            azim_end=88.0,
            plot_only=True,
        )
        sleep_snap_composite_3d(
            buf,
            page,
            base,
            1.0,
            zoom_start=1.42,
            zoom_end=1.5,
            elev_start=34.0,
            elev_end=30.0,
            azim_start=88.0,
            azim_end=96.0,
            plot_only=True,
        )

    with _timed_step(
        log, "Pair plot: load page and wait for initial PNG", slow_note=True
    ):
        page.goto(base + "/pairplot", wait_until="domcontentloaded", timeout=120_000)
        _wait_pairplot_ready(page, 180_000)
    buf.sleep_snap(page, 0.45)

    radios = page.locator('input[name="color_cov"]')
    rc = radios.count()
    if rc >= 2:
        with _timed_step(
            log,
            "Pair plot: first covariate change and wait for PNG",
            slow_note=True,
        ):
            r = radios.nth(min(1, rc - 1))
            r.scroll_into_view_if_needed()
            r.click(force=True)
            _snap_until_pairplot_ready(buf, page, max_steps=22)
            _wait_pairplot_ready(page, 120_000)
        buf.sleep_snap(page, 1.05)

    if rc >= 3:
        with _timed_step(
            log,
            "Pair plot: second covariate change and wait for PNG",
            slow_note=True,
        ):
            r2 = radios.nth(min(2, rc - 1))
            r2.scroll_into_view_if_needed()
            r2.click(force=True)
            _snap_until_pairplot_ready(buf, page, max_steps=22)
            _wait_pairplot_ready(page, 120_000)
        buf.sleep_snap(page, 1.25)


def record_pairplot_after_selection_sequence(
    page, base: str, buf: FrameBuffer, *, log
) -> None:
    """Pair-plot only: covariate radio ``selection`` → wait for PNG → long hold on finished grid."""
    with _timed_step(log, "Companion: pair plot load + initial PNG", slow_note=True):
        page.goto(base + "/pairplot", wait_until="domcontentloaded", timeout=120_000)
        _wait_pairplot_ready(page, 180_000)
    buf.sleep_snap(page, 0.55)

    radios = page.locator('input[name="color_cov"]')
    rc = radios.count()
    if rc >= 2:
        with _timed_step(
            log,
            "Companion: pair plot covariate 1 + long hold segment",
            slow_note=True,
        ):
            r = radios.nth(1)
            r.scroll_into_view_if_needed()
            r.click(force=True)
            _snap_until_pairplot_ready(buf, page, max_steps=120)
        buf.sleep_snap(page, 2.35)
    if rc >= 3:
        with _timed_step(
            log,
            "Companion: pair plot covariate 2 + long hold segment",
            slow_note=True,
        ):
            r2 = radios.nth(2)
            r2.scroll_into_view_if_needed()
            r2.click(force=True)
            _snap_until_pairplot_ready(buf, page, max_steps=120)
        buf.sleep_snap(page, 3.85)
    elif rc < 2:
        buf.snap(page)


def _save_buf_to_gif(
    buf: FrameBuffer,
    output: Path,
    *,
    frame_ms: int,
    max_width: int,
    max_wall_ms: int,
    log,
    label: str,
) -> tuple[int, float]:
    from PIL import Image

    if not buf.pngs:
        raise RuntimeError("no frames to save")
    with _timed_step(
        log,
        f"Encode GIF ({label}): resize {len(buf.pngs)} frames, write {output.name}",
        slow_note=True,
    ):
        max_frames = max(1, max_wall_ms // max(1, frame_ms))
        pngs = buf.pngs[:max_frames]
        pil_frames: list[Image.Image] = []
        for raw in pngs:
            im = Image.open(io.BytesIO(raw)).convert("RGB")
            if max_width and im.width > max_width:
                ratio = max_width / im.width
                pil_frames.append(
                    im.resize(
                        (max_width, max(1, int(im.height * ratio))),
                        Image.Resampling.LANCZOS,
                    )
                )
            else:
                pil_frames.append(im)
        durations = [frame_ms] * len(pil_frames)
        output.parent.mkdir(parents=True, exist_ok=True)
        pil_frames[0].save(
            output,
            save_all=True,
            append_images=pil_frames[1:],
            duration=durations,
            loop=0,
            optimize=True,
        )
    wall_s = len(pil_frames) * frame_ms / 1000.0
    return len(pil_frames), wall_s


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", type=Path)
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=REPO_ROOT / "cryodrgn/dashboard/static/dashboard_interactions_demo.gif",
    )
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument("--conda-env", type=str, default=None)
    parser.add_argument(
        "--frame-ms",
        type=int,
        default=FRAME_MS_30FPS,
        help="Frame duration (33 ≈ 30 fps)",
    )
    parser.add_argument(
        "--max-width",
        type=int,
        default=1008,
        help="Resize GIF to this width (0 = no resize)",
    )
    parser.add_argument(
        "--headed",
        action="store_true",
        help="Run Chromium with a visible window (optional; 3D GIF frames use server PNG composite, not WebGL)",
    )
    parser.add_argument(
        "--pairplot-companion",
        action="store_true",
        help="Also write dashboard_interactions_pairplot_after_selection.gif (long holds, separate file)",
    )
    parser.add_argument(
        "--pairplot-companion-output",
        type=Path,
        default=DEFAULT_PAIRPLOT_COMPANION,
        help="Output path for the pair-plot-after-selection companion GIF",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Minimal output (only final 'Wrote …' lines and errors)",
    )
    args = parser.parse_args()
    outdir = args.outdir.resolve()
    if not outdir.is_dir():
        print(f"error: not a directory: {outdir}", file=sys.stderr)
        return 1

    try:
        import PIL.Image  # noqa: F401
    except ImportError:
        print("error: install pillow", file=sys.stderr)
        return 1
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        print("error: install playwright + chromium", file=sys.stderr)
        return 1

    port = args.port or _free_port()
    base = f"http://127.0.0.1:{port}"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO_ROOT) + os.pathsep + env.get("PYTHONPATH", "")

    dash = [
        "cryodrgn",
        "dashboard",
        str(outdir),
        "--no-browser",
        "--host",
        "127.0.0.1",
        "--port",
        str(port),
        "--filter-max",
        "40000",
        "--cpus",
        "8",
    ]
    if args.conda_env:
        cmd = ["conda", "run", "-n", args.conda_env, "--no-capture-output", *dash]
    else:
        cryo = shutil.which("cryodrgn")
        if not cryo:
            print("error: `cryodrgn` not on PATH", file=sys.stderr)
            return 1
        cmd = [cryo, *dash[1:]]
    log = _make_logger(quiet=args.quiet)
    proc = subprocess.Popen(
        cmd,
        env=env,
        cwd=str(REPO_ROOT),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        text=True,
    )
    try:
        log(f"Starting dashboard subprocess on port {port} (outdir={outdir})")
        _wait_http(base + "/", timeout_s=240.0, log=log)
        buf = FrameBuffer(args.frame_ms)
        ph, launch_args = _chromium_launch_options(headed=args.headed)
        with sync_playwright() as p:
            with _timed_step(
                log,
                f"Launch Chromium (headless={ph})",
                slow_note=False,
            ):
                browser = p.chromium.launch(headless=ph, args=launch_args)
            context = browser.new_context(
                viewport={"width": 1440, "height": 900},
                device_scale_factor=1,
            )
            page = context.new_page()
            with _timed_step(
                log, "Main capture sequence (explorer → 3D → pair plot)", slow_note=True
            ):
                record_sequence(page, base, buf, log=log)
            if args.pairplot_companion:
                buf_pf = FrameBuffer(args.frame_ms)
                with _timed_step(log, "Companion pair-plot sequence", slow_note=True):
                    record_pairplot_after_selection_sequence(
                        page, base, buf_pf, log=log
                    )
                nf, wall_s_pf = _save_buf_to_gif(
                    buf_pf,
                    args.pairplot_companion_output,
                    frame_ms=args.frame_ms,
                    max_width=args.max_width,
                    max_wall_ms=COMPANION_MAX_WALL_MS,
                    log=log,
                    label="companion",
                )
                print(
                    f"Wrote {args.pairplot_companion_output} ({nf} frames × {args.frame_ms} ms "
                    f"≈ {wall_s_pf:.2f} s @ {1000/args.frame_ms:.0f} fps)"
                )
            browser.close()

        if not buf.pngs:
            raise RuntimeError("no frames captured")

        nf, wall_s = _save_buf_to_gif(
            buf,
            args.output,
            frame_ms=args.frame_ms,
            max_width=args.max_width,
            max_wall_ms=MAX_GIF_WALL_MS,
            log=log,
            label="main",
        )
        print(
            f"Wrote {args.output} ({nf} frames × {args.frame_ms} ms ≈ {wall_s:.2f} s @ {1000/args.frame_ms:.0f} fps)"
        )
        return 0
    except Exception as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=15)
        except subprocess.TimeoutExpired:
            proc.kill()


if __name__ == "__main__":
    raise SystemExit(main())
