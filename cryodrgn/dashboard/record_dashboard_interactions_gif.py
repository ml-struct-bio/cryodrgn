#!/usr/bin/env python3
"""Headless capture: particle explorer (lasso + image cache), 3D latent view, pair-plot covariates.

Uses a **50 ms** frame interval by default (~20 fps) for a slower, easier-to-follow GIF. Everything is
written to **one** output GIF (default ``dashboard_interactions_demo.gif``), capped at **30 s** wall
time. Pair-plot segments wait for the PNG grid to finish, then hold on the finished grid for the same
duration as the preceding ``Rendering…`` clip (so the two read as equal-length halves of the update).

Before each major interface segment, the recorder shows a short clip on the **launch hub** (hover +
click the corresponding card) so the transition from hub → tool is visible in the GIF. Hub dwell
times are scaled by ``HUB_CLIP_TIME_SCALE``. For pair-plot and montage **loading** clips, polling uses
``PAIRPLOT_SNAP_MAX_STEPS`` / ``MONTAGE_PRELOAD_MAX_FRAMES`` with ``PAIRPLOT_RENDER_POLL_SNAPSHOT_STRIDE``
(twice the stride used elsewhere so the pair-plot ``Rendering…`` strip is half as long) plus
``RENDER_POLL_FRAME_DURATION_MULT``. Each kept frame uses ``frame_ms`` × ``RENDER_POLL_FRAME_DURATION_MULT`` in
the encoded GIF so those segments stay shorter overall but play a bit slower per frame.

**3D plot:** The GIF drives Plotly's ``scene.camera`` (orbit + zoom) and eases ``scene.camera.center``
along the data **z** axis between segments for a slow vertical pan layered on the existing motion.
The in-plot Plotly colour bar is hidden; the aside covariate legend (violin + threshold UI) is hidden
for the 3-D segment so the locked base frame stays pixel-stable next to the orbit.

Depends on the same setup as ``record_dashboard_demo_gif.py`` (playwright, pillow, chromium).

Example::

    PYTHONPATH=/path/to/cryodrgn_beta conda run -n cdrgn-4.2.0_py-3.13_beta \\
        python -m cryodrgn.dashboard.record_dashboard_interactions_gif \\
        /scratch/.../004_train-vae_1gpu_dim.1024/

Or from the repo root (same entrypoint)::

    python scripts/record_dashboard_interactions_gif.py /scratch/.../outdir/

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

# Repo root: …/cryodrgn/dashboard/this_file.py → parents[2] == checkout containing ``cryodrgn/``.
REPO_ROOT = Path(__file__).resolve().parents[2]
FRAME_MS_DEFAULT = 50
MAX_GIF_WALL_MS = 30_000
# Pair-plot: poll while overlay/skeleton is up, then block until PNG is ready.
PAIRPLOT_SNAP_MAX_STEPS = 75
PAIRPLOT_WAIT_MS = 300_000
# Montage preload: max screenshot polls before giving up (still waits below).
MONTAGE_PRELOAD_MAX_FRAMES = 20
# Launch hub: scale all dwell times in ``_hub_launch_clip`` (20% slower than 2.25×).
HUB_CLIP_TIME_SCALE = 2.8
# During render polling, keep every Nth screenshot (coarser = shorter loading segment in the GIF).
RENDER_POLL_SNAPSHOT_STRIDE = 4
# Pair-plot loading in this GIF: coarser snapshots so the ``Rendering…`` segment is half as long.
PAIRPLOT_RENDER_POLL_SNAPSHOT_STRIDE = RENDER_POLL_SNAPSHOT_STRIDE * 4
# GIF display time per kept render-poll frame (multiple of ``frame_ms``, e.g. 2 = half playback speed).
RENDER_POLL_FRAME_DURATION_MULT = 4
# Preview PNG is pasted over ``#latent3d``; keep this close to 1 so in-page Plotly colour bars
# do not show in the gutters around the matplotlib composite. (The aside covariate legend is Plotly —
# violin + threshold UI — see ``color_covariate_legend.js``; we hide it during the 3-D segment.)
LATENT3D_PREVIEW_PANEL_SCALE = 0.996
# ``base.html`` loads Plotly from this CDN URL. On compute nodes without outbound
# internet the script fails to load and the explorer shows an "Plotly is not
# defined" error in place of the scatter. We pre-cache the bundle once (typically
# on a login node where the CDN is reachable) and serve it via Playwright's route
# interception so subsequent runs work even when the compute node is offline.
PLOTLY_CDN_VERSION = "2.35.2"
PLOTLY_CDN_URL = f"https://cdn.plot.ly/plotly-{PLOTLY_CDN_VERSION}.min.js"
# ``**/plotly-<version>.min.js`` glob so the interception keys off the file name
# and version regardless of scheme/host (Playwright's ``route`` takes a glob or regex).
PLOTLY_CDN_ROUTE_GLOB = f"**/plotly-{PLOTLY_CDN_VERSION}.min.js"

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


def _plotly_cache_path() -> Path:
    """Return the per-user cache path for ``plotly-<version>.min.js``.

    Respects ``$XDG_CACHE_HOME`` if set, else falls back to ``~/.cache/``.
    """
    xdg = os.environ.get("XDG_CACHE_HOME")
    base = Path(xdg).expanduser() if xdg else Path.home() / ".cache"
    return base / "cryodrgn" / "dashboard_gif" / f"plotly-{PLOTLY_CDN_VERSION}.min.js"


def _ensure_plotly_cached(log) -> Path | None:
    """Return a path to a cached ``plotly.min.js``, downloading on demand.

    Returns ``None`` if the bundle is neither on disk nor fetchable — in that case
    the recorder lets Playwright hit the real CDN (and the GIF shows the explorer's
    "Plotly is not defined" error if the compute node has no outbound internet).
    """
    dst = _plotly_cache_path()
    if dst.exists() and dst.stat().st_size > 1024:
        return dst
    dst.parent.mkdir(parents=True, exist_ok=True)
    with _timed_step(log, f"Cache Plotly bundle → {dst}", slow_note=False):
        try:
            with urllib.request.urlopen(PLOTLY_CDN_URL, timeout=30) as r:
                data = r.read()
        except (urllib.error.URLError, OSError) as err:
            log(f"Plotly CDN unreachable ({err}); continuing without local cache")
            return None
        tmp = dst.with_suffix(dst.suffix + ".part")
        tmp.write_bytes(data)
        tmp.replace(dst)
    return dst


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
        # Parallel to ``pngs``: per-frame duration in the encoded GIF (ms).
        self.durations_ms: list[int] = []

    def snap(self, page, *, duration_ms: int | None = None) -> None:
        self.pngs.append(page.screenshot(type="png", full_page=False))
        self.durations_ms.append(
            int(duration_ms) if duration_ms is not None else self.frame_ms
        )

    def add_png(self, data: bytes, *, duration_ms: int | None = None) -> None:
        self.pngs.append(data)
        self.durations_ms.append(
            int(duration_ms) if duration_ms is not None else self.frame_ms
        )

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


def _snap_until_pairplot_ready(buf: FrameBuffer, page, max_steps: int = 120) -> int:
    """Poll until the pair-plot PNG is ready; return the total GIF duration (ms)
    of the ``Rendering…`` frames appended to ``buf`` so callers can mirror the
    length with a matching hold on the finished plot."""
    render_ms = max(1, int(buf.frame_ms * RENDER_POLL_FRAME_DURATION_MULT))
    rendering_ms = 0
    for step in range(max_steps):
        if step % PAIRPLOT_RENDER_POLL_SNAPSHOT_STRIDE == 0:
            buf.snap(page, duration_ms=render_ms)
            rendering_ms += render_ms
        if page.evaluate(_pairplot_ready_js()):
            return rendering_ms
        time.sleep(buf.frame_ms / 1000.0)
    _wait_pairplot_ready(page, 120_000)
    return rendering_ms


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
    locked_query: dict[str, str] | None = None,
) -> None:
    from PIL import Image

    params = (
        locked_query if locked_query is not None else _latent3d_query_from_dom(page)
    )
    if not params:
        buf.snap(page)
        return
    req_params = dict(params)
    req_params["elev"] = f"{float(elev):.4f}"
    req_params["azim"] = f"{float(azim):.4f}"
    req_params["colorbar"] = "0"
    # Keep PNG dimensions deterministic across elev/azim so the composited sweep
    # does not jitter along the plot edges (visible as flicker at the bright end
    # of the colour gradient for palettes like Cividis / Viridis).
    req_params["stable_size"] = "1"
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
            buf.add_png(shot)
            return
        box = host.bounding_box()
        if not box:
            buf.add_png(shot)
            return
        bx, by = int(box["x"]), int(box["y"])
        bw, bh = int(box["width"]), int(box["height"])
    pw, ph = preview.size
    if bw < 8 or bh < 8:
        buf.snap(page)
        return
    panel = Image.new("RGBA", (bw, bh), (250, 248, 244, 255))
    scale = min(bw / pw, bh / ph) * LATENT3D_PREVIEW_PANEL_SCALE * max(0.5, zoom)
    nw, nh = max(1, int(pw * scale)), max(1, int(ph * scale))
    pr = preview.resize((nw, nh), Image.Resampling.LANCZOS)
    ox = (bw - nw) // 2 + int(pan_x * bw * 0.14)
    oy = (bh - nh) // 2 + int(pan_y * bh * 0.14)
    panel.paste(pr, (ox, oy), pr)
    canvas.paste(panel, (bx, by), panel)
    out = io.BytesIO()
    canvas.convert("RGB").save(out, format="PNG")
    buf.add_png(out.getvalue())


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
    # Lock axis / colour / palette for the whole sweep. Re-reading the DOM every frame can
    # briefly miss a checked palette radio during Plotly refresh, which flips the preview
    # API default (Viridis) vs the real selection and reads as palette popping in the GIF.
    locked = _latent3d_query_from_dom(page)
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
            locked_query=locked,
        )
        time.sleep(buf.frame_ms / 1000.0)


def _plotly_camera_from_spherical(
    elev_deg: float,
    azim_deg: float,
    *,
    r: float = 2.1,
    center_xyz: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> dict:
    """Plotly ``scene.camera`` dict from matplotlib-style spherical coords.

    ``r`` is the camera distance in normalised data-space units (Plotly's default
    eye is ``(1.25, 1.25, 1.25)`` → r ≈ 2.17). ``center_xyz`` translates both eye
    and center for a pan.
    """
    elev = math.radians(elev_deg)
    azim = math.radians(azim_deg)
    cx, cy, cz = center_xyz
    return {
        "eye": {
            "x": cx + r * math.cos(elev) * math.cos(azim),
            "y": cy + r * math.cos(elev) * math.sin(azim),
            "z": cz + r * math.sin(elev),
        },
        "center": {"x": cx, "y": cy, "z": cz},
        "up": {"x": 0.0, "y": 0.0, "z": 1.0},
    }


def _set_latent3d_camera(page, camera: dict) -> None:
    """Drive the 3D visualizer's camera via ``Plotly.relayout`` (no transition)."""
    page.evaluate(
        """(cam) => {
          const gd = document.getElementById('latent3d');
          if (!gd || !window.Plotly) return false;
          Plotly.relayout(gd, {'scene.camera': cam});
          return true;
        }""",
        camera,
    )


def _hide_latent3d_plotly_colorbar(page) -> None:
    """Hide Plotly's in-plot colour bar on the 3D page.

    The dashboard already renders a covariate histogram / filter legend in the
    right-hand aside, so Plotly's default bar is redundant and clutters the GIF.
    Safe no-op on discrete covariates or when Plotly is not ready.
    """
    page.evaluate(
        """() => {
          const gd = document.getElementById('latent3d');
          if (!gd || !window.Plotly || !gd.data || !gd.data.length) return false;
          try { Plotly.restyle(gd, {'marker.showscale': false}, [0]); } catch (e) {}
          return true;
        }"""
    )
    # Give Plotly a moment to reflow after the restyle so the next base screenshot
    # does not capture an intermediate state of the plot SVG layer.
    time.sleep(0.05)


def _hide_latent3d_color_legend(page) -> None:
    """Hide the 3-D aside covariate legend panel during capture.

    The panel (``CryoColorCovariateLegend`` — Plotly violin, threshold chrome,
    collapsible palette) sits beside the GL plot; hiding it keeps the locked base
    frame stable across camera sweeps and avoids any aside reflow beside the orbit.
    """
    page.evaluate(
        """() => {
          const el = document.getElementById('latent3d-color-legend-panel');
          if (el) el.style.display = 'none';
          return !!el;
        }"""
    )


def _capture_latent3d_locked_base(
    page,
) -> tuple[bytes, tuple[int, int, int, int]] | None:
    """Snapshot the whole page plus the ``#latent3d`` bounding box for re-use.

    Intended to be called once before a sequence of :func:`sleep_snap_latent3d_pan`
    calls so that the non-plot chrome (aside legend is hidden for capture) stays
    byte-for-byte identical across every sweep segment, not
    just within a single call. Returns ``None`` if the plot host is missing or
    degenerate, in which case callers should fall back to per-sweep bases.
    """
    shot = page.screenshot(type="png", full_page=False)
    host = page.locator("#latent3d")
    if not host.count():
        return None
    box = host.bounding_box()
    if not box:
        return None
    bx = int(round(box["x"]))
    by = int(round(box["y"]))
    bw = int(round(box["width"]))
    bh = int(round(box["height"]))
    if bw < 8 or bh < 8:
        return None
    return shot, (bx, by, bw, bh)


def sleep_snap_latent3d_pan(
    buf: FrameBuffer,
    page,
    seconds: float,
    *,
    elev_start: float,
    elev_end: float,
    azim_start: float,
    azim_end: float,
    r_start: float = 2.1,
    r_end: float = 2.1,
    center_start: tuple[float, float, float] = (0.0, 0.0, 0.0),
    center_end: tuple[float, float, float] = (0.0, 0.0, 0.0),
    locked_base: tuple[bytes, tuple[int, int, int, int]] | None = None,
) -> None:
    """Plotly-driven camera pan on ``#latent3d``.

    A full-page base screenshot is captured up-front (or supplied via
    ``locked_base`` so that it can be shared across back-to-back sweep segments);
    each subsequent frame re-screenshots the page, crops the ``#latent3d``
    region, and pastes it onto the locked base. Everything outside the plot
    (including the aside covariate legend, which the recorder hides) is therefore
    pixel-identical across all frames.
    """
    from PIL import Image

    if locked_base is None:
        captured = _capture_latent3d_locked_base(page)
        if captured is None:
            buf.snap(page)
            return
        base_shot, (bx, by, bw, bh) = captured
    else:
        base_shot, (bx, by, bw, bh) = locked_base
    base_canvas = Image.open(io.BytesIO(base_shot)).convert("RGBA")
    n = max(1, int(round(seconds * 1000 / buf.frame_ms)))
    for i in range(n):
        t = i / max(1, n - 1)
        tw = 0.5 - 0.5 * math.cos(math.pi * t)
        elev = elev_start + (elev_end - elev_start) * tw
        azim = azim_start + (azim_end - azim_start) * tw
        r = r_start + (r_end - r_start) * tw
        cx = center_start[0] + (center_end[0] - center_start[0]) * tw
        cy = center_start[1] + (center_end[1] - center_start[1]) * tw
        cz = center_start[2] + (center_end[2] - center_start[2]) * tw
        cam = _plotly_camera_from_spherical(elev, azim, r=r, center_xyz=(cx, cy, cz))
        _set_latent3d_camera(page, cam)
        time.sleep(buf.frame_ms / 1000.0)
        shot = page.screenshot(type="png", full_page=False)
        frame = Image.open(io.BytesIO(shot)).convert("RGBA")
        plot_crop = frame.crop((bx, by, bx + bw, by + bh))
        out = base_canvas.copy()
        out.paste(plot_crop, (bx, by))
        buf_out = io.BytesIO()
        out.convert("RGB").save(buf_out, format="PNG")
        buf.add_png(buf_out.getvalue())


_HUB_CARD_TITLES: dict[str, str] = {
    "explorer": "Particle explorer",
    "latent_3d": "3D visualizer",
    "pairplot": "Pair-plot generator",
}
_HUB_DIRECT_PATHS: dict[str, str] = {
    "explorer": "/explorer",
    "latent_3d": "/latent-3d",
    "pairplot": "/pairplot",
}


def _hub_launch_clip(
    buf: FrameBuffer,
    page,
    base: str,
    interface_key: str,
    *,
    log,
) -> None:
    """On the dashboard hub: show cards, hover the right interface, click to navigate."""
    s = HUB_CLIP_TIME_SCALE
    title = _HUB_CARD_TITLES[interface_key]
    direct = _HUB_DIRECT_PATHS[interface_key]
    with _timed_step(
        log,
        f"Launch hub: highlight and open “{title}” (~{1.0 * s:.1f}s clip)",
        slow_note=False,
    ):
        page.goto(base + "/", wait_until="domcontentloaded", timeout=120_000)
        buf.sleep_snap(page, 0.35 * s)
        link = page.locator("a.landing-card-link").filter(
            has=page.locator("h2.landing-card-title", has_text=title)
        )
        if link.count() == 0:
            log(
                f"No active hub link for “{title}” (interface inactive?); navigating to {direct}"
            )
            page.goto(base + direct, wait_until="domcontentloaded", timeout=120_000)
            buf.sleep_snap(page, 0.65 * s)
            return
        card = link.first
        card.scroll_into_view_if_needed()
        buf.sleep_snap(page, 0.2 * s)
        card.hover()
        buf.sleep_snap(page, 0.12 * s)
        with page.expect_navigation(wait_until="domcontentloaded", timeout=120_000):
            card.click()
        buf.sleep_snap(page, 0.33 * s)


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
    _hub_launch_clip(buf, page, base, "explorer", log=log)
    with _timed_step(
        log, "Particle explorer: wait for scatter plot after hub navigation"
    ):
        _wait_scatter_ready(page, 180_000)
    buf.sleep_snap(page, 0.7)

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
        buf.sleep_snap(page, 0.7)

    with _timed_step(
        log,
        "Particle explorer: preload image montage after View images (server decode + grid)",
        slow_note=True,
    ):
        page.locator("#btn-view-images").click()
        render_ms = max(1, int(buf.frame_ms * RENDER_POLL_FRAME_DURATION_MULT))
        for step in range(MONTAGE_PRELOAD_MAX_FRAMES):
            if page.evaluate(_montage_cache_ready_js()):
                break
            if step % RENDER_POLL_SNAPSHOT_STRIDE == 0:
                buf.snap(page, duration_ms=render_ms)
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
    buf.sleep_snap(page, 1.6)

    # --- 3D visualizer: use Plotly throughout (intro + camera pan) ---
    _hub_launch_clip(buf, page, base, "latent_3d", log=log)
    with _timed_step(
        log, "3D latent page: wait for Plotly overlay after hub navigation"
    ):
        _wait_latent3d_ready(page, 180_000)
    _hide_latent3d_plotly_colorbar(page)
    # Hide the aside covariate legend for the entire 3-D portion (user sees it again
    # after navigating away). Keeps the locked base frame unchanged beside the orbit.
    _hide_latent3d_color_legend(page)
    buf.sleep_snap(page, 0.45)

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
        _hide_latent3d_plotly_colorbar(page)

    with _timed_step(log, "3D: Cividis palette for GIF capture", slow_note=False):
        # Palette lives in a collapsible control; radios stay in DOM. Setting
        # ``checked`` + dispatching ``change`` matches clicking an option when the
        # dropdown is closed (listeners call ``loadPlot`` / ``refreshHistogramColors``).
        page.evaluate(
            """() => {
              const r = document.querySelector('input[name="latent3d_palette"][value="Cividis"]');
              if (!r) return;
              r.checked = true;
              r.dispatchEvent(new Event('change', { bubbles: true }));
            }"""
        )
        _wait_latent3d_ready(page, 180_000)
        _hide_latent3d_plotly_colorbar(page)

    # Capture the base once and share it across every sweep segment so chrome
    # outside ``#latent3d`` (headers, side controls, hidden legend slot) stays
    # byte-identical for the entire 3D pan.
    locked_base = _capture_latent3d_locked_base(page)

    # Camera sweep uses Plotly directly (same renderer as the intro clip) so the
    # 3D part never switches plotting backends mid-animation. The sweep is split
    # into contiguous segments with matched start/end cameras so there are no
    # jumps between segments.
    sleep_snap_latent3d_pan(
        buf,
        page,
        0.65,
        elev_start=20.0,
        elev_end=24.0,
        azim_start=-72.0,
        azim_end=-62.0,
        r_start=2.10,
        r_end=2.06,
        center_start=(0.0, 0.0, 0.0),
        center_end=(0.0, 0.0, 0.04),
        locked_base=locked_base,
    )
    sleep_snap_latent3d_pan(
        buf,
        page,
        0.75,
        elev_start=24.0,
        elev_end=28.0,
        azim_start=-62.0,
        azim_end=-48.0,
        r_start=2.06,
        r_end=2.02,
        center_start=(0.0, 0.0, 0.04),
        center_end=(0.0, 0.0, 0.10),
        locked_base=locked_base,
    )
    sleep_snap_latent3d_pan(
        buf,
        page,
        0.65,
        elev_start=28.0,
        elev_end=26.0,
        azim_start=-48.0,
        azim_end=-35.0,
        r_start=2.02,
        r_end=1.98,
        center_start=(0.0, 0.0, 0.10),
        center_end=(0.0, 0.0, 0.16),
        locked_base=locked_base,
    )
    with _timed_step(
        log,
        "3D: in-page orbit + gentle zoom (Plotly camera sweep)",
        slow_note=True,
    ):
        # Duration trimmed to the first 80% of the originally-planned 3.2 s
        # orbit (i.e., last 20% of the whole 3D pan removed). End-state picked
        # to match where the original cos-eased trajectory would have been at
        # t = 0.672 * 3.2 s (tw ≈ 0.757), so the motion still eases out to a
        # soft stop rather than reaching the same endpoint faster.
        sleep_snap_latent3d_pan(
            buf,
            page,
            2.15,
            elev_start=26.0,
            elev_end=20.0,
            azim_start=-35.0,
            azim_end=13.0,
            r_start=1.98,
            r_end=1.87,
            center_start=(0.0, 0.0, 0.16),
            center_end=(0.0, 0.0, 0.0),
            locked_base=locked_base,
        )

    _hub_launch_clip(buf, page, base, "pairplot", log=log)
    with _timed_step(
        log,
        "Pair plot: initial grid render to completion (sample frames + wait)",
        slow_note=True,
    ):
        rendering_ms = _snap_until_pairplot_ready(
            buf, page, max_steps=PAIRPLOT_SNAP_MAX_STEPS
        )
        _wait_pairplot_ready(page, PAIRPLOT_WAIT_MS)
    buf.sleep_snap(page, max(0.001, rendering_ms / 1000.0))

    radios = page.locator('input[name="color_cov"]')
    rc = radios.count()
    if rc >= 2:
        with _timed_step(
            log,
            "Pair plot: first covariate change — render to completion",
            slow_note=True,
        ):
            r = radios.nth(min(1, rc - 1))
            r.scroll_into_view_if_needed()
            r.click(force=True)
            rendering_ms = _snap_until_pairplot_ready(
                buf, page, max_steps=PAIRPLOT_SNAP_MAX_STEPS
            )
            _wait_pairplot_ready(page, PAIRPLOT_WAIT_MS)
        buf.sleep_snap(page, max(0.001, rendering_ms / 1000.0))

    if rc >= 3:
        with _timed_step(
            log,
            "Pair plot: second covariate change — render to completion",
            slow_note=True,
        ):
            r2 = radios.nth(min(2, rc - 1))
            r2.scroll_into_view_if_needed()
            r2.click(force=True)
            rendering_ms = _snap_until_pairplot_ready(
                buf, page, max_steps=PAIRPLOT_SNAP_MAX_STEPS
            )
            _wait_pairplot_ready(page, PAIRPLOT_WAIT_MS)
        buf.sleep_snap(page, max(0.001, rendering_ms / 1000.0))
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
    label: str = "GIF",
) -> tuple[int, float]:
    from PIL import Image

    if not buf.pngs:
        raise RuntimeError("no frames to save")
    if len(buf.durations_ms) != len(buf.pngs):
        raise RuntimeError("internal error: durations_ms out of sync with pngs")
    with _timed_step(
        log,
        f"Encode GIF ({label}): resize {len(buf.pngs)} frames, write {output.name}",
        slow_note=True,
    ):
        pngs_out: list[bytes] = []
        durs_out: list[int] = []
        used_ms = 0
        for raw, d in zip(buf.pngs, buf.durations_ms):
            d = max(1, int(d))
            if used_ms + d > max_wall_ms:
                break
            pngs_out.append(raw)
            durs_out.append(d)
            used_ms += d
        if not pngs_out:
            pngs_out = [buf.pngs[0]]
            durs_out = [max(1, min(buf.durations_ms[0], max_wall_ms))]
        pil_frames: list[Image.Image] = []
        for raw in pngs_out:
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
        output.parent.mkdir(parents=True, exist_ok=True)
        pil_frames[0].save(
            output,
            save_all=True,
            append_images=pil_frames[1:],
            duration=durs_out,
            loop=0,
            optimize=True,
        )
    wall_s = sum(durs_out) / 1000.0
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
        default=FRAME_MS_DEFAULT,
        help=f"Frame duration in ms (default {FRAME_MS_DEFAULT} ≈ 20 fps)",
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
        plotly_cached = _ensure_plotly_cached(log)
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
            if plotly_cached is not None:
                # Serve the cached bundle in place of any request for
                # ``plotly-<version>.min.js`` so explorer/3D pages have
                # ``window.Plotly`` even on compute nodes without outbound
                # internet. The CDN URL in ``base.html`` still resolves to
                # the local file via ``context.route``.
                plotly_bytes = plotly_cached.read_bytes()

                def _fulfill_plotly(route):
                    route.fulfill(
                        status=200,
                        headers={
                            "content-type": "application/javascript; charset=utf-8",
                            "cache-control": "public, max-age=86400",
                        },
                        body=plotly_bytes,
                    )

                context.route(PLOTLY_CDN_ROUTE_GLOB, _fulfill_plotly)
            page = context.new_page()
            with _timed_step(
                log,
                "Full capture (explorer → 3D → pair plot, single GIF)",
                slow_note=True,
            ):
                record_sequence(page, base, buf, log=log)
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
        )
        mean_ms = (wall_s * 1000.0 / nf) if nf else 0.0
        print(
            f"Wrote {args.output} ({nf} frames, ≈ {wall_s:.2f} s playback, "
            f"mean {mean_ms:.1f} ms/frame)"
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
