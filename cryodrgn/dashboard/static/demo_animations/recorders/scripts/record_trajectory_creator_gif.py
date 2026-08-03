#!/usr/bin/env python3
"""Headless capture: trajectory creator VTK / ChimeraX / direct-trace demo.

Opens the trajectory creator directly (same destination as
``cryodrgn dashboard OUTDIR -c N --traj``) and records:

1. the default PC1 path (1 s)
2. decode/hydrate volumes, switch to VTK, orbit, cycle the slider (8 s)
3. colour covariate → ``znorm`` when available (1 s)
4. switch to ChimeraX and show / re-render volumes (8 s)
5. direct-trace mode: shrink to 3 points, move them onto three k-means
   centroids, densify to 8 points, then re-show in ChimeraX (12 s)

Long decode / ChimeraX waits are sampled for at most ``RENDER_CLIP_MAX_MS``
each, then polled silently so the GIF wall budget stays near 30 s.
Indeterminate progress bars are frozen (static fill + loading labels) so the
GIF does not flicker a sliding bar.

On GPU nodes, default Chromium flags use **ANGLE + GL/EGL** so the VTK viewer
can paint under headless Chromium against the NVIDIA device (no ``xvfb`` needed).
Vulkan ANGLE can enumerate the GPU yet leave WebGL unreadable; ``gl-egl`` both
draws and appears in Playwright screenshots. Pass ``--swiftshader`` to fall
back to software WebGL, or ``--headed`` under ``xvfb-run`` if preferred.

Example::

    conda run -n cdrgn_beta --no-capture-output \\
        python -m cryodrgn.dashboard.static.demo_animations.recorders.scripts.record_trajectory_creator_gif \\
        /scratch/.../run-1 --cpus 10
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from cryodrgn.dashboard.static.demo_animations.recorders.scripts.record_dashboard_interactions_gif import (
    DEMO_ANIMATIONS_RECORDERS_DIR,
    PLOTLY_CDN_ROUTE_GLOB,
    REPO_ROOT,
    FrameBuffer,
    _chromium_args_headless_webgl,
    _ensure_plotly_cached,
    _free_port,
    _install_explorer_scatter_recording_routes,
    _save_buf_to_gif_under_budget,
    _timed_step,
    _wait_http,
)


def _chromium_args_nvidia_egl_webgl() -> list[str]:
    """Headless Chromium flags that bind WebGL to an NVIDIA GPU via ANGLE GL/EGL.

    On this head node, ``--use-angle=vulkan`` reports the H100 in ``SystemInfo`` but
    yields an empty WebGL renderer and black ``readPixels`` / screenshots.
    ``--use-angle=gl-egl`` draws on the H100 and composites into Playwright captures
    without ``--disable-gpu-compositing``.
    """
    return [
        "--no-sandbox",
        "--disable-setuid-sandbox",
        "--disable-dev-shm-usage",
        "--ignore-gpu-blocklist",
        "--disable-gpu-sandbox",
        # Do not pass --disable-software-rasterizer: after a large volume hydrate
        # the NVIDIA GPU process can die, and without software fallback every
        # subsequent getContext('webgl') fails with vtk.js "no webgl context".
        "--enable-unsafe-swiftshader",
        "--enable-webgl",
        "--enable-webgl2",
        "--max-active-webgl-contexts=32",
        "--use-gl=angle",
        "--use-angle=gl-egl",
        "--headless=new",
    ]


def _trajectory_chromium_launch_options(
    *, headed: bool, swiftshader: bool
) -> tuple[bool, list[str]]:
    """Return ``(headless, args)`` for trajectory / VTK capture."""
    if swiftshader:
        return (not headed), _chromium_args_headless_webgl()
    args = list(_chromium_args_nvidia_egl_webgl())
    if headed:
        # Keep NVIDIA ANGLE GL/EGL flags under xvfb; only drop headless.
        args = [
            a
            for a in args
            if a != "--headless=new" and not str(a).startswith("--headless=")
        ]
        return False, args
    return True, args


_PRESERVE_WEBGL_DRAWING_BUFFER_JS = """
(() => {
  if (window.__cryoGifPreserveDrawingBuffer) return;
  window.__cryoGifPreserveDrawingBuffer = true;
  const orig = HTMLCanvasElement.prototype.getContext;
  HTMLCanvasElement.prototype.getContext = function (type, attrs) {
    if (type === 'webgl' || type === 'webgl2' || type === 'experimental-webgl') {
      attrs = Object.assign({}, attrs || {}, { preserveDrawingBuffer: true });
    }
    return orig.call(this, type, attrs);
  };
})();
"""

FRAME_MS_DEFAULT = 80  # ~12.5 fps (was 50 ms / 20 fps)
MAX_GIF_WALL_MS = 30_000
MAX_GIF_BYTES = 10 * 1024 * 1024
GIF_MAX_WIDTH_DEFAULT = 720
RENDER_CLIP_MAX_MS = 1_200
WAIT_TIMEOUT_MS = 600_000

SEG_PC1_S = 1.0
SEG_ZNORM_S = 1.0
SEG_VTK_S = 8.0
# Slow, short orbit with a gentle yaw / tip / twist.
SEG_VTK_ORBIT_S = 2.2
SEG_CHIMERAX_S = 8.0
SEG_DIRECT_S = 12.0
# Same dwell on every volume tick for VTK / ChimeraX / direct-trace scrubs.
# Prior dwell was (0.30 / 1.3) s/tick; another 45% speed-up → ÷ 1.45.
VOLUME_SLIDER_SEC_PER_TICK = (0.30 / 1.3) / 1.45
# Short calm arc (|x|,|z| ≪ 90°) plus a mild Z twist; both sliders use default view.
VTK_ORBIT_START_TURN_Y_DEG = 12.0
VTK_ORBIT_START_TURN_X_DEG = -4.0
VTK_ORBIT_START_TURN_Z_DEG = -10.0
VTK_ORBIT_END_TURN_Y_DEG = -10.0
VTK_ORBIT_END_TURN_X_DEG = 4.0
VTK_ORBIT_END_TURN_Z_DEG = 10.0
# Keep 40% of the post-slider segment pad on the final volume (60% shorter linger).
FINAL_VOLUME_LINGER_SCALE = 0.4
# Direct-trace: k-means waypoints then densified ChimeraX path.
DIRECT_KMEANS_K = 3
DIRECT_DENSIFY_N = 8

_LOG_PREFIX = "[record-trajectory-creator-gif]"

# Freeze sliding indeterminate bars; keep the track + loading label visible.
_FREEZE_PROGRESS_CSS = """
.cryo-plot-rendering-overlay__indeterminate,
.pairplot-skel-indeterminate {
  animation: none !important;
  left: 31% !important;
  transform: none !important;
}
"""

# Recording layout: keep the square volume viewport sized, hide empty ChimeraX
# <img> broken-icons, and suppress the VTK/slice flash while Decode forces VTK
# underneath a ChimeraX job (body[data-gif-chimerax-decode="1"]).
_GIF_VOLUME_VIEWER_CSS = """
body.cryo-gif-recording .cryo-traj-vol-panel-host .cryo-vslice-viewport {
  min-width: 300px !important;
  min-height: 300px !important;
  width: 100% !important;
  aspect-ratio: 1 / 1 !important;
  flex: 0 0 auto !important;
}
body.cryo-gif-recording #vslice-chimerax-preview:not([src]),
body.cryo-gif-recording #vslice-chimerax-preview[src=""],
body.cryo-gif-recording #vslice-chimerax-preview[src="about:blank"] {
  opacity: 0 !important;
  visibility: hidden !important;
  pointer-events: none !important;
}
.cryo-vslice-panel-shell[data-vol-backend="chimerax"] #vslice-canvas,
.cryo-vslice-panel-shell[data-vol-backend="chimerax"] #vslice-vtk-container,
body[data-gif-chimerax-decode="1"] #vslice-canvas,
body[data-gif-chimerax-decode="1"] #vslice-vtk-container,
body[data-gif-chimerax-decode="1"] #vslice-chimerax-preview {
  opacity: 0 !important;
  visibility: hidden !important;
  pointer-events: none !important;
}
.cryo-vslice-panel-shell[data-vol-backend="chimerax"]
  #vslice-rendering-overlay.cryo-plot-rendering-overlay--nonblocking.cryo-plot-rendering-overlay--show,
body[data-gif-chimerax-decode="1"]
  #vslice-rendering-overlay.cryo-plot-rendering-overlay--show {
  inset: 0 !important;
  top: 0 !important;
  right: 0 !important;
  bottom: 0 !important;
  left: 0 !important;
  width: auto !important;
  max-width: none !important;
  min-width: 0 !important;
  height: auto !important;
  align-items: center !important;
  justify-content: center !important;
  background: rgba(250, 248, 244, 0.96) !important;
  box-shadow: none !important;
  border: none !important;
  border-radius: 0 !important;
  pointer-events: auto !important;
  display: flex !important;
}
"""


def _make_logger(quiet: bool):
    def log(msg: str) -> None:
        if not quiet:
            print(f"{_LOG_PREFIX} {msg}", flush=True)

    return log


def _buf_wall_s(buf: FrameBuffer) -> float:
    return sum(buf.durations_ms) / 1000.0


def _hold_until(
    buf: FrameBuffer, page, target_wall_s: float, *, scale: float = 1.0
) -> None:
    """Sleep-snap toward ``target_wall_s`` (optionally a fraction of the remaining gap)."""
    remaining = (target_wall_s - _buf_wall_s(buf)) * float(scale)
    if remaining > 0.02:
        buf.sleep_snap(page, remaining)


def _snap_until_ready_capped(
    buf: FrameBuffer,
    page,
    ready_js: str,
    *,
    max_clip_ms: int = RENDER_CLIP_MAX_MS,
    timeout_ms: int = WAIT_TIMEOUT_MS,
) -> None:
    """Capture at most ``max_clip_ms`` while polling, then wait silently."""
    deadline = time.monotonic() + timeout_ms / 1000.0
    captured_ms = 0
    while time.monotonic() < deadline:
        if page.evaluate(ready_js):
            return
        if captured_ms < max_clip_ms:
            buf.snap(page)
            captured_ms += buf.frame_ms
        time.sleep(buf.frame_ms / 1000.0)
    page.wait_for_function(ready_js, timeout=5_000)


def _scatter_ready_js() -> str:
    return """() => {
      const overlay = document.getElementById('scatter-rendering-overlay');
      const overlayHidden = overlay
        && !overlay.classList.contains('cryo-plot-rendering-overlay--show');
      const gd = document.getElementById('scatter');
      const x = gd && gd.data && gd.data[0] ? gd.data[0].x : null;
      const P = window.CryoPlotlyArrays;
      const n = (P && x != null) ? P.length(x) : (x && x.length) || 0;
      const plotlyOk = typeof Plotly !== 'undefined';
      return !!(overlayHidden && plotlyOk && n > 0);
    }"""


def _pc1_path_ready_js() -> str:
    return """() => {
      const overlay = document.getElementById('traj-glyph-overlay');
      const hasGlyph = !!(overlay && overlay.querySelector('.cryo-traj-glyph-path'));
      const markers = overlay
        ? overlay.querySelectorAll('.cryo-traj-glyph-marker').length
        : 0;
      const volBtns = document.querySelectorAll('.cryo-vslice-vol-btn').length;
      return hasGlyph && markers >= 2 && volBtns >= 2;
    }"""


def _volumes_busy_js() -> str:
    """True while Decode/render is in flight (button label or volume overlay)."""
    return """() => {
      const btn = document.getElementById('btn-generate-volumes');
      const txt = btn ? (btn.textContent || '') : '';
      if (/decoding|rendering/i.test(txt)) return true;
      const overlay = document.getElementById('vslice-rendering-overlay');
      return !!(overlay && !overlay.hidden);
    }"""


def _volumes_ready_js() -> str:
    """True when decode/render finished (ChimeraX preview or interactive volumes)."""
    return """() => {
      const btn = document.getElementById('btn-generate-volumes');
      const txt = btn ? (btn.textContent || '') : '';
      if (/decoding|rendering/i.test(txt)) return false;
      const overlay = document.getElementById('vslice-rendering-overlay');
      if (overlay && !overlay.hidden) return false;
      const preview = document.getElementById('vslice-chimerax-preview');
      if (preview && preview.src && !preview.hidden) return true;
      const d = window.__cryoTrajVolDisplay;
      const vols = d && Array.isArray(d.volumes) ? d.volumes : [];
      if (vols.some(function (v) { return v && v.volume_b64; })) return true;
      const title = btn ? (btn.title || '') : '';
      return !!(btn && btn.disabled
        && /already (available|match)|volumes and chimeraX images are already available|already loaded/i.test(title));
    }"""


def _vtk_volume_ready_js() -> str:
    """True when the VTK raycast view holds a volume actor."""
    return """() => {
      const d = window.__cryoTrajVolDisplay;
      return !!(d && d.raycastView && d.raycastView.volume);
    }"""


def _click_generate_and_await(
    buf: FrameBuffer,
    page,
    *,
    expect_ticks: int | None = None,
    log=None,
) -> None:
    """Click Decode/render and wait for completion without trusting a stale preview.

    When a ChimeraX preview already exists, ``_volumes_ready_js`` can return true
    before the new job flips the button into a busy state. Require busy first in
    that case, then poll for ready (and optionally ``expect_ticks`` slider cells).
    """
    page.wait_for_function(
        """() => {
          const btn = document.getElementById('btn-generate-volumes');
          return !!(btn && !btn.hidden && !btn.disabled);
        }""",
        timeout=WAIT_TIMEOUT_MS,
    )
    had_preview = bool(
        page.evaluate(
            """() => {
              const preview = document.getElementById('vslice-chimerax-preview');
              return !!(preview && preview.src);
            }"""
        )
    )
    page.click("#btn-generate-volumes")
    if had_preview:
        page.wait_for_function(_volumes_busy_js(), timeout=120_000)
        if log:
            log("Decode/render busy after Generate (cleared stale-preview race)")
    else:
        try:
            page.wait_for_function(_volumes_busy_js(), timeout=15_000)
        except Exception:
            pass
    _snap_until_ready_capped(buf, page, _volumes_ready_js())
    if expect_ticks is not None:
        n = int(expect_ticks)
        page.wait_for_function(
            f"""() => {{
              const host = document.getElementById('traj-vol-volume-slider-ticks');
              if (!host) return false;
              return host.querySelectorAll('[data-vol-index]').length >= {n};
            }}""",
            timeout=WAIT_TIMEOUT_MS,
        )


def _vtk_active_js() -> str:
    return """() => {
      const vtk = document.getElementById('traj-vol-backend-vtk');
      const host = document.getElementById('vslice-vtk-container');
      return !!(vtk && vtk.checked && !vtk.disabled && host && !host.hidden);
    }"""


def _vtk_canvas_painted_js() -> str:
    """True when VTK host is visible with a sized canvas and a volume actor.

    vtk.js does not set ``preserveDrawingBuffer``, so ``readPixels`` after
    compositing often returns zeros even when Playwright screenshots show the
    volume. Treat a visible, non-trivial canvas plus ``raycastView.volume`` as
    painted.
    """
    return """() => {
      const d = window.__cryoTrajVolDisplay;
      if (!(d && d.raycastView && d.raycastView.volume)) return false;
      const host = document.getElementById('vslice-vtk-container');
      if (!host || host.hidden) return false;
      if (host.clientWidth < 32 || host.clientHeight < 32) return false;
      const canvases = Array.from(host.querySelectorAll('canvas'));
      for (const c of canvases) {
        if (c.width >= 64 && c.height >= 64) return true;
      }
      return false;
    }"""


def _webgl_renderer_info(page, *, lose: bool = True) -> dict:
    return page.evaluate(
        """(lose) => {
          const c = document.createElement('canvas');
          c.width = 64; c.height = 64;
          const gl = c.getContext('webgl2') || c.getContext('webgl');
          if (!gl) return { ok: false };
          const dbg = gl.getExtension('WEBGL_debug_renderer_info');
          const info = {
            ok: true,
            vendor: dbg ? gl.getParameter(dbg.UNMASKED_VENDOR_WEBGL) : gl.getParameter(gl.VENDOR),
            renderer: dbg ? gl.getParameter(dbg.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER),
          };
          if (lose) {
            const ext = gl.getExtension('WEBGL_lose_context');
            if (ext) ext.loseContext();
          }
          return info;
        }""",
        bool(lose),
    )


def _vtk_debug_state(page) -> dict:
    """Snapshot VTK / WebGL state for recorder diagnostics."""
    return page.evaluate(
        """() => {
          const d = window.__cryoTrajVolDisplay;
          const host = document.getElementById('vslice-vtk-container');
          const canvases = host
            ? Array.from(host.querySelectorAll('canvas')).map(function (c) {
                return { w: c.width, h: c.height };
              })
            : [];
          const vols = d && Array.isArray(d.volumes) ? d.volumes : [];
          let webglOk = false;
          try {
            const c = document.createElement('canvas');
            c.width = 8; c.height = 8;
            webglOk = !!(c.getContext('webgl2') || c.getContext('webgl'));
          } catch (err) { webglOk = false; }
          return {
            backend: d ? d.backend : null,
            hasRaycast: !!(d && d.raycastView),
            hasVolume: !!(d && d.raycastView && d.raycastView.volume),
            nB64: vols.filter(function (v) { return v && v.volume_b64; }).length,
            hostHidden: !!(host && host.hidden),
            hostSize: host ? [host.clientWidth, host.clientHeight] : null,
            canvases: canvases,
            webglOk: webglOk,
            lastErr: d && d._lastVtkError ? String(d._lastVtkError) : null,
          };
        }"""
    )


def _chimerax_active_js() -> str:
    return """() => {
      const cx = document.getElementById('traj-vol-backend-chimerax');
      const preview = document.getElementById('vslice-chimerax-preview');
      return !!(cx && cx.checked && !cx.disabled
        && preview && !preview.hidden && preview.src);
    }"""


def _ticks_active_js() -> str:
    return """() => {
      const host = document.getElementById('traj-vol-volume-slider-ticks');
      if (!host) return false;
      const ticks = Array.from(host.querySelectorAll('[data-vol-index]'));
      return ticks.length > 0 && ticks.every(function (t) {
        return !t.disabled
          && !t.classList.contains('cryo-vslice-volume-slider-tick--inactive');
      });
    }"""


def _set_select_value(page, select_id: str, value: str) -> None:
    page.evaluate(
        """([id, val]) => {
          const el = document.getElementById(id);
          if (!el) throw new Error('missing select #' + id);
          el.value = val;
          el.dispatchEvent(new Event('change', { bubbles: true }));
        }""",
        [select_id, value],
    )


def _choose_znorm(page) -> str:
    opts = page.eval_on_selector(
        "#sc",
        """(sel) => Array.from(sel.options).map((o) => ({
          value: o.value,
          label: (o.textContent || '').trim(),
          disabled: !!o.disabled,
        }))""",
    )
    for o in opts:
        if o.get("disabled"):
            continue
        value = str(o.get("value") or "")
        label = str(o.get("label") or "")
        if value == "znorm" or label.lower().startswith("znorm"):
            _set_select_value(page, "sc", value)
            return value
    # Prefer any non-empty colour covariate when znorm is absent.
    for o in opts:
        if o.get("disabled"):
            continue
        value = str(o.get("value") or "")
        if value and value.lower() not in ("none", ""):
            _set_select_value(page, "sc", value)
            return value
    raise RuntimeError(
        "No usable colour covariate in #sc; options: "
        + ", ".join(f"{o.get('value')} ({o.get('label')})" for o in opts)
    )


def _select_backend(page, backend: str) -> None:
    radio_id = f"traj-vol-backend-{backend}"
    page.evaluate(
        """(rid) => {
          const el = document.getElementById(rid);
          if (!el || el.disabled) throw new Error('backend unavailable: ' + rid);
          el.checked = true;
          el.dispatchEvent(new Event('change', { bubbles: true }));
          // Decode temporarily forces trajVolDisplay.backend to "vtk"; keep the
          // display object aligned with the radio after a combined ChimeraX job.
          const d = window.__cryoTrajVolDisplay;
          const name = rid.replace('traj-vol-backend-', '');
          if (d && typeof d.setBackend === 'function') d.setBackend(name);
        }""",
        radio_id,
    )


def _chimerax_images_ready_js(n: int) -> str:
    """True when at least ``n`` nonempty ChimeraX preview slots are present."""
    n = int(n)
    return f"""() => {{
      const d = window.__cryoTrajVolDisplay;
      const imgs = d && Array.isArray(d.chimeraxImages) ? d.chimeraxImages : [];
      return imgs.filter(Boolean).length >= {n};
    }}"""


def _freeze_progress_bar_motion(page) -> None:
    """Keep progress tracks + labels, but stop the sliding indeterminate fill."""
    page.add_style_tag(content=_FREEZE_PROGRESS_CSS)
    page.add_style_tag(content=_GIF_VOLUME_VIEWER_CSS)
    page.evaluate("""() => { document.body.classList.add('cryo-gif-recording'); }""")
    # vtk.js defaults preserveDrawingBuffer=false, so Playwright screenshots of
    # the VTK canvas are often empty even when the volume actor is live.
    page.evaluate(_PRESERVE_WEBGL_DRAWING_BUFFER_JS)


def _set_chimerax_decode_veil(page, enabled: bool) -> None:
    """Hide VTK/slice/empty preview while Decode forces an intermediate VTK paint."""
    page.evaluate(
        """(on) => {
          if (on) document.body.setAttribute('data-gif-chimerax-decode', '1');
          else document.body.removeAttribute('data-gif-chimerax-decode');
        }""",
        bool(enabled),
    )


def _zoom_scatter_to_cloud(
    page, *, pad_frac: float = 0.04, blend: float = 0.5, log=None
) -> None:
    """Zoom the scatter between the original autorange and a tight cloud fit.

    ``blend=0`` keeps the pre-zoom (original) ranges; ``blend=1`` uses the tight
    cloud framing (``pad_frac``). The demo default ``blend=0.5`` sits halfway.
    """
    page.wait_for_function(
        """() => {
          const gd = document.getElementById('scatter');
          return !!(gd && window.Plotly && gd._fullLayout
            && gd._fullLayout.xaxis && gd._fullLayout.yaxis
            && Array.isArray(gd.data) && gd.data.length
            && Array.isArray(gd._fullLayout.xaxis.range)
            && Array.isArray(gd._fullLayout.yaxis.range));
        }""",
        timeout=180_000,
    )
    ok = page.evaluate(
        """({ pad, blend }) => {
          const gd = document.getElementById('scatter');
          if (!gd || !window.Plotly || !Array.isArray(gd.data)) return false;
          const xa = gd._fullLayout.xaxis;
          const ya = gd._fullLayout.yaxis;
          const origX = xa.range.slice();
          const origY = ya.range.slice();
          const xs = [];
          const ys = [];
          gd.data.forEach((tr) => {
            if (!tr || !tr.x || !tr.y) return;
            const n = Math.min(tr.x.length, tr.y.length);
            for (let i = 0; i < n; i++) {
              const x = Number(tr.x[i]);
              const y = Number(tr.y[i]);
              if (Number.isFinite(x) && Number.isFinite(y)) {
                xs.push(x);
                ys.push(y);
              }
            }
          });
          if (xs.length < 8) return false;
          xs.sort((a, b) => a - b);
          ys.sort((a, b) => a - b);
          const q = (arr, t) => {
            const i = Math.max(0, Math.min(arr.length - 1, Math.floor(t * (arr.length - 1))));
            return arr[i];
          };
          // Tight cloud: trim outliers, pad, then square-ish framing.
          let xmin = q(xs, 0.01);
          let xmax = q(xs, 0.99);
          let ymin = q(ys, 0.01);
          let ymax = q(ys, 0.99);
          const dx = Math.max(1e-3, xmax - xmin);
          const dy = Math.max(1e-3, ymax - ymin);
          xmin -= dx * pad;
          xmax += dx * pad;
          ymin -= dy * pad;
          ymax += dy * pad;
          const cx = 0.5 * (xmin + xmax);
          const cy = 0.5 * (ymin + ymax);
          const half = 0.5 * Math.max(xmax - xmin, ymax - ymin);
          xmin = cx - half;
          xmax = cx + half;
          ymin = cy - half;
          ymax = cy + half;
          const b = Math.max(0, Math.min(1, Number(blend)));
          const lerp = (a, c) => a + (c - a) * b;
          Plotly.relayout(gd, {
            'xaxis.autorange': false,
            'yaxis.autorange': false,
            'xaxis.range': [lerp(origX[0], xmin), lerp(origX[1], xmax)],
            'yaxis.range': [lerp(origY[0], ymin), lerp(origY[1], ymax)],
          });
          return true;
        }""",
        {"pad": float(pad_frac), "blend": float(blend)},
    )
    if log:
        log(
            f"Zoomed scatter (blend={blend:.2f} toward cloud)"
            if ok
            else "WARNING: scatter zoom skipped"
        )
    page.wait_for_timeout(250)


def _volume_slot_ready_js(index: int, *, backend: str | None = None) -> str:
    """True when slider slot ``index`` has paintable material for the backend."""
    idx = int(index)
    be = "null" if backend is None else json.dumps(str(backend))
    return f"""() => {{
      const d = window.__cryoTrajVolDisplay;
      if (!d) return false;
      const backend = {be} || d.backend || '';
      if (backend === 'chimerax') {{
        const imgs = Array.isArray(d.chimeraxImages) ? d.chimeraxImages : [];
        return !!(imgs[{idx}]);
      }}
      if (typeof d._volumeReadyAt === 'function') return !!d._volumeReadyAt({idx});
      const vols = Array.isArray(d.volumes) ? d.volumes : [];
      return !!(vols[{idx}] && vols[{idx}].volume_b64);
    }}"""


def _wait_volume_slot_ready(
    page, index: int, *, backend: str | None = None, timeout_ms: int = 180_000
) -> bool:
    try:
        page.wait_for_function(
            _volume_slot_ready_js(index, backend=backend), timeout=timeout_ms
        )
        return True
    except Exception:
        return False


def _count_ready_volumes(page, *, backend: str | None = None) -> int:
    be = backend
    return int(
        page.evaluate(
            """(backend) => {
              const d = window.__cryoTrajVolDisplay;
              if (!d) return 0;
              const b = backend || d.backend || '';
              if (b === 'chimerax') {
                const imgs = Array.isArray(d.chimeraxImages) ? d.chimeraxImages : [];
                return imgs.filter(Boolean).length;
              }
              const vols = Array.isArray(d.volumes) ? d.volumes : [];
              if (typeof d._volumeReadyAt === 'function') {
                let n = 0;
                const total = typeof d._volumeNavCount === 'function'
                  ? d._volumeNavCount() : vols.length;
                for (let i = 0; i < total; i++) if (d._volumeReadyAt(i)) n++;
                return n;
              }
              return vols.filter(function (v) { return v && v.volume_b64; }).length;
            }""",
            be,
        )
        or 0
    )


def _await_all_volume_slots(
    page,
    n: int,
    *,
    backend: str | None = None,
    timeout_ms: int = WAIT_TIMEOUT_MS,
    log=None,
    label: str = "volumes",
) -> bool:
    """Wait until ``n`` slots are ready; return False if the budget expires."""
    n = max(1, int(n))
    deadline = time.monotonic() + timeout_ms / 1000.0
    while time.monotonic() < deadline:
        ready = _count_ready_volumes(page, backend=backend)
        if ready >= n:
            if log:
                log(f"{label}: {ready}/{n} slots ready")
            return True
        time.sleep(0.5)
    ready = _count_ready_volumes(page, backend=backend)
    if log:
        log(f"WARNING: {label} incomplete ({ready}/{n}) after wait")
    return ready >= n


def _ensure_chimerax_images(
    buf: FrameBuffer,
    page,
    n: int,
    *,
    log,
    label: str,
    timeout_ms: int = WAIT_TIMEOUT_MS,
) -> None:
    """Generate / wait until ``n`` ChimeraX preview slots exist; retry once."""
    n = max(1, int(n))
    if _count_ready_volumes(page, backend="chimerax") >= n:
        return
    _select_backend(page, "chimerax")
    _set_chimerax_decode_veil(page, True)
    try:
        for attempt in range(1, 3):
            if _count_ready_volumes(page, backend="chimerax") >= n:
                break
            if log:
                log(f"{label}: ChimeraX Generate attempt {attempt}/2 (need {n})")
            _click_generate_and_await(buf, page, expect_ticks=n, log=log)
            if _await_all_volume_slots(
                page,
                n,
                backend="chimerax",
                timeout_ms=timeout_ms,
                log=log,
                label=label,
            ):
                break
        _select_backend(page, "chimerax")
        try:
            page.wait_for_function(_chimerax_active_js(), timeout=120_000)
        except Exception:
            if log:
                log(f"{label}: ChimeraX preview not immediately active; continuing")
        try:
            page.wait_for_function(_ticks_active_js(), timeout=120_000)
        except Exception:
            if log:
                log(f"{label}: ChimeraX ticks not all active yet; continuing")
    finally:
        _set_chimerax_decode_veil(page, False)
        _select_backend(page, "chimerax")


def _reset_vtk_view(page) -> None:
    """Destroy any prior vtk.js viewer so the next activate gets a fresh canvas."""
    page.evaluate(
        """() => {
          const d = window.__cryoTrajVolDisplay;
          if (d && d.raycastView && typeof d.raycastView.destroy === 'function') {
            try { d.raycastView.destroy(); } catch (err) { /* ignore */ }
          }
          if (d) {
            d.raycastView = null;
            d.raycastVolIndex = null;
          }
          const host = document.getElementById('vslice-vtk-container');
          if (host) {
            host.querySelectorAll('canvas').forEach(function (c) {
              try {
                const gl = c.getContext('webgl2') || c.getContext('webgl')
                  || c.getContext('experimental-webgl');
                if (!gl) return;
                const ext = gl.getExtension('WEBGL_lose_context');
                if (ext) ext.loseContext();
              } catch (err) { /* ignore */ }
            });
            host.innerHTML = '';
          }
        }"""
    )


def _install_light_vtk_hydrate_route(context, *, max_volumes: int = 1, log=None):
    """Truncate ``volumes_from_cache`` payloads so VTK WebGL is not OOM'd.

    Full GTPase cache hydrates (~10×128³) often kill headless ANGLE before
    vtk.js can create a context. Keep only the first ``max_volumes`` blobs for
    the initial VTK orbit; later slider hops can fetch the rest one at a time.
    """
    max_volumes = max(1, int(max_volumes))

    def _route(route):
        try:
            post = route.request.post_data or ""
        except Exception:
            post = ""
        if "volumes_from_cache" not in post:
            route.continue_()
            return
        try:
            response = route.fetch()
            data = response.json()
        except Exception:
            route.continue_()
            return
        vols = data.get("volumes") if isinstance(data, dict) else None
        if not isinstance(vols, list) or len(vols) <= max_volumes:
            route.fulfill(status=response.status, json=data)
            return
        data = dict(data)
        data["volumes"] = vols[:max_volumes]
        slots = data.get("slot_indices")
        if isinstance(slots, list) and slots:
            data["slot_indices"] = slots[:max_volumes]
        if log:
            log(
                f"Light VTK hydrate: truncated cache payload "
                f"{len(vols)} → {max_volumes} volume(s)"
            )
        route.fulfill(status=response.status, json=data)

    context.route("**/api/trajectory_volumes", _route)
    return lambda: context.unroute("**/api/trajectory_volumes", _route)


def _activate_vtk(page, *, already_hydrated: bool = False) -> None:
    """Switch to VTK after ``volume_b64`` is present (avoid heavy context resets)."""
    page.wait_for_function(
        """() => {
          const vtk = document.getElementById('traj-vol-backend-vtk');
          return !!(vtk && !vtk.disabled);
        }""",
        timeout=WAIT_TIMEOUT_MS,
    )
    if not already_hydrated:
        page.click("label[for='traj-vol-backend-slice']")
        page.wait_for_function(
            """() => {
              const d = window.__cryoTrajVolDisplay;
              const vols = d && Array.isArray(d.volumes) ? d.volumes : [];
              return vols.some(function (v) { return v && v.volume_b64; });
            }""",
            timeout=WAIT_TIMEOUT_MS,
        )
    # Only destroy an existing raycast view; do not touch unrelated canvases.
    page.evaluate(
        """() => {
          const d = window.__cryoTrajVolDisplay;
          if (d && d.raycastView && typeof d.raycastView.destroy === 'function') {
            try { d.raycastView.destroy(); } catch (err) { /* ignore */ }
            d.raycastView = null;
            d.raycastVolIndex = null;
          }
          const preview = document.getElementById('vslice-chimerax-preview');
          if (preview) {
            preview.hidden = true;
            preview.removeAttribute('src');
            preview.src = '';
          }
        }"""
    )
    time.sleep(0.2)
    page.click("label[for='traj-vol-backend-vtk']")
    page.wait_for_function(_vtk_active_js(), timeout=WAIT_TIMEOUT_MS)


def _lose_webgl_contexts(page) -> None:
    """Best-effort release of abandoned VTK WebGL state before a retry."""
    _reset_vtk_view(page)


def _force_vtk_render(page) -> None:
    page.evaluate(
        """() => {
          const d = window.__cryoTrajVolDisplay;
          if (!d) return;
          // Hide any stale ChimeraX <img> (empty src shows a broken-image icon).
          const preview = document.getElementById('vslice-chimerax-preview');
          if (preview) {
            preview.hidden = true;
            preview.removeAttribute('src');
            preview.src = '';
          }
          const host = document.getElementById('vslice-vtk-container');
          if (host) host.hidden = false;
          // _renderVtk refuses to paint while vtkFocusIndex points at an unready
          // slot (even if a neighbour already has volume_b64). Snap to the first
          // ready volume so the initial orbit can start.
          if (typeof d._volumeReadyAt === 'function'
              && typeof d._volumeNavCount === 'function') {
            const total = d._volumeNavCount() || 0;
            let focus = Math.floor(Number(d.vtkFocusIndex)) || 0;
            if (total > 0 && !d._volumeReadyAt(focus)) {
              for (let i = 0; i < total; i++) {
                if (d._volumeReadyAt(i)) { focus = i; break; }
              }
              d.vtkFocusIndex = focus;
            }
          }
          // Never bump _vtkPaintGen while a module load / paint is already in
          // flight — SwiftShader vtk.js can take >60 s, and a second _renderVtk
          // aborts the almost-ready promise.
          if (d.vtkModulePromise && !d.raycastView) return;
          if (d.backend !== 'vtk' && typeof d.setBackend === 'function') {
            d.setBackend('vtk');
          } else if (typeof d._renderVtk === 'function') {
            d._renderVtk();
          } else if (typeof d._renderCurrent === 'function') {
            d._renderCurrent();
          }
          const view = d.raycastView;
          if (view && view.renderWindow
              && typeof view.renderWindow.render === 'function') {
            view.renderWindow.render();
          }
          if (typeof d._scheduleVtkResize === 'function') d._scheduleVtkResize();
        }"""
    )


def _ensure_vtk_painted(page, *, log, attempts: int = 2) -> bool:
    """Wait for a painted VTK volume, retrying after WebGL context loss."""
    for attempt in range(1, attempts + 1):
        try:
            page.wait_for_function(
                """() => {
                  const d = window.__cryoTrajVolDisplay;
                  const vols = d && Array.isArray(d.volumes) ? d.volumes : [];
                  return vols.some(function (v) { return v && v.volume_b64; });
                }""",
                timeout=180_000,
            )
            if log and attempt == 1:
                log(f"VTK pre-paint state: {_vtk_debug_state(page)}")
            # Kick one paint, then poll without re-entering _renderVtk (paintGen).
            # Do not create a throwaway WebGL probe context here — loseContext can
            # kill the NVIDIA ANGLE process before vtk.js gets a canvas.
            _force_vtk_render(page)
            # Software WebGL (SwiftShader) often needs well over a minute for the
            # first vtk.js volume raycast; keep this deadline generous.
            deadline = time.monotonic() + 180.0
            while time.monotonic() < deadline:
                if page.evaluate(_vtk_volume_ready_js()):
                    break
                # Re-kick only if the async paint never started a raycast view.
                started = page.evaluate(
                    """() => {
                      const d = window.__cryoTrajVolDisplay;
                      return !!(d && (d.raycastView || d.vtkModulePromise));
                    }"""
                )
                if not started:
                    _force_vtk_render(page)
                time.sleep(0.5)
            else:
                raise TimeoutError("raycastView.volume not ready within 180s")
            page.wait_for_function(_vtk_canvas_painted_js(), timeout=60_000)
            # Give the compositor a beat, then confirm pixels (not just canvas size).
            time.sleep(0.35)
            page.evaluate(
                """() => {
                  const d = window.__cryoTrajVolDisplay;
                  const view = d && d.raycastView;
                  if (view && view.renderWindow
                      && typeof view.renderWindow.render === 'function') {
                    view.renderWindow.render();
                  }
                }"""
            )
            time.sleep(0.2)
            if _vtk_host_has_volume_pixels(page):
                if log:
                    log(f"VTK canvas appears painted (attempt {attempt})")
                return True
            # SwiftShader / some ANGLE paths composite to an empty Playwright
            # screenshot even when vtk.js holds a live volume actor. Prefer the
            # actor+canvas readiness signal over a destructive retry that would
            # bump _vtkPaintGen and lose the context.
            if page.evaluate(_vtk_canvas_painted_js()):
                if log:
                    log(
                        f"VTK canvas appears painted (attempt {attempt}; "
                        "screenshot empty — trusting raycast actor)"
                    )
                return True
            if log:
                log(
                    f"WARNING: VTK canvas present but screenshot is empty "
                    f"(attempt {attempt}); state={_vtk_debug_state(page)}; retrying"
                )
        except Exception as err:
            if log:
                # lose=False: a diagnostic probe must not tear down an in-flight
                # SwiftShader / ANGLE context before the next attempt.
                log(
                    f"WARNING: VTK paint attempt {attempt}/{attempts} failed ({err}); "
                    f"state={_vtk_debug_state(page)}; "
                    f"webgl={_webgl_renderer_info(page, lose=False)}"
                )
            if attempt >= attempts:
                break
        if attempt >= attempts:
            break
        _lose_webgl_contexts(page)
        time.sleep(0.5)
        _select_backend(page, "vtk")
        time.sleep(0.35)
    return False


def _vtk_host_has_volume_pixels(page) -> bool:
    """True when a Playwright clip of the VTK host shows non-background content."""
    try:
        from PIL import Image
        import io
    except ImportError:
        return False
    box = page.locator("#vslice-vtk-container").bounding_box()
    if not box or box["width"] < 32 or box["height"] < 32:
        return False
    png = page.screenshot(
        type="png",
        clip={
            "x": box["x"],
            "y": box["y"],
            "width": box["width"],
            "height": box["height"],
        },
    )
    im = Image.open(io.BytesIO(png)).convert("RGB")
    px = list(im.getdata())
    if not px:
        return False
    # VTK clear colour is ~beige (250,248,244); volume isosurface is mid-gray.
    darkish = sum(1 for r, g, b in px if (r + g + b) < 620)
    return darkish >= max(80, int(0.01 * len(px)))


def _orbit_vtk_volume(
    buf: FrameBuffer,
    page,
    *,
    budget_s: float,
    log,
    start_turn_y_deg: float = VTK_ORBIT_START_TURN_Y_DEG,
    start_turn_x_deg: float = VTK_ORBIT_START_TURN_X_DEG,
    start_turn_z_deg: float = VTK_ORBIT_START_TURN_Z_DEG,
    end_turn_y_deg: float = VTK_ORBIT_END_TURN_Y_DEG,
    end_turn_x_deg: float = VTK_ORBIT_END_TURN_X_DEG,
    end_turn_z_deg: float = VTK_ORBIT_END_TURN_Z_DEG,
) -> None:
    """Slow short VTK reframe: yaw, tip, and a mild Z twist.

    Uses ``applyChimeraxViewTurns`` each frame. Angles stay well below 90° so
    the motion is a calm twist, not a polar flip or multi-spin.
    """
    has_view = page.evaluate(
        """() => {
          const disp = window.__cryoTrajVolDisplay;
          return !!(disp && disp.raycastView
            && typeof disp.raycastView.applyChimeraxViewTurns === 'function'
            && disp.raycastView.volume);
        }"""
    )
    if not has_view:
        log("WARNING: VTK raycast view unavailable for orbit; trying mouse drag")
        box = page.locator("#vslice-vtk-container").bounding_box()
        if not box or box["width"] < 32 or box["height"] < 32:
            buf.sleep_snap(page, max(0.2, budget_s))
            return
        cx = box["x"] + box["width"] * 0.5
        cy = box["y"] + box["height"] * 0.5
        steps = max(8, int(round(budget_s * 1000 / buf.frame_ms)))
        page.mouse.move(cx, cy)
        page.mouse.down()
        for s in range(1, steps + 1):
            t = s / steps
            ease = t * t * (3.0 - 2.0 * t)
            # Mild diagonal drag (fallback when ChimeraX-turn API is unavailable).
            page.mouse.move(cx + 90.0 * ease, cy + 40.0 * ease)
            buf.snap(page)
            time.sleep(buf.frame_ms / 1000.0)
        page.mouse.up()
        buf.sleep_snap(page, 0.25)
        return

    steps = max(8, int(round(budget_s * 1000 / buf.frame_ms)))
    log(
        f"VTK: ChimeraX-turn sweep ({steps} frames, "
        f"start y={start_turn_y_deg:.0f}° x={start_turn_x_deg:.0f}° "
        f"z={start_turn_z_deg:.0f}° → "
        f"end y={end_turn_y_deg:.0f}° x={end_turn_x_deg:.0f}° "
        f"z={end_turn_z_deg:.0f}°)"
    )
    page.evaluate(
        """() => {
          const disp = window.__cryoTrajVolDisplay;
          if (disp && typeof disp._markVtkCameraUserAdjusted === 'function') {
            disp._markVtkCameraUserAdjusted();
          }
        }"""
    )

    def _apply_turns(turn_y: float, turn_x: float, turn_z: float) -> None:
        turns = []
        if abs(turn_y) > 1e-6:
            turns.append({"axis": "y", "degrees": float(turn_y)})
        if abs(turn_x) > 1e-6:
            turns.append({"axis": "x", "degrees": float(turn_x)})
        if abs(turn_z) > 1e-6:
            turns.append({"axis": "z", "degrees": float(turn_z)})
        page.evaluate(
            """(turns) => {
              const d = window.__cryoTrajVolDisplay;
              const view = d && d.raycastView;
              if (!view || typeof view.applyChimeraxViewTurns !== 'function') return;
              if (!turns.length) {
                if (typeof view.resetCamera === 'function') view.resetCamera();
              } else {
                view.applyChimeraxViewTurns(turns);
              }
              if (view.renderWindow && typeof view.renderWindow.render === 'function') {
                view.renderWindow.render();
              }
            }""",
            turns,
        )

    _apply_turns(start_turn_y_deg, start_turn_x_deg, start_turn_z_deg)
    page.evaluate(
        """() => new Promise((resolve) => {
          requestAnimationFrame(() => requestAnimationFrame(resolve));
        })"""
    )
    buf.sleep_snap(page, 0.2)

    d_y = end_turn_y_deg - start_turn_y_deg
    d_x = end_turn_x_deg - start_turn_x_deg
    d_z = end_turn_z_deg - start_turn_z_deg
    for s in range(1, steps + 1):
        t = s / steps
        ease = t * t * (3.0 - 2.0 * t)
        _apply_turns(
            start_turn_y_deg + d_y * ease,
            start_turn_x_deg + d_x * ease,
            start_turn_z_deg + d_z * ease,
        )
        page.evaluate(
            """() => new Promise((resolve) => {
              requestAnimationFrame(() => requestAnimationFrame(resolve));
            })"""
        )
        buf.snap(page)
        time.sleep(buf.frame_ms / 1000.0)
    _apply_turns(end_turn_y_deg, end_turn_x_deg, end_turn_z_deg)
    buf.sleep_snap(page, 0.3)


def _reset_volume_view_default(page, *, log=None) -> None:
    """Reset VTK / ChimeraX view state to the default orient pose."""
    page.evaluate(
        """() => {
          const d = window.__cryoTrajVolDisplay;
          if (d && typeof d.resetInteractiveView === 'function') {
            d.resetInteractiveView();
          } else if (d && d.raycastView
              && typeof d.raycastView.resetCamera === 'function') {
            d.raycastView.resetCamera();
            if (d.raycastView.renderWindow
                && typeof d.raycastView.renderWindow.render === 'function') {
              d.raycastView.renderWindow.render();
            }
          }
          // Clear any ChimeraX turn / matrix state so Generate uses identity.
          if (d && typeof d.onSyncVtkViewToChimerax === 'function') {
            d.onSyncVtkViewToChimerax({ view_matrix: '', view_turns: [] });
          }
          if (d) {
            d.vtkCameraUserAdjusted = false;
            if (typeof d.setSharedViewMatrix === 'function') d.setSharedViewMatrix('');
          }
        }"""
    )
    if log:
        log("Reset volume view to default orient")


def _slider_max(page) -> int:
    return int(
        page.evaluate(
            """() => {
              const el = document.getElementById('traj-vol-volume-slider');
              if (!el) return 0;
              const mx = parseInt(el.max, 10);
              return Number.isFinite(mx) ? mx : 0;
            }"""
        )
        or 0
    )


def _set_slider_index(page, index: int) -> None:
    page.evaluate(
        """(idx) => {
          const el = document.getElementById('traj-vol-volume-slider');
          if (!el) throw new Error('missing #traj-vol-volume-slider');
          el.value = String(idx);
          el.dispatchEvent(new Event('input', { bubbles: true }));
          el.dispatchEvent(new Event('change', { bubbles: true }));
        }""",
        int(index),
    )


def _cycle_volume_slider(
    buf: FrameBuffer,
    page,
    *,
    log,
    label: str,
    sec_per_tick: float = VOLUME_SLIDER_SEC_PER_TICK,
    budget_s: float | None = None,
    wait_ready: bool = True,
    backend: str | None = None,
) -> None:
    """Step the volume slider with a fixed dwell on each tick.

    ``sec_per_tick`` is the GIF-time spent on every volume (same for VTK and
    ChimeraX). When ``wait_ready`` is true, each tick waits for that slot's
    volume / ChimeraX preview before dwelling so unfinished slots are not shown.
    """
    del budget_s  # dwell is per-tick, not derived from a segment budget
    mx = _slider_max(page)
    n_steps = max(1, mx + 1)
    per = max(buf.frame_ms / 1000.0, float(sec_per_tick))
    log(f"{label}: cycle slider 0…{mx} (~{per:.2f}s/step, wait_ready={wait_ready})")
    for i in range(n_steps):
        _set_slider_index(page, i)
        if wait_ready:
            ok = _wait_volume_slot_ready(page, i, backend=backend, timeout_ms=180_000)
            if not ok and log:
                log(f"WARNING: {label} slot {i} not ready; continuing")
            # Give VTK / ChimeraX a beat to paint the newly focused slot.
            page.evaluate(
                """() => new Promise((resolve) => {
                  requestAnimationFrame(() => requestAnimationFrame(resolve));
                })"""
            )
            time.sleep(0.08)
        else:
            time.sleep(0.05)
        buf.sleep_snap(page, per)


def _traj_marker_client_xy(page, index: int) -> tuple[float, float]:
    pt = page.evaluate(
        """(i) => {
          const overlay = document.getElementById('traj-glyph-overlay');
          if (!overlay) return null;
          const markers = overlay.querySelectorAll('.cryo-traj-glyph-marker');
          const m = markers[i];
          if (!m) return null;
          const r = m.getBoundingClientRect();
          return [r.left + r.width / 2, r.top + r.height / 2];
        }""",
        int(index),
    )
    if not pt or len(pt) != 2:
        raise RuntimeError(f"No trajectory glyph marker at index {index}")
    return float(pt[0]), float(pt[1])


def _data_to_client_xy(page, x: float, y: float) -> tuple[float, float]:
    """Map scatter data coordinates to viewport client pixels via Plotly axes."""
    pt = page.evaluate(
        """([x, y]) => {
          const gd = document.getElementById('scatter');
          if (!gd || !gd._fullLayout) return null;
          const xa = gd._fullLayout.xaxis;
          const ya = gd._fullLayout.yaxis;
          if (!xa || !ya || typeof xa.l2p !== 'function' || typeof ya.l2p !== 'function') {
            return null;
          }
          const rect = gd.getBoundingClientRect();
          return [rect.left + xa._offset + xa.l2p(x), rect.top + ya._offset + ya.l2p(y)];
        }""",
        [float(x), float(y)],
    )
    if not pt or len(pt) != 2:
        raise RuntimeError(f"Could not map data ({x}, {y}) to client pixels")
    return float(pt[0]), float(pt[1])


def _scatter_axis_cols(page) -> tuple[str, str]:
    cols = page.evaluate(
        """() => {
          const sx = document.getElementById('sx');
          const sy = document.getElementById('sy');
          return [
            (sx && sx.value) ? sx.value : 'PC1',
            (sy && sy.value) ? sy.value : 'PC2',
          ];
        }"""
    )
    return str(cols[0]), str(cols[1])


def _kmeans_centroid_xy(page, *, k: int = DIRECT_KMEANS_K) -> list[list[float]]:
    """Return ``k`` well-spaced k-means centre positions in current scatter axes.

    Uses ``/api/trajectory_kmeans_centers`` (analyze ``kmeansK`` folder). When more
    than ``k`` centres exist, farthest-point samples ``k`` of them and orders the
    path by angle around their mean so the densified polyline reads cleanly.
    """
    xcol, ycol = _scatter_axis_cols(page)
    payload = page.evaluate(
        """async ([xcol, ycol]) => {
          const r = await fetch('/api/trajectory_kmeans_centers', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ x: xcol, y: ycol, mode: 'direct', n_points: 0 }),
          });
          if (!r.ok) {
            const t = await r.text();
            throw new Error('kmeans centers HTTP ' + r.status + ': ' + t.slice(0, 200));
          }
          return await r.json();
        }""",
        [xcol, ycol],
    )
    if not isinstance(payload, dict) or not payload.get("ok", True):
        raise RuntimeError(f"k-means centres request failed: {payload!r}")
    anchors = payload.get("anchor_indices") or []
    traj_xy = payload.get("traj_xy") or []
    n_anchors = len(anchors)
    if n_anchors < 2:
        raise RuntimeError(f"Need ≥2 k-means centres; got {n_anchors}")
    # With interpolation_points=0, traj_xy should be one row per anchor.
    if len(traj_xy) == n_anchors:
        centers = [[float(p[0]), float(p[1])] for p in traj_xy]
    else:
        # Fallback: stride through a densified path (n_points>0 responses).
        stride = max(1, len(traj_xy) // max(1, n_anchors))
        centers = []
        for i in range(n_anchors):
            row = traj_xy[min(i * stride, len(traj_xy) - 1)]
            centers.append([float(row[0]), float(row[1])])
    if len(centers) <= k:
        chosen = centers
    else:
        # Farthest-point sampling in plot space.
        chosen_idx = [0]
        for _ in range(k - 1):
            best_i, best_d = -1, -1.0
            for i, c in enumerate(centers):
                if i in chosen_idx:
                    continue
                d = min(
                    (c[0] - centers[j][0]) ** 2 + (c[1] - centers[j][1]) ** 2
                    for j in chosen_idx
                )
                if d > best_d:
                    best_d, best_i = d, i
            if best_i < 0:
                break
            chosen_idx.append(best_i)
        chosen = [centers[i] for i in chosen_idx]
    # Angular order around the mean for a non-crossing polyline.
    mx = sum(p[0] for p in chosen) / len(chosen)
    my = sum(p[1] for p in chosen) / len(chosen)
    chosen = sorted(chosen, key=lambda p: math.atan2(p[1] - my, p[0] - mx))
    return chosen


def _wait_marker_count(page, n: int, *, timeout_ms: int = 180_000) -> None:
    page.wait_for_function(
        """(n) => {
          const overlay = document.getElementById('traj-glyph-overlay');
          const markers = overlay
            ? overlay.querySelectorAll('.cryo-traj-glyph-marker').length
            : 0;
          return markers === n;
        }""",
        arg=int(n),
        timeout=timeout_ms,
    )


def _wait_traj_controls_idle(page, *, timeout_ms: int = WAIT_TIMEOUT_MS) -> None:
    """Wait until Points / path mutations are accepted (not volume/path-busy)."""
    page.wait_for_function(
        """() => {
          if (typeof trajectoryAddPointsBusy === 'function') {
            return !trajectoryAddPointsBusy();
          }
          return true;
        }""",
        timeout=timeout_ms,
    )


def _set_n_points_and_wait(page, n: int, *, log=None) -> None:
    """Set ``#n-points`` and wait for glyph markers; retry if a busy gate no-ops."""
    n = int(n)
    _wait_traj_controls_idle(page)
    _set_select_value(page, "n-points", str(n))
    try:
        _wait_marker_count(page, n, timeout_ms=45_000)
        return
    except Exception:
        if log:
            log(
                f"Points={n} did not apply (likely busy gate); "
                f"markers={_marker_count(page)}; retrying"
            )
    _wait_traj_controls_idle(page)
    page.evaluate(
        """(val) => {
          const el = document.getElementById('n-points');
          if (!el) throw new Error('missing #n-points');
          el.value = String(val);
          el.dispatchEvent(new Event('change', { bubbles: true }));
        }""",
        n,
    )
    _wait_marker_count(page, n, timeout_ms=WAIT_TIMEOUT_MS)


def _drag_traj_point(
    buf: FrameBuffer,
    page,
    index: int,
    dx: float,
    dy: float,
    *,
    steps: int = 14,
) -> None:
    x, y = _traj_marker_client_xy(page, index)
    page.mouse.move(x, y)
    page.mouse.down()
    for s in range(1, steps + 1):
        t = s / steps
        page.mouse.move(x + dx * t, y + dy * t)
        buf.snap(page)
        time.sleep(buf.frame_ms / 1000.0)
    page.mouse.up()
    buf.sleep_snap(page, 0.25)


def _move_traj_point_to_data(
    buf: FrameBuffer,
    page,
    index: int,
    data_xy: list[float] | tuple[float, float],
    *,
    steps: int = 16,
) -> None:
    """Drag glyph marker ``index`` to scatter data coordinates ``(x, y)``."""
    x0, y0 = _traj_marker_client_xy(page, index)
    x1, y1 = _data_to_client_xy(page, float(data_xy[0]), float(data_xy[1]))
    page.mouse.move(x0, y0)
    page.mouse.down()
    for s in range(1, steps + 1):
        t = s / steps
        ease = t * t * (3.0 - 2.0 * t)
        page.mouse.move(x0 + (x1 - x0) * ease, y0 + (y1 - y0) * ease)
        buf.snap(page)
        time.sleep(buf.frame_ms / 1000.0)
    page.mouse.up()
    buf.sleep_snap(page, 0.2)


def _marker_count(page) -> int:
    return int(
        page.evaluate(
            """() => {
              const overlay = document.getElementById('traj-glyph-overlay');
              return overlay
                ? overlay.querySelectorAll('.cryo-traj-glyph-marker').length
                : 0;
            }"""
        )
        or 0
    )


def _click_generate_if_enabled(page) -> bool:
    enabled = page.evaluate(
        """() => {
          const btn = document.getElementById('btn-generate-volumes');
          return !!(btn && !btn.hidden && !btn.disabled);
        }"""
    )
    if enabled:
        page.click("#btn-generate-volumes")
    return bool(enabled)


def record_sequence(page, base: str, buf: FrameBuffer, *, log) -> None:
    t0 = time.monotonic()

    with _timed_step(log, "Open trajectory creator", slow_note=True):
        page.goto(base + "/trajectory", wait_until="domcontentloaded", timeout=180_000)
        _freeze_progress_bar_motion(page)
        page.wait_for_function(_scatter_ready_js(), timeout=WAIT_TIMEOUT_MS)
        page.wait_for_function(_pc1_path_ready_js(), timeout=WAIT_TIMEOUT_MS)
        _zoom_scatter_to_cloud(page, pad_frac=0.04, blend=0.5, log=log)
    log(
        f"PC1 path ready with {_marker_count(page)} glyph markers "
        f"(wall {time.monotonic() - t0:.1f}s)"
    )
    # --- 1. Initial PC1 path -------------------------------------------------
    seg_end = SEG_PC1_S
    _hold_until(buf, page, seg_end)

    # --- 2. Hydrate → VTK (before Plotly colour churn / ChimeraX) ------------
    # NVIDIA ANGLE can lose its GPU process after a full multi-volume hydrate;
    # load a single cache blob first, paint VTK, then allow per-tick fetches.
    vtk_budget_start = _buf_wall_s(buf)
    unroute_light = None
    with _timed_step(log, "Hydrate volumes for VTK", slow_note=True):
        hydrated = False
        try:
            unroute_light = _install_light_vtk_hydrate_route(
                page.context, max_volumes=1, log=log
            )
            page.wait_for_function(
                """() => {
                  const slice = document.getElementById('traj-vol-backend-slice');
                  return !!(slice && !slice.disabled);
                }""",
                timeout=30_000,
            )
            page.click("label[for='traj-vol-backend-slice']")
            page.wait_for_function(
                """() => {
                  const d = window.__cryoTrajVolDisplay;
                  const vols = d && Array.isArray(d.volumes) ? d.volumes : [];
                  return vols.some(function (v) { return v && v.volume_b64; });
                }""",
                timeout=180_000,
            )
            hydrated = True
            n_b64 = int(
                page.evaluate(
                    """() => {
                      const d = window.__cryoTrajVolDisplay;
                      const vols = d && Array.isArray(d.volumes) ? d.volumes : [];
                      return vols.filter(function (v) { return v && v.volume_b64; }).length;
                    }"""
                )
                or 0
            )
            log(
                f"Hydrated volume_b64 via slice backend ({n_b64} slot(s), pre-ChimeraX)"
            )
        except Exception as err:
            log(f"Slice pre-hydrate unavailable ({err}); using ChimeraX Generate")
        if not hydrated:
            if unroute_light is not None:
                try:
                    unroute_light()
                except Exception:
                    pass
                unroute_light = None
            _select_backend(page, "chimerax")
            buf.sleep_snap(page, 0.2)
            _click_generate_and_await(buf, page, log=log)
    with _timed_step(log, "Activate VTK backend", slow_note=True):
        _activate_vtk(page, already_hydrated=hydrated)
        if not _ensure_vtk_painted(page, log=log):
            log("WARNING: VTK canvas did not report a painted framebuffer")
        # Allow subsequent slider focus-fetches to load the remaining volumes.
        if unroute_light is not None:
            try:
                unroute_light()
            except Exception:
                pass
            unroute_light = None
        buf.sleep_snap(page, 0.35)
        n_vtk = max(1, _slider_max(page) + 1)
        # Decode any remaining path volumes so every slider tick has volume_b64.
        if _count_ready_volumes(page, backend="vtk") < n_vtk:
            if _click_generate_if_enabled(page):
                log(
                    f"VTK: Generate remaining volumes (have {_count_ready_volumes(page, backend='vtk')}/{n_vtk})"
                )
                _snap_until_ready_capped(
                    buf, page, _volumes_ready_js(), max_clip_ms=RENDER_CLIP_MAX_MS
                )
        # Walk every slot once so interactive VTK hydrates the full path before
        # the recorded scrub (avoids blank unfinished volumes mid-slider).
        for i in range(n_vtk):
            _set_slider_index(page, i)
            _wait_volume_slot_ready(page, i, backend="vtk", timeout_ms=180_000)
        _set_slider_index(page, 0)
        _await_all_volume_slots(
            page, n_vtk, backend="vtk", timeout_ms=WAIT_TIMEOUT_MS, log=log, label="VTK"
        )
        try:
            page.wait_for_function(_ticks_active_js(), timeout=120_000)
        except Exception:
            log("VTK ticks not all active yet; continuing")
    orbit_budget = min(
        SEG_VTK_ORBIT_S,
        max(0.6, SEG_VTK_S - (_buf_wall_s(buf) - vtk_budget_start) - 2.0),
    )
    with _timed_step(log, "Orbit VTK volume", slow_note=False):
        _orbit_vtk_volume(buf, page, budget_s=orbit_budget, log=log)
    with _timed_step(log, "Reset VTK view to default", slow_note=False):
        _reset_volume_view_default(page, log=log)
        buf.sleep_snap(page, 0.25)
    _cycle_volume_slider(buf, page, log=log, label="VTK", backend="vtk")
    seg_end += SEG_VTK_S
    _hold_until(buf, page, seg_end, scale=FINAL_VOLUME_LINGER_SCALE)
    seg_end = _buf_wall_s(buf)

    # --- 3. Colour → znorm (after VTK so Plotly updates cannot steal WebGL) --
    with _timed_step(log, "Colour scatter by znorm", slow_note=False):
        chosen = _choose_znorm(page)
        log(f"Selected colour covariate {chosen!r}")
        _snap_until_ready_capped(
            buf, page, _scatter_ready_js(), max_clip_ms=400, timeout_ms=180_000
        )
        # Colour restyle can restore autorange; keep the blended framing.
        _zoom_scatter_to_cloud(page, pad_frac=0.04, blend=0.5, log=log)
    seg_end += SEG_ZNORM_S
    _hold_until(buf, page, seg_end)

    # --- 4. ChimeraX mode + rendering ----------------------------------------
    # Decode temporarily forces the volume display onto VTK even when ChimeraX
    # is selected; veil that intermediate so the GIF does not show a second VTK
    # paint after the orbit / slider scrub.
    n_cx = max(1, _slider_max(page) + 1)
    with _timed_step(log, "Switch to ChimeraX backend", slow_note=True):
        # Keep the same default orient used for the VTK slider scrub.
        _reset_volume_view_default(page, log=log)
        _ensure_chimerax_images(
            buf, page, n_cx, log=log, label="ChimeraX", timeout_ms=WAIT_TIMEOUT_MS
        )
        buf.sleep_snap(page, 0.25)
    _cycle_volume_slider(buf, page, log=log, label="ChimeraX", backend="chimerax")
    seg_end += SEG_CHIMERAX_S
    _hold_until(buf, page, seg_end, scale=FINAL_VOLUME_LINGER_SCALE)
    seg_end = _buf_wall_s(buf)

    # --- 5. Direct-trace: k centroids → densify → ChimeraX ------------------
    with _timed_step(log, "Switch to direct-trace mode", slow_note=False):
        _wait_traj_controls_idle(page)
        page.click("label[for='traj-mode-direct']")
        page.wait_for_function(
            """() => {
              const el = document.getElementById('traj-mode-direct');
              return !!(el && el.checked);
            }""",
            timeout=60_000,
        )
        page.wait_for_function(
            """() => {
              const overlay = document.getElementById('traj-glyph-overlay');
              return !!(overlay
                && overlay.querySelectorAll('.cryo-traj-glyph-marker').length >= 2);
            }""",
            timeout=120_000,
        )
        _wait_traj_controls_idle(page)
        buf.sleep_snap(page, 0.35)

    with _timed_step(
        log, f"Set Trajectory controls Points={DIRECT_KMEANS_K}", slow_note=False
    ):
        _set_n_points_and_wait(page, DIRECT_KMEANS_K, log=log)
        buf.sleep_snap(page, 0.4)

    with _timed_step(
        log, f"Move path onto {DIRECT_KMEANS_K} k-means centroids", slow_note=True
    ):
        centers = _kmeans_centroid_xy(page, k=DIRECT_KMEANS_K)
        if len(centers) != DIRECT_KMEANS_K:
            raise RuntimeError(
                f"Expected {DIRECT_KMEANS_K} k-means centres; got {len(centers)}"
            )
        log(
            "Centroid targets: "
            + ", ".join(f"({c[0]:.2f},{c[1]:.2f})" for c in centers)
        )
        for i, xy in enumerate(centers):
            _move_traj_point_to_data(buf, page, i, xy, steps=8)
        buf.sleep_snap(page, 0.35)

    with _timed_step(
        log, f"Densify path to Points={DIRECT_DENSIFY_N}", slow_note=False
    ):
        _set_n_points_and_wait(page, DIRECT_DENSIFY_N, log=log)
        buf.sleep_snap(page, 0.45)

    # Decode forces trajVolDisplay onto VTK; ChimeraX PNGs still land in
    # chimeraxImages. Veil that intermediate and flip back to ChimeraX once the
    # densified frames exist so the GIF does not linger on a second VTK paint.
    with _timed_step(log, "Render ChimeraX after densify", slow_note=True):
        _ensure_chimerax_images(
            buf,
            page,
            DIRECT_DENSIFY_N,
            log=log,
            label="Direct ChimeraX",
            timeout_ms=WAIT_TIMEOUT_MS,
        )
        buf.sleep_snap(page, 0.35)

    _cycle_volume_slider(
        buf, page, log=log, label="Direct ChimeraX", backend="chimerax"
    )
    seg_end += SEG_DIRECT_S
    _hold_until(buf, page, seg_end, scale=FINAL_VOLUME_LINGER_SCALE)
    seg_end = _buf_wall_s(buf)
    log(
        f"Sequence complete: {_buf_wall_s(buf):.2f}s GIF wall, "
        f"{time.monotonic() - t0:.1f}s wall-clock"
    )


def _dashboard_command(args: argparse.Namespace, outdir: Path, port: int) -> list[str]:
    dash = [
        "cryodrgn",
        "dashboard",
        str(outdir),
        "--no-browser",
        "--host",
        "127.0.0.1",
        "--port",
        str(port),
        "--traj",
        "--cpus",
        str(args.cpus),
    ]
    if args.conda_prefix:
        return [
            "conda",
            "run",
            "-p",
            str(Path(args.conda_prefix).expanduser()),
            "--no-capture-output",
            *dash,
        ]
    if args.conda_env:
        return ["conda", "run", "-n", args.conda_env, "--no-capture-output", *dash]
    cryo = shutil.which("cryodrgn")
    if not cryo:
        raise RuntimeError(
            "`cryodrgn` is not on PATH; pass --conda-prefix or --conda-env."
        )
    return [cryo, *dash[1:]]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", type=Path)
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=DEMO_ANIMATIONS_RECORDERS_DIR / "trajectory_creator_demo.gif",
    )
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument("--conda-env", type=str, default=None)
    parser.add_argument("--conda-prefix", type=Path, default=None)
    parser.add_argument("--cpus", type=int, default=10)
    parser.add_argument(
        "--frame-ms",
        type=int,
        default=FRAME_MS_DEFAULT,
        help=(
            f"Frame duration in ms (default {FRAME_MS_DEFAULT}, "
            f"about {1000 / FRAME_MS_DEFAULT:.1f} fps)"
        ),
    )
    parser.add_argument(
        "--max-width",
        type=int,
        default=GIF_MAX_WIDTH_DEFAULT,
        help=f"GIF encode max width (default {GIF_MAX_WIDTH_DEFAULT})",
    )
    parser.add_argument(
        "--headed",
        action="store_true",
        help="Run Chromium headed (use with xvfb-run if no DISPLAY)",
    )
    parser.add_argument(
        "--swiftshader",
        action="store_true",
        help="Use software WebGL (SwiftShader) instead of NVIDIA ANGLE GL/EGL",
    )
    parser.add_argument("-q", "--quiet", action="store_true")
    args = parser.parse_args()

    outdir = args.outdir.resolve()
    if not outdir.is_dir():
        print(f"error: not a directory: {outdir}", file=sys.stderr)
        return 1

    try:
        import PIL.Image  # noqa: F401
        from playwright.sync_api import sync_playwright
    except ImportError as err:
        print(f"error: install dashboard-gif dependencies ({err})", file=sys.stderr)
        return 1

    log = _make_logger(quiet=args.quiet)
    port = args.port or _free_port()
    base = f"http://127.0.0.1:{port}"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO_ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    # Force ChimeraX onto Mesa software GLX. Default GLX can bind NVIDIA and then
    # fail with BadAlloc after Chromium has opened an ANGLE GL/EGL WebGL context.
    env["__GLX_VENDOR_LIBRARY_NAME"] = "mesa"
    env["LIBGL_ALWAYS_SOFTWARE"] = "1"
    env["GALLIUM_DRIVER"] = "llvmpipe"
    env.setdefault("CRYODRGN_CHIMERAX_XVFB", "1")
    browser_env = os.environ.copy()
    # Mesa software GL is for the ChimeraX server process only. Leaving those
    # vars in Chromium breaks NVIDIA ANGLE GL/EGL WebGL (VTK paints blank).
    for _gl_key in (
        "__GLX_VENDOR_LIBRARY_NAME",
        "LIBGL_ALWAYS_SOFTWARE",
        "GALLIUM_DRIVER",
    ):
        browser_env.pop(_gl_key, None)
    if not args.swiftshader:
        browser_env.setdefault(
            "__EGL_VENDOR_LIBRARY_FILENAMES",
            "/usr/share/glvnd/egl_vendor.d/10_nvidia.json",
        )
        browser_env.setdefault("NVIDIA_DRIVER_CAPABILITIES", "all")
        # Do not set __GLX_VENDOR_LIBRARY_NAME=nvidia here: it can leave WebGL
        # dead after /trajectory loads even though a blank-page probe succeeds.
    else:
        browser_env = None
    if not env.get("CHIMERAX_PATH"):
        default_cx = "/projects/MOLBIO/local/src/chimerax-1.9/usr/bin/chimerax"
        if Path(default_cx).is_file():
            env["CHIMERAX_PATH"] = default_cx
            log(f"Using CHIMERAX_PATH={default_cx}")

    try:
        cmd = _dashboard_command(args, outdir, port)
    except RuntimeError as err:
        print(f"error: {err}", file=sys.stderr)
        return 1

    log(f"Dashboard command: {' '.join(cmd)}")
    dash_log = Path("/tmp/cryodrgn_traj_gif_dashboard.log")
    dash_log_fh = open(dash_log, "w", encoding="utf-8")
    proc = subprocess.Popen(
        cmd,
        env=env,
        cwd=str(REPO_ROOT),
        stdout=dash_log_fh,
        stderr=subprocess.STDOUT,
        text=True,
    )
    try:
        log(f"Starting dashboard subprocess on port {port} (outdir={outdir})")
        log(f"Dashboard log: {dash_log}")
        _wait_http(base + "/trajectory", timeout_s=300.0, log=log)
        plotly_cached = _ensure_plotly_cached(log)
        buf = FrameBuffer(args.frame_ms)
        ph, launch_args = _trajectory_chromium_launch_options(
            headed=args.headed, swiftshader=args.swiftshader
        )
        with sync_playwright() as p:
            with _timed_step(log, f"Launch Chromium (headless={ph})", slow_note=False):
                launch_kwargs = {"headless": ph, "args": launch_args}
                if browser_env is not None:
                    launch_kwargs["env"] = browser_env
                browser = p.chromium.launch(**launch_kwargs)
            context = browser.new_context(
                # Full UI viewport so side panels / volume square are not clipped;
                # GIF encode still downscales via --max-width.
                viewport={"width": 1440, "height": 900},
                device_scale_factor=1,
            )
            if plotly_cached is not None:
                plotly_bytes = plotly_cached.read_bytes()

                def _fulfill_plotly(route):
                    route.fulfill(
                        status=200,
                        headers={
                            "content-type": "application/javascript; charset=utf-8"
                        },
                        body=plotly_bytes,
                    )

                context.route(PLOTLY_CDN_ROUTE_GLOB, _fulfill_plotly)
            _install_explorer_scatter_recording_routes(context, log=log)
            page = context.new_page()
            with _timed_step(log, "Full trajectory creator capture", slow_note=True):
                record_sequence(page, base, buf, log=log)
            # Probe after capture only for diagnostics (blank-page probes can
            # leave a half-dead GPU context before the trajectory page loads).
            try:
                log(f"WebGL after capture: {_webgl_renderer_info(page)}")
            except Exception:
                pass
            browser.close()

        nf, wall_s, used_w, size_b = _save_buf_to_gif_under_budget(
            buf,
            args.output,
            frame_ms=args.frame_ms,
            primary_max_width=args.max_width,
            max_wall_ms=MAX_GIF_WALL_MS,
            max_bytes=MAX_GIF_BYTES,
            log=log,
        )
        mean_ms = (wall_s * 1000.0 / nf) if nf else 0.0
        print(
            f"Wrote {args.output} ({nf} frames, approx {wall_s:.2f} s playback, "
            f"mean {mean_ms:.1f} ms/frame, width={used_w}, "
            f"{size_b / (1024 * 1024):.2f} MiB)"
        )
        return 0
    except Exception as err:
        print(f"error: {err}", file=sys.stderr)
        try:
            dash_log_fh.flush()
            tail = dash_log.read_text(encoding="utf-8", errors="replace")[-4000:]
            if tail.strip():
                print(tail, file=sys.stderr)
        except Exception:
            pass
        return 1
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=20)
        except subprocess.TimeoutExpired:
            proc.kill()
        try:
            dash_log_fh.close()
        except Exception:
            pass


if __name__ == "__main__":
    raise SystemExit(main())
