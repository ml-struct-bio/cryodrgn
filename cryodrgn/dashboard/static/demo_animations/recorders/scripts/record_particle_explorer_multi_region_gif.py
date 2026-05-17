#!/usr/bin/env python3
"""Headless GIF: particle explorer — high‑K k-means on the scatter, **hull-shaped** lassos, then image grid.

Fits **k-means with ``k_fit``** (default **37**) on the same Plotly trace points shown in the explorer,
then **lassos ``pick_n``** separated clusters using **irregular polygons** from each cluster's
**2‑D convex hull** (slightly expanded), not axis-aligned boxes. Lasso cursor motion is **~70% slower**
than the baseline recorder speed so selections read clearly on the GIF.

Then **image cache** + **montage** as before.

Same Playwright / Pillow stack; default output ``demo_animations/recorders/particle_explorer_multi_region_demo.gif``.
Default encoded length is longer (**32 s**) to cover cache + grid work.

Example::

    PYTHONPATH=/path/to/cryodrgn_beta \\
      conda run -p …/cdrgn_beta --no-capture-output \\
        python -m cryodrgn.dashboard.static.demo_animations.recorders.scripts.record_particle_explorer_multi_region_gif \\
        /path/to/train-vae_out/ \\
        --k-means-fit 37 --lasso-regions 3 \\
        --duration 36
"""

from __future__ import annotations

import argparse
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
    _chromium_launch_options,
    _ensure_plotly_cached,
    _free_port,
    _montage_cache_ready_js,
    _scatter_plot_box,
    _timed_step,
    _wait_http,
    _wait_scatter_ready,
)
from cryodrgn.dashboard.static.demo_animations.recorders.scripts.record_particle_explorer_actions_gif import (
    CaptionedExplorerBuffer,
    _densify_lasso_rel,
    _snap_montage_load,
    _wait_preload_overlay_clear,
    extend_caption_segments_min_duration,
    normalize_frames_fixed_fps,
    save_pngs_duration_gif,
    _stamp_step,
)

DEFAULT_CONDA_PREFIX = Path("/projects/CRYOEM/zhonglab/mg2332/conda_envs/cdrgn_beta")
DEFAULT_OUTPUT = (
    DEMO_ANIMATIONS_RECORDERS_DIR / "particle_explorer_multi_region_demo.gif"
)
_LOG_PREFIX = "[record-particle-explorer-multi-region-gif]"


def _make_mr_logger(quiet: bool):
    def log(msg: str) -> None:
        if not quiet:
            print(f"{_LOG_PREFIX} {msg}", flush=True)

    return log


# Plotly selection handler debounce (~200 ms) plus margin for chip overlay sync.
# Kept tight so longer montage/grid holds do not balloon total raw capture.
_POST_LASSO_SETTLE_S = 0.24

# Cursor moves along each lasso: 70% slower → speed × (1 − 0.70); stretch delays by 1/0.3.
LASSO_DRAW_SLOWDOWN = 1.0 / (1.0 - 0.70)

# Injected: k-means with ``kFit`` on trace (x,y), pick ``pickN`` clusters → irregular lasso rings
# (convex hull of cluster points in data space, lightly expanded, subsampled for Playwright).
_KMEANS_LASSO_PATHS_JS = """
(opts) => {
  opts = opts || {};
  const kFit = Math.max(4, Math.min(64, Number(opts.kFit) || 37));
  const pickN = Math.max(1, Math.min(8, Number(opts.pickN) || 3));
  const gd = document.getElementById('scatter');
  if (!gd || !gd.data || !gd.data[0]) return { ok: false, reason: 'no trace' };
  const xs = gd.data[0].x;
  const ys = gd.data[0].y;
  const n = Math.min(xs.length, ys.length);
  const pts = [];
  for (let i = 0; i < n; i++) {
    const x = +xs[i];
    const y = +ys[i];
    if (isFinite(x) && isFinite(y)) pts.push([x, y]);
  }
  const minPts = Math.max(3, Math.floor(pts.length / (kFit * 4)));
  if (pts.length < kFit * minPts) return { ok: false, reason: 'too few points for kFit' };

  function dist2(a, b) {
    const dx = a[0] - b[0];
    const dy = a[1] - b[1];
    return dx * dx + dy * dy;
  }

  const centers = [];
  for (let j = 0; j < kFit; j++) {
    const idx = Math.min(pts.length - 1, Math.floor(((j + 0.5) / kFit) * pts.length));
    centers.push(pts[idx].slice());
  }

  const labels = new Int32Array(pts.length);
  for (let it = 0; it < 28; it++) {
    for (let i = 0; i < pts.length; i++) {
      let bi = 0;
      let bd = dist2(pts[i], centers[0]);
      for (let j = 1; j < kFit; j++) {
        const d = dist2(pts[i], centers[j]);
        if (d < bd) {
          bd = d;
          bi = j;
        }
      }
      labels[i] = bi;
    }
    const sums = [];
    for (let j = 0; j < kFit; j++) sums.push([0, 0, 0]);
    for (let i = 0; i < pts.length; i++) {
      const b = labels[i];
      sums[b][0] += pts[i][0];
      sums[b][1] += pts[i][1];
      sums[b][2]++;
    }
    let moved = false;
    for (let j = 0; j < kFit; j++) {
      const cnt = sums[j][2];
      const nx = cnt ? sums[j][0] / cnt : centers[j][0];
      const ny = cnt ? sums[j][1] / cnt : centers[j][1];
      if (
        Math.abs(nx - centers[j][0]) > 1e-12 * (1 + Math.abs(nx)) ||
        Math.abs(ny - centers[j][1]) > 1e-12 * (1 + Math.abs(ny))
      ) {
        moved = true;
      }
      centers[j] = [nx, ny];
    }
    if (!moved && it > 3) break;
  }

  function convexHull2D(cps) {
    if (cps.length === 0) return [];
    if (cps.length === 1) return [cps[0].slice()];
    if (cps.length === 2) return [cps[0].slice(), cps[1].slice()];
    const sorted = cps.slice().sort(function (a, b) {
      return a[0] === b[0] ? a[1] - b[1] : a[0] - b[0];
    });
    function cross(o, a, b) {
      return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    }
    const lower = [];
    for (let i = 0; i < sorted.length; i++) {
      while (
        lower.length >= 2 &&
        cross(lower[lower.length - 2], lower[lower.length - 1], sorted[i]) <= 0
      ) {
        lower.pop();
      }
      lower.push(sorted[i]);
    }
    const upper = [];
    for (let i = sorted.length - 1; i >= 0; i--) {
      while (
        upper.length >= 2 &&
        cross(upper[upper.length - 2], upper[upper.length - 1], sorted[i]) <= 0
      ) {
        upper.pop();
      }
      upper.push(sorted[i]);
    }
    upper.pop();
    lower.pop();
    return lower.concat(upper);
  }

  function scaleHullFromCentroid(h, scale) {
    if (!h.length) return h;
    let sx = 0;
    let sy = 0;
    for (let i = 0; i < h.length; i++) {
      sx += h[i][0];
      sy += h[i][1];
    }
    const cx = sx / h.length;
    const cy = sy / h.length;
    return h.map(function (p) {
      return [cx + (p[0] - cx) * scale, cy + (p[1] - cy) * scale];
    });
  }

  function subsampleHullRing(h, maxV) {
    if (h.length <= maxV) return h.slice();
    const out = [];
    const L = h.length;
    for (let i = 0; i < maxV; i++) {
      const idx = Math.min(L - 1, Math.floor((i * L) / maxV));
      out.push(h[idx].slice());
    }
    return out;
  }

  const clusters = [];
  for (let j = 0; j < kFit; j++) {
    const cps = [];
    for (let i = 0; i < pts.length; i++) {
      if (labels[i] === j) cps.push(pts[i].slice());
    }
    if (cps.length < minPts) continue;
    let sx = 0;
    let sy = 0;
    for (let t = 0; t < cps.length; t++) {
      sx += cps[t][0];
      sy += cps[t][1];
    }
    const cx = sx / cps.length;
    const cy = sy / cps.length;
    clusters.push({ id: j, pts: cps, cx: cx, cy: cy, n: cps.length });
  }
  if (clusters.length < 1) return { ok: false, reason: 'no nonempty clusters' };

  let cminx = Infinity;
  let cmaxx = -Infinity;
  let cminy = Infinity;
  let cmaxy = -Infinity;
  for (let c = 0; c < clusters.length; c++) {
    const cl = clusters[c];
    if (cl.cx < cminx) cminx = cl.cx;
    if (cl.cx > cmaxx) cmaxx = cl.cx;
    if (cl.cy < cminy) cminy = cl.cy;
    if (cl.cy > cmaxy) cmaxy = cl.cy;
  }
  const spanx = Math.max(1e-9, cmaxx - cminx);
  const spany = Math.max(1e-9, cmaxy - cminy);
  for (let c = 0; c < clusters.length; c++) {
    const cl = clusters[c];
    cl.ncx = (cl.cx - cminx) / spanx;
    cl.ncy = (cl.cy - cminy) / spany;
  }

  const sorted = clusters.slice().sort(function (a, b) {
    return b.n - a.n;
  });
  const picked = [];
  const minSep2 = 0.07 * 0.07;
  for (let i = 0; i < sorted.length && picked.length < pickN; i++) {
    const c = sorted[i];
    let ok = true;
    for (let p = 0; p < picked.length; p++) {
      const dx = c.ncx - picked[p].ncx;
      const dy = c.ncy - picked[p].ncy;
      if (dx * dx + dy * dy < minSep2) {
        ok = false;
        break;
      }
    }
    if (ok) picked.push(c);
  }
  for (let i = 0; i < sorted.length && picked.length < pickN; i++) {
    const c = sorted[i];
    let dup = false;
    for (let p = 0; p < picked.length; p++) {
      if (picked[p].id === c.id) {
        dup = true;
        break;
      }
    }
    if (!dup) picked.push(c);
  }
  picked.sort(function (a, b) {
    return a.cx - b.cx;
  });

  const fl = gd._fullLayout;
  const xa = fl.xaxis;
  const ya = fl.yaxis;
  const gw = fl.width || gd.clientWidth || 720;
  const gh = fl.height || gd.clientHeight || 480;

  function dataToFrac(xd, yd) {
    if (xa && ya && typeof xa.l2p === 'function' && typeof ya.l2p === 'function') {
      const px = (xa._offset || 0) + xa.l2p(xd);
      const py = (ya._offset || 0) + ya.l2p(yd);
      return [px / gw, py / gh];
    }
    const xr = xa.range;
    const yr = ya.range;
    const dx = xr[1] - xr[0];
    const dy = yr[1] - yr[0];
    const ml = (fl.margin && fl.margin.l) || 0;
    const mr = (fl.margin && fl.margin.r) || 0;
    const mt = (fl.margin && fl.margin.t) || 0;
    const mb = (fl.margin && fl.margin.b) || 0;
    const plotW = Math.max(1e-6, gw - ml - mr);
    const plotH = Math.max(1e-6, gh - mt - mb);
    const fx = (ml + (plotW * (xd - xr[0])) / dx) / gw;
    const fy = (mt + (plotH * (yr[1] - yd)) / dy) / gh;
    return [fx, fy];
  }

  function clamp01(t) {
    return Math.max(0.03, Math.min(0.97, t));
  }

  const maxHullVerts = 20;
  const expandScale = 1.065;
  const paths = [];
  const sizes = [];
  const ids = [];
  for (let c = 0; c < picked.length; c++) {
    const meta = picked[c];
    let h = convexHull2D(meta.pts);
    if (h.length < 3) {
      let xmin = Infinity;
      let xmax = -Infinity;
      let ymin = Infinity;
      let ymax = -Infinity;
      for (let t = 0; t < meta.pts.length; t++) {
        const x = meta.pts[t][0];
        const y = meta.pts[t][1];
        if (x < xmin) xmin = x;
        if (x > xmax) xmax = x;
        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
      }
      const pwx = Math.max(1e-9, xmax - xmin);
      const phy = Math.max(1e-9, ymax - ymin);
      const px = 0.05 * pwx;
      const py = 0.05 * phy;
      h = [
        [xmin - px, ymin - py],
        [xmax + px, ymin - py],
        [xmax + px, ymax + py],
        [xmin - px, ymax + py]
      ];
    }
    h = scaleHullFromCentroid(h, expandScale);
    h = subsampleHullRing(h, maxHullVerts);
    h.push(h[0].slice());
    const ring = h.map(function (p) {
      const q = dataToFrac(p[0], p[1]);
      return [clamp01(q[0]), clamp01(q[1])];
    });
    paths.push(ring);
    sizes.push(meta.n);
    ids.push(meta.id);
  }
  return { ok: true, paths: paths, sizes: sizes, ids: ids, kFit: kFit, pickN: pickN };
}
"""


def _kmeans_lasso_paths_from_page(
    page,
    *,
    k_fit: int,
    pick_n: int,
    log,
) -> list[list[tuple[float, float]]] | None:
    """Return ``pick_n`` closed polygons (fractions of ``#scatter``), or ``None`` for fallback."""
    raw = page.evaluate(
        _KMEANS_LASSO_PATHS_JS,
        {"kFit": int(k_fit), "pickN": int(pick_n)},
    )
    if not isinstance(raw, dict) or not raw.get("ok"):
        reason = raw.get("reason") if isinstance(raw, dict) else "bad payload"
        log(f"k-means lasso path compute failed ({reason}); using fallback polygons")
        return None
    paths_raw = raw.get("paths")
    if not isinstance(paths_raw, list) or len(paths_raw) < 1:
        log("k-means returned no paths; using fallback polygons")
        return None
    out: list[list[tuple[float, float]]] = []
    for poly in paths_raw:
        if not isinstance(poly, list) or len(poly) < 3:
            continue
        ring: list[tuple[float, float]] = []
        for pt in poly:
            if (
                isinstance(pt, (list, tuple))
                and len(pt) == 2
                and isinstance(pt[0], (int, float))
                and isinstance(pt[1], (int, float))
            ):
                ring.append((float(pt[0]), float(pt[1])))
        if len(ring) >= 3:
            out.append(ring)
    if len(out) < 1:
        log("parsed no k-means polygons; using fallback")
        return None
    sizes = raw.get("sizes")
    ids = raw.get("ids")
    if isinstance(sizes, list) and isinstance(ids, list) and len(sizes) == len(ids):
        log(
            f"k-means fit K={raw.get('kFit', '?')}: lasso {len(out)} region(s); "
            "cluster ids "
            + ", ".join(str(int(i)) for i in ids if isinstance(i, (int, float)))
            + " — sizes "
            + ", ".join(str(int(s)) for s in sizes if isinstance(s, (int, float)))
        )
    return out


def _fallback_lasso_demo_paths() -> list[list[tuple[float, float]]]:
    """Hand-tuned relative polys when k-means / axis mapping is unavailable."""
    return [
        [
            (0.32, 0.48),
            (0.44, 0.38),
            (0.58, 0.42),
            (0.64, 0.55),
            (0.52, 0.64),
            (0.38, 0.58),
            (0.30, 0.50),
        ],
        [
            (0.66, 0.18),
            (0.80, 0.16),
            (0.88, 0.28),
            (0.82, 0.38),
            (0.68, 0.34),
            (0.60, 0.24),
        ],
        [
            (0.14, 0.72),
            (0.26, 0.68),
            (0.28, 0.82),
            (0.16, 0.86),
        ],
    ]


def _lasso_polygon(
    page,
    box: dict,
    rel: list[tuple[float, float]],
    buf: CaptionedExplorerBuffer,
    *,
    frame_ms: int,
    subdivisions_between: int = 8,
    intro_snaps: int = 0,
    slowdown: float = LASSO_DRAW_SLOWDOWN,
) -> None:
    """Trace one closed lasso in plot fraction coordinates; caller sets dragmode lasso.

    ``slowdown`` > 1 stretches pauses between moves and adds intermediate waypoints
    (default ``LASSO_DRAW_SLOWDOWN`` ≈ 3.33 for **70% slower** tracing than baseline).
    """
    w, h = box["width"], box["height"]
    sub_n = max(2, int(round(float(subdivisions_between) * float(slowdown))))
    rel_dense = _densify_lasso_rel(rel, sub_n)
    step_s = max(frame_ms / 1000.0 * 2.2, 0.048) * float(slowdown)
    for _ in range(max(0, intro_snaps)):
        buf.snap(page)
        time.sleep(frame_ms / 1000.0)
    page.mouse.move(box["x"] + w * rel_dense[0][0], box["y"] + h * rel_dense[0][1])
    page.mouse.down()
    buf.snap(page)
    time.sleep(step_s)
    for i in range(1, len(rel_dense)):
        page.mouse.move(box["x"] + w * rel_dense[i][0], box["y"] + h * rel_dense[i][1])
        time.sleep(step_s)
        buf.snap(page)
    page.mouse.up()
    buf.snap(page)
    time.sleep(_POST_LASSO_SETTLE_S)


# Montage / grid dwells: first grid, menu, mid re-layout, final grid were scaled up earlier; this
# pass multiplies those holds by **1.3** again. Intro, lasso pauses, post-lasso settle, union, cache
# holds, montage poll budget, and tail are tightened to limit raw capture growth before normalization.
def _image_cache_and_montage_clip(
    buf: CaptionedExplorerBuffer,
    page,
    *,
    frame_ms: int,
    log,
) -> None:
    """Expand cache, preload selection, open montage — mirrors ``record_particle_explorer_actions_gif``."""
    _stamp_step(
        buf,
        "Image cache — widen the preload window so more particle thumbnails stay in RAM.",
    )
    with _timed_step(log, "Build image cache (expand window)", slow_note=True):
        expand = page.locator("#btn-expand-cache")
        if expand.count() and expand.is_enabled():
            expand.click()
            try:
                _wait_preload_overlay_clear(page, 360_000)
            except Exception:
                pass
        for _ in range(max(1, int(round(0.10 * 1000 / frame_ms)))):
            buf.snap(page)
            time.sleep(frame_ms / 1000.0)

    _stamp_step(
        buf,
        "Selection preload — add uncached particles from the multi-region pick into the cache.",
    )
    with _timed_step(log, "Cache selection thumbnails (uncached)", slow_note=True):
        sel_uncached = page.locator("#btn-cache-selection-uncached")
        if sel_uncached.count() and sel_uncached.is_enabled():
            sel_uncached.click()
            try:
                _wait_preload_overlay_clear(page, 360_000)
            except Exception:
                pass
        for _ in range(max(1, int(round(0.09 * 1000 / frame_ms)))):
            buf.snap(page)
            time.sleep(frame_ms / 1000.0)

    with _timed_step(log, "Open image grid from cache", slow_note=True):
        view_btn = page.locator("#btn-view-images")
        if view_btn.count() and view_btn.is_enabled():
            _stamp_step(
                buf,
                "Image montage — grid of thumbnails for the combined k-means selection.",
            )
            view_btn.click()
            try:
                page.wait_for_function(_montage_cache_ready_js(), timeout=180_000)
            except Exception:
                pass
            _snap_montage_load(buf, page, frame_ms=frame_ms, max_polls=110)
            # Grid dwells vs baseline actions recorder: earlier ×2 ×1.6, then ×1.3 again (this pass).
            for _ in range(max(1, int(round(3.536 * 1000 / frame_ms)))):
                buf.snap(page)
                time.sleep(frame_ms / 1000.0)
            _stamp_step(
                buf,
                "Grid layout — adjust rows/columns in the image grid menu.",
            )
            gs = page.locator("#grid-size")
            toggle = page.locator("#image-grid-menu-toggle")
            if toggle.count():
                aria = toggle.get_attribute("aria-expanded") or "false"
                if aria.strip().lower() == "false":
                    toggle.click()
                    for _ in range(max(1, int(round(1.456 * 1000 / frame_ms)))):
                        buf.snap(page)
                        time.sleep(frame_ms / 1000.0)
            if gs.count() and gs.is_enabled():
                try:
                    cur = gs.input_value() or "3"
                except Exception:
                    cur = "3"
                pick = {"2": "5", "3": "5", "4": "6", "5": "3", "6": "4"}.get(cur, "5")
                try:
                    gs.select_option(value=pick)
                except Exception:
                    try:
                        gs.select_option(label=f"{pick} × {pick}")
                    except Exception:
                        pass
                try:
                    page.wait_for_function(_montage_cache_ready_js(), timeout=120_000)
                except Exception:
                    pass
                for _ in range(max(1, int(round(1.872 * 1000 / frame_ms)))):
                    buf.snap(page)
                    time.sleep(frame_ms / 1000.0)
                _snap_montage_load(buf, page, frame_ms=frame_ms, max_polls=42)
                for _ in range(max(1, int(round(4.16 * 1000 / frame_ms)))):
                    buf.snap(page)
                    time.sleep(frame_ms / 1000.0)
        else:
            log(
                "#btn-view-images stayed disabled (build cache first); "
                "montage segment skipped — selection clip still captured"
            )
            for _ in range(max(1, int(round(0.5 * 1000 / frame_ms)))):
                buf.snap(page)
                time.sleep(frame_ms / 1000.0)


def record_multi_region_sequence(
    page,
    base: str,
    buf: CaptionedExplorerBuffer,
    *,
    frame_ms: int,
    log,
    k_fit: int = 37,
    pick_n: int = 3,
) -> None:
    page.goto(base + "/explorer", wait_until="domcontentloaded", timeout=180_000)
    with _timed_step(log, "Wait for scatter render", slow_note=True):
        _wait_scatter_ready(page, 300_000)

    _stamp_step(
        buf,
        f"Particle explorer — k-means with K={k_fit} on the visible scatter; "
        f"lasso {pick_n} separated clusters, then image cache + montage.",
    )
    for _ in range(max(1, int(round(0.06 * 1000 / frame_ms)))):
        buf.snap(page)
        time.sleep(frame_ms / 1000.0)

    box = _scatter_plot_box(page)
    if not box:
        raise RuntimeError("could not find #scatter bounding box")

    page.evaluate(
        """() => {
          const gd = document.getElementById('scatter');
          if (!gd || !window.Plotly) return false;
          Plotly.relayout(gd, { dragmode: 'lasso' });
          return true;
        }"""
    )

    paths = _kmeans_lasso_paths_from_page(page, k_fit=k_fit, pick_n=pick_n, log=log)
    if paths is None:
        paths = _fallback_lasso_demo_paths()

    captions = [
        f"Lasso {i + 1} — convex hull of one of {pick_n} k-means groups (K={k_fit}), traced slowly."
        for i in range(len(paths))
    ]

    for idx, rel in enumerate(paths):
        box = _scatter_plot_box(page) or box
        _stamp_step(buf, captions[idx])
        _lasso_polygon(page, box, rel, buf, frame_ms=frame_ms, intro_snaps=0)
        hold = 0.24 if idx < len(paths) - 1 else 0.12
        for _ in range(max(1, int(round(hold * 1000 / frame_ms)))):
            buf.snap(page)
            time.sleep(frame_ms / 1000.0)

    _stamp_step(
        buf,
        "Union selection — region chips for each lasso; combined particles drive cache + grid below.",
    )
    for _ in range(max(1, int(round(0.12 * 1000 / frame_ms)))):
        buf.snap(page)
        time.sleep(frame_ms / 1000.0)

    _image_cache_and_montage_clip(buf, page, frame_ms=frame_ms, log=log)

    _stamp_step(buf, "Montage ready — thumbnails from the multi-cluster selection.")
    for _ in range(max(1, int(round(0.07 * 1000 / frame_ms)))):
        buf.snap(page)
        time.sleep(frame_ms / 1000.0)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "outdir",
        type=Path,
        help="experiment output dir passed to cryodrgn dashboard",
    )
    parser.add_argument("-o", "--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument(
        "--conda-prefix",
        type=Path,
        default=DEFAULT_CONDA_PREFIX,
        help="conda env prefix for conda run -p …",
    )
    parser.add_argument("--fps", type=float, default=40.0)
    parser.add_argument(
        "--duration",
        type=float,
        default=32.0,
        help="encoded playback duration in seconds (default 32 for cache + montage)",
    )
    parser.add_argument("--filter-max", type=int, default=40_000)
    parser.add_argument("--cpus", type=int, default=8)
    parser.add_argument("--max-width", type=int, default=512)
    parser.add_argument("--annotate-caption", action="store_true")
    parser.add_argument(
        "--high-resolution-output",
        type=Path,
        default=None,
        metavar="PATH",
    )
    parser.add_argument(
        "--high-resolution-max-width",
        type=int,
        default=0,
        metavar="PX",
    )
    parser.add_argument("--caption-font-size", type=int, default=21, metavar="PX")
    parser.add_argument(
        "--annotate-min-segment-seconds",
        type=float,
        default=3.5,
        metavar="SEC",
    )
    parser.add_argument(
        "--k-means-fit",
        type=int,
        default=37,
        metavar="K",
        help="k-means cluster count on the rendered scatter trace (4–64, default 37)",
    )
    parser.add_argument(
        "--lasso-regions",
        type=int,
        default=3,
        metavar="N",
        help="how many disjoint k-means groups to lasso (1–8, default 3)",
    )
    parser.add_argument("--headed", action="store_true")
    parser.add_argument("-q", "--quiet", action="store_true")
    args = parser.parse_args()

    conda_prefix = Path(args.conda_prefix)
    fps = float(args.fps)
    wall_s = float(args.duration)
    if fps <= 0 or wall_s <= 0:
        print("error: fps and duration must be positive", file=sys.stderr)
        return 1
    frame_ms = max(1, int(round(1000.0 / fps)))
    k_fit = max(4, min(64, int(args.k_means_fit)))
    pick_n = max(1, min(8, int(args.lasso_regions)))
    log = _make_mr_logger(quiet=args.quiet)

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
        str(conda_prefix),
        "--no-capture-output",
        "cryodrgn",
        "dashboard",
        str(outdir),
        "--particle",
        "--no-browser",
        "--host",
        "127.0.0.1",
        "--port",
        str(port),
        "--filter-max",
        str(int(args.filter_max)),
        "--cpus",
        str(int(args.cpus)),
    ]
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
        log(f"Dashboard subprocess conda -p {conda_prefix} …")
        _wait_http(base + "/", timeout_s=420.0, log=log)
        plotly_cached = _ensure_plotly_cached(log)

        ph, launch_args = _chromium_launch_options(headed=args.headed)
        with sync_playwright() as p:
            with _timed_step(log, "Launch Chromium", slow_note=False):
                browser = p.chromium.launch(headless=ph, args=launch_args)
            context = browser.new_context(
                viewport={"width": 1440, "height": 900},
                device_scale_factor=1,
            )
            if plotly_cached is not None:
                plotly_bytes = plotly_cached.read_bytes()

                def _fulfill(route):
                    route.fulfill(
                        status=200,
                        headers={
                            "content-type": "application/javascript; charset=utf-8",
                            "cache-control": "public, max-age=86400",
                        },
                        body=plotly_bytes,
                    )

                context.route(PLOTLY_CDN_ROUTE_GLOB, _fulfill)
            page = context.new_page()
            with _timed_step(
                log, "Record multi-region particle explorer clip", slow_note=True
            ):
                record_multi_region_sequence(
                    page,
                    base,
                    buf,
                    frame_ms=frame_ms,
                    log=log,
                    k_fit=k_fit,
                    pick_n=pick_n,
                )
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
                log(
                    f"Annotated captions: extended {n_seg_ext} span(s); playback ~{new_play_s:.1f}s"
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
            f"Wrote {std_path} ({nf} frames, ~{fps:g} fps avg, {play_s:.2f} s playback)"
        )
        hr_out = args.high_resolution_output
        if hr_out is not None:
            hr_path = Path(hr_out).expanduser().resolve()
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
