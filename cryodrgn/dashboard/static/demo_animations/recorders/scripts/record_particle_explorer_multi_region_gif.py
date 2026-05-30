#!/usr/bin/env python3
"""Headless GIF: particle explorer — high‑K k-means on the scatter, **hull-shaped** lassos, then image grid.

Fits **k-means with ``k_fit``** (default **37**) on the same Plotly trace points shown in the explorer,
then **lassos ``pick_n``** separated clusters using **irregular polygons** from each cluster's
**2‑D convex hull** (slightly expanded), not axis-aligned boxes. One cluster is chosen for **low
mean ``znorm``**, one for **high mean ``znorm``**, and a third well-separated cluster. Region chips use
fixed **red / blue / yellow** primaries. Lasso cursor motion is much faster than the original slow demo rate.

Then **image cache** + **montage** as before.

Same Playwright / Pillow stack. By default writes two GIFs beside ``recorders/``:

* ``particle_explorer_multi_region_demo.gif`` — **512 px** wide, **16 fps × 20 s** (targets **< 10 MB**).
* ``particle_explorer_multi_region_demo_hr.gif`` — **native** screenshot width, **~40 fps × 36 s** (targets **> 50 MB**).

Use ``--no-high-resolution-output`` to skip the HR file; tune ``--fps`` / ``--duration`` and
``--high-resolution-fps`` / ``--high-resolution-duration`` if a run lands outside those size bands.

Example::

    PYTHONPATH=/path/to/cryodrgn_beta \\
      conda run -p /scratch/gpfs/ZHONGE/mg2332/conda-envs/cdrgn_gifs --no-capture-output \\
        python -m cryodrgn.dashboard.static.demo_animations.recorders.scripts.record_particle_explorer_multi_region_gif \\
        /path/to/train-vae_out/ \\
        --k-means-fit 37 --lasso-regions 3
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
    _montage_gif_ready_count_js,
    _scatter_plot_box,
    _set_gif_suppress_scatter_grid_letters,
    _timed_step,
    _wait_http,
    _install_explorer_scatter_recording_routes,
    _wait_explorer_scatter_ready,
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

DEFAULT_CONDA_PREFIX = Path("/scratch/gpfs/ZHONGE/mg2332/conda-envs/cdrgn_gifs")
DEFAULT_OUTPUT = (
    DEMO_ANIMATIONS_RECORDERS_DIR / "particle_explorer_multi_region_demo.gif"
)
DEFAULT_HIGH_RES_OUTPUT = (
    DEMO_ANIMATIONS_RECORDERS_DIR / "particle_explorer_multi_region_demo_hr.gif"
)
# Encoded playback targets; tune via CLI if a run lands outside the size band.
DEFAULT_FPS = 16.0
DEFAULT_DURATION_S = 20.0
DEFAULT_HIGH_RES_FPS = 40.0
DEFAULT_HIGH_RES_DURATION_S = 36.0
_LOG_PREFIX = "[record-particle-explorer-multi-region-gif]"
EXPLORER_GIF_RECORDER_QUERY = "gif_recorder=1"


def _make_mr_logger(quiet: bool):
    def log(msg: str) -> None:
        if not quiet:
            print(f"{_LOG_PREFIX} {msg}", flush=True)

    return log


def _explorer_gif_recorder_url(base: str) -> str:
    return f"{base.rstrip('/')}/explorer?{EXPLORER_GIF_RECORDER_QUERY}"


def _wait_explorer_gif_hooks(page, timeout_ms: float = 30_000) -> None:
    page.wait_for_function(
        """() => {
      const api = window.__cdrgnExplorerGif;
      return !!(api && typeof api.commitRegionFromDataRing === 'function');
    }""",
        timeout=timeout_ms,
    )


def _gif_file_size_mib(path: Path) -> float:
    try:
        return path.stat().st_size / (1024.0 * 1024.0)
    except OSError:
        return 0.0


# Plotly selection handler debounce (~200 ms) plus margin for chip overlay sync.
# Kept tight so longer montage/grid holds do not balloon total raw capture.
_POST_LASSO_SETTLE_S = 0.24
# Scales k-means intro, lasso traces, per-region holds, and union clip (not montage/cache).
_SELECTION_CLIP_SCALE = 0.5

# Cursor moves along each lasso: baseline was 70% slower; further speed-ups for GIF pacing.
LASSO_DRAW_SLOWDOWN = (1.0 / (1.0 - 0.70)) * 0.7 * 0.7 * 0.5


def _selection_clip_seconds(seconds: float) -> float:
    return seconds * _SELECTION_CLIP_SCALE


# Fixed RGB primary colours for the three region chips / overlays (red, blue, yellow).
EXPLORER_GIF_REGION_COLORS = ("#e41a1c", "#1e88e5", "#fdd835")

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

  function znormAtTraceIdx(ti) {
    const cd = gd.data[0].customdata;
    if (!cd || !cd[ti]) return NaN;
    const row = cd[ti];
    if (row.length > 2) {
      const v = Number(row[2]);
      if (isFinite(v)) return v;
    }
    return NaN;
  }

  const clusters = [];
  for (let j = 0; j < kFit; j++) {
    const cps = [];
    let zsum = 0;
    let zcnt = 0;
    for (let i = 0; i < pts.length; i++) {
      if (labels[i] !== j) continue;
      cps.push(pts[i].slice());
      const z = znormAtTraceIdx(i);
      if (isFinite(z)) {
        zsum += z;
        zcnt++;
      }
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
    const meanZ = zcnt ? zsum / zcnt : NaN;
    clusters.push({ id: j, pts: cps, cx: cx, cy: cy, n: cps.length, meanZ: meanZ });
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

  function alreadyPicked(c) {
    for (let p = 0; p < picked.length; p++) {
      if (picked[p].id === c.id) return true;
    }
    return false;
  }
  function spatiallySeparated(c) {
    const minSep2 = 0.07 * 0.07;
    for (let p = 0; p < picked.length; p++) {
      const dx = c.ncx - picked[p].ncx;
      const dy = c.ncy - picked[p].ncy;
      if (dx * dx + dy * dy < minSep2) return false;
    }
    return true;
  }

  const picked = [];
  const withZ = clusters.filter(function (c) { return isFinite(c.meanZ); });
  if (withZ.length >= 1 && pickN >= 1) {
    const byZ = withZ.slice().sort(function (a, b) { return a.meanZ - b.meanZ; });
    const lowCl = byZ[0];
    const highCl = byZ[byZ.length - 1];
    picked.push(lowCl);
    if (pickN >= 2 && highCl.id !== lowCl.id) picked.push(highCl);
  }
  const sorted = clusters.slice().sort(function (a, b) { return b.n - a.n; });
  for (let i = 0; i < sorted.length && picked.length < pickN; i++) {
    const c = sorted[i];
    if (alreadyPicked(c)) continue;
    if (!spatiallySeparated(c)) continue;
    picked.push(c);
  }
  for (let i = 0; i < sorted.length && picked.length < pickN; i++) {
    const c = sorted[i];
    if (!alreadyPicked(c)) picked.push(c);
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
  const dataPaths = [];
  const sizes = [];
  const ids = [];
  const meanZs = [];
  const roles = [];
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
    dataPaths.push(h.map(function (p) { return [p[0], p[1]]; }));
    const ring = h.map(function (p) {
      const q = dataToFrac(p[0], p[1]);
      return [clamp01(q[0]), clamp01(q[1])];
    });
    paths.push(ring);
    sizes.push(meta.n);
    ids.push(meta.id);
    meanZs.push(meta.meanZ);
  }
  for (let pi = 0; pi < picked.length; pi++) {
    const p = picked[pi];
    let role = 'mid';
    if (isFinite(p.meanZ) && withZ.length >= 2) {
      const byZ = withZ.slice().sort(function (a, b) { return a.meanZ - b.meanZ; });
      if (p.id === byZ[0].id) role = 'low-z';
      else if (p.id === byZ[byZ.length - 1].id) role = 'high-z';
    }
    roles[pi] = role;
  }
  return {
    ok: true,
    paths: paths,
    dataPaths: dataPaths,
    sizes: sizes,
    ids: ids,
    meanZs: meanZs,
    roles: roles,
    kFit: kFit,
    pickN: pickN
  };
}
"""

# Headless recorders call ``window.__cdrgnExplorerGif`` (``_particle_explorer_gif_recorder_hooks.html``,
# included when opening ``/explorer?gif_recorder=1``).
_COMMIT_EXPLORER_REGION_JS = """
(opts) => {
  opts = opts || {};
  const api = window.__cdrgnExplorerGif;
  if (!api || typeof api.commitRegionFromDataRing !== 'function') {
    return { ok: false, reason: 'explorer GIF hooks missing' };
  }
  let ring = opts.dataRing;
  if ((!ring || ring.length < 3) && opts.relRing && opts.relRing.length >= 3) {
    const gd0 = document.getElementById('scatter');
    const fl = gd0 && gd0._fullLayout;
    const xa = fl && fl.xaxis;
    const ya = fl && fl.yaxis;
    const gw = (fl && fl.width) || (gd0 && gd0.clientWidth) || 720;
    const gh = (fl && fl.height) || (gd0 && gd0.clientHeight) || 480;
    function fracToData(fx, fy) {
      if (xa && ya && typeof xa.p2l === 'function' && typeof ya.p2l === 'function') {
        const px = fx * gw;
        const py = fy * gh;
        return [
          xa.p2l(px - (xa._offset || 0)),
          ya.p2l(py - (ya._offset || 0))
        ];
      }
      const xr = xa.range;
      const yr = ya.range;
      const ml = (fl.margin && fl.margin.l) || 0;
      const mr = (fl.margin && fl.margin.r) || 0;
      const mt = (fl.margin && fl.margin.t) || 0;
      const mb = (fl.margin && fl.margin.b) || 0;
      const plotW = Math.max(1e-6, gw - ml - mr);
      const plotH = Math.max(1e-6, gh - mt - mb);
      const xd = xr[0] + ((fx * gw - ml) / plotW) * (xr[1] - xr[0]);
      const yd = yr[1] - ((fy * gh - mt) / plotH) * (yr[1] - yr[0]);
      return [xd, yd];
    }
    ring = opts.relRing.map(function (p) { return fracToData(p[0], p[1]); });
  }
  return api.commitRegionFromDataRing(ring, opts.colorHex || null);
}
"""

_ENSURE_SCATTER_COLOR_ZNORM_JS = """
() => {
  const sc = document.getElementById('sc');
  if (!sc) return false;
  let zv = null;
  for (let i = 0; i < sc.options.length; i++) {
    const o = sc.options[i];
    if (o.value === 'znorm' || /^znorm\\b/i.test(String(o.value || ''))) {
      zv = o.value;
      break;
    }
  }
  if (!zv) return false;
  if (sc.value !== zv) {
    sc.value = zv;
    sc.dispatchEvent(new Event('change', { bubbles: true }));
  }
  return true;
}
"""

_CLEAR_EXPLORER_REGIONS_JS = """
() => {
  const api = window.__cdrgnExplorerGif;
  if (api && typeof api.clearRegions === 'function') api.clearRegions();
  return true;
}
"""

_SET_SUPPRESS_PLOTLY_SELECTION_JS = """
(on) => {
  const api = window.__cdrgnExplorerGif;
  if (!api || typeof api.setSuppressPlotlySelection !== 'function') return false;
  api.setSuppressPlotlySelection(!!on);
  return true;
}
"""

_CLEAR_PLOTLY_TRANSIENT_SELECTIONS_JS = """
() => {
  const api = window.__cdrgnExplorerGif;
  if (api && typeof api.clearPlotlyTransientSelections === 'function') {
    api.clearPlotlyTransientSelections();
    return true;
  }
  return false;
}
"""

# Must exceed explorer ``LASSO_SELECTION_DEBOUNCE_MS`` (200) after the last synthetic lasso.
_LASSO_DEBOUNCE_QUIET_S = 0.45  # multiplied by ``_SELECTION_CLIP_SCALE`` at use sites

# Main image-grid hold multipliers (×2 linger on first open and after layout change).
_MONTAGE_GRID_DWELL_MULT = 2.0


def _parse_kmeans_path_rings(
    paths_raw: object,
) -> list[list[tuple[float, float]]]:
    out: list[list[tuple[float, float]]] = []
    if not isinstance(paths_raw, list):
        return out
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
    return out


def _kmeans_lasso_paths_from_page(
    page,
    *,
    k_fit: int,
    pick_n: int,
    log,
) -> tuple[list[list[tuple[float, float]]], list[list[list[float]]]] | None:
    """Return fraction paths (for cursor motion) and data-space hulls (for selection commit)."""
    raw = page.evaluate(
        _KMEANS_LASSO_PATHS_JS,
        {"kFit": int(k_fit), "pickN": int(pick_n)},
    )
    if not isinstance(raw, dict) or not raw.get("ok"):
        reason = raw.get("reason") if isinstance(raw, dict) else "bad payload"
        log(f"k-means lasso path compute failed ({reason}); using fallback polygons")
        return None
    rel_paths = _parse_kmeans_path_rings(raw.get("paths"))
    data_paths = _parse_kmeans_path_rings(raw.get("dataPaths"))
    if len(rel_paths) < 1:
        log("k-means returned no paths; using fallback polygons")
        return None
    if len(data_paths) != len(rel_paths):
        data_paths = []
    sizes = raw.get("sizes")
    ids = raw.get("ids")
    roles = raw.get("roles")
    mean_zs = raw.get("meanZs")
    if isinstance(sizes, list) and isinstance(ids, list) and len(sizes) == len(ids):
        parts = [
            f"k-means fit K={raw.get('kFit', '?')}: lasso {len(rel_paths)} region(s); "
            "cluster ids "
            + ", ".join(str(int(i)) for i in ids if isinstance(i, (int, float)))
            + " — sizes "
            + ", ".join(str(int(s)) for s in sizes if isinstance(s, (int, float)))
        ]
        if isinstance(roles, list) and roles:
            parts.append("roles " + ", ".join(str(r) for r in roles))
        if isinstance(mean_zs, list) and mean_zs:
            zfmt = ", ".join(
                f"{float(z):.3g}" if isinstance(z, (int, float)) else "?"
                for z in mean_zs
            )
            parts.append(f"mean znorm [{zfmt}]")
        log("; ".join(parts))
    return rel_paths, data_paths


def _commit_explorer_region(
    page,
    data_ring: list[list[float]] | None,
    *,
    rel_ring: list[tuple[float, float]] | None,
    color_hex: str | None,
    log,
) -> None:
    payload: dict = {}
    if data_ring:
        payload["dataRing"] = data_ring
    elif rel_ring:
        payload["relRing"] = [[float(x), float(y)] for x, y in rel_ring]
    else:
        raise RuntimeError("commit region needs data_ring or rel_ring")
    if color_hex:
        payload["colorHex"] = str(color_hex)
    result = page.evaluate(_COMMIT_EXPLORER_REGION_JS, payload)
    if not isinstance(result, dict) or not result.get("ok"):
        reason = result.get("reason") if isinstance(result, dict) else "bad payload"
        raise RuntimeError(f"could not commit explorer region ({reason})")
    log(
        f"committed region {result.get('nRegions', '?')}: "
        f"{result.get('regionRows', '?')} rows in hull, "
        f"{result.get('nSelected', '?')} selected total"
    )


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
    post_settle_s: float = _POST_LASSO_SETTLE_S,
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
    time.sleep(post_settle_s)


# Montage / grid dwells: first grid, menu, mid re-layout, final grid were scaled up earlier; this
# pass multiplies those holds by **1.3** again. Intro, lasso pauses, post-lasso settle, union, cache
# holds, montage poll budget, and tail are tightened to limit raw capture growth before normalization.
def _volume_animate_montage_clip(
    buf: CaptionedExplorerBuffer,
    page,
    *,
    frame_ms: int,
    log,
) -> None:
    """Generate volume stills, then snapshot each cell as its rotating GIF finishes (ChimeraX)."""
    gen = page.locator("#btn-volume-generate")
    if gen.count() == 0:
        return
    vol_toggle = page.locator("#volumes-panel-toggle")
    if vol_toggle.count():
        expanded = vol_toggle.get_attribute("aria-expanded") or "false"
        if expanded.strip().lower() == "false":
            vol_toggle.click()
            for _ in range(max(1, int(round(0.35 * 1000 / frame_ms)))):
                buf.snap(page)
                time.sleep(frame_ms / 1000.0)
    if not gen.is_enabled():
        log("Generate volumes stayed disabled; skipping volume animate clip")
        return
    with _timed_step(log, "Generate volume stills for montage grid", slow_note=True):
        gen.click()
        try:
            page.wait_for_function(
                """() => {
                  const b = document.getElementById('btn-volume-animate');
                  return b && !b.disabled;
                }""",
                timeout=600_000,
            )
        except Exception:
            log("volume stills did not finish in time; skipping animate clip")
            return
        for _ in range(max(1, int(round(0.25 * 1000 / frame_ms)))):
            buf.snap(page)
            time.sleep(frame_ms / 1000.0)
    _stamp_step(
        buf,
        "Animate volumes — each grid cell starts rotating as soon as its GIF is ready.",
    )
    started = page.evaluate(
        """() => {
      const api = window.__cdrgnExplorerGif;
      if (!api || typeof api.startVolumeAnimateSequential !== 'function') {
        return { ok: false };
      }
      return api.startVolumeAnimateSequential();
    }"""
    )
    if not isinstance(started, dict) or not started.get("ok"):
        anim = page.locator("#btn-volume-animate")
        if anim.count() and anim.is_enabled():
            anim.click()
        else:
            log("Animate volumes unavailable; skipping per-cell GIF clip")
            return
    with _timed_step(
        log, "Capture per-cell volume GIFs as they finish", slow_note=True
    ):
        _snap_montage_load(
            buf,
            page,
            frame_ms=frame_ms,
            max_polls=600,
            count_js=_montage_gif_ready_count_js(),
        )
        for _ in range(max(1, int(round(1.0 * 1000 / frame_ms)))):
            buf.snap(page)
            time.sleep(frame_ms / 1000.0)


def _image_cache_and_montage_clip(
    buf: CaptionedExplorerBuffer,
    page,
    *,
    frame_ms: int,
    log,
) -> None:
    """Expand cache, preload selection, open montage — mirrors ``record_particle_explorer_actions_gif``."""
    _set_gif_suppress_scatter_grid_letters(page, True)
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
            dwell_s = 3.536 * _MONTAGE_GRID_DWELL_MULT
            for _ in range(max(1, int(round(dwell_s * 1000 / frame_ms)))):
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
                dwell2_s = 4.16 * _MONTAGE_GRID_DWELL_MULT
                for _ in range(max(1, int(round(dwell2_s * 1000 / frame_ms)))):
                    buf.snap(page)
                    time.sleep(frame_ms / 1000.0)
            _volume_animate_montage_clip(buf, page, frame_ms=frame_ms, log=log)
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
    page.goto(
        _explorer_gif_recorder_url(base),
        wait_until="domcontentloaded",
        timeout=180_000,
    )
    with _timed_step(log, "Wait for GIF recorder hooks", slow_note=False):
        _wait_explorer_gif_hooks(page, 30_000)
    with _timed_step(log, "Wait for scatter render", slow_note=True):
        _wait_explorer_scatter_ready(page, 300_000)

    _stamp_step(
        buf,
        f"Particle explorer — k-means with K={k_fit} on the visible scatter; "
        f"lasso {pick_n} separated clusters, then image cache + montage.",
    )
    for _ in range(max(1, int(round(_selection_clip_seconds(0.06) * 1000 / frame_ms)))):
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

    if page.evaluate(_ENSURE_SCATTER_COLOR_ZNORM_JS):
        with _timed_step(
            log, "Colour scatter by znorm for cluster picking", slow_note=True
        ):
            _wait_explorer_scatter_ready(page, 300_000)
    else:
        log("znorm not in colour dropdown; cluster pick uses size/separation only")

    km = _kmeans_lasso_paths_from_page(page, k_fit=k_fit, pick_n=pick_n, log=log)
    if km is None:
        rel_paths = _fallback_lasso_demo_paths()
        data_paths: list[list[list[float]]] = []
    else:
        rel_paths, data_paths = km

    page.evaluate(_CLEAR_EXPLORER_REGIONS_JS)
    page.evaluate(_SET_SUPPRESS_PLOTLY_SELECTION_JS, True)
    time.sleep(_selection_clip_seconds(0.12))

    lasso_slowdown = LASSO_DRAW_SLOWDOWN * _SELECTION_CLIP_SCALE
    lasso_settle_s = _selection_clip_seconds(_POST_LASSO_SETTLE_S)
    debounce_quiet_s = _selection_clip_seconds(_LASSO_DEBOUNCE_QUIET_S)

    captions = [
        f"Lasso {i + 1} — convex hull of one of {pick_n} k-means groups (K={k_fit}), traced slowly."
        for i in range(len(rel_paths))
    ]

    for idx, rel in enumerate(rel_paths):
        box = _scatter_plot_box(page) or box
        _stamp_step(buf, captions[idx])
        _lasso_polygon(
            page,
            box,
            rel,
            buf,
            frame_ms=frame_ms,
            intro_snaps=0,
            slowdown=lasso_slowdown,
            post_settle_s=lasso_settle_s,
        )
        page.evaluate(_CLEAR_PLOTLY_TRANSIENT_SELECTIONS_JS)
        data_ring = data_paths[idx] if idx < len(data_paths) else None
        color_hex = EXPLORER_GIF_REGION_COLORS[idx % len(EXPLORER_GIF_REGION_COLORS)]
        with _timed_step(log, f"Commit lasso region {idx + 1}", slow_note=False):
            _commit_explorer_region(
                page,
                data_ring,
                rel_ring=rel,
                color_hex=color_hex,
                log=log,
            )
        time.sleep(lasso_settle_s)
        hold = _selection_clip_seconds(0.24 if idx < len(rel_paths) - 1 else 0.12)
        for _ in range(max(1, int(round(hold * 1000 / frame_ms)))):
            buf.snap(page)
            time.sleep(frame_ms / 1000.0)

    time.sleep(debounce_quiet_s)
    page.evaluate(_CLEAR_PLOTLY_TRANSIENT_SELECTIONS_JS)
    page.evaluate(_SET_SUPPRESS_PLOTLY_SELECTION_JS, False)
    time.sleep(debounce_quiet_s)

    _stamp_step(
        buf,
        "Union selection — region chips for each lasso; combined particles drive cache + grid below.",
    )
    for _ in range(max(1, int(round(_selection_clip_seconds(0.12) * 1000 / frame_ms)))):
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
    parser.add_argument(
        "--fps",
        type=float,
        default=DEFAULT_FPS,
        help=f"standard GIF frame rate (default {DEFAULT_FPS:g}, targets < 10 MB at 512 px)",
    )
    parser.add_argument(
        "--duration",
        type=float,
        default=DEFAULT_DURATION_S,
        help=(
            f"standard GIF playback seconds (default {DEFAULT_DURATION_S:g} for cache + montage)"
        ),
    )
    parser.add_argument(
        "--filter-max",
        type=int,
        default=28_000,
        help="``cryodrgn dashboard --filter-max`` (lower default for SVG scatter recording)",
    )
    parser.add_argument("--cpus", type=int, default=8)
    parser.add_argument("--max-width", type=int, default=512)
    parser.add_argument("--annotate-caption", action="store_true")
    parser.add_argument(
        "--high-resolution-output",
        type=Path,
        default=DEFAULT_HIGH_RES_OUTPUT,
        metavar="PATH",
        help=(
            "second GIF at native width (default "
            f"{DEFAULT_HIGH_RES_OUTPUT.name}; use --no-high-resolution-output to skip)"
        ),
    )
    parser.add_argument(
        "--no-high-resolution-output",
        action="store_true",
        help="do not write the high-resolution GIF",
    )
    parser.add_argument(
        "--high-resolution-fps",
        type=float,
        default=DEFAULT_HIGH_RES_FPS,
        help=f"HR GIF frame rate (default {DEFAULT_HIGH_RES_FPS:g}, targets > 50 MB)",
    )
    parser.add_argument(
        "--high-resolution-duration",
        type=float,
        default=DEFAULT_HIGH_RES_DURATION_S,
        help=(f"HR GIF playback seconds (default {DEFAULT_HIGH_RES_DURATION_S:g})"),
    )
    parser.add_argument(
        "--high-resolution-max-width",
        type=int,
        default=0,
        metavar="PX",
        help="max width for HR GIF (0 = native screenshot width)",
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
    hr_fps = float(args.high_resolution_fps)
    hr_wall_s = float(args.high_resolution_duration)
    if fps <= 0 or wall_s <= 0:
        print("error: fps and duration must be positive", file=sys.stderr)
        return 1
    if hr_fps <= 0 or hr_wall_s <= 0:
        print(
            "error: high-resolution fps and duration must be positive", file=sys.stderr
        )
        return 1
    frame_ms = max(1, int(round(1000.0 / fps)))
    write_hr = not args.no_high_resolution_output
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
            _install_explorer_scatter_recording_routes(context, log=log)
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
        std_mib = _gif_file_size_mib(std_path)
        print(
            f"Wrote {std_path} ({nf} frames, ~{fps:g} fps avg, {play_s:.2f} s playback, "
            f"{std_mib:.2f} MiB)"
        )
        if std_mib >= 10.0:
            log(
                f"standard GIF is {std_mib:.1f} MiB (target < 10 MiB); try lower --fps or --duration"
            )
        if write_hr:
            hr_path = Path(args.high_resolution_output).expanduser().resolve()
            hr_mw = max(0, int(args.high_resolution_max_width))
            pngs_hr, durs_hr, caps_hr = normalize_frames_fixed_fps(
                buf.pngs,
                fps=hr_fps,
                wall_seconds=hr_wall_s,
                captions=(
                    list(buf.captions_per_frame) if args.annotate_caption else None
                ),
            )
            cap_lines_hr = caps_hr if annotate else None
            if annotate and cap_lines_hr is not None and min_seg > 0:
                (
                    pngs_hr,
                    durs_hr,
                    cap_lines_hr,
                    n_seg_ext_hr,
                ) = extend_caption_segments_min_duration(
                    pngs_hr,
                    durs_hr,
                    cap_lines_hr,
                    min_segment_s=min_seg,
                    fps=hr_fps,
                )
                if n_seg_ext_hr:
                    log(
                        f"HR captions: extended {n_seg_ext_hr} span(s); "
                        f"playback ~{sum(durs_hr) / 1000.0:.1f}s"
                    )
            save_pngs_duration_gif(
                pngs_hr,
                durs_hr,
                hr_path,
                max_width=hr_mw,
                log=log,
                annotate=annotate,
                step_captions=cap_lines_hr,
                caption_font_px=caption_font_px,
            )
            hr_mib = _gif_file_size_mib(hr_path)
            tag = "full native width" if hr_mw == 0 else f"max width {hr_mw}px"
            print(
                f"Wrote {hr_path} ({len(pngs_hr)} frames, ~{hr_fps:g} fps, "
                f"{sum(durs_hr) / 1000.0:.2f} s playback, {hr_mib:.2f} MiB, {tag})"
            )
            if hr_mib <= 50.0:
                log(
                    f"high-resolution GIF is {hr_mib:.1f} MiB (target > 50 MiB); "
                    "try higher --high-resolution-fps / --high-resolution-duration"
                )
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
