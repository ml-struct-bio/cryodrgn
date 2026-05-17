#!/usr/bin/env python3
"""Headless capture: volume sketched landscape colour changes, lasso, and cycle GIFs.

The recorder opens the landscape explorer directly, so the launch page is not
included. It captures a short sequence for the volume PCA landscape interface:

1. choose the agglomerative state colour covariate
2. make a lasso selection
3. choose the znorm colour covariate

Loading and ChimeraX rendering clips are sampled for at most 0.4 seconds each;
the script waits silently for the real work to finish before recording the next
finished state. It leaves the animation output in cycle mode, so rotating volume
GIF previews are never shown.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from cryodrgn.dashboard.static.demo_animations.recorders.scripts.record_dashboard_interactions_gif import (
    DEMO_ANIMATIONS_RECORDERS_DIR,
    PLOTLY_CDN_ROUTE_GLOB,
    REPO_ROOT,
    FrameBuffer,
    _chromium_launch_options,
    _ensure_plotly_cached,
    _free_port,
    _save_buf_to_gif,
    _timed_step,
    _wait_http,
)

FRAME_MS_DEFAULT = 50
MAX_GIF_WALL_MS = 6_000
RENDER_CLIP_MAX_MS = 400

_LOG_PREFIX = "[record-landscape-volpca-gif]"


def _make_logger(quiet: bool):
    def log(msg: str) -> None:
        if not quiet:
            print(f"{_LOG_PREFIX} {msg}", flush=True)

    return log


def _route_landscape_scatter_as_svg(context) -> None:
    """Rewrite capture-only Plotly traces from WebGL ``scattergl`` to SVG ``scatter``."""

    def _fulfill_svg_scatter(route):
        try:
            resp = route.fetch()
            body = resp.text()
            fig = json.loads(body)
            for trace in fig.get("data", []):
                if trace.get("type") == "scattergl":
                    trace["type"] = "scatter"
            headers = dict(resp.headers)
            headers["content-type"] = "application/json"
            headers.pop("content-length", None)
            route.fulfill(
                status=resp.status,
                headers=headers,
                body=json.dumps(fig, separators=(",", ":")),
            )
        except Exception:
            route.continue_()

    context.route("**/api/landscape_volpca/scatter**", _fulfill_svg_scatter)


def _scatter_ready_js() -> str:
    return """() => {
      const gd = document.getElementById('volsketch');
      const status = document.getElementById('volsketch-plot-status');
      const text = status ? (status.textContent || '') : '';
      return !!(
        window.Plotly
        && gd
        && gd.data
        && gd.data.length
        && gd.data[0]
        && gd.data[0].x
        && gd.data[0].x.length
        && text.indexOf('Loading scatter') < 0
      );
    }"""


def _animation_ready_js() -> str:
    return """() => {
      const grid = document.getElementById('volsketch-preview-grid');
      const progress = document.getElementById('volsketch-animate-progress');
      const status = document.getElementById('volsketch-animate-status');
      const text = status ? (status.textContent || '') : '';
      if (!grid) return false;
      const imgs = Array.from(grid.querySelectorAll('img'));
      if (!imgs.length) return false;
      const busy = progress && !progress.hidden;
      const rendering = /Rendering|Re-rendering/.test(text);
      return !busy && !rendering && imgs.every((im) => im.complete && im.naturalWidth > 0);
    }"""


def _direct_render_done_js() -> str:
    return """() => !!window.__cdrgnLandscapeGifDone"""


def _make_volume_points_visible(page) -> None:
    """Bump marker/text styling for readability in the recorded GIF only."""

    page.evaluate(
        """() => {
          const gd = document.getElementById('volsketch');
          if (!gd || !window.Plotly || !gd.data || !gd.data[0]) return false;
          const tr = gd.data[0];
          const n = tr.x ? tr.x.length : 0;
          if (!n) return false;
          const selected = new Set(Array.isArray(tr.selectedpoints) ? tr.selectedpoints : []);
          const hasSelection = selected.size > 0;
          const sizes = new Array(n);
          const opacities = new Array(n);
          for (let i = 0; i < n; i++) {
            const on = selected.has(i);
            sizes[i] = on ? 18 : 14;
            opacities[i] = on ? 0.96 : (hasSelection ? 0.78 : 0.90);
          }
          const update = {
            'marker.size': [sizes],
            'marker.opacity': [opacities],
            'marker.line.width': [1.4],
            'marker.line.color': ['rgba(255,255,255,0.95)'],
            textfont: [{
              size: 16,
              color: '#111827',
              family: 'system-ui, Segoe UI, sans-serif',
            }],
          };
          Plotly.restyle(gd, update, [0]).catch(() => {});
          return true;
        }"""
    )


def _wait_scatter_ready(page, timeout_ms: float = 180_000) -> None:
    page.wait_for_function(_scatter_ready_js(), timeout=timeout_ms)


def _wait_animation_ready(page, timeout_ms: float = 300_000) -> None:
    page.wait_for_function(_animation_ready_js(), timeout=timeout_ms)


def _snap_until_ready_capped(
    buf: FrameBuffer,
    page,
    ready_js: str,
    *,
    max_clip_ms: int = RENDER_CLIP_MAX_MS,
    timeout_ms: int = 180_000,
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
    page.wait_for_function(ready_js, timeout=1)


def _select_option_by_token(page, selector: str, token: str) -> str:
    """Select the first option whose value or label matches ``token``."""

    opts = page.eval_on_selector(
        selector,
        """(sel) => Array.from(sel.options).map((o) => ({
          value: o.value,
          label: o.textContent || '',
          disabled: !!o.disabled,
        }))""",
    )
    needle = token.strip().lower()
    candidates: list[dict[str, Any]] = [o for o in opts if not o.get("disabled")]

    def score(opt: dict[str, Any]) -> int:
        value = str(opt.get("value", "")).lower()
        label = str(opt.get("label", "")).lower()
        if value == needle:
            return 0
        if label == needle:
            return 1
        if needle in value:
            return 2
        if needle in label:
            return 3
        return 99

    ranked = sorted(candidates, key=score)
    if not ranked or score(ranked[0]) == 99:
        labels = ", ".join(
            f"{o.get('value')} ({str(o.get('label', '')).strip()})" for o in candidates
        )
        raise RuntimeError(
            f"No selectable colour option matching {token!r}; options: {labels}"
        )
    chosen = str(ranked[0]["value"])
    page.locator(selector).select_option(value=chosen)
    return chosen


def _choose_color(buf: FrameBuffer, page, token: str, *, log) -> str:
    with _timed_step(log, f"Choose colour covariate {token!r}", slow_note=False):
        sel = page.locator("#volsketch-color")
        sel.scroll_into_view_if_needed()
        box = sel.bounding_box()
        if box:
            page.mouse.move(
                box["x"] + box["width"] * 0.75, box["y"] + box["height"] * 0.5
            )
            buf.sleep_snap(page, 0.12)
            page.mouse.click(
                box["x"] + box["width"] * 0.75, box["y"] + box["height"] * 0.5
            )
            buf.sleep_snap(page, 0.08)
        chosen = _select_option_by_token(page, "#volsketch-color", token)
        _snap_until_ready_capped(buf, page, _scatter_ready_js(), timeout_ms=180_000)
        _wait_scatter_ready(page)
        _make_volume_points_visible(page)
        return chosen


def _force_cycle_mode(page) -> None:
    page.evaluate(
        """() => {
          const cyc = document.getElementById('volsketch-mode-cycle');
          if (!cyc) return false;
          if (!cyc.checked) {
            cyc.checked = true;
            cyc.dispatchEvent(new Event('change', { bubbles: true }));
          }
          return true;
        }"""
    )


def _lasso_path_from_plot(page) -> list[tuple[float, float]]:
    """Return viewport coordinates for a compact lasso around central points."""

    path = page.evaluate(
        """() => {
          const gd = document.getElementById('volsketch');
          if (!gd || !gd.data || !gd.data[0] || !gd._fullLayout) return null;
          const tr = gd.data[0];
          const xa = gd._fullLayout.xaxis;
          const ya = gd._fullLayout.yaxis;
          const bb = gd.getBoundingClientRect();
          if (!tr.x || !tr.y || !tr.x.length || !xa || !ya) return null;
          const pts = [];
          for (let i = 0; i < tr.x.length; i++) {
            const x = Number(tr.x[i]);
            const y = Number(tr.y[i]);
            if (Number.isFinite(x) && Number.isFinite(y)) pts.push([x, y]);
          }
          if (pts.length < 3) return null;
          function q(vals, p) {
            const a = vals.slice().sort((u, v) => u - v);
            const idx = Math.max(0, Math.min(a.length - 1, Math.round((a.length - 1) * p)));
            return a[idx];
          }
          const xs = pts.map((p) => p[0]);
          const ys = pts.map((p) => p[1]);
          const cx = q(xs, 0.50);
          const cy = q(ys, 0.50);
          const rx = Math.max((q(xs, 0.70) - q(xs, 0.30)) * 0.70, (q(xs, 0.90) - q(xs, 0.10)) * 0.12);
          const ry = Math.max((q(ys, 0.70) - q(ys, 0.30)) * 0.70, (q(ys, 0.90) - q(ys, 0.10)) * 0.12);
          const rel = [
            [-0.95, -0.05],
            [-0.52, -0.72],
            [ 0.28, -0.82],
            [ 0.92, -0.30],
            [ 0.74,  0.48],
            [ 0.05,  0.88],
            [-0.70,  0.52],
            [-1.00,  0.00],
          ];
          return rel.map(([dx, dy]) => [
            bb.left + xa._offset + xa.l2p(cx + dx * rx),
            bb.top + ya._offset + ya.l2p(cy + dy * ry),
          ]);
        }"""
    )
    if path and len(path) >= 3:
        return [(float(x), float(y)) for x, y in path]

    box = page.locator("#volsketch").bounding_box()
    if not box:
        raise RuntimeError("Could not locate #volsketch for lasso selection.")
    w, h = box["width"], box["height"]
    rel = [
        (0.34, 0.48),
        (0.44, 0.38),
        (0.58, 0.42),
        (0.64, 0.55),
        (0.52, 0.64),
        (0.38, 0.58),
        (0.32, 0.50),
    ]
    return [(box["x"] + w * x, box["y"] + h * y) for x, y in rel]


def _selected_count(page) -> int:
    text = page.locator("#volsketch-sel-summary").inner_text(timeout=5_000)
    try:
        return int(text.split("/", 1)[0].strip())
    except (ValueError, IndexError):
        return 0


def _selected_vol_ids(page) -> list[int]:
    vals = page.evaluate(
        """() => {
          const gd = document.getElementById('volsketch');
          const fallback = Array.isArray(window.__cdrgnLandscapeSelectedVols)
            ? window.__cdrgnLandscapeSelectedVols
            : [];
          if (!gd || !gd.data || !gd.data[0]) return fallback;
          const tr = gd.data[0];
          const idxs = Array.isArray(tr.selectedpoints) ? tr.selectedpoints.slice() : [];
          if (!idxs.length && tr.text && tr.text.length) {
            for (let i = 0; i < tr.text.length; i++) {
              if (tr.text[i] != null && String(tr.text[i]).trim()) idxs.push(i);
            }
          }
          if (!idxs.length) return fallback;
          function volAt(i) {
            if (tr.ids && tr.ids[i] != null && String(tr.ids[i]) !== '') {
              const v = parseInt(String(tr.ids[i]), 10);
              if (!Number.isNaN(v)) return v;
            }
            const cd = tr.customdata;
            if (cd && cd[i] != null) {
              const row = Array.isArray(cd[i]) ? cd[i] : Array.from(cd[i] || []);
              if (row.length) {
                const v = parseInt(String(row[0]), 10);
                if (!Number.isNaN(v)) return v;
              }
            }
            return NaN;
          }
          return Array.from(new Set(idxs.map(volAt).filter((v) => !Number.isNaN(v))))
            .sort((a, b) => a - b);
        }"""
    )
    return [int(v) for v in vals]


def _ensure_visible_selection(page) -> int:
    """Fallback for headless Plotly lasso misses: select central points visibly."""

    return int(
        page.evaluate(
            """() => {
              const gd = document.getElementById('volsketch');
              if (!gd || !window.Plotly || !gd.data || !gd.data[0]) return 0;
              const tr = gd.data[0];
              const n = tr.x ? tr.x.length : 0;
              if (!n) return 0;
              const existing = Array.isArray(tr.selectedpoints) ? tr.selectedpoints.slice() : [];
              if (existing.length) return existing.length;
              const pts = [];
              for (let i = 0; i < n; i++) {
                const x = Number(tr.x[i]);
                const y = Number(tr.y[i]);
                if (Number.isFinite(x) && Number.isFinite(y)) pts.push({i, x, y});
              }
              if (!pts.length) return 0;
              function q(vals, p) {
                const a = vals.slice().sort((u, v) => u - v);
                const idx = Math.max(0, Math.min(a.length - 1, Math.round((a.length - 1) * p)));
                return a[idx];
              }
              const xs = pts.map((p) => p.x);
              const ys = pts.map((p) => p.y);
              const cx = q(xs, 0.50);
              const cy = q(ys, 0.50);
              const sx = Math.max(q(xs, 0.75) - q(xs, 0.25), 1e-9);
              const sy = Math.max(q(ys, 0.75) - q(ys, 0.25), 1e-9);
              pts.forEach((p) => {
                p.d = Math.pow((p.x - cx) / sx, 2) + Math.pow((p.y - cy) / sy, 2);
              });
              pts.sort((a, b) => a.d - b.d);
              const idxs = pts.slice(0, Math.min(8, pts.length)).map((p) => p.i).sort((a, b) => a - b);
              const selected = new Set(idxs);
              const selectedVols = [];
              function volAt(i) {
                if (tr.ids && tr.ids[i] != null && String(tr.ids[i]) !== '') {
                  const v = parseInt(String(tr.ids[i]), 10);
                  if (!Number.isNaN(v)) return v;
                }
                const cd = tr.customdata;
                if (cd && cd[i] != null) {
                  const row = Array.isArray(cd[i]) ? cd[i] : Array.from(cd[i] || []);
                  if (row.length) {
                    const v = parseInt(String(row[0]), 10);
                    if (!Number.isNaN(v)) return v;
                  }
                }
                return NaN;
              }
              idxs.forEach((idx) => {
                const v = volAt(idx);
                if (!Number.isNaN(v)) selectedVols.push(v);
              });
              window.__cdrgnLandscapeSelectedVols = selectedVols;
              const safeLetters = 'ABCDEFGHJKLMNPQRSTVWXYZ';
              const labels = new Array(n).fill('');
              idxs.forEach((idx, j) => { labels[idx] = safeLetters[j] || String(j + 1); });
              const sizes = new Array(n);
              const opacities = new Array(n);
              for (let i = 0; i < n; i++) {
                const on = selected.has(i);
                sizes[i] = on ? 18 : 14;
                opacities[i] = on ? 0.96 : 0.78;
              }
              const update = {
                selectedpoints: [idxs],
                mode: ['markers+text'],
                text: [labels],
                textposition: ['top center'],
                textfont: [{
                  size: 16,
                  color: '#111827',
                  family: 'system-ui, Segoe UI, sans-serif',
                }],
                'marker.size': [sizes],
                'marker.opacity': [opacities],
                'marker.line.width': [1.4],
                'marker.line.color': ['rgba(255,255,255,0.95)'],
              };
              Plotly.restyle(gd, update, [0]).catch(() => {});
              const summary = document.getElementById('volsketch-sel-summary');
              if (summary) {
                const total = tr.x ? tr.x.length : n;
                summary.textContent = idxs.length + '/' + total + ' sketched volumes selected';
              }
              return idxs.length;
            }"""
        )
    )


def _restore_visible_selection_for_vols(page, vols: list[int]) -> None:
    page.evaluate(
        """(vols) => {
          window.__cdrgnLandscapeSelectedVols = vols.slice();
          const gd = document.getElementById('volsketch');
          if (!gd || !window.Plotly || !gd.data || !gd.data[0]) return false;
          const tr = gd.data[0];
          const n = tr.x ? tr.x.length : 0;
          const volSet = new Set(vols.map((v) => Number(v)));
          function volAt(i) {
            if (tr.ids && tr.ids[i] != null && String(tr.ids[i]) !== '') {
              const v = parseInt(String(tr.ids[i]), 10);
              if (!Number.isNaN(v)) return v;
            }
            const cd = tr.customdata;
            if (cd && cd[i] != null) {
              const row = Array.isArray(cd[i]) ? cd[i] : Array.from(cd[i] || []);
              if (row.length) {
                const v = parseInt(String(row[0]), 10);
                if (!Number.isNaN(v)) return v;
              }
            }
            return NaN;
          }
          const idxs = [];
          for (let i = 0; i < n; i++) {
            const v = volAt(i);
            if (!Number.isNaN(v) && volSet.has(v)) idxs.push(i);
          }
          if (!idxs.length) return false;
          const safeLetters = 'ABCDEFGHJKLMNPQRSTVWXYZ';
          const labels = new Array(n).fill('');
          idxs.forEach((idx, j) => { labels[idx] = safeLetters[j] || String(j + 1); });
          const selected = new Set(idxs);
          const sizes = new Array(n);
          const opacities = new Array(n);
          for (let i = 0; i < n; i++) {
            const on = selected.has(i);
            sizes[i] = on ? 18 : 14;
            opacities[i] = on ? 0.96 : 0.78;
          }
          Plotly.restyle(gd, {
            selectedpoints: [idxs],
            mode: ['markers+text'],
            text: [labels],
            textposition: ['top center'],
            textfont: [{
              size: 16,
              color: '#111827',
              family: 'system-ui, Segoe UI, sans-serif',
            }],
            'marker.size': [sizes],
            'marker.opacity': [opacities],
            'marker.line.width': [1.4],
            'marker.line.color': ['rgba(255,255,255,0.95)'],
          }, [0]).catch(() => {});
          const summary = document.getElementById('volsketch-sel-summary');
          if (summary) summary.textContent = idxs.length + '/' + n + ' sketched volumes selected';
          return true;
        }""",
        vols,
    )


def _render_cycle_preview_from_selection(
    buf: FrameBuffer,
    page,
    *,
    keep_previews: bool,
    log,
) -> None:
    vols = _selected_vol_ids(page)
    if not vols:
        raise RuntimeError(
            "No selected sketch volume IDs were available for cycle rendering."
        )
    _restore_visible_selection_for_vols(page, vols)
    with _timed_step(log, "Render cycle preview through landscape API", slow_note=True):
        page.evaluate(
            """({vols, keepPreviews}) => {
              window.__cdrgnLandscapeGifDone = false;
              window.__cdrgnLandscapeGifError = '';
              const colorSel = document.getElementById('volsketch-color');
              const status = document.getElementById('volsketch-animate-status');
              const progress = document.getElementById('volsketch-animate-progress');
              const grid = document.getElementById('volsketch-preview-grid');
              if (status) {
                status.textContent = keepPreviews
                  ? 'Re-rendering volume previews with the new plot colours (ChimeraX)...'
                  : 'Rendering GIFs with ChimeraX...';
                status.style.color = '';
              }
              if (progress) {
                progress.hidden = false;
                progress.setAttribute('aria-busy', 'true');
              }
              if (grid && !keepPreviews) grid.innerHTML = '';
              function applyBadgeBackground(badge, colour) {
                badge.style.color = '#1a1a1a';
                badge.classList.remove('volsketch-gif-overlay--on-tint');
                if (typeof colour !== 'string' || !colour.trim()) {
                  badge.style.backgroundColor = '';
                  badge.style.borderColor = '';
                  return;
                }
                const m = /^#([0-9a-fA-F]{6})$/.exec(colour.trim());
                if (!m) return;
                const x = parseInt(m[1], 16);
                const rgba = 'rgba(' + ((x >> 16) & 255) + ',' + ((x >> 8) & 255)
                  + ',' + (x & 255) + ',0.88)';
                badge.style.backgroundColor = rgba;
                badge.style.borderColor = 'rgba(27, 31, 36, 0.28)';
                badge.classList.add('volsketch-gif-overlay--on-tint');
              }
              function mountOverlay(wrap, img, spec) {
                spec = spec || {};
                const badge = document.createElement('span');
                badge.className = 'volsketch-gif-overlay';
                const covBadge = document.createElement('span');
                covBadge.className = 'volsketch-gif-overlay volsketch-gif-overlay--cov';
                const covVar = document.createElement('span');
                covVar.className = 'volsketch-gif-cov-var';
                const covVal = document.createElement('span');
                covVal.className = 'volsketch-gif-cov-val';
                covBadge.appendChild(covVar);
                covBadge.appendChild(covVal);
                wrap.appendChild(badge);
                wrap.appendChild(covBadge);
                function hideCov() {
                  covVar.textContent = '';
                  covVal.textContent = '';
                  covVar.hidden = true;
                  covBadge.hidden = true;
                }
                function showCov(label, value, bg) {
                  covVar.textContent = label || '';
                  covVar.hidden = !covVar.textContent;
                  covVal.textContent = value || '';
                  covBadge.hidden = !covVal.textContent;
                  applyBadgeBackground(covBadge, bg || '');
                }
                function applySegment(idx) {
                  const labels = spec.segment_labels || [];
                  const bgs = spec.segment_backgrounds || [];
                  const covs = spec.segment_covariate_texts || [];
                  idx = Math.max(0, Math.min(labels.length - 1, idx | 0));
                  badge.textContent = labels[idx] || '';
                  badge.hidden = !badge.textContent;
                  applyBadgeBackground(badge, bgs[idx] || '');
                  if (covs.length && covs[idx]) {
                    showCov(spec.covariate_variable_label || '', String(covs[idx]), bgs[idx] || '');
                  } else {
                    hideCov();
                  }
                }
                if (spec.style === 'static') {
                  badge.textContent = spec.text || '';
                  badge.hidden = !badge.textContent;
                  applyBadgeBackground(badge, spec.badge_background || '');
                  if (spec.covariate_text) {
                    showCov(
                      spec.covariate_variable_label || '',
                      String(spec.covariate_text),
                      spec.covariate_badge_background || spec.badge_background || ''
                    );
                  } else {
                    hideCov();
                  }
                  return;
                }
                if (spec.style !== 'cycle_segments' || !(spec.segment_labels || []).length) {
                  badge.hidden = true;
                  hideCov();
                  return;
                }
                const fpv = Math.max(1, parseInt(spec.frames_per_segment, 10) || 8);
                const fd = Math.max(1, parseInt(spec.frame_duration_ms, 10) || 100);
                const totalF = Math.max(
                  1,
                  parseInt(spec.total_frames, 10) || (spec.segment_labels.length * fpv)
                );
                let raf = 0;
                let start = null;
                function tick(now) {
                  if (!wrap.isConnected) {
                    if (raf) cancelAnimationFrame(raf);
                    return;
                  }
                  if (start === null) start = now;
                  const frameIdx = Math.floor(((now - start) % (totalF * fd)) / fd);
                  applySegment(Math.floor(frameIdx / fpv));
                  raf = requestAnimationFrame(tick);
                }
                applySegment(0);
                const begin = () => { start = null; raf = requestAnimationFrame(tick); };
                if (img.decode) {
                  img.decode().then(begin).catch(begin);
                } else if (img.complete) {
                  begin();
                } else {
                  img.addEventListener('load', begin, { once: true });
                }
              }
              const payload = {
                vol_indices: vols,
                mode: 'cycle',
                chimerax_cpus: 8,
                cycle_frames_per_vol: 8,
                color_mode: colorSel ? colorSel.value : 'none',
              };
              fetch('/api/landscape_volpca/generate_animations', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(payload),
              })
              .then((r) => r.json().then((j) => {
                if (!r.ok) throw new Error(j.error || r.status);
                return j;
              }))
              .then((j) => {
                if (grid) {
                  grid.innerHTML = '';
                  grid.classList.add('volsketch-preview-grid--cycle');
                  (j.items || []).forEach((it) => {
                    const fig = document.createElement('figure');
                    const wrap = document.createElement('div');
                    wrap.className = 'volsketch-anim-wrap';
                    const img = document.createElement('img');
                    img.alt = '';
                    img.src = 'data:image/gif;base64,' + it.gif_b64;
                    wrap.appendChild(img);
                    mountOverlay(wrap, img, it.preview_overlay || {});
                    fig.appendChild(wrap);
                    grid.appendChild(fig);
                  });
                }
                if (status) {
                  const ds = j.duration_s != null ? Number(j.duration_s) : NaN;
                  status.textContent = Number.isNaN(ds)
                    ? 'Finished animating.'
                    : 'Finished animating in ' + Math.round(ds) + ' s.';
                }
              })
              .catch((e) => {
                window.__cdrgnLandscapeGifError = String(e);
                if (status) {
                  status.textContent = String(e);
                  status.style.color = 'var(--error, #b42318)';
                }
              })
              .finally(() => {
                if (progress) {
                  progress.hidden = true;
                  progress.removeAttribute('aria-busy');
                }
                window.__cdrgnLandscapeGifDone = true;
              });
            }""",
            {"vols": vols, "keepPreviews": bool(keep_previews)},
        )
        _snap_until_ready_capped(
            buf, page, _direct_render_done_js(), timeout_ms=300_000
        )
        page.wait_for_function(_direct_render_done_js(), timeout=300_000)
        err = page.evaluate("() => window.__cdrgnLandscapeGifError || ''")
        if err:
            raise RuntimeError(f"Landscape animation API failed: {err}")
        _wait_animation_ready(page, timeout_ms=30_000)


def _make_lasso_selection(buf: FrameBuffer, page, *, log) -> None:
    with _timed_step(log, "Draw lasso selection", slow_note=False):
        _force_cycle_mode(page)
        page.evaluate(
            """() => {
              const gd = document.getElementById('volsketch');
              if (gd && window.Plotly) Plotly.relayout(gd, { dragmode: 'lasso' });
            }"""
        )
        path = _lasso_path_from_plot(page)
        page.mouse.move(*path[0])
        page.mouse.down()
        for xy in path[1:]:
            page.mouse.move(*xy)
            time.sleep(buf.frame_ms / 1000.0)
            buf.snap(page)
        page.mouse.up()
        buf.sleep_snap(page, 0.25)

        if _selected_count(page) <= 0 and _ensure_visible_selection(page) <= 0:
            raise RuntimeError("The lasso did not select any sketch volumes.")
        _make_volume_points_visible(page)

    _render_cycle_preview_from_selection(buf, page, keep_previews=False, log=log)


def _capture_after_color_with_previews(
    buf: FrameBuffer, page, token: str, *, log
) -> str:
    chosen = _choose_color(buf, page, token, log=log)
    _render_cycle_preview_from_selection(buf, page, keep_previews=True, log=log)
    return chosen


def record_sequence(page, base: str, buf: FrameBuffer, *, log) -> None:
    with _timed_step(log, "Open landscape explorer directly", slow_note=True):
        page.goto(
            base + "/landscape-volpca", wait_until="domcontentloaded", timeout=120_000
        )
        _wait_scatter_ready(page, timeout_ms=180_000)
        _force_cycle_mode(page)
        _make_volume_points_visible(page)
    buf.sleep_snap(page, 0.35)

    _choose_color(buf, page, "state", log=log)
    buf.sleep_snap(page, 0.45)

    _make_lasso_selection(buf, page, log=log)
    buf.sleep_snap(page, 1.55)

    _capture_after_color_with_previews(buf, page, "znorm", log=log)
    buf.sleep_snap(page, 1.75)


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
        "--landscape-explorer",
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
        default=DEMO_ANIMATIONS_RECORDERS_DIR / "landscape_volpca_demo.gif",
    )
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument("--conda-env", type=str, default=None)
    parser.add_argument("--conda-prefix", type=Path, default=None)
    parser.add_argument("--cpus", type=int, default=8)
    parser.add_argument(
        "--frame-ms",
        type=int,
        default=FRAME_MS_DEFAULT,
        help=f"Frame duration in ms (default {FRAME_MS_DEFAULT}, about 20 fps)",
    )
    parser.add_argument("--max-width", type=int, default=1008)
    parser.add_argument("--headed", action="store_true")
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

    try:
        cmd = _dashboard_command(args, outdir, port)
    except RuntimeError as err:
        print(f"error: {err}", file=sys.stderr)
        return 1

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
        _wait_http(base + "/landscape-volpca", timeout_s=240.0, log=log)
        plotly_cached = _ensure_plotly_cached(log)
        buf = FrameBuffer(args.frame_ms)
        ph, launch_args = _chromium_launch_options(headed=args.headed)
        with sync_playwright() as p:
            with _timed_step(log, f"Launch Chromium (headless={ph})", slow_note=False):
                browser = p.chromium.launch(headless=ph, args=launch_args)
            context = browser.new_context(
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
            _route_landscape_scatter_as_svg(context)
            page = context.new_page()
            with _timed_step(log, "Full landscape capture", slow_note=True):
                record_sequence(page, base, buf, log=log)
            browser.close()

        nf, wall_s = _save_buf_to_gif(
            buf,
            args.output,
            frame_ms=args.frame_ms,
            max_width=args.max_width,
            max_wall_ms=MAX_GIF_WALL_MS,
            log=log,
            label="landscape GIF",
        )
        mean_ms = (wall_s * 1000.0 / nf) if nf else 0.0
        print(
            f"Wrote {args.output} ({nf} frames, approx {wall_s:.2f} s playback, "
            f"mean {mean_ms:.1f} ms/frame)"
        )
        return 0
    except Exception as err:
        print(f"error: {err}", file=sys.stderr)
        return 1
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=15)
        except subprocess.TimeoutExpired:
            proc.kill()


if __name__ == "__main__":
    raise SystemExit(main())
