#!/usr/bin/env python3
"""Headless GIF: particle explorer (lasso → cache → k-means → legend toggle → ``znorm`` histogram).

Encodes ``duration`` × ``fps`` frames at ~``fps`` (default **16** s × **40** ⇒ **640** frames). GIF delay is centisecond-quantized, so
per-frame delay alternates (e.g. **20** / **30** ms for 40 fps) so the average rate matches ``fps``
and total playback matches ``duration``.

Optional ``--annotate-caption`` composites a beige caption band under each frame (consistent height).
With ``--annotate-caption``, each contiguous run of the same caption is extended to at least
``--annotate-min-segment-seconds`` wall time (default **3.5** s); the GIF gets longer accordingly.
Use ``--high-resolution-output PATH`` plus ``--high-resolution-max-width 0`` to keep native browser
resolution (typically 1440 px wide screenshots) alongside the resized ``--output``.

Depends on Pillow, Playwright/Chromium — same assumptions as ``record_dashboard_interactions_gif``.

Example::

    PYTHONPATH=/path/to/cryodrgn_beta \\
      conda run -p /projects/CRYOEM/zhonglab/mg2332/conda_envs/cdrgn_beta \\
        --no-capture-output \\
        python -m cryodrgn.dashboard.record_particle_explorer_actions_gif \\
        /scratch/gpfs/.../004_train-vae_1gpu_dim.1024/ \\
        -o particle_explorer_actions.gif \\
        --annotate-caption \\
        --high-resolution-output particle_explorer_actions_hr.gif \\
        --high-resolution-max-width 0
"""

from __future__ import annotations

import argparse
import contextlib
import io
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from cryodrgn.dashboard.record_dashboard_interactions_gif import (
    DEMO_ANIMATIONS_DIR,
    PLOTLY_CDN_ROUTE_GLOB,
    REPO_ROOT,
    RENDER_POLL_FRAME_DURATION_MULT,
    RENDER_POLL_SNAPSHOT_STRIDE,
    FrameBuffer,
    _chromium_launch_options,
    _ensure_plotly_cached,
    _free_port,
    _montage_cache_ready_js,
    _scatter_plot_box,
    _timed_step,
    _wait_http,
    _wait_scatter_ready,
)

DEFAULT_CONDA_PREFIX = Path("/projects/CRYOEM/zhonglab/mg2332/conda_envs/cdrgn_beta")
DEFAULT_OUTPUT = DEMO_ANIMATIONS_DIR / "particle_explorer_actions_demo.gif"
_LOG_PREFIX = "[record-particle-explorer-actions-gif]"


class CaptionedExplorerBuffer(FrameBuffer):
    """Like :class:`FrameBuffer` but records the current ``step_label`` for each frame."""

    def __init__(self, frame_ms: int) -> None:
        super().__init__(frame_ms)
        self.step_label = ""
        self.captions_per_frame: list[str] = []

    def snap(self, page, *, duration_ms: int | None = None) -> None:
        self.captions_per_frame.append(self.step_label.strip())
        super().snap(page, duration_ms=duration_ms)

    def add_png(self, data: bytes, *, duration_ms: int | None = None) -> None:
        self.captions_per_frame.append(self.step_label.strip())
        super().add_png(data, duration_ms=duration_ms)


def _stamp_step(buf: FrameBuffer, caption: str) -> None:
    if hasattr(buf, "step_label"):
        setattr(buf, "step_label", caption)


def _pil_ui_font(size_px: int):
    from PIL import ImageFont

    candidates = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
        "/usr/share/fonts/truetype/freefont/FreeSans.ttf",
    ]
    for path in candidates:
        try:
            return ImageFont.truetype(path, int(size_px))
        except OSError:
            continue
    return ImageFont.load_default()


def _wrap_words_to_fit(
    words: list[str],
    *,
    joiner: str,
    font,
    draw,
    max_px: int,
) -> list[str]:
    lines: list[str] = []
    cur: list[str] = []
    for w in words:
        cand = joiner.join(cur + [w]) if cur else w
        bbox = draw.textbbox((0, 0), cand, font=font)
        if bbox[2] - bbox[0] <= max_px or not cur:
            cur.append(w)
        else:
            lines.append(joiner.join(cur))
            cur = [w]
    if cur:
        lines.append(joiner.join(cur))
    return lines


def _caption_lines_wrapped(text: str, *, font, draw, max_line_px: int) -> list[str]:
    paragraphs = [p.strip() for p in text.replace("\r", "").split("\n") if p.strip()]
    if not paragraphs:
        return []
    lines: list[str] = []
    for pi, pg in enumerate(paragraphs):
        lines.extend(
            _wrap_words_to_fit(
                pg.split(), joiner=" ", font=font, draw=draw, max_px=max_line_px
            )
        )
        if pi < len(paragraphs) - 1:
            lines.append("")
    return lines


def _max_footer_height_for_captions(
    captions_seq: list[str],
    *,
    content_width_px: int,
    font_px: int,
    h_pad: int = 26,
    line_gap: int = 4,
    cap_max_px: int = 200,
) -> int:
    from PIL import ImageDraw, Image

    if not captions_seq:
        return 0
    font = _pil_ui_font(font_px)
    im = Image.new("RGB", (max(240, content_width_px), 32))
    draw = ImageDraw.Draw(im)
    max_line_px = max(120, content_width_px - 32)
    need = 8
    for cap in captions_seq:
        if not cap.strip():
            need = max(need, h_pad)
            continue
        ls = _caption_lines_wrapped(cap, font=font, draw=draw, max_line_px=max_line_px)
        h = h_pad
        for ln in ls:
            bbox = draw.textbbox((0, 0), ln or " ", font=font)
            h += (bbox[3] - bbox[1]) + line_gap
        need = max(need, min(cap_max_px, h))
    return int(need)


def _composite_fixed_footer_below(
    im_rgb,
    *,
    caption: str,
    footer_h: int,
    font_px: int = 21,
):
    """Return a taller RGB image with a caption band (`footer_h`) under ``im_rgb``."""
    from PIL import ImageDraw, Image

    w, h = im_rgb.size
    if footer_h <= 0:
        return im_rgb
    top_rgb = (250, 248, 244)
    foot_rgb = (245, 243, 238)
    border_rgb = (180, 192, 202)
    text_rgb = (36, 59, 83)

    canvas = Image.new("RGB", (w, h + footer_h), top_rgb)
    canvas.paste(im_rgb, (0, 0))
    draw = ImageDraw.Draw(canvas)
    draw.rectangle([0, h, w, h + footer_h], fill=foot_rgb)
    draw.rectangle([0, h - 1, w, h + 1], outline=border_rgb, fill=border_rgb)

    font = _pil_ui_font(font_px)
    inner_w = max(40, w - 28)
    cap = caption.strip()
    if not cap:
        return canvas
    lines = _caption_lines_wrapped(cap, font=font, draw=draw, max_line_px=inner_w)
    cx = 14
    y = h + 10
    for ln in lines:
        draw.text((cx, y), ln, fill=text_rgb, font=font)
        bbox = draw.textbbox((cx, y), ln or " ", font=font)
        y += (bbox[3] - bbox[1]) + (5 if ln else 2)

    return canvas


def _make_pe_logger(quiet: bool):
    def log(msg: str) -> None:
        if not quiet:
            print(f"{_LOG_PREFIX} {msg}", flush=True)

    return log


def _wait_preload_overlay_clear(page, timeout_ms: float = 360_000) -> None:
    page.wait_for_function(
        """() => {
          const o = document.getElementById('montage-preload-overlay');
          return !o || !o.classList.contains('cryo-plot-rendering-overlay--show');
        }""",
        timeout=int(timeout_ms),
    )


def _znorm_option_value(page) -> str | None:
    """Return the ``#sc`` option value for ``znorm`` (exact value or label match), else ``None``."""
    return page.evaluate(
        """() => {
          const sel = document.getElementById('sc');
          if (!sel) return null;
          const opts = sel.querySelectorAll('option');
          for (let i = 0; i < opts.length; i++) {
            const o = opts[i];
            if (!o.value || o.value === 'none') continue;
            if (o.value === 'znorm') return o.value;
          }
          for (let j = 0; j < opts.length; j++) {
            const o = opts[j];
            const t = ('' + (o.textContent || '') + ' ' + (o.value || '')).trim();
            if (/^znorm\\b/i.test(o.value || '') || /^znorm\\b/i.test(t)) return o.value || null;
          }
          return null;
        }"""
    )


def _kmeans_covariate_value(page) -> str | None:
    """Return `<option value=…>` for a k-means column, else ``None``."""
    return page.evaluate(
        """() => {
          const sel = document.getElementById('sc');
          if (!sel) return null;
          const opts = sel.querySelectorAll('option');
          for (let i = 0; i < opts.length; i++) {
            const o = opts[i];
            const t = (o.textContent || '') + ' ' + (o.value || '');
            if (/k[-\\s_]*means|kmeans|analyze[^\\s\\/]*kmeans/i.test(t)) return o.value || null;
          }
          return null;
        }"""
    )


def _densify_lasso_rel(
    rel: list[tuple[float, float]], subdivisions_between: int
) -> list[tuple[float, float]]:
    """Insert linear steps between consecutive waypoints so the lasso traces more slowly."""
    if subdivisions_between < 2 or len(rel) < 2:
        return list(rel)
    out: list[tuple[float, float]] = []
    for i in range(len(rel) - 1):
        ax, ay = rel[i]
        bx, by = rel[i + 1]
        for k in range(subdivisions_between):
            t = k / subdivisions_between
            out.append((ax + (bx - ax) * t, ay + (by - ay) * t))
    out.append(rel[-1])
    return out


def _snap_montage_load(
    buf: FrameBuffer,
    page,
    *,
    frame_ms: int,
    max_polls: int = 400,
) -> None:
    render_ms = max(1, int(frame_ms * RENDER_POLL_FRAME_DURATION_MULT))
    for step in range(max_polls):
        if page.evaluate(_montage_cache_ready_js()):
            break
        if step % RENDER_POLL_SNAPSHOT_STRIDE == 0:
            buf.snap(page, duration_ms=render_ms)
        time.sleep(frame_ms / 1000.0)


def normalize_frames_fixed_fps(
    pngs: list[bytes],
    *,
    fps: float,
    wall_seconds: float,
    captions: list[str] | None = None,
) -> tuple[list[bytes], list[int], list[str] | None]:
    if not pngs:
        raise RuntimeError("no frames captured")
    if captions is not None and len(captions) != len(pngs):
        raise RuntimeError("captions length must match pngs")
    target_n = max(1, int(round(float(wall_seconds) * float(fps))))
    # GIF delay is centiseconds; ~25 ms (40 fps) must be approximated (e.g. 20 + 30 ms).
    target_ms = 1000.0 / float(fps)
    short_cs = max(1, int(target_ms // 10))
    long_ms = (short_cs + 1) * 10
    short_ms = short_cs * 10
    durs = [(long_ms if i % 2 else short_ms) for i in range(target_n)]
    expected_ms = float(wall_seconds) * 1000.0
    drift = int(round(expected_ms - float(sum(durs))))
    if drift != 0 and durs:
        durs[-1] = max(10, durs[-1] + drift)
    n = len(pngs)
    out_caps: list[str] | None = None
    if n == target_n:
        out_pngs = list(pngs)
        if captions is not None:
            out_caps = list(captions)
    elif n > target_n:
        if target_n == 1:
            out_pngs = [pngs[-1]]
            if captions is not None:
                out_caps = [captions[-1]]
        else:
            idxs = [int(round(i * (n - 1) / (target_n - 1))) for i in range(target_n)]
            out_pngs = [pngs[j] for j in idxs]
            if captions is not None:
                out_caps = [captions[j] for j in idxs]
    else:
        out_pngs = list(pngs) + [pngs[-1]] * (target_n - n)
        if captions is not None:
            out_caps = list(captions) + [captions[-1]] * (target_n - n)
    assert len(out_pngs) == target_n
    assert len(durs) == target_n
    if out_caps is not None:
        assert len(out_caps) == target_n
    return out_pngs, durs, out_caps


def extend_caption_segments_min_duration(
    pngs: list[bytes],
    durations_ms: list[int],
    captions: list[str],
    *,
    min_segment_s: float,
    fps: float,
) -> tuple[list[bytes], list[int], list[str], int]:
    """For each contiguous block of identical caption text, ensure total display time ≥ ``min_segment_s``.

    Appends duplicate PNGs (with short delays) at the end of under-long segments. Empty captions
    are left unchanged. Returns new lists and how many non-empty segments were extended.
    """

    n = len(pngs)
    if n == 0 or min_segment_s <= 0:
        return pngs, durations_ms, captions, 0
    if len(durations_ms) != n or len(captions) != n:
        raise RuntimeError("pngs, durations_ms, and captions must align")
    min_ms = float(min_segment_s) * 1000.0
    per = max(10, int(round(1000.0 / max(float(fps), 1e-6))))

    out_p: list[bytes] = []
    out_d: list[int] = []
    out_c: list[str] = []
    n_extended = 0
    i = 0
    while i < n:
        j = i + 1
        while j < n and captions[j] == captions[i]:
            j += 1
        block_p = list(pngs[i:j])
        block_d = list(durations_ms[i:j])
        cap_key = captions[i]
        cap_strip = cap_key.strip()
        total_ms = float(sum(block_d))
        if cap_strip and total_ms < min_ms:
            n_extended += 1
            need_ms = min_ms - total_ms
            last_png = block_p[-1]
            while need_ms > 0.5:
                d = max(10, min(per, int(round(need_ms))))
                block_p.append(last_png)
                block_d.append(d)
                need_ms -= d
        out_p.extend(block_p)
        out_d.extend(block_d)
        out_c.extend([cap_key] * len(block_p))
        assert len(block_p) == len(block_d)
        i = j
    assert len(out_p) == len(out_d) == len(out_c)
    return out_p, out_d, out_c, n_extended


def save_pngs_duration_gif(
    pngs: list[bytes],
    durations_ms: list[int],
    output: Path,
    *,
    max_width: int,
    log,
    annotate: bool = False,
    step_captions: list[str] | None = None,
    caption_font_px: int = 21,
    caption_band_max_height: int = 200,
) -> None:
    from PIL import Image

    with _timed_step(
        log,
        f"Encode GIF: {len(pngs)} frames → {output.name}",
        slow_note=True,
    ):
        if annotate and step_captions is not None and len(step_captions) == len(pngs):
            probe = Image.open(io.BytesIO(pngs[0]))
            content_w = probe.width
            probe.close()
            footer_h = _max_footer_height_for_captions(
                step_captions,
                content_width_px=content_w,
                font_px=caption_font_px,
                cap_max_px=caption_band_max_height,
            )
        else:
            footer_h = 0

        pil_frames: list[Image.Image] = []
        for i, raw in enumerate(pngs):
            im = Image.open(io.BytesIO(raw)).convert("RGB")
            if annotate and step_captions is not None and footer_h > 0:
                im = _composite_fixed_footer_below(
                    im,
                    caption=step_captions[i] if i < len(step_captions) else "",
                    footer_h=footer_h,
                    font_px=caption_font_px,
                )
            if max_width and im.width > max_width:
                ratio = max_width / im.width
                im = im.resize(
                    (max_width, max(1, int(im.height * ratio))),
                    Image.Resampling.LANCZOS,
                )
            # Pillow drops byte-identical frames when building animated GIFs; nudge one
            # corner pixel per frame index so the full fps × duration frame count is kept.
            px = im.load()
            w, h = im.size
            if w > 0 and h > 0:
                u = max(1, min(i + 1, (1 << 24) - 1))
                px[w - 1, h - 1] = (
                    u & 0xFF,
                    (u >> 8) & 0xFF,
                    (u >> 16) & 0xFF,
                )
            pil_frames.append(im)
        output.parent.mkdir(parents=True, exist_ok=True)
        pil_frames[0].save(
            output,
            save_all=True,
            append_images=pil_frames[1:],
            duration=durations_ms,
            loop=0,
            optimize=False,
        )
        for im in pil_frames:
            with contextlib.suppress(Exception):
                im.close()


def record_sequence(
    page,
    base: str,
    buf: FrameBuffer,
    *,
    frame_ms: int,
    log,
) -> None:
    """Capture explorer actions; every ``buf.snap`` uses ``frame_ms`` display time."""

    page.goto(base + "/explorer", wait_until="domcontentloaded", timeout=180_000)
    with _timed_step(log, "Wait for scatter render", slow_note=True):
        _wait_scatter_ready(page, 300_000)

    _stamp_step(
        buf,
        "Latent embedding scatter — particle thumbnails load as you browse the map.",
    )
    for _ in range(max(1, int(round(0.35 * 1000 / frame_ms)))):
        buf.snap(page)
        time.sleep(frame_ms / 1000.0)

    box = _scatter_plot_box(page)
    if not box:
        raise RuntimeError("could not find #scatter bounding box")
    w, h = box["width"], box["height"]
    page.evaluate(
        """() => {
          const gd = document.getElementById('scatter');
          if (!gd || !window.Plotly) return false;
          Plotly.relayout(gd, { dragmode: 'lasso' });
          return true;
        }"""
    )
    rel = [
        (0.32, 0.48),
        (0.44, 0.38),
        (0.58, 0.42),
        (0.64, 0.55),
        (0.52, 0.64),
        (0.38, 0.58),
        (0.30, 0.50),
    ]
    rel_dense = _densify_lasso_rel(rel, subdivisions_between=10)
    lasso_step_s = max(frame_ms / 1000.0 * 2.4, 0.055)
    _stamp_step(
        buf,
        "Lasso tool — tracing a polygon to grab particles in one neighbourhood.",
    )
    page.mouse.move(box["x"] + w * rel_dense[0][0], box["y"] + h * rel_dense[0][1])
    for _ in range(max(1, int(round(0.2 * 1000 / frame_ms)))):
        buf.snap(page)
        time.sleep(frame_ms / 1000.0)
    page.mouse.down()
    buf.snap(page)
    time.sleep(lasso_step_s)
    for i in range(1, len(rel_dense)):
        page.mouse.move(box["x"] + w * rel_dense[i][0], box["y"] + h * rel_dense[i][1])
        time.sleep(lasso_step_s)
        buf.snap(page)
    page.mouse.up()
    buf.snap(page)

    for _ in range(max(1, int(round(0.45 * 1000 / frame_ms)))):
        buf.snap(page)
        time.sleep(frame_ms / 1000.0)

    _stamp_step(
        buf,
        "Preload panel — widen the sliding window so more thumbnails stay hot in RAM.",
    )
    with _timed_step(log, "Build image cache (first chunk)", slow_note=True):
        expand = page.locator("#btn-expand-cache")
        if expand.count() and expand.is_enabled():
            expand.click()
            try:
                _wait_preload_overlay_clear(page, 360_000)
            except Exception:
                pass
        for _ in range(max(1, int(round(0.5 * 1000 / frame_ms)))):
            buf.snap(page)
            time.sleep(frame_ms / 1000.0)

    _stamp_step(
        buf,
        "Selection preload — enqueue uncached picks from your lasso into the window.",
    )
    with _timed_step(
        log, "Add selection images to cache (if available)", slow_note=True
    ):
        sel_uncached = page.locator("#btn-cache-selection-uncached")
        if sel_uncached.count() and sel_uncached.is_enabled():
            sel_uncached.click()
            try:
                _wait_preload_overlay_clear(page, 360_000)
            except Exception:
                pass
        for _ in range(max(1, int(round(0.35 * 1000 / frame_ms)))):
            buf.snap(page)
            time.sleep(frame_ms / 1000.0)

    with _timed_step(log, "Open image grid from cache", slow_note=True):
        view_btn = page.locator("#btn-view-images")
        if view_btn.count() and view_btn.is_enabled():
            _stamp_step(
                buf,
                "Image montage — open a grid of thumbnails from what's cached.",
            )
            view_btn.click()
            try:
                page.wait_for_function(_montage_cache_ready_js(), timeout=180_000)
            except Exception:
                pass
            _snap_montage_load(buf, page, frame_ms=frame_ms)
            for _ in range(max(1, int(round(0.85 * 1000 / frame_ms)))):
                buf.snap(page)
                time.sleep(frame_ms / 1000.0)
            _stamp_step(
                buf,
                "Grid layout menu — tweak how many images appear per row and column.",
            )
            gs = page.locator("#grid-size")
            toggle = page.locator("#image-grid-menu-toggle")
            if toggle.count():
                aria = toggle.get_attribute("aria-expanded") or "false"
                if aria.strip().lower() == "false":
                    toggle.click()
                    for _ in range(max(1, int(round(0.35 * 1000 / frame_ms)))):
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
                for _ in range(max(1, int(round(0.45 * 1000 / frame_ms)))):
                    buf.snap(page)
                    time.sleep(frame_ms / 1000.0)
                _snap_montage_load(buf, page, frame_ms=frame_ms, max_polls=80)
                for _ in range(max(1, int(round(1.1 * 1000 / frame_ms)))):
                    buf.snap(page)
                    time.sleep(frame_ms / 1000.0)
        else:
            for _ in range(max(1, int(round(0.6 * 1000 / frame_ms)))):
                buf.snap(page)
                time.sleep(frame_ms / 1000.0)

    km = _kmeans_covariate_value(page)
    if km:
        _stamp_step(
            buf,
            "Discrete colours — shade the scatter by k-means (or similar) clusters.",
        )
        with _timed_step(log, "Colour by k-means covariate", slow_note=False):
            page.locator("#sc").select_option(value=str(km))
            _wait_scatter_ready(page, 300_000)
            for _ in range(max(1, int(round(0.85 * 1000 / frame_ms)))):
                buf.snap(page)
                time.sleep(frame_ms / 1000.0)
    else:
        log("No k-means-like column in colour dropdown; skipping recolour step.")

    discrete_ready = False
    try:
        page.wait_for_function(
            """() => {
              const inputs = document.querySelectorAll('#color-discrete-switches input[type=\"checkbox\"]');
              return inputs.length >= 2;
            }""",
            timeout=60_000,
        )
        discrete_ready = True
    except Exception:
        log("Discrete legend switches never appeared (need k-means-style column?).")

    if discrete_ready:
        _stamp_step(
            buf,
            "Legend switches — mute or highlight individual clusters without redrawing.",
        )
        with _timed_step(log, "Toggle first legend cluster switch", slow_note=False):
            page.locator(
                "#color-discrete-switches label.cryo-color-discrete-switch"
            ).first.click(force=True)
            try:
                _wait_scatter_ready(page, 120_000)
            except Exception:
                pass
            for _ in range(max(1, int(round(1.2 * 1000 / frame_ms)))):
                buf.snap(page)
                time.sleep(frame_ms / 1000.0)

    zn = _znorm_option_value(page)
    if zn:
        _stamp_step(
            buf,
            "Continuous covariate + histogram — z-scores and click-drag on the histogram to threshold.",
        )
        with _timed_step(
            log, "Colour by znorm + histogram threshold selection", slow_note=True
        ):
            page.locator("#sc").select_option(value=str(zn))
            _wait_scatter_ready(page, 300_000)
            page.wait_for_function(
                """() => {
                  const p = document.getElementById('color-hist-panel');
                  const gd = document.getElementById('color-hist');
                  const ax = gd && gd._fullLayout && gd._fullLayout.xaxis;
                  return !!(
                    p && !p.hidden && gd && gd.data && gd.data.length
                    && ax && typeof ax.p2c === 'function'
                  );
                }""",
                timeout=120_000,
            )
            for _ in range(max(1, int(round(0.55 * 1000 / frame_ms)))):
                buf.snap(page)
                time.sleep(frame_ms / 1000.0)
            hist = page.locator("#color-hist")
            hist.scroll_into_view_if_needed()
            frac_seq = (0.30, 0.52, 0.74, 0.42)
            for fi, frac_x in enumerate(frac_seq):
                box = hist.bounding_box()
                if box:
                    click_x = max(6.0, min(box["width"] - 6.0, box["width"] * frac_x))
                    click_y = max(6.0, min(box["height"] - 6.0, box["height"] * 0.45))
                    hist.click(position={"x": click_x, "y": click_y})
                try:
                    _wait_scatter_ready(page, 180_000)
                except Exception:
                    pass
                hold_ratio = 0.85 if fi < len(frac_seq) - 1 else 1.05
                for _ in range(max(1, int(round(hold_ratio * 1000 / frame_ms)))):
                    buf.snap(page)
                    time.sleep(frame_ms / 1000.0)
    else:
        log("No znorm covariate in colour dropdown; skipping histogram step.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "outdir",
        type=Path,
        help="experiment output dir passed to cryodrgn dashboard",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
    )
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument(
        "--conda-prefix",
        type=Path,
        default=DEFAULT_CONDA_PREFIX,
        help="conda env prefix for conda run -p … (default cdrgn_beta path on zhonglab fs)",
    )
    parser.add_argument(
        "--fps",
        type=float,
        default=40.0,
        help="output GIF playback frame rate (default 40)",
    )
    parser.add_argument(
        "--duration",
        type=float,
        default=16.0,
        help=(
            "encoded playback duration in seconds (default is ~60 percent longer "
            "than the original 10 s demo so lasso/grid/histogram have more room)"
        ),
    )
    parser.add_argument("--filter-max", type=int, default=40_000)
    parser.add_argument("--cpus", type=int, default=8)
    parser.add_argument("--max-width", type=int, default=512)
    parser.add_argument(
        "--annotate-caption",
        action="store_true",
        help="draw a caption band below each frame (uses step text recorded during recording)",
    )
    parser.add_argument(
        "--high-resolution-output",
        type=Path,
        default=None,
        metavar="PATH",
        help=(
            "write a second GIF (e.g. full-width); use "
            "`--high-resolution-max-width 0` to skip downscaling after capture"
        ),
    )
    parser.add_argument(
        "--high-resolution-max-width",
        type=int,
        default=0,
        metavar="PX",
        help="max width for the high-resolution file (default 0 = native screenshot width)",
    )
    parser.add_argument(
        "--caption-font-size",
        type=int,
        default=21,
        metavar="PX",
        help="pixel size for caption band text when --annotate-caption is set (default 21)",
    )
    parser.add_argument(
        "--annotate-min-segment-seconds",
        type=float,
        default=3.5,
        metavar="SEC",
        help=(
            "when using --annotate-caption, each distinct on-screen caption span "
            "is held at least this long (default 3.5); lengthens GIF if needed"
        ),
    )
    parser.add_argument(
        "--headed",
        action="store_true",
        help="headed Chromium",
    )
    parser.add_argument("-q", "--quiet", action="store_true")
    args = parser.parse_args()

    conda_prefix = Path(args.conda_prefix)
    fps = float(args.fps)
    wall_s = float(args.duration)
    if fps <= 0 or wall_s <= 0:
        print("error: fps and duration must be positive", file=sys.stderr)
        return 1
    frame_ms = max(1, int(round(1000.0 / fps)))
    log = _make_pe_logger(quiet=args.quiet)

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
            with _timed_step(log, "Record particle explorer clip", slow_note=True):
                record_sequence(page, base, buf, frame_ms=frame_ms, log=log)
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
