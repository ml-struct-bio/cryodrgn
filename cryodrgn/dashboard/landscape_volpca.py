"""Volume-space PCA scatter and ChimeraX animations for ``analyze_landscape`` outputs."""

from __future__ import annotations

import atexit
import base64
import os
import random
import re
import secrets
import shutil
import tempfile
import threading
import time
from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.colors import sample_colorscale

from cryodrgn import utils
from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.explorer_volumes import (
    DEFAULT_CHIMERAX_PARALLEL,
    DEFAULT_GIF_FRAMES,
    montage_cell_label,
    mrc_to_rotating_gif,
    mrc_to_static_png,
)
from cryodrgn.dashboard.plots import (
    _DASHBOARD_CREAM,
    _PLOTLY_FONT,
    _labels_colors_and_legend_items,
    _lower_legend_entry_label,
    _plotly_to_json,
    normalize_continuous_palette,
)

_LANDSCAPE_ANIM_LOCK = threading.Lock()
_LANDSCAPE_ANIM_ROOT: str | None = None
# token -> job_dir, files, landscape_dir, kmeans_k, t0, mode, color_mode,
# optional rotate_keyframes: { vol_index: png_path } (first frame of each rotate GIF)
_LANDSCAPE_ANIM_ENTRIES: dict[str, dict[str, Any]] = {}

_LANDSCAPE_VOLPCA_KMEANS_RE = re.compile(r"^kmeans(\d+)$")
_LANDSCAPE_VOLPCA_PKL_RE = re.compile(r"^vol_pca_(\d+)\.pkl$")
# ``vol_mean.mrc`` and similar match ``vol_*.mrc`` but are not k-means centroids.
_KMEANS_SKETCH_VOL_MRC_RE = re.compile(r"^vol_(\d+)\.mrc$")

# ChimeraX cost: cap how many volumes go into each animation style (extras subsampled).
LANDSCAPE_ANIM_MAX_ROTATE = 5
LANDSCAPE_ANIM_MAX_CYCLE = 50


def _sample_landscape_vols(vols: list[int], k: int, rng: random.Random) -> list[int]:
    """Up to ``k`` volumes without replacement; sorted ascending."""
    if len(vols) <= k:
        return list(vols)
    return sorted(rng.sample(vols, k))


def list_landscape_epochs(workdir: str) -> list[int]:
    """Epochs ``N`` with a directory ``landscape.N`` under ``workdir``."""
    if not os.path.isdir(workdir):
        return []
    out: list[int] = []
    with os.scandir(workdir) as it:
        for entry in it:
            if not entry.is_dir():
                continue
            m = re.fullmatch(r"landscape\.([0-9]+)", entry.name)
            if m:
                out.append(int(m.group(1)))
    return sorted(out)


def landscape_dir_for_epoch(workdir: str, epoch: int) -> str:
    """Path to ``landscape.{epoch}`` under a training workdir."""
    return os.path.join(workdir, f"landscape.{int(epoch)}")


def resolve_kmeans_sketch_bundle(landscape_dir: str) -> tuple[int, str]:
    """Return ``(K, kmeans_dir)`` where ``vol_pca_K.pkl`` matches ``kmeansK``."""
    if not os.path.isdir(landscape_dir):
        raise FileNotFoundError(landscape_dir)
    for name in sorted(os.listdir(landscape_dir)):
        m = _LANDSCAPE_VOLPCA_PKL_RE.fullmatch(name)
        if not m:
            continue
        k = int(m.group(1))
        kd = os.path.join(landscape_dir, f"kmeans{k}")
        if os.path.isdir(kd):
            pkl = os.path.join(landscape_dir, f"vol_pca_{k}.pkl")
            if os.path.isfile(pkl):
                return k, kd
    raise FileNotFoundError(
        f"No vol_pca_<K>.pkl + kmeans<K> pair under {landscape_dir!r}."
    )


def landscape_analysis_ready(workdir: str, epoch: int) -> bool:
    """True when ``landscape.{epoch}`` exists and has a matching ``vol_pca_*`` + ``kmeans*`` bundle."""
    try:
        resolve_kmeans_sketch_bundle(landscape_dir_for_epoch(workdir, int(epoch)))
    except FileNotFoundError:
        return False
    return True


def kmeans_sorted_vol_indices(kmeans_dir: str) -> list[int]:
    """Sorted numeric indices for ``vol_NNN.mrc`` sketch maps only (excludes ``vol_mean.mrc``, etc.)."""
    if not os.path.isdir(kmeans_dir):
        return []
    idx: list[int] = []
    with os.scandir(kmeans_dir) as it:
        for entry in it:
            if not entry.is_file():
                continue
            m = _KMEANS_SKETCH_VOL_MRC_RE.fullmatch(entry.name)
            if m:
                idx.append(int(m.group(1)))
    return sorted(idx)


def load_vol_pca_matrix(landscape_dir: str, k_sketch: int) -> np.ndarray:
    """2-D PCA coordinates matrix from ``vol_pca_{k_sketch}.pkl``."""
    pkl = os.path.join(landscape_dir, f"vol_pca_{k_sketch}.pkl")
    pc = np.asarray(utils.load_pkl(pkl), dtype=np.float64)
    if pc.ndim != 2:
        raise ValueError(f"Expected 2-D vol PCA array in {pkl}")
    return pc


def load_pca_explained_variance(landscape_dir: str) -> np.ndarray | None:
    obj_path = os.path.join(landscape_dir, "vol_pca_obj.pkl")
    if not os.path.isfile(obj_path):
        return None
    pca = utils.load_pkl(obj_path)
    evr = getattr(pca, "explained_variance_ratio_", None)
    if evr is None:
        return None
    return np.asarray(evr, dtype=np.float64)


def _first_sketch_clustering_dir(landscape_dir: str) -> str | None:
    cands = [
        d
        for d in os.listdir(landscape_dir)
        if d.startswith("sketch_clustering_")
        and os.path.isdir(os.path.join(landscape_dir, d))
    ]
    if not cands:
        return None
    cands.sort()
    return os.path.join(landscape_dir, cands[0])


def load_sketch_state_labels(landscape_dir: str) -> np.ndarray | None:
    sub = _first_sketch_clustering_dir(landscape_dir)
    if not sub:
        return None
    lab_path = os.path.join(sub, "state_labels.pkl")
    if not os.path.isfile(lab_path):
        return None
    return np.asarray(utils.load_pkl(lab_path), dtype=np.int64)


def landscape_vol_state_hex_by_vol_index(
    landscape_dir: str,
    kmeans_dir: str,
    k_sketch: int,
) -> dict[int, str] | None:
    """Vol index (``11`` for ``vol_011.mrc``) → fill hex, matching state scatter colours."""
    states = load_sketch_state_labels(landscape_dir)
    if states is None:
        return None
    pc = load_vol_pca_matrix(landscape_dir, k_sketch)
    n = int(pc.shape[0])
    if len(states) != n:
        return None
    sorted_vols = kmeans_sorted_vol_indices(kmeans_dir)
    vol_ids = sorted_vols[:n] if len(sorted_vols) >= n else sorted_vols
    if len(vol_ids) < n:
        return None
    ser = pd.Series(states)
    colors, _ = _labels_colors_and_legend_items(ser)
    return {int(vol_ids[i]): str(colors[i]) for i in range(n)}


def vol_mrc_path(kmeans_dir: str, vol_index: int) -> str:
    p = os.path.join(kmeans_dir, f"vol_{int(vol_index):03d}.mrc")
    if not os.path.isfile(p):
        raise FileNotFoundError(
            f"Sketched k-means volume not found: {p}",
        )
    return p


def _volsketch_covariate_display(name: str) -> str:
    """Match particle explorer naming for covariate selectors."""
    if name == "labels":
        return "k-means labels"
    return str(name)


def landscape_color_options(exp: DashboardExperiment) -> list[dict[str, str]]:
    opts: list[dict[str, str]] = [
        {"value": "none", "label": "None"},
        {"value": "state", "label": "Agglomerative state"},
    ]
    for c in exp.numeric_columns:
        opts.append({"value": c, "label": _volsketch_covariate_display(c)})
    return opts


def load_sketch_centroid_plot_df_rows(kmeans_dir: str, n_expected: int) -> np.ndarray:
    path = os.path.join(kmeans_dir, "centers_ind.txt")
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Missing {path} (expected from analyze_landscape k-means export)."
        )
    ind = np.loadtxt(path)
    ind = np.atleast_1d(ind).astype(np.int64, copy=False)
    if ind.shape[0] != n_expected:
        raise ValueError(
            f"centers_ind.txt has {ind.shape[0]} rows but vol_pca has {n_expected}."
        )
    return ind


def _numeric_array_to_plotly_hex(vals: np.ndarray, plotly_cs: str) -> list[str]:
    """Map finite floats through a Plotly colorscale name; NaN → neutral gray."""
    vals = np.asarray(vals, dtype=float)
    n = int(vals.shape[0])
    finite = np.isfinite(vals)
    if not finite.any():
        return ["#9ca3af"] * n
    vmin = float(np.nanmin(vals))
    vmax = float(np.nanmax(vals))
    span = vmax - vmin
    out: list[str] = []
    for i in range(n):
        if not finite[i]:
            out.append("#9ca3af")
            continue
        t = 0.5 if span <= 0 else (float(vals[i]) - vmin) / span
        t = max(0.0, min(1.0, t))
        out.append(sample_colorscale(plotly_cs, [t])[0])
    return out


def sketch_vol_marker_hex_by_vol_index(
    landscape_dir: str,
    exp: DashboardExperiment,
    *,
    color_mode: str,
    continuous_palette: str | None = None,
) -> dict[int, str]:
    """Sketch volume id → marker fill hex (same logic as the scatter); used for preview badge fills."""
    k_sketch, kmeans_dir = resolve_kmeans_sketch_bundle(landscape_dir)
    pc = load_vol_pca_matrix(landscape_dir, k_sketch)
    n = int(pc.shape[0])
    sorted_vols = kmeans_sorted_vol_indices(kmeans_dir)
    vol_ids = sorted_vols[:n] if len(sorted_vols) >= n else sorted_vols
    if len(vol_ids) < n:
        raise ValueError(
            f"Found {len(vol_ids)} vol_*.mrc files but vol_pca has {n} rows."
        )
    centers_plot_rows = load_sketch_centroid_plot_df_rows(kmeans_dir, n)
    states = load_sketch_state_labels(landscape_dir)
    if states is not None and len(states) != n:
        states = None
    plotly_cs = normalize_continuous_palette(continuous_palette)
    color_mode = (color_mode or "none").strip().lower()
    df = exp.plot_df

    if color_mode == "state" and states is not None:
        ser = pd.Series(states)
        colors, _ = _labels_colors_and_legend_items(ser)
        return {int(vol_ids[i]): str(colors[i]) for i in range(n)}
    if (
        color_mode not in ("none", "state")
        and color_mode in df.columns
        and color_mode in exp.numeric_columns
    ):
        cov_vals: list[Any] = []
        for i in range(n):
            r = int(centers_plot_rows[i])
            if r < 0 or r >= len(df):
                cov_vals.append(np.nan)
            else:
                cov_vals.append(df.iloc[r][color_mode])
        ser = pd.Series(cov_vals)
        if color_mode == "labels":
            colors, _ = _labels_colors_and_legend_items(ser)
            return {int(vol_ids[i]): str(colors[i]) for i in range(n)}
        color_num = pd.to_numeric(ser, errors="coerce")
        hexes = _numeric_array_to_plotly_hex(
            color_num.to_numpy(dtype=float),
            plotly_cs,
        )
        return {int(vol_ids[i]): hexes[i] for i in range(n)}
    if color_mode == "none":
        return {}
    return {int(vol_ids[i]): "#4a5568" for i in range(n)}


def landscape_volpca_scatter_figure(
    landscape_dir: str,
    exp: DashboardExperiment,
    *,
    pc_x: int,
    pc_y: int,
    color_mode: str,
    continuous_palette: str | None = None,
) -> go.Figure:
    k_sketch, kmeans_dir = resolve_kmeans_sketch_bundle(landscape_dir)
    pc = load_vol_pca_matrix(landscape_dir, k_sketch)
    n, dim = pc.shape
    if pc_x < 0 or pc_x >= dim or pc_y < 0 or pc_y >= dim:
        raise ValueError(f"PCA axis out of range (0..{dim - 1}).")
    if pc_x == pc_y:
        raise ValueError("Choose two distinct PCA axes.")

    sorted_vols = kmeans_sorted_vol_indices(kmeans_dir)
    if len(sorted_vols) < n:
        # Still plot first n rows aligned with analyze_landscape ordering
        pass
    vol_ids = sorted_vols[:n] if len(sorted_vols) >= n else sorted_vols
    if len(vol_ids) < n:
        raise ValueError(
            f"Found {len(vol_ids)} vol_*.mrc files but vol_pca has {n} rows."
        )

    centers_plot_rows = load_sketch_centroid_plot_df_rows(kmeans_dir, n)
    _ai = exp.all_indices
    train_image_idx_str: list[str] = []
    for i in range(n):
        r = int(centers_plot_rows[i])
        if 0 <= r < len(_ai):
            train_image_idx_str.append(str(int(_ai[r])))
        else:
            train_image_idx_str.append("—")
    train_idx_col = np.asarray(train_image_idx_str, dtype=object)

    states = load_sketch_state_labels(landscape_dir)
    if states is not None and len(states) != n:
        states = None

    evr = load_pca_explained_variance(landscape_dir)
    xtitle = f"Volume PC{pc_x + 1}"
    ytitle = f"Volume PC{pc_y + 1}"
    if evr is not None and pc_x < len(evr):
        xtitle += f" (EV: {float(evr[pc_x]):.4f})"
    if evr is not None and pc_y < len(evr):
        ytitle += f" (EV: {float(evr[pc_y]):.4f})"

    plotly_cs = normalize_continuous_palette(continuous_palette)
    color_mode = (color_mode or "none").strip().lower()
    df = exp.plot_df

    if color_mode == "state" and states is not None:
        ser = pd.Series(states)
        colors, _ = _labels_colors_and_legend_items(ser)
        marker = dict(size=9, opacity=0.75, color=colors)
        customdata = np.column_stack([vol_ids, states, train_idx_col])
        hover = (
            "vol %{customdata[0]} · training image %{customdata[2]} "
            "· state %{customdata[1]}<extra></extra>"
        )
    elif (
        color_mode not in ("none", "state")
        and color_mode in df.columns
        and color_mode in exp.numeric_columns
    ):
        cov_vals: list[Any] = []
        for i in range(n):
            r = int(centers_plot_rows[i])
            if r < 0 or r >= len(df):
                cov_vals.append(np.nan)
            else:
                cov_vals.append(df.iloc[r][color_mode])
        ser = pd.Series(cov_vals)
        if color_mode == "labels":
            colors, _ = _labels_colors_and_legend_items(ser)
            marker = dict(size=9, opacity=0.75, color=colors)
        else:
            marker = dict(
                size=9,
                opacity=0.75,
                color=ser,
                colorscale=plotly_cs,
            )
        color_num = pd.to_numeric(ser, errors="coerce")
        disp: list[Any] = []
        for i in range(n):
            if color_mode == "labels":
                disp.append(_lower_legend_entry_label("labels", ser.iloc[i]))
            else:
                v = color_num.iloc[i]
                if pd.isna(v):
                    disp.append(None)
                else:
                    fv = float(v)
                    disp.append(fv if np.isfinite(fv) else None)
        customdata = np.column_stack(
            [
                np.asarray(vol_ids, dtype=np.int64),
                np.asarray(disp, dtype=object),
                train_idx_col,
            ],
        )
        hover = (
            "vol %{customdata[0]} · training image %{customdata[2]} "
            "· %{customdata[1]}<extra></extra>"
        )
    else:
        marker = dict(size=9, opacity=0.75, color="#4a5568")
        customdata = np.column_stack(
            [np.asarray(vol_ids, dtype=np.int64), train_idx_col],
        )
        hover = "vol %{customdata[0]} · training image %{customdata[1]}<extra></extra>"

    sketch_ids = [str(int(v)) for v in vol_ids]

    sc = go.Scattergl(
        x=pc[:, pc_x],
        y=pc[:, pc_y],
        mode="markers",
        ids=sketch_ids,
        customdata=customdata,
        hovertemplate=hover,
        marker=marker,
    )
    fig = go.Figure(sc)
    fig.update_layout(
        template="plotly_white",
        paper_bgcolor=_DASHBOARD_CREAM,
        margin=dict(l=50, r=20, t=16, b=50),
        title=None,
        xaxis=dict(title=xtitle),
        yaxis=dict(title=ytitle),
        dragmode="lasso",
        uirevision="vol_sketch_landscape",
        font=_PLOTLY_FONT,
        showlegend=False,
    )
    return fig


def landscape_volpca_scatter_json(
    landscape_dir: str,
    exp: DashboardExperiment,
    *,
    pc_x: int,
    pc_y: int,
    color_mode: str,
    continuous_palette: str | None = None,
) -> str:
    fig = landscape_volpca_scatter_figure(
        landscape_dir,
        exp,
        pc_x=pc_x,
        pc_y=pc_y,
        color_mode=color_mode,
        continuous_palette=continuous_palette,
    )
    return _plotly_to_json(fig)


def _landscape_anim_root_dir() -> str:
    global _LANDSCAPE_ANIM_ROOT
    if _LANDSCAPE_ANIM_ROOT is None:
        _LANDSCAPE_ANIM_ROOT = tempfile.mkdtemp(
            prefix="cryodrgn_dash_landscape_anim_",
        )
    return _LANDSCAPE_ANIM_ROOT


def clear_landscape_animation_cache() -> None:
    """Remove all generated landscape animation files and the temp root directory."""
    global _LANDSCAPE_ANIM_ROOT
    with _LANDSCAPE_ANIM_LOCK:
        _LANDSCAPE_ANIM_ENTRIES.clear()
        if _LANDSCAPE_ANIM_ROOT and os.path.isdir(_LANDSCAPE_ANIM_ROOT):
            shutil.rmtree(_LANDSCAPE_ANIM_ROOT, ignore_errors=True)
        _LANDSCAPE_ANIM_ROOT = None


atexit.register(clear_landscape_animation_cache)


def _landscape_cycle_gif_timing(frames_per_vol: int) -> tuple[int, int]:
    """Clamp ``frames_per_vol`` and return ``(frames_per_vol, duration_ms_per_frame)`` for cycle GIFs."""
    fpv = max(2, min(int(frames_per_vol), 30))
    duration_ms = max(30, int(800 / fpv))
    return fpv, duration_ms


def _gif_first_frame_to_png(gif_path: str, out_png: str) -> None:
    """Save frame 0 of an animated GIF as a PNG (e.g. rotate preview → cycle keyframe)."""
    from PIL import Image

    with Image.open(gif_path) as im:
        im.seek(0)
        im.convert("RGBA").save(out_png)


def _cycle_gif_from_png_paths(
    png_paths: list[str],
    out_gif: str,
    *,
    cycle_frames_per_vol: int,
) -> None:
    """Hold each static PNG for ``frames_per_vol`` GIF frames (one ChimeraX view per volume)."""
    from PIL import Image

    fpv, duration_ms = _landscape_cycle_gif_timing(cycle_frames_per_vol)
    if not png_paths:
        raise ValueError("No PNG paths for cycle GIF.")
    frames: list[Any] = []
    try:
        for pth in png_paths:
            with Image.open(pth) as im:
                rgb = im.convert("RGB")
                for _ in range(fpv):
                    frames.append(rgb.copy())
        frames[0].save(
            out_gif,
            save_all=True,
            append_images=frames[1:],
            duration=duration_ms,
            loop=0,
        )
    finally:
        for f in frames:
            try:
                f.close()
            except Exception:
                pass


def _lookup_rotate_keyframe_pngs(
    reuse_rotate_keyframes_token: str | None,
    landscape_dir: str,
    color_mode: str,
    plot_color_mode: str,
    continuous_palette: str | None,
) -> dict[int, str] | None:
    """Return ``vol → png`` from a prior rotate batch if valid for reuse, else ``None``."""
    if (
        not reuse_rotate_keyframes_token
        or not str(reuse_rotate_keyframes_token).strip()
    ):
        return None
    tok = str(reuse_rotate_keyframes_token).strip()
    with _LANDSCAPE_ANIM_LOCK:
        prev = _LANDSCAPE_ANIM_ENTRIES.get(tok)
    if not prev:
        return None
    if os.path.abspath(prev["landscape_dir"]) != os.path.abspath(landscape_dir):
        return None
    if prev.get("mode") != "rotate_each":
        return None
    if prev.get("color_mode") != color_mode:
        return None
    if "plot_color_mode" in prev and prev.get("plot_color_mode") != plot_color_mode:
        return None
    if "continuous_palette" in prev:
        if (prev.get("continuous_palette") or None) != (continuous_palette or None):
            return None
    rk = prev.get("rotate_keyframes")
    if not isinstance(rk, dict) or not rk:
        return None
    out: dict[int, str] = {}
    for k, v in rk.items():
        p = str(v)
        if os.path.isfile(p):
            out[int(k)] = p
    return out or None


def generate_landscape_volume_animations(
    landscape_dir: str,
    vol_indices: list[int],
    *,
    exp: DashboardExperiment,
    mode: str,
    gif_frames: int = DEFAULT_GIF_FRAMES,
    chimerax_cpus: int = DEFAULT_CHIMERAX_PARALLEL,
    cycle_frames_per_vol: int = 8,
    color_mode: str = "none",
    plot_color_mode: str = "none",
    continuous_palette: str | None = None,
    reuse_rotate_keyframes_token: str | None = None,
) -> tuple[str, list[dict[str, Any]], list[int]]:
    """Write GIF(s) under a new token directory.

    Returns ``(token, file_metadata, rendered_vol_indices)`` — the last list is the
    sketch volume indices actually rendered (after subsampling caps), sorted unique.
    """
    vol_indices = sorted({int(v) for v in vol_indices})
    if not vol_indices:
        raise ValueError("Select at least one k-means volume.")

    mode = (mode or "cycle").strip().lower()
    if mode == "both":
        raise ValueError(
            'mode "both" is no longer supported; use "cycle" or "rotate_each".'
        )
    if mode not in ("rotate_each", "cycle"):
        raise ValueError('mode must be "cycle" or "rotate_each".')

    k_sketch, kmeans_dir = resolve_kmeans_sketch_bundle(landscape_dir)
    valid_sketch = set(kmeans_sorted_vol_indices(kmeans_dir))
    vol_indices = sorted(set(vol_indices) & valid_sketch)
    if not vol_indices:
        raise ValueError(
            "No valid k-means sketch volumes in the selection — indices must match "
            "vol_NNN.mrc files exported under the landscape k-means folder."
        )

    label_for_vol: dict[int, str] = {
        int(v): montage_cell_label(i) for i, v in enumerate(vol_indices)
    }

    pcm = (plot_color_mode or "none").strip().lower()
    if pcm == "state" and load_sketch_state_labels(landscape_dir) is None:
        pcm = "none"
    hex_by_vol_anim = sketch_vol_marker_hex_by_vol_index(
        landscape_dir,
        exp,
        color_mode=pcm,
        continuous_palette=continuous_palette,
    )

    def _label_hex(vol: int) -> str:
        if pcm == "none":
            return ""
        return hex_by_vol_anim.get(int(vol), "#4a5568")

    rng = random.Random()
    rot_vols: list[int] = []
    cyc_vols: list[int] = []

    if mode == "rotate_each":
        rot_vols = _sample_landscape_vols(vol_indices, LANDSCAPE_ANIM_MAX_ROTATE, rng)
    else:
        cyc_vols = _sample_landscape_vols(vol_indices, LANDSCAPE_ANIM_MAX_CYCLE, rng)

    color_mode = (color_mode or "none").strip().lower()
    vol_state_hex: dict[int, str] | None = None
    if color_mode == "state":
        vol_state_hex = landscape_vol_state_hex_by_vol_index(
            landscape_dir, kmeans_dir, k_sketch
        )

    def _vol_fill(vol: int) -> str | None:
        if vol_state_hex is None:
            return None
        return vol_state_hex.get(int(vol))

    gif_frames = max(4, min(int(gif_frames), 120))
    chimerax_cpus = max(1, min(int(chimerax_cpus), 32))

    token = secrets.token_urlsafe(20)
    job_dir = os.path.join(_landscape_anim_root_dir(), token)
    os.makedirs(job_dir, exist_ok=True)

    out_files: list[dict[str, Any]] = []
    rotate_keyframes: dict[int, str] = {}
    keyframe_dir = os.path.join(job_dir, "keyframes")

    if mode == "rotate_each" and rot_vols:
        os.makedirs(keyframe_dir, exist_ok=True)
        for v in rot_vols:
            mp = vol_mrc_path(kmeans_dir, v)
            out_gif = os.path.join(job_dir, f"vol_{v:03d}_rotate.gif")
            mrc_to_rotating_gif(
                mp,
                out_gif,
                gif_frames=gif_frames,
                ncpus=chimerax_cpus,
                volume_color=_vol_fill(v),
            )
            kpng = os.path.join(keyframe_dir, f"vol_{int(v):03d}.png")
            _gif_first_frame_to_png(out_gif, kpng)
            rotate_keyframes[int(v)] = kpng
            out_files.append(
                {
                    "vol": v,
                    "kind": "rotate",
                    "filename": os.path.basename(out_gif),
                    "path": out_gif,
                    "label": label_for_vol[int(v)],
                    "preview_overlay": {
                        "style": "static",
                        "text": label_for_vol[int(v)],
                        "badge_background": _label_hex(int(v)),
                    },
                }
            )

    reuse_pngs = _lookup_rotate_keyframe_pngs(
        reuse_rotate_keyframes_token,
        landscape_dir,
        color_mode,
        pcm,
        continuous_palette,
    )

    if mode == "cycle" and cyc_vols:
        fpv_c, frame_ms = _landscape_cycle_gif_timing(cycle_frames_per_vol)
        png_sequence: list[str] = []
        tmp_cycle_pngs: list[str] = []
        try:
            for v in cyc_vols:
                src = reuse_pngs.get(int(v)) if reuse_pngs else None
                if src:
                    png_sequence.append(src)
                else:
                    mp = vol_mrc_path(kmeans_dir, v)
                    tmp = os.path.join(job_dir, f"_cyc_vol_{int(v):03d}.png")
                    mrc_to_static_png(
                        mp,
                        tmp,
                        dpi=100,
                        volume_color=_vol_fill(v),
                    )
                    png_sequence.append(tmp)
                    tmp_cycle_pngs.append(tmp)

            if len(cyc_vols) == 1:
                v0 = cyc_vols[0]
                out_gif = os.path.join(job_dir, f"vol_{v0:03d}_cycle.gif")
                _cycle_gif_from_png_paths(
                    png_sequence,
                    out_gif,
                    cycle_frames_per_vol=cycle_frames_per_vol,
                )
                out_files.append(
                    {
                        "vol": v0,
                        "kind": "cycle",
                        "filename": os.path.basename(out_gif),
                        "path": out_gif,
                        "label": label_for_vol[int(v0)],
                        "preview_overlay": {
                            "style": "static",
                            "text": label_for_vol[int(v0)],
                            "badge_background": _label_hex(int(v0)),
                        },
                    }
                )
            else:
                out_gif = os.path.join(
                    job_dir,
                    "cycle_" + "_".join(f"{v:03d}" for v in cyc_vols) + ".gif",
                )
                _cycle_gif_from_png_paths(
                    png_sequence,
                    out_gif,
                    cycle_frames_per_vol=cycle_frames_per_vol,
                )
                out_files.append(
                    {
                        "vol": None,
                        "kind": "cycle",
                        "filename": os.path.basename(out_gif),
                        "path": out_gif,
                        "label": ", ".join(label_for_vol[int(v)] for v in cyc_vols),
                        "preview_overlay": {
                            "style": "cycle_segments",
                            "segment_labels": [label_for_vol[int(v)] for v in cyc_vols],
                            "segment_backgrounds": [
                                _label_hex(int(v)) for v in cyc_vols
                            ],
                            "frames_per_segment": fpv_c,
                            "frame_duration_ms": frame_ms,
                            "total_frames": int(len(cyc_vols) * fpv_c),
                        },
                    }
                )
        finally:
            for tpath in tmp_cycle_pngs:
                try:
                    os.remove(tpath)
                except OSError:
                    pass

    with _LANDSCAPE_ANIM_LOCK:
        entry: dict[str, Any] = {
            "job_dir": job_dir,
            "files": out_files,
            "landscape_dir": os.path.abspath(landscape_dir),
            "kmeans_k": k_sketch,
            "t0": time.monotonic(),
            "mode": mode,
            "color_mode": color_mode,
        }
        if rotate_keyframes:
            entry["rotate_keyframes"] = rotate_keyframes
        entry["plot_color_mode"] = pcm
        entry["continuous_palette"] = continuous_palette
        _LANDSCAPE_ANIM_ENTRIES[token] = entry

    rendered_vol_indices = sorted(
        {int(v) for v in rot_vols} | {int(v) for v in cyc_vols}
    )
    return token, out_files, rendered_vol_indices


def animation_payload_b64(token: str) -> list[dict[str, Any]]:
    with _LANDSCAPE_ANIM_LOCK:
        meta = _LANDSCAPE_ANIM_ENTRIES.get(token)
        if not meta:
            raise ValueError("Unknown or expired animation batch id.")
    items: list[dict[str, Any]] = []
    for ent in meta["files"]:
        p = ent["path"]
        with open(p, "rb") as fh:
            b64 = base64.standard_b64encode(fh.read()).decode("ascii")
        item: dict[str, Any] = {
            "vol": ent["vol"],
            "kind": ent["kind"],
            "filename": ent["filename"],
            "gif_b64": b64,
        }
        if "label" in ent:
            item["label"] = ent["label"]
        if "preview_overlay" in ent:
            item["preview_overlay"] = ent["preview_overlay"]
        items.append(item)
    return items


def save_landscape_animations(
    token: str,
    out_dir: str | None,
    *,
    exp: DashboardExperiment,
    landscape_epoch: int,
) -> list[str]:
    """Copy GIFs from the temp batch to ``out_dir`` (default: landscape kmeans folder)."""
    landscape_dir = landscape_dir_for_epoch(exp.workdir, landscape_epoch)
    k_sketch, kmeans_dir = resolve_kmeans_sketch_bundle(landscape_dir)

    if out_dir is None or not str(out_dir).strip():
        target = kmeans_dir
    else:
        target = os.path.abspath(str(out_dir).strip())
        if not os.path.isdir(target):
            raise ValueError(f"Output directory does not exist: {target}")

    with _LANDSCAPE_ANIM_LOCK:
        meta = _LANDSCAPE_ANIM_ENTRIES.get(token)
        if not meta:
            raise ValueError("Unknown or expired animation batch id.")
        if os.path.abspath(meta["landscape_dir"]) != os.path.abspath(landscape_dir):
            raise ValueError("Animation token does not match this landscape directory.")
        if int(meta["kmeans_k"]) != int(k_sketch):
            raise ValueError("Animation token does not match k-means sketch size.")
        files = list(meta["files"])

    os.makedirs(target, exist_ok=True)
    saved: list[str] = []
    for ent in files:
        src = ent["path"]
        dst_name = ent["filename"]
        dst = os.path.join(target, dst_name)
        if os.path.isfile(dst):
            root, ext = os.path.splitext(dst_name)
            dst = os.path.join(target, f"{root}_saved{ext}")
        shutil.copy2(src, dst)
        saved.append(dst)
    return saved


def meta_for_api(exp: DashboardExperiment) -> dict[str, Any]:
    """JSON-serializable summary for the volume sketched landscape explorer (current epoch)."""
    epochs = list_landscape_epochs(exp.workdir)
    if not epochs:
        return {"ok": False, "error": "No landscape.N directories found.", "epochs": []}

    le = int(exp.epoch)
    if le not in epochs:
        return {
            "ok": False,
            "error": (
                f"No landscape.{le} folder for the current dashboard epoch. "
                "Switch epoch in the nav or run analyze_landscape for this epoch."
            ),
            "epochs": epochs,
        }

    landscape_dir = landscape_dir_for_epoch(exp.workdir, le)
    try:
        k_sketch, kmeans_dir = resolve_kmeans_sketch_bundle(landscape_dir)
    except FileNotFoundError as err:
        return {"ok": False, "error": str(err), "epochs": epochs}

    pc = load_vol_pca_matrix(landscape_dir, k_sketch)
    _, dim = pc.shape
    evr = load_pca_explained_variance(landscape_dir)
    evr_list = [float(x) for x in evr] if evr is not None else None
    states = load_sketch_state_labels(landscape_dir)
    chimerax_cpus = max(1, min(int(DEFAULT_CHIMERAX_PARALLEL), 32))

    return {
        "ok": True,
        "landscape_epoch": le,
        "landscape_dir": landscape_dir,
        "kmeans_k": k_sketch,
        "kmeans_dir": kmeans_dir,
        "n_volumes": int(pc.shape[0]),
        "n_pc": int(dim),
        "explained_variance_ratio": evr_list,
        "has_state_color": states is not None,
        "color_options": landscape_color_options(exp),
        "default_save_dir": kmeans_dir,
        "chimerax_cpus": int(chimerax_cpus),
    }
