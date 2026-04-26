"""Volume-space PCA scatter and ChimeraX animations for ``analyze_landscape`` outputs."""

from __future__ import annotations

import atexit
import base64
import glob
import os
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

from cryodrgn import utils
from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.explorer_volumes import (
    DEFAULT_CHIMERAX_PARALLEL,
    DEFAULT_GIF_FRAMES,
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
# token -> { "job_dir", "files": [{vol, rel_path}], "landscape_dir", "kmeans_k", "t0", "mode" }
_LANDSCAPE_ANIM_ENTRIES: dict[str, dict[str, Any]] = {}

_LANDSCAPE_VOLPCA_KMEANS_RE = re.compile(r"^kmeans(\d+)$")
_LANDSCAPE_VOLPCA_PKL_RE = re.compile(r"^vol_pca_(\d+)\.pkl$")


def list_landscape_epochs(workdir: str) -> list[int]:
    """Epochs ``N`` with a directory ``landscape.N`` under ``workdir``."""
    if not os.path.isdir(workdir):
        return []
    out: list[int] = []
    for name in os.listdir(workdir):
        m = re.fullmatch(r"landscape\.([0-9]+)", name)
        if m and os.path.isdir(os.path.join(workdir, name)):
            out.append(int(m.group(1)))
    return sorted(out)


def landscape_dir_for_epoch(workdir: str, epoch: int) -> str:
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
    paths = glob.glob(os.path.join(kmeans_dir, "vol_*.mrc"))

    def _idx(p: str) -> int:
        mm = re.search(r"vol_(\d+)", os.path.basename(p))
        return int(mm.group(1)) if mm else 0

    return sorted(_idx(p) for p in paths)


def load_vol_pca_matrix(landscape_dir: str, k_sketch: int) -> np.ndarray:
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


def vol_mrc_path(kmeans_dir: str, vol_index: int) -> str:
    p = os.path.join(kmeans_dir, f"vol_{int(vol_index):03d}.mrc")
    if not os.path.isfile(p):
        raise FileNotFoundError(p)
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

    sc = go.Scattergl(
        x=pc[:, pc_x],
        y=pc[:, pc_y],
        mode="markers",
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


def _mrc_cycle_gif(
    mrc_paths: list[str],
    out_gif: str,
    *,
    frames_per_vol: int = 8,
    dpi: int = 100,
) -> None:
    """GIF that holds each volume at a fixed view for ``frames_per_vol`` frames."""
    from PIL import Image

    frames_per_vol = max(2, min(int(frames_per_vol), 30))
    tmpdir = tempfile.mkdtemp(prefix="cryodrgn_volpca_cycle_")
    try:
        pil_frames: list[Any] = []
        for mrc_path in mrc_paths:
            for _ in range(frames_per_vol):
                png = os.path.join(tmpdir, f"f_{len(pil_frames)}.png")
                mrc_to_static_png(mrc_path, png, dpi=dpi)
                pil_frames.append(Image.open(png))
        if not pil_frames:
            raise ValueError("No frames for cycle GIF.")
        duration_ms = max(30, int(800 / frames_per_vol))
        pil_frames[0].save(
            out_gif,
            save_all=True,
            append_images=pil_frames[1:],
            duration=duration_ms,
            loop=0,
        )
        for im in pil_frames:
            im.close()
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _mrc_single_frame_cycle_gif(mrc_path: str, out_gif: str, *, dpi: int = 100) -> None:
    """One fixed ChimeraX view, single GIF frame (no rotation)."""
    from PIL import Image

    tmpdir = tempfile.mkdtemp(prefix="cryodrgn_volpca_cycle1_")
    try:
        png = os.path.join(tmpdir, "f0.png")
        mrc_to_static_png(mrc_path, png, dpi=dpi)
        im = Image.open(png)
        im.save(out_gif, duration=500, loop=0)
        im.close()
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def generate_landscape_volume_animations(
    landscape_dir: str,
    vol_indices: list[int],
    *,
    mode: str,
    gif_frames: int = DEFAULT_GIF_FRAMES,
    chimerax_cpus: int = DEFAULT_CHIMERAX_PARALLEL,
    cycle_frames_per_vol: int = 8,
) -> tuple[str, list[dict[str, Any]]]:
    """Write GIF(s) under a new token directory; return ``(token, file_metadata)``."""
    vol_indices = sorted({int(v) for v in vol_indices})
    if not vol_indices:
        raise ValueError("Select at least one k-means volume.")
    if len(vol_indices) > 10:
        raise ValueError("At most 10 volumes per animation batch.")

    k_sketch, kmeans_dir = resolve_kmeans_sketch_bundle(landscape_dir)
    mrc_paths = [vol_mrc_path(kmeans_dir, v) for v in vol_indices]

    mode = (mode or "cycle").strip().lower()
    if mode not in ("rotate_each", "cycle", "both"):
        raise ValueError('mode must be "rotate_each", "cycle", or "both".')

    gif_frames = max(4, min(int(gif_frames), 120))
    chimerax_cpus = max(1, min(int(chimerax_cpus), 32))

    token = secrets.token_urlsafe(20)
    job_dir = os.path.join(_landscape_anim_root_dir(), token)
    os.makedirs(job_dir, exist_ok=True)

    out_files: list[dict[str, Any]] = []

    # Rotating GIFs only when the user asks for rotation and (if "both") there
    # are at least two volumes — a single selection stays fixed-view until then.
    do_rotate_each = mode == "rotate_each" or (mode == "both" and len(mrc_paths) >= 2)
    if do_rotate_each:
        for v, mp in zip(vol_indices, mrc_paths):
            out_gif = os.path.join(job_dir, f"vol_{v:03d}_rotate.gif")
            mrc_to_rotating_gif(
                mp,
                out_gif,
                gif_frames=gif_frames,
                ncpus=chimerax_cpus,
            )
            out_files.append(
                {
                    "vol": v,
                    "kind": "rotate",
                    "filename": os.path.basename(out_gif),
                    "path": out_gif,
                }
            )

    if mode in ("cycle", "both"):
        if len(mrc_paths) == 1:
            v0, mp0 = vol_indices[0], mrc_paths[0]
            out_gif = os.path.join(job_dir, f"vol_{v0:03d}_cycle.gif")
            _mrc_single_frame_cycle_gif(mp0, out_gif)
            out_files.append(
                {
                    "vol": v0,
                    "kind": "cycle",
                    "filename": os.path.basename(out_gif),
                    "path": out_gif,
                }
            )
        else:
            out_gif = os.path.join(
                job_dir,
                "cycle_" + "_".join(f"{v:03d}" for v in vol_indices) + ".gif",
            )
            _mrc_cycle_gif(
                mrc_paths,
                out_gif,
                frames_per_vol=cycle_frames_per_vol,
            )
            out_files.append(
                {
                    "vol": None,
                    "kind": "cycle",
                    "filename": os.path.basename(out_gif),
                    "path": out_gif,
                }
            )

    with _LANDSCAPE_ANIM_LOCK:
        _LANDSCAPE_ANIM_ENTRIES[token] = {
            "job_dir": job_dir,
            "files": out_files,
            "landscape_dir": os.path.abspath(landscape_dir),
            "kmeans_k": k_sketch,
            "t0": time.monotonic(),
            "mode": mode,
        }

    return token, out_files


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
        items.append(
            {
                "vol": ent["vol"],
                "kind": ent["kind"],
                "filename": ent["filename"],
                "gif_b64": b64,
            }
        )
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
    }
