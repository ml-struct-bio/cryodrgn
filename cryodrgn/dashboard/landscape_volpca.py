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
from cryodrgn import utils
from cryodrgn.dashboard.covariate_labels import covariate_display_name
from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.plots_color_covariate import numeric_array_to_plotly_hex
from cryodrgn.dashboard.particle_explorer import (
    ChimeraxViewTurn,
    DEFAULT_CHIMERAX_PARALLEL,
    DEFAULT_GIF_FRAMES,
    chimerax_view_matrix_camera_arg,
    format_chimerax_view_matrix_display,
    montage_cell_label,
    mrc_to_rotating_gif,
    mrc_to_static_png,
)
from cryodrgn.dashboard.palette_config import normalize_continuous_palette
from cryodrgn.dashboard.plots_color_covariate import (
    _labels_colors_and_legend_items,
    _lower_legend_entry_label,
)
from cryodrgn.dashboard.plots_color_covariate import _continuous_series_stats
from cryodrgn.dashboard.plots_figure_utils import (
    _DASHBOARD_CREAM,
    _PLOTLY_FONT,
    DASHBOARD_SCATTER_HOVERLABEL_FONT_SIZE,
    _plotly_to_json,
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
# Particle PCA columns from ``analysis.load_dataframe`` (``PC1``, ``PC2``, …).
_PARTICLE_PC_COV_COL_RE = re.compile(r"^PC(\d+)$", re.IGNORECASE)

# ChimeraX cost: cap how many volumes go into each animation style (extras subsampled).
LANDSCAPE_ANIM_MAX_ROTATE = 5
LANDSCAPE_ANIM_MAX_CYCLE = 50
# Default ChimeraX rotation-frame count (sketch GIFs + vol-PCA ``generate_animations``).
LANDSCAPE_SKETCH_GIF_FRAMES_DEFAULT = 20
_LANDSCAPE_VIEW_ROTATION_AXES: tuple[str, ...] = ("x", "y", "z")

# Vol-PCA scatter: legacy 9 → −13%; match dashboard 3D (+30% size, −20% opacity vs that base).
_VOLSKETCH_SCATTER_MARKER = float(9 * (1.0 - 0.13) * 1.3)
_VOLSKETCH_SCATTER_OPACITY = float(min(0.98, 0.75 * 0.8))


def _sample_landscape_vols(vols: list[int], k: int, rng: random.Random) -> list[int]:
    """Up to ``k`` volumes without replacement; sorted ascending."""
    if len(vols) <= k:
        return list(vols)
    return sorted(rng.sample(vols, k))


def normalize_landscape_view_rotations(
    view_rotations: dict[str, Any] | None,
) -> list[ChimeraxViewTurn]:
    """Return validated base-view rotations for ChimeraX animation previews."""
    if not view_rotations:
        return []
    if not isinstance(view_rotations, dict):
        raise ValueError("view_rotations must be an object with x, y, and z degrees.")
    out: list[ChimeraxViewTurn] = []
    for axis in _LANDSCAPE_VIEW_ROTATION_AXES:
        raw = view_rotations.get(axis, 0.0)
        try:
            degrees = float(raw)
        except (TypeError, ValueError) as err:
            raise ValueError(f"view_rotations.{axis} must be a number.") from err
        if not np.isfinite(degrees):
            raise ValueError(f"view_rotations.{axis} must be finite.")
        if abs(degrees) > 1e-9:
            out.append((axis, degrees))
    return out


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
    """Vol index (``11`` for ``vol_011.mrc``) → fill hex, matching state scatter colors."""
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


def sketch_plot_color_covariate_variable_label(plot_color_mode: str) -> str | None:
    """Human-readable color dimension name for GIF preview (matches Color by column)."""
    raw = (plot_color_mode or "").strip()
    pcm = raw.lower()
    if pcm in ("none", ""):
        return None
    if pcm == "state":
        return "Agglomerative state"
    m = _PARTICLE_PC_COV_COL_RE.match(raw)
    if m:
        return f"latent-space PC{m.group(1)}"
    return _volsketch_covariate_display(pcm)


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


def load_landscape_sketch_umap_coords(
    landscape_dir: str,
    kmeans_dir: str,
    n: int,
) -> np.ndarray | None:
    """Sketch-centroid rows in the particle UMAP from ``landscape.*/umap.pkl``.

    Row ``i`` aligns with ``vol_pca`` row ``i`` via ``centers_ind.txt`` (same as
    ``analyze_landscape`` plots). Returns ``(n, d)`` float64 or ``None`` if UMAP
    is missing, mis-shaped, or indices are out of range.
    """
    umap_path = os.path.join(landscape_dir, "umap.pkl")
    if not os.path.isfile(umap_path):
        return None
    try:
        umap_full = np.asarray(utils.load_pkl(umap_path), dtype=np.float64)
    except Exception:
        return None
    if umap_full.ndim != 2 or umap_full.shape[0] < 1 or umap_full.shape[1] < 1:
        return None
    try:
        centers = load_sketch_centroid_plot_df_rows(kmeans_dir, n)
    except (FileNotFoundError, ValueError):
        return None
    centers_i = centers.astype(np.intp, copy=False)
    if np.any(centers_i < 0) or np.any(centers_i >= umap_full.shape[0]):
        return None
    return umap_full[centers_i].astype(np.float64, copy=False)


def parse_landscape_axis_spec(spec: str) -> tuple[str, int]:
    """Parse ``pc:0`` / ``umap:1`` style axis strings from the scatter API."""
    s = (spec or "").strip().lower().replace(" ", "")
    if ":" not in s:
        raise ValueError(
            f"Invalid axis {spec!r}; expected a form like pc:0 or umap:1.",
        )
    kind, idx_s = s.split(":", 1)
    if kind not in ("pc", "umap"):
        raise ValueError(
            f"Unknown axis kind in {spec!r} (expected pc or umap).",
        )
    try:
        idx = int(idx_s)
    except ValueError as err:
        raise ValueError(f"Invalid axis index in {spec!r}.") from err
    if idx < 0:
        raise ValueError("Axis index must be non-negative.")
    return (kind, idx)


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
        hexes = numeric_array_to_plotly_hex(
            color_num.to_numpy(dtype=float),
            plotly_cs,
        )
        return {int(vol_ids[i]): hexes[i] for i in range(n)}
    if color_mode == "none":
        return {}
    return {int(vol_ids[i]): "#4a5568" for i in range(n)}


def sketch_vol_color_covariate_overlay_text(
    landscape_dir: str,
    kmeans_dir: str,
    k_sketch: int,
    exp: DashboardExperiment,
    plot_color_mode: str,
    vol_index: int,
) -> str | None:
    """Display string for the plot color column at this sketch volume (preview badge)."""
    pcm = (plot_color_mode or "none").strip().lower()
    if pcm == "none":
        return None

    n = int(load_vol_pca_matrix(landscape_dir, k_sketch).shape[0])
    sorted_vols = kmeans_sorted_vol_indices(kmeans_dir)
    vol_ids = sorted_vols[:n] if len(sorted_vols) >= n else sorted_vols
    if len(vol_ids) < n:
        return None
    try:
        row_idx = vol_ids.index(int(vol_index))
    except ValueError:
        return None

    centers_plot_rows = load_sketch_centroid_plot_df_rows(kmeans_dir, n)
    states = load_sketch_state_labels(landscape_dir)
    if states is not None and len(states) != n:
        states = None
    df = exp.plot_df
    r_plot = int(centers_plot_rows[row_idx])

    if pcm == "state" and states is not None:
        return str(int(states[row_idx]))

    if (
        pcm not in ("none", "state")
        and pcm in df.columns
        and pcm in exp.numeric_columns
    ):
        if r_plot < 0 or r_plot >= len(df):
            return "—"
        raw = df.iloc[r_plot][pcm]
        if pcm == "labels":
            return _lower_legend_entry_label("labels", raw)
        color_num = pd.to_numeric(pd.Series([raw]), errors="coerce").iloc[0]
        if pd.isna(color_num):
            return "—"
        fv = float(color_num)
        if not np.isfinite(fv):
            return "—"
        return f"{fv:.6g}"

    return None


def _plot_color_mode_is_continuous_numeric(exp: DashboardExperiment, pcm: str) -> bool:
    pcm = (pcm or "none").strip().lower()
    return bool(
        pcm not in ("none", "state", "labels")
        and pcm in exp.plot_df.columns
        and pcm in exp.numeric_columns
    )


def sketch_vol_rotate_frames_covariate_overlay(
    exp: DashboardExperiment,
    landscape_dir: str,
    kmeans_dir: str,
    k_sketch: int,
    plot_color_mode: str,
    vol_index: int,
    gif_frames: int,
    continuous_palette: str | None,
) -> dict[str, Any] | None:
    """Per-frame covariate text + palette hex for rotate-each GIF previews (continuous colour only)."""
    pcm = (plot_color_mode or "none").strip().lower()
    if not _plot_color_mode_is_continuous_numeric(exp, pcm):
        return None

    from cryodrgn.dashboard.landscape_full_3d import (
        VOL_LANDSCAPE_NEAREST_SKETCH_VOL,
        attach_landscape_nearest_sketch_vol_column,
    )

    df = exp.plot_df
    if VOL_LANDSCAPE_NEAREST_SKETCH_VOL not in df.columns:
        df = attach_landscape_nearest_sketch_vol_column(exp, df)
    if VOL_LANDSCAPE_NEAREST_SKETCH_VOL not in df.columns:
        return None

    mask = df[VOL_LANDSCAPE_NEAREST_SKETCH_VOL].to_numpy(dtype=np.int64) == int(
        vol_index
    )
    if not mask.any():
        return None
    ser = pd.to_numeric(df.loc[mask, pcm], errors="coerce").dropna()
    if ser.empty:
        return None

    plotly_cs = normalize_continuous_palette(continuous_palette)
    _, cmin, cmax = _continuous_series_stats(exp.plot_df[pcm])
    sorted_vals = np.sort(ser.to_numpy(dtype=float))
    n = int(sorted_vals.shape[0])
    gif_frames = max(4, int(gif_frames))
    texts: list[str] = []
    bgs: list[str] = []
    for i in range(gif_frames):
        idx = min(n - 1, int(i * n / gif_frames))
        fv = float(sorted_vals[idx])
        texts.append(f"{fv:.6g}")
        bgs.append(
            numeric_array_to_plotly_hex(
                np.array([fv], dtype=float),
                plotly_cs,
                vmin=cmin,
                vmax=cmax,
            )[0]
        )
    from cryodrgn.dashboard.particle_explorer import GIF_DURATION_S

    frame_ms = max(20, int(GIF_DURATION_S * 1000 / gif_frames))
    out: dict[str, Any] = {
        "style": "rotate_frames",
        "frame_covariate_texts": texts,
        "frame_badge_backgrounds": bgs,
        "frame_duration_ms": frame_ms,
        "total_frames": gif_frames,
    }
    vlab = sketch_plot_color_covariate_variable_label(pcm)
    if vlab:
        out["covariate_variable_label"] = vlab
    return out


def landscape_volpca_scatter_figure(
    landscape_dir: str,
    exp: DashboardExperiment,
    *,
    axis_x: tuple[str, int],
    axis_y: tuple[str, int],
    color_mode: str,
    continuous_palette: str | None = None,
) -> go.Figure:
    k_sketch, kmeans_dir = resolve_kmeans_sketch_bundle(landscape_dir)
    pc = load_vol_pca_matrix(landscape_dir, k_sketch)
    n, dim = pc.shape
    kx, ix = axis_x
    ky, iy = axis_y
    if (kx, ix) == (ky, iy):
        raise ValueError("Choose two distinct axes.")

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
    umap_sk = load_landscape_sketch_umap_coords(landscape_dir, kmeans_dir, n)

    def _axis_coord_and_title(
        kind: str,
        index: int,
        *,
        evr_local: np.ndarray | None,
        umap_local: np.ndarray | None,
        pc_local: np.ndarray,
        dim_local: int,
    ) -> tuple[np.ndarray, str]:
        if kind == "pc":
            if index < 0 or index >= dim_local:
                raise ValueError(f"PCA axis out of range (0..{dim_local - 1}).")
            title = f"Vol PC{index + 1}"
            if evr_local is not None and index < len(evr_local):
                title += f" ({100.0 * float(evr_local[index]):.1f}%)"
            return pc_local[:, index], title
        if kind == "umap":
            if umap_local is None:
                raise ValueError(
                    "UMAP axes need landscape umap.pkl (particle embedding) "
                    "and valid centroid indices in centers_ind.txt.",
                )
            if index < 0 or index >= umap_local.shape[1]:
                raise ValueError(
                    f"UMAP axis out of range (0..{umap_local.shape[1] - 1}).",
                )
            return umap_local[:, index], f"Volume UMAP{index + 1}"
        raise ValueError(f"Unknown axis kind {kind!r}.")

    xv, xtitle = _axis_coord_and_title(
        kx,
        ix,
        evr_local=evr,
        umap_local=umap_sk,
        pc_local=pc,
        dim_local=dim,
    )
    yv, ytitle = _axis_coord_and_title(
        ky,
        iy,
        evr_local=evr,
        umap_local=umap_sk,
        pc_local=pc,
        dim_local=dim,
    )

    plotly_cs = normalize_continuous_palette(continuous_palette)
    color_mode = (color_mode or "none").strip().lower()
    df = exp.plot_df

    if color_mode == "state" and states is not None:
        ser = pd.Series(states)
        colors, _ = _labels_colors_and_legend_items(ser)
        marker = dict(
            size=_VOLSKETCH_SCATTER_MARKER,
            opacity=_VOLSKETCH_SCATTER_OPACITY,
            color=colors,
        )
        customdata = np.column_stack([vol_ids, states, train_idx_col])
        hover = "volume: %{customdata[0]}<br>%{customdata[1]}<extra></extra>"
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
            marker = dict(
                size=_VOLSKETCH_SCATTER_MARKER,
                opacity=_VOLSKETCH_SCATTER_OPACITY,
                color=colors,
            )
        else:
            marker = dict(
                size=_VOLSKETCH_SCATTER_MARKER,
                opacity=_VOLSKETCH_SCATTER_OPACITY,
                color=ser,
                colorscale=plotly_cs,
            )
        color_num = pd.to_numeric(ser, errors="coerce")
        disp: list[Any] = []
        for i in range(n):
            if color_mode == "labels":
                raw = ser.iloc[i]
                if pd.isna(raw):
                    disp.append("(missing)")
                else:
                    disp.append(f"kmeans={_lower_legend_entry_label('labels', raw)}")
            else:
                v = color_num.iloc[i]
                lab_nm = covariate_display_name(color_mode)
                if pd.isna(v):
                    disp.append(f"{lab_nm}: —")
                else:
                    fv = float(v)
                    if np.isfinite(fv):
                        disp.append(f"{lab_nm}: {fv:.5g}")
                    else:
                        disp.append(f"{lab_nm}: —")
        customdata = np.column_stack(
            [
                np.asarray(vol_ids, dtype=np.int64),
                np.asarray(disp, dtype=object),
                train_idx_col,
            ],
        )
        hover = "volume: %{customdata[0]}<br>%{customdata[1]}<extra></extra>"
    else:
        marker = dict(
            size=_VOLSKETCH_SCATTER_MARKER,
            opacity=_VOLSKETCH_SCATTER_OPACITY,
            color="#4a5568",
        )
        customdata = np.column_stack(
            [np.asarray(vol_ids, dtype=np.int64), train_idx_col],
        )
        hover = "volume: %{customdata[0]}<extra></extra>"

    sketch_ids = [str(int(v)) for v in vol_ids]

    sc = go.Scattergl(
        x=xv,
        y=yv,
        mode="markers",
        ids=sketch_ids,
        customdata=customdata,
        hovertemplate=hover,
        marker=marker,
    )
    fig = go.Figure(sc)
    hov_font = dict(_PLOTLY_FONT)
    hov_font["size"] = DASHBOARD_SCATTER_HOVERLABEL_FONT_SIZE
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
        hoverlabel=dict(font=hov_font, align="left"),
    )
    return fig


def landscape_volpca_scatter_json(
    landscape_dir: str,
    exp: DashboardExperiment,
    *,
    axis_x: str | tuple[str, int],
    axis_y: str | tuple[str, int],
    color_mode: str,
    continuous_palette: str | None = None,
) -> str:
    ax_x = parse_landscape_axis_spec(axis_x) if isinstance(axis_x, str) else axis_x
    ax_y = parse_landscape_axis_spec(axis_y) if isinstance(axis_y, str) else axis_y
    fig = landscape_volpca_scatter_figure(
        landscape_dir,
        exp,
        axis_x=ax_x,
        axis_y=ax_y,
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
    view_turns: list[ChimeraxViewTurn],
    view_matrix_camera: str | None,
) -> tuple[dict[int, str], str | None] | None:
    """Return reusable keyframes and their ChimeraX view matrix, if valid."""
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
    prev_camera = prev.get("view_matrix_camera")
    if view_matrix_camera:
        if str(prev_camera or "") != str(view_matrix_camera):
            return None
    elif prev_camera or tuple(prev.get("view_turns") or ()) != tuple(view_turns):
        return None
    rk = prev.get("rotate_keyframes")
    if not isinstance(rk, dict) or not rk:
        return None
    out: dict[int, str] = {}
    for k, v in rk.items():
        p = str(v)
        if os.path.isfile(p):
            out[int(k)] = p
    if not out:
        return None
    vm = prev.get("view_matrix")
    return out, str(vm) if vm else None


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
    view_rotations: dict[str, Any] | None = None,
    view_matrix: str | None = None,
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
    view_matrix_camera = chimerax_view_matrix_camera_arg(view_matrix)
    if view_matrix_camera:
        view_turns: list[ChimeraxViewTurn] = []
    else:
        view_turns = normalize_landscape_view_rotations(view_rotations)

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

    def _cov_preview_fields(vol: int) -> dict[str, Any]:
        s = sketch_vol_color_covariate_overlay_text(
            landscape_dir,
            kmeans_dir,
            k_sketch,
            exp,
            pcm,
            int(vol),
        )
        if not s:
            return {}
        cov_bg = _label_hex(int(vol)) if pcm != "none" else ""
        out: dict[str, Any] = {
            "covariate_text": s,
            "covariate_badge_background": cov_bg,
        }
        vlab = sketch_plot_color_covariate_variable_label(pcm)
        if vlab:
            out["covariate_variable_label"] = vlab
        return out

    def _rotate_each_preview_overlay(vol: int) -> dict[str, Any]:
        rot = sketch_vol_rotate_frames_covariate_overlay(
            exp,
            landscape_dir,
            kmeans_dir,
            k_sketch,
            pcm,
            int(vol),
            gif_frames,
            continuous_palette,
        )
        if rot:
            rot["text"] = label_for_vol[int(vol)]
            rot["badge_background"] = _label_hex(int(vol))
            return rot
        overlay: dict[str, Any] = {
            "style": "static",
            "text": label_for_vol[int(vol)],
            "badge_background": _label_hex(int(vol)),
        }
        overlay.update(_cov_preview_fields(int(vol)))
        return overlay

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

    def _vol_chimerax_color(vol: int) -> str:
        """Solid colour for ChimeraX ``volume color`` — align with scatter / preview badges."""
        from cryodrgn.dashboard.particle_explorer import chimerax_volume_color_spec

        if vol_state_hex is not None:
            h = vol_state_hex.get(int(vol))
            if h:
                return chimerax_volume_color_spec(str(h))
        if pcm != "none":
            return chimerax_volume_color_spec(
                str(hex_by_vol_anim.get(int(vol), "#4a5568"))
            )
        return chimerax_volume_color_spec("#4a5568")

    gif_frames = max(4, min(int(gif_frames), 120))
    chimerax_cpus = max(1, min(int(chimerax_cpus), 32))

    token = secrets.token_urlsafe(20)
    job_dir = os.path.join(_landscape_anim_root_dir(), token)
    os.makedirs(job_dir, exist_ok=True)

    out_files: list[dict[str, Any]] = []
    rotate_keyframes: dict[int, str] = {}
    view_matrix_text: str | None = None
    keyframe_dir = os.path.join(job_dir, "keyframes")

    if mode == "rotate_each" and rot_vols:
        os.makedirs(keyframe_dir, exist_ok=True)
        for v in rot_vols:
            mp = vol_mrc_path(kmeans_dir, v)
            out_gif = os.path.join(job_dir, f"vol_{v:03d}_rotate.gif")
            rot_overlay = sketch_vol_rotate_frames_covariate_overlay(
                exp,
                landscape_dir,
                kmeans_dir,
                k_sketch,
                pcm,
                int(v),
                gif_frames,
                continuous_palette,
            )
            frame_volume_colors = None
            if rot_overlay and rot_overlay.get("frame_badge_backgrounds"):
                from cryodrgn.dashboard.particle_explorer import (
                    chimerax_volume_color_spec,
                )

                frame_volume_colors = [
                    chimerax_volume_color_spec(str(c))
                    for c in rot_overlay["frame_badge_backgrounds"]
                ]
            vm = mrc_to_rotating_gif(
                mp,
                out_gif,
                gif_frames=gif_frames,
                ncpus=chimerax_cpus,
                volume_color=_vol_chimerax_color(v),
                volume_colors=frame_volume_colors,
                view_turns=view_turns,
                view_matrix_camera=view_matrix_camera,
                report_view_matrix=view_matrix_text is None,
            )
            if vm and view_matrix_text is None:
                view_matrix_text = vm
            elif view_matrix_camera and view_matrix_text is None:
                view_matrix_text = format_chimerax_view_matrix_display(
                    view_matrix_camera
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
                    "preview_overlay": _rotate_each_preview_overlay(int(v)),
                }
            )

    reuse = _lookup_rotate_keyframe_pngs(
        reuse_rotate_keyframes_token,
        landscape_dir,
        color_mode,
        pcm,
        continuous_palette,
        view_turns,
        view_matrix_camera,
    )
    reuse_pngs = reuse[0] if reuse else None
    if reuse and view_matrix_text is None:
        view_matrix_text = reuse[1]

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
                    vm = mrc_to_static_png(
                        mp,
                        tmp,
                        dpi=100,
                        volume_color=_vol_chimerax_color(v),
                        view_turns=view_turns,
                        view_matrix_camera=view_matrix_camera,
                        report_view_matrix=view_matrix_text is None,
                    )
                    if vm and view_matrix_text is None:
                        view_matrix_text = vm
                    elif view_matrix_camera and view_matrix_text is None:
                        view_matrix_text = format_chimerax_view_matrix_display(
                            view_matrix_camera
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
                            **_cov_preview_fields(int(v0)),
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
                cyc_preview: dict[str, Any] = {
                    "style": "cycle_segments",
                    "segment_labels": [label_for_vol[int(v)] for v in cyc_vols],
                    "segment_backgrounds": [_label_hex(int(v)) for v in cyc_vols],
                    "frames_per_segment": fpv_c,
                    "frame_duration_ms": frame_ms,
                    "total_frames": int(len(cyc_vols) * fpv_c),
                }
                if pcm != "none":
                    cyc_preview["segment_covariate_texts"] = [
                        sketch_vol_color_covariate_overlay_text(
                            landscape_dir,
                            kmeans_dir,
                            k_sketch,
                            exp,
                            pcm,
                            int(vx),
                        )
                        or ""
                        for vx in cyc_vols
                    ]
                    if _plot_color_mode_is_continuous_numeric(exp, pcm):
                        cyc_preview["segment_covariate_backgrounds"] = [
                            _label_hex(int(vx)) for vx in cyc_vols
                        ]
                    vlab_c = sketch_plot_color_covariate_variable_label(pcm)
                    if vlab_c:
                        cyc_preview["covariate_variable_label"] = vlab_c
                out_files.append(
                    {
                        "vol": None,
                        "kind": "cycle",
                        "filename": os.path.basename(out_gif),
                        "path": out_gif,
                        "label": ", ".join(label_for_vol[int(v)] for v in cyc_vols),
                        "preview_overlay": cyc_preview,
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
        entry["view_turns"] = tuple(view_turns)
        if view_matrix_camera:
            entry["view_matrix_camera"] = view_matrix_camera
        if view_matrix_text:
            entry["view_matrix"] = view_matrix_text
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
        view_matrix = meta.get("view_matrix")
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
        if view_matrix:
            item["view_matrix"] = str(view_matrix)
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
    rotate_seq = 0
    for ent in files:
        src = ent["path"]
        kind = str(ent.get("kind") or "")
        if kind == "cycle":
            dst_name = "cycle.gif"
        elif kind == "rotate":
            rotate_seq += 1
            dst_name = f"rotate_{rotate_seq:02d}.gif"
        else:
            dst_name = str(ent.get("filename") or os.path.basename(src))
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
    n_vol, dim = pc.shape
    umap_sk = load_landscape_sketch_umap_coords(landscape_dir, kmeans_dir, n_vol)
    n_umap = int(umap_sk.shape[1]) if umap_sk is not None else 0
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
        "n_umap": n_umap,
        "explained_variance_ratio": evr_list,
        "has_state_color": states is not None,
        "color_options": landscape_color_options(exp),
        "default_save_dir": kmeans_dir,
        "chimerax_cpus": int(chimerax_cpus),
    }
