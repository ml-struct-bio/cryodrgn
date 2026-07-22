"""Volume decode and analyze-volume catalog for the dashboard trajectory creator.

Decodes latent ``z`` vectors to 3-D volumes (reusing the particle explorer decoder)
and returns float32 arrays for client-side slicing in the browser canvas.
"""

from __future__ import annotations

import base64
import os
import re
import shutil
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional

import numpy as np

from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.chimerax_animation import (
    LandscapeStaticView,
    chimerax_iso_response_fields,
    normalize_chimerax_view_turns,
    render_landscape_cycle_static_views,
    resolve_chimerax_volume_level,
)
from cryodrgn.dashboard.particle_explorer import (
    DEFAULT_CHIMERAX_PARALLEL,
    _decode_z_values_to_vol_paths,
)
from cryodrgn.mrcfile import parse_mrc
from cryodrgn import analysis as cryo_analysis

_VOL_MRC_RE = re.compile(r"^vol_(\d+)\.mrc$", re.IGNORECASE)
_PC_DIR_RE = re.compile(r"^pc(\d+)$", re.IGNORECASE)
_RECON_WINDOW_OUT_RAD = 0.99
# Matches ``CryoVolume3dUtils.PLOT3D_TARGET_D`` / client box-average downsample.
PLOT3D_TARGET_D = 128


def spherical_window_mask_3d(
    vol: Optional[np.ndarray] = None,
    *,
    D: Optional[int] = None,
    in_rad: float = 1.0,
    out_rad: float = 1.0,
) -> np.ndarray:
    """Create a 3-D radial mask for cubic volumes (soft or hard edge).

    Uses normalized radial coordinates with ``in_rad`` / ``out_rad`` semantics
    matching training-time ``spherical_window_mask``, extended to
    ``sqrt(x^2 + y^2 + z^2)``.
    """
    if (vol is None) == (D is None):
        raise ValueError("Either `vol` or `D` must be specified!")
    if vol is not None:
        D = int(vol.shape[0])

    assert in_rad <= out_rad
    ax = np.linspace(-1.0, 1.0, int(D) + 1, dtype=np.float32)[:-1]
    x0, x1, x2 = np.meshgrid(ax, ax, ax, indexing="ij")
    dists = np.sqrt(x0**2 + x1**2 + x2**2)

    if in_rad == out_rad:
        return (dists <= out_rad).astype(np.float32)

    return np.clip(
        1.0 - (dists - in_rad) / (out_rad - in_rad),
        0.0,
        1.0,
    ).astype(np.float32)


def apply_spherical_window_mask_3d(
    vol: np.ndarray,
    *,
    in_rad: float = 0.85,
    out_rad: float = 0.99,
) -> np.ndarray:
    """Multiply a cubic volume by a 3-D spherical window."""
    vol = np.asarray(vol, dtype=np.float32)
    mask = spherical_window_mask_3d(D=int(vol.shape[0]), in_rad=in_rad, out_rad=out_rad)
    return (vol * mask).astype(np.float32)


def reconstruction_window_params(exp: DashboardExperiment) -> tuple[bool, float, float]:
    """Window settings from the training config (matches ``ImageDataset`` defaults)."""
    cfg = exp.train_configs or {}
    ds = cfg.get("dataset_args") or {}
    if "window" in ds:
        enabled = bool(ds["window"])
    else:
        enabled = True
    if not enabled:
        return False, 0.85, _RECON_WINDOW_OUT_RAD
    if ds.get("window_r") is not None:
        in_rad = float(ds["window_r"])
    elif cfg.get("window_radius_gt_real") is not None:
        in_rad = float(cfg["window_radius_gt_real"])
    else:
        in_rad = 0.85
    return True, in_rad, _RECON_WINDOW_OUT_RAD


def apply_reconstruction_window(
    vol: np.ndarray, exp: DashboardExperiment
) -> np.ndarray:
    """Apply the same real-space spherical window used during reconstruction."""
    enabled, in_rad, out_rad = reconstruction_window_params(exp)
    if not enabled:
        return np.asarray(vol, dtype=np.float32)
    return apply_spherical_window_mask_3d(vol, in_rad=in_rad, out_rad=out_rad)


def _decode_volume_array(exp: DashboardExperiment, row: int) -> np.ndarray:
    row = int(row)
    n_df = len(exp.plot_df)
    if row < 0 or row >= n_df:
        raise ValueError("Row index out of range.")
    mrc_dir = tempfile.mkdtemp(prefix="cryodrgn_vslice_mrc_")
    try:
        zsel = exp.z[np.asarray([row], dtype=int)]
        paths = _decode_z_values_to_vol_paths(exp, zsel, mrc_dir)
        vol, _ = parse_mrc(paths[0])
        return apply_reconstruction_window(np.asarray(vol, dtype=np.float32), exp)
    finally:
        shutil.rmtree(mrc_dir, ignore_errors=True)


def volume_array_b64(vol: np.ndarray) -> str:
    """Encode a cubic volume as base64 float32 bytes for the browser canvas."""
    vol = np.asarray(vol, dtype=np.float32)
    if vol.ndim != 3:
        raise ValueError("Volume must be 3-D.")
    return base64.standard_b64encode(vol.tobytes()).decode("ascii")


def downsample_volume_box_average(
    vol: np.ndarray,
    target_d: int | None = None,
) -> np.ndarray:
    """Box-average downsample to ``target_d``³ (matches ``volume_3d_utils.js``)."""
    vol = np.asarray(vol, dtype=np.float32)
    src_d = int(vol.shape[0])
    if target_d is None:
        target_d = PLOT3D_TARGET_D
    target_d = int(target_d)
    if src_d == target_d:
        return vol
    if target_d < 1:
        raise ValueError("target_d must be positive.")
    scale = src_d / target_d
    out = np.zeros((target_d, target_d, target_d), dtype=np.float32)
    for iz in range(target_d):
        z0 = int(np.floor(iz * scale))
        z1 = int(min(src_d, np.ceil((iz + 1) * scale)))
        for iy in range(target_d):
            y0 = int(np.floor(iy * scale))
            y1 = int(min(src_d, np.ceil((iy + 1) * scale)))
            for ix in range(target_d):
                x0 = int(np.floor(ix * scale))
                x1 = int(min(src_d, np.ceil((ix + 1) * scale)))
                block = vol[x0:x1, y0:y1, z0:z1]
                out[ix, iy, iz] = float(block.mean()) if block.size else 0.0
    return out


def prepare_vtk_transfer_volume(
    vol: np.ndarray,
    *,
    target_d: int | None = None,
) -> tuple[np.ndarray, dict]:
    """Downsample large cubes before browser transfer (matches client VTK prefilter)."""
    vol = np.asarray(vol, dtype=np.float32)
    source_d = int(vol.shape[0])
    if target_d is None:
        target_d = PLOT3D_TARGET_D
    target_d = int(target_d)
    if source_d > target_d:
        vol = downsample_volume_box_average(vol, target_d)
    d = int(vol.shape[0])
    meta = {
        "D": d,
        "source_D": source_d,
        "downsample": "box_average" if source_d > d else "none",
    }
    return vol, meta


def vtk_transfer_volume_payload(
    vol: np.ndarray,
    *,
    target_d: int | None = None,
) -> dict:
    """Float32 volume blob metadata for client VTK / slice backends."""
    vol, meta = prepare_vtk_transfer_volume(vol, target_d=target_d)
    return {
        **meta,
        "volume_b64": volume_array_b64(vol),
        "volume_dtype": "float32",
    }


def decode_volume_payload(exp: DashboardExperiment, row: int) -> dict:
    """Decode one particle volume and return array data for client-side slicing."""
    vol = _decode_volume_array(exp, row)
    transfer = vtk_transfer_volume_payload(vol)
    ds_idx = int(exp.all_indices[int(row)])
    return {
        "ok": True,
        "row": int(row),
        **transfer,
        "dataset_index": ds_idx,
        "default_slice_z": int(transfer["D"]) // 2,
    }


_ANALYZE_VOL_CACHE: dict[tuple[str, int, int, int], dict] = {}
_ANALYZE_VOL_ENTRY_CACHE: dict[tuple[str, int, int, str], dict] = {}
_ANALYZE_MARKER_CACHE: dict[tuple[str, int, int], list[dict]] = {}
_ANALYZE_CATALOG_CACHE: dict[tuple[str, int, int], list[dict]] = {}
_ANALYZE_VOL_CACHE_LOCK = threading.Lock()


def _analyze_dir(exp: DashboardExperiment) -> str:
    return os.path.join(exp.workdir, f"analyze.{int(exp.epoch)}")


def _sorted_vol_mrc_paths(directory: str) -> list[tuple[int, str]]:
    """Return ``(vol_index, path)`` for ``vol_NNN.mrc`` files under ``directory``."""
    if not os.path.isdir(directory):
        return []
    indexed: list[tuple[int, str]] = []
    with os.scandir(directory) as it:
        for entry in it:
            if not entry.is_file():
                continue
            m = _VOL_MRC_RE.match(entry.name)
            if m:
                indexed.append((int(m.group(1)), entry.path))
    indexed.sort(key=lambda t: t[0])
    return indexed


def discover_analyze_volume_catalog(exp: DashboardExperiment) -> list[dict]:
    """Catalog ``kmeans*`` and ``pc*`` volumes from ``analyze.{epoch}/``."""
    cache_key = (exp.workdir, int(exp.epoch), int(exp.kmeans_folder_id))
    with _ANALYZE_VOL_CACHE_LOCK:
        cached = _ANALYZE_CATALOG_CACHE.get(cache_key)
        if cached is not None:
            return [dict(e) for e in cached]

    anlz = _analyze_dir(exp)
    if not os.path.isdir(anlz):
        return []

    catalog: list[dict] = []
    km_dir = os.path.join(anlz, f"kmeans{int(exp.kmeans_folder_id)}")
    for vol_index, path in _sorted_vol_mrc_paths(km_dir):
        cluster_label = vol_index - 1
        display_cluster = cluster_label + 1
        catalog.append(
            {
                "id": f"kmeans:{cluster_label}",
                "kind": "kmeans",
                "title": f"k-means cluster {display_cluster}",
                "cluster_label": int(cluster_label),
                "vol_index": int(vol_index),
                "path": path,
            }
        )

    with os.scandir(anlz) as it:
        pc_dirs = sorted(
            (
                entry.name
                for entry in it
                if entry.is_dir() and _PC_DIR_RE.fullmatch(entry.name)
            ),
            key=lambda name: int(_PC_DIR_RE.fullmatch(name).group(1)),  # type: ignore[union-attr]
        )
    for pc_name in pc_dirs:
        m_pc = _PC_DIR_RE.fullmatch(pc_name)
        assert m_pc is not None
        pc_num = int(m_pc.group(1))
        pc_dir = os.path.join(anlz, pc_name)
        paths = _sorted_vol_mrc_paths(pc_dir)
        if not paths:
            continue
        vol_indices = [vi for vi, _ in paths]
        vol_start = min(vol_indices)
        for vol_index, path in paths:
            sample_idx = vol_index - vol_start
            catalog.append(
                {
                    "id": f"pc{pc_num}:{sample_idx}",
                    "kind": "pc",
                    "title": f"PC{pc_num} sample {sample_idx + 1}/{len(paths)}",
                    "pc": pc_num,
                    "sample_index": int(sample_idx),
                    "vol_index": int(vol_index),
                    "path": path,
                }
            )
    with _ANALYZE_VOL_CACHE_LOCK:
        _ANALYZE_CATALOG_CACHE[cache_key] = [dict(e) for e in catalog]
    return catalog


def _load_pc_trajectory_z_values(pc_dir: str) -> np.ndarray | None:
    """Load ``pc*/z_values.txt`` as ``(n, zdim)``, or ``None`` if missing."""
    path = os.path.join(pc_dir, "z_values.txt")
    if not os.path.isfile(path):
        return None
    try:
        arr = np.loadtxt(path)
    except (OSError, ValueError):
        return None
    arr = np.asarray(arr, dtype=np.float64)
    if arr.size == 0:
        return None
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.ndim != 2:
        return None
    return arr


def _pc_trajectory_nearest_plot_rows(
    exp: DashboardExperiment,
    z_pc: np.ndarray,
) -> list[int]:
    """Nearest on-data plot rows for PC trajectory ``z`` samples (UMAP fallback)."""
    z_pc = np.asarray(z_pc, dtype=np.float32)
    if z_pc.ndim != 2 or z_pc.shape[0] < 1:
        return []
    z = np.asarray(exp.z, dtype=np.float32)
    if z_pc.shape[1] != z.shape[1]:
        return []
    _, pc_ind = cryo_analysis.get_nearest_point(z, z_pc)
    return [int(x) for x in np.atleast_1d(pc_ind).tolist()]


def _plot_axis_values_from_latent_z(
    exp: DashboardExperiment,
    z_mat: np.ndarray,
    col: str,
    *,
    pca=None,
) -> np.ndarray | None:
    """Map latent ``z`` rows onto one scatter axis (``z*`` or ``PC*`` only)."""
    z_mat = np.atleast_2d(np.asarray(z_mat, dtype=np.float64))
    if z_mat.ndim != 2 or z_mat.shape[0] < 1:
        return None
    zdim = int(exp.z.shape[1])
    if z_mat.shape[1] != zdim:
        return None
    col = str(col)
    if re.fullmatch(r"z[0-9]+", col):
        idx = int(col[1:])
        if 0 <= idx < zdim:
            return z_mat[:, idx].astype(np.float64, copy=False)
        return None
    m_pc = re.fullmatch(r"PC([0-9]+)", col)
    if m_pc:
        pc_idx = int(m_pc.group(1)) - 1
        if pca is None:
            _, pca = cryo_analysis.run_pca(np.asarray(exp.z, dtype=np.float32))
        pc_coords = np.asarray(pca.transform(z_mat), dtype=np.float64)
        if 0 <= pc_idx < pc_coords.shape[1]:
            return pc_coords[:, pc_idx]
        return None
    return None


def project_latent_z_to_plot_xy(
    exp: DashboardExperiment,
    z_mat: np.ndarray,
    xcol: str,
    ycol: str,
) -> np.ndarray | None:
    """Project latent ``z`` rows into ``(xcol, ycol)`` plot space.

    Supports ``z*`` and ``PC*`` axes (matching analyze PCA traversal plots).
    Returns ``None`` for axes that require out-of-sample embedding (e.g. UMAP).
    """
    xcol = str(xcol)
    ycol = str(ycol)
    pca = None
    if re.fullmatch(r"PC([0-9]+)", xcol) or re.fullmatch(r"PC([0-9]+)", ycol):
        _, pca = cryo_analysis.run_pca(np.asarray(exp.z, dtype=np.float32))
    xs = _plot_axis_values_from_latent_z(exp, z_mat, xcol, pca=pca)
    ys = _plot_axis_values_from_latent_z(exp, z_mat, ycol, pca=pca)
    if xs is None or ys is None:
        return None
    return np.column_stack([xs, ys])


def enrich_analyze_volume_markers_plot_xy(
    exp: DashboardExperiment,
    markers: list[dict],
    xcol: str,
    ycol: str,
) -> list[dict]:
    """Attach ``xy`` plot coordinates for the current scatter axes (in-place copies)."""
    xcol = str(xcol)
    ycol = str(ycol)
    coords = exp.plot_df[[xcol, ycol]].values.astype(np.float64)
    out: list[dict] = []
    z_rows: list[np.ndarray] = []
    z_marker_indices: list[int] = []
    for marker in markers:
        m = dict(marker)
        m.pop("xy", None)
        z_raw = m.get("z")
        if z_raw is not None:
            try:
                z_arr = np.asarray(z_raw, dtype=np.float64).reshape(-1)
            except (TypeError, ValueError):
                z_arr = np.asarray([], dtype=np.float64)
            if z_arr.size == int(exp.z.shape[1]) and np.all(np.isfinite(z_arr)):
                z_marker_indices.append(len(out))
                z_rows.append(z_arr)
        plot_row = m.get("plot_row")
        if plot_row is not None:
            ri = int(plot_row)
            if 0 <= ri < len(coords) and np.all(np.isfinite(coords[ri])):
                m["xy"] = [float(coords[ri, 0]), float(coords[ri, 1])]
        out.append(m)

    if z_rows:
        z_mat = np.vstack(z_rows)
        projected = project_latent_z_to_plot_xy(exp, z_mat, xcol, ycol)
        if projected is not None:
            for local_i, marker_i in enumerate(z_marker_indices):
                pt = projected[local_i]
                if np.all(np.isfinite(pt)):
                    # Prefer true PC/z trajectory coordinates over nearest-particle xy.
                    out[marker_i]["xy"] = [float(pt[0]), float(pt[1])]
    return out


def _load_kmeans_center_plot_rows(km_dir: str) -> list[int] | None:
    path = os.path.join(km_dir, "centers_ind.txt")
    if not os.path.isfile(path):
        return None
    ind = np.loadtxt(path)
    return [int(x) for x in np.atleast_1d(ind).tolist()]


def discover_analyze_volume_markers(
    exp: DashboardExperiment,
    catalog: list[dict] | None = None,
) -> list[dict]:
    """Scatterplot anchor metadata for each analyze volume.

    PC markers carry latent ``z`` from ``analyze.*/pc*/z_values.txt`` so the
    trajectory creator can place them at the true PC-traversal coordinates
    (instead of snapping to a nearest particle). Volumes still load from the
    corresponding ``vol_*.mrc`` files. ``plot_row`` remains the nearest particle
    for UMAP axes and colour fallbacks.
    """
    cache_key = (exp.workdir, int(exp.epoch), int(exp.kmeans_folder_id))
    with _ANALYZE_VOL_CACHE_LOCK:
        cached = _ANALYZE_MARKER_CACHE.get(cache_key)
        if cached is not None:
            return [dict(m) for m in cached]

    if catalog is None:
        catalog = discover_analyze_volume_catalog(exp)

    km_dir = os.path.join(_analyze_dir(exp), f"kmeans{int(exp.kmeans_folder_id)}")
    km_rows = _load_kmeans_center_plot_rows(km_dir)

    pc_nums = sorted(
        {int(entry["pc"]) for entry in catalog if entry.get("kind") == "pc"}
    )
    pc_z_cache: dict[int, np.ndarray] = {}
    pc_row_cache: dict[int, list[int]] = {}
    anlz = _analyze_dir(exp)
    pca = None
    for pc_num in pc_nums:
        pc_dir = os.path.join(anlz, f"pc{pc_num}")
        z_pc = _load_pc_trajectory_z_values(pc_dir)
        n_samples = sum(
            1
            for entry in catalog
            if entry.get("kind") == "pc" and int(entry["pc"]) == pc_num
        )
        if z_pc is None or z_pc.shape[1] != int(exp.z.shape[1]):
            # Older analyze folders without z_values.txt: regenerate the PC
            # trajectory and snap to nearest particles.
            if n_samples < 1:
                continue
            if pca is None:
                _, pca = cryo_analysis.run_pca(np.asarray(exp.z, dtype=np.float32))
            pc = np.asarray(exp.pc, dtype=np.float64)
            col = int(pc_num) - 1
            if col < 0 or col >= pc.shape[1]:
                continue
            lo, hi = np.percentile(pc[:, col], (5, 95))
            z_pc = cryo_analysis.get_pc_traj(
                pca, int(exp.z.shape[1]), n_samples, pc_num, lo, hi
            )
        pc_z_cache[pc_num] = np.asarray(z_pc, dtype=np.float64)
        pc_row_cache[pc_num] = _pc_trajectory_nearest_plot_rows(exp, z_pc)

    markers: list[dict] = []
    for entry in catalog:
        vol_id = str(entry["id"])
        kind = entry.get("kind")
        plot_row: int | None = None
        label = vol_id
        z_list: list[float] | None = None
        if kind == "kmeans":
            cl = int(entry["cluster_label"])
            label = f"K{cl + 1}"
            if km_rows is not None and 0 <= cl < len(km_rows):
                plot_row = int(km_rows[cl])
        elif kind == "pc":
            pc_num = int(entry["pc"])
            si = int(entry["sample_index"])
            label = f"PC{pc_num}·{si + 1}"
            z_pc = pc_z_cache.get(pc_num)
            if z_pc is not None and 0 <= si < len(z_pc):
                z_list = [float(v) for v in z_pc[si].tolist()]
            rows = pc_row_cache.get(pc_num) or []
            if 0 <= si < len(rows):
                plot_row = int(rows[si])
        marker: dict = {
            "vol_id": vol_id,
            "plot_row": plot_row,
            "label": label,
            "kind": kind,
        }
        if z_list is not None:
            marker["z"] = z_list
        markers.append(marker)

    with _ANALYZE_VOL_CACHE_LOCK:
        _ANALYZE_MARKER_CACHE[cache_key] = markers
    return [dict(m) for m in markers]


def kmeans_volume_ids_for_anchor_indices(
    exp: DashboardExperiment, anchor_indices: list[int]
) -> list[str] | None:
    """Map path-ordered anchor rows to pre-decoded ``kmeans:*`` catalog ids."""
    if len(anchor_indices) < 2:
        return None
    markers = discover_analyze_volume_markers(exp)
    row_to_vol: dict[int, str] = {}
    for marker in markers:
        if marker.get("kind") != "kmeans":
            continue
        plot_row = marker.get("plot_row")
        if plot_row is None:
            continue
        row_to_vol[int(plot_row)] = str(marker["vol_id"])
    if not row_to_vol:
        return None
    vol_ids: list[str] = []
    for anchor in anchor_indices:
        vol_id = row_to_vol.get(int(anchor))
        if vol_id is None:
            return None
        vol_ids.append(vol_id)
    return vol_ids


def analyze_volume_marker_index(
    exp: DashboardExperiment,
    xcol: str,
    ycol: str,
    allowed_vol_ids: set[str] | frozenset[str] | None = None,
) -> tuple[np.ndarray, list[str], list[int]]:
    """Scatter coordinates, catalog ids, and plot rows for analyze volume markers."""
    markers = enrich_analyze_volume_markers_plot_xy(
        exp, discover_analyze_volume_markers(exp), xcol, ycol
    )
    xy_parts: list[np.ndarray] = []
    vol_ids: list[str] = []
    plot_rows: list[int] = []
    for marker in markers:
        vol_id = str(marker["vol_id"])
        if allowed_vol_ids is not None and vol_id not in allowed_vol_ids:
            continue
        xy = marker.get("xy")
        if (
            not isinstance(xy, (list, tuple))
            or len(xy) != 2
            or not np.all(np.isfinite(np.asarray(xy, dtype=np.float64)))
        ):
            continue
        plot_row = marker.get("plot_row")
        # Free PC latent points may lack a particle row; keep index alignment with
        # a non-negative sentinel only when a real row exists.
        if plot_row is None:
            continue
        ri = int(plot_row)
        xy_parts.append(np.asarray(xy, dtype=np.float64))
        vol_ids.append(vol_id)
        plot_rows.append(ri)
    if not xy_parts:
        raise ValueError(
            "No analyze volume markers with plot coordinates for the current axes."
        )
    return np.vstack(xy_parts), vol_ids, plot_rows


def snap_xy_to_analyze_volume_markers(
    exp: DashboardExperiment,
    xy_pts: np.ndarray,
    xcol: str,
    ycol: str,
    allowed_vol_ids: set[str] | frozenset[str] | None = None,
) -> tuple[list[int], np.ndarray, list[str]]:
    """Nearest analyze-volume marker for each ``(x, y)`` query in plot space."""
    marker_xy, vol_ids, plot_rows = analyze_volume_marker_index(
        exp, xcol, ycol, allowed_vol_ids=allowed_vol_ids
    )
    xy_pts = np.atleast_2d(np.asarray(xy_pts, dtype=np.float64))
    diff = xy_pts[:, None, :] - marker_xy[None, :, :]
    dist2 = np.sum(diff * diff, axis=2)
    nearest = np.argmin(dist2, axis=1)
    out_rows = [plot_rows[int(i)] for i in nearest]
    out_xy = marker_xy[nearest]
    out_vols = [vol_ids[int(i)] for i in nearest]
    return out_rows, out_xy, out_vols


def analyze_volume_ids_for_plot_rows(
    exp: DashboardExperiment, plot_rows: list[int]
) -> list[str | None]:
    """Map plot rows to catalog volume ids where each row is an analyze marker."""
    markers = discover_analyze_volume_markers(exp)
    row_to_vol: dict[int, str] = {}
    for marker in markers:
        plot_row = marker.get("plot_row")
        if plot_row is not None:
            row_to_vol[int(plot_row)] = str(marker["vol_id"])
    return [row_to_vol.get(int(row)) for row in plot_rows]


def _catalog_public(catalog: list[dict]) -> list[dict]:
    return [{k: v for k, v in entry.items() if k != "path"} for entry in catalog]


def _enrich_kmeans_catalog_znorm(exp: DashboardExperiment, catalog: list[dict]) -> None:
    """Attach ``znorm`` at each k-means center's nearest plot row (in-place)."""
    if "znorm" not in exp.plot_df.columns:
        return
    km_dir = os.path.join(_analyze_dir(exp), f"kmeans{int(exp.kmeans_folder_id)}")
    km_rows = _load_kmeans_center_plot_rows(km_dir)
    if km_rows is None:
        return
    znorm = np.asarray(exp.plot_df["znorm"], dtype=np.float64)
    for entry in catalog:
        if entry.get("kind") != "kmeans":
            continue
        cl = int(entry["cluster_label"])
        if not (0 <= cl < len(km_rows)):
            continue
        plot_row = int(km_rows[cl])
        if not (0 <= plot_row < len(znorm)):
            continue
        zn = float(znorm[plot_row])
        if np.isfinite(zn):
            entry["znorm"] = zn


def default_analyze_volume_id(catalog: list[dict]) -> str | None:
    """Return the k-means catalog id with the lowest ``znorm``, if any."""
    best_id: str | None = None
    best_znorm: float | None = None
    for entry in catalog:
        if entry.get("kind") != "kmeans":
            continue
        zn = entry.get("znorm")
        if zn is None:
            continue
        zn_f = float(zn)
        if not np.isfinite(zn_f):
            continue
        if best_znorm is None or zn_f < best_znorm:
            best_znorm = zn_f
            best_id = str(entry["id"])
    if best_id is not None:
        return best_id
    for entry in catalog:
        if entry.get("kind") == "kmeans":
            return str(entry["id"])
    return str(catalog[0]["id"]) if catalog else None


def _catalog_entry_by_id(catalog: list[dict], vol_id: str) -> dict | None:
    for entry in catalog:
        if entry["id"] == vol_id:
            return entry
    return None


def analyze_volumes_catalog_payload(
    exp: DashboardExperiment,
    *,
    include_markers: bool = True,
) -> dict:
    """Catalog of analyze k-means / PC volumes (no volume arrays)."""
    cache_key = (
        exp.workdir,
        int(exp.epoch),
        int(exp.kmeans_folder_id),
        -1 if include_markers else -2,
    )
    with _ANALYZE_VOL_CACHE_LOCK:
        cached = _ANALYZE_VOL_CACHE.get(cache_key)
        if cached is not None:
            return dict(cached["payload"])

    catalog = discover_analyze_volume_catalog(exp)
    catalog_public = _catalog_public(catalog)
    _enrich_kmeans_catalog_znorm(exp, catalog_public)
    markers: list[dict] = []
    if include_markers:
        markers = discover_analyze_volume_markers(exp, catalog)

    payload = {
        "ok": True,
        "catalog": catalog_public,
        "default_vol_id": default_analyze_volume_id(catalog_public),
        "markers": markers,
        "D": None,
    }
    with _ANALYZE_VOL_CACHE_LOCK:
        _ANALYZE_VOL_CACHE[cache_key] = {"payload": payload, "t0": time.monotonic()}
    return payload


def analyze_volume_markers_payload(
    exp: DashboardExperiment,
    *,
    xcol: str | None = None,
    ycol: str | None = None,
) -> dict:
    """Scatterplot anchor rows and labels for analyze volumes."""
    catalog = discover_analyze_volume_catalog(exp)
    markers = discover_analyze_volume_markers(exp, catalog)
    if xcol and ycol and xcol in exp.plot_df.columns and ycol in exp.plot_df.columns:
        markers = enrich_analyze_volume_markers_plot_xy(exp, markers, xcol, ycol)
    return {"ok": True, "markers": markers}


def analyze_volume_by_id_payload(
    exp: DashboardExperiment,
    vol_id: str,
    *,
    target_d: int | None = None,
) -> dict:
    """Load one analyze volume by catalog id (e.g. ``kmeans:0``, ``pc1:3``)."""
    vol_id = str(vol_id)
    if target_d is None:
        target_d = PLOT3D_TARGET_D
    cache_key = (
        exp.workdir,
        int(exp.epoch),
        int(exp.kmeans_folder_id),
        vol_id,
        int(target_d),
    )
    with _ANALYZE_VOL_CACHE_LOCK:
        cached = _ANALYZE_VOL_ENTRY_CACHE.get(cache_key)
        if cached is not None:
            return dict(cached["payload"])

    catalog = discover_analyze_volume_catalog(exp)
    entry = _catalog_entry_by_id(catalog, vol_id)
    if entry is None:
        raise ValueError(f"Unknown analyze volume id: {vol_id}")

    _, vol_payload = _load_catalog_volume(entry, exp, target_d=target_d)
    meta = next(e for e in _catalog_public(catalog) if e["id"] == vol_id)
    payload = {
        "ok": True,
        "id": vol_id,
        "meta": meta,
        **vol_payload,
    }
    with _ANALYZE_VOL_CACHE_LOCK:
        _ANALYZE_VOL_ENTRY_CACHE[cache_key] = {
            "payload": payload,
            "t0": time.monotonic(),
        }
    return payload


def analyze_volumes_batch_payload(
    exp: DashboardExperiment,
    vol_ids: list[str],
    *,
    n_cpus: int = 1,
    target_d: int | None = None,
) -> dict:
    """Load several analyze volumes in parallel (for background prefetch)."""
    if target_d is None:
        target_d = PLOT3D_TARGET_D
    catalog = discover_analyze_volume_catalog(exp)
    by_id = {entry["id"]: entry for entry in catalog}
    entries = []
    for vol_id in vol_ids:
        entry = by_id.get(str(vol_id))
        if entry is not None:
            entries.append(entry)
    volumes = _load_catalog_volumes_parallel(
        entries, exp=exp, n_cpus=n_cpus, target_d=target_d
    )
    return {"ok": True, "volumes": volumes, "target_d": int(target_d)}


def _png_is_mostly_blank(path: str, *, min_std: float = 0.5) -> bool:
    """True when a PNG is essentially uniform (failed/off-screen ChimeraX render)."""
    try:
        from PIL import Image

        with Image.open(path) as im:
            gray = np.asarray(im.convert("L"), dtype=np.float32)
        if gray.size == 0:
            return True
        return float(gray.std()) < min_std
    except Exception:
        return True


def parse_iso_level_from_request(data: dict) -> float | None:
    """Parse optional ``iso_level`` (map data units) from dashboard JSON."""
    if data.get("iso_level") is None:
        return None
    try:
        level = float(data["iso_level"])
    except (TypeError, ValueError) as err:
        raise ValueError("iso_level must be a finite number.") from err
    if not np.isfinite(level):
        raise ValueError("iso_level must be a finite number.")
    return level


def parse_chimerax_view_turns_from_request(raw_turns) -> list[tuple[str, float]] | None:
    """Parse ``view_turns`` JSON from dashboard API requests."""
    if raw_turns is None:
        return None
    if not isinstance(raw_turns, list):
        raise ValueError("view_turns must be a list of {axis, degrees} objects.")
    parsed: list[tuple[str, float]] = []
    for item in raw_turns:
        if not isinstance(item, dict):
            raise ValueError(
                "view_turns entries must be objects with axis and degrees."
            )
        parsed.append((str(item.get("axis", "")), float(item.get("degrees", 0))))
    return normalize_chimerax_view_turns(parsed)


def analyze_volumes_chimerax_batch_payload(
    exp: DashboardExperiment,
    vol_ids: list[str],
    *,
    chimerax_cpus: int = DEFAULT_CHIMERAX_PARALLEL,
    view_matrix_camera: str | None = None,
    view_turns: list[tuple[str, float]] | None = None,
    volume_level: float | None = None,
) -> dict:
    """Render ChimeraX PNGs from pre-generated analyze ``.mrc`` files (no GPU decode)."""
    catalog = discover_analyze_volume_catalog(exp)
    by_id = {entry["id"]: entry for entry in catalog}
    entries = []
    for vol_id in vol_ids:
        entry = by_id.get(str(vol_id))
        if entry is None:
            raise ValueError(f"Unknown analyze volume id: {vol_id}")
        entries.append(entry)
    if not entries:
        raise ValueError("ids must be a non-empty list of volume ids.")

    level = resolve_chimerax_volume_level(entries[0]["path"], volume_level)
    cc = max(1, min(int(chimerax_cpus), 32))
    attempts: list[tuple[str | None, list[tuple[str, float]] | None]] = []
    if view_turns:
        attempts.append((None, list(view_turns)))
    elif view_matrix_camera:
        attempts.append((view_matrix_camera, None))
        attempts.append((None, None))
    else:
        attempts.append((None, None))
    last_blank_err = (
        "ChimeraX produced blank volume images. "
        "Check CHIMERAX_PATH and that analyze .mrc files are valid."
    )
    for attempt_idx, (attempt_vm, attempt_turns) in enumerate(attempts):
        view_matrix_text: str | None = None
        with tempfile.TemporaryDirectory(prefix="cryodrgn_analyze_cx_png_") as png_dir:
            views = [
                LandscapeStaticView(
                    mrc_path=entry["path"],
                    out_png=os.path.join(png_dir, f"cell_{i}.png"),
                    volume_level=level,
                    view_turns=attempt_turns,
                    view_matrix_camera=attempt_vm,
                    report_view_matrix=(i == 0),
                )
                for i, entry in enumerate(entries)
            ]
            paths, view_matrix_text = render_landscape_cycle_static_views(
                views, chimerax_cpus=cc
            )
            images: list[str] = []
            blank_count = 0
            for pth in paths:
                if _png_is_mostly_blank(pth):
                    blank_count += 1
                with open(pth, "rb") as fh:
                    images.append(base64.standard_b64encode(fh.read()).decode("ascii"))
            if blank_count == len(paths):
                if attempt_idx + 1 < len(attempts):
                    continue
                raise ValueError(last_blank_err)
        out: dict = {
            "ok": True,
            "images": images,
            "ids": [str(v) for v in vol_ids],
        }
        out.update(
            chimerax_iso_response_fields(entries[0]["path"], iso_level=volume_level)
        )
        if view_matrix_text:
            out["view_matrix"] = view_matrix_text
        if attempt_turns:
            out["view_turns"] = [
                {"axis": axis, "degrees": degrees} for axis, degrees in attempt_turns
            ]
        return out
    raise ValueError(last_blank_err)


def _load_catalog_volume(
    entry: dict,
    exp: DashboardExperiment,
    *,
    target_d: int | None = None,
) -> tuple[str, dict]:
    """Load one analyze ``.mrc`` volume (module-level for process pools)."""
    vol, _ = parse_mrc(entry["path"])
    vol = apply_reconstruction_window(np.asarray(vol, dtype=np.float32), exp)
    return entry["id"], vtk_transfer_volume_payload(vol, target_d=target_d)


def _load_catalog_volumes_parallel(
    catalog: list[dict],
    *,
    exp: DashboardExperiment,
    n_cpus: int,
    target_d: int | None = None,
) -> dict[str, dict]:
    """Load catalog volumes using up to ``n_cpus`` worker threads."""
    n_cpus = max(1, int(n_cpus))
    if len(catalog) <= 1:
        volumes: dict[str, dict] = {}
        for entry in catalog:
            vol_id, payload = _load_catalog_volume(entry, exp, target_d=target_d)
            volumes[vol_id] = payload
        return volumes

    max_workers = min(n_cpus, len(catalog))
    volumes = {}
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = [
            pool.submit(_load_catalog_volume, entry, exp, target_d=target_d)
            for entry in catalog
        ]
        for fut in as_completed(futures):
            vol_id, payload = fut.result()
            volumes[vol_id] = payload
    return volumes
