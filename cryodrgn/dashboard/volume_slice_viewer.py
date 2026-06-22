"""Orthogonal volume slice previews for the dashboard volume slice viewer.

Decodes latent ``z`` vectors to 3-D volumes (reusing the particle explorer decoder)
and renders matplotlib slice panels similar to
:func:`cryodrgn.commands.analyze_landscape.view_slices`.
"""

from __future__ import annotations

import base64
import io
import os
import re
import secrets
import shutil
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import matplotlib
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np
from scipy.ndimage import affine_transform
from scipy.spatial.transform import Rotation

from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.particle_explorer import _decode_z_values_to_vol_paths
from cryodrgn.mrcfile import parse_mrc
from cryodrgn import analysis as cryo_analysis

matplotlib.use("Agg")

_VOL_MRC_RE = re.compile(r"^vol_(\d+)\.mrc$", re.IGNORECASE)
_PC_DIR_RE = re.compile(r"^pc(\d+)$", re.IGNORECASE)

_SLICE_CACHE: dict[str, dict] = {}
_SLICE_CACHE_LOCK = threading.Lock()
_SLICE_CACHE_MAX_ENTRIES = 16
_SLICE_CACHE_TTL_S = 7200.0

_SLICE_TITLES = ("Y–Z (fixed X)", "X–Z (fixed Y)", "X–Y (fixed Z)")


def rotate_volume_array(
    vol: np.ndarray,
    rot_x_deg: float,
    rot_y_deg: float,
    rot_z_deg: float,
) -> np.ndarray:
    """Rotate a cubic volume by Euler angles (degrees, intrinsic XYZ)."""
    vol = np.asarray(vol, dtype=np.float32)
    if vol.ndim != 3 or vol.shape[0] != vol.shape[1] or vol.shape[1] != vol.shape[2]:
        raise ValueError("Volume must be a cubic 3-D array.")
    rot = Rotation.from_euler(
        "xyz",
        [float(rot_x_deg), float(rot_y_deg), float(rot_z_deg)],
        degrees=True,
    )
    matrix = rot.as_matrix()
    center = (np.array(vol.shape, dtype=np.float64) - 1.0) / 2.0
    offset = center - matrix @ center
    fill = float(np.min(vol))
    return affine_transform(
        vol,
        matrix,
        offset=offset,
        order=1,
        mode="constant",
        cval=fill,
    ).astype(np.float32, copy=False)


def _slice_indices(
    d: int,
    slice_ix: int | None,
    slice_iy: int | None,
    slice_iz: int | None,
) -> tuple[int, int, int]:
    d = int(d)
    if d < 1:
        raise ValueError("Volume dimension must be positive.")
    mid = d // 2

    def _one(raw: int | None, name: str) -> int:
        if raw is None:
            return mid
        idx = int(raw)
        if idx < 0 or idx >= d:
            raise ValueError(f"{name} slice index must be in [0, {d - 1}].")
        return idx

    return _one(slice_ix, "X"), _one(slice_iy, "Y"), _one(slice_iz, "Z")


def orthogonal_slice_png_b64_list(
    vol: np.ndarray,
    *,
    slice_ix: int | None = None,
    slice_iy: int | None = None,
    slice_iz: int | None = None,
) -> list[str]:
    """Return three separate PNG panels (one per orthogonal axis) as base64 strings."""
    vol = np.asarray(vol)
    if vol.ndim != 3:
        raise ValueError("Volume must be 3-D.")
    d = int(vol.shape[0])
    ix, iy, iz = _slice_indices(d, slice_ix, slice_iy, slice_iz)
    panels = (vol[ix, :, :], vol[:, iy, :], vol[:, :, iz])
    out: list[str] = []
    for data, title in zip(panels, _SLICE_TITLES, strict=True):
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        ax.imshow(data, cmap="gray", aspect="equal")
        ax.set_title(title, fontsize=9)
        ax.axis("off")
        fig.tight_layout(pad=0.2)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=110, facecolor="white")
        plt.close(fig)
        out.append(base64.standard_b64encode(buf.getvalue()).decode("ascii"))
    return out


def _slice_cache_evict_unlocked(token: str) -> None:
    _SLICE_CACHE.pop(token, None)


def _slice_cache_prune_unlocked() -> None:
    now = time.monotonic()
    dead = [
        tok
        for tok, meta in _SLICE_CACHE.items()
        if now - meta["t0"] > _SLICE_CACHE_TTL_S
    ]
    for tok in dead:
        _slice_cache_evict_unlocked(tok)
    while len(_SLICE_CACHE) >= _SLICE_CACHE_MAX_ENTRIES:
        oldest = min(_SLICE_CACHE.items(), key=lambda kv: kv[1]["t0"])[0]
        _slice_cache_evict_unlocked(oldest)


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
        return np.asarray(vol, dtype=np.float32)
    finally:
        shutil.rmtree(mrc_dir, ignore_errors=True)


def decode_volume_to_cache(exp: DashboardExperiment, row: int) -> tuple[str, int]:
    """Decode one particle volume and register it in the in-memory slice cache."""
    vol = _decode_volume_array(exp, row)
    with _SLICE_CACHE_LOCK:
        _slice_cache_prune_unlocked()
        token = secrets.token_urlsafe(24)
        _SLICE_CACHE[token] = {
            "vol": vol,
            "row": int(row),
            "t0": time.monotonic(),
        }
        return token, int(vol.shape[0])


def slices_from_cache_payload(
    token: str,
    row_expected: int,
    *,
    rot_x_deg: float = 0.0,
    rot_y_deg: float = 0.0,
    rot_z_deg: float = 0.0,
    slice_ix: int | None = None,
    slice_iy: int | None = None,
    slice_iz: int | None = None,
) -> dict:
    """Slice images (base64 PNG) for a cached decoded volume."""
    with _SLICE_CACHE_LOCK:
        meta = _SLICE_CACHE.get(token)
        if not meta:
            raise ValueError("Unknown or expired volume cache id.")
        if time.monotonic() - meta["t0"] > _SLICE_CACHE_TTL_S:
            _slice_cache_evict_unlocked(token)
            raise ValueError("Volume cache expired. Decode the volume again.")
        if int(meta["row"]) != int(row_expected):
            raise ValueError("Cached volume row does not match request.")
        vol = np.array(meta["vol"], copy=True)

    if rot_x_deg or rot_y_deg or rot_z_deg:
        vol = rotate_volume_array(vol, rot_x_deg, rot_y_deg, rot_z_deg)

    images = orthogonal_slice_png_b64_list(
        vol,
        slice_ix=slice_ix,
        slice_iy=slice_iy,
        slice_iz=slice_iz,
    )
    d = int(vol.shape[0])
    ix, iy, iz = _slice_indices(d, slice_ix, slice_iy, slice_iz)
    return {
        "ok": True,
        "row": int(row_expected),
        "D": d,
        "slice_ix": ix,
        "slice_iy": iy,
        "slice_iz": iz,
        "rot_x_deg": float(rot_x_deg),
        "rot_y_deg": float(rot_y_deg),
        "rot_z_deg": float(rot_z_deg),
        "images": images,
        "labels": list(_SLICE_TITLES),
    }


def volume_array_b64(vol: np.ndarray) -> str:
    """Encode a cubic volume as base64 float32 bytes for the browser canvas."""
    vol = np.asarray(vol, dtype=np.float32)
    if vol.ndim != 3:
        raise ValueError("Volume must be 3-D.")
    return base64.standard_b64encode(vol.tobytes()).decode("ascii")


def decode_volume_payload(exp: DashboardExperiment, row: int) -> dict:
    """Decode one particle volume and return array data for client-side slicing."""
    vol = _decode_volume_array(exp, row)
    d = int(vol.shape[0])
    ds_idx = int(exp.all_indices[int(row)])
    return {
        "ok": True,
        "row": int(row),
        "D": d,
        "volume_b64": volume_array_b64(vol),
        "volume_dtype": "float32",
        "dataset_index": ds_idx,
        "default_slice_z": d // 2,
    }


def decode_and_initial_slices(exp: DashboardExperiment, row: int) -> dict:
    """Decode a volume for the interactive slice canvas (legacy name)."""
    return decode_volume_payload(exp, row)


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


def _pc_trajectory_plot_rows_with_pca(
    exp: DashboardExperiment,
    pca,
    pc_num: int,
    n_samples: int,
) -> list[int]:
    """Nearest on-data plot rows for ``pc{pc_num}/`` trajectory samples."""
    if n_samples < 1:
        return []
    z = np.asarray(exp.z, dtype=np.float32)
    pc = np.asarray(exp.pc, dtype=np.float64)
    col = int(pc_num) - 1
    if col < 0 or col >= pc.shape[1]:
        return []
    lo, hi = np.percentile(pc[:, col], (5, 95))
    z_pc = cryo_analysis.get_pc_traj(pca, z.shape[1], n_samples, pc_num, lo, hi)
    _, pc_ind = cryo_analysis.get_nearest_point(z, z_pc)
    return [int(x) for x in np.atleast_1d(pc_ind).tolist()]


def _pc_trajectory_plot_rows(
    exp: DashboardExperiment, pc_num: int, n_samples: int
) -> list[int]:
    z = np.asarray(exp.z, dtype=np.float32)
    _, pca = cryo_analysis.run_pca(z)
    return _pc_trajectory_plot_rows_with_pca(exp, pca, pc_num, n_samples)


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
    """Scatterplot anchor row and short label for each analyze volume."""
    cache_key = (exp.workdir, int(exp.epoch), int(exp.kmeans_folder_id))
    with _ANALYZE_VOL_CACHE_LOCK:
        cached = _ANALYZE_MARKER_CACHE.get(cache_key)
        if cached is not None:
            return [dict(m) for m in cached]

    if catalog is None:
        catalog = discover_analyze_volume_catalog(exp)

    km_dir = os.path.join(_analyze_dir(exp), f"kmeans{int(exp.kmeans_folder_id)}")
    km_rows = _load_kmeans_center_plot_rows(km_dir)

    pc_counts: dict[int, int] = {}
    for entry in catalog:
        if entry.get("kind") == "pc":
            pc_counts[int(entry["pc"])] = pc_counts.get(int(entry["pc"]), 0) + 1
    pc_row_cache: dict[int, list[int]] = {}
    if pc_counts:
        z = np.asarray(exp.z, dtype=np.float32)
        _, pca = cryo_analysis.run_pca(z)
        for pc_num, n_samples in pc_counts.items():
            pc_row_cache[pc_num] = _pc_trajectory_plot_rows_with_pca(
                exp, pca, pc_num, n_samples
            )

    markers: list[dict] = []
    for entry in catalog:
        vol_id = str(entry["id"])
        kind = entry.get("kind")
        plot_row: int | None = None
        label = vol_id
        if kind == "kmeans":
            cl = int(entry["cluster_label"])
            label = f"K{cl + 1}"
            if km_rows is not None and 0 <= cl < len(km_rows):
                plot_row = int(km_rows[cl])
        elif kind == "pc":
            pc_num = int(entry["pc"])
            si = int(entry["sample_index"])
            label = f"PC{pc_num}·{si + 1}"
            rows = pc_row_cache.get(pc_num) or []
            if 0 <= si < len(rows):
                plot_row = int(rows[si])
        markers.append(
            {
                "vol_id": vol_id,
                "plot_row": plot_row,
                "label": label,
                "kind": kind,
            }
        )

    with _ANALYZE_VOL_CACHE_LOCK:
        _ANALYZE_MARKER_CACHE[cache_key] = markers
    return [dict(m) for m in markers]


def _catalog_public(catalog: list[dict]) -> list[dict]:
    return [{k: v for k, v in entry.items() if k != "path"} for entry in catalog]


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
    markers: list[dict] = []
    if include_markers:
        markers = discover_analyze_volume_markers(exp, catalog)

    payload = {
        "ok": True,
        "catalog": _catalog_public(catalog),
        "markers": markers,
        "D": None,
    }
    with _ANALYZE_VOL_CACHE_LOCK:
        _ANALYZE_VOL_CACHE[cache_key] = {"payload": payload, "t0": time.monotonic()}
    return payload


def analyze_volume_markers_payload(exp: DashboardExperiment) -> dict:
    """Scatterplot anchor rows and labels for analyze volumes."""
    catalog = discover_analyze_volume_catalog(exp)
    markers = discover_analyze_volume_markers(exp, catalog)
    return {"ok": True, "markers": markers}


def analyze_volume_by_id_payload(exp: DashboardExperiment, vol_id: str) -> dict:
    """Load one analyze volume by catalog id (e.g. ``kmeans:0``, ``pc1:3``)."""
    vol_id = str(vol_id)
    cache_key = (exp.workdir, int(exp.epoch), int(exp.kmeans_folder_id), vol_id)
    with _ANALYZE_VOL_CACHE_LOCK:
        cached = _ANALYZE_VOL_ENTRY_CACHE.get(cache_key)
        if cached is not None:
            return dict(cached["payload"])

    catalog = discover_analyze_volume_catalog(exp)
    entry = _catalog_entry_by_id(catalog, vol_id)
    if entry is None:
        raise ValueError(f"Unknown analyze volume id: {vol_id}")

    _, vol_payload = _load_catalog_volume(entry)
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
) -> dict:
    """Load several analyze volumes in parallel (for background prefetch)."""
    catalog = discover_analyze_volume_catalog(exp)
    by_id = {entry["id"]: entry for entry in catalog}
    entries = []
    for vol_id in vol_ids:
        entry = by_id.get(str(vol_id))
        if entry is not None:
            entries.append(entry)
    volumes = _load_catalog_volumes_parallel(entries, n_cpus=n_cpus)
    return {"ok": True, "volumes": volumes}


def analyze_volumes_payload(
    exp: DashboardExperiment,
    *,
    n_cpus: int = 1,
) -> dict:
    """Backward-compatible alias: catalog only (volumes load on demand)."""
    del n_cpus
    payload = analyze_volumes_catalog_payload(exp)
    payload["volumes"] = {}
    return payload


def _load_catalog_volume(entry: dict) -> tuple[str, dict]:
    """Load one analyze ``.mrc`` volume (module-level for process pools)."""
    vol, _ = parse_mrc(entry["path"])
    vol = np.asarray(vol, dtype=np.float32)
    return entry["id"], {
        "D": int(vol.shape[0]),
        "volume_b64": volume_array_b64(vol),
    }


def _load_catalog_volumes_parallel(
    catalog: list[dict],
    *,
    n_cpus: int,
) -> dict[str, dict]:
    """Load catalog volumes using up to ``n_cpus`` worker threads."""
    n_cpus = max(1, int(n_cpus))
    if len(catalog) <= 1:
        volumes: dict[str, dict] = {}
        for entry in catalog:
            vol_id, payload = _load_catalog_volume(entry)
            volumes[vol_id] = payload
        return volumes

    max_workers = min(n_cpus, len(catalog))
    volumes = {}
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = [pool.submit(_load_catalog_volume, entry) for entry in catalog]
        for fut in as_completed(futures):
            vol_id, payload = fut.result()
            volumes[vol_id] = payload
    return volumes
