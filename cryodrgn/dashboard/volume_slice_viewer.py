"""Volume decode and analyze-volume catalog for the dashboard volume slice viewer.

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

import numpy as np

from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.particle_explorer import _decode_z_values_to_vol_paths
from cryodrgn.mrcfile import parse_mrc
from cryodrgn import analysis as cryo_analysis

_VOL_MRC_RE = re.compile(r"^vol_(\d+)\.mrc$", re.IGNORECASE)
_PC_DIR_RE = re.compile(r"^pc(\d+)$", re.IGNORECASE)


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
