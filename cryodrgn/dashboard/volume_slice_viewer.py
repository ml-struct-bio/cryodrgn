"""Orthogonal volume slice previews for the dashboard volume slice viewer.

Decodes latent ``z`` vectors to 3-D volumes (reusing the particle explorer decoder)
and renders matplotlib slice panels similar to
:func:`cryodrgn.commands.analyze_landscape.view_slices`.
"""

from __future__ import annotations

import base64
import io
import secrets
import shutil
import tempfile
import threading
import time
import numpy as np
import matplotlib

import matplotlib.pyplot as plt  # noqa: E402
from scipy.ndimage import affine_transform
from scipy.spatial.transform import Rotation

from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.particle_explorer import _decode_z_values_to_vol_paths
from cryodrgn.mrcfile import parse_mrc

matplotlib.use("Agg")

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
