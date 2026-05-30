"""On-demand decoder volumes + ChimeraX static PNGs and rotating GIFs.

For the particle explorer. ChimeraX rendering lives in
:mod:`cryodrgn.dashboard.chimerax_animation` (shared with landscape vol-PCA and
trajectory creator).
"""

from __future__ import annotations

import os
import re
import secrets
import shutil
import tempfile
import threading
import time
import numpy as np
from cryodrgn.dashboard import chimerax_animation as cx
from cryodrgn.dashboard.data import DashboardExperiment

# Re-export unified ChimeraX animation API (backward compatible).
DEFAULT_GIF_FRAMES = cx.DEFAULT_GIF_FRAMES
DEFAULT_CHIMERAX_PARALLEL = cx.DEFAULT_CHIMERAX_PARALLEL
GIF_DURATION_S = cx.GIF_DURATION_S
EXPLORER_GIF_DURATION_S = cx.EXPLORER_GIF_DURATION_S
explorer_rotation_gif_frames = cx.explorer_rotation_gif_frames
rotation_turn_y_increments = cx.rotation_turn_y_increments
LANDSCAPE_ROTATE_BATCH_MIN_TASKS = cx.LANDSCAPE_ROTATE_BATCH_MIN_TASKS
LANDSCAPE_ROTATE_MAX_VOL_PARALLEL = cx.LANDSCAPE_ROTATE_MAX_VOL_PARALLEL
PARALLEL_VOLUME_JOB_CAP = cx.PARALLEL_VOLUME_JOB_CAP
CHIMERAX_MAX_CONCURRENT_SLOTS = cx.CHIMERAX_MAX_CONCURRENT_SLOTS
ChimeraxViewTurn = cx.ChimeraxViewTurn
ChimeraxPngTask = cx.ChimeraxPngTask
LandscapeRotateSpec = cx.LandscapeRotateSpec
LandscapeStaticView = cx.LandscapeStaticView
chimerax_view_matrix_camera_arg = cx.chimerax_view_matrix_camera_arg
format_chimerax_view_matrix_display = cx.format_chimerax_view_matrix_display
chimerax_volume_color_spec = cx.chimerax_volume_color_spec
normalize_chimerax_view_turns = cx.normalize_chimerax_view_turns
chimerax_path = cx.chimerax_path
run_chimerax_cmds = cx.run_chimerax_cmds
parallel_jobs = cx.parallel_jobs
chimerax_animation_meta = cx.chimerax_animation_meta
mrc_to_static_png = cx.mrc_to_static_png
parallel_chimerax_static_pngs = cx.parallel_chimerax_static_pngs
mrc_to_rotating_gif = cx.mrc_to_rotating_gif
mrc_to_rotating_gif_single_session = cx.mrc_to_rotating_gif_single_session
batch_landscape_rotate_gifs = cx.batch_landscape_rotate_gifs
render_landscape_rotate_gifs = cx.render_landscape_rotate_gifs
render_landscape_cycle_static_views = cx.render_landscape_cycle_static_views
render_static_png = cx.render_static_png
render_static_pngs_parallel = cx.render_static_pngs_parallel
render_rotating_gif = cx.render_rotating_gif
landscape_rotate_render_strategy = cx.landscape_rotate_render_strategy
landscape_rotate_parallel_jobs = cx.landscape_rotate_parallel_jobs
_chimerax_render_cmds = cx._chimerax_render_cmds
_chimerax_rotation_session_cmds = cx._chimerax_rotation_session_cmds
_mpl_retrim_png = cx._mpl_retrim_png

_MONTAGE_SAFE_LETTERS: tuple[str, ...] = tuple(
    chr(c)
    for c in range(ord("A"), ord("Z") + 1)
    if c not in (ord("I"), ord("O"), ord("U"))
)


def montage_cell_label(idx: int) -> str:
    """Montage-style label for linear index ``idx`` (0-based).

    Matches the particle explorer grid.
    """
    idx = int(idx)
    if idx < 0:
        raise ValueError("montage_cell_label idx must be non-negative")
    letters = _MONTAGE_SAFE_LETTERS
    n = len(letters)
    if idx < n:
        return letters[idx]
    j = idx - n
    return letters[j // n] + letters[j % n]


_VOL_MRC_CACHE: dict[str, dict] = {}
_VOL_CACHE_LOCK = threading.Lock()
_VOL_CACHE_MAX_ENTRIES = 32
_VOL_CACHE_TTL_S = 7200.0


def _vol_cache_evict_unlocked(token: str) -> None:
    meta = _VOL_MRC_CACHE.pop(token, None)
    if meta and meta.get("mrc_dir"):
        shutil.rmtree(meta["mrc_dir"], ignore_errors=True)


def _vol_cache_prune_unlocked() -> None:
    now = time.monotonic()
    dead = [
        tok
        for tok, meta in _VOL_MRC_CACHE.items()
        if now - meta["t0"] > _VOL_CACHE_TTL_S
    ]
    for tok in dead:
        _vol_cache_evict_unlocked(tok)
    while len(_VOL_MRC_CACHE) >= _VOL_CACHE_MAX_ENTRIES:
        oldest = min(_VOL_MRC_CACHE.items(), key=lambda kv: kv[1]["t0"])[0]
        _vol_cache_evict_unlocked(oldest)


def _register_vol_mrc_cache(
    mrc_dir: str, vol_files: list[str], rows: tuple[int, ...]
) -> str:
    with _VOL_CACHE_LOCK:
        _vol_cache_prune_unlocked()
        token = secrets.token_urlsafe(24)
        _VOL_MRC_CACHE[token] = {
            "mrc_dir": mrc_dir,
            "vol_files": list(vol_files),
            "rows": rows,
            "t0": time.monotonic(),
        }
        return token


def volume_cell_gif_from_cache(
    token: str,
    cell_index: int,
    *,
    rows_expected: tuple[int, ...],
    gif_frames: int | None = None,
    chimerax_cpus: int = DEFAULT_CHIMERAX_PARALLEL,
) -> bytes:
    """One rotating GIF from a .mrc path kept after montage PNG generation.

    See :func:`generate_montage_volume_pngs`.
    """
    chimerax_cpus = max(1, min(int(chimerax_cpus), 32))
    if gif_frames is None:
        gif_frames = explorer_rotation_gif_frames(chimerax_cpus)
    else:
        gif_frames = max(4, min(int(gif_frames), 120))
    with _VOL_CACHE_LOCK:
        meta = _VOL_MRC_CACHE.get(token)
        if not meta:
            raise ValueError("Unknown or expired volume cache id.")
        if time.monotonic() - meta["t0"] > _VOL_CACHE_TTL_S:
            _vol_cache_evict_unlocked(token)
            raise ValueError("Volume cache expired. Generate volumes again.")
        if meta["rows"] != rows_expected:
            raise ValueError("Montage rows do not match cached volumes.")
        vfs = meta["vol_files"]
        if cell_index < 0 or cell_index >= len(vfs):
            raise ValueError("cell_index out of range for cached volumes.")
        mrc_path = vfs[cell_index]
    with tempfile.TemporaryDirectory(prefix="cryodrgn_explorer_gif_") as gif_dir:
        out_gif = os.path.join(gif_dir, "cell.gif")
        cx.render_rotating_gif(
            mrc_path,
            out_gif,
            gif_frames=gif_frames,
            ncpus=chimerax_cpus,
            gif_duration_s=EXPLORER_GIF_DURATION_S,
        )
        with open(out_gif, "rb") as fh:
            return fh.read()


def save_cached_volumes_to_dir(
    token: str,
    out_dir: str,
    *,
    filename_prefix: str = "volume",
) -> list[str]:
    """Copy cached decoded ``.mrc`` files to a user-selected folder.

    Returns absolute output paths in save order.
    """
    if not out_dir:
        raise ValueError("Choose an output folder.")
    out_dir = os.path.abspath(out_dir)
    with _VOL_CACHE_LOCK:
        meta = _VOL_MRC_CACHE.get(token)
        if not meta:
            raise ValueError("Unknown or expired volume cache id.")
        if time.monotonic() - meta["t0"] > _VOL_CACHE_TTL_S:
            _vol_cache_evict_unlocked(token)
            raise ValueError("Volume cache expired. Generate volumes again.")
        vol_files = list(meta["vol_files"])
    if not vol_files:
        raise ValueError("No cached volumes available to save.")
    os.makedirs(out_dir, exist_ok=True)
    saved_paths: list[str] = []
    for i, src in enumerate(vol_files, start=1):
        if not os.path.isfile(src):
            raise ValueError(
                "Cached volume files are no longer available. Regenerate first."
            )
        dst = os.path.join(out_dir, f"{filename_prefix}_{i:03d}.mrc")
        shutil.copy2(src, dst)
        saved_paths.append(dst)
    return saved_paths


def torch_cuda_available() -> bool:
    """True if ``torch.cuda.is_available()`` (import errors count as False)."""
    try:
        import torch

        return bool(torch.cuda.is_available())
    except Exception:
        return False


def explorer_volumes_eligible(exp: DashboardExperiment) -> bool:
    """GPU present, SPA (not tilt), and weights file exists for this epoch."""
    if not exp.can_preview_particles:
        return False
    if not torch_cuda_available():
        return False
    w = os.path.join(exp.workdir, f"weights.{exp.epoch}.pkl")
    return os.path.isfile(w)


def _config_yaml_path(workdir: str) -> str:
    y = os.path.join(workdir, "config.yaml")
    if os.path.isfile(y):
        return y
    p = os.path.join(workdir, "config.pkl")
    if os.path.isfile(p):
        return p
    raise FileNotFoundError(f"No config.yaml or config.pkl under {workdir}")


def _is_drgnai_config(train_configs: dict) -> bool:
    return "data_norm_mean" in train_configs


def _decode_z_values_classic(
    exp: DashboardExperiment,
    z_values: np.ndarray,
    out_dir: str,
    device: int = 0,
) -> None:
    """Decode arbitrary z-values to ``.mrc`` volumes using classic cryoDRGN."""
    from cryodrgn import analysis

    os.makedirs(out_dir, exist_ok=True)
    zfile = os.path.join(out_dir, "z_values.txt")
    np.savetxt(zfile, z_values)
    weights = os.path.join(exp.workdir, f"weights.{exp.epoch}.pkl")
    cfg = _config_yaml_path(exp.workdir)
    analysis.gen_volumes(
        weights,
        cfg,
        zfile,
        out_dir,
        device=device,
        Apix=1.0,
        vol_start_index=1,
    )


def _drgnai_volume_generator(exp: DashboardExperiment):
    """Build a DRGN-AI ``VolumeGenerator`` from checkpoint + train config."""
    import torch

    from cryodrgn import models_ai as models
    from cryodrgn.analysis_drgnai import VolumeGenerator
    from cryodrgn.lattice import Lattice

    ckpt_path = os.path.join(exp.workdir, f"weights.{exp.epoch}.pkl")
    try:
        checkpoint = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    except TypeError:
        checkpoint = torch.load(ckpt_path, map_location="cpu")
    hypervolume_params = checkpoint["hypervolume_params"]
    hypervolume = models.HyperVolume(**hypervolume_params)
    hypervolume.load_state_dict(checkpoint["hypervolume_state_dict"])
    hypervolume.eval()
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    hypervolume.to(device)

    lattice = Lattice(
        checkpoint["hypervolume_params"]["resolution"],
        extent=0.5,
        device=device,
    )
    zdim = int(checkpoint["hypervolume_params"]["z_dim"])
    radius_mask = checkpoint.get("output_mask_radius")
    tc = exp.train_configs
    data_norm = (
        float(tc.get("data_norm_mean", 0.0)),
        float(tc.get("data_norm_std", 1.0)),
    )
    vg = VolumeGenerator(
        hypervolume,
        lattice,
        zdim,
        invert=False,
        radius_mask=radius_mask,
        data_norm=data_norm,
        vol_start_index=1,
        apix=1.0,
    )
    return vg


def _decode_z_values_drgnai(
    exp: DashboardExperiment,
    z_values: np.ndarray,
    out_dir: str,
) -> None:
    """Decode arbitrary z-values to ``.mrc`` volumes using DRGN-AI."""
    os.makedirs(out_dir, exist_ok=True)
    vg = _drgnai_volume_generator(exp)
    vg.gen_volumes(out_dir, z_values)


def _sorted_vol_mrc_paths(mrc_dir: str, n_take: int) -> list[str]:
    indexed: list[tuple[int, str]] = []
    with os.scandir(mrc_dir) as it:
        for entry in it:
            if not entry.is_file():
                continue
            m = re.search(r"vol_(\d+)", entry.name)
            if m:
                indexed.append((int(m.group(1)), entry.path))
    indexed.sort(key=lambda t: t[0])
    paths = [p for _, p in indexed]
    if len(paths) < n_take:
        raise RuntimeError(
            f"Expected {n_take} volumes, found {len(paths)} under {mrc_dir}."
        )
    return paths[:n_take]


def _decode_z_values_to_vol_paths(
    exp: DashboardExperiment, z_values: np.ndarray, mrc_dir: str
) -> list[str]:
    if _is_drgnai_config(exp.train_configs):
        _decode_z_values_drgnai(exp, z_values, mrc_dir)
    else:
        _decode_z_values_classic(exp, z_values, mrc_dir, device=0)
    return _sorted_vol_mrc_paths(mrc_dir, len(z_values))


def generate_trajectory_volume_pngs(
    exp: DashboardExperiment,
    z_values: np.ndarray,
    *,
    chimerax_cpus: int = DEFAULT_CHIMERAX_PARALLEL,
) -> tuple[list[bytes], str]:
    """Decode volumes along a z-space trajectory and render ChimeraX static PNGs.

    ``z_values`` is an ``(n_points, zdim)`` array of z-latent-space coordinates
    (e.g. from direct interpolation or nearest-neighbor lookup).
    Returns ``(png_bytes_list, cache_token)``.
    """
    z_values = np.asarray(z_values, dtype=np.float64)
    if z_values.ndim != 2 or z_values.shape[1] != exp.z.shape[1]:
        raise ValueError(
            f"z_values must be (n, {exp.z.shape[1]}); got shape {z_values.shape}"
        )

    mrc_dir = tempfile.mkdtemp(prefix="cryodrgn_trajectory_mrc_")
    try:
        vol_files = _decode_z_values_to_vol_paths(exp, z_values, mrc_dir)
        cc = max(1, min(int(chimerax_cpus), 32))
        with tempfile.TemporaryDirectory(prefix="cryodrgn_trajectory_png_") as png_dir:
            tasks = [
                (i, vf, os.path.join(png_dir, f"cell_{i}.png"), 100)
                for i, vf in enumerate(vol_files)
            ]
            paths = parallel_chimerax_static_pngs(tasks, chimerax_cpus=cc)
            png_bytes_list: list[bytes] = []
            for pth in paths:
                with open(pth, "rb") as fh:
                    png_bytes_list.append(fh.read())
        token = _register_vol_mrc_cache(mrc_dir, vol_files, ())
        return png_bytes_list, token
    except Exception:
        shutil.rmtree(mrc_dir, ignore_errors=True)
        raise


def generate_montage_volume_pngs(
    exp: DashboardExperiment,
    rows: list[int],
    *,
    chimerax_cpus: int = DEFAULT_CHIMERAX_PARALLEL,
) -> tuple[list[bytes], str]:
    """Decode volumes for ``rows``, static ChimeraX PNG per cell, and cache .mrc paths.

    Returns ``(png_bytes_list, cache_token)``. The token is required for
    :func:`volume_cell_gif_from_cache` so animations reuse the same decoded volumes.
    """
    rows = [int(r) for r in rows]
    n = len(exp.z)
    if not rows or any(r < 0 or r >= n for r in rows):
        raise ValueError("Invalid plot row indices for volume generation.")

    mrc_dir = tempfile.mkdtemp(prefix="cryodrgn_explorer_mrc_")
    try:
        zsel = exp.z[np.asarray(rows, dtype=int)]
        vol_files = _decode_z_values_to_vol_paths(exp, zsel, mrc_dir)
        cc = max(1, min(int(chimerax_cpus), 32))
        with tempfile.TemporaryDirectory(prefix="cryodrgn_explorer_png_") as png_dir:
            tasks = [
                (i, vf, os.path.join(png_dir, f"cell_{i}.png"), 100)
                for i, vf in enumerate(vol_files)
            ]
            paths = parallel_chimerax_static_pngs(tasks, chimerax_cpus=cc)
            png_bytes_list: list[bytes] = []
            for pth in paths:
                with open(pth, "rb") as fh:
                    png_bytes_list.append(fh.read())
        token = _register_vol_mrc_cache(mrc_dir, vol_files, tuple(rows))
        return png_bytes_list, token
    except Exception:
        shutil.rmtree(mrc_dir, ignore_errors=True)
        raise
