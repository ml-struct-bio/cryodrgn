"""On-demand decoder volumes + ChimeraX static PNGs and rotating GIFs.

For the particle explorer. ChimeraX rendering lives in
:mod:`cryodrgn.dashboard.chimerax_animation` (shared with landscape vol-PCA and
trajectory creator).
"""

from __future__ import annotations

import os
import re
import shutil
import tempfile
import threading
import time
from concurrent.futures import Future, ThreadPoolExecutor, wait
import numpy as np
from cryodrgn.dashboard import chimerax_animation as cx
from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.volume_caches import VolumeJobStore, VolumeMrcCache

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


_mrc_cache = VolumeMrcCache()
_job_store = VolumeJobStore()

# Test / diagnostic aliases (backward compatible).
_VOL_MRC_CACHE = _mrc_cache.entries
_VOL_CACHE_LOCK = _mrc_cache.lock
_VOL_CACHE_MAX_ENTRIES = _mrc_cache._max_entries


def _vol_cache_evict_unlocked(token: str) -> None:
    _mrc_cache.evict_unlocked(token)


def _vol_cache_prune_unlocked() -> None:
    _mrc_cache.prune_unlocked()


def _register_vol_mrc_cache(
    mrc_dir: str, vol_files: list[str], rows: tuple[int, ...]
) -> str:
    return _mrc_cache.register(mrc_dir, vol_files, rows)


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
    # ``require_meta`` already holds ``_mrc_cache.lock`` (``_VOL_CACHE_LOCK``);
    # do not nest another acquire — ``threading.Lock`` is not re-entrant.
    meta = _mrc_cache.require_meta(token, rows_expected=rows_expected)
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
    meta = _mrc_cache.require_meta(token)
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
    return os.path.isfile(exp.weights_path)


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


_DECODE_PROGRESS_LOCK = _job_store._progress_lock
_VOLUME_JOB_PARTIAL_LOCK = _job_store._partial_lock


def volume_job_progress_register(
    token: str,
    total: int,
    workers: int,
    phase: str,
    *,
    rerender: bool = False,
) -> None:
    _job_store.progress_register(token, total, workers, phase, rerender=rerender)


def volume_job_progress_set_done(token: str, done: int) -> None:
    _job_store.progress_set_done(token, done)


def volume_job_progress_snapshot(token: str) -> dict[str, object] | None:
    return _job_store.progress_snapshot(token)


def volume_job_progress_unregister(token: str) -> None:
    _job_store.progress_unregister(token)


def cuda_gpu_count_for_decode() -> int:
    """CUDA devices available for parallel trajectory volume decoding."""
    try:
        import torch

        if torch.cuda.is_available():
            return max(1, int(torch.cuda.device_count()))
    except Exception:
        pass
    return 1


def trajectory_volume_pipeline_enabled() -> bool:
    """Whether decode and ChimeraX render overlap for trajectory volumes."""
    raw = os.environ.get("CRYODRGN_DASHBOARD_VOLUME_PIPELINE", "1").strip().lower()
    return raw not in ("0", "false", "no", "off")


def volume_job_progress_register_pipeline(
    token: str,
    total: int,
    *,
    n_gpus: int,
    n_cpus: int,
    display_total: int | None = None,
) -> None:
    _job_store.progress_register_pipeline(
        token, total, n_gpus=n_gpus, n_cpus=n_cpus, display_total=display_total
    )


def volume_job_progress_update_pipeline(
    token: str,
    *,
    decode_done: int | None = None,
    render_done: int | None = None,
) -> None:
    _job_store.progress_update_pipeline(
        token, decode_done=decode_done, render_done=render_done
    )


def volume_job_partial_register(token: str, total: int) -> None:
    _job_store.partial_register(token, total)


def volume_job_partial_set_image(token: str, index: int, png_bytes: bytes) -> None:
    _job_store.partial_set_image(token, index, png_bytes)


def volume_job_partial_update_decode(token: str, decode_done: int) -> None:
    _job_store.partial_update_decode(token, decode_done)


def volume_job_partial_mark_complete(
    token: str, view_matrix: str | None = None
) -> None:
    _job_store.partial_mark_complete(token, view_matrix)


def volume_job_partial_snapshot(token: str) -> dict[str, object] | None:
    return _job_store.partial_snapshot(token)


def volume_job_partial_unregister(token: str) -> None:
    _job_store.partial_unregister(token)


def _vol_mrc_index_from_name(name: str) -> int | None:
    m = re.search(r"vol_(\d+)", name)
    if not m:
        return None
    return int(m.group(1)) - 1


def _scan_decoded_mrc_paths(mrc_dir: str) -> dict[int, str]:
    out: dict[int, str] = {}
    if not os.path.isdir(mrc_dir):
        return out
    with os.scandir(mrc_dir) as it:
        for entry in it:
            if not entry.is_file():
                continue
            idx = _vol_mrc_index_from_name(entry.name)
            if idx is not None:
                out[idx] = entry.path
    return out


_MRC_HEADER_BYTES = 1024


def _ready_decoded_mrc_paths(
    mrc_map: dict[int, str],
    stable_sizes: dict[str, int],
) -> tuple[dict[int, str], dict[str, int]]:
    """Return MRC paths whose on-disk size has stabilized (decode finished writing)."""
    ready: dict[int, str] = {}
    next_sizes: dict[str, int] = {}
    for idx, path in mrc_map.items():
        try:
            size = os.path.getsize(path)
        except OSError:
            continue
        next_sizes[path] = size
        if size >= _MRC_HEADER_BYTES and stable_sizes.get(path) == size:
            ready[idx] = path
    return ready, next_sizes


def _count_vol_mrc_in_dir(mrc_dir: str) -> int:
    if not os.path.isdir(mrc_dir):
        return 0
    n = 0
    with os.scandir(mrc_dir) as it:
        for entry in it:
            if entry.is_file() and re.search(r"vol_\d+", entry.name):
                n += 1
    return n


def _run_decode_with_mrc_progress(
    mrc_dir: str,
    n_total: int,
    progress_token: str | None,
    decode_fn,
) -> None:
    if not progress_token:
        decode_fn()
        return

    stop = threading.Event()

    def _watch() -> None:
        while not stop.is_set():
            volume_job_progress_set_done(progress_token, _count_vol_mrc_in_dir(mrc_dir))
            if _count_vol_mrc_in_dir(mrc_dir) >= n_total:
                break
            time.sleep(0.15)
        volume_job_progress_set_done(progress_token, n_total)

    watcher = threading.Thread(target=_watch, daemon=True)
    watcher.start()
    try:
        decode_fn()
    finally:
        stop.set()
        watcher.join(timeout=2.0)
        volume_job_progress_set_done(progress_token, n_total)


def _count_png_in_dir(png_dir: str) -> int:
    if not os.path.isdir(png_dir):
        return 0
    n = 0
    with os.scandir(png_dir) as it:
        for entry in it:
            if entry.is_file() and entry.name.lower().endswith(".png"):
                n += 1
    return n


def _run_chimerax_render_with_png_progress(
    png_dir: str,
    n_total: int,
    progress_token: str | None,
    *,
    n_cpus: int,
    rerender: bool,
    render_fn,
) -> None:
    if not progress_token:
        render_fn()
        return

    volume_job_progress_register(
        progress_token, n_total, n_cpus, "chimerax", rerender=rerender
    )
    stop = threading.Event()

    def _watch() -> None:
        while not stop.is_set():
            volume_job_progress_set_done(progress_token, _count_png_in_dir(png_dir))
            if _count_png_in_dir(png_dir) >= n_total:
                break
            time.sleep(0.15)
        volume_job_progress_set_done(progress_token, n_total)

    watcher = threading.Thread(target=_watch, daemon=True)
    watcher.start()
    try:
        render_fn()
    finally:
        stop.set()
        watcher.join(timeout=2.0)
        volume_job_progress_set_done(progress_token, n_total)


def _coerce_z_decode_matrix(z_values: np.ndarray, zdim: int) -> np.ndarray:
    """Return ``z_values`` as ``(n, zdim)`` for ``np.savetxt`` / ``eval_vol``."""
    z_values = np.asarray(z_values, dtype=np.float64)
    zdim = int(zdim)
    if zdim < 1:
        raise ValueError(f"zdim must be positive; got {zdim}.")
    if z_values.ndim == 1:
        if z_values.size == zdim:
            z_values = z_values.reshape(1, zdim)
        elif z_values.size % zdim != 0:
            raise ValueError(
                f"z_values length {z_values.size} is not a multiple of zdim {zdim}."
            )
        else:
            z_values = z_values.reshape(-1, zdim)
    elif z_values.ndim == 2:
        if z_values.shape[1] != zdim:
            raise ValueError(
                f"z_values must have zdim {zdim}; got shape {z_values.shape}."
            )
    else:
        raise ValueError(f"z_values must be 1D or 2D; got shape {z_values.shape}.")
    return z_values


def _decode_z_values_classic(
    exp: DashboardExperiment,
    z_values: np.ndarray,
    out_dir: str,
    device: int = 0,
    vol_start_index: int = 1,
) -> None:
    """Decode arbitrary z-values to ``.mrc`` volumes using classic cryoDRGN."""
    from cryodrgn import analysis

    os.makedirs(out_dir, exist_ok=True)
    zdim = int(exp.z.shape[1])
    z_values = _coerce_z_decode_matrix(z_values, zdim)
    zfile = os.path.join(out_dir, f"z_values_{int(vol_start_index):04d}.txt")
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
        vol_start_index=vol_start_index,
    )


def _drgnai_volume_generator(exp: DashboardExperiment, device_id: int = 0):
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
    if torch.cuda.is_available():
        device = torch.device(f"cuda:{int(device_id)}")
    else:
        device = torch.device("cpu")
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
    *,
    device: int = 0,
    vol_start_index: int = 1,
) -> None:
    """Decode arbitrary z-values to ``.mrc`` volumes using DRGN-AI."""
    from cryodrgn import models_ai as models
    from cryodrgn.mrcfile import write_mrc

    os.makedirs(out_dir, exist_ok=True)
    zdim = int(exp.z.shape[1])
    z_values = _coerce_z_decode_matrix(z_values, zdim)
    vg = _drgnai_volume_generator(exp, device_id=device)
    vol_start_index = int(vol_start_index)
    for i, z in enumerate(z_values):
        out_mrc = "{}/{}{:03d}.mrc".format(out_dir, "vol_", i + vol_start_index)
        vol = models.eval_volume_method(
            vg.hypervolume,
            vg.lattice,
            vg.zdim,
            vg.data_norm,
            zval=z,
            radius=vg.radius_mask,
        )
        if vg.invert:
            vol *= -1
        write_mrc(out_mrc, vol.cpu().numpy().astype(np.float32), Apix=vg.apix)


def _decode_z_values_worker(
    exp: DashboardExperiment,
    z_values: np.ndarray,
    out_dir: str,
    *,
    device: int = 0,
    vol_start_index: int = 1,
) -> None:
    if _is_drgnai_config(exp.train_configs):
        _decode_z_values_drgnai(
            exp,
            z_values,
            out_dir,
            device=device,
            vol_start_index=vol_start_index,
        )
    else:
        _decode_z_values_classic(
            exp,
            z_values,
            out_dir,
            device=device,
            vol_start_index=vol_start_index,
        )


def _decode_z_values_parallel_impl(
    exp: DashboardExperiment,
    z_values: np.ndarray,
    mrc_dir: str,
    n_gpus: int,
) -> None:
    import joblib

    zdim = int(exp.z.shape[1])
    z_values = _coerce_z_decode_matrix(z_values, zdim)
    chunks = [c for c in np.array_split(z_values, n_gpus) if len(c)]
    if len(chunks) <= 1:
        _decode_z_values_worker(exp, z_values, mrc_dir, device=0, vol_start_index=1)
        return

    os.makedirs(mrc_dir, exist_ok=True)
    vol_start = 1
    tasks: list[tuple[int, np.ndarray, int]] = []
    for device_id, chunk in enumerate(chunks):
        tasks.append((device_id, chunk, vol_start))
        vol_start += len(chunk)

    def _run_parallel_decode(
        device_id: int, chunk: np.ndarray, vol_start_index: int
    ) -> None:
        work_dir = tempfile.mkdtemp(prefix="cryodrgn_decode_", dir=mrc_dir)
        try:
            _decode_z_values_worker(
                exp,
                chunk,
                work_dir,
                device=device_id,
                vol_start_index=vol_start_index,
            )
            for name in os.listdir(work_dir):
                if not name.endswith(".mrc"):
                    continue
                src = os.path.join(work_dir, name)
                dst = os.path.join(mrc_dir, name)
                if os.path.exists(dst):
                    os.remove(dst)
                shutil.move(src, dst)
        finally:
            shutil.rmtree(work_dir, ignore_errors=True)

    joblib.Parallel(n_jobs=len(tasks))(
        joblib.delayed(_run_parallel_decode)(device_id, chunk, vol_start_index)
        for device_id, chunk, vol_start_index in tasks
    )


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
    exp: DashboardExperiment,
    z_values: np.ndarray,
    mrc_dir: str,
    *,
    progress_token: str | None = None,
) -> list[str]:
    z_values = np.asarray(z_values)
    n_total = len(z_values)
    n_gpus = cuda_gpu_count_for_decode()

    def _decode() -> None:
        if n_gpus <= 1:
            _decode_z_values_worker(exp, z_values, mrc_dir, device=0, vol_start_index=1)
        else:
            _decode_z_values_parallel_impl(exp, z_values, mrc_dir, n_gpus)

    if progress_token:
        volume_job_progress_register(progress_token, n_total, n_gpus, "decode")
        _run_decode_with_mrc_progress(mrc_dir, n_total, progress_token, _decode)
    else:
        _decode()
    return _sorted_vol_mrc_paths(mrc_dir, n_total)


def generate_trajectory_volume_b64_list(
    exp: DashboardExperiment,
    z_values: np.ndarray,
    *,
    progress_token: str | None = None,
    trajectory_slot_indices: list[int] | None = None,
    include_transfer_payloads: bool = True,
) -> tuple[list[dict[str, object]], str]:
    """Decode trajectory z values and optionally return VTK/slice transfer blobs.

    Returns ``(volume_payloads, cache_token)``. When ``include_transfer_payloads``
    is true each payload has ``volume_b64``, ``D``, and ``index``. When false
    (Generate / cache-only), payloads are lightweight ``{index, decoded: true}``
    markers so the client can arm Render without shipping multi-hundred-MB
    base64 volumes after decode.
    """
    z_values = np.asarray(z_values, dtype=np.float64)
    if z_values.ndim != 2 or z_values.shape[1] != exp.z.shape[1]:
        raise ValueError(
            f"z_values must be (n, {exp.z.shape[1]}); got shape {z_values.shape}"
        )

    n_total = int(z_values.shape[0])
    if trajectory_slot_indices is not None:
        rows = tuple(int(i) for i in trajectory_slot_indices)
        if len(rows) != n_total:
            raise ValueError("trajectory_slot_indices length must match z_values rows.")
    else:
        rows = tuple(range(n_total))

    mrc_dir = tempfile.mkdtemp(prefix="cryodrgn_trajectory_mrc_")
    try:
        vol_files = _decode_z_values_to_vol_paths(
            exp, z_values, mrc_dir, progress_token=progress_token
        )
        if not include_transfer_payloads:
            token = _register_vol_mrc_cache(mrc_dir, vol_files, rows)
            payloads = [
                {"index": int(rows[i]), "decoded": True} for i in range(len(vol_files))
            ]
            return payloads, token

        from cryodrgn.dashboard.volume_slice_viewer import (
            apply_reconstruction_window,
            vtk_transfer_volume_payload,
        )
        from cryodrgn.mrcfile import parse_mrc

        payloads = []
        for i, vf in enumerate(vol_files):
            vol, _ = parse_mrc(vf)
            vol = apply_reconstruction_window(np.asarray(vol, dtype=np.float32), exp)
            transfer = vtk_transfer_volume_payload(vol)
            payloads.append(
                {
                    "index": int(rows[i]),
                    **transfer,
                }
            )
        token = _register_vol_mrc_cache(mrc_dir, vol_files, rows)
        return payloads, token
    except Exception:
        shutil.rmtree(mrc_dir, ignore_errors=True)
        raise


def _chimerax_png_bytes_from_mrc_paths(
    vol_files: list[str],
    *,
    chimerax_cpus: int = DEFAULT_CHIMERAX_PARALLEL,
    view_matrix_camera: str | None = None,
    view_turns: list[tuple[str, float]] | None = None,
    volume_level: float | None = None,
    progress_token: str | None = None,
    rerender: bool = False,
) -> tuple[list[bytes], str | None]:
    """Render ChimeraX PNG bytes from cached ``.mrc`` paths (parallel when ``cpus > 1``)."""
    if not vol_files:
        return [], None
    from cryodrgn.dashboard.chimerax_animation import resolve_chimerax_volume_level

    level = resolve_chimerax_volume_level(vol_files[0], volume_level)
    cc = max(1, min(int(chimerax_cpus), 32))
    n_total = len(vol_files)

    def _render(png_dir: str) -> tuple[list[str], str | None]:
        views = [
            LandscapeStaticView(
                mrc_path=vf,
                out_png=os.path.join(png_dir, f"cell_{i}.png"),
                volume_level=level,
                view_turns=view_turns,
                view_matrix_camera=view_matrix_camera,
                report_view_matrix=(i == 0),
            )
            for i, vf in enumerate(vol_files)
        ]
        return render_landscape_cycle_static_views(views, chimerax_cpus=cc)

    with tempfile.TemporaryDirectory(prefix="cryodrgn_trajectory_png_") as png_dir:
        paths: list[str] = []
        view_matrix: str | None = None

        def _do_render() -> None:
            nonlocal paths, view_matrix
            paths, view_matrix = _render(png_dir)

        _run_chimerax_render_with_png_progress(
            png_dir,
            n_total,
            progress_token,
            n_cpus=cc,
            rerender=rerender,
            render_fn=_do_render,
        )
        png_bytes_list: list[bytes] = []
        for pth in paths:
            with open(pth, "rb") as fh:
                png_bytes_list.append(fh.read())
    return png_bytes_list, view_matrix


def rerender_chimerax_pngs_from_volume_cache(
    token: str,
    *,
    chimerax_cpus: int = DEFAULT_CHIMERAX_PARALLEL,
    view_matrix_camera: str | None = None,
    view_turns: list[tuple[str, float]] | None = None,
    volume_level: float | None = None,
    progress_token: str | None = None,
) -> tuple[list[bytes], str | None]:
    """Re-render ChimeraX PNGs from a prior trajectory/montage decode cache."""
    meta = _mrc_cache.require_meta(token)
    vol_files = list(meta["vol_files"])
    return _chimerax_png_bytes_from_mrc_paths(
        vol_files,
        chimerax_cpus=chimerax_cpus,
        view_matrix_camera=view_matrix_camera,
        view_turns=view_turns,
        volume_level=volume_level,
        progress_token=progress_token,
        rerender=True,
    )


def primary_mrc_path_from_volume_cache(token: str) -> str | None:
    """First cached ``.mrc`` path for a volume cache token, if any."""
    meta = _mrc_cache.get_meta(token)
    if not meta:
        return None
    vol_files = meta.get("vol_files") or []
    return str(vol_files[0]) if vol_files else None


def volume_cache_slot_indices(token: str) -> tuple[int, ...]:
    """Trajectory slot indices for each cached volume in ``images`` order.

    For trajectory PNG cache tokens created by
    :func:`generate_trajectory_volume_pngs`, the stored ``rows`` metadata
    aligns each cached PNG back onto the corresponding slider index.
    """
    meta = _mrc_cache.require_meta(token)
    rows = meta.get("rows") or ()
    if not rows:
        return ()
    return tuple(int(r) for r in rows)


def trajectory_volume_b64_list_from_cache(
    token: str,
    exp: DashboardExperiment,
) -> list[dict[str, object]]:
    """Build VTK/slice ``volume_b64`` payloads from a prior decode cache."""
    from cryodrgn.dashboard.volume_slice_viewer import (
        apply_reconstruction_window,
        vtk_transfer_volume_payload,
    )
    from cryodrgn.mrcfile import parse_mrc

    meta = _mrc_cache.require_meta(token)
    vol_files = list(meta["vol_files"])
    # Align dense cache order onto trajectory slider slots (partial decode of
    # interiors leaves endpoints out of the MRC cache).
    try:
        slot_indices = volume_cache_slot_indices(token)
    except ValueError:
        slot_indices = ()
    if slot_indices and len(slot_indices) != len(vol_files):
        slot_indices = ()

    payloads: list[dict[str, object]] = []
    for i, vf in enumerate(vol_files):
        if not os.path.isfile(vf):
            raise ValueError(
                "Cached volume files are no longer available. Regenerate first."
            )
        vol, _ = parse_mrc(vf)
        vol = apply_reconstruction_window(np.asarray(vol, dtype=np.float32), exp)
        transfer = vtk_transfer_volume_payload(vol)
        slot_idx = int(slot_indices[i]) if slot_indices else int(i)
        payloads.append({"index": slot_idx, **transfer})
    return payloads


def _pipelined_decode_and_chimerax_pngs(
    exp: DashboardExperiment,
    z_values: np.ndarray,
    mrc_dir: str,
    png_dir: str,
    *,
    chimerax_cpus: int = DEFAULT_CHIMERAX_PARALLEL,
    view_matrix_camera: str | None = None,
    view_turns: list[tuple[str, float]] | None = None,
    volume_level: float | None = None,
    progress_token: str | None = None,
    trajectory_slot_indices: list[int] | None = None,
    n_trajectory_total: int | None = None,
    progress_display_total: int | None = None,
) -> tuple[list[bytes], str | None, list[str]]:
    """Decode trajectory z values while rendering ChimeraX PNGs as each .mrc appears."""
    from cryodrgn.dashboard.chimerax_animation import resolve_chimerax_volume_level

    z_values = np.asarray(z_values, dtype=np.float64)
    n_total = len(z_values)
    slot_indices = (
        [int(i) for i in trajectory_slot_indices]
        if trajectory_slot_indices is not None
        else list(range(n_total))
    )
    if len(slot_indices) != n_total:
        raise ValueError("trajectory_slot_indices length must match z_values rows.")
    n_partial_ui = (
        int(n_trajectory_total) if n_trajectory_total is not None else n_total
    )
    cc = max(1, min(int(chimerax_cpus), 32))
    n_jobs = parallel_jobs(cc, n_total)
    n_gpus = cuda_gpu_count_for_decode()
    os.makedirs(mrc_dir, exist_ok=True)
    os.makedirs(png_dir, exist_ok=True)

    if progress_token:
        volume_job_progress_register_pipeline(
            progress_token,
            n_total,
            n_gpus=n_gpus,
            n_cpus=cc,
            display_total=progress_display_total,
        )
        volume_job_partial_register(progress_token, n_partial_ui)

    decode_exc: list[BaseException] = []
    decode_done = threading.Event()

    def _decode() -> None:
        try:
            if n_gpus <= 1:
                _decode_z_values_worker(
                    exp, z_values, mrc_dir, device=0, vol_start_index=1
                )
            else:
                _decode_z_values_parallel_impl(exp, z_values, mrc_dir, n_gpus)
        except BaseException as err:
            decode_exc.append(err)
        finally:
            decode_done.set()

    decode_thread = threading.Thread(target=_decode, daemon=True)
    decode_thread.start()

    scheduled: set[int] = set()
    rendered: dict[int, bytes] = {}
    view_matrix: str | None = None
    iso_level: float | None = None
    render_lock = threading.Lock()

    def _render_one(idx: int, mrc_path: str) -> tuple[int, bytes, str | None]:
        out_png = os.path.join(png_dir, f"cell_{idx}.png")
        vm = render_static_png(
            mrc_path,
            out_png,
            dpi=100,
            volume_level=iso_level,
            view_turns=view_turns,
            view_matrix_camera=view_matrix_camera,
            report_view_matrix=(idx == 0),
        )
        with open(out_png, "rb") as fh:
            data = fh.read()
        return idx, data, vm

    pending: dict[Future[tuple[int, bytes, str | None]], int] = {}
    mrc_stable_sizes: dict[str, int] = {}

    def _collect_finished(*, block: bool = False) -> None:
        nonlocal view_matrix
        if not pending:
            return

        done, _ = wait(
            tuple(pending),
            timeout=None if block else 0.1,
            return_when="FIRST_COMPLETED",
        )

        for fut in done:
            i, data, vm = fut.result()
            pending.pop(fut, None)
            with render_lock:
                rendered[i] = data
                if vm and view_matrix is None:
                    view_matrix = vm

            if progress_token:
                volume_job_partial_set_image(progress_token, slot_indices[i], data)

    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        while len(rendered) < n_total:
            if decode_exc:
                raise decode_exc[0]

            mrc_map = _scan_decoded_mrc_paths(mrc_dir)
            ready_map, mrc_stable_sizes = _ready_decoded_mrc_paths(
                mrc_map, mrc_stable_sizes
            )
            decode_count = len(ready_map)
            if progress_token:
                volume_job_progress_update_pipeline(
                    progress_token,
                    decode_done=decode_count,
                    render_done=len(rendered),
                )
                volume_job_partial_update_decode(progress_token, decode_count)

            if iso_level is None and volume_level is not None:
                iso_level = resolve_chimerax_volume_level(None, volume_level)

            # ``iso_level is None`` → ChimeraX ``sdLevel 2`` per volume; render as
            # each .mrc becomes ready (no need to wait on a shared absolute level).
            for idx, mrc_path in ready_map.items():
                if idx in scheduled or idx < 0 or idx >= n_total:
                    continue
                scheduled.add(idx)
                fut = executor.submit(_render_one, idx, mrc_path)
                pending[fut] = idx

            _collect_finished()

            if len(rendered) >= n_total:
                break

            if decode_done.is_set() and not pending:
                all_mrc = _scan_decoded_mrc_paths(mrc_dir)
                if len(all_mrc) < n_total:
                    missing_files = [i for i in range(n_total) if i not in all_mrc]
                    if missing_files:
                        raise RuntimeError(
                            "Decode finished but volumes missing at indices: "
                            f"{missing_files}"
                        )

            if not pending:
                time.sleep(0.05)

        decode_thread.join()
        if decode_exc:
            raise decode_exc[0]

        while pending:
            _collect_finished(block=True)
            if progress_token:
                mrc_map = _scan_decoded_mrc_paths(mrc_dir)
                volume_job_progress_update_pipeline(
                    progress_token,
                    decode_done=len(mrc_map),
                    render_done=len(rendered),
                )

    vol_files = _sorted_vol_mrc_paths(mrc_dir, n_total)
    png_bytes_list = [rendered[i] for i in range(n_total)]

    if progress_token:
        volume_job_partial_mark_complete(progress_token, view_matrix)
        volume_job_progress_update_pipeline(
            progress_token,
            decode_done=n_total,
            render_done=n_total,
        )
        volume_job_progress_set_done(progress_token, n_total)

    return png_bytes_list, view_matrix, vol_files


def generate_trajectory_volume_pngs(
    exp: DashboardExperiment,
    z_values: np.ndarray,
    *,
    chimerax_cpus: int = DEFAULT_CHIMERAX_PARALLEL,
    view_matrix_camera: str | None = None,
    view_turns: list[tuple[str, float]] | None = None,
    volume_level: float | None = None,
    progress_token: str | None = None,
    pipeline: bool | None = None,
    trajectory_slot_indices: list[int] | None = None,
    n_trajectory_total: int | None = None,
    progress_display_total: int | None = None,
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
    use_pipeline = (
        trajectory_volume_pipeline_enabled() if pipeline is None else bool(pipeline)
    )
    try:
        if use_pipeline:
            with tempfile.TemporaryDirectory(
                prefix="cryodrgn_trajectory_png_"
            ) as png_dir:
                (
                    png_bytes_list,
                    _view_matrix,
                    vol_files,
                ) = _pipelined_decode_and_chimerax_pngs(
                    exp,
                    z_values,
                    mrc_dir,
                    png_dir,
                    chimerax_cpus=chimerax_cpus,
                    view_matrix_camera=view_matrix_camera,
                    view_turns=view_turns,
                    volume_level=volume_level,
                    progress_token=progress_token,
                    trajectory_slot_indices=trajectory_slot_indices,
                    n_trajectory_total=n_trajectory_total,
                    progress_display_total=progress_display_total,
                )
        else:
            vol_files = _decode_z_values_to_vol_paths(
                exp, z_values, mrc_dir, progress_token=progress_token
            )
            if progress_token:
                cc = max(1, min(int(chimerax_cpus), 32))
                volume_job_progress_register(
                    progress_token,
                    len(vol_files),
                    cc,
                    "chimerax",
                    rerender=False,
                )
            png_bytes_list, _view_matrix = _chimerax_png_bytes_from_mrc_paths(
                vol_files,
                chimerax_cpus=chimerax_cpus,
                view_matrix_camera=view_matrix_camera,
                view_turns=view_turns,
                volume_level=volume_level,
                progress_token=progress_token,
                rerender=False,
            )
        # Store the mapping between the returned PNG order and the trajectory
        # slot indices that each PNG belongs to. This is required so the
        # frontend can re-align cached ChimeraX renders back onto the correct
        # slider indices after partial decoding.
        rows = (
            tuple(int(i) for i in trajectory_slot_indices)
            if trajectory_slot_indices is not None
            else tuple(range(int(z_values.shape[0])))
        )
        token = _register_vol_mrc_cache(mrc_dir, vol_files, rows)
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
            paths = render_static_pngs_parallel(tasks, chimerax_cpus=cc)
            png_bytes_list: list[bytes] = []
            for pth in paths:
                with open(pth, "rb") as fh:
                    png_bytes_list.append(fh.read())
        token = _register_vol_mrc_cache(mrc_dir, vol_files, tuple(rows))
        return png_bytes_list, token
    except Exception:
        shutil.rmtree(mrc_dir, ignore_errors=True)
        raise
