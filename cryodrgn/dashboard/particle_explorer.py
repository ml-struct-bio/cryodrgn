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


_DECODE_PROGRESS_LOCK = threading.Lock()
_VOLUME_JOB_PROGRESS: dict[str, dict[str, object]] = {}
_VOLUME_JOB_PROGRESS_TTL_S = 600.0


def cuda_gpu_count_for_decode() -> int:
    """CUDA devices available for parallel trajectory volume decoding."""
    try:
        import torch

        if torch.cuda.is_available():
            return max(1, int(torch.cuda.device_count()))
    except Exception:
        pass
    return 1


def volume_job_progress_register(
    token: str,
    total: int,
    workers: int,
    phase: str,
    *,
    rerender: bool = False,
) -> None:
    with _DECODE_PROGRESS_LOCK:
        _VOLUME_JOB_PROGRESS[token] = {
            "total": max(0, int(total)),
            "done": 0,
            "workers": max(1, int(workers)),
            "phase": str(phase),
            "rerender": bool(rerender),
            "t0": time.monotonic(),
        }


def volume_job_progress_set_done(token: str, done: int) -> None:
    with _DECODE_PROGRESS_LOCK:
        entry = _VOLUME_JOB_PROGRESS.get(token)
        if not entry:
            return
        total = int(entry["total"])
        entry["done"] = max(0, min(total, int(done)))


def volume_job_progress_snapshot(token: str) -> dict[str, object] | None:
    with _DECODE_PROGRESS_LOCK:
        entry = _VOLUME_JOB_PROGRESS.get(token)
        if not entry:
            return None
        if time.monotonic() - float(entry["t0"]) > _VOLUME_JOB_PROGRESS_TTL_S:
            _VOLUME_JOB_PROGRESS.pop(token, None)
            return None
        total = int(entry["total"])
        done = int(entry["done"])
        pct = (100.0 * done / total) if total else 0.0
        phase = str(entry.get("phase") or "decode")
        workers = int(entry["workers"])
        snap: dict[str, object] = {
            "total": total,
            "done": done,
            "workers": workers,
            "percent": round(pct, 1),
            "phase": phase,
            "rerender": bool(entry.get("rerender")),
        }
        if phase == "decode":
            snap["n_gpus"] = workers
        else:
            snap["n_cpus"] = workers
        return snap


def volume_job_progress_unregister(token: str) -> None:
    with _DECODE_PROGRESS_LOCK:
        _VOLUME_JOB_PROGRESS.pop(token, None)


def decode_progress_register(token: str, total: int, n_gpus: int) -> None:
    volume_job_progress_register(token, total, n_gpus, "decode")


def decode_progress_set_done(token: str, done: int) -> None:
    volume_job_progress_set_done(token, done)


def decode_progress_snapshot(token: str) -> dict[str, object] | None:
    return volume_job_progress_snapshot(token)


def decode_progress_unregister(token: str) -> None:
    volume_job_progress_unregister(token)


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
            decode_progress_set_done(progress_token, _count_vol_mrc_in_dir(mrc_dir))
            if _count_vol_mrc_in_dir(mrc_dir) >= n_total:
                break
            time.sleep(0.15)
        decode_progress_set_done(progress_token, n_total)

    watcher = threading.Thread(target=_watch, daemon=True)
    watcher.start()
    try:
        decode_fn()
    finally:
        stop.set()
        watcher.join(timeout=2.0)
        decode_progress_set_done(progress_token, n_total)


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

    z_values = np.asarray(z_values)
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

    joblib.Parallel(n_jobs=len(tasks))(
        joblib.delayed(_decode_z_values_worker)(
            exp,
            chunk,
            mrc_dir,
            device=device_id,
            vol_start_index=vol_start_index,
        )
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
        decode_progress_register(progress_token, n_total, n_gpus)
        _run_decode_with_mrc_progress(mrc_dir, n_total, progress_token, _decode)
    else:
        _decode()
    return _sorted_vol_mrc_paths(mrc_dir, n_total)


def generate_trajectory_volume_b64_list(
    exp: DashboardExperiment,
    z_values: np.ndarray,
    *,
    progress_token: str | None = None,
) -> tuple[list[dict[str, object]], str]:
    """Decode trajectory z values and return float32 volume blobs for client VTK/slice.

    Returns ``(volume_payloads, cache_token)`` where each payload has
    ``volume_b64``, ``D``, and ``index``.
    """
    from cryodrgn.dashboard.volume_slice_viewer import (
        apply_reconstruction_window,
        vtk_transfer_volume_payload,
    )
    from cryodrgn.mrcfile import parse_mrc

    z_values = np.asarray(z_values, dtype=np.float64)
    if z_values.ndim != 2 or z_values.shape[1] != exp.z.shape[1]:
        raise ValueError(
            f"z_values must be (n, {exp.z.shape[1]}); got shape {z_values.shape}"
        )

    mrc_dir = tempfile.mkdtemp(prefix="cryodrgn_trajectory_mrc_")
    try:
        vol_files = _decode_z_values_to_vol_paths(
            exp, z_values, mrc_dir, progress_token=progress_token
        )
        payloads: list[dict[str, object]] = []
        for i, vf in enumerate(vol_files):
            vol, _ = parse_mrc(vf)
            vol = apply_reconstruction_window(np.asarray(vol, dtype=np.float32), exp)
            transfer = vtk_transfer_volume_payload(vol)
            payloads.append(
                {
                    "index": int(i),
                    **transfer,
                }
            )
        token = _register_vol_mrc_cache(mrc_dir, vol_files, ())
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
    with _VOL_CACHE_LOCK:
        meta = _VOL_MRC_CACHE.get(token)
        if not meta:
            raise ValueError("Unknown or expired volume cache id.")
        if time.monotonic() - meta["t0"] > _VOL_CACHE_TTL_S:
            _vol_cache_evict_unlocked(token)
            raise ValueError("Volume cache expired. Generate volumes again.")
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
    with _VOL_CACHE_LOCK:
        meta = _VOL_MRC_CACHE.get(token)
        if not meta:
            return None
        vol_files = meta.get("vol_files") or []
        return str(vol_files[0]) if vol_files else None


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

    with _VOL_CACHE_LOCK:
        meta = _VOL_MRC_CACHE.get(token)
        if not meta:
            raise ValueError("Unknown or expired volume cache id.")
        if time.monotonic() - meta["t0"] > _VOL_CACHE_TTL_S:
            _vol_cache_evict_unlocked(token)
            raise ValueError("Volume cache expired. Generate volumes again.")
        vol_files = list(meta["vol_files"])

    payloads: list[dict[str, object]] = []
    for i, vf in enumerate(vol_files):
        if not os.path.isfile(vf):
            raise ValueError(
                "Cached volume files are no longer available. Regenerate first."
            )
        vol, _ = parse_mrc(vf)
        vol = apply_reconstruction_window(np.asarray(vol, dtype=np.float32), exp)
        transfer = vtk_transfer_volume_payload(vol)
        payloads.append({"index": int(i), **transfer})
    return payloads


def generate_trajectory_volume_pngs(
    exp: DashboardExperiment,
    z_values: np.ndarray,
    *,
    chimerax_cpus: int = DEFAULT_CHIMERAX_PARALLEL,
    view_matrix_camera: str | None = None,
    view_turns: list[tuple[str, float]] | None = None,
    volume_level: float | None = None,
    progress_token: str | None = None,
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
        vol_files = _decode_z_values_to_vol_paths(
            exp, z_values, mrc_dir, progress_token=progress_token
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
