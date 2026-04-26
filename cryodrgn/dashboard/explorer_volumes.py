"""On-demand decoder volumes + ChimeraX static PNGs and rotating GIFs for the particle explorer."""

from __future__ import annotations

from collections.abc import Sequence
import glob
import os
import re
import secrets
import shlex
import shutil
import subprocess
import tempfile
import threading
import time

import numpy as np

from cryodrgn.dashboard.data import DashboardExperiment

DEFAULT_GIF_FRAMES = 40
DEFAULT_CHIMERAX_PARALLEL = 16
GIF_DURATION_S = 3.0

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
    gif_frames: int = DEFAULT_GIF_FRAMES,
    chimerax_cpus: int = DEFAULT_CHIMERAX_PARALLEL,
) -> bytes:
    """One rotating GIF from a .mrc path kept after :func:`generate_montage_volume_pngs`."""
    gif_frames = max(4, min(int(gif_frames), 120))
    chimerax_cpus = max(1, min(int(chimerax_cpus), 32))
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
        mrc_to_rotating_gif(
            mrc_path,
            out_gif,
            gif_frames=gif_frames,
            ncpus=chimerax_cpus,
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
    try:
        import torch

        return bool(torch.cuda.is_available())
    except Exception:
        return False


def chimerax_path() -> str:
    p = os.environ.get("CHIMERAX_PATH", "").strip()
    if not p:
        raise EnvironmentError(
            "CHIMERAX_PATH is not set. Set it to your ChimeraX executable "
            "(same as for cryodrgn-experiments pipelines) and restart the dashboard."
        )
    return p


def run_chimerax_cmds(
    cmds: list[str], *, catch_errors: bool = False
) -> tuple[str, str]:
    """Run ChimeraX offscreen with semicolon-separated commands (see pipelines/tile.py)."""
    cx = chimerax_path()
    joined = " ; ".join(cmds)
    proc = subprocess.run(
        f"{shlex.quote(cx)} --offscreen --cmd {shlex.quote(joined)}",
        shell=True,
        capture_output=True,
        text=True,
    )
    out, err = proc.stdout or "", proc.stderr or ""
    if catch_errors and err.strip():
        raise RuntimeError(f"ChimeraX stderr:\n{err}\nstdout:\n{out}")
    if proc.returncode != 0:
        raise RuntimeError(f"ChimeraX exited with {proc.returncode}.\n{err}\n{out}")
    return out, err


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
    vol_files = glob.glob(os.path.join(mrc_dir, "vol_*.mrc"))

    def _vol_index(p: str) -> int:
        m = re.search(r"vol_(\d+)", os.path.basename(p))
        return int(m.group(1)) if m else 0

    vol_files.sort(key=_vol_index)
    if len(vol_files) < n_take:
        raise RuntimeError(
            f"Expected {n_take} volumes, found {len(vol_files)} under {mrc_dir}."
        )
    return vol_files[:n_take]


def _mpl_retrim_png(out_png: str, dpi: int) -> None:
    """Re-save ``out_png`` with matplotlib for uniform framing.

    Pads to a centered square and uses a fixed figure bbox so GIF frames stay
    aligned across rotation angles and across volumes of different extents.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    img = np.asarray(plt.imread(out_png))
    if img.ndim == 2:
        img = np.stack([img, img, img], axis=-1)
    elif img.shape[2] == 1:
        img = np.repeat(img, 3, axis=2)

    is_float = np.issubdtype(img.dtype, np.floating)
    h, w, c = int(img.shape[0]), int(img.shape[1]), int(img.shape[2])
    side = max(h, w)
    pad_y = (side - h) // 2
    pad_x = (side - w) // 2

    canvas = np.empty((side, side, c), dtype=img.dtype)
    if is_float:
        canvas[:] = 1.0
        if c == 4:
            canvas[..., 3] = 1.0
    else:
        canvas[:] = 255
        if c == 4:
            canvas[..., 3] = 255
    canvas[pad_y : pad_y + h, pad_x : pad_x + w] = img

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.imshow(canvas, interpolation="nearest", aspect="equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlim(-0.5, side - 0.5)
    ax.set_ylim(side - 0.5, -0.5)
    fig.savefig(out_png, dpi=dpi, bbox_inches=None, pad_inches=0)
    plt.close(fig)


def _chimerax_render_cmds(
    mrc_path: str,
    out_png: str,
    dpi: int,
    *,
    vol_name: str,
    turn_y: float | None,
    volume_color: str | None = None,
) -> list[str]:
    """ChimeraX cmd list rendering a single view (optional Y-axis rotation)."""
    qvol = shlex.quote(mrc_path)
    qpng = shlex.quote(out_png)
    vc = (volume_color or "cornflowerblue").strip()
    if not vc:
        vc = "cornflowerblue"
    cmds = [
        f"open {qvol} name {vol_name} ",
        "set bgColor white ",
        "volume center #1",
        f"volume color {vc} ",
        # Standard orientation + zoom-to-fit for consistent framing.  Do not use
        # ``camera ortho`` here: ``--offscreen`` uses OffScreenRenderingContext,
        # which lacks attributes the ortho camera path expects (e.g. stereo).
        "view #1 orient ",
    ]
    if turn_y is not None:
        cmds.append(f"turn y {turn_y} ")
    cmds += [f"save {qpng} width {dpi * 5} height {dpi * 5} ", "exit"]
    return cmds


def mrc_to_static_png(
    mrc_path: str,
    out_png: str,
    dpi: int = 100,
    *,
    volume_color: str | None = None,
) -> None:
    """Single ChimeraX view (default camera after ``volume center``) — no rotation."""
    cmds = _chimerax_render_cmds(
        mrc_path,
        out_png,
        dpi,
        vol_name="vol000",
        turn_y=None,
        volume_color=volume_color,
    )
    run_chimerax_cmds(cmds, catch_errors=False)
    _mpl_retrim_png(out_png, dpi)


ChimeraxPngTask = tuple[int, str, str, int] | tuple[int, str, str, int, str | None]


def parallel_chimerax_static_pngs(
    tasks: Sequence[ChimeraxPngTask],
    *,
    chimerax_cpus: int,
) -> list[str]:
    """Run :func:`mrc_to_static_png` for each task in parallel.

    Each task is ``(sort_index, mrc_path, out_png_path, dpi)`` or the same with an
    optional fifth element ``volume_color`` (hex or ChimeraX colour name).
    Returns ``out_png`` paths sorted by ``sort_index``.
    """
    n = len(tasks)
    if n == 0:
        return []
    n_jobs = max(1, min(int(chimerax_cpus), 32, n))

    def _one(task: ChimeraxPngTask) -> tuple[int, str]:
        if len(task) == 5:
            idx, mrc_path, out_png, dpi, vcol = task
            mrc_to_static_png(mrc_path, out_png, dpi=dpi, volume_color=vcol)
        else:
            idx, mrc_path, out_png, dpi = task
            mrc_to_static_png(mrc_path, out_png, dpi=dpi)
        return idx, out_png

    if n_jobs <= 1:
        pairs = [_one(t) for t in tasks]
    else:
        try:
            from joblib import Parallel, delayed

            pairs = Parallel(n_jobs=n_jobs)(delayed(_one)(t) for t in tasks)
        except ImportError:
            pairs = [_one(t) for t in tasks]
    pairs.sort(key=lambda x: x[0])
    return [p for _, p in pairs]


def mrc_to_rotating_gif(
    mrc_path: str,
    out_gif: str,
    *,
    gif_frames: int = DEFAULT_GIF_FRAMES,
    ncpus: int = DEFAULT_CHIMERAX_PARALLEL,
    dpi: int = 100,
    volume_color: str | None = None,
) -> None:
    """ChimeraX renders each rotation frame; assemble an animated GIF (cf. pipelines/tile.py)."""
    from PIL import Image

    gif_frames = max(4, int(gif_frames))
    ncpus = max(1, min(int(ncpus), 32))
    frame_rots = np.linspace(0, 360, gif_frames, endpoint=False)
    tmpdir = tempfile.mkdtemp(prefix="cryodrgn_explorer_vol_")
    try:

        def one_frame(frame_rot: float) -> tuple[float, str]:
            png = os.path.join(tmpdir, f"f{frame_rot:.6f}.png")
            vol_name = f"vol000_frame{frame_rot:.4g}"
            cmds = _chimerax_render_cmds(
                mrc_path,
                png,
                dpi,
                vol_name=vol_name,
                turn_y=float(frame_rot),
                volume_color=volume_color,
            )
            run_chimerax_cmds(cmds, catch_errors=False)
            _mpl_retrim_png(png, dpi)
            return frame_rot, png

        try:
            from joblib import Parallel, delayed

            pairs = Parallel(n_jobs=min(ncpus, gif_frames))(
                delayed(one_frame)(rot) for rot in frame_rots
            )
        except ImportError:
            pairs = [one_frame(float(rot)) for rot in frame_rots]
        pairs.sort(key=lambda x: x[0])
        paths = [p for _, p in pairs]
        images = [Image.open(p) for p in paths]
        frame_ms = max(20, int(GIF_DURATION_S * 1000 / gif_frames))
        images[0].save(
            out_gif,
            save_all=True,
            append_images=images[1:],
            duration=frame_ms,
            loop=0,
        )
        for im in images:
            im.close()
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


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
