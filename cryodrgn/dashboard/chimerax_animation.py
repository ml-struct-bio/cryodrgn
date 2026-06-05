"""Unified ChimeraX volume rendering for dashboard animations.

Shared by the particle explorer montage, 3D landscape vol-PCA animations,
and trajectory creator volume strips.

- One ChimeraX session per volume for rotation GIFs (not one process per frame).
- Multi-volume sessions pack several landscape volumes per ChimeraX process.
- Parallel landscape retrim via a thread pool after each session.
- Hybrid parallel across volumes when ``n_vols * gif_frames`` is large enough.
- Parallel worker cap of :data:`PARALLEL_VOLUME_JOB_CAP` (8) on shared nodes.
- Global ChimeraX concurrency limit via semaphore.
"""

from __future__ import annotations

import html
import os
import re
import shlex
import shutil
import subprocess
import tempfile
import threading
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

DEFAULT_GIF_FRAMES = 40
DEFAULT_CHIMERAX_PARALLEL = 16
GIF_DURATION_S = 3.0
# Particle explorer montage cell GIFs: fixed loop length, constant deg/frame.
EXPLORER_GIF_DURATION_S = 4.0


def explorer_rotation_gif_frames(chimerax_cpus: int) -> int:
    """Default explorer rotation frame count: ``2 * chimerax_cpus`` (clamped 4–120)."""
    cpus = max(1, min(int(chimerax_cpus), 32))
    return max(4, min(cpus * 2, 120))


# Parallelize across volumes when total frame work is large enough.
LANDSCAPE_ROTATE_BATCH_MIN_TASKS = 48
# Cap parallel volume workers — ChimeraX oversubscribes above ~8 on shared nodes.
PARALLEL_VOLUME_JOB_CAP = 8
LANDSCAPE_ROTATE_MAX_VOL_PARALLEL = PARALLEL_VOLUME_JOB_CAP
# Global ChimeraX process slots on shared nodes.
CHIMERAX_MAX_CONCURRENT_SLOTS = max(
    1, int(os.environ.get("CRYODRGN_CHIMERAX_MAX_CONCURRENT", "4") or "4")
)
_CHIMERAX_SLOT_SEM = threading.Semaphore(CHIMERAX_MAX_CONCURRENT_SLOTS)

ChimeraxViewTurn = tuple[str, float]


def rotation_turn_y_increments(gif_frames: int) -> list[float]:
    """Per-frame ``turn y`` degrees for one full 360° loop at constant angular speed.

    ChimeraX ``turn`` is relative to the current orientation, so each frame after the
    first uses the same increment (``360 / n`` degrees).
    """
    n = max(4, int(gif_frames))
    step = 360.0 / n
    return [0.0 if i == 0 else step for i in range(n)]


def chimerax_path() -> str:
    """ChimeraX executable from ``CHIMERAX_PATH`` (process environment) or raise.

    When the dashboard was started without ``CHIMERAX_PATH``, the UI can save a path
    in the Flask session; the server mirrors that into ``os.environ`` on each request
    so ChimeraX subprocesses (including joblib workers) see the same executable.
    """
    p = os.environ.get("CHIMERAX_PATH", "").strip()
    if not p:
        raise EnvironmentError(
            "CHIMERAX_PATH is not set. Set it to your ChimeraX executable "
            "(same as for cryodrgn-experiments pipelines), use the dashboard header "
            '"ChimeraX…" control if shown, or restart after exporting CHIMERAX_PATH.'
        )
    return p


def run_chimerax_cmds(
    cmds: list[str], *, catch_errors: bool = False
) -> tuple[str, str]:
    """Run ChimeraX offscreen with semicolon-separated commands.

    See ``pipelines/tile.py``. Limited by :data:`CHIMERAX_MAX_CONCURRENT_SLOTS`.
    """
    cx = chimerax_path()
    joined = " ; ".join(cmds)
    with _CHIMERAX_SLOT_SEM:
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


_CHIMERAX_MATRIX_NUM_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


def chimerax_view_matrix_camera_arg(text: str | None) -> str | None:
    """Parse UI/API matrix text into comma-separated floats.

    For ``view matrix camera``.
    """
    if text is None:
        return None
    raw = str(text).strip()
    if not raw:
        return None
    if raw.lower().startswith("camera"):
        raw = raw[6:].strip()
    nums = _CHIMERAX_MATRIX_NUM_RE.findall(raw)
    if len(nums) < 12:
        raise ValueError(
            "ChimeraX view matrix must contain 12 numbers "
            "(optionally prefixed with 'camera')."
        )
    floats = [float(n) for n in nums[:12]]
    if not all(np.isfinite(f) for f in floats):
        raise ValueError("ChimeraX view matrix numbers must be finite.")
    return ",".join(f"{f:g}" for f in floats)


def format_chimerax_view_matrix_display(camera_arg: str | None) -> str:
    """Format comma-separated camera floats for the dashboard textarea."""
    if not camera_arg or not str(camera_arg).strip():
        return ""
    return "camera " + str(camera_arg).strip()


def chimerax_volume_color_spec(color: str | None) -> str:
    """Normalize Plotly hex or ChimeraX colour names for ``volume color``."""
    s = (color or "").strip()
    if not s:
        return "cornflowerblue"
    if s.startswith("#"):
        return s
    if re.fullmatch(r"[0-9a-fA-F]{6}", s):
        return "#" + s
    return s


def normalize_chimerax_view_turns(
    view_turns: Sequence[tuple[str, float]] | None,
) -> list[ChimeraxViewTurn]:
    """Validated ChimeraX ``turn <axis> <degrees>`` commands for the base view."""
    out: list[ChimeraxViewTurn] = []
    if not view_turns:
        return out
    for axis_raw, degrees_raw in view_turns:
        axis = str(axis_raw).strip().lower()
        if axis not in ("x", "y", "z"):
            raise ValueError(f"Invalid ChimeraX view rotation axis: {axis_raw!r}.")
        degrees = float(degrees_raw)
        if not np.isfinite(degrees):
            raise ValueError("ChimeraX view rotation degrees must be finite.")
        if abs(degrees) > 1e-9:
            out.append((axis, degrees))
    return out


def _extract_chimerax_view_matrix_text(
    stdout: str,
    stderr: str,
    log_path: str | None = None,
) -> str | None:
    """Best-effort extraction of ChimeraX ``view matrix`` output from logs."""
    chunks = [x for x in (stdout, stderr) if x]
    if log_path and os.path.isfile(log_path):
        try:
            with open(log_path, encoding="utf-8", errors="replace") as fh:
                chunks.append(fh.read())
        except OSError:
            pass
    text = "\n".join(chunks)
    if not text.strip():
        return None
    text = html.unescape(re.sub(r"<[^>]+>", "\n", text))
    lines = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        low = line.lower()
        if low.startswith("executing:") or low == "view matrix":
            continue
        lines.append(line)
    num_re = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"

    def _matrix_from_text(label: str, txt: str) -> str | None:
        nums = re.findall(num_re, txt)
        if len(nums) < 12:
            return None
        return label + " " + ",".join(nums[:12])

    camera_lines = [ln for ln in lines if "camera" in ln.lower()]
    for line in camera_lines:
        got = _matrix_from_text("camera", line)
        if got:
            return got
    joined = "\n".join(lines)
    m = re.search(r"camera\b(.{0,800})", joined, flags=re.IGNORECASE | re.DOTALL)
    if m:
        got = _matrix_from_text("camera", m.group(1))
        if got:
            return got
    return None


def _mpl_retrim_png(out_png: str, dpi: int, *, corner_label: str | None = None) -> None:
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
    if corner_label:
        fs = max(11, int(dpi * 0.14))
        ax.text(
            0.02,
            0.98,
            corner_label,
            transform=ax.transAxes,
            fontsize=fs,
            fontweight="bold",
            color="#1a1a1a",
            ha="left",
            va="top",
            bbox=dict(
                boxstyle="round,pad=0.28",
                facecolor="white",
                edgecolor="#374151",
                linewidth=0.8,
                alpha=0.92,
            ),
        )
    fig.savefig(out_png, dpi=dpi, bbox_inches=None, pad_inches=0)
    plt.close(fig)


def chimerax_render_cmds(
    mrc_path: str,
    out_png: str,
    dpi: int,
    *,
    vol_name: str,
    turn_y: float | None,
    volume_color: str | None = None,
    view_turns: Sequence[tuple[str, float]] | None = None,
    view_matrix_camera: str | None = None,
    report_view_matrix: bool = False,
    view_matrix_log_path: str | None = None,
) -> list[str]:
    """ChimeraX cmd list rendering a single view with optional camera turns."""
    qvol = shlex.quote(mrc_path)
    qpng = shlex.quote(out_png)
    vc = chimerax_volume_color_spec(volume_color)
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
    if view_matrix_camera:
        cmds.append(f"view matrix camera {view_matrix_camera} ")
    else:
        for axis, degrees in normalize_chimerax_view_turns(view_turns):
            cmds.append(f"turn {axis} {degrees} ")
    if turn_y is not None:
        cmds.append(f"turn y {turn_y} ")
        cmds.append("volume center #1")
    if report_view_matrix:
        cmds.append("view matrix ")
        if view_matrix_log_path:
            cmds.append(f"log save {shlex.quote(view_matrix_log_path)} ")
    cmds += [f"save {qpng} width {dpi * 5} height {dpi * 5} ", "exit"]
    return cmds


def render_static_png(
    mrc_path: str,
    out_png: str,
    dpi: int = 100,
    *,
    volume_color: str | None = None,
    corner_label: str | None = None,
    view_turns: Sequence[tuple[str, float]] | None = None,
    view_matrix_camera: str | None = None,
    report_view_matrix: bool = False,
) -> str | None:
    """Single ChimeraX view after ``volume center`` and optional base-view turns."""
    matrix_log = os.path.splitext(out_png)[0] + "_view_matrix.html"
    cmds = _chimerax_render_cmds(
        mrc_path,
        out_png,
        dpi,
        vol_name="vol000",
        turn_y=None,
        volume_color=volume_color,
        view_turns=view_turns,
        view_matrix_camera=view_matrix_camera,
        report_view_matrix=report_view_matrix,
        view_matrix_log_path=matrix_log if report_view_matrix else None,
    )
    out, err = run_chimerax_cmds(cmds, catch_errors=False)
    _mpl_retrim_png(out_png, dpi, corner_label=corner_label)
    if report_view_matrix:
        try:
            return _extract_chimerax_view_matrix_text(out, err, matrix_log)
        finally:
            try:
                os.remove(matrix_log)
            except OSError:
                pass
    return None


ChimeraxPngTask = (
    tuple[int, str, str, int]
    | tuple[int, str, str, int, str | None]
    | tuple[int, str, str, int, str | None, str | None]
)


def render_static_pngs_parallel(
    tasks: Sequence[ChimeraxPngTask],
    *,
    chimerax_cpus: int,
) -> list[str]:
    """Run :func:`mrc_to_static_png` for each task in parallel.

    Each task is ``(sort_index, mrc_path, out_png_path, dpi)`` or with optional
    ``volume_color`` and/or ``corner_label`` (hex or ChimeraX color name).
    Returns ``out_png`` paths sorted by ``sort_index``.
    """
    n = len(tasks)
    if n == 0:
        return []
    n_jobs = parallel_jobs(chimerax_cpus, n)

    def _one(task: ChimeraxPngTask) -> tuple[int, str]:
        if len(task) == 6:
            idx, mrc_path, out_png, dpi, vcol, clab = task
            mrc_to_static_png(
                mrc_path,
                out_png,
                dpi=dpi,
                volume_color=vcol,
                corner_label=clab,
            )
        elif len(task) == 5:
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


LandscapeRotateSpec = dict[str, Any]


def _assemble_rotation_pngs_to_gif(
    frame_paths: list[str],
    out_gif: str,
    gif_frames: int,
    *,
    gif_duration_s: float = GIF_DURATION_S,
) -> None:
    from PIL import Image

    images = [Image.open(p) for p in frame_paths]
    duration_s = max(0.5, float(gif_duration_s))
    frame_ms = max(20, int(duration_s * 1000 / gif_frames))
    images[0].save(
        out_gif,
        save_all=True,
        append_images=images[1:],
        duration=frame_ms,
        loop=0,
    )
    for im in images:
        im.close()


def chimerax_rotation_session_cmds(
    mrc_path: str,
    frame_pngs: Sequence[tuple[str, float, str | None]],
    dpi: int,
    *,
    volume_color: str | None = None,
    view_turns: Sequence[ChimeraxViewTurn] | None = None,
    view_matrix_camera: str | None = None,
    report_view_matrix: bool = False,
    vol_name: str = "vol000",
) -> tuple[list[str], str | None]:
    """Build one ChimeraX script: open once, save each rotation, single exit."""
    qvol = shlex.quote(mrc_path)
    vc = chimerax_volume_color_spec(volume_color)
    cmds: list[str] = [
        f"open {qvol} name {vol_name} ",
        "set bgColor white ",
        "volume center #1",
        f"volume color {vc} ",
        "view #1 orient ",
    ]
    matrix_log: str | None = None
    if view_matrix_camera:
        cmds.append(f"view matrix camera {view_matrix_camera} ")
    else:
        for axis, degrees in normalize_chimerax_view_turns(view_turns):
            cmds.append(f"turn {axis} {degrees} ")
    for i, (png, rot, frame_color) in enumerate(frame_pngs):
        if frame_color:
            cmds.append(f"volume color {chimerax_volume_color_spec(frame_color)} ")
        if report_view_matrix and i == 0:
            matrix_log = os.path.splitext(png)[0] + "_view_matrix.html"
            cmds.append("view matrix ")
            if matrix_log:
                cmds.append(f"log save {shlex.quote(matrix_log)} ")
        cmds.append(f"turn y {float(rot)} ")
        cmds.append("volume center #1")
        qpng = shlex.quote(png)
        cmds.append(f"save {qpng} width {dpi * 5} height {dpi * 5} ")
    cmds.append("exit")
    return cmds, matrix_log


@dataclass(frozen=True)
class RotationVolumeInSession:
    """One volume's rotation frames inside a shared ChimeraX session."""

    mrc_path: str
    vol_name: str
    frame_pngs: tuple[tuple[str, float, str | None], ...]
    volume_color: str | None = None
    view_turns: Sequence[ChimeraxViewTurn] | None = None
    view_matrix_camera: str | None = None
    report_view_matrix: bool = False


def chimerax_multi_volume_rotation_session_cmds(
    volumes: Sequence[RotationVolumeInSession],
    dpi: int,
) -> tuple[list[str], str | None]:
    """One ChimeraX script: open each volume in turn, save frames, single exit."""
    if not volumes:
        raise ValueError("chimerax_multi_volume_rotation_session_cmds requires volumes")
    cmds = ["set bgColor white "]
    matrix_log: str | None = None
    for vi, vol in enumerate(volumes):
        qvol = shlex.quote(vol.mrc_path)
        vc = chimerax_volume_color_spec(vol.volume_color)
        cmds.extend(
            [
                f"open {qvol} name {vol.vol_name} ",
                "volume center #1",
                f"volume color {vc} ",
                "view #1 orient ",
            ]
        )
        if vol.view_matrix_camera:
            cmds.append(f"view matrix camera {vol.view_matrix_camera} ")
        else:
            for axis, degrees in normalize_chimerax_view_turns(vol.view_turns):
                cmds.append(f"turn {axis} {degrees} ")
        for fi, (png, rot, frame_color) in enumerate(vol.frame_pngs):
            if frame_color:
                cmds.append(f"volume color {chimerax_volume_color_spec(frame_color)} ")
            if vol.report_view_matrix and fi == 0 and matrix_log is None:
                matrix_log = os.path.splitext(png)[0] + "_view_matrix.html"
                cmds.append("view matrix ")
                cmds.append(f"log save {shlex.quote(matrix_log)} ")
            cmds.append(f"turn y {float(rot)} ")
            cmds.append("volume center #1")
            qpng = shlex.quote(png)
            cmds.append(f"save {qpng} width {dpi * 5} height {dpi * 5} ")
        if vi < len(volumes) - 1:
            cmds.append("close #1 ")
    cmds.append("exit")
    return cmds, matrix_log


_RETRIM_PIPELINE_WORKERS = 8


def _pipeline_worker_count(n_pngs: int) -> int:
    return max(1, min(_RETRIM_PIPELINE_WORKERS, max(1, n_pngs)))


def _retrim_rotation_pngs(
    pngs: Sequence[str],
    dpi: int,
    *,
    corner_label: str | None = None,
    landscape: bool = False,
) -> None:
    """Retrim rotation frames (thread pool for large landscape batches)."""
    if corner_label is not None:
        for png in pngs:
            _mpl_retrim_png(png, dpi, corner_label=corner_label)
        return
    if landscape and len(pngs) > 2:
        from concurrent.futures import ThreadPoolExecutor

        workers = _pipeline_worker_count(len(pngs))
        with ThreadPoolExecutor(max_workers=workers) as ex:
            list(ex.map(lambda p: _mpl_retrim_png(p, dpi), pngs))
        return
    for png in pngs:
        _mpl_retrim_png(png, dpi)


def landscape_rotate_session_groups(
    specs: Sequence[LandscapeRotateSpec],
    chimerax_cpus: int,
) -> list[list[LandscapeRotateSpec]]:
    """Partition rotate specs into groups for multi-volume sessions."""
    specs = list(specs)
    n_vol = len(specs)
    if n_vol <= 1:
        return [specs]
    n_jobs = landscape_rotate_parallel_jobs(chimerax_cpus, n_vol)
    if n_jobs <= 1:
        return [specs]
    per = (n_vol + n_jobs - 1) // n_jobs
    return [specs[i : i + per] for i in range(0, n_vol, per)]


def render_rotating_gif_single_session(
    mrc_path: str,
    out_gif: str,
    *,
    gif_frames: int = DEFAULT_GIF_FRAMES,
    dpi: int = 100,
    volume_color: str | None = None,
    volume_colors: Sequence[str] | None = None,
    corner_label: str | None = None,
    view_turns: Sequence[ChimeraxViewTurn] | None = None,
    view_matrix_camera: str | None = None,
    report_view_matrix: bool = False,
    landscape_rotate: bool = False,
    gif_duration_s: float | None = None,
) -> str | None:
    """One ChimeraX subprocess per volume (all rotation frames in one session)."""

    gif_frames = max(4, int(gif_frames))
    frame_rots = rotation_turn_y_increments(gif_frames)
    loop_duration_s = (
        GIF_DURATION_S if gif_duration_s is None else float(gif_duration_s)
    )
    per_frame_colors: list[str | None] | None = None
    if volume_colors is not None:
        per_frame_colors = [
            str(c).strip() if c is not None and str(c).strip() else None
            for c in volume_colors
        ]
        if len(per_frame_colors) < gif_frames:
            per_frame_colors.extend([None] * (gif_frames - len(per_frame_colors)))
        else:
            per_frame_colors = per_frame_colors[:gif_frames]

    tmpdir = tempfile.mkdtemp(prefix="cryodrgn_rot_session_")
    try:
        frame_pngs: list[tuple[str, float, str | None]] = []
        for fi, rot in enumerate(frame_rots):
            png = os.path.join(tmpdir, f"f{fi:03d}.png")
            fc = None
            if per_frame_colors is not None and fi < len(per_frame_colors):
                fc = per_frame_colors[fi]
            frame_pngs.append((png, float(rot), fc))
        cmds, matrix_log = chimerax_rotation_session_cmds(
            mrc_path,
            frame_pngs,
            dpi,
            volume_color=volume_color,
            view_turns=view_turns,
            view_matrix_camera=view_matrix_camera,
            report_view_matrix=report_view_matrix,
        )
        out, err = run_chimerax_cmds(cmds, catch_errors=False)
        paths = [p for p, _, _ in frame_pngs]
        _retrim_rotation_pngs(
            paths,
            dpi,
            corner_label=corner_label,
            landscape=landscape_rotate,
        )
        view_matrix = (
            _extract_chimerax_view_matrix_text(out, err, matrix_log)
            if report_view_matrix and matrix_log
            else None
        )
        _assemble_rotation_pngs_to_gif(
            paths, out_gif, gif_frames, gif_duration_s=loop_duration_s
        )
        return view_matrix
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def landscape_rotate_render_strategy(
    n_volumes: int, gif_frames: int, chimerax_cpus: int
) -> str:
    """Choose serial or parallel one-session-per-volume rendering."""
    n_volumes = max(0, int(n_volumes))
    gif_frames = max(1, int(gif_frames))
    cpus = max(1, int(chimerax_cpus))
    total = n_volumes * gif_frames
    min_tasks = max(LANDSCAPE_ROTATE_BATCH_MIN_TASKS, cpus * 4)

    if n_volumes <= 1 or total < min_tasks:
        return "serial_session"

    return "parallel_session"


def landscape_rotate_parallel_jobs(chimerax_cpus: int, n_volumes: int) -> int:
    """Worker count for parallel one-session-per-volume rendering."""
    cpus = max(1, min(int(chimerax_cpus), 32))
    n_volumes = max(0, int(n_volumes))
    return max(1, min(cpus, n_volumes, LANDSCAPE_ROTATE_MAX_VOL_PARALLEL))


def _render_landscape_rotate_spec(
    spec: LandscapeRotateSpec,
    *,
    gif_frames: int,
    dpi: int,
) -> tuple[int, str | None]:
    """Render one landscape rotate GIF via a single ChimeraX session."""
    vol_key = int(spec["vol_key"])
    per_frame = spec.get("frame_volume_colors")
    vm = render_rotating_gif_single_session(
        str(spec["mrc_path"]),
        str(spec["out_gif"]),
        gif_frames=gif_frames,
        dpi=dpi,
        volume_color=spec.get("volume_color"),
        volume_colors=per_frame,
        view_turns=spec.get("view_turns"),
        view_matrix_camera=spec.get("view_matrix_camera"),
        report_view_matrix=bool(spec.get("report_view_matrix")),
        landscape_rotate=True,
    )
    return vol_key, vm


def _rotation_frame_plan(
    spec: LandscapeRotateSpec,
    *,
    gif_frames: int,
    tmpdir: str,
) -> tuple[int, str, list[tuple[str, float, str | None]], list[str | None] | None]:
    """PNG paths and per-frame colours for one landscape rotate volume."""
    vol_key = int(spec["vol_key"])
    frame_rots = rotation_turn_y_increments(gif_frames)
    per_frame_colors: list[str | None] | None = None
    volume_colors = spec.get("frame_volume_colors")
    if volume_colors is not None:
        per_frame_colors = [
            str(c).strip() if c is not None and str(c).strip() else None
            for c in volume_colors
        ]
        if len(per_frame_colors) < gif_frames:
            per_frame_colors.extend([None] * (gif_frames - len(per_frame_colors)))
        else:
            per_frame_colors = per_frame_colors[:gif_frames]
    frame_pngs: list[tuple[str, float, str | None]] = []
    for fi, rot in enumerate(frame_rots):
        fc = chimerax_volume_color_spec(spec.get("volume_color"))
        if per_frame_colors is not None:
            raw = per_frame_colors[fi]
            if raw:
                fc = chimerax_volume_color_spec(raw)
        png = os.path.join(tmpdir, f"v{vol_key:03d}_f{fi:03d}.png")
        frame_pngs.append((png, float(rot), fc))
    return vol_key, str(spec["out_gif"]), frame_pngs, per_frame_colors


def render_multi_volume_rotate_session(
    group_specs: Sequence[LandscapeRotateSpec],
    *,
    gif_frames: int,
    dpi: int,
) -> list[tuple[int, str | None]]:
    """Render a group of landscape rotate GIFs in one ChimeraX subprocess."""
    if not group_specs:
        return []
    gif_frames = max(4, int(gif_frames))
    tmpdir = tempfile.mkdtemp(prefix="cryodrgn_landscape_multivol_")
    view_matrix: str | None = None
    try:
        volumes: list[RotationVolumeInSession] = []
        plans: list[tuple[int, str, list[tuple[str, float, str | None]]]] = []
        for si, spec in enumerate(group_specs):
            vol_key, out_gif, frame_pngs, _ = _rotation_frame_plan(
                spec, gif_frames=gif_frames, tmpdir=tmpdir
            )
            plans.append((vol_key, out_gif, frame_pngs))
            volumes.append(
                RotationVolumeInSession(
                    mrc_path=str(spec["mrc_path"]),
                    vol_name=f"vol{si:03d}",
                    frame_pngs=tuple(frame_pngs),
                    volume_color=spec.get("volume_color"),
                    view_turns=spec.get("view_turns"),
                    view_matrix_camera=spec.get("view_matrix_camera"),
                    report_view_matrix=bool(spec.get("report_view_matrix")),
                )
            )
        cmds, matrix_log = chimerax_multi_volume_rotation_session_cmds(volumes, dpi)
        out, err = run_chimerax_cmds(cmds, catch_errors=False)
        if matrix_log:
            view_matrix = _extract_chimerax_view_matrix_text(out, err, matrix_log)
        all_pngs = [png for _, _, fps in plans for png, _, _ in fps]
        _retrim_rotation_pngs(all_pngs, dpi, landscape=True)
        results: list[tuple[int, str | None]] = []
        for i, (vol_key, out_gif, frame_pngs) in enumerate(plans):
            paths = [p for p, _, _ in frame_pngs]
            _assemble_rotation_pngs_to_gif(paths, out_gif, gif_frames)
            vm_out: str | None = None
            if volumes[i].report_view_matrix and view_matrix:
                vm_out = view_matrix
            results.append((vol_key, vm_out))
        return results
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def render_landscape_rotate_gifs(
    vol_specs: Sequence[LandscapeRotateSpec],
    *,
    gif_frames: int,
    chimerax_cpus: int,
    dpi: int = 100,
) -> list[tuple[int, str | None]]:
    """Hybrid scheduler with multi-volume ChimeraX sessions."""
    vol_specs = list(vol_specs)
    gif_frames = max(4, int(gif_frames))
    ncpus = max(1, min(int(chimerax_cpus), 32))
    n_vol = len(vol_specs)
    if n_vol == 0:
        return []
    strategy = landscape_rotate_render_strategy(n_vol, gif_frames, ncpus)
    if strategy in ("serial_session", "parallel_session"):
        specs = list(vol_specs)
        use_multivol = n_vol >= 1

        def _one(s: LandscapeRotateSpec) -> tuple[int, str | None]:
            return _render_landscape_rotate_spec(s, gif_frames=gif_frames, dpi=dpi)

        def _render_group(
            group: list[LandscapeRotateSpec],
        ) -> list[tuple[int, str | None]]:
            return render_multi_volume_rotate_session(
                group, gif_frames=gif_frames, dpi=dpi
            )

        if use_multivol:
            groups = landscape_rotate_session_groups(specs, ncpus)
            if len(groups) == 1:
                pairs = _render_group(groups[0])
            else:
                try:
                    from joblib import Parallel, delayed

                    batch_results = Parallel(n_jobs=len(groups))(
                        delayed(_render_group)(g) for g in groups
                    )
                except ImportError:
                    batch_results = [_render_group(g) for g in groups]
                pairs = [p for batch in batch_results for p in batch]
            pairs.sort(key=lambda x: x[0])
            return pairs

        n_jobs = landscape_rotate_parallel_jobs(ncpus, n_vol)
        if strategy == "serial_session" or n_jobs <= 1:
            return [_one(s) for s in specs]
        try:
            from joblib import Parallel, delayed

            return Parallel(n_jobs=n_jobs)(delayed(_one)(s) for s in specs)
        except ImportError:
            return [_one(s) for s in specs]
    return []


def parallel_jobs(chimerax_cpus: int, n_tasks: int) -> int:
    """Effective parallel workers for volume-level ChimeraX jobs."""
    cpus = max(1, min(int(chimerax_cpus), 32))
    n_tasks = max(0, int(n_tasks))
    if n_tasks <= 0:
        return 1
    return max(1, min(cpus, n_tasks, PARALLEL_VOLUME_JOB_CAP))


def chimerax_animation_meta(chimerax_cpus: int | None = None) -> dict[str, int]:
    """Meta fields for API responses (default / effective CPUs and concurrency cap)."""
    default = max(1, min(int(DEFAULT_CHIMERAX_PARALLEL), 32))
    if chimerax_cpus is None:
        effective = default
    else:
        effective = max(1, min(int(chimerax_cpus), 32))
    return {
        "chimerax_cpus_default": default,
        "chimerax_cpus_effective": effective,
        "chimerax_max_concurrent": CHIMERAX_MAX_CONCURRENT_SLOTS,
        "chimerax_parallel_volume_cap": PARALLEL_VOLUME_JOB_CAP,
    }


@dataclass(frozen=True)
class LandscapeStaticView:
    """One static ChimeraX PNG for landscape cycle mode."""

    mrc_path: str
    out_png: str
    volume_color: str | None = None
    view_turns: Sequence[ChimeraxViewTurn] | None = None
    view_matrix_camera: str | None = None
    report_view_matrix: bool = False


def render_landscape_cycle_static_views(
    views: Sequence[LandscapeStaticView],
    *,
    chimerax_cpus: int,
) -> tuple[list[str], str | None]:
    """Render cycle-mode static PNGs in parallel (one ChimeraX process per volume)."""
    if not views:
        return [], None
    n_jobs = parallel_jobs(chimerax_cpus, len(views))

    def _one(i: int, view: LandscapeStaticView) -> tuple[int, str, str | None]:
        vm = render_static_png(
            view.mrc_path,
            view.out_png,
            dpi=100,
            volume_color=view.volume_color,
            view_turns=view.view_turns,
            view_matrix_camera=view.view_matrix_camera,
            report_view_matrix=view.report_view_matrix,
        )
        return i, view.out_png, vm

    if n_jobs <= 1:
        pairs = [_one(i, v) for i, v in enumerate(views)]
    else:
        try:
            from joblib import Parallel, delayed

            pairs = Parallel(n_jobs=n_jobs)(
                delayed(_one)(i, v) for i, v in enumerate(views)
            )
        except ImportError:
            pairs = [_one(i, v) for i, v in enumerate(views)]
    pairs.sort(key=lambda x: x[0])
    view_matrix = next((m for _, _, m in pairs if m), None)
    return [p for _, p, _ in pairs], view_matrix


def render_rotating_gif(
    mrc_path: str,
    out_gif: str,
    *,
    gif_frames: int = DEFAULT_GIF_FRAMES,
    ncpus: int = DEFAULT_CHIMERAX_PARALLEL,
    dpi: int = 100,
    volume_color: str | None = None,
    volume_colors: Sequence[str] | None = None,
    corner_label: str | None = None,
    view_turns: Sequence[ChimeraxViewTurn] | None = None,
    view_matrix_camera: str | None = None,
    report_view_matrix: bool = False,
    gif_duration_s: float | None = None,
) -> str | None:
    """Rotating GIF via one ChimeraX session per volume (preferred over per-frame processes)."""
    _ = ncpus  # kept for API compatibility; session path does not fan out frames
    return render_rotating_gif_single_session(
        mrc_path,
        out_gif,
        gif_frames=gif_frames,
        dpi=dpi,
        volume_color=volume_color,
        volume_colors=volume_colors,
        corner_label=corner_label,
        view_turns=view_turns,
        view_matrix_camera=view_matrix_camera,
        report_view_matrix=report_view_matrix,
        gif_duration_s=gif_duration_s,
    )


# Backward-compatible aliases used across the dashboard codebase.
mrc_to_static_png = render_static_png
parallel_chimerax_static_pngs = render_static_pngs_parallel
mrc_to_rotating_gif_single_session = render_rotating_gif_single_session
mrc_to_rotating_gif = render_rotating_gif
batch_landscape_rotate_gifs = render_landscape_rotate_gifs
_chimerax_render_cmds = chimerax_render_cmds
_chimerax_rotation_session_cmds = chimerax_rotation_session_cmds
