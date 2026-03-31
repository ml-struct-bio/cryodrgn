"""Flask app for the cryoDRGN analysis dashboard."""

from __future__ import annotations

import base64
import io
import logging
import os
import re
import pickle
import shlex
import time
import uuid
from heapq import heappop, heappush
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml
from flask import (
    Flask,
    Response,
    current_app,
    g,
    jsonify,
    redirect,
    render_template,
    request,
    session,
    url_for,
)

from cryodrgn import utils
from cryodrgn.dashboard.data import (
    DashboardExperiment,
    list_z_epochs,
    load_experiment,
    particle_image_array,
)
from cryodrgn.dashboard.command_builder_data import (
    COMMAND_BUILDER_REQUIRED_FIELD_TITLES,
    COMMAND_BUILDER_SCHEMA,
)
from cryodrgn.dashboard.explorer_volumes import (
    DEFAULT_CHIMERAX_PARALLEL,
    DEFAULT_GIF_FRAMES,
    explorer_volumes_eligible,
    generate_montage_volume_pngs,
    generate_trajectory_volume_pngs,
    save_cached_volumes_to_dir,
    volume_cell_gif_from_cache,
)
from cryodrgn.dashboard.mpl_style import ezlab_matplotlib_rc
from cryodrgn.dashboard.plots import (
    normalize_continuous_palette,
    pair_grid_png,
    pair_grid_skeleton_placeholder_layout,
    scatter3d_z_json,
    scatter3d_z_preview_png,
    scatter_json,
)

logger = logging.getLogger(__name__)

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_TEMPLATE_DIR = os.path.join(_THIS_DIR, "templates")
_STATIC_DIR = os.path.join(_THIS_DIR, "static")


def _encode_particle_batch(
    mrcfile: str,
    datadir: str | None,
    global_indices: list[int],
    max_px: int,
) -> list[str]:
    """Load and encode a batch of particle images as base64 JPEG strings.

    Designed to run in a ``ProcessPoolExecutor`` worker — all imports
    are local so the function is self-contained after a fork.
    """
    import base64 as _b64
    import io as _io

    import numpy as _np
    from PIL import Image as PILImage

    from cryodrgn.dataset import ImageDataset

    ds = ImageDataset(mrcfile=mrcfile, lazy=True, datadir=datadir)
    out: list[str] = []
    for gidx in global_indices:
        raw = ds.src.images(gidx, as_numpy=True)
        if raw.ndim == 3:
            raw = raw[0]
        arr = _np.asarray(raw, dtype=_np.float32)
        lo, hi = _np.percentile(arr, (2, 98))
        u8 = (_np.clip((arr - lo) / (hi - lo + 1e-9), 0, 1) * 255).astype(_np.uint8)
        pil = PILImage.fromarray(u8, mode="L")
        if max(pil.size) > max_px:
            pil = pil.resize((max_px, max_px), PILImage.LANCZOS)
        buf = _io.BytesIO()
        pil.save(buf, format="JPEG", quality=85)
        out.append(_b64.b64encode(buf.getvalue()).decode("ascii"))
    return out


# (workdir, epoch, kmeans) -> loaded experiment (invalidated when workdir/epoch changes).
_EXP_CACHE: dict[tuple[str, int, int], DashboardExperiment] = {}
# (epoch, kmeans, xcol, ycol, selection_tuple_or_None) -> (rows, images_b64, elapsed_s).
_PRELOAD_CACHE: dict[
    tuple[int, int, str, str, tuple[int, ...] | None],
    tuple[list[int], list[str], float],
] = {}
_EPOCHS_BY_WORKDIR_CACHE: dict[str, list[int]] = {}
# (workdir, epoch, max_neighbors, avg_neighbors) -> (neighbor_ids, neighbor_dists)
_TRAJ_GRAPH_NEIGHBOR_CACHE: dict[
    tuple[str, int, int, int], tuple[np.ndarray, np.ndarray]
] = {}
_EXP_REQUIRED_ENDPOINTS = {
    "abinit_builder_redirect",
    "filter_page_redirect",
    "api_save_selection",
    "explorer",
    "api_explorer_volume_media",
    "api_scatter",
    "latent_3d_page",
    "api_scatter3d_z",
    "api_latent3d_preview_png",
    "api_preview_montage",
    "api_preload_images",
    "pairplot_page",
    "api_pairplot",
    "api_save_pairplot_png",
    "trajectory_creator_page",
    "api_trajectory_volumes",
    "api_trajectory_coords",
    "api_trajectory_save_volumes",
    "api_trajectory_save_zpath",
    "api_trajectory_import_anchors",
    "api_list_server_files",
    "api_trajectory_kmeans_centers",
    "api_trajectory_random_indices",
    "api_default_trajectory_endpoints",
}


def _config_has_cryodrgn_cmd(config: object) -> bool:
    if not isinstance(config, dict):
        return False
    cmd = config.get("cmd")
    if isinstance(cmd, str):
        return "cryodrgn" in cmd.lower()
    if isinstance(cmd, list):
        return any("cryodrgn" in str(part).lower() for part in cmd)
    return False


def _discover_cryodrgn_workdirs(cwd: str) -> list[str]:
    """Direct subfolders of cwd with config.yaml containing a cryodrgn command."""
    out: list[str] = []
    base = Path(cwd)
    if not base.is_dir():
        return out
    for child in sorted(base.iterdir(), key=lambda p: p.name.lower()):
        if not child.is_dir():
            continue
        cfg = child / "config.yaml"
        if not cfg.is_file():
            continue
        try:
            with cfg.open("r", encoding="utf-8") as fh:
                parsed = yaml.safe_load(fh)
        except Exception:
            continue
        if _config_has_cryodrgn_cmd(parsed):
            out.append(str(child.resolve()))
    return out


def _workdir_options(abs_paths: list[str], base_cwd: str) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for wd in abs_paths:
        try:
            label = os.path.relpath(wd, base_cwd)
        except Exception:
            label = wd
        out.append({"value": wd, "label": label})
    return out


def _filter_ui_scatter_max_points() -> int:
    """Max points when ``filter_ui=1`` on ``/api/scatter`` (env ``CRYODRGN_DASHBOARD_FILTER_MAX_POINTS``)."""
    try:
        v = int(os.environ.get("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", "500000"))
    except ValueError:
        v = 500_000
    return max(50_000, min(v, 2_000_000))


def _default_xy_cols(cols: list[str]) -> tuple[str, str]:
    """Pick sensible default X/Y axes (UMAP if available, else first two)."""
    x = "UMAP1" if "UMAP1" in cols else cols[0]
    y = "UMAP2" if "UMAP2" in cols else cols[min(1, len(cols) - 1)]
    return x, y


def _trajectory_default_xy_cols(cols: list[str], zdim: int) -> tuple[str, str]:
    """Pick default X/Y from the trajectory-allowed axis list only."""
    if zdim > 2 and "PC1" in cols and "PC2" in cols:
        return "PC1", "PC2"
    if "UMAP1" in cols and "UMAP2" in cols:
        return "UMAP1", "UMAP2"
    if len(cols) >= 2:
        return cols[0], cols[1]
    if len(cols) == 1:
        return cols[0], cols[0]
    return cols[0] if cols else "z0", cols[1] if len(cols) > 1 else "z1"


def _trajectory_plot_axis_columns(e: DashboardExperiment) -> list[str]:
    """Axes allowed: z0..z_{D-1}, PC1..PC_D when D>2, UMAP* when UMAP embedding exists."""
    zdim = int(e.z.shape[1])
    seen: set[str] = set()
    ordered: list[str] = []

    def add(name: str) -> None:
        if name in e.plot_df.columns and name not in seen:
            seen.add(name)
            ordered.append(name)

    for i in range(zdim):
        add(f"z{i}")
    if zdim > 2:
        for i in range(1, zdim + 1):
            add(f"PC{i}")
    if _has_umap_columns(e):
        umap_cols = [
            str(c) for c in e.plot_df.columns if re.fullmatch(r"UMAP[0-9]+", str(c))
        ]
        umap_cols.sort(key=lambda c: int(c.replace("UMAP", "")))
        for c in umap_cols:
            add(c)
    return ordered


def _trajectory_axis_column_set(e: DashboardExperiment) -> frozenset[str]:
    return frozenset(_trajectory_plot_axis_columns(e))


def _trajectory_xy_ok_for_direct(xcol: str, ycol: str) -> bool:
    """Direct traversal linearly interpolates z between NN endpoints — only valid in z or PC plot space."""

    def ok_one(c: str) -> bool:
        return bool(re.fullmatch(r"z[0-9]+", c) or re.fullmatch(r"PC[0-9]+", c))

    return ok_one(xcol) and ok_one(ycol)


def _validate_trajectory_plot_axes(
    e: DashboardExperiment, xcol: str, ycol: str
) -> None:
    allowed = _trajectory_axis_column_set(e)
    if xcol not in allowed or ycol not in allowed:
        raise ValueError(
            "That axis combination is not available in trajectory creator."
        )


def _z_traj_to_savetxt_str(z_traj: np.ndarray) -> str:
    """Same text layout as ``cryodrgn graph_traversal -o`` / ``np.savetxt`` on the z rows."""
    buf = io.StringIO()
    np.savetxt(buf, z_traj)
    return buf.getvalue()


def _covariate_display_name(name: str) -> str:
    """Human-friendly covariate names in dashboard selectors."""
    if name == "labels":
        return "k-means labels"
    return name


def _montage_bytes(exp: DashboardExperiment, row_indices: list[int]) -> bytes:
    rows = [int(r) for r in row_indices[:25]]
    with ezlab_matplotlib_rc():
        if not rows:
            fig, ax = plt.subplots(figsize=(4, 4))
            ax.text(
                0.5, 0.5, "Hover points to preview particles", ha="center", va="center"
            )
            ax.axis("off")
            buf = io.BytesIO()
            fig.savefig(buf, format="png", bbox_inches="tight")
            plt.close(fig)
            buf.seek(0)
            return buf.getvalue()

        imgs = [particle_image_array(exp, r) for r in rows]
        n = len(imgs)
        ncol = int(np.ceil(n**0.5))
        nrow = int(np.ceil(n / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(1.8 * ncol, 1.8 * nrow))
        axes_flat = np.atleast_1d(axes).ravel()
        for ax in axes_flat[n:]:
            ax.axis("off")
        for i, img in enumerate(imgs):
            ax = axes_flat[i]
            lo, hi = np.percentile(img, (2, 98))
            disp = np.clip((img - lo) / (hi - lo + 1e-9), 0, 1)
            ax.imshow(disp, cmap="gray")
            ax.set_title(f"idx {int(exp.plot_df.iloc[rows[i]]['index'])}", fontsize=8)
            ax.axis("off")
        plt.tight_layout()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight", dpi=120)
        plt.close(fig)
        buf.seek(0)
        return buf.getvalue()


def _particle_thumbnail_b64_from_row(
    exp: DashboardExperiment, row_index: int, max_side: int = 72
) -> str | None:
    """Small greyscale JPEG (base64) for trajectory NN inset previews."""
    if not exp.can_preview_particles:
        return None
    try:
        from PIL import Image

        img = particle_image_array(exp, int(row_index))
        lo, hi = np.percentile(img, (2, 98))
        disp = np.clip((img - lo) / (hi - lo + 1e-9), 0, 1)
        u8 = (disp * 255).astype(np.uint8)
        pil = Image.fromarray(u8, mode="L")
        w, h = pil.size
        s = max(w, h)
        if s > max_side:
            scale = max_side / float(s)
            pil = pil.resize(
                (max(1, int(w * scale)), max(1, int(h * scale))),
                Image.LANCZOS,
            )
        buf = io.BytesIO()
        pil.save(buf, format="JPEG", quality=88)
        return base64.standard_b64encode(buf.getvalue()).decode("ascii")
    except Exception:
        logger.debug("particle thumbnail encode failed", exc_info=True)
        return None


def _plot_row_particle_index(exp: DashboardExperiment, row_index: int) -> int:
    ri = int(row_index)
    if "index" in exp.plot_df.columns:
        return int(exp.plot_df.iloc[ri]["index"])
    return ri


def _round_direct_mode_traj_xy(traj_xy: np.ndarray) -> np.ndarray:
    """Snap axis coordinates to a short decimal representation (2 or 3 places per axis)."""
    out = np.asarray(traj_xy, dtype=np.float64).copy()
    for j in range(out.shape[1]):
        col = out[:, j]
        finite = col[np.isfinite(col)]
        if finite.size == 0:
            continue
        mx = float(np.nanmax(np.abs(finite)))
        decimals = 2 if mx >= 100.0 else 3
        out[:, j] = np.round(col, decimals)
    return out


def _parse_anchor_indices_txt(raw: bytes) -> list[int]:
    """Parse anchor indices from a text file of whitespace-delimited integers."""
    try:
        txt = raw.decode("utf-8")
    except UnicodeDecodeError as err:
        raise ValueError(
            "Anchor file must be UTF-8 text with space-delimited indices."
        ) from err
    toks = [t for t in re.split(r"[\s,;]+", txt.strip()) if t]
    if len(toks) < 2:
        raise ValueError("Need at least two anchor indices")
    out: list[int] = []
    for t in toks:
        if not re.fullmatch(r"-?\d+", t):
            raise ValueError(f"Invalid anchor index token: {t!r}")
        out.append(int(t))
    return out


def _compute_trajectory_from_anchor_indices(
    e: DashboardExperiment,
    anchor_indices: list[int],
    xcol: str,
    ycol: str,
) -> tuple[np.ndarray, list[int], np.ndarray]:
    """One latent point per anchor (rows in ``z.N.pkl``), in order."""
    n = int(e.z.shape[0])
    rows: list[int] = []
    for a in anchor_indices:
        ai = int(a)
        if ai < 0 or ai >= n:
            raise ValueError(
                f"Anchor index {ai} out of range for z embeddings [0, {n})."
            )
        rows.append(ai)
    zs = np.asarray(rows, dtype=int)
    z_traj = e.z[zs]
    coords = e.plot_df[[xcol, ycol]].values.astype(np.float64)
    traj_xy = coords[zs]
    return z_traj, rows, traj_xy


def _parse_traj_points_value(data: dict, default: int = 4) -> int:
    raw_n = data.get("n_points", default)
    try:
        n_points = max(2, min(int(raw_n), 20))
    except (TypeError, ValueError):
        n_points = default
    return n_points


def _parse_traj_interpolation_value(
    data: dict, key: str = "n_points", default: int = 4
) -> int:
    raw_n = data.get(key, default)
    try:
        n_points = int(raw_n)
    except (TypeError, ValueError):
        n_points = default
    return max(0, min(n_points, 20))


def _parse_traj_neighbor_value(data: dict, key: str, default: int) -> int:
    raw_v = data.get(key, default)
    try:
        v = int(raw_v)
    except (TypeError, ValueError):
        v = default
    return max(2, min(v, 200))


def _compute_direct_anchor_trajectory(
    e: DashboardExperiment,
    anchor_indices: list[int],
    xcol: str,
    ycol: str,
    interpolation_points: int,
) -> tuple[np.ndarray, list[int] | None, np.ndarray]:
    """Match ``cryodrgn direct_traversal`` for anchor interpolation in z-space."""
    z_anchor, rows, anchor_xy = _compute_trajectory_from_anchor_indices(
        e, anchor_indices, xcol, ycol
    )
    if len(rows) < 2:
        raise ValueError("Need at least two anchor indices")
    if interpolation_points < 0:
        raise ValueError("Interpolation points must be >= 0")
    n_points = int(interpolation_points) + 2

    z_parts: list[np.ndarray] = []
    xy_parts: list[np.ndarray] = []
    for i in range(len(rows) - 1):
        z_start = z_anchor[i]
        z_end = z_anchor[i + 1]
        xy_start = anchor_xy[i]
        xy_end = anchor_xy[i + 1]
        t = np.linspace(0.0, 1.0, n_points, dtype=np.float64)
        z_seg = np.outer(1.0 - t, z_start) + np.outer(t, z_end)
        xy_seg = np.outer(1.0 - t, xy_start) + np.outer(t, xy_end)
        z_parts.append(z_seg[:-1])
        xy_parts.append(xy_seg[:-1])

    z_parts.append(z_anchor[-1:].copy())
    xy_parts.append(anchor_xy[-1:].copy())
    return np.concatenate(z_parts, axis=0), None, np.concatenate(xy_parts, axis=0)


def _graph_neighbor_arrays(
    e: DashboardExperiment,
    *,
    max_neighbors: int,
    avg_neighbors: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Nearest-neighbor arrays for graph traversal, cached by workdir+epoch."""
    k = max(2, int(max_neighbors))
    avg = max(1, int(avg_neighbors))
    key = (e.workdir, int(e.epoch), k, avg)
    cached = _TRAJ_GRAPH_NEIGHBOR_CACHE.get(key)
    if cached is not None:
        return cached

    from scipy import spatial

    z = np.asarray(e.z, dtype=np.float64)
    n = int(z.shape[0])
    if n < 2:
        raise ValueError("Need at least two latent points for graph traversal.")
    k_eff = min(k + 1, n)
    tree = spatial.KDTree(z)
    q_dist, q_neighbors = tree.query(z, k=k_eff)
    neighbors = np.atleast_2d(np.asarray(q_neighbors, dtype=np.int64))
    ndist = np.atleast_2d(np.asarray(q_dist, dtype=np.float64))

    # Drop self-neighbor (distance 0); keep up to ``k`` true neighbors.
    neighbors = np.asarray(neighbors[:, 1:], dtype=np.int64)
    ndist = np.asarray(ndist[:, 1:], dtype=np.float64)
    if neighbors.shape[1] == 0:
        raise ValueError("Could not build nearest-neighbor graph from latent points.")

    # Match graph_traversal's average-neighbor thresholding behaviour.
    total_neighbors = max(1, min(int(n * avg), int(ndist.size)))
    flat = np.sort(ndist.reshape(-1))
    max_dist = float(flat[total_neighbors - 1])
    keep = ndist <= max_dist
    neighbors = np.where(keep, neighbors, -1)
    ndist = np.where(keep, ndist, np.inf)

    _TRAJ_GRAPH_NEIGHBOR_CACHE[key] = (neighbors, ndist)
    return neighbors, ndist


def _dijkstra_path_from_neighbors(
    neighbors: np.ndarray, dists: np.ndarray, src: int, dest: int
) -> list[int] | None:
    """Shortest path in a sparse directed neighbor graph."""
    if src == dest:
        return [int(src)]
    n = int(neighbors.shape[0])
    inf = float("inf")
    dist = np.full(n, inf, dtype=np.float64)
    pred = np.full(n, -1, dtype=np.int64)
    visited = np.zeros(n, dtype=bool)
    dist[src] = 0.0
    q: list[tuple[float, int]] = [(0.0, int(src))]

    while q:
        cur_d, v = heappop(q)
        if visited[v]:
            continue
        visited[v] = True
        if v == dest:
            break
        n_ids = neighbors[v]
        n_ds = dists[v]
        for j in range(n_ids.shape[0]):
            u = int(n_ids[j])
            w = float(n_ds[j])
            if u < 0 or not np.isfinite(w):
                continue
            nd = cur_d + w
            if nd < dist[u]:
                dist[u] = nd
                pred[u] = v
                heappush(q, (nd, u))

    if not np.isfinite(dist[dest]):
        return None
    path = [int(dest)]
    cur = int(dest)
    while cur != int(src):
        cur = int(pred[cur])
        if cur < 0:
            return None
        path.append(cur)
    path.reverse()
    return path


def _compute_graph_anchor_trajectory(
    e: DashboardExperiment,
    anchor_indices: list[int],
    xcol: str,
    ycol: str,
    max_neighbors: int,
    avg_neighbors: int,
) -> tuple[np.ndarray, list[int], np.ndarray]:
    """Match ``cryodrgn graph_traversal`` semantics for anchor-to-anchor shortest paths."""
    _z_anchor, rows, _anchor_xy = _compute_trajectory_from_anchor_indices(
        e, anchor_indices, xcol, ycol
    )
    neighbors, ndist = _graph_neighbor_arrays(
        e, max_neighbors=max_neighbors, avg_neighbors=avg_neighbors
    )
    full_path: list[int] = []
    for i in range(len(rows) - 1):
        src = int(rows[i])
        dest = int(rows[i + 1])
        path = _dijkstra_path_from_neighbors(neighbors, ndist, src, dest)
        if not path:
            raise ValueError(
                f"Could not find graph path between anchors {src} and {dest}."
            )
        if full_path and full_path[-1] == path[0]:
            full_path.extend(path[1:])
        else:
            full_path.extend(path)

    if not full_path:
        raise ValueError("Could not construct graph traversal path from given anchors.")

    pidx = np.asarray(full_path, dtype=int)
    z_traj = e.z[pidx]
    coords = e.plot_df[[xcol, ycol]].values.astype(np.float64)
    traj_xy = coords[pidx]
    return z_traj, full_path, traj_xy


def _parse_trajectory_request_body(e: DashboardExperiment, data: dict) -> dict:
    """Validate trajectory POST JSON (shared by coords-only and volume decode)."""
    raw_anchors = data.get("anchor_indices")
    if isinstance(raw_anchors, list) and len(raw_anchors) >= 2:
        try:
            anchor_indices = [int(x) for x in raw_anchors]
        except (TypeError, ValueError) as err:
            raise ValueError("anchor_indices must be a list of integers") from err
        xcol = str(data.get("x", ""))
        ycol = str(data.get("y", ""))
        if xcol not in e.plot_df.columns or ycol not in e.plot_df.columns:
            raise ValueError("bad axis column")
        _validate_trajectory_plot_axes(e, xcol, ycol)
        mode = str(data.get("mode", "direct")).strip().lower()
        if mode not in ("direct", "graph"):
            raise ValueError('anchor mode must be "direct" or "graph".')
        if mode == "direct":
            n_points = _parse_traj_interpolation_value(data, key="n_points", default=0)
        else:
            n_points = _parse_traj_points_value(data, default=4)
        max_neighbors = _parse_traj_neighbor_value(
            data, "max_neighbors", default=max(2, n_points)
        )
        avg_neighbors = _parse_traj_neighbor_value(
            data, "avg_neighbors", default=max(2, n_points)
        )
        return {
            "use_anchors": True,
            "anchor_indices": anchor_indices,
            "xcol": xcol,
            "ycol": ycol,
            "mode": mode,
            "n_points": n_points,
            "max_neighbors": max_neighbors,
            "avg_neighbors": avg_neighbors,
        }

    mode = str(data.get("mode", "direct")).strip().lower()
    if mode not in ("direct", "nearest"):
        raise ValueError('mode must be "direct" or "nearest".')

    raw_traj_xy = data.get("traj_xy")
    traj_xy_custom = None
    sx0 = sy0 = ex0 = ey0 = 0.0
    if isinstance(raw_traj_xy, list) and len(raw_traj_xy) >= 2:
        pts: list[list[float]] = []
        for i, p in enumerate(raw_traj_xy[:200]):
            if not isinstance(p, list) or len(p) != 2:
                raise ValueError(f"traj_xy[{i}] must be [x, y].")
            try:
                px = float(p[0])
                py = float(p[1])
            except (TypeError, ValueError) as err:
                raise ValueError(f"traj_xy[{i}] has non-numeric coordinates.") from err
            pts.append([px, py])
        traj_xy_custom = pts
        sx0, sy0 = pts[0]
        ex0, ey0 = pts[-1]
    else:
        start_xy = data.get("start")
        end_xy = data.get("end")
        if (
            not isinstance(start_xy, list)
            or len(start_xy) != 2
            or not isinstance(end_xy, list)
            or len(end_xy) != 2
        ):
            raise ValueError("start and end must be [x, y] coordinate pairs.")
        try:
            sx0, sy0 = float(start_xy[0]), float(start_xy[1])
            ex0, ey0 = float(end_xy[0]), float(end_xy[1])
        except (TypeError, ValueError) as err:
            raise ValueError("start/end coordinates must be numeric.") from err

    xcol = str(data.get("x", ""))
    ycol = str(data.get("y", ""))
    if xcol not in e.plot_df.columns or ycol not in e.plot_df.columns:
        raise ValueError("bad axis column")
    _validate_trajectory_plot_axes(e, xcol, ycol)
    if mode == "direct" and not _trajectory_xy_ok_for_direct(xcol, ycol):
        raise ValueError(
            "Direct traversal is only available for principal-component or z latent axes."
        )

    n_points = _parse_traj_points_value(data, default=4)

    return {
        "use_anchors": False,
        "mode": mode,
        "sx0": sx0,
        "sy0": sy0,
        "ex0": ex0,
        "ey0": ey0,
        "xcol": xcol,
        "ycol": ycol,
        "n_points": n_points,
        "traj_xy_custom": traj_xy_custom,
    }


def _compute_trajectory_latent_path(
    e: DashboardExperiment, p: dict
) -> tuple[np.ndarray, list[int] | None, np.ndarray]:
    """Return (z_traj, traj_rows or None, traj_xy) for the trajectory line and settings."""
    if p.get("use_anchors"):
        if p["mode"] == "direct":
            return _compute_direct_anchor_trajectory(
                e,
                p["anchor_indices"],
                p["xcol"],
                p["ycol"],
                int(p["n_points"]),
            )
        return _compute_graph_anchor_trajectory(
            e,
            p["anchor_indices"],
            p["xcol"],
            p["ycol"],
            max_neighbors=int(p.get("max_neighbors", p["n_points"])),
            avg_neighbors=int(p.get("avg_neighbors", p["n_points"])),
        )

    mode = p["mode"]
    sx0, sy0 = p["sx0"], p["sy0"]
    ex0, ey0 = p["ex0"], p["ey0"]
    xcol, ycol = p["xcol"], p["ycol"]
    n_points = p["n_points"]
    coords = e.plot_df[[xcol, ycol]].values.astype(np.float64)

    custom_xy = p.get("traj_xy_custom")
    if custom_xy is not None:
        pts = np.asarray(custom_xy, dtype=np.float64)
        traj_rows = []
        for pt in pts:
            dists = np.sum((coords - pt) ** 2, axis=1)
            traj_rows.append(int(np.argmin(dists)))
        z_traj = e.z[np.asarray(traj_rows, dtype=int)]
        if mode == "nearest":
            traj_xy = coords[np.asarray(traj_rows, dtype=int)]
            return z_traj, traj_rows, traj_xy
        traj_xy = _round_direct_mode_traj_xy(pts)
        return z_traj, None, traj_xy

    if mode == "direct":
        start_pt = np.array([sx0, sy0])
        end_pt = np.array([ex0, ey0])
        dists_start = np.sum((coords - start_pt) ** 2, axis=1)
        dists_end = np.sum((coords - end_pt) ** 2, axis=1)
        start_row = int(np.argmin(dists_start))
        end_row = int(np.argmin(dists_end))
        z_start = e.z[start_row]
        z_end = e.z[end_row]
        t = np.linspace(0.0, 1.0, n_points, dtype=np.float64)
        z_traj = np.outer(1.0 - t, z_start) + np.outer(t, z_end)
        traj_rows = None
        traj_xy = np.outer(1.0 - t, start_pt) + np.outer(t, end_pt)
        traj_xy = _round_direct_mode_traj_xy(traj_xy)
    else:
        t = np.linspace(0.0, 1.0, n_points, dtype=np.float64)
        line_xy = np.outer(1.0 - t, np.array([sx0, sy0])) + np.outer(
            t, np.array([ex0, ey0])
        )
        traj_rows = []
        for pt in line_xy:
            dists = np.sum((coords - pt) ** 2, axis=1)
            traj_rows.append(int(np.argmin(dists)))
        z_traj = e.z[np.array(traj_rows)]
        traj_xy = coords[np.array(traj_rows)]

    return z_traj, traj_rows, traj_xy


def _trajectory_shared_json_payload(
    e: DashboardExperiment,
    z_traj: np.ndarray,
    traj_rows: list[int] | None,
    traj_xy: np.ndarray,
    *,
    mode: str,
    n_points: int,
    xcol: str,
    ycol: str,
) -> dict:
    payload: dict = {
        "ok": True,
        "n_points": n_points,
        "mode": mode,
        "traj_rows": traj_rows,
        "traj_xy": traj_xy.tolist(),
        "z_traj": z_traj.tolist(),
        "z_path_txt": _z_traj_to_savetxt_str(z_traj),
        "xcol": xcol,
        "ycol": ycol,
    }
    if traj_rows is not None and mode in ("nearest", "graph"):
        payload["traj_particle_indices"] = [
            _plot_row_particle_index(e, r) for r in traj_rows
        ]
    return payload


def _trajectory_anchor_payload_from_indices(
    e: DashboardExperiment,
    anchor_indices: list[int],
    xcol: str,
    ycol: str,
    *,
    mode: str = "direct",
    n_points: int = 4,
    max_neighbors: int | None = None,
    avg_neighbors: int | None = None,
) -> dict:
    """Build coords/volume JSON for a list of dataset indices (``z.N.pkl`` rows)."""
    mode = str(mode).strip().lower()
    if mode not in ("direct", "graph"):
        mode = "direct"
    if mode == "direct":
        n_points = max(0, min(int(n_points), 20))
    else:
        n_points = max(2, min(int(n_points), 20))
    max_neighbors = (
        n_points if max_neighbors is None else max(2, min(int(max_neighbors), 200))
    )
    avg_neighbors = (
        n_points if avg_neighbors is None else max(2, min(int(avg_neighbors), 200))
    )
    p = {
        "use_anchors": True,
        "anchor_indices": anchor_indices,
        "xcol": xcol,
        "ycol": ycol,
        "mode": mode,
        "n_points": n_points,
        "max_neighbors": max_neighbors,
        "avg_neighbors": avg_neighbors,
    }
    z_traj, traj_rows, traj_xy = _compute_trajectory_latent_path(e, p)
    payload = _trajectory_shared_json_payload(
        e,
        z_traj,
        traj_rows,
        traj_xy,
        mode=mode,
        n_points=n_points,
        xcol=xcol,
        ycol=ycol,
    )
    if "traj_particle_indices" not in payload and mode == "direct":
        n_total = int(np.asarray(z_traj).shape[0])
        pidx = _direct_anchor_particle_indices_payload(
            anchor_indices=anchor_indices,
            interpolation_points=n_points,
            n_total=n_total,
        )
        if pidx is not None:
            payload["traj_particle_indices"] = pidx
    payload["anchor_indices"] = anchor_indices
    return payload


def _direct_anchor_particle_indices_payload(
    *, anchor_indices: list[int], interpolation_points: int, n_total: int
) -> list[int | None] | None:
    """Return row-aligned particle IDs for direct-anchor trajectories.

    Anchor rows carry particle indices; interpolated rows are ``None``.
    """
    if len(anchor_indices) < 2 or n_total <= 0:
        return None
    if interpolation_points <= 0:
        if len(anchor_indices) != n_total:
            return None
        return [int(a) for a in anchor_indices]
    per_seg = int(interpolation_points) + 1
    expected = (len(anchor_indices) - 1) * per_seg + 1
    if expected != n_total:
        return None
    out: list[int | None] = [None] * n_total
    out[0] = int(anchor_indices[0])
    for seg_i in range(1, len(anchor_indices)):
        row_i = seg_i * per_seg
        if 0 <= row_i < n_total:
            out[row_i] = int(anchor_indices[seg_i])
    return out


def _kmeans_centers_ind_path(e: DashboardExperiment) -> str:
    return os.path.join(
        e.workdir,
        f"analyze.{e.epoch}",
        f"kmeans{e.kmeans_folder_id}",
        "centers_ind.txt",
    )


def _load_kmeans_center_indices(e: DashboardExperiment) -> list[int]:
    path = _kmeans_centers_ind_path(e)
    if not os.path.isfile(path):
        raise ValueError(
            f"No k-means centers_ind.txt at {path}. Run cryodrgn analyze with k-means first."
        )
    raw = np.loadtxt(path)
    arr = np.atleast_1d(raw).astype(int).ravel()
    if arr.size < 2:
        raise ValueError("Need at least two k-means center indices")
    return arr.tolist()


def _random_dataset_indices(e: DashboardExperiment, k: int = 10) -> list[int]:
    """Up to ``k`` distinct random row indices into ``z`` (fewer if the dataset is smaller)."""
    n = int(e.z.shape[0])
    if n < 2:
        raise ValueError("Dataset has fewer than 2 particles")
    take = min(int(k), n)
    rng = np.random.default_rng()
    inds = rng.choice(n, size=take, replace=False).astype(int)
    return inds.tolist()


def _plot_df_rows_for_dataset_indices(
    exp: DashboardExperiment, dataset_indices: np.ndarray
) -> list[int]:
    """Map saved dataset indices (values in ``indices.pkl``) to ``plot_df`` row positions."""
    want_arr = np.asarray(dataset_indices).ravel().astype(int, copy=False)
    if want_arr.size == 0:
        return []
    ai = np.asarray(exp.all_indices)
    mask = np.isin(ai, want_arr, assume_unique=False)
    return np.nonzero(mask)[0].astype(int).tolist()


def _load_plot_df_rows_from_plot_inds_file(
    exp: DashboardExperiment, plot_inds_path: str | None
) -> list[int]:
    """Map ``--plot-inds`` pickle dataset indices to ``plot_df`` rows; [] on missing path or error."""
    path = (plot_inds_path or "").strip()
    if not path or not os.path.isfile(path):
        return []
    try:
        pre = utils.load_pkl(path)
        return _plot_df_rows_for_dataset_indices(exp, np.asarray(pre))
    except Exception as err:
        logger.warning("Could not load plot-inds from %s: %s", path, err)
        return []


def _has_umap_columns(exp: DashboardExperiment) -> bool:
    return (
        exp.umap is not None
        and "UMAP1" in exp.plot_df.columns
        and "UMAP2" in exp.plot_df.columns
    )


def _has_pc_columns(exp: DashboardExperiment) -> bool:
    return "PC1" in exp.plot_df.columns and "PC2" in exp.plot_df.columns


def _parse_preselect_rows_param(raw: str | None) -> tuple[list[int] | None, str | None]:
    """Parse ``preselect_rows`` query (comma-separated ints). Returns (rows, error_message)."""
    s = (raw or "").strip()
    if not s:
        return None, None
    try:
        return [int(p) for p in s.split(",") if p.strip()][:5000], None
    except ValueError as exc:
        return None, f"invalid preselect_rows: {exc}"


def _command_builder_template_kwargs(
    exp: DashboardExperiment | None,
) -> dict[str, object]:
    """Template variables for ``command_builder.html`` from experiment config."""
    if exp is None:
        return {
            "default_particles": "",
            "default_ctf": "",
            "default_zdim": 8,
            "default_outdir_abinit": "abinit_run",
            "default_outdir_train": "train_next",
            "default_poses": "",
            "command_builder_schema": COMMAND_BUILDER_SCHEMA,
            "command_builder_required_field_titles": COMMAND_BUILDER_REQUIRED_FIELD_TITLES,
        }

    cfg = exp.train_configs
    da = cfg.get("dataset_args", {}) or {}
    ma = cfg.get("model_args", {}) or {}
    raw_p = da.get("particles") or ""
    default_particles = raw_p if isinstance(raw_p, str) else str(raw_p)
    raw_c = da.get("ctf")
    default_ctf = raw_c if isinstance(raw_c, str) else (str(raw_c) if raw_c else "")
    zd = ma.get("zdim", 8)
    try:
        default_zdim = int(zd)
    except (TypeError, ValueError):
        default_zdim = 8
    raw_poses = da.get("poses")
    if isinstance(raw_poses, str) and raw_poses.strip():
        default_poses = raw_poses.strip()
    else:
        default_poses = os.path.join(exp.workdir, f"pose.{exp.epoch}.pkl")
    return {
        "default_particles": default_particles,
        "default_ctf": default_ctf,
        "default_zdim": default_zdim,
        "default_outdir_abinit": os.path.join(exp.workdir, "abinit_run"),
        "default_outdir_train": os.path.join(exp.workdir, "train_next"),
        "default_poses": default_poses,
        "command_builder_schema": COMMAND_BUILDER_SCHEMA,
        "command_builder_required_field_titles": COMMAND_BUILDER_REQUIRED_FIELD_TITLES,
    }


def _redirect(endpoint: str):
    return redirect(url_for(endpoint), code=302)


def _sync_discovery_session_boot() -> None:
    """Drop stale session workdir when restarting ``cryodrgn dashboard`` with CWD discovery.

    Only applies when the server was started without an outdir and at least one output
    folder was discovered — avoids auto-resuming a previous run from the cookie.
    """
    boot = current_app.config.get("DASHBOARD_DISCOVERY_BOOT_ID")
    if not boot:
        return
    if session.get("dashboard_discovery_boot") == boot:
        return
    session.pop("dashboard_workdir", None)
    session.pop("dashboard_epoch", None)
    session["dashboard_discovery_boot"] = boot


def _active_workdir(app: Flask) -> str | None:
    default_wd = app.config.get("DASHBOARD_WORKDIR")
    candidates = set(app.config.get("DASHBOARD_DISCOVERED_WORKDIRS", []))
    if default_wd:
        candidates.add(default_wd)
    selected = session.get("dashboard_workdir")
    if selected and selected in candidates:
        return str(selected)
    return default_wd


def _epochs_for_workdir(workdir: str) -> list[int]:
    cached = _EPOCHS_BY_WORKDIR_CACHE.get(workdir)
    if cached is not None:
        return cached
    epochs = list_z_epochs(workdir)
    _EPOCHS_BY_WORKDIR_CACHE[workdir] = epochs
    return epochs


def _resolve_epoch(app: Flask) -> int:
    wd = _active_workdir(app)
    if not wd:
        return 0
    epochs = _epochs_for_workdir(wd)
    if not epochs:
        raise RuntimeError("No z.N.pkl epochs in workdir.")
    sess = session.get("dashboard_epoch")
    if sess is None:
        return max(epochs)
    try:
        ep = int(sess)
    except (TypeError, ValueError):
        return max(epochs)
    if ep not in epochs:
        return max(epochs)
    return ep


def _get_dashboard_exp(app: Flask) -> DashboardExperiment:
    wd = _active_workdir(app)
    if not wd:
        raise RuntimeError("No output directory selected.")
    ep = _resolve_epoch(app)
    km = int(app.config["DASHBOARD_KMEANS"])
    key = (wd, ep, km)
    if key not in _EXP_CACHE:
        _EXP_CACHE[key] = load_experiment(
            wd,
            epoch=ep,
            kmeans=km,
        )
    return _EXP_CACHE[key]


def _stratified_xy_row_indices(
    coords: np.ndarray,
    rng: np.random.Generator,
    total_k: int,
) -> set[int]:
    """Local row indices 0..m-1: random + high mean dist to refs + high NN dist (cf. preload)."""
    from scipy.spatial import cKDTree
    from scipy.spatial.distance import cdist as _cdist

    m = int(coords.shape[0])
    if m == 0:
        return set()
    total_k = min(m, total_k)
    k_random = total_k * 2 // 3
    k_mean = total_k // 6
    random_rows = set(
        rng.choice(m, size=min(k_random, m), replace=False).tolist(),
    )

    n_ref = min(m, 500)
    ref_coords = coords[rng.choice(m, size=n_ref, replace=False)]
    avg_dists = np.zeros(m)
    batch = 20_000
    for start in range(0, m, batch):
        avg_dists[start : start + batch] = _cdist(
            coords[start : start + batch],
            ref_coords,
        ).mean(axis=1)

    mean_rows: set[int] = set()
    for idx in np.argsort(-avg_dists):
        if len(mean_rows) >= k_mean:
            break
        r = int(idx)
        if r not in random_rows:
            mean_rows.add(r)

    tree = cKDTree(coords)
    nn_dists = tree.query(coords, k=2)[0][:, 1]

    k_nn = max(0, total_k - k_random - len(mean_rows))
    excluded = random_rows | mean_rows
    nn_rows: set[int] = set()
    for idx in np.argsort(-nn_dists):
        if len(nn_rows) >= k_nn:
            break
        r = int(idx)
        if r not in excluded:
            nn_rows.add(r)

    return random_rows | mean_rows | nn_rows


def _sample_plot_df_rows_for_preload(
    exp: DashboardExperiment,
    xcol: str,
    ycol: str,
    *,
    restrict_to_rows: list[int] | None = None,
) -> tuple[list[int], list[int]]:
    """Pick up to 5k plot_df rows for thumbnail preload (random + spread-out in XY).

    If ``restrict_to_rows`` is set, only those plot_df indices are eligible.
    """
    n = len(exp.plot_df)
    coords_full = exp.plot_df[[xcol, ycol]].values.astype(np.float64)
    rng = np.random.default_rng(42)

    if restrict_to_rows is not None:
        sub_idx = sorted({int(r) for r in restrict_to_rows if 0 <= int(r) < n})
        if not sub_idx:
            return [], []
        coords = coords_full[np.array(sub_idx, dtype=int)]
        inv_map = sub_idx
    else:
        coords = coords_full
        inv_map = list(range(n))

    m = len(coords)
    total_k = min(m, 5000)
    local_sel = _stratified_xy_row_indices(coords, rng, total_k)
    plot_rows_set = {inv_map[i] for i in local_sel}
    pairs = sorted(
        ((r, int(exp.all_indices[r])) for r in plot_rows_set),
        key=lambda p: p[1],
    )
    plot_rows = [p[0] for p in pairs]
    global_indices = [p[1] for p in pairs]
    return plot_rows, global_indices


def _preload_cache_time_estimate_bounds(cpus: int) -> tuple[int, int]:
    """Rough wall-clock range (seconds) for caching up to 5k 96px JPEG thumbs.

    Scales inversely with ``cpus`` (preload worker count), anchored at ~20–30s
    with 4 cores (empirical). I/O and fixed overhead keep a modest floor at high
    core counts.
    """
    cpus = max(1, int(cpus))
    low = max(8, int(round(20 * 4 / cpus)))
    high_parallel = int(round(30 * 4 / cpus))
    high = max(max(18, low + 6), high_parallel)
    return low, high


def _format_preload_cache_time_hint(cpus: int) -> str:
    lo, hi = _preload_cache_time_estimate_bounds(cpus)
    cpus = max(1, int(cpus))
    if lo == hi:
        span = f"~{lo}s"
    else:
        span = f"~{lo}\u2013{hi}s"
    cword = "core" if cpus == 1 else "cores"
    return f"{span} for {cpus} {cword}"


def abbrev_middle(text: object, maxlen: int = 30) -> str:
    """Shorten a string with a middle Unicode ellipsis when longer than ``maxlen``."""
    s = "" if text is None else str(text)
    if len(s) <= maxlen:
        return s
    if maxlen < 4:
        return s[:maxlen]
    ell = "\u2026"
    inner = maxlen - len(ell)
    left = inner // 2
    right = inner - left
    return s[:left] + ell + s[-right:]


def _bind_dashboard_exp() -> None:
    _sync_discovery_session_boot()
    wd = _active_workdir(current_app)
    if not wd:
        if request.endpoint in _EXP_REQUIRED_ENDPOINTS:
            return _redirect("index")
        return
    g.dashboard_exp = _get_dashboard_exp(current_app)


def _cmd_argv_for_nav_display(cmd_parts: list[str]) -> list[str]:
    """Drop the filesystem path to the cryodrgn entrypoint; show ``cryodrgn <args>``."""
    if not cmd_parts:
        return []
    parts = [str(x) for x in cmd_parts]
    if len(parts) >= 3 and parts[1] == "-m":
        mod = parts[2]
        if mod == "cryodrgn" or mod.startswith("cryodrgn."):
            return ["cryodrgn", *parts[3:]]
    base0 = os.path.basename(parts[0])
    if base0 == "cryodrgn" or base0.lower().startswith("cryodrgn."):
        return ["cryodrgn", *parts[1:]]
    if len(parts) >= 2:
        base1 = os.path.basename(parts[1])
        if base1 == "cryodrgn" or base1.lower().startswith("cryodrgn."):
            return ["cryodrgn", *parts[2:]]
    return parts


def _argv_four_command_lines(argv: list[str]) -> list[str]:
    """Format argv as at most four lines.

    Line 1 is argv[0:2] (so it shows e.g. ``cryodrgn abinit``).
    Args longer than 120 characters are abbreviated with an ellipsis-in-the-middle.
    """

    def _abbrev_middle_token(text: str, maxlen: int = 120) -> str:
        s = "" if text is None else str(text)
        if len(s) <= maxlen:
            return s
        if maxlen < 4:
            return s[:maxlen]
        ell = "…"
        inner = maxlen - len(ell)
        left = inner // 2
        right = inner - left
        return s[:left] + ell + s[-right:]

    def _display_join(tokens: list[str]) -> str:
        return " ".join(_abbrev_middle_token(t) for t in tokens)

    if not argv:
        return []
    if len(argv) == 1:
        return [_abbrev_middle_token(argv[0])]
    if len(argv) == 2:
        return [_display_join(argv)]

    head_tokens = argv[0:2]
    head = _display_join(head_tokens)
    rest = argv[2:]

    def _can_break_after(token: str) -> bool:
        """Allow a line break only after argument values or key=value pairs."""
        if "=" in token:
            return True
        # Argument "values" are tokens not starting with "-" / "--".
        return not token.startswith("-")

    def chunk_weight(chunk: list[str]) -> int:
        # Keep line-splitting based on original token lengths (not abbreviated display),
        # so we still get "similar amounts of space".
        if not chunk:
            return 0
        return sum(len(x) for x in chunk) + max(0, len(chunk) - 1)

    # If the tail is tiny, keep it compact (fewer than 4 lines).
    if len(rest) == 1:
        return [head, _abbrev_middle_token(rest[0])]
    if len(rest) == 2:
        # Only split into 3 lines if the first tail token is a valid break point.
        if _can_break_after(rest[0]):
            return [
                head,
                _abbrev_middle_token(rest[0]),
                _abbrev_middle_token(rest[1]),
            ]
        # Otherwise keep both tail tokens on the same line.
        return [head, _display_join(rest)]

    # Split rest into three chunks (lines 2-4) with similar token "weight".
    # Brute force two cut points; `rest` is small (typically CLI args), so O(n^2) is fine.
    n = len(rest)
    avg = chunk_weight(rest) / 3.0
    best_score: float | None = None
    best_i = 1
    best_j = n - 1

    # Cut points i and j are boundaries:
    # - boundary at i breaks after rest[i-1] (so rest[i-1] must be breakable)
    # - boundary at j breaks after rest[j-1] (so rest[j-1] must be breakable)
    for i in range(1, n - 1):
        if not _can_break_after(rest[i - 1]):
            continue
        for j in range(i + 1, n):
            if not _can_break_after(rest[j - 1]):
                continue
            c1 = rest[:i]
            c2 = rest[i:j]
            c3 = rest[j:]
            if not c1 or not c2 or not c3:
                continue
            w1 = chunk_weight(c1)
            w2 = chunk_weight(c2)
            w3 = chunk_weight(c3)
            score = (w1 - avg) ** 2 + (w2 - avg) ** 2 + (w3 - avg) ** 2
            if best_score is None or score < best_score:
                best_score = score
                best_i = i
                best_j = j

    # If no feasible cut points exist (e.g. tail is all flags), fall back to the
    # unconstrained split so we still get a readable display.
    if best_score is None:
        best_score = 0.0
        best_i = 1
        best_j = n - 1

    mid1 = _display_join(rest[:best_i])
    mid2 = _display_join(rest[best_i:best_j])
    last = _display_join(rest[best_j:])
    return [head, mid1, mid2, last]


def inject_meta():
    e: DashboardExperiment = g.dashboard_exp
    discovered_workdirs = current_app.config.get("DASHBOARD_DISCOVERED_WORKDIRS", [])
    discovery_cwd = current_app.config.get("DASHBOARD_DISCOVERY_CWD", os.getcwd())
    epochs = _epochs_for_workdir(e.workdir)
    cfg = e.train_configs
    cmd_list = cfg.get("cmd", [])
    zdim = cfg.get("model_args", {}).get("zdim", "?")
    raw_parts = [str(x) for x in cmd_list] if isinstance(cmd_list, list) else []
    cmd_parts = _cmd_argv_for_nav_display(raw_parts)
    if len(cmd_parts) > 1:
        model_type = cmd_parts[1]
    elif len(cmd_list) > 1:
        model_type = cmd_list[1]
    else:
        model_type = "unknown"
    cfg_train_command = shlex.join(cmd_parts) if cmd_parts else ""
    cfg_cmd_display_lines = (
        _argv_four_command_lines(cmd_parts) if cmd_parts else [f"{model_type} z{zdim}"]
    )
    return {
        "exp_workdir": e.workdir,
        "exp_epoch": e.epoch,
        "exp_kmeans": e.kmeans_folder_id,
        "filter_plot_inds_default": current_app.config["FILTER_PLOT_INDS"] or "",
        "dashboard_epochs": epochs if epochs else [e.epoch],
        "cfg_model_type": model_type,
        "cfg_zdim": zdim,
        "cfg_train_command": cfg_train_command,
        "cfg_cmd_display_lines": cfg_cmd_display_lines,
        "command_builder_only": False,
        "discovered_workdirs": discovered_workdirs,
        "discovered_workdir_options": _workdir_options(
            discovered_workdirs, discovery_cwd
        ),
        "selected_workdir": e.workdir,
    }


def inject_meta_command_builder_only():
    if _active_workdir(current_app) and hasattr(g, "dashboard_exp"):
        return inject_meta()
    discovered_workdirs = current_app.config.get("DASHBOARD_DISCOVERED_WORKDIRS", [])
    discovery_cwd = current_app.config.get("DASHBOARD_DISCOVERY_CWD", os.getcwd())
    active_wd = _active_workdir(current_app)
    if active_wd:
        epochs = _epochs_for_workdir(active_wd)
        exp_epoch = _resolve_epoch(current_app) if epochs else 0
    else:
        epochs = []
        exp_epoch = 0
    return {
        "exp_workdir": "",
        "exp_epoch": exp_epoch,
        "exp_kmeans": -1,
        "filter_plot_inds_default": "",
        "dashboard_epochs": epochs if epochs else [0],
        "cfg_model_type": "cryodrgn",
        "cfg_zdim": "",
        "cfg_train_command": "",
        "cfg_cmd_display_lines": ["No experiment loaded"],
        "command_builder_only": not bool(active_wd),
        "discovered_workdirs": discovered_workdirs,
        "discovered_workdir_options": _workdir_options(
            discovered_workdirs, discovery_cwd
        ),
        "selected_workdir": active_wd or "",
    }


def api_set_epoch():
    wd = _active_workdir(current_app)
    if not wd:
        return jsonify(error="Select an output folder first."), 400
    data = request.get_json(force=True, silent=True) or {}
    raw_epoch = data.get("epoch")
    if raw_epoch is None:
        return jsonify(error="Invalid epoch."), 400
    try:
        ep = int(raw_epoch)
    except (TypeError, ValueError):
        return jsonify(error="Invalid epoch."), 400
    if ep not in _epochs_for_workdir(wd):
        return jsonify(error="Epoch not available for this output folder."), 400
    session["dashboard_epoch"] = ep
    _EXP_CACHE.clear()
    _PRELOAD_CACHE.clear()
    _TRAJ_GRAPH_NEIGHBOR_CACHE.clear()
    return jsonify(ok=True, epoch=ep)


def api_set_workdir():
    data = request.get_json(force=True, silent=True) or {}
    raw = data.get("workdir")
    candidates = set(current_app.config.get("DASHBOARD_DISCOVERED_WORKDIRS", []))

    if raw is None or (isinstance(raw, str) and not raw.strip()):
        if not current_app.config.get("COMMAND_BUILDER_ONLY"):
            return jsonify(error="Cannot clear output folder in this mode."), 400
        if not candidates:
            return jsonify(error="No output folders available."), 400
        session.pop("dashboard_workdir", None)
        session.pop("dashboard_epoch", None)
        boot = current_app.config.get("DASHBOARD_DISCOVERY_BOOT_ID")
        if boot:
            session["dashboard_discovery_boot"] = boot
        _EXP_CACHE.clear()
        _PRELOAD_CACHE.clear()
        _TRAJ_GRAPH_NEIGHBOR_CACHE.clear()
        return jsonify(ok=True, workdir=None)

    requested = str(raw).strip()
    if requested not in candidates:
        return jsonify(error="Invalid output folder."), 400
    epochs = _epochs_for_workdir(requested)
    if not epochs:
        return jsonify(error="No analyzed epochs found in selected output folder."), 400
    session["dashboard_workdir"] = requested
    session["dashboard_epoch"] = max(epochs)
    boot = current_app.config.get("DASHBOARD_DISCOVERY_BOOT_ID")
    if boot:
        session["dashboard_discovery_boot"] = boot
    _EXP_CACHE.clear()
    _PRELOAD_CACHE.clear()
    _TRAJ_GRAPH_NEIGHBOR_CACHE.clear()
    return jsonify(ok=True, workdir=requested, epoch=max(epochs))


def index():
    if not _active_workdir(current_app):
        return render_template(
            "index.html",
            can_images=False,
            zdim=0,
            show_trajectory_creator=False,
            command_builder_only=True,
        )
    e: DashboardExperiment = g.dashboard_exp
    return render_template(
        "index.html",
        can_images=e.can_preview_particles,
        zdim=int(e.z.shape[1]),
        show_trajectory_creator=explorer_volumes_eligible(e),
        command_builder_only=False,
    )


def command_builder_page():
    e: DashboardExperiment | None = None
    if _active_workdir(current_app):
        e = g.dashboard_exp
    return render_template(
        "command_builder.html",
        **_command_builder_template_kwargs(e),
    )


def abinit_builder_redirect():
    return _redirect("command_builder_page")


def filter_page_redirect():
    return _redirect("explorer")


def api_save_selection():
    e: DashboardExperiment = g.dashboard_exp
    data = request.get_json(force=True, silent=True) or {}
    rows_raw = data.get("rows")
    if not isinstance(rows_raw, list) or len(rows_raw) == 0:
        return jsonify(error="No particles selected."), 400
    try:
        rows = sorted({int(r) for r in rows_raw})
    except (TypeError, ValueError):
        return jsonify(error="Invalid row list."), 400
    n_df = len(e.plot_df)
    if any(r < 0 or r >= n_df for r in rows):
        return jsonify(error="Row index out of range."), 400
    selected_ds = np.asarray(e.all_indices[rows], dtype=int)
    force = bool(data.get("force"))
    save_inverse = bool(data.get("save_inverse"))
    basename = (data.get("basename") or "indices").strip() or "indices"
    basename = os.path.basename(basename)
    if force:
        out_base = os.path.join(e.workdir, "indices")
    else:
        sel_dir = (data.get("sel_dir") or "").strip()
        dir_abs = os.path.abspath(sel_dir) if sel_dir else os.path.abspath(e.workdir)
        if not os.path.isdir(dir_abs):
            return jsonify(error=f"Directory does not exist: {dir_abs}"), 400
        out_base = os.path.join(dir_abs, basename)
    path_main = out_base + ".pkl"
    path_inv = out_base + "_inverse.pkl"
    parent = os.path.dirname(path_main)
    if parent:
        os.makedirs(parent, exist_ok=True)
    inverse_ds = np.setdiff1d(np.asarray(e.all_indices, dtype=int), selected_ds)
    try:
        with open(path_main, "wb") as fh:
            pickle.dump(selected_ds, fh)
        if save_inverse:
            with open(path_inv, "wb") as fh:
                pickle.dump(inverse_ds, fh)
    except OSError as err:
        logger.exception("save selection failed")
        return jsonify(error=str(err)), 500
    payload = {
        "ok": True,
        "path": path_main,
        "n_selected": int(selected_ds.size),
    }
    if save_inverse:
        payload["inverse_path"] = path_inv
    return jsonify(payload)


def explorer():
    e: DashboardExperiment = g.dashboard_exp
    if not e.can_preview_particles:
        return (
            render_template(
                "no_images.html",
                reason="Tilt-series data does not include single-particle thumbnails in this view.",
            ),
            200,
        )
    cols = e.numeric_columns
    dx, dy = _default_xy_cols(cols)
    initial_rows = _load_plot_df_rows_from_plot_inds_file(
        e,
        current_app.config.get("FILTER_PLOT_INDS"),
    )
    pc = int(current_app.config["PRELOAD_CPUS"])
    return render_template(
        "scatter_explorer.html",
        numeric_cols=cols,
        covariate_display_map={c: _covariate_display_name(c) for c in cols},
        default_x=dx,
        default_y=dy,
        initial_rows=initial_rows,
        total_particles=int(len(e.all_indices)),
        workdir=e.workdir,
        preload_cache_time_hint=_format_preload_cache_time_hint(pc),
        show_volume_explorer=explorer_volumes_eligible(e),
    )


def api_explorer_volume_media():
    e: DashboardExperiment = g.dashboard_exp
    if not explorer_volumes_eligible(e):
        return (
            jsonify(
                error="Volume explorer needs a CUDA GPU, single-particle data, "
                "and weights for the current epoch.",
            ),
            400,
        )
    data = request.get_json(force=True, silent=True) or {}
    rows_raw = data.get("rows")
    if not isinstance(rows_raw, list) or len(rows_raw) == 0:
        return jsonify(error="No montage rows supplied."), 400
    try:
        rows = [int(r) for r in rows_raw]
    except (TypeError, ValueError):
        return jsonify(error="Invalid row indices."), 400
    if len(rows) > 64:
        return jsonify(error="At most 64 montage cells are supported."), 400
    mode = str(data.get("mode") or "static").strip().lower()
    try:
        if mode == "static":
            blobs, cache_token = generate_montage_volume_pngs(e, rows)
            fmt = "png"
            b64s = [base64.standard_b64encode(b).decode("ascii") for b in blobs]
            return jsonify(
                ok=True,
                format=fmt,
                images=b64s,
                rows=rows,
                volume_cache_id=cache_token,
            )
        if mode in ("animate", "gif"):
            cache_id = data.get("volume_cache_id")
            if not cache_id or not isinstance(cache_id, str):
                return (
                    jsonify(
                        error="Animate requires volume_cache_id from Generate volumes.",
                    ),
                    400,
                )
            raw_ci = data.get("cell_index")
            if raw_ci is None:
                return jsonify(error="cell_index is required for animate."), 400
            try:
                cell_index = int(raw_ci)
            except (TypeError, ValueError):
                return jsonify(error="cell_index must be an integer."), 400
            gf = int(data.get("gif_frames", DEFAULT_GIF_FRAMES))
            cc = int(data.get("chimerax_cpus", DEFAULT_CHIMERAX_PARALLEL))
            gif_bytes = volume_cell_gif_from_cache(
                cache_id,
                cell_index,
                rows_expected=tuple(rows),
                gif_frames=gf,
                chimerax_cpus=cc,
            )
            b64_one = base64.standard_b64encode(gif_bytes).decode("ascii")
            return jsonify(
                ok=True,
                format="gif",
                image=b64_one,
                cell_index=cell_index,
                rows=rows,
            )
        return jsonify(error='mode must be "static" or "animate".'), 400
    except EnvironmentError as err:
        return jsonify(error=str(err), need_chimerax=True), 503
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("explorer volume media failed")
        return jsonify(error=str(err)), 500


def api_scatter():
    e: DashboardExperiment = g.dashboard_exp
    xcol = request.args.get("x", e.numeric_columns[0])
    ycol = request.args.get("y", e.numeric_columns[0])
    ccol = request.args.get("color") or "none"
    filter_ui = request.args.get("filter_ui") == "1"
    full = request.args.get("full") == "1"
    if xcol not in e.plot_df.columns or ycol not in e.plot_df.columns:
        return jsonify(error="bad axis column"), 400
    if ccol != "none" and ccol not in e.plot_df.columns:
        return jsonify(error="bad color column"), 400
    if filter_ui:
        max_pts = _filter_ui_scatter_max_points()
    elif full:
        max_pts = None
    else:
        max_pts = 200_000
    preselect_rows, pre_err = _parse_preselect_rows_param(
        request.args.get("preselect_rows"),
    )
    if pre_err:
        return jsonify(error=pre_err), 400
    marker_size = 4.0
    raw_ms = request.args.get("marker_size")
    if raw_ms:
        try:
            marker_size = max(0.5, min(float(raw_ms), 20))
        except ValueError:
            pass
    try:
        js = scatter_json(
            e,
            xcol,
            ycol,
            None if ccol == "none" else ccol,
            max_points=max_pts,
            preselect_plot_df_rows=preselect_rows,
            use_webgl=False,
            marker_size=marker_size,
            continuous_palette=request.args.get("palette"),
        )
    except Exception as err:
        logger.exception("scatter plot failed")
        return jsonify(error=str(err)), 500
    return Response(js, mimetype="application/json")


def latent_3d_page():
    e: DashboardExperiment = g.dashboard_exp
    zdim = int(e.z.shape[1])
    if zdim < 3:
        return (
            render_template(
                "pair_grid_need_more_cols.html",
                kind="z3",
                n=zdim,
            ),
            200,
        )
    z_cols = [f"z{i}" for i in range(zdim)]
    cols = e.numeric_columns
    return render_template(
        "latent_3d.html",
        z_cols=z_cols,
        numeric_cols=cols,
        covariate_display_map={c: _covariate_display_name(c) for c in cols},
        default_x="z0",
        default_y="z1",
        default_z="z2",
    )


def api_scatter3d_z():
    e: DashboardExperiment = g.dashboard_exp
    xcol = request.args.get("x", "z0")
    ycol = request.args.get("y", "z1")
    zcol = request.args.get("z", "z2")
    ccol = request.args.get("color") or "none"
    if ccol != "none" and ccol not in e.plot_df.columns:
        return jsonify(error="bad color column"), 400
    try:
        js = scatter3d_z_json(
            e,
            xcol,
            ycol,
            zcol,
            None if ccol == "none" else ccol,
            continuous_palette=request.args.get("palette"),
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("3-D Latent Space Visualizer failed")
        return jsonify(error=str(err)), 500
    return Response(js, mimetype="application/json")


def api_latent3d_preview_png():
    """PNG snapshot of the 3D latent scatter (matplotlib), same data rules as ``api_scatter3d_z``."""
    e: DashboardExperiment = g.dashboard_exp
    xcol = request.args.get("x", "z0")
    ycol = request.args.get("y", "z1")
    zcol = request.args.get("z", "z2")
    ccol = request.args.get("color") or "none"
    try:
        elev = float(request.args.get("elev", "22"))
        azim = float(request.args.get("azim", "-65"))
    except ValueError:
        return jsonify(error="elev/azim must be numeric"), 400
    if ccol != "none" and ccol not in e.plot_df.columns:
        return jsonify(error="bad color column"), 400
    try:
        png = scatter3d_z_preview_png(
            e,
            xcol,
            ycol,
            zcol,
            None if ccol == "none" else ccol,
            continuous_palette=request.args.get("palette"),
            elev=elev,
            azim=azim,
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("3-D latent preview PNG failed")
        return jsonify(error=str(err)), 500
    return Response(png, mimetype="image/png")


def api_preview_montage():
    e: DashboardExperiment = g.dashboard_exp
    raw = request.args.get("rows", "")
    if not raw.strip():
        data = _montage_bytes(e, [])
    else:
        parts = [p.strip() for p in raw.split(",") if p.strip()]
        try:
            idxs = [int(p) for p in parts]
        except ValueError:
            return jsonify(error="rows must be integers"), 400
        data = _montage_bytes(e, idxs)
    return Response(data, mimetype="image/png")


def api_preload_images():
    """Return a subsample of particle images as base64 JPEGs for the viewer.

    Selection: ~2/3 random, ~1/6 high mean distance to reference points,
    ~1/6 high nearest-neighbor distance (see ``_sample_plot_df_rows_for_preload``).

    Use **POST** with a JSON body when ``selected_rows`` is large (lasso selections);
    query strings hit proxy/server URI length limits (414).
    """
    e: DashboardExperiment = g.dashboard_exp
    if not e.can_preview_particles:
        return jsonify(rows=[], images=[], elapsed=0)

    cols = e.plot_df.columns
    restrict_list: list[int] | None = None

    if request.method == "POST":
        data = request.get_json(force=True, silent=True) or {}
        xcol = str(data.get("x") or "")
        ycol = str(data.get("y") or "")
        raw_list = data.get("selected_rows")
        if raw_list is None:
            restrict_list = None
        elif not isinstance(raw_list, list):
            return jsonify(error="selected_rows must be a JSON array of integers."), 400
        else:
            try:
                restrict_list = [int(r) for r in raw_list]
            except (TypeError, ValueError):
                return (
                    jsonify(error="selected_rows must be a JSON array of integers."),
                    400,
                )
            if not restrict_list:
                restrict_list = None
    else:
        xcol = request.args.get("x", "")
        ycol = request.args.get("y", "")
        raw_sel = (request.args.get("selected_rows") or "").strip()
        if raw_sel:
            try:
                restrict_list = [
                    int(x.strip()) for x in raw_sel.split(",") if x.strip()
                ]
            except ValueError:
                return (
                    jsonify(error="selected_rows must be comma-separated integers."),
                    400,
                )
            if not restrict_list:
                restrict_list = None

    if xcol not in cols:
        xcol = str(cols[0])
    if ycol not in cols:
        ycol = str(cols[min(1, len(cols) - 1)])
    sel_key: tuple[int, ...] | None = (
        tuple(sorted(set(restrict_list))) if restrict_list else None
    )

    key = (e.epoch, e.kmeans_folder_id, xcol, ycol, sel_key)
    if key in _PRELOAD_CACHE:
        rows, imgs, elapsed = _PRELOAD_CACHE[key]
        return jsonify(rows=rows, images=imgs, elapsed=elapsed)

    t0 = time.monotonic()
    cpus = int(current_app.config.get("PRELOAD_CPUS") or 4)
    rows, global_indices = _sample_plot_df_rows_for_preload(
        e, xcol, ycol, restrict_to_rows=restrict_list
    )
    if not global_indices:
        return jsonify(rows=[], images=[], elapsed=0.0)

    if cpus > 1 and len(global_indices) > cpus:
        from concurrent.futures import ProcessPoolExecutor

        chunk_sz = -(-len(global_indices) // cpus)
        chunks = [
            global_indices[i : i + chunk_sz]
            for i in range(0, len(global_indices), chunk_sz)
        ]
        with ProcessPoolExecutor(max_workers=len(chunks)) as pool:
            futures = [
                pool.submit(
                    _encode_particle_batch,
                    e.particles_path,
                    e.datadir,
                    ch,
                    96,
                )
                for ch in chunks
            ]
            imgs: list[str] = []
            for f in futures:
                imgs.extend(f.result())
    else:
        imgs = _encode_particle_batch(
            e.particles_path,
            e.datadir,
            global_indices,
            96,
        )

    elapsed = round(time.monotonic() - t0, 1)
    _PRELOAD_CACHE[key] = (rows, imgs, elapsed)
    return jsonify(rows=rows, images=imgs, elapsed=elapsed)


def pairplot_page():
    e: DashboardExperiment = g.dashboard_exp
    zdim = int(e.z.shape[1])
    if zdim < 2:
        return (
            render_template(
                "pair_grid_need_more_cols.html",
                kind="zdim",
                n=zdim,
            ),
            200,
        )
    has_pc = _has_pc_columns(e)
    has_umap = _has_umap_columns(e)
    if not has_pc:
        return (
            render_template(
                "pair_grid_need_more_cols.html",
                kind="pca",
                n=zdim,
            ),
            200,
        )
    z_names = {f"z{i}" for i in range(zdim)}
    color_choices = [c for c in e.numeric_columns if c not in z_names]
    if not color_choices:
        return (
            render_template(
                "pair_grid_need_more_cols.html",
                kind="numeric",
                n=0,
            ),
            200,
        )
    default_color = next(
        (c for c in ("labels", "znorm", "UMAP1", "PC1") if c in color_choices),
        color_choices[0],
    )
    skeleton_cells = pair_grid_skeleton_placeholder_layout(zdim)
    return render_template(
        "pair_grid.html",
        color_choices=color_choices,
        covariate_display_map={c: _covariate_display_name(c) for c in color_choices},
        default_color=default_color,
        has_umap=has_umap,
        zdim=zdim,
        skeleton_placeholder_cells=skeleton_cells,
        pairplot_save_default_name="zdim_pairplot.png",
        pairplot_save_default_dir=os.path.join(e.workdir, f"analyze.{e.epoch}"),
    )


def api_pairplot():
    e: DashboardExperiment = g.dashboard_exp
    payload = request.get_json(force=True, silent=True) or {}
    color_col = payload.get("color_col") or payload.get("lower_color_col")
    if not color_col or not isinstance(color_col, str):
        return jsonify(error="Choose a color covariate."), 400
    if color_col not in e.plot_df.columns:
        return jsonify(error="Invalid color column."), 400
    if color_col in {f"z{i}" for i in range(int(e.z.shape[1]))}:
        return (
            jsonify(error="Latent z columns cannot be used as the color covariate."),
            400,
        )
    raw_diag = payload.get("diagonal_emb")
    if raw_diag is None or (isinstance(raw_diag, str) and raw_diag.strip() == ""):
        diagonal_emb = "umap" if _has_umap_columns(e) else "pc"
    else:
        diagonal_emb = str(raw_diag).lower()
    upper_style = (payload.get("upper_style") or "scatter").lower()
    raw_palette = payload.get("palette")
    pair_palette = normalize_continuous_palette(
        str(raw_palette) if raw_palette is not None else None,
    )
    if diagonal_emb not in ("pc", "umap"):
        return jsonify(error="diagonal_emb must be pc or umap."), 400
    if upper_style not in ("scatter", "hex"):
        return jsonify(error="upper_style must be scatter or hex."), 400
    if diagonal_emb == "umap" and not _has_umap_columns(e):
        return jsonify(error="UMAP is not available for this run."), 400
    if diagonal_emb == "pc" and not _has_pc_columns(e):
        return jsonify(error="PCA components are not available."), 400
    try:
        png, cells = pair_grid_png(
            e,
            lower_color_col=color_col,
            diagonal_emb=diagonal_emb,
            upper_style=upper_style,
            continuous_palette=pair_palette,
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("pair grid failed")
        return jsonify(error=str(err)), 500
    b64 = base64.standard_b64encode(png).decode("ascii")
    return jsonify(png_b64=b64, cells=cells)


def api_save_pairplot_png():
    e: DashboardExperiment = g.dashboard_exp
    payload = request.get_json(force=True, silent=True) or {}
    color_col = payload.get("color_col") or payload.get("lower_color_col")
    if not color_col or not isinstance(color_col, str):
        return jsonify(error="Choose a color covariate."), 400
    if color_col not in e.plot_df.columns:
        return jsonify(error="Invalid color column."), 400
    if color_col in {f"z{i}" for i in range(int(e.z.shape[1]))}:
        return (
            jsonify(error="Latent z columns cannot be used as the color covariate."),
            400,
        )
    raw_diag = payload.get("diagonal_emb")
    if raw_diag is None or (isinstance(raw_diag, str) and raw_diag.strip() == ""):
        diagonal_emb = "umap" if _has_umap_columns(e) else "pc"
    else:
        diagonal_emb = str(raw_diag).lower()
    upper_style = (payload.get("upper_style") or "scatter").lower()
    if diagonal_emb not in ("pc", "umap"):
        return jsonify(error="diagonal_emb must be pc or umap."), 400
    if upper_style not in ("scatter", "hex"):
        return jsonify(error="upper_style must be scatter or hex."), 400
    if diagonal_emb == "umap" and not _has_umap_columns(e):
        return jsonify(error="UMAP is not available for this run."), 400
    if diagonal_emb == "pc" and not _has_pc_columns(e):
        return jsonify(error="PCA components are not available."), 400

    raw_palette_save = payload.get("palette")
    pair_palette_save = normalize_continuous_palette(
        str(raw_palette_save) if raw_palette_save is not None else None,
    )

    raw_name = str(payload.get("filename") or "zdim_pairplot.png").strip()
    filename = os.path.basename(raw_name) or "zdim_pairplot.png"
    if not filename.lower().endswith(".png"):
        filename += ".png"
    out_dir = os.path.join(e.workdir, f"analyze.{e.epoch}")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, filename)
    try:
        png, _ = pair_grid_png(
            e,
            lower_color_col=color_col,
            diagonal_emb=diagonal_emb,
            upper_style=upper_style,
            continuous_palette=pair_palette_save,
        )
        with open(out_path, "wb") as fh:
            fh.write(png)
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except OSError as err:
        logger.exception("save pairplot failed")
        return jsonify(error=str(err)), 500
    except Exception as err:
        logger.exception("save pairplot failed")
        return jsonify(error=str(err)), 500
    return jsonify(ok=True, path=out_path, filename=filename)


def api_default_trajectory_endpoints():
    """Return start/end in plot space along the long axis of the point cloud (2D PCA).

    The segment passes through the centroid and spans from the minimum to the maximum
    projection of points onto the first principal direction — a line through the
    middle of the mass with endpoints at opposite extents of the occupied region.
    """
    e: DashboardExperiment = g.dashboard_exp
    if not explorer_volumes_eligible(e):
        return (
            jsonify(
                error="Trajectory creator needs a CUDA GPU, single-particle data, "
                "and weights for the current epoch.",
            ),
            400,
        )
    xcol = request.args.get("x", "")
    ycol = request.args.get("y", "")
    if xcol not in e.plot_df.columns or ycol not in e.plot_df.columns:
        return jsonify(error="bad axis column"), 400
    try:
        _validate_trajectory_plot_axes(e, xcol, ycol)
    except ValueError as err:
        return jsonify(error=str(err)), 400

    sub = e.plot_df[[xcol, ycol]].dropna()
    if len(sub) < 2:
        return jsonify(error="not enough finite points for default trajectory"), 400
    xy = sub.values.astype(np.float64)
    mu = xy.mean(axis=0)
    xc = xy - mu
    span = xc.max(axis=0) - xc.min(axis=0)
    if not np.any(np.isfinite(span)) or float(np.nanmax(span)) < 1e-15:
        return jsonify(
            ok=True,
            start=mu.tolist(),
            end=(mu + np.array([1e-6, 0.0])).tolist(),
        )

    if float(span[0]) >= float(span[1]):
        v = np.array([1.0, 0.0])
    else:
        v = np.array([0.0, 1.0])

    try:
        _u, _s, vt = np.linalg.svd(xc, full_matrices=False)
        if vt.shape[0] >= 1:
            cand = vt[0].astype(np.float64)
            nrm = float(np.linalg.norm(cand))
            if nrm > 1e-15:
                v = cand / nrm
    except np.linalg.LinAlgError:
        pass

    t = xc @ v
    t_min = float(np.min(t))
    t_max = float(np.max(t))
    if np.isclose(t_min, t_max):
        bump = max(float(np.nanmax(span)), 1.0) * 0.05
        start = (mu + (t_min - bump) * v).tolist()
        end = (mu + (t_max + bump) * v).tolist()
    else:
        start = (mu + t_min * v).tolist()
        end = (mu + t_max * v).tolist()
    return jsonify(ok=True, start=start, end=end)


def trajectory_creator_page():
    e: DashboardExperiment = g.dashboard_exp
    if not explorer_volumes_eligible(e):
        return (
            render_template(
                "no_images.html",
                reason="Trajectory creator needs a CUDA GPU and model weights for the current epoch.",
            ),
            200,
        )
    zdim = int(e.z.shape[1])
    traj_cols = _trajectory_plot_axis_columns(e)
    color_cols = e.numeric_columns
    dx, dy = _trajectory_default_xy_cols(traj_cols, zdim)
    cov_keys = list(dict.fromkeys(traj_cols + color_cols))
    return render_template(
        "trajectory_creator.html",
        traj_axis_cols=traj_cols,
        numeric_cols=color_cols,
        covariate_display_map={c: _covariate_display_name(c) for c in cov_keys},
        default_x=dx,
        default_y=dy,
        zdim=zdim,
        exp_workdir=e.workdir,
    )


def api_trajectory_save_zpath():
    """Write z-path text to a server-side path (defaults to ``workdir/z-path.txt``)."""
    e: DashboardExperiment = g.dashboard_exp
    if not explorer_volumes_eligible(e):
        return (
            jsonify(
                error="Trajectory creator needs a CUDA GPU, single-particle data, "
                "and weights for the current epoch.",
            ),
            400,
        )
    data = request.get_json(force=True, silent=True) or {}
    txt = data.get("z_path_txt")
    if not isinstance(txt, str):
        return jsonify(error="z_path_txt must be a string"), 400
    raw_out_path = data.get("out_path")
    if raw_out_path is None or str(raw_out_path).strip() == "":
        out_path = os.path.join(os.path.abspath(e.workdir), "z-path.txt")
    else:
        if not isinstance(raw_out_path, str):
            return jsonify(error="out_path must be a string path"), 400
        req_path = raw_out_path.strip()
        if not req_path:
            return jsonify(error="out_path must not be empty"), 400
        if not req_path.lower().endswith(".txt"):
            req_path = req_path + ".txt"
        if os.path.isabs(req_path):
            out_path = os.path.abspath(req_path)
        else:
            out_path = os.path.abspath(os.path.join(e.workdir, req_path))
    out_dir = os.path.dirname(out_path) or os.path.abspath(e.workdir)
    try:
        os.makedirs(out_dir, exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(txt)
    except OSError as err:
        return jsonify(error=str(err)), 500
    return jsonify(ok=True, path=out_path)


def api_trajectory_save_volumes():
    """Save cached trajectory ``.mrc`` files into a chosen server-side folder."""
    e: DashboardExperiment = g.dashboard_exp
    if not explorer_volumes_eligible(e):
        return (
            jsonify(
                error="Trajectory creator needs a CUDA GPU, single-particle data, "
                "and weights for the current epoch.",
            ),
            400,
        )
    data = request.get_json(force=True, silent=True) or {}
    token = str(data.get("volume_cache_id", "") or "").strip()
    out_dir = str(data.get("out_dir", "") or "").strip()
    if not token:
        return jsonify(error="Missing volume_cache_id. Generate volumes first."), 400
    if not out_dir:
        return jsonify(error="Choose an output folder."), 400
    try:
        saved = save_cached_volumes_to_dir(
            token,
            out_dir,
            filename_prefix="trajectory_volume",
        )
        return jsonify(
            ok=True, out_dir=os.path.abspath(out_dir), files=saved, n_saved=len(saved)
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except OSError as err:
        return jsonify(error=str(err)), 500


def api_trajectory_import_anchors():
    """Import anchor indices from a server-side ``.txt`` path."""
    e: DashboardExperiment = g.dashboard_exp
    if not explorer_volumes_eligible(e):
        return (
            jsonify(
                error="Trajectory creator needs a CUDA GPU, single-particle data, "
                "and weights for the current epoch.",
            ),
            400,
        )
    data = request.get_json(force=True, silent=True) or {}
    server_path = str(data.get("server_path", "") or "").strip()
    if not server_path:
        return jsonify(error="no file path provided"), 400
    if not server_path.lower().endswith(".txt"):
        return jsonify(error="select a .txt file"), 400
    abs_path = os.path.abspath(server_path)
    if not os.path.isfile(abs_path):
        return jsonify(error="file not found on server"), 400
    xcol = str(data.get("x", "") or "")
    ycol = str(data.get("y", "") or "")
    if xcol not in e.plot_df.columns or ycol not in e.plot_df.columns:
        return jsonify(error="bad axis column"), 400
    try:
        _validate_trajectory_plot_axes(e, xcol, ycol)
    except ValueError as err:
        return jsonify(error=str(err)), 400
    mode = str(data.get("mode", "direct") or "direct")
    if mode.strip().lower() == "direct":
        n_points = _parse_traj_interpolation_value(data, key="n_points", default=0)
    else:
        n_points = _parse_traj_points_value(data, default=4)
    max_neighbors = _parse_traj_neighbor_value(data, "max_neighbors", default=n_points)
    avg_neighbors = _parse_traj_neighbor_value(data, "avg_neighbors", default=n_points)
    try:
        raw = Path(abs_path).read_bytes()
        anchor_indices = _parse_anchor_indices_txt(raw)
        return jsonify(
            _trajectory_anchor_payload_from_indices(
                e,
                anchor_indices,
                xcol,
                ycol,
                mode=mode,
                n_points=n_points,
                max_neighbors=max_neighbors,
                avg_neighbors=avg_neighbors,
            )
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("trajectory anchor import failed")
        return jsonify(error=str(err)), 500


def api_list_server_files():
    """List directories and ``.txt`` files for the server-side file browser."""
    e: DashboardExperiment = g.dashboard_exp
    req_dir = request.args.get("dir", "").strip()
    root = os.path.abspath(e.workdir)
    if req_dir:
        browse = os.path.abspath(req_dir)
    else:
        analyze_dir = os.path.join(root, f"analyze.{e.epoch}")
        browse = analyze_dir if os.path.isdir(analyze_dir) else root
    if not os.path.isdir(browse):
        return jsonify(error="directory not found"), 400
    entries: list[dict] = []
    try:
        for name in sorted(os.listdir(browse)):
            full = os.path.join(browse, name)
            if os.path.isdir(full):
                entries.append({"name": name, "type": "dir"})
            elif name.lower().endswith(".txt"):
                entries.append({"name": name, "type": "file"})
    except PermissionError:
        return jsonify(error="permission denied"), 403
    parent = os.path.dirname(browse) if browse != "/" else None
    return jsonify(ok=True, dir=browse, parent=parent, entries=entries)


def api_trajectory_kmeans_centers():
    """Load ``kmeansK/centers_ind.txt`` for the current ``analyze.N`` / k-means folder."""
    e: DashboardExperiment = g.dashboard_exp
    if not explorer_volumes_eligible(e):
        return (
            jsonify(
                error="Trajectory creator needs a CUDA GPU, single-particle data, "
                "and weights for the current epoch.",
            ),
            400,
        )
    data = request.get_json(force=True, silent=True) or {}
    xcol = str(data.get("x", "") or "")
    ycol = str(data.get("y", "") or "")
    if xcol not in e.plot_df.columns or ycol not in e.plot_df.columns:
        return jsonify(error="bad axis column"), 400
    try:
        _validate_trajectory_plot_axes(e, xcol, ycol)
    except ValueError as err:
        return jsonify(error=str(err)), 400
    mode = str(data.get("mode", "direct") or "direct")
    if mode.strip().lower() == "direct":
        n_points = _parse_traj_interpolation_value(data, key="n_points", default=0)
    else:
        n_points = _parse_traj_points_value(data, default=4)
    max_neighbors = _parse_traj_neighbor_value(data, "max_neighbors", default=n_points)
    avg_neighbors = _parse_traj_neighbor_value(data, "avg_neighbors", default=n_points)
    try:
        anchor_indices = _load_kmeans_center_indices(e)
        return jsonify(
            _trajectory_anchor_payload_from_indices(
                e,
                anchor_indices,
                xcol,
                ycol,
                mode=mode,
                n_points=n_points,
                max_neighbors=max_neighbors,
                avg_neighbors=avg_neighbors,
            )
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("trajectory k-means centers failed")
        return jsonify(error=str(err)), 500


def api_trajectory_random_indices():
    """Choose up to 10 random dataset indices (fewer if the stack is smaller)."""
    e: DashboardExperiment = g.dashboard_exp
    if not explorer_volumes_eligible(e):
        return (
            jsonify(
                error="Trajectory creator needs a CUDA GPU, single-particle data, "
                "and weights for the current epoch.",
            ),
            400,
        )
    data = request.get_json(force=True, silent=True) or {}
    xcol = str(data.get("x", "") or "")
    ycol = str(data.get("y", "") or "")
    if xcol not in e.plot_df.columns or ycol not in e.plot_df.columns:
        return jsonify(error="bad axis column"), 400
    try:
        _validate_trajectory_plot_axes(e, xcol, ycol)
    except ValueError as err:
        return jsonify(error=str(err)), 400
    mode = str(data.get("mode", "direct") or "direct")
    if mode.strip().lower() == "direct":
        n_points = _parse_traj_interpolation_value(data, key="n_points", default=0)
    else:
        n_points = _parse_traj_points_value(data, default=4)
    max_neighbors = _parse_traj_neighbor_value(data, "max_neighbors", default=n_points)
    avg_neighbors = _parse_traj_neighbor_value(data, "avg_neighbors", default=n_points)
    try:
        anchor_indices = _random_dataset_indices(e, 10)
        return jsonify(
            _trajectory_anchor_payload_from_indices(
                e,
                anchor_indices,
                xcol,
                ycol,
                mode=mode,
                n_points=n_points,
                max_neighbors=max_neighbors,
                avg_neighbors=avg_neighbors,
            )
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("trajectory random indices failed")
        return jsonify(error=str(err)), 500


def api_trajectory_coords():
    """Latent z along the trajectory line (no ChimeraX / volume decode)."""
    e: DashboardExperiment = g.dashboard_exp
    if not explorer_volumes_eligible(e):
        return (
            jsonify(
                error="Trajectory creator needs a CUDA GPU, single-particle data, "
                "and weights for the current epoch.",
            ),
            400,
        )
    data = request.get_json(force=True, silent=True) or {}
    try:
        p = _parse_trajectory_request_body(e, data)
    except ValueError as err:
        return jsonify(error=str(err)), 400
    try:
        z_traj, traj_rows, traj_xy = _compute_trajectory_latent_path(e, p)
        payload = _trajectory_shared_json_payload(
            e,
            z_traj,
            traj_rows,
            traj_xy,
            mode=p["mode"],
            n_points=p["n_points"],
            xcol=p["xcol"],
            ycol=p["ycol"],
        )
        if p.get("use_anchors") and p["mode"] == "direct":
            pidx = _direct_anchor_particle_indices_payload(
                anchor_indices=p["anchor_indices"],
                interpolation_points=p["n_points"],
                n_total=int(np.asarray(z_traj).shape[0]),
            )
            if pidx is not None:
                payload["traj_particle_indices"] = pidx
        return jsonify(payload)
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("trajectory coords failed")
        return jsonify(error=str(err)), 500


def api_trajectory_volumes():
    e: DashboardExperiment = g.dashboard_exp
    if not explorer_volumes_eligible(e):
        return (
            jsonify(
                error="Trajectory creator needs a CUDA GPU, single-particle data, "
                "and weights for the current epoch.",
            ),
            400,
        )
    data = request.get_json(force=True, silent=True) or {}
    try:
        p = _parse_trajectory_request_body(e, data)
    except ValueError as err:
        return jsonify(error=str(err)), 400

    try:
        z_traj, traj_rows, traj_xy = _compute_trajectory_latent_path(e, p)
        mode = p["mode"]
        n_points = p["n_points"]
        xcol = p["xcol"]
        ycol = p["ycol"]

        blobs, cache_token = generate_trajectory_volume_pngs(e, z_traj)
        b64s = [base64.standard_b64encode(b).decode("ascii") for b in blobs]
        payload = _trajectory_shared_json_payload(
            e,
            z_traj,
            traj_rows,
            traj_xy,
            mode=mode,
            n_points=n_points,
            xcol=xcol,
            ycol=ycol,
        )
        if p.get("use_anchors") and mode == "direct":
            pidx = _direct_anchor_particle_indices_payload(
                anchor_indices=p["anchor_indices"],
                interpolation_points=p["n_points"],
                n_total=int(np.asarray(z_traj).shape[0]),
            )
            if pidx is not None:
                payload["traj_particle_indices"] = pidx
        payload["images"] = b64s
        payload["volume_cache_id"] = cache_token
        if traj_rows is not None and mode in ("nearest", "graph"):
            payload["particle_thumbs"] = [
                _particle_thumbnail_b64_from_row(e, int(r)) for r in traj_rows
            ]
        return jsonify(payload)
    except EnvironmentError as err:
        return jsonify(error=str(err), need_chimerax=True), 503
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("trajectory volume generation failed")
        return jsonify(error=str(err)), 500


def create_app(
    workdir: str | None,
    epoch: int = -1,
    kmeans: int = -1,
    filter_plot_inds: str | None = None,
    cpus: int = 4,
) -> Flask:
    app = Flask(
        __name__,
        template_folder=_TEMPLATE_DIR,
        static_folder=_STATIC_DIR,
        static_url_path="/static",
    )
    app.jinja_env.filters["abbrev_middle"] = abbrev_middle
    app.config["PRELOAD_CPUS"] = max(1, cpus)
    app.secret_key = os.environ.get(
        "CRYODRGN_DASHBOARD_SECRET", "cryodrgn-dashboard-dev-key"
    )
    command_builder_only = workdir is None
    app.config["COMMAND_BUILDER_ONLY"] = command_builder_only
    app.config["DASHBOARD_DISCOVERY_CWD"] = os.getcwd()
    discovered = (
        _discover_cryodrgn_workdirs(os.getcwd()) if command_builder_only else []
    )
    app.config["DASHBOARD_DISCOVERED_WORKDIRS"] = discovered
    app.config["DASHBOARD_DISCOVERY_BOOT_ID"] = (
        str(uuid.uuid4()) if command_builder_only and discovered else None
    )
    app.before_request(_bind_dashboard_exp)
    if command_builder_only:
        app.config["DASHBOARD_WORKDIR"] = None
        app.config["DASHBOARD_KMEANS"] = -1
        app.config["DASHBOARD_EPOCHS"] = [0]
        app.config["DASHBOARD_START_EPOCH"] = 0
        app.config["FILTER_PLOT_INDS"] = None
        app.context_processor(inject_meta_command_builder_only)
    else:
        workdir = os.path.abspath(workdir)
        epochs = list_z_epochs(workdir)
        if not epochs:
            raise ValueError(
                f"No analyzed epochs under {workdir!r} — need z.N.pkl and analyze.N/ "
                "(run `cryodrgn analyze` first)."
            )
        start_epoch = epoch if epoch != -1 else max(epochs)
        if start_epoch not in epochs:
            start_epoch = max(epochs)
        app.config["DASHBOARD_WORKDIR"] = workdir
        app.config["DASHBOARD_KMEANS"] = kmeans
        app.config["DASHBOARD_EPOCHS"] = epochs
        app.config["DASHBOARD_START_EPOCH"] = start_epoch
        app.config["FILTER_PLOT_INDS"] = filter_plot_inds
        app.context_processor(inject_meta)
    app.add_url_rule("/api/set_epoch", view_func=api_set_epoch, methods=["POST"])
    app.add_url_rule("/api/set_workdir", view_func=api_set_workdir, methods=["POST"])
    app.add_url_rule("/", view_func=index, methods=["GET"])
    app.add_url_rule(
        "/command-builder", view_func=command_builder_page, methods=["GET"]
    )
    app.add_url_rule(
        "/abinit-builder", view_func=abinit_builder_redirect, methods=["GET"]
    )
    app.add_url_rule("/filter", view_func=filter_page_redirect, methods=["GET"])
    app.add_url_rule(
        "/api/save_selection", view_func=api_save_selection, methods=["POST"]
    )
    app.add_url_rule("/explorer", view_func=explorer, methods=["GET"])
    app.add_url_rule(
        "/api/explorer_volume_media",
        view_func=api_explorer_volume_media,
        methods=["POST"],
    )
    app.add_url_rule("/api/scatter", view_func=api_scatter, methods=["GET"])
    app.add_url_rule("/latent-3d", view_func=latent_3d_page, methods=["GET"])
    app.add_url_rule("/api/scatter3d_z", view_func=api_scatter3d_z, methods=["GET"])
    app.add_url_rule(
        "/api/latent3d_preview.png",
        view_func=api_latent3d_preview_png,
        methods=["GET"],
    )
    app.add_url_rule(
        "/api/preview_montage", view_func=api_preview_montage, methods=["GET"]
    )
    app.add_url_rule(
        "/api/preload_images", view_func=api_preload_images, methods=["GET", "POST"]
    )
    app.add_url_rule("/pairplot", view_func=pairplot_page, methods=["GET"])
    app.add_url_rule("/api/pairplot", view_func=api_pairplot, methods=["POST"])
    app.add_url_rule(
        "/api/save_pairplot_png", view_func=api_save_pairplot_png, methods=["POST"]
    )
    app.add_url_rule("/trajectory", view_func=trajectory_creator_page, methods=["GET"])
    app.add_url_rule(
        "/api/trajectory_volumes",
        view_func=api_trajectory_volumes,
        methods=["POST"],
    )
    app.add_url_rule(
        "/api/trajectory_coords",
        view_func=api_trajectory_coords,
        methods=["POST"],
    )
    app.add_url_rule(
        "/api/trajectory_save_zpath",
        view_func=api_trajectory_save_zpath,
        methods=["POST"],
    )
    app.add_url_rule(
        "/api/trajectory_save_volumes",
        view_func=api_trajectory_save_volumes,
        methods=["POST"],
    )
    app.add_url_rule(
        "/api/trajectory_import_anchors",
        view_func=api_trajectory_import_anchors,
        methods=["POST"],
    )
    app.add_url_rule(
        "/api/list_server_files",
        view_func=api_list_server_files,
        methods=["GET"],
    )
    app.add_url_rule(
        "/api/trajectory_kmeans_centers",
        view_func=api_trajectory_kmeans_centers,
        methods=["POST"],
    )
    app.add_url_rule(
        "/api/trajectory_random_indices",
        view_func=api_trajectory_random_indices,
        methods=["POST"],
    )
    app.add_url_rule(
        "/api/default_trajectory_endpoints",
        view_func=api_default_trajectory_endpoints,
        methods=["GET"],
    )

    return app


def run_server(
    workdir: str | None,
    epoch: int = -1,
    kmeans: int = -1,
    plot_inds: str | None = None,
    host: str = "127.0.0.1",
    port: int = 5050,
    debug: bool = False,
    cpus: int = 4,
) -> None:
    app = create_app(
        workdir=workdir,
        epoch=epoch,
        kmeans=kmeans,
        filter_plot_inds=plot_inds,
        cpus=cpus,
    )

    app.run(host=host, port=port, debug=debug, threaded=True)
