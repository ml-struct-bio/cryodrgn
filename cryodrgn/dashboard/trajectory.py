"""Pure-logic helpers for the trajectory-creator view.

Nothing in this module imports Flask — the functions take a
:class:`DashboardExperiment` and plain Python dicts so they can be reused from
tests, notebooks, and API routes alike. The routes in :mod:`app` translate
``ValueError`` into HTTP 400 responses.
"""

from __future__ import annotations

import io
import os
import re
from heapq import heappop, heappush

import numpy as np

from cryodrgn.dashboard.data import DashboardExperiment

# (workdir, epoch, max_neighbors, avg_neighbors) -> (neighbor_ids, neighbor_dists)
_TRAJ_GRAPH_NEIGHBOR_CACHE: dict[
    tuple[str, int, int, int], tuple[np.ndarray, np.ndarray]
] = {}


# ---------------------------------------------------------------------------
# Axis selection / validation
# ---------------------------------------------------------------------------


def has_umap_columns(exp: DashboardExperiment) -> bool:
    return (
        exp.umap is not None
        and "UMAP1" in exp.plot_df.columns
        and "UMAP2" in exp.plot_df.columns
    )


def has_pc_columns(exp: DashboardExperiment) -> bool:
    return "PC1" in exp.plot_df.columns and "PC2" in exp.plot_df.columns


def trajectory_default_xy_cols(cols: list[str], zdim: int) -> tuple[str, str]:
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


def trajectory_plot_axis_columns(e: DashboardExperiment) -> list[str]:
    """Axes allowed: z0..z_{D-1}, PC1..PC_D when D>2, UMAP* when UMAP exists."""
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
    if has_umap_columns(e):
        umap_cols = [
            str(c) for c in e.plot_df.columns if re.fullmatch(r"UMAP[0-9]+", str(c))
        ]
        umap_cols.sort(key=lambda c: int(c.replace("UMAP", "")))
        for c in umap_cols:
            add(c)
    return ordered


def _trajectory_xy_ok_for_direct(xcol: str, ycol: str) -> bool:
    """Direct mode linearly interpolates z — only valid in z or PC plot space."""

    def ok_one(c: str) -> bool:
        return bool(re.fullmatch(r"z[0-9]+", c) or re.fullmatch(r"PC[0-9]+", c))

    return ok_one(xcol) and ok_one(ycol)


def validate_trajectory_plot_axes(e: DashboardExperiment, xcol: str, ycol: str) -> None:
    allowed = frozenset(trajectory_plot_axis_columns(e))
    if xcol not in allowed or ycol not in allowed:
        raise ValueError(
            "That axis combination is not available in trajectory creator."
        )


def trajectory_axes_from_payload(e: DashboardExperiment, data: dict) -> tuple[str, str]:
    xcol = str(data.get("x", "") or "")
    ycol = str(data.get("y", "") or "")
    if xcol not in e.plot_df.columns or ycol not in e.plot_df.columns:
        raise ValueError("bad axis column")
    validate_trajectory_plot_axes(e, xcol, ycol)
    return xcol, ycol


# ---------------------------------------------------------------------------
# Integer parsing (clamped) — shared across routes
# ---------------------------------------------------------------------------


def parse_int_from_dict(data: dict, key: str, *, default: int, lo: int, hi: int) -> int:
    """Coerce ``data[key]`` to an int; clamp to ``[lo, hi]``; fall back to ``default``."""
    try:
        v = int(data.get(key, default))
    except (TypeError, ValueError):
        return default
    return max(lo, min(v, hi))


def parse_traj_points_value(data: dict, default: int = 4) -> int:
    """Anchor-count for graph / nearest trajectories (``[2, 20]``)."""
    return parse_int_from_dict(data, "n_points", default=default, lo=2, hi=20)


def parse_traj_interpolation_value(
    data: dict, key: str = "n_points", default: int = 4
) -> int:
    """Interpolation-count for direct trajectories (``[0, 20]``)."""
    return parse_int_from_dict(data, key, default=default, lo=0, hi=20)


def parse_traj_neighbor_value(data: dict, key: str, default: int) -> int:
    """Neighbor-count for graph trajectories (``[2, 200]``)."""
    return parse_int_from_dict(data, key, default=default, lo=2, hi=200)


def trajectory_anchor_mode_params(data: dict) -> tuple[str, int, int, int]:
    mode = str(data.get("mode", "direct") or "direct")
    if mode.strip().lower() == "direct":
        n_points = parse_traj_interpolation_value(data, key="n_points", default=0)
    else:
        n_points = parse_traj_points_value(data, default=4)
    max_neighbors = parse_traj_neighbor_value(data, "max_neighbors", default=n_points)
    avg_neighbors = parse_traj_neighbor_value(data, "avg_neighbors", default=n_points)
    return mode, n_points, max_neighbors, avg_neighbors


# ---------------------------------------------------------------------------
# Small formatting / parsing helpers
# ---------------------------------------------------------------------------


def z_traj_to_savetxt_str(z_traj: np.ndarray) -> str:
    """Same text layout as ``cryodrgn graph_traversal -o`` / ``np.savetxt``."""
    buf = io.StringIO()
    np.savetxt(buf, z_traj)
    return buf.getvalue()


def _round_direct_mode_traj_xy(traj_xy: np.ndarray) -> np.ndarray:
    """Snap axis coordinates to 2 or 3 decimal places per axis (depending on scale)."""
    out = np.asarray(traj_xy, dtype=np.float64).copy()
    for j in range(out.shape[1]):
        col = out[:, j]
        finite = col[np.isfinite(col)]
        if finite.size == 0:
            continue
        decimals = 2 if float(np.nanmax(np.abs(finite))) >= 100.0 else 3
        out[:, j] = np.round(col, decimals)
    return out


def parse_anchor_indices_txt(raw: bytes) -> list[int]:
    """Parse whitespace-delimited integer indices from a UTF-8 text blob."""
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


def _plot_row_particle_index(exp: DashboardExperiment, row_index: int) -> int:
    ri = int(row_index)
    if "index" in exp.plot_df.columns:
        return int(exp.plot_df.iloc[ri]["index"])
    return ri


# ---------------------------------------------------------------------------
# Trajectory computation
# ---------------------------------------------------------------------------


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
        z_parts.append((np.outer(1.0 - t, z_start) + np.outer(t, z_end))[:-1])
        xy_parts.append((np.outer(1.0 - t, xy_start) + np.outer(t, xy_end))[:-1])
    z_parts.append(z_anchor[-1:].copy())
    xy_parts.append(anchor_xy[-1:].copy())
    return np.concatenate(z_parts, axis=0), None, np.concatenate(xy_parts, axis=0)


def _graph_neighbor_arrays(
    e: DashboardExperiment,
    *,
    max_neighbors: int,
    avg_neighbors: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Nearest-neighbor arrays for graph traversal (cached by workdir+epoch)."""
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
    """Shortest path in a sparse directed neighbor graph (``-1``/``inf`` = no edge)."""
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
    """Match ``cryodrgn graph_traversal`` for anchor-to-anchor shortest paths."""
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


def parse_trajectory_request_body(e: DashboardExperiment, data: dict) -> dict:
    """Validate a trajectory POST JSON body (shared by coords-only and volume decode)."""
    raw_anchors = data.get("anchor_indices")
    if isinstance(raw_anchors, list) and len(raw_anchors) >= 2:
        try:
            anchor_indices = [int(x) for x in raw_anchors]
        except (TypeError, ValueError) as err:
            raise ValueError("anchor_indices must be a list of integers") from err
        xcol, ycol = trajectory_axes_from_payload(e, data)
        mode = str(data.get("mode", "direct")).strip().lower()
        if mode not in ("direct", "graph"):
            raise ValueError('anchor mode must be "direct" or "graph".')
        if mode == "direct":
            n_points = parse_traj_interpolation_value(data, key="n_points", default=0)
        else:
            n_points = parse_traj_points_value(data, default=4)
        max_neighbors = parse_traj_neighbor_value(
            data, "max_neighbors", default=max(2, n_points)
        )
        avg_neighbors = parse_traj_neighbor_value(
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
    traj_xy_custom = sx0 = sy0 = ex0 = ey0 = None
    if isinstance(raw_traj_xy, list) and len(raw_traj_xy) >= 2:
        pts: list[list[float]] = []
        for i, p in enumerate(raw_traj_xy[:200]):
            if not isinstance(p, list) or len(p) != 2:
                raise ValueError(f"traj_xy[{i}] must be [x, y].")
            try:
                pts.append([float(p[0]), float(p[1])])
            except (TypeError, ValueError) as err:
                raise ValueError(f"traj_xy[{i}] has non-numeric coordinates.") from err
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

    xcol, ycol = trajectory_axes_from_payload(e, data)
    if mode == "direct" and not _trajectory_xy_ok_for_direct(xcol, ycol):
        raise ValueError(
            "Direct traversal is only available for principal-component or z latent axes."
        )
    n_points = parse_traj_points_value(data, default=4)
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


def compute_trajectory_latent_path(
    e: DashboardExperiment, p: dict
) -> tuple[np.ndarray, list[int] | None, np.ndarray]:
    """Return ``(z_traj, traj_rows_or_None, traj_xy)`` for the parsed body ``p``."""
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
        traj_rows = [int(np.argmin(np.sum((coords - pt) ** 2, axis=1))) for pt in pts]
        z_traj = e.z[np.asarray(traj_rows, dtype=int)]
        if mode == "nearest":
            return z_traj, traj_rows, coords[np.asarray(traj_rows, dtype=int)]
        return z_traj, None, _round_direct_mode_traj_xy(pts)

    if mode == "direct":
        start_pt = np.array([sx0, sy0])
        end_pt = np.array([ex0, ey0])
        start_row = int(np.argmin(np.sum((coords - start_pt) ** 2, axis=1)))
        end_row = int(np.argmin(np.sum((coords - end_pt) ** 2, axis=1)))
        t = np.linspace(0.0, 1.0, n_points, dtype=np.float64)
        z_traj = np.outer(1.0 - t, e.z[start_row]) + np.outer(t, e.z[end_row])
        traj_xy = _round_direct_mode_traj_xy(
            np.outer(1.0 - t, start_pt) + np.outer(t, end_pt)
        )
        return z_traj, None, traj_xy

    t = np.linspace(0.0, 1.0, n_points, dtype=np.float64)
    line_xy = np.outer(1.0 - t, np.array([sx0, sy0])) + np.outer(
        t, np.array([ex0, ey0])
    )
    traj_rows = [int(np.argmin(np.sum((coords - pt) ** 2, axis=1))) for pt in line_xy]
    z_traj = e.z[np.array(traj_rows)]
    return z_traj, traj_rows, coords[np.array(traj_rows)]


# ---------------------------------------------------------------------------
# JSON payload assembly
# ---------------------------------------------------------------------------


def trajectory_shared_json_payload(
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
        "z_path_txt": z_traj_to_savetxt_str(z_traj),
        "xcol": xcol,
        "ycol": ycol,
    }
    if traj_rows is not None and mode in ("nearest", "graph"):
        payload["traj_particle_indices"] = [
            _plot_row_particle_index(e, r) for r in traj_rows
        ]
    return payload


def direct_anchor_particle_indices_payload(
    *, anchor_indices: list[int], interpolation_points: int, n_total: int
) -> list[int | None] | None:
    """Row-aligned particle IDs for direct-anchor trajectories.

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


def trajectory_anchor_payload_from_indices(
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
    z_traj, traj_rows, traj_xy = compute_trajectory_latent_path(e, p)
    payload = trajectory_shared_json_payload(
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
        pidx = direct_anchor_particle_indices_payload(
            anchor_indices=anchor_indices,
            interpolation_points=n_points,
            n_total=int(np.asarray(z_traj).shape[0]),
        )
        if pidx is not None:
            payload["traj_particle_indices"] = pidx
    payload["anchor_indices"] = anchor_indices
    return payload


# ---------------------------------------------------------------------------
# Dataset-index lookups (k-means centers, random, filter-inds)
# ---------------------------------------------------------------------------


def kmeans_centers_ind_path(e: DashboardExperiment) -> str:
    return os.path.join(
        e.workdir,
        f"analyze.{e.epoch}",
        f"kmeans{e.kmeans_folder_id}",
        "centers_ind.txt",
    )


def load_kmeans_center_indices(e: DashboardExperiment) -> list[int]:
    path = kmeans_centers_ind_path(e)
    if not os.path.isfile(path):
        raise ValueError(
            f"No k-means centers_ind.txt at {path}. "
            "Run cryodrgn analyze with k-means first."
        )
    raw = np.loadtxt(path)
    arr = np.atleast_1d(raw).astype(int).ravel()
    if arr.size < 2:
        raise ValueError("Need at least two k-means center indices")
    return arr.tolist()


def random_dataset_indices(e: DashboardExperiment, k: int = 10) -> list[int]:
    """Up to ``k`` distinct random row indices into ``z``."""
    n = int(e.z.shape[0])
    if n < 2:
        raise ValueError("Dataset has fewer than 2 particles")
    rng = np.random.default_rng()
    return rng.choice(n, size=min(int(k), n), replace=False).astype(int).tolist()


def plot_df_rows_for_dataset_indices(
    exp: DashboardExperiment, dataset_indices: np.ndarray
) -> list[int]:
    """Map saved dataset indices to ``plot_df`` row positions (``[]`` if empty)."""
    want_arr = np.asarray(dataset_indices).ravel().astype(int, copy=False)
    if want_arr.size == 0:
        return []
    mask = np.isin(np.asarray(exp.all_indices), want_arr, assume_unique=False)
    return np.nonzero(mask)[0].astype(int).tolist()


def default_trajectory_endpoints_xy(
    e: DashboardExperiment, xcol: str, ycol: str
) -> tuple[list[float], list[float]]:
    """Start/end in plot space along the long axis of the 2-D point cloud.

    The segment passes through the centroid and spans from the minimum to the
    maximum projection of points onto the first principal direction (SVD of
    the centred coordinates, with an axis-aligned fallback).
    """
    sub = e.plot_df[[xcol, ycol]].dropna()
    if len(sub) < 2:
        raise ValueError("not enough finite points for default trajectory")
    xy = sub.values.astype(np.float64)
    mu = xy.mean(axis=0)
    xc = xy - mu
    span = xc.max(axis=0) - xc.min(axis=0)
    if not np.any(np.isfinite(span)) or float(np.nanmax(span)) < 1e-15:
        return mu.tolist(), (mu + np.array([1e-6, 0.0])).tolist()

    v = (
        np.array([1.0, 0.0])
        if float(span[0]) >= float(span[1])
        else np.array([0.0, 1.0])
    )
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
        return (
            (mu + (t_min - bump) * v).tolist(),
            (mu + (t_max + bump) * v).tolist(),
        )
    return (mu + t_min * v).tolist(), (mu + t_max * v).tolist()
