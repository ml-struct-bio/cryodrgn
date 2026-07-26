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
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.palette_config import normalize_continuous_palette
from cryodrgn.dashboard.plots_color_covariate import (
    _continuous_series_stats,
    _lower_color_series_is_discrete,
    _scatter_discrete_marker_arrays,
    numeric_array_to_plotly_hex,
)

if TYPE_CHECKING:
    from scipy.sparse import csr_matrix as _csr_matrix

from cryodrgn.dashboard.experiment_store import EXPERIMENT_STORE

# (workdir, epoch, max_neighbors, avg_neighbors) -> (neighbor_ids, dists, csr_graph)
_TRAJ_GRAPH_NEIGHBOR_CACHE = EXPERIMENT_STORE.traj_graph_neighbors

_CSGRAPH_MISSING_PREDECESSOR = -9999


# ---------------------------------------------------------------------------
# Axis selection / validation
# ---------------------------------------------------------------------------


def has_umap_columns(exp: DashboardExperiment) -> bool:
    """True when UMAP coordinates were loaded into ``plot_df``."""
    return (
        exp.umap is not None
        and "UMAP1" in exp.plot_df.columns
        and "UMAP2" in exp.plot_df.columns
    )


def has_pc_columns(exp: DashboardExperiment) -> bool:
    """True when at least PC1 and PC2 exist (PCA latent embedding in the dataframe)."""
    return "PC1" in exp.plot_df.columns and "PC2" in exp.plot_df.columns


def trajectory_default_xy_cols(cols: list[str], zdim: int) -> tuple[str, str]:
    """Pick default X/Y from the trajectory-allowed axis list only."""
    from cryodrgn.dashboard.route_helpers import default_embedding_xy_cols

    return default_embedding_xy_cols(cols, zdim=zdim, prefer_pc=True)


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
    """Coerce ``data[key]`` to an int.

    Clamp to ``[lo, hi]``; fall back to ``default``.
    """
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


def resolve_trajectory_decode_plan(
    data: dict, z_traj: np.ndarray
) -> tuple[np.ndarray, list[int], int]:
    """Choose which trajectory slots to decode for a volume request.

    Returns ``(z_decode, slot_indices, n_traj)`` where ``slot_indices[j]`` is the
    full-trajectory index for ``z_decode[j]``. When ``decode_indices`` is omitted,
    all points are decoded in order.
    """
    z_traj = np.asarray(z_traj, dtype=np.float64)
    n_traj = int(z_traj.shape[0])
    raw = data.get("decode_indices")
    if raw is None:
        return z_traj, list(range(n_traj)), n_traj
    if not isinstance(raw, list) or not raw:
        raise ValueError("decode_indices must be a non-empty list of integers.")
    slot_indices = [int(i) for i in raw]
    for idx in slot_indices:
        if idx < 0 or idx >= n_traj:
            raise ValueError(f"decode_indices out of range: {idx}")
    z_decode = z_traj[np.asarray(slot_indices, dtype=np.int64)]
    if z_decode.ndim == 1:
        zdim = int(z_traj.shape[1]) if z_traj.ndim == 2 else 1
        if z_decode.size == zdim:
            z_decode = z_decode.reshape(1, zdim)
        else:
            z_decode = z_decode.reshape(-1, zdim)
    return z_decode, slot_indices, n_traj


def parse_anchor_path_order(data: dict) -> str:
    """``preserve`` (default), ``heuristic``, or ``exact``."""
    raw = str(data.get("anchor_path_order", "preserve") or "preserve").strip().lower()
    if raw in ("exact", "held_karp", "held-karp", "optimal"):
        return "exact"
    if raw in ("heuristic", "inexact"):
        return "heuristic"
    if raw in ("preserve", "none", "selection", "natural", "original"):
        return "preserve"
    return "preserve"


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


def _direct_mode_xy_atol(*pts: np.ndarray) -> float:
    """Half-unit tolerance matching :func:`_round_direct_mode_traj_xy` decimals."""
    max_abs = 0.0
    for pt in pts:
        arr = np.asarray(pt, dtype=np.float64).ravel()
        finite = arr[np.isfinite(arr)]
        if finite.size:
            max_abs = max(max_abs, float(np.nanmax(np.abs(finite))))
    decimals = 2 if max_abs >= 100.0 else 3
    return 0.5 * (10.0**-decimals)


def _xy_already_on_particle(pt: np.ndarray, particle_xy: np.ndarray) -> bool:
    """True when ``pt`` already coincides with a particle under direct-mode rounding."""
    a = np.asarray(pt, dtype=np.float64).ravel()
    b = np.asarray(particle_xy, dtype=np.float64).ravel()
    if a.size < 2 or b.size < 2:
        return False
    atol = _direct_mode_xy_atol(a[:2], b[:2])
    return bool(np.allclose(a[:2], b[:2], rtol=0.0, atol=atol))


def _snap_points_to_nearest_particles(
    coords: np.ndarray, pts: np.ndarray
) -> tuple[list[int], np.ndarray]:
    """Nearest-particle rows / XY for each sample.

    Samples that already lie on their nearest particle (within direct-mode
    rounding) keep their input XY unchanged so snap has no visual effect on
    points that already correspond to real particles.
    """
    coords = np.asarray(coords, dtype=np.float64)
    pts = np.asarray(pts, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] < 2:
        raise ValueError("pts must be an (N, 2) array")
    rows: list[int] = []
    out_xy = np.empty((pts.shape[0], 2), dtype=np.float64)
    for i, pt in enumerate(pts):
        row = _nearest_plot_row_for_xy(coords, pt)
        particle_xy = np.asarray(coords[row], dtype=np.float64)
        rows.append(row)
        if _xy_already_on_particle(pt, particle_xy):
            out_xy[i] = pt[:2]
        else:
            out_xy[i] = particle_xy[:2]
    return rows, out_xy


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
# Shortest non-crossing path through anchor scatter coordinates
# ---------------------------------------------------------------------------


def _orientation_2d(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    return float((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def _on_segment_2d(a: np.ndarray, b: np.ndarray, c: np.ndarray, eps: float) -> bool:
    return (
        min(a[0], b[0]) - eps <= c[0] <= max(a[0], b[0]) + eps
        and min(a[1], b[1]) - eps <= c[1] <= max(a[1], b[1]) + eps
    )


def _segments_cross(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> bool:
    """True when segments ``ab`` and ``cd`` intersect (including collinear overlap)."""
    eps = 1e-12
    o1 = _orientation_2d(a, b, c)
    o2 = _orientation_2d(a, b, d)
    o3 = _orientation_2d(c, d, a)
    o4 = _orientation_2d(c, d, b)
    if o1 * o2 < -eps and o3 * o4 < -eps:
        return True
    if abs(o1) <= eps and _on_segment_2d(a, c, b, eps):
        return True
    if abs(o2) <= eps and _on_segment_2d(a, d, b, eps):
        return True
    if abs(o3) <= eps and _on_segment_2d(c, a, d, eps):
        return True
    if abs(o4) <= eps and _on_segment_2d(c, b, d, eps):
        return True
    return False


def _segments_cross_proper(
    a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray
) -> bool:
    """Crossing excluding shared endpoints (adjacent path edges)."""
    eps = 1e-12
    for p in (a, b):
        for q in (c, d):
            if float(np.linalg.norm(p - q)) <= eps:
                return False
    return _segments_cross(a, b, c, d)


def _path_has_crossings(order: list[int], points: np.ndarray) -> bool:
    n = len(order)
    if n < 4:
        return False
    for i in range(n - 1):
        a = points[order[i]]
        b = points[order[i + 1]]
        for j in range(i + 2, n - 1):
            c = points[order[j]]
            d = points[order[j + 1]]
            if _segments_cross_proper(a, b, c, d):
                return True
    return False


def _path_length(order: list[int], points: np.ndarray) -> float:
    total = 0.0
    for i in range(len(order) - 1):
        total += float(np.linalg.norm(points[order[i]] - points[order[i + 1]]))
    return total


def _uncross_path_order(order: list[int], points: np.ndarray) -> list[int]:
    """Remove segment crossings by reversing subpaths (2-opt uncross)."""
    order = list(order)
    if len(order) < 4:
        return order
    changed = True
    guard = 0
    max_guard = max(64, len(order) * len(order))
    seen: set[str] = set()
    while changed:
        changed = False
        guard += 1
        if guard > max_guard:
            break
        sig = ",".join(str(i) for i in order)
        if sig in seen:
            break
        seen.add(sig)
        n = len(order)
        for i in range(n - 1):
            for j in range(i + 2, n - 1):
                a = points[order[i]]
                b = points[order[i + 1]]
                c = points[order[j]]
                d = points[order[j + 1]]
                if _segments_cross_proper(a, b, c, d):
                    order[i + 1 : j + 1] = reversed(order[i + 1 : j + 1])
                    changed = True
                    break
            if changed:
                break
    return order


def _held_karp_open_path(dist: np.ndarray) -> list[int]:
    """Shortest open Hamiltonian path (exact, ``O(n^2 2^n)``)."""
    n = int(dist.shape[0])
    if n == 0:
        return []
    if n == 1:
        return [0]
    size = 1 << n
    inf = float("inf")
    dp = np.full((size, n), inf, dtype=np.float64)
    parent = np.full((size, n), -1, dtype=np.int32)
    for i in range(n):
        dp[1 << i, i] = 0.0
    for mask in range(size):
        for j in range(n):
            if not (mask & (1 << j)):
                continue
            cost_j = dp[mask, j]
            if cost_j == inf:
                continue
            rem = (~mask) & (size - 1)
            k = rem
            while k:
                lsb = k & -k
                kk = lsb.bit_length() - 1
                nmask = mask | lsb
                nd = cost_j + dist[j, kk]
                if nd < dp[nmask, kk]:
                    dp[nmask, kk] = nd
                    parent[nmask, kk] = j
                k -= lsb
    full = size - 1
    end = int(np.argmin(dp[full]))
    order = [end]
    mask = full
    cur = end
    while len(order) < n:
        prev = int(parent[mask, cur])
        order.append(prev)
        mask ^= 1 << cur
        cur = prev
    order.reverse()
    return order


def _nearest_neighbor_open_path(dist: np.ndarray, points: np.ndarray) -> list[int]:
    n = int(dist.shape[0])
    best_order: list[int] | None = None
    best_dist = float("inf")
    for start in range(n):
        unvisited = set(range(n))
        order = [start]
        unvisited.remove(start)
        cur = start
        while unvisited:
            nxt = min(unvisited, key=lambda j: dist[cur, j])
            order.append(nxt)
            unvisited.remove(nxt)
            cur = nxt
        order = _uncross_path_order(order, points)
        d = _path_length(order, points)
        if d < best_dist:
            best_dist = d
            best_order = order
    return best_order if best_order is not None else list(range(n))


def _improve_noncrossing_path(order: list[int], points: np.ndarray) -> list[int]:
    order = list(order)
    n = len(order)
    improved = True
    while improved:
        improved = False
        for i in range(n - 1):
            for j in range(i + 2, n):
                if j == n - 1 and i == 0:
                    continue
                candidate = (
                    order[: i + 1]
                    + list(reversed(order[i + 1 : j + 1]))
                    + order[j + 1 :]
                )
                if _path_has_crossings(candidate, points):
                    continue
                if _path_length(candidate, points) + 1e-12 < _path_length(
                    order, points
                ):
                    order = candidate
                    improved = True
    return order


def order_points_noncrossing_path_heuristic(points: np.ndarray) -> list[int]:
    """Fast multi-start nearest-neighbour tour, uncrossed (non-optimal length)."""
    pts = np.asarray(points, dtype=np.float64).reshape(-1, 2)
    n = int(pts.shape[0])
    if n <= 1:
        return list(range(n))
    if n == 2:
        return [0, 1]
    dist = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=2)
    return _nearest_neighbor_open_path(dist, pts)


def order_points_noncrossing_path_exact(points: np.ndarray) -> list[int]:
    """Exact minimum-length open Hamiltonian path among non-crossing permutations."""
    import itertools

    pts = np.asarray(points, dtype=np.float64).reshape(-1, 2)
    n = int(pts.shape[0])
    if n <= 1:
        return list(range(n))
    if n == 2:
        return [0, 1]

    dist = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=2)
    order: list[int] | None = None

    if n <= 9:
        best_dist = float("inf")
        for perm in itertools.permutations(range(n)):
            cand = list(perm)
            if _path_has_crossings(cand, pts):
                continue
            d = _path_length(cand, pts)
            if d < best_dist:
                best_dist = d
                order = cand

    if order is None and n <= 20:
        order = _held_karp_open_path(dist)

    if order is None:
        order = _nearest_neighbor_open_path(dist, pts)

    order = _uncross_path_order(order, pts)
    order = _improve_noncrossing_path(order, pts)
    return order


def order_points_shortest_noncrossing_path(
    points: np.ndarray, *, path_order: str = "heuristic"
) -> list[int]:
    """Permutation visiting each point once with no crossing segments."""
    if path_order == "exact":
        return order_points_noncrossing_path_exact(points)
    return order_points_noncrossing_path_heuristic(points)


def order_anchor_indices_for_direct_path(
    e: DashboardExperiment,
    anchor_indices: list[int],
    xcol: str,
    ycol: str,
    *,
    path_order: str = "heuristic",
) -> list[int]:
    """Reorder anchors for a direct traversal with no crossing scatter segments."""
    if path_order == "preserve" or len(anchor_indices) <= 2:
        return list(anchor_indices)
    coords = e.plot_df[[xcol, ycol]].values.astype(np.float64)
    points = np.vstack([coords[int(a)] for a in anchor_indices])
    order = order_points_shortest_noncrossing_path(points, path_order=path_order)
    return [anchor_indices[i] for i in order]


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
    # Path-order changes with no interpolation must remain a pure permutation of
    # the anchors (same particle set / point count), not a denser polyline.
    if int(interpolation_points) == 0:
        return z_anchor, rows, anchor_xy
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


def _nearest_plot_row_for_xy(coords: np.ndarray, pt: np.ndarray) -> int:
    """Dataset row nearest to ``pt`` in the current plot axis space."""
    return int(np.argmin(np.sum((coords - pt) ** 2, axis=1)))


def _direct_path_anchor_positions(
    n_anchors: int, interpolation_points: int, path_len: int
) -> dict[int, int]:
    """Map direct-interpolation path index -> anchor ordinal (0..n_anchors-1)."""
    if n_anchors < 2 or path_len < 2:
        return {}
    if interpolation_points <= 0:
        if path_len != n_anchors:
            return {}
        return {i: i for i in range(n_anchors)}
    per_seg = int(interpolation_points) + 1
    expected = (n_anchors - 1) * per_seg + 1
    if expected != path_len:
        return {}
    out: dict[int, int] = {0: 0}
    for seg_i in range(1, n_anchors):
        out[seg_i * per_seg] = seg_i
    return out


def _dedupe_consecutive_trajectory_points(
    rows: list[int], z_parts: list[np.ndarray], xy: np.ndarray
) -> tuple[list[int], np.ndarray, np.ndarray]:
    """Drop consecutive duplicate rows while keeping aligned z and xy."""
    if not rows:
        zdim = int(z_parts[0].shape[0]) if z_parts else 0
        return (
            [],
            np.zeros((0, zdim), dtype=np.float64),
            np.zeros((0, 2), dtype=np.float64),
        )
    keep_rows = [rows[0]]
    keep_z = [np.asarray(z_parts[0], dtype=np.float64)]
    keep_xy = [np.asarray(xy[0], dtype=np.float64)]
    for i in range(1, len(rows)):
        if rows[i] == keep_rows[-1]:
            continue
        keep_rows.append(rows[i])
        keep_z.append(np.asarray(z_parts[i], dtype=np.float64))
        keep_xy.append(np.asarray(xy[i], dtype=np.float64))
    return keep_rows, np.stack(keep_z, axis=0), np.asarray(keep_xy, dtype=np.float64)


def _dedupe_consecutive_rows(
    rows: list[int], vol_ids: list[str], xy: np.ndarray
) -> tuple[list[int], list[str], np.ndarray]:
    """Drop consecutive duplicate plot rows while keeping aligned vol ids and xy."""
    if not rows:
        return [], [], np.zeros((0, 2), dtype=np.float64)
    keep_rows = [rows[0]]
    keep_vols = [vol_ids[0]]
    keep_xy = [xy[0]]
    for i in range(1, len(rows)):
        if rows[i] == keep_rows[-1]:
            continue
        keep_rows.append(rows[i])
        keep_vols.append(vol_ids[i])
        keep_xy.append(xy[i])
    return keep_rows, keep_vols, np.asarray(keep_xy, dtype=np.float64)


def _compute_direct_anchor_trajectory_snap_volume_markers(
    e: DashboardExperiment,
    anchor_indices: list[int],
    xcol: str,
    ycol: str,
    interpolation_points: int,
    allowed_vol_ids: set[str] | frozenset[str] | None = None,
) -> tuple[np.ndarray, list[int], np.ndarray, list[str]]:
    """Direct anchor interpolation with each sample snapped to analyze volume markers."""
    from cryodrgn.dashboard.volume_slice_viewer import snap_xy_to_analyze_volume_markers

    _z_anchor, rows, anchor_xy = _compute_trajectory_from_anchor_indices(
        e, anchor_indices, xcol, ycol
    )
    if len(rows) < 2:
        raise ValueError("Need at least two anchor indices")
    if interpolation_points < 0:
        raise ValueError("Interpolation points must be >= 0")

    if interpolation_points == 0:
        snap_rows, snap_xy, snap_vols = snap_xy_to_analyze_volume_markers(
            e, anchor_xy, xcol, ycol, allowed_vol_ids=allowed_vol_ids
        )
        deduped_rows, deduped_vols, deduped_xy = _dedupe_consecutive_rows(
            snap_rows, snap_vols, snap_xy
        )
        zs = np.asarray(deduped_rows, dtype=int)
        return e.z[zs], deduped_rows, deduped_xy, deduped_vols

    n_points = int(interpolation_points) + 2
    full_rows: list[int] = []
    full_vols: list[str] = []
    full_xy: list[np.ndarray] = []
    for i in range(len(rows) - 1):
        xy_start = anchor_xy[i]
        xy_end = anchor_xy[i + 1]
        t = np.linspace(0.0, 1.0, n_points, dtype=np.float64)
        seg_xy = (np.outer(1.0 - t, xy_start) + np.outer(t, xy_end))[:-1]
        seg_rows, seg_snap_xy, seg_vols = snap_xy_to_analyze_volume_markers(
            e, seg_xy, xcol, ycol, allowed_vol_ids=allowed_vol_ids
        )
        full_rows.extend(seg_rows)
        full_vols.extend(seg_vols)
        full_xy.extend(seg_snap_xy)
    last_rows, last_xy, last_vols = snap_xy_to_analyze_volume_markers(
        e, anchor_xy[-1:], xcol, ycol, allowed_vol_ids=allowed_vol_ids
    )
    full_rows.extend(last_rows)
    full_vols.extend(last_vols)
    full_xy.extend(last_xy)
    full_xy_arr = np.asarray(full_xy, dtype=np.float64)
    if allowed_vol_ids is not None and interpolation_points > 0:
        deduped_rows, deduped_vols, deduped_xy = full_rows, full_vols, full_xy_arr
    else:
        deduped_rows, deduped_vols, deduped_xy = _dedupe_consecutive_rows(
            full_rows, full_vols, full_xy_arr
        )
    zs = np.asarray(deduped_rows, dtype=int)
    return e.z[zs], deduped_rows, deduped_xy, deduped_vols


def _compute_direct_anchor_trajectory_snap_nearest_particles(
    e: DashboardExperiment,
    anchor_indices: list[int],
    xcol: str,
    ycol: str,
    interpolation_points: int,
) -> tuple[np.ndarray, list[int], np.ndarray]:
    """Snap midpoint-split anchor trajectories to the nearest particles.

    The frontend inserts ``interpolation_points`` interior points between each
    consecutive anchor pair. Each inserted point is snapped to the nearest
    particle to the midpoint between its two adjacent trajectory points in
    *currently-plotted* coordinates.
    """
    if interpolation_points < 0:
        raise ValueError("Interpolation points must be >= 0")

    _z_anchor, anchor_rows, anchor_xy = _compute_trajectory_from_anchor_indices(
        e, anchor_indices, xcol, ycol
    )
    if len(anchor_rows) < 2:
        raise ValueError("Need at least two anchor indices")

    coords = e.plot_df[[xcol, ycol]].values.astype(np.float64)

    from collections import deque

    class _MidNode:
        __slots__ = ("row", "xy", "prev", "next")

        def __init__(self, row: int, xy: np.ndarray):
            self.row = int(row)
            self.xy = xy
            self.prev: "_MidNode | None" = None
            self.next: "_MidNode | None" = None

    interior = int(interpolation_points)

    full_rows: list[int] = []
    full_z_parts: list[np.ndarray] = []
    full_xy_parts: list[np.ndarray] = []

    seg_count = len(anchor_rows) - 1
    for seg_i in range(seg_count):
        left_row = int(anchor_rows[seg_i])
        right_row = int(anchor_rows[seg_i + 1])
        left_xy = np.asarray(anchor_xy[seg_i], dtype=np.float64)
        right_xy = np.asarray(anchor_xy[seg_i + 1], dtype=np.float64)

        left_node = _MidNode(left_row, left_xy)
        right_node = _MidNode(right_row, right_xy)
        left_node.next = right_node
        right_node.prev = left_node

        q = deque([(left_node, right_node)])
        inserted = 0

        while inserted < interior and q:
            l, r = q.popleft()
            # Endpoints must still be consecutive for this queued segment.
            if l.next is not r:
                continue

            mid_xy = (l.xy + r.xy) / 2.0
            mid_row = _nearest_plot_row_for_xy(coords, mid_xy)
            mid_xy_snapped = np.asarray(coords[mid_row], dtype=np.float64)

            m_node = _MidNode(mid_row, mid_xy_snapped)
            # Insert between l and r.
            m_node.prev = l
            m_node.next = r
            l.next = m_node
            r.prev = m_node

            inserted += 1
            q.append((l, m_node))
            q.append((m_node, r))

        # Collect segment points from left -> right.
        seg_rows: list[int] = []
        seg_xy: list[np.ndarray] = []
        node = left_node
        while node is not None:
            seg_rows.append(node.row)
            seg_xy.append(np.asarray(node.xy, dtype=np.float64))
            if node is right_node:
                break
            node = node.next

        if len(seg_rows) < 2:
            raise ValueError("Could not construct snapped trajectory segment.")

        # Avoid duplicating the right endpoint across adjacent segments.
        if seg_i < seg_count - 1:
            seg_rows = seg_rows[:-1]
            seg_xy = seg_xy[:-1]

        full_rows.extend(seg_rows)
        full_xy_parts.extend(seg_xy)
        full_z_parts.extend(
            [np.asarray(e.z[row], dtype=np.float64) for row in seg_rows]
        )

    # Preserve consecutive duplicate snapped particle rows.
    # This keeps trajectory point indices aligned with the volume-slider slots.
    full_xy_arr = np.asarray(full_xy_parts, dtype=np.float64)
    full_z_arr = np.stack(
        [np.asarray(z, dtype=np.float64) for z in full_z_parts], axis=0
    )
    return full_z_arr, full_rows, full_xy_arr


def _graph_neighbor_max_dist(ndist: np.ndarray, n: int, avg_neighbors: int) -> float:
    """Distance cutoff matching ``cryodrgn graph_traversal`` (top-k order statistic).

    Uses ``np.partition`` instead of a full sort over all query distances.
    """
    total_neighbors = max(1, min(int(n * avg_neighbors), int(ndist.size)))
    flat = ndist.reshape(-1)
    k_idx = total_neighbors - 1
    if k_idx <= 0:
        return float(flat[0])
    # Use the array returned by ``partition`` (in-place view can leave wrong value at k_idx).
    return float(np.partition(flat.copy(), k_idx)[k_idx])


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
        return cached[0], cached[1]

    from scipy.spatial import cKDTree

    z = np.asarray(e.z, dtype=np.float64)
    n = int(z.shape[0])
    if n < 2:
        raise ValueError("Need at least two latent points for graph traversal.")
    k_eff = min(k + 1, n)
    tree = cKDTree(z)
    q_dist, q_neighbors = tree.query(z, k=k_eff, workers=-1)
    neighbors = np.atleast_2d(np.asarray(q_neighbors, dtype=np.int64))
    ndist = np.atleast_2d(np.asarray(q_dist, dtype=np.float64))

    # Drop self-neighbor (distance 0); keep up to ``k`` true neighbors.
    neighbors = np.asarray(neighbors[:, 1:], dtype=np.int64)
    ndist = np.asarray(ndist[:, 1:], dtype=np.float64)
    if neighbors.shape[1] == 0:
        raise ValueError("Could not build nearest-neighbor graph from latent points.")

    max_dist = _graph_neighbor_max_dist(ndist, n, avg)
    keep = ndist <= max_dist
    neighbors = np.where(keep, neighbors, -1)
    ndist = np.where(keep, ndist, np.inf)

    csr = _neighbors_to_csr(neighbors, ndist)
    _TRAJ_GRAPH_NEIGHBOR_CACHE[key] = (neighbors, ndist, csr)

    return neighbors, ndist


def _neighbors_to_csr(neighbors: np.ndarray, dists: np.ndarray) -> "_csr_matrix":
    """Directed CSR graph from per-row neighbour lists (``-1``/``inf`` = no edge)."""
    from scipy.sparse import coo_matrix

    n, k = neighbors.shape
    row_rep = np.broadcast_to(np.arange(n, dtype=np.int64)[:, None], (n, k))
    valid = (neighbors >= 0) & np.isfinite(dists)
    rows = row_rep[valid]
    cols = neighbors[valid].astype(np.int64, copy=False)
    data = dists[valid].astype(np.float64, copy=False)
    return coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()


def _graph_csr_for_neighbors(neighbors: np.ndarray, dists: np.ndarray) -> "_csr_matrix":
    """Return a cached CSR graph when possible, else build one."""
    for _neighbors, _dists, csr in _TRAJ_GRAPH_NEIGHBOR_CACHE.values():
        if _neighbors is neighbors and _dists is dists:
            return csr
    return _neighbors_to_csr(neighbors, dists)


def _path_from_csgraph_predecessors(
    predecessors: np.ndarray, src: int, dest: int, dist_row: np.ndarray
) -> list[int] | None:
    if src == dest:
        return [int(src)]
    if not np.isfinite(float(dist_row[dest])):
        return None
    path = [int(dest)]
    cur = int(dest)
    while cur != int(src):
        cur = int(predecessors[cur])
        if cur < 0 or cur == _CSGRAPH_MISSING_PREDECESSOR:
            return None
        path.append(cur)
    path.reverse()
    return path


def _dijkstra_path_from_neighbors(
    neighbors: np.ndarray, dists: np.ndarray, src: int, dest: int
) -> list[int] | None:
    """Shortest path in a sparse directed neighbor graph (``-1``/``inf`` = no edge)."""
    from scipy.sparse.csgraph import dijkstra

    if src == dest:
        return [int(src)]
    csr = _graph_csr_for_neighbors(neighbors, dists)
    dist_matrix, predecessors = dijkstra(
        csr,
        directed=True,
        indices=int(src),
        return_predecessors=True,
    )
    dist_row = np.asarray(dist_matrix, dtype=np.float64).reshape(-1)
    pred_row = np.asarray(predecessors, dtype=np.int64).reshape(-1)
    return _path_from_csgraph_predecessors(pred_row, int(src), int(dest), dist_row)


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


def _parse_optional_traj_xy_custom(data: dict) -> list[list[float]] | None:
    """Parse optional client-supplied traj_xy polyline from a request body."""
    raw_traj_xy = data.get("traj_xy")
    if not isinstance(raw_traj_xy, list) or len(raw_traj_xy) < 2:
        return None
    pts: list[list[float]] = []
    for i, p in enumerate(raw_traj_xy[:200]):
        if not isinstance(p, list) or len(p) != 2:
            raise ValueError(f"traj_xy[{i}] must be [x, y].")
        try:
            pts.append([float(p[0]), float(p[1])])
        except (TypeError, ValueError) as err:
            raise ValueError(f"traj_xy[{i}] has non-numeric coordinates.") from err
    return pts


def parse_trajectory_request_body(e: DashboardExperiment, data: dict) -> dict:
    """Validate a trajectory POST JSON body.

    Shared by coords-only and volume decode.
    """
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
        anchor_path_order = parse_anchor_path_order(data)
        # Backward compatibility: older clients could request volume-marker
        # snapping via `snap_to_volume_points`; treat that as snap-to-particles.
        snap_to_nearest_particles = bool(
            data.get("snap_to_nearest_particles") or data.get("snap_to_volume_points")
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
            "anchor_path_order": anchor_path_order,
            "snap_to_nearest_particles": snap_to_nearest_particles,
            # Direct-line densify may send the PC/catalog polyline for display.
            "traj_xy_custom": _parse_optional_traj_xy_custom(data),
        }

    mode = str(data.get("mode", "direct")).strip().lower()
    if mode not in ("direct", "nearest"):
        raise ValueError('mode must be "direct" or "nearest".')

    traj_xy_custom = _parse_optional_traj_xy_custom(data)
    sx0 = sy0 = ex0 = ey0 = None
    if traj_xy_custom is not None:
        sx0, sy0 = traj_xy_custom[0]
        ex0, ey0 = traj_xy_custom[-1]
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
            "Direct traversal is only available for principal-component "
            "or z latent axes."
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
    """Return ``(z_traj, traj_rows_or_None, traj_xy)``."""
    if p.get("use_anchors"):
        if p["mode"] == "direct":
            path_order = str(p.get("anchor_path_order", "preserve"))
            ordered = order_anchor_indices_for_direct_path(
                e,
                p["anchor_indices"],
                p["xcol"],
                p["ycol"],
                path_order=path_order,
            )
            p["anchor_indices"] = ordered
            if p.get("snap_to_nearest_particles"):
                (
                    z_traj,
                    traj_rows,
                    traj_xy,
                ) = _compute_direct_anchor_trajectory_snap_nearest_particles(
                    e,
                    ordered,
                    p["xcol"],
                    p["ycol"],
                    int(p["n_points"]),
                )
                return z_traj, traj_rows, traj_xy
            z_traj, traj_rows, traj_xy = _compute_direct_anchor_trajectory(
                e,
                ordered,
                p["xcol"],
                p["ycol"],
                int(p["n_points"]),
            )
            # Client direct-line densify sends the PC/catalog polyline as traj_xy.
            # Keep that geometry for display while still decoding z from anchors.
            custom_xy = p.get("traj_xy_custom")
            if custom_xy is not None:
                custom = np.asarray(custom_xy, dtype=np.float64)
                if custom.ndim == 2 and custom.shape[0] == int(z_traj.shape[0]):
                    return z_traj, None, _round_direct_mode_traj_xy(custom)
            return z_traj, traj_rows, traj_xy
        z_traj, traj_rows, traj_xy = _compute_graph_anchor_trajectory(
            e,
            p["anchor_indices"],
            p["xcol"],
            p["ycol"],
            max_neighbors=int(p.get("max_neighbors", p["n_points"])),
            avg_neighbors=int(p.get("avg_neighbors", p["n_points"])),
        )
        return z_traj, traj_rows, traj_xy

    mode = p["mode"]
    sx0, sy0 = p["sx0"], p["sy0"]
    ex0, ey0 = p["ex0"], p["ey0"]
    xcol, ycol = p["xcol"], p["ycol"]
    n_points = p["n_points"]
    coords = e.plot_df[[xcol, ycol]].values.astype(np.float64)

    custom_xy = p.get("traj_xy_custom")
    if custom_xy is not None:
        pts = np.asarray(custom_xy, dtype=np.float64)
        if mode == "nearest":
            traj_rows, traj_xy = _snap_points_to_nearest_particles(coords, pts)
            z_traj = e.z[np.asarray(traj_rows, dtype=int)]
            return z_traj, traj_rows, traj_xy
        traj_rows = [int(np.argmin(np.sum((coords - pt) ** 2, axis=1))) for pt in pts]
        z_traj = e.z[np.asarray(traj_rows, dtype=int)]
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
    traj_rows, traj_xy = _snap_points_to_nearest_particles(coords, line_xy)
    z_traj = e.z[np.asarray(traj_rows, dtype=int)]
    return z_traj, traj_rows, traj_xy


# ---------------------------------------------------------------------------
# JSON payload assembly
# ---------------------------------------------------------------------------


def _resolve_traj_ref_to_plot_row(exp: DashboardExperiment, ref: int) -> int | None:
    """Map a trajectory row reference to a ``plot_df`` iloc position."""
    ri = int(ref)
    mapped = plot_df_rows_for_dataset_indices(exp, np.asarray([ri], dtype=int))
    if mapped:
        return int(mapped[0])
    if 0 <= ri < len(exp.plot_df):
        return ri
    return None


def _marker_colors_for_plot_rows(
    e: DashboardExperiment,
    plot_rows: list[int],
    color_col: str,
    *,
    continuous_palette: str | None = None,
    discrete_label_colors: dict[str, str] | None = None,
) -> list[str]:
    df = e.plot_df
    sub = df.iloc[plot_rows]
    if _lower_color_series_is_discrete(df[color_col]):
        hex_colors, _fk = _scatter_discrete_marker_arrays(
            df, sub, color_col, discrete_label_colors
        )
        return hex_colors
    _cvals, cmin, cmax = _continuous_series_stats(df[color_col])
    plotly_cs = normalize_continuous_palette(continuous_palette)
    vals = pd.to_numeric(sub[color_col], errors="coerce").to_numpy(dtype=np.float64)
    return numeric_array_to_plotly_hex(vals, plotly_cs, vmin=cmin, vmax=cmax)


def trajectory_marker_colors_for_rows(
    e: DashboardExperiment,
    traj_rows: list[int],
    color_col: str | None,
    *,
    continuous_palette: str | None = None,
    discrete_label_colors: dict[str, str] | None = None,
) -> list[str | None] | None:
    """Per-trajectory-point hex colours matching the scatter covariate scale."""
    if not traj_rows or not color_col or color_col == "none":
        return None
    if color_col not in e.plot_df.columns:
        return None
    resolved: list[int | None] = [
        _resolve_traj_ref_to_plot_row(e, int(r)) for r in traj_rows
    ]
    if not any(r is not None for r in resolved):
        return None
    batch_rows: list[int] = []
    batch_pos: list[int] = []
    for i, pr in enumerate(resolved):
        if pr is not None:
            batch_rows.append(int(pr))
            batch_pos.append(i)
    batch_colors = _marker_colors_for_plot_rows(
        e,
        batch_rows,
        color_col,
        continuous_palette=continuous_palette,
        discrete_label_colors=discrete_label_colors,
    )
    out: list[str | None] = [None] * len(traj_rows)
    for j, pos in enumerate(batch_pos):
        out[pos] = batch_colors[j]
    return out


def trajectory_marker_colors_for_particle_indices(
    e: DashboardExperiment,
    particle_indices: list[int | None],
    color_col: str | None,
    *,
    continuous_palette: str | None = None,
    discrete_label_colors: dict[str, str] | None = None,
) -> list[str | None] | None:
    """Sparse per-point colours when only some trajectory samples are particles."""
    if not particle_indices or not color_col or color_col == "none":
        return None
    if color_col not in e.plot_df.columns:
        return None
    out: list[str | None] = [None] * len(particle_indices)
    batch_rows: list[int] = []
    batch_pos: list[int] = []
    for i, pidx in enumerate(particle_indices):
        if pidx is None:
            continue
        pr = _resolve_traj_ref_to_plot_row(e, int(pidx))
        if pr is None:
            continue
        batch_rows.append(int(pr))
        batch_pos.append(i)
    if not batch_rows:
        return out
    batch_colors = _marker_colors_for_plot_rows(
        e,
        batch_rows,
        color_col,
        continuous_palette=continuous_palette,
        discrete_label_colors=discrete_label_colors,
    )
    for j, pos in enumerate(batch_pos):
        out[pos] = batch_colors[j]
    return out


def attach_trajectory_marker_colors(
    e: DashboardExperiment,
    payload: dict,
    color_col: str | None,
    *,
    continuous_palette: str | None = None,
    discrete_label_colors: dict[str, str] | None = None,
) -> None:
    """Add ``traj_marker_colors`` to a trajectory coords/volumes JSON payload."""
    if payload.get("traj_marker_colors"):
        return
    traj_rows = payload.get("traj_rows")
    if traj_rows:
        marker_colors = trajectory_marker_colors_for_rows(
            e,
            traj_rows,
            color_col,
            continuous_palette=continuous_palette,
            discrete_label_colors=discrete_label_colors,
        )
        if marker_colors:
            payload["traj_marker_colors"] = marker_colors
        return
    pidx = payload.get("traj_particle_indices")
    if pidx:
        marker_colors = trajectory_marker_colors_for_particle_indices(
            e,
            pidx,
            color_col,
            continuous_palette=continuous_palette,
            discrete_label_colors=discrete_label_colors,
        )
        if marker_colors:
            payload["traj_marker_colors"] = marker_colors


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
    color_col: str | None = None,
    continuous_palette: str | None = None,
    discrete_label_colors: dict[str, str] | None = None,
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
    if traj_rows is not None and mode in ("nearest", "graph", "direct"):
        payload["traj_particle_indices"] = [
            _plot_row_particle_index(e, r) for r in traj_rows
        ]
    attach_trajectory_marker_colors(
        e,
        payload,
        color_col,
        continuous_palette=continuous_palette,
        discrete_label_colors=discrete_label_colors,
    )
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
    anchor_path_order: str = "preserve",
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
        "anchor_path_order": parse_anchor_path_order(
            {"anchor_path_order": anchor_path_order}
        ),
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
            anchor_indices=p["anchor_indices"],
            interpolation_points=n_points,
            n_total=int(np.asarray(z_traj).shape[0]),
        )
        if pidx is not None:
            payload["traj_particle_indices"] = pidx
    payload["anchor_indices"] = p["anchor_indices"]
    payload["anchor_path_order"] = p.get("anchor_path_order", "preserve")
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


def random_dataset_indices(
    e: DashboardExperiment,
    k: int = 10,
    exclude: list[int] | tuple[int, ...] | None = None,
) -> list[int]:
    """Up to ``k`` distinct random row indices into ``z``.

    When ``exclude`` is provided, those indices are omitted so callers can
    append unused particles to an existing trajectory.

    Sampling avoids materializing an ``O(n)`` candidate list for typical
    large stacks (rejection sampling). A dense fallback is used only when
    most indices are excluded.
    """
    n = int(e.z.shape[0])
    if n < 2:
        raise ValueError("Dataset has fewer than 2 particles")
    exclude_set = {int(i) for i in (exclude or []) if 0 <= int(i) < n}
    available = n - len(exclude_set)
    if available < 1:
        raise ValueError("No unused particle indices available to sample")
    take = min(max(0, int(k)), available)
    if take < 1:
        return []
    rng = np.random.default_rng()

    # Fast path: reject draws that hit exclude / duplicates. Expected work is
    # O(take) when exclude is small relative to n (the common dashboard case).
    if available > max(take * 20, 1024):
        chosen: set[int] = set()
        guard = 0
        max_guard = max(10_000, take * 200)
        while len(chosen) < take and guard < max_guard:
            batch = int(min(max(take * 4, 32), take * 8, max_guard - guard))
            draws = rng.integers(0, n, size=batch)
            guard += batch
            for raw in draws:
                idx = int(raw)
                if idx in exclude_set or idx in chosen:
                    continue
                chosen.add(idx)
                if len(chosen) >= take:
                    break
        if len(chosen) >= take:
            return list(chosen)

    # Dense / mostly-excluded: build the available index array once.
    mask = np.ones(n, dtype=bool)
    if exclude_set:
        mask[list(exclude_set)] = False
    candidates = np.flatnonzero(mask)
    if candidates.size < 1:
        raise ValueError("No unused particle indices available to sample")
    take = min(take, int(candidates.size))
    return rng.choice(candidates, size=take, replace=False).astype(int).tolist()


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
    the centered coordinates, with an axis-aligned fallback).
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
        # Ill-conditioned SVD; keep the axis-aligned default direction.
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
