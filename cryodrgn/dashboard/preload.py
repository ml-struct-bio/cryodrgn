"""Particle-thumbnail sampling, encoding, and montage helpers.

These feed the explorer's preview montage, the hover pre-load cache, and the
trajectory nearest-neighbor insets. Kept separate from :mod:`app` so each
helper is a pure function over a :class:`DashboardExperiment`.
"""

from __future__ import annotations

import base64
import io
import logging
import math
import os
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np

from cryodrgn.dashboard.data import DashboardExperiment, particle_image_array
from cryodrgn.dashboard.plots_figure_utils import ezlab_matplotlib_rc

logger = logging.getLogger(__name__)

DEFAULT_PRELOAD_IMAGE_LIMIT = 1000

# Max thumbnails encoded in one ``api_preload_images`` response (JSON + base64 size).
MAX_PRELOAD_IMAGES_PER_HTTP_RESPONSE = 400


def explorer_cache_size_power10_step(scatter_plotted_point_count: int) -> int:
    """Largest power of ten not greater than 5% of scatter-plotted point count.

    Used as the HTML ``step`` for the particle-explorer cache size field and as
    the scale for the default initial cache size (clamped to the plotted cap).
    """
    n = max(0, int(scatter_plotted_point_count))
    if n == 0:
        return 1
    x = 0.05 * n
    if x < 1:
        return 1
    k = int(math.floor(math.log10(x)))
    step = 10**k
    return max(1, int(step))


def explorer_initial_preload_image_limit(scatter_plotted_point_count: int) -> int:
    """Default first-build cache size.

    :func:`explorer_cache_size_power10_step` capped by plotted *n*.
    """
    n = max(0, int(scatter_plotted_point_count))
    if n == 0:
        return 1
    step = explorer_cache_size_power10_step(n)
    return min(n, step)


def encode_particle_batch(
    mrcfile: str,
    datadir: str | None,
    global_indices: Iterable[int],
    max_px: int,
) -> list[str]:
    """Load and encode raw particles as base64 JPEGs for thumbnail preloading.

    Use :class:`ImageSource` directly instead of :class:`ImageDataset`; dashboard
    thumbnails do their own percentile scaling and do not need ImageDataset's
    costly normalization estimates.
    """
    import base64 as _b64
    import io as _io

    import numpy as _np
    from PIL import Image as PILImage

    from cryodrgn.source import ImageSource

    src = ImageSource.from_file(mrcfile, lazy=True, datadir=datadir or "")
    out: list[str] = []
    for gidx in global_indices:
        raw = src.images(gidx, as_numpy=True)
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


def montage_bytes(exp: DashboardExperiment, row_indices: list[int]) -> bytes:
    """Up-to-5x5 PNG montage of particles for ``row_indices`` (empty = hint image)."""
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


def particle_thumbnail_b64_from_row(
    exp: DashboardExperiment, row_index: int, max_side: int = 72
) -> str | None:
    """Small grayscale JPEG (base64) for trajectory NN inset previews."""
    if not exp.can_preview_particles:
        return None
    try:
        from PIL import Image

        img = particle_image_array(exp, int(row_index))
        lo, hi = np.percentile(img, (2, 98))
        u8 = (np.clip((img - lo) / (hi - lo + 1e-9), 0, 1) * 255).astype(np.uint8)
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


# k-th NN distance for outlier scoring; grid sizing uses
# sqrt(n_outlier * factor) bins per axis.
_PRELOAD_KNN_OUTLIER_K = 10
_PRELOAD_SPACING_BIN_FACTOR = 2.0


def _kth_nn_distances(coords: np.ndarray, k: int) -> np.ndarray:
    """Per-point distance to the *k*-th nearest neighbor in 2-D.

    Requires ``k < n_points``.
    """
    from scipy.spatial import cKDTree

    m = int(coords.shape[0])
    d = np.zeros(m, dtype=np.float64)
    if m <= 1 or k < 1:
        return d
    kk = min(int(k), m - 1)
    k_query = kk + 1
    tree = cKDTree(coords)
    dists, _ = tree.query(coords, k=k_query)
    if dists.ndim == 1:
        dists = dists.reshape(m, k_query)
    return dists[:, -1].astype(np.float64, copy=False)


def _grid_spaced_outlier_pick(
    coords: np.ndarray,
    order_by_score_desc: np.ndarray,
    want: int,
    exclude: set[int],
    *,
    n_bins: int,
) -> list[int]:
    """Take up to ``want`` local indices.

    Prefer high scores, at most one per grid cell.
    """
    m = int(coords.shape[0])
    if want <= 0 or m == 0:
        return []
    x0, y0 = coords.min(axis=0)
    x1, y1 = coords.max(axis=0)
    rx = float(x1 - x0) + 1e-12
    ry = float(y1 - y0) + 1e-12
    nx = (coords[:, 0] - x0) / rx
    ny = (coords[:, 1] - y0) / ry

    picked: list[int] = []
    used_cells: set[tuple[int, int]] = set()
    nb = max(2, int(n_bins))

    for li in order_by_score_desc:
        if len(picked) >= want:
            break
        i = int(li)
        if i in exclude:
            continue
        cx = min(int(nx[i] * nb), nb - 1)
        cy = min(int(ny[i] * nb), nb - 1)
        cell = (cx, cy)
        if cell in used_cells:
            continue
        used_cells.add(cell)
        picked.append(i)

    if len(picked) < want:
        for li in order_by_score_desc:
            if len(picked) >= want:
                break
            i = int(li)
            if i in exclude or i in picked:
                continue
            picked.append(i)
    return picked[:want]


def _hybrid_random_knn_spaced_local_indices(
    coords: np.ndarray,
    rng: np.random.Generator,
    total_k: int,
    *,
    k_nn: int = _PRELOAD_KNN_OUTLIER_K,
    spacing_bin_factor: float = _PRELOAD_SPACING_BIN_FACTOR,
) -> set[int]:
    """``total_k`` local indices: half uniform random, half high k-NN distance.

    XY grid spacing applied to the k-NN half.
    """
    m = int(coords.shape[0])
    if m == 0:
        return set()
    total_k = min(m, int(total_k))

    n_rand = total_k // 2
    n_out = total_k - n_rand

    rand_local = rng.choice(m, size=n_rand, replace=False)
    rand_set = {int(x) for x in rand_local}

    kk = min(max(1, int(k_nn)), m - 1)
    scores = _kth_nn_distances(coords, kk)
    order = np.argsort(-scores)
    n_bins = max(2, int(np.ceil(np.sqrt(max(n_out, 1) * float(spacing_bin_factor)))))
    outlier_list = _grid_spaced_outlier_pick(
        coords,
        order,
        n_out,
        rand_set,
        n_bins=n_bins,
    )
    out_set = {int(x) for x in outlier_list}

    return rand_set | out_set


def sample_plot_df_rows_for_preload(
    exp: DashboardExperiment,
    xcol: str,
    ycol: str,
    *,
    restrict_to_rows: list[int] | None = None,
    exclude_rows: list[int] | None = None,
    max_images: int = DEFAULT_PRELOAD_IMAGE_LIMIT,
) -> tuple[list[int], list[int]]:
    """Pick plot_df rows for thumbnail preload.

    Uses a hybrid rule on the current scatter axes: half of the picks are uniform
    random; the other half favor large k-th nearest-neighbor distances in the
    2-D embedding, with a coarse XY grid so high-scoring outliers do not stack
    in the same visual bin.

    If ``restrict_to_rows`` is set, only those plot_df indices are eligible.
    Returns ``(plot_df_rows, dataset_indices)`` sorted by dataset index.
    """
    max_images = max(1, int(max_images))
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

    if exclude_rows:
        excluded = {int(r) for r in exclude_rows}
        keep = [i for i, r in enumerate(inv_map) if r not in excluded]
        if not keep:
            return [], []
        coords = coords[np.array(keep, dtype=int)]
        inv_map = [inv_map[i] for i in keep]

    local_sel = _hybrid_random_knn_spaced_local_indices(
        coords, rng, min(len(coords), max_images)
    )
    plot_rows_set = {inv_map[i] for i in local_sel}
    pairs = sorted(
        ((r, int(exp.all_indices[r])) for r in plot_rows_set),
        key=lambda p: p[1],
    )
    return [p[0] for p in pairs], [p[1] for p in pairs]


def _preload_cache_time_estimate_bounds(cpus: int) -> tuple[int, int]:
    """Rough wall-clock range (seconds) for caching up to 5k 96px JPEG thumbs.

    Scales inversely with ``cpus`` (preload worker count), anchored at ~20–30s
    with 4 cores (empirical). I/O and fixed overhead keep a modest floor at
    high core counts.
    """
    cpus = max(1, int(cpus))
    low = max(8, int(round(20 * 4 / cpus)))
    high = max(max(18, low + 6), int(round(30 * 4 / cpus)))
    return low, high


def format_preload_cache_time_hint(cpus: int) -> str:
    """Human-readable wall-time hint for thumbnail preload (e.g. UI copy)."""
    lo, hi = _preload_cache_time_estimate_bounds(cpus)
    cpus = max(1, int(cpus))
    span = f"~{lo}s" if lo == hi else f"~{lo}\u2013{hi}s"
    cword = "core" if cpus == 1 else "cores"
    return f"{span} for {cpus} {cword}"


def load_plot_df_rows_from_plot_inds_file(
    exp: DashboardExperiment, plot_inds_path: str | None
) -> list[int]:
    """Map ``--plot-inds`` pickle dataset indices to ``plot_df`` rows.

    Returns ``[]`` when the path is empty, missing, or errors on load.
    """
    from cryodrgn import utils

    from cryodrgn.dashboard.trajectory import plot_df_rows_for_dataset_indices

    path = (plot_inds_path or "").strip()
    if not path or not os.path.isfile(path):
        return []
    try:
        pre = utils.load_pkl(path)
        return plot_df_rows_for_dataset_indices(exp, np.asarray(pre))
    except Exception as err:
        logger.warning("Could not load plot-inds from %s: %s", path, err)
        return []
