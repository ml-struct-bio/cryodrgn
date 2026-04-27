"""Particle-thumbnail sampling, encoding, and montage helpers.

These feed the explorer's preview montage, the hover pre-load cache, and the
trajectory nearest-neighbor insets. Kept separate from :mod:`app` so each
helper is a pure function over a :class:`DashboardExperiment`.
"""

from __future__ import annotations

import base64
import io
import logging
import os
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np

from cryodrgn.dashboard.data import DashboardExperiment, particle_image_array
from cryodrgn.dashboard.mpl_style import ezlab_matplotlib_rc

logger = logging.getLogger(__name__)


def encode_particle_batch(
    mrcfile: str,
    datadir: str | None,
    global_indices: Iterable[int],
    max_px: int,
) -> list[str]:
    """Load and encode particles as base64 JPEGs for a ``ProcessPoolExecutor`` worker.

    All imports are local so the function is self-contained after fork.
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
    """Small greyscale JPEG (base64) for trajectory NN inset previews."""
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


def _stratified_xy_row_indices(
    coords: np.ndarray,
    rng: np.random.Generator,
    total_k: int,
) -> set[int]:
    """Mix: random + far-from-reference-centroids + large NN distance (``total_k`` picks)."""
    from scipy.spatial import cKDTree
    from scipy.spatial.distance import cdist as _cdist

    m = int(coords.shape[0])
    if m == 0:
        return set()
    total_k = min(m, total_k)
    k_random = total_k * 2 // 3
    k_mean = total_k // 6
    random_rows = set(rng.choice(m, size=min(k_random, m), replace=False).tolist())

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


def sample_plot_df_rows_for_preload(
    exp: DashboardExperiment,
    xcol: str,
    ycol: str,
    *,
    restrict_to_rows: list[int] | None = None,
) -> tuple[list[int], list[int]]:
    """Pick up to 5k plot_df rows for thumbnail preload (random + spread-out in XY).

    If ``restrict_to_rows`` is set, only those plot_df indices are eligible.
    Returns ``(plot_df_rows, dataset_indices)`` sorted by dataset index.
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

    local_sel = _stratified_xy_row_indices(coords, rng, min(len(coords), 5000))
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
