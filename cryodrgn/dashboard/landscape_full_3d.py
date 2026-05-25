"""Dashboard support for the ``analyze_landscape_full`` 3D volume-landscapes view.

Loads ``landscape.{epoch}/landscape_full`` pickles (``z.sampled.pkl``,
``ind.sampled.pkl``, ``vol_pca_all.pkl``, …), builds the sampled ``plot_df``
used for colouring, and produces the Plotly JSON for the scatter (axes
restricted to ``landscape_vol_PC*`` from ``vol_pca_all.pkl``).
"""

from __future__ import annotations

import os
import re
from typing import Any

import numpy as np
import pandas as pd

from cryodrgn import utils
from cryodrgn.dashboard.column_names import (
    VOL_LANDSCAPE_3D_PLOT_DF_ROW,
    VOL_LANDSCAPE_IS_SKETCH_CENTROID,
    VOL_LANDSCAPE_NEAREST_SKETCH_VOL,
)
from cryodrgn.dashboard.covariate_labels import covariate_display_map
from cryodrgn.dashboard.data import DashboardExperiment

_LANDSCAPE_FULL_SUBDIR = "landscape_full"
_LANDSCAPE_VOL_PC_COL_RE = re.compile(r"^landscape_vol_PC(\d+)$")
_SKETCH_CENTERS_IND_CACHE: dict[tuple[str, int], np.ndarray] = {}
_LANDSCAPE_FULL_SAMPLED_DF_CACHE: dict[tuple, pd.DataFrame] = {}


def _landscape_full_sampled_df_cache_key(exp: DashboardExperiment) -> tuple:
    """L3-E1 cache key from output mtimes + table size."""
    outdir = landscape_full_outdir(exp.workdir, exp.epoch)
    ind_path = os.path.join(outdir, "ind.sampled.pkl")
    z_path = os.path.join(outdir, "z.sampled.pkl")
    vol_path = os.path.join(outdir, "vol_pca_all.pkl")

    def _mt(p: str) -> float:
        return os.path.getmtime(p) if os.path.isfile(p) else 0.0

    return (
        os.path.abspath(exp.workdir),
        int(exp.epoch),
        _mt(ind_path),
        _mt(z_path),
        _mt(vol_path),
        len(exp.plot_df),
    )


def _sketch_centers_ind_for_epoch(workdir: str, epoch: int) -> np.ndarray | None:
    """Centroid particle row-indices (`centers_ind.txt`) for this epoch's sketch."""
    key = (workdir, int(epoch))
    if key in _SKETCH_CENTERS_IND_CACHE:
        return _SKETCH_CENTERS_IND_CACHE[key]
    try:
        from cryodrgn.dashboard.landscape_volpca import (
            load_sketch_centroid_plot_df_rows,
            load_vol_pca_matrix,
            resolve_kmeans_sketch_bundle,
            landscape_dir_for_epoch,
        )

        landscape_dir = landscape_dir_for_epoch(workdir, int(epoch))
        k_sketch, kmeans_dir = resolve_kmeans_sketch_bundle(landscape_dir)
        centers = load_vol_pca_matrix(landscape_dir, k_sketch)
        n_expected = int(centers.shape[0])
        centers_ind = load_sketch_centroid_plot_df_rows(kmeans_dir, n_expected)
    except Exception:
        centers_ind = None
    if centers_ind is None:
        _SKETCH_CENTERS_IND_CACHE[key] = None  # type: ignore[assignment]
    else:
        _SKETCH_CENTERS_IND_CACHE[key] = centers_ind
    return centers_ind


LANDSCAPE_FULL_3D_LEAD_HTML = (
    '<p class="cryo-dash-lead">Double-click a point to add or remove the nearest '
    "k-means sketch volume to the animation, or choose volumes to animate at random."
    "<br /><br />"
    "Use toggle or histogram selection to filter volumes based on color covariates.</p>"
)

LANDSCAPE_FULL_3D_LEGEND_CONTEXT_EXTRA: dict[str, str] = {
    "scope": "landscape_full_sampled"
}


def landscape_full_outdir(workdir: str, epoch: int) -> str:
    """``analyze_landscape_full`` output dir (``landscape.N/landscape_full``)."""
    return os.path.join(workdir, f"landscape.{int(epoch)}", _LANDSCAPE_FULL_SUBDIR)


def landscape_full_ready(workdir: str, epoch: int) -> bool:
    """True when ``z.sampled.pkl`` and ``ind.sampled.pkl`` exist for this epoch."""
    d = landscape_full_outdir(workdir, int(epoch))
    return os.path.isfile(os.path.join(d, "z.sampled.pkl")) and os.path.isfile(
        os.path.join(d, "ind.sampled.pkl")
    )


def landscape_full_vol_pca_all_path(workdir: str, epoch: int) -> str:
    """Path to ``vol_pca_all.pkl`` under ``landscape.{epoch}/landscape_full``."""
    return os.path.join(landscape_full_outdir(workdir, int(epoch)), "vol_pca_all.pkl")


def landscape_full_vol_pca_dim(workdir: str, epoch: int) -> int | None:
    """Number of volume PCA columns in ``vol_pca_all.pkl`` (``None`` if missing)."""
    p = landscape_full_vol_pca_all_path(workdir, int(epoch))
    if not os.path.isfile(p):
        return None
    try:
        v = np.asarray(utils.load_pkl(p), dtype=np.float64)
    except Exception:
        return None
    if v.ndim != 2 or v.shape[1] < 1:
        return None
    return int(v.shape[1])


def landscape_full_3d_ready(workdir: str, epoch: int) -> bool:
    """True when sampled latents exist and ``vol_pca_all.pkl`` has >= 3 columns."""
    if not landscape_full_ready(workdir, int(epoch)):
        return False
    n_pc = landscape_full_vol_pca_dim(workdir, int(epoch))
    return n_pc is not None and n_pc >= 3


def landscape_full_vol_pca_axis_columns(df: pd.DataFrame) -> list[str]:
    """Sorted ``landscape_vol_PC*`` columns in ``df`` (from ``vol_pca_all.pkl``)."""
    found: list[tuple[int, str]] = []
    for name in df.columns:
        m = _LANDSCAPE_VOL_PC_COL_RE.fullmatch(str(name))
        if m:
            found.append((int(m.group(1)), str(name)))
    found.sort(key=lambda t: t[0])
    return [t[1] for t in found]


def landscape_full_sampled_numeric_covariates(df: pd.DataFrame) -> list[str]:
    """Numeric columns for the colour selector (excludes helpers)."""
    skip = frozenset(
        {
            "index",
            VOL_LANDSCAPE_3D_PLOT_DF_ROW,
            VOL_LANDSCAPE_NEAREST_SKETCH_VOL,
            VOL_LANDSCAPE_IS_SKETCH_CENTROID,
        }
    )
    return [c for c in df.select_dtypes(include=[np.number]).columns if c not in skip]


def attach_landscape_nearest_sketch_vol_column(
    exp: DashboardExperiment, base: pd.DataFrame
) -> pd.DataFrame:
    """Add ``VOL_LANDSCAPE_NEAREST_SKETCH_VOL`` when sketch volumes exist (for GIFs)."""
    if VOL_LANDSCAPE_NEAREST_SKETCH_VOL in base.columns:
        return base
    from cryodrgn.dashboard.landscape_volpca import (
        kmeans_sorted_vol_indices,
        landscape_analysis_ready,
        landscape_dir_for_epoch,
        load_sketch_centroid_plot_df_rows,
        load_vol_pca_matrix,
        resolve_kmeans_sketch_bundle,
    )

    if not landscape_analysis_ready(exp.workdir, exp.epoch):
        return base
    landscape_dir = landscape_dir_for_epoch(exp.workdir, exp.epoch)
    try:
        k_sketch, kmeans_dir = resolve_kmeans_sketch_bundle(landscape_dir)
    except FileNotFoundError:
        return base
    vol_pc_cols = landscape_full_vol_pca_axis_columns(base)
    if not vol_pc_cols:
        return base
    centers = load_vol_pca_matrix(landscape_dir, k_sketch)
    sorted_vols = kmeans_sorted_vol_indices(kmeans_dir)
    n = int(centers.shape[0])
    vol_ids = sorted_vols[:n] if len(sorted_vols) >= n else sorted_vols
    if len(vol_ids) < n:
        return base
    d_plot = len(vol_pc_cols)
    d_ctr = int(centers.shape[1])
    d = min(d_plot, d_ctr)
    if d < 1:
        return base
    use_cols = vol_pc_cols[:d]
    x = base[use_cols].to_numpy(dtype=np.float64)
    c = centers[:, :d].astype(np.float64)
    xx = np.sum(x * x, axis=1, keepdims=True)
    cc = np.sum(c * c, axis=1, keepdims=True).T
    dist2 = xx + cc - 2.0 * (x @ c.T)
    nearest_rows = np.argmin(dist2, axis=1).astype(np.int64)
    nearest_vol = np.array([vol_ids[int(r)] for r in nearest_rows], dtype=np.int64)
    out = base.copy()
    out[VOL_LANDSCAPE_NEAREST_SKETCH_VOL] = nearest_vol
    is_cent = np.zeros(len(out), dtype=np.int64)
    try:
        centers_ind = load_sketch_centroid_plot_df_rows(kmeans_dir, n)
    except (FileNotFoundError, OSError, ValueError):
        pass
    else:
        plot_rows = out[VOL_LANDSCAPE_3D_PLOT_DF_ROW].to_numpy(dtype=np.int64)
        # Mark sketch-centroid particles directly by their plot-row indices.
        # This avoids relying on any particular ordering alignment between
        # `centers_ind.txt` and the k-means volume-id mapping used for
        # `nearest_vol`.
        is_cent = np.isin(plot_rows, centers_ind).astype(np.int64, copy=False)
    out[VOL_LANDSCAPE_IS_SKETCH_CENTROID] = is_cent
    return out


def _build_landscape_full_sampled_plot_df(exp: DashboardExperiment) -> pd.DataFrame:
    """Merge ``ind.sampled.pkl`` / ``z.sampled.pkl`` with the analysis table.

    Adds volume covariates from ``vol_pca_all.pkl`` to the sampled rows.
    """
    outdir = landscape_full_outdir(exp.workdir, exp.epoch)
    ind_path = os.path.join(outdir, "ind.sampled.pkl")
    z_path = os.path.join(outdir, "z.sampled.pkl")
    if not (os.path.isfile(ind_path) and os.path.isfile(z_path)):
        raise FileNotFoundError(
            f"Landscape full outputs missing under {outdir} (expected "
            "ind.sampled.pkl and z.sampled.pkl from analyze_landscape_full)."
        )

    ind = np.asarray(utils.load_pkl(ind_path), dtype=np.int64)
    z_s = np.asarray(utils.load_pkl(z_path), dtype=np.float64)
    if z_s.ndim != 2:
        raise ValueError("z.sampled.pkl must contain a 2-D array.")
    if len(ind) != len(z_s):
        raise ValueError("ind.sampled.pkl and z.sampled.pkl have different lengths.")
    zdim = int(exp.z.shape[1])
    if z_s.shape[1] != zdim:
        raise ValueError(f"z.sampled.pkl width {z_s.shape[1]} != run z dim {zdim}.")

    n_full = len(exp.plot_df)
    if ind.min() < 0 or ind.max() >= n_full:
        raise ValueError(
            "ind.sampled.pkl contains row indices outside the current analysis table."
        )

    # Ensure sketch-centroid particle rows are present.
    # `ind.sampled.pkl` is a performance sample; if it drops centroid rows, the
    # volume-landscape UI cannot show italics for those volumes.
    # We union the sketch centroid rows (`centers_ind.txt`) into the sampled indices.
    centers_ind = _sketch_centers_ind_for_epoch(exp.workdir, exp.epoch)
    if centers_ind is not None and len(centers_ind):
        ind_used = np.unique(
            np.concatenate([ind, centers_ind.astype(np.int64, copy=False)])
        )
        ind_used = ind_used[(ind_used >= 0) & (ind_used < n_full)]
    else:
        ind_used = ind

    base = exp.plot_df.iloc[ind_used].copy().reset_index(drop=True)
    zdim = int(exp.z.shape[1])
    z_used = np.asarray(exp.z[ind_used], dtype=np.float64)
    for i in range(zdim):
        base[f"z{i}"] = z_used[:, i]

    base[VOL_LANDSCAPE_3D_PLOT_DF_ROW] = ind_used.astype(np.int64, copy=False)

    n_particles = n_full
    vol_all_path = os.path.join(outdir, "vol_pca_all.pkl")
    if os.path.isfile(vol_all_path):
        v = np.asarray(utils.load_pkl(vol_all_path), dtype=np.float64)
        if v.ndim == 2 and v.shape[0] == n_particles:
            for j in range(v.shape[1]):
                base[f"landscape_vol_PC{j + 1}"] = v[ind_used, j]

    umap_path = os.path.join(outdir, "umap_vol_pca.pkl")
    if os.path.isfile(umap_path):
        u = np.asarray(utils.load_pkl(umap_path), dtype=np.float64)
        if u.ndim == 2 and u.shape[0] == n_particles and u.shape[1] >= 2:
            base["landscape_vol_UMAP1"] = u[ind_used, 0]
            base["landscape_vol_UMAP2"] = u[ind_used, 1]

    cl_path = os.path.join(outdir, "full_clustering", "cluster_labels.pkl")
    if os.path.isfile(cl_path):
        cl = np.asarray(utils.load_pkl(cl_path), dtype=np.int64)
        if cl.shape[0] == n_particles:
            base["landscape_vol_cluster"] = cl[ind_used]

    return attach_landscape_nearest_sketch_vol_column(exp, base)


def landscape_full_sampled_plot_df(exp: DashboardExperiment) -> pd.DataFrame:
    """L3-E1: process-level cache for sampled plot table."""
    key = _landscape_full_sampled_df_cache_key(exp)
    hit = _LANDSCAPE_FULL_SAMPLED_DF_CACHE.get(key)
    if hit is not None:
        return hit
    built = _build_landscape_full_sampled_plot_df(exp)
    _LANDSCAPE_FULL_SAMPLED_DF_CACHE[key] = built
    return built


def landscape_full_3d_scatter_plotly_json(
    exp: DashboardExperiment,
    *,
    xcol: str,
    ycol: str,
    zcol: str,
    color_col: str | None,
    continuous_palette: str | None = None,
    color_filter: dict[str, Any] | None = None,
    discrete_label_colors: dict[str, str] | None = None,
) -> str:
    """Plotly figure JSON for the volume-landscape 3D scatter (volume PCA axes only)."""
    from cryodrgn.dashboard.plots_scatter import scatter3d_z_json

    sampled = landscape_full_sampled_plot_df(exp)
    vol_axes = landscape_full_vol_pca_axis_columns(sampled)
    if len(vol_axes) < 3:
        raise ValueError(
            "vol_pca_all.pkl must yield at least three landscape_vol_PC* columns "
            "in the sampled table."
        )
    ax_allow = frozenset(vol_axes)
    return scatter3d_z_json(
        exp,
        xcol,
        ycol,
        zcol,
        color_col,
        continuous_palette=continuous_palette,
        color_filter=color_filter,
        discrete_label_colors=discrete_label_colors,
        plot_df=sampled,
        uirevision="scatter3d_z_landscape_full",
        xyz_axes_allowed=ax_allow,
        volume_landscape_3d_style=True,
    )


def landscape_full_3d_latent_3d_template_kwargs(
    exp: DashboardExperiment,
    *,
    scatter3d_url: str,
) -> dict[str, Any]:
    """Keyword arguments for ``latent_3d.html`` when the interface is fully available.

    Raises ``FileNotFoundError``, ``ValueError``, or ``OSError`` if the sampled
    table cannot be built. Raises ``ValueError`` if fewer than three volume
    PCA columns are present.
    """
    sampled = landscape_full_sampled_plot_df(exp)
    vol_axes = landscape_full_vol_pca_axis_columns(sampled)
    if len(vol_axes) < 3:
        raise ValueError(
            "Expected at least three columns in vol_pca_all.pkl (landscape_vol_PC1 …) "
            "for this epoch."
        )
    dx, dy, dz = vol_axes[0], vol_axes[1], vol_axes[2]
    cols = landscape_full_sampled_numeric_covariates(sampled)
    return {
        "page_title": "3D volume landscapes · cryoDRGN",
        "nav_bar_title": "3D volume landscapes",
        "lead_html": LANDSCAPE_FULL_3D_LEAD_HTML,
        "axis_cols": vol_axes,
        "axes_fieldset_legend": "Volume PCA axes",
        "numeric_cols": cols,
        "covariate_display_map": covariate_display_map(cols),
        "default_x": dx,
        "default_y": dy,
        "default_z": dz,
        "scatter3d_url": scatter3d_url,
        "legend_context_body_extra": LANDSCAPE_FULL_3D_LEGEND_CONTEXT_EXTRA,
        "show_vol_landscape_quick_actions": True,
    }


def landscape_full_3d_not_ready_template_kwargs(
    exp_epoch: int, error_message: str = ""
) -> dict[str, Any]:
    """Keyword arguments for ``volume_landscape_3d_need_outputs.html``."""
    return {"exp_epoch": int(exp_epoch), "error_message": error_message}
