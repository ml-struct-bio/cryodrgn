"""Low-level helpers shared by dashboard Plotly and Matplotlib figures."""

from __future__ import annotations

from contextlib import contextmanager

import math
import re
import numpy as np
import pandas as pd
import plotly.graph_objects as go

# Pastel blue for upper-triangle scatter (matplotlib).
_UPPER_SCATTER_BLUE = "#9ec5e8"
# Matches ``base.html`` ``--cream`` (page background).
_DASHBOARD_CREAM = "#faf8f4"

_UMAP_RE = re.compile(r"(?i)^umap")

_PLOTLY_FONT = dict(
    family="Barlow, Roboto, -apple-system, BlinkMacSystemFont, sans-serif",
    size=13,
)

# Dashboard-only discrete palette for labels/k-means (kept out of non-dashboard tools).
_DASHBOARD_CHIMERAX_COLORS: tuple[str, ...] = (
    "#949494",
    "#ffd84d",
    "#81bfd1",
    "#929ad1",
    "#c39bc3",
    "#d1a0a0",
    "#92c29a",
    "#b5a086",
    "#829db1",
    "#b5b547",
)


def _dashboard_chimerax_colors(k: int) -> list[str]:
    if k <= 0:
        return []
    return [
        _DASHBOARD_CHIMERAX_COLORS[i % len(_DASHBOARD_CHIMERAX_COLORS)]
        for i in range(k)
    ]


def _subsample(
    df: pd.DataFrame,
    max_points: int | None,
    *,
    seed: int = 0,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Random subsample of *df*; returns ``(sub_df, original_row_indices)``."""
    n = len(df)
    if max_points is not None and n > max_points:
        pick = np.random.default_rng(seed).choice(
            n,
            size=max_points,
            replace=False,
        )
        return (
            df.iloc[pick].reset_index(drop=True),
            pick.astype(np.int64, copy=False),
        )
    return df, np.arange(n, dtype=np.int64)


def _scatter3d_marker_size_opacity(
    n_points: int, *, point_cap: int
) -> tuple[float, float]:
    """Scale 3D point glyphs for readability — larger / more opaque when fewer points."""
    n = max(1, int(n_points))
    hi_n = float(max(2, int(point_cap)))
    if n >= hi_n:
        return 0.75, 0.4
    t = math.log10(float(n)) / math.log10(hi_n)
    t = max(0.0, min(1.0, t))
    size = 6.0 + t * (0.75 - 6.0)
    opacity = 0.95 + t * (0.4 - 0.95)
    return float(size), float(opacity)


def _mpl_pair_grid_marker_s_alpha(
    n_points: int,
    *,
    point_cap: int,
    dense_s: float,
    dense_a: float,
    sparse_s: float,
    sparse_a: float,
) -> tuple[float, float]:
    """Matplotlib scatter ``s`` / ``alpha`` vs row count (pair-grid cells)."""
    n = max(1, int(n_points))
    hi_n = float(max(2, int(point_cap)))
    if n >= hi_n:
        return dense_s, dense_a
    t = math.log10(float(n)) / math.log10(hi_n)
    t = max(0.0, min(1.0, t))
    s = sparse_s + t * (dense_s - sparse_s)
    a = sparse_a + t * (dense_a - sparse_a)
    return float(s), float(a)


def plot_df_row_indices_for_explorer_scatter(
    plot_df: pd.DataFrame,
    max_points: int = 200_000,
    *,
    seed: int = 0,
) -> np.ndarray:
    """``plot_df`` row positions shown in the particle explorer scatter (``api_scatter`` defaults).

    Matches :func:`scatter_json` / ``api_scatter`` when the request does not set
    ``filter_ui`` or ``full`` — same ``max_points`` cap and fixed RNG ``seed``.

    Returns
    -------
    np.ndarray
        Integer indices into ``plot_df`` (same as ``customdata[:, 1]`` for those points).
    """
    _, row_indices = _subsample(plot_df, max_points, seed=seed)
    return row_indices


def _axes_cell_bboxes(axes: np.ndarray) -> list[dict[str, float]]:
    """Per-cell axis bounding boxes in figure coordinates (row-major)."""
    zdim = int(axes.shape[0])
    out: list[dict[str, float]] = []
    for i in range(zdim):
        for j in range(zdim):
            bb = axes[i, j].get_position()
            out.append(
                {
                    "x0": float(bb.x0),
                    "y0": float(bb.y0),
                    "x1": float(bb.x1),
                    "y1": float(bb.y1),
                }
            )
    return out


def _plotly_to_json(fig: go.Figure) -> str:
    sj = fig.to_json()
    if sj is None:
        raise RuntimeError("Plotly failed to serialize figure.")
    return sj


# ---------------------------------------------------------------------------
# Matplotlib rc settings (used by matplotlib figures)
# ---------------------------------------------------------------------------

# Matplotlib rc typography aligned with https://ezlab.princeton.edu/ (Barlow / Roboto).
EZLAB_SANS_SERIF = [
    "Barlow",
    "Roboto",
    "DejaVu Sans",
    "Helvetica",
    "Arial",
    "sans-serif",
]


@contextmanager
def ezlab_matplotlib_rc():
    """Matplotlib rc context for consistent typography (lazy pyplot import)."""
    import matplotlib.pyplot as plt

    with plt.rc_context(
        {
            "font.family": "sans-serif",
            "font.sans-serif": EZLAB_SANS_SERIF,
        }
    ):
        yield


# ---------------------------------------------------------------------------
# Pair-grid axis limit helpers (stable and symmetric)
# ---------------------------------------------------------------------------


def _pair_grid_finite_axis_lim(values: pd.Series) -> tuple[float, float]:
    """Min/max of finite entries (stable pair-grid axes across colour filters)."""
    v = values.to_numpy(dtype=np.float64)
    m = np.isfinite(v)
    if not np.any(m):
        return -1.0, 1.0
    lo = float(np.min(v[m]))
    hi = float(np.max(v[m]))
    if hi <= lo:
        hi = lo + 1e-9
    return lo, hi


# Fraction of the axis span added as padding on each side (pair grid + pinned 3-D axes).
_PAIR_GRID_AXIS_PAD_FRAC = 0.045


def _pair_grid_axis_lim_padded(
    values: pd.Series,
    *,
    pad_frac: float | None = None,
) -> tuple[float, float]:
    """Finite axis bounds widened symmetrically so points clear panel edges."""
    pf = _PAIR_GRID_AXIS_PAD_FRAC if pad_frac is None else pad_frac
    pf = max(0.0, float(pf))
    lo, hi = _pair_grid_finite_axis_lim(values)
    span = hi - lo
    pad = span * pf if span > 0 else max(abs(lo), abs(hi), 1.0) * pf
    return lo - pad, hi + pad


def _pair_grid_square_xy_limits(
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Align x/y data ranges so span matches on square ``equal``-aspect axes.

    Matplotlib ≥3.10 logs when ``set_xlim`` / ``set_ylim`` are later overridden by
    ``set_aspect(..., adjustable="datalim")``. Matching spans avoids that adjustment.
    """
    lo_x, hi_x = float(xlim[0]), float(xlim[1])
    lo_y, hi_y = float(ylim[0]), float(ylim[1])
    span_x = hi_x - lo_x
    span_y = hi_y - lo_y
    if not math.isfinite(span_x) or span_x <= 0:
        span_x = max(abs(lo_x), abs(hi_x), 1.0) * 1e-9 + 1e-15
        hi_x = lo_x + span_x
    if not math.isfinite(span_y) or span_y <= 0:
        span_y = max(abs(lo_y), abs(hi_y), 1.0) * 1e-9 + 1e-15
        hi_y = lo_y + span_y
    cx = 0.5 * (lo_x + hi_x)
    cy = 0.5 * (lo_y + hi_y)
    if span_x >= span_y:
        half = 0.5 * span_x
        return (lo_x, hi_x), (cy - half, cy + half)
    half = 0.5 * span_y
    return (cx - half, cx + half), (lo_y, hi_y)
