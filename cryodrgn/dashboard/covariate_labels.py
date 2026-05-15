"""Human-readable labels for covariate columns (selectors, legends)."""

from __future__ import annotations

import re
from collections.abc import Iterable, Sequence

import numpy as np


def landscape_vol_pc_pretty_label(
    pc_index_1_based: int,
    explained_variance_ratio: np.ndarray | Sequence[float] | None,
) -> str:
    """Label for a volume PCA component, e.g. ``Vol PC2 (14.2%)`` when ratios are known."""
    i = int(pc_index_1_based)
    base = f"Vol PC{i}"
    if explained_variance_ratio is None or i < 1:
        return base
    evr = np.asarray(explained_variance_ratio, dtype=np.float64).ravel()
    j = i - 1
    if j < 0 or j >= evr.size:
        return base
    val = float(evr[j])
    if not np.isfinite(val):
        return base
    return f"{base} ({100.0 * val:.1f}%)"


def landscape_vol_pc_column_pretty_label(
    column: str,
    explained_variance_ratio: np.ndarray | Sequence[float] | None,
) -> str:
    """Pretty label for ``landscape_vol_PC*`` columns; other names returned unchanged."""
    m = re.fullmatch(r"landscape_vol_PC(\d+)", str(column))
    if not m:
        return str(column)
    return landscape_vol_pc_pretty_label(int(m.group(1)), explained_variance_ratio)


def covariate_display_name(name: str) -> str:
    """Human-friendly covariate names in dashboard selectors."""
    if name == "labels":
        return "k-means labels"
    if name == "landscape_vol_cluster":
        return "Vol cluster"
    m = re.fullmatch(r"landscape_vol_PC(\d+)", name)
    if m:
        return landscape_vol_pc_pretty_label(int(m.group(1)), None)
    return name


def covariate_display_map(
    names: Iterable[str],
    *,
    vol_pc_explained_variance_ratio: np.ndarray | Sequence[float] | None = None,
) -> dict[str, str]:
    """Map column names to display strings for template covariate dropdowns."""
    out: dict[str, str] = {}
    for c in names:
        m = re.fullmatch(r"landscape_vol_PC(\d+)", str(c))
        if m:
            out[c] = landscape_vol_pc_pretty_label(
                int(m.group(1)),
                vol_pc_explained_variance_ratio,
            )
        else:
            out[c] = covariate_display_name(c)
    return out
