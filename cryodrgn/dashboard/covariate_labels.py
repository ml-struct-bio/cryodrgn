"""Human-readable labels for covariate columns (selectors, legends)."""

from __future__ import annotations

import re
from collections.abc import Iterable


def covariate_display_name(name: str) -> str:
    """Human-friendly covariate names in dashboard selectors."""
    if name == "labels":
        return "k-means labels"
    m = re.fullmatch(r"landscape_vol_PC(\d+)", name)
    if m:
        return f"Volume PCA {m.group(1)}"
    return name


def covariate_display_map(names: Iterable[str]) -> dict[str, str]:
    """Map column names to display strings for template covariate dropdowns."""
    return {c: covariate_display_name(c) for c in names}
