"""Plotly and Matplotlib continuous palette names for the dashboard."""

from __future__ import annotations

DEFAULT_DASHBOARD_CONTINUOUS_PALETTE = "Viridis"
_DASHBOARD_PALETTE_TO_MPL: dict[str, str] = {
    "Viridis": "viridis",
    "Plasma": "plasma",
    "Inferno": "inferno",
    "Magma": "magma",
    "Cividis": "cividis",
    "Turbo": "turbo",
    "Blues": "Blues",
    "Greens": "Greens",
    "Greys": "Greys",
    "Oranges": "Oranges",
    "Purples": "Purples",
    "Reds": "Reds",
    "YlGnBu": "YlGnBu",
    "YlOrRd": "YlOrRd",
    "RdBu": "RdBu",
    "Portland": "turbo",
    "Jet": "jet",
    "Hot": "hot",
    "Blackbody": "afmhot",
    "Electric": "magma",
    "Rainbow": "rainbow",
    "Earth": "gist_earth",
}


def normalize_continuous_palette(raw: str | None) -> str:
    """Map a user/API string to a supported Plotly colorscale name."""
    if not raw or not isinstance(raw, str):
        return DEFAULT_DASHBOARD_CONTINUOUS_PALETTE
    s = raw.strip()
    for pk in _DASHBOARD_PALETTE_TO_MPL:
        if pk.lower() == s.lower():
            return pk
    return DEFAULT_DASHBOARD_CONTINUOUS_PALETTE


def mpl_cmap_for_palette(plotly_palette: str) -> str:
    """Matplotlib colormap name for a Plotly palette (falls back to Viridis)."""
    return _DASHBOARD_PALETTE_TO_MPL.get(
        plotly_palette,
        _DASHBOARD_PALETTE_TO_MPL[DEFAULT_DASHBOARD_CONTINUOUS_PALETTE],
    )
