"""Serve Plotly.js from the installed Python ``plotly`` package (offline-safe)."""

from __future__ import annotations

from functools import lru_cache
from importlib.resources import files as pkg_files

from flask import Response


@lru_cache(maxsize=1)
def bundled_plotly_js_version() -> str:
    """Plotly.js version shipped with the installed ``plotly`` wheel."""
    from plotly.offline.offline import get_plotlyjs_version

    return get_plotlyjs_version()


@lru_cache(maxsize=1)
def bundled_plotly_min_js_bytes() -> bytes:
    """Minified Plotly.js bytes from ``plotly.package_data``."""
    path = pkg_files("plotly").joinpath("package_data", "plotly.min.js")
    return path.read_bytes()


def bundled_plotly_js() -> Response:
    """Flask view: serve bundled Plotly.js for dashboard pages."""
    body = bundled_plotly_min_js_bytes()
    return Response(
        body,
        mimetype="application/javascript",
        headers={
            "Cache-Control": "public, max-age=86400",
            "X-Plotly-Version": bundled_plotly_js_version(),
        },
    )
