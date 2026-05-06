"""HTTP helper utilities shared by :mod:`cryodrgn.dashboard.app` route handlers."""

from __future__ import annotations

import os
from collections.abc import Iterable
from typing import Any

import numpy as np
from flask import jsonify, redirect, url_for

from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.particle_explorer import explorer_volumes_eligible
from cryodrgn.dashboard.plots import normalize_continuous_palette
from cryodrgn.dashboard.preload import DEFAULT_PRELOAD_IMAGE_LIMIT
from cryodrgn.dashboard.trajectory import (
    direct_anchor_particle_indices_payload,
    has_pc_columns,
    has_umap_columns,
)

_TRAJECTORY_INELIGIBLE_MSG = (
    "Trajectory creator needs a CUDA GPU, single-particle data, "
    "and weights for the current epoch."
)
_EXPLORER_VOLUMES_INELIGIBLE_MSG = (
    "Volume explorer needs a CUDA GPU, single-particle data, "
    "and weights for the current epoch."
)


def _filter_ui_scatter_max_points() -> int:
    """Cap for ``/api/scatter`` when ``filter_ui=1`` (env override available)."""
    try:
        v = int(os.environ.get("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", "500000"))
    except ValueError:
        v = 500_000
    return max(50_000, min(v, 2_000_000))


def _particle_explorer_scatter_max_points() -> int:
    """Cap for particle explorer scatter + thumbnail restrict pool.

    When ``CRYODRGN_DASHBOARD_FILTER_MAX_POINTS`` is set (e.g. ``cryodrgn dashboard
    … --filter-max N``), use that clamped cap. Otherwise keep the historical
    200k default for ``/api/scatter`` without ``filter_ui``.
    """
    raw = (os.environ.get("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS") or "").strip()
    if raw:
        try:
            return max(50_000, min(int(raw), 2_000_000))
        except ValueError:
            pass
    return 200_000


def _particle_explorer_scatter_cap_from_env() -> bool:
    """True when the explorer scatter cap comes from FILTER_MAX_POINTS (CLI or env)."""
    raw = (os.environ.get("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS") or "").strip()
    if not raw:
        return False
    try:
        int(raw)
    except ValueError:
        return False
    return True


def _default_xy_cols(cols: list[str]) -> tuple[str, str]:
    """Pick sensible default X/Y axes (UMAP if available, else first two)."""
    x = "UMAP1" if "UMAP1" in cols else cols[0]
    y = "UMAP2" if "UMAP2" in cols else cols[min(1, len(cols) - 1)]
    return x, y


def _covariate_display_name(name: str) -> str:
    """Human-friendly covariate names in dashboard selectors."""
    if name == "labels":
        return "k-means labels"
    return name


def _covariate_display_map(names: Iterable[str]) -> dict[str, str]:
    """Map column names to display strings for template covariate dropdowns."""
    return {c: _covariate_display_name(c) for c in names}


def _parse_preselect_rows_param(raw: str | None) -> tuple[list[int] | None, str | None]:
    """Parse ``preselect_rows`` query (comma-separated ints). Returns ``(rows, err)``."""
    s = (raw or "").strip()
    if not s:
        return None, None
    try:
        return [int(p) for p in s.split(",") if p.strip()][:5000], None
    except ValueError as exc:
        return None, f"invalid preselect_rows: {exc}"


def _parse_preload_image_limit(raw: object) -> int:
    """Positive thumbnail cache size from a request value, defaulting to 5000."""
    if raw is None or raw == "":
        return DEFAULT_PRELOAD_IMAGE_LIMIT
    if isinstance(raw, bool) or not isinstance(raw, (int, float, str)):
        raise ValueError("cache_size must be a positive integer.")
    if isinstance(raw, float) and not raw.is_integer():
        raise ValueError("cache_size must be a positive integer.")
    try:
        value = int(raw)
    except (TypeError, ValueError) as exc:
        raise ValueError("cache_size must be a positive integer.") from exc
    if value < 1:
        raise ValueError("cache_size must be a positive integer.")
    return value


def _redirect(endpoint: str):
    """302 redirect to a named view (same-origin)."""
    return redirect(url_for(endpoint), code=302)


def _trajectory_eligibility_error(e: DashboardExperiment):
    """JSON 400 when the trajectory creator is not available, else ``None``."""
    if explorer_volumes_eligible(e):
        return None
    return jsonify(error=_TRAJECTORY_INELIGIBLE_MSG), 400


def _parse_pairplot_request(
    e: DashboardExperiment, payload: dict
) -> tuple[str, str, str, str]:
    color_col = payload.get("color_col") or payload.get("lower_color_col")
    if not color_col or not isinstance(color_col, str):
        raise ValueError("Choose a color covariate.")
    if color_col not in e.plot_df.columns:
        raise ValueError("Invalid color column.")
    if color_col in {f"z{i}" for i in range(int(e.z.shape[1]))}:
        raise ValueError("Latent z columns cannot be used as the color covariate.")
    raw_diag = payload.get("diagonal_emb")
    if raw_diag is None or (isinstance(raw_diag, str) and raw_diag.strip() == ""):
        diagonal_emb = "umap" if has_umap_columns(e) else "pc"
    else:
        diagonal_emb = str(raw_diag).lower()
    upper_style = (payload.get("upper_style") or "scatter").lower()
    if diagonal_emb not in ("pc", "umap"):
        raise ValueError("diagonal_emb must be pc or umap.")
    if upper_style not in ("scatter", "hex"):
        raise ValueError("upper_style must be scatter or hex.")
    if diagonal_emb == "umap" and not has_umap_columns(e):
        raise ValueError("UMAP is not available for this run.")
    if diagonal_emb == "pc" and not has_pc_columns(e):
        raise ValueError("PCA components are not available.")
    raw_palette = payload.get("palette")
    pair_palette = normalize_continuous_palette(
        str(raw_palette) if raw_palette is not None else None,
    )
    return color_col, diagonal_emb, upper_style, pair_palette


def _parse_color_filter_for_column(
    color_col: str, payload: dict
) -> dict[str, Any] | None:
    """Parse optional ``color_filter`` from a JSON body (threshold or discrete keys)."""
    raw = payload.get("color_filter")
    if raw in (None, "", {}):
        return None
    if not isinstance(raw, dict):
        raise ValueError("color_filter must be an object.")
    kind = raw.get("kind")
    if kind == "threshold":
        try:
            level = float(raw["level"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("Invalid colour threshold level.") from exc
        use_max = bool(raw.get("use_max"))
        return {"kind": "threshold", "level": level, "use_max": use_max}
    if kind == "discrete":
        keys = raw.get("keys")
        if not isinstance(keys, list):
            raise ValueError("Discrete colour filter requires a keys array.")
        return {"kind": "discrete", "keys": [str(k) for k in keys]}
    raise ValueError("color_filter kind must be threshold or discrete.")


def _add_direct_anchor_pidx(payload: dict, p: dict, z_traj: np.ndarray) -> None:
    """Merge direct-anchor particle IDs into ``payload`` when applicable."""
    if not (p.get("use_anchors") and p["mode"] == "direct"):
        return
    pidx = direct_anchor_particle_indices_payload(
        anchor_indices=p["anchor_indices"],
        interpolation_points=p["n_points"],
        n_total=int(np.asarray(z_traj).shape[0]),
    )
    if pidx is not None:
        payload["traj_particle_indices"] = pidx
