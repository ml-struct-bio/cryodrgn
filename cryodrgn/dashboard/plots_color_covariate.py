"""Covariate colouring, filters, and legend payloads for dashboard plots."""

from __future__ import annotations

import re
from typing import Any, cast

import numpy as np
import pandas as pd

from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.plots_figure_utils import (
    _DASHBOARD_CREAM,
    _dashboard_chimerax_colors,
    _mpl_pair_grid_marker_s_alpha,
)


def _lower_color_series_is_discrete(s: pd.Series) -> bool:
    """Treat integer / categorical / few whole-valued floats as discrete coloring."""
    s = s.dropna()
    if len(s) == 0:
        return False
    if pd.api.types.is_bool_dtype(s):
        return True
    if isinstance(s.dtype, pd.CategoricalDtype):
        return True
    if pd.api.types.is_integer_dtype(s) and not pd.api.types.is_float_dtype(s):
        return True
    if s.dtype == object:
        return True
    if pd.api.types.is_float_dtype(s):
        u = np.unique(s.to_numpy())
        if len(u) > 40:
            return False
        if len(u) == 0:
            return False
        return bool(
            np.all(np.isfinite(u)) and np.allclose(u, np.round(u), rtol=0, atol=1e-8)
        )
    return False


def _continuous_series_stats(values: pd.Series) -> tuple[np.ndarray, float, float]:
    """Return numeric array plus robust (min,max) bounds for continuous coloring."""
    color_num = cast(pd.Series, pd.to_numeric(values, errors="coerce"))
    cvals = np.asarray(color_num, dtype=np.float64)
    cfinite = cvals[np.isfinite(cvals)]
    if cfinite.size == 0:
        cmin, cmax = 0.0, 1.0
    elif np.isclose(cfinite.min(), cfinite.max()):
        cmin = float(cfinite.min()) - 0.5
        cmax = float(cfinite.max()) + 0.5
    else:
        cmin = float(np.min(cfinite))
        cmax = float(np.max(cfinite))
    return cvals, cmin, cmax


def _lower_legend_entry_label(lower_color_col: str, u: Any) -> str:
    """1-based cluster indices for the ``labels`` (k-means) column."""
    if pd.isna(u):
        return str(u)
    if lower_color_col == "labels":
        if isinstance(u, (bool, np.bool_)):
            return str(u)
        if isinstance(u, (int, np.integer)):
            return str(int(u) + 1)
        if isinstance(u, float) and np.isfinite(u) and u == int(u):
            return str(int(u) + 1)
        try:
            fu = float(u)
            if np.isfinite(fu) and fu == int(fu):
                return str(int(fu) + 1)
        except (TypeError, ValueError):
            pass
    return str(u)


def _labels_colors_and_legend_items(
    values: pd.Series,
) -> tuple[list[str], list[dict[str, str]]]:
    """Point colors and legend items for k-means-style integer labels."""
    codes_arr, uniques = pd.factorize(values, sort=True)
    codes = np.asarray(codes_arr, dtype=np.int64)
    pal = _dashboard_chimerax_colors(max(len(uniques), 1))
    colors: list[str] = []
    for code in codes:
        if int(code) < 0:
            colors.append("#aab4bf")
        else:
            colors.append(pal[int(code) % len(pal)])
    items: list[dict[str, str]] = []
    for idx, u in enumerate(uniques):
        if pd.isna(u):
            continue
        items.append(
            {
                "label": _lower_legend_entry_label("labels", u),
                "color": pal[idx % len(pal)],
            }
        )
    if bool(values.isna().any()):
        items.append({"label": "(missing)", "color": "#aab4bf"})
    return colors, items


def covariate_row_filter_key(color_col: str, u: Any) -> str:
    """Stable string key for a covariate value (matches dashboard scatter customdata)."""
    if pd.isna(u):
        return "__na__"
    if color_col == "labels":
        return str(_lower_legend_entry_label(color_col, u))
    num = pd.to_numeric(pd.Series([u], dtype=object), errors="coerce").iloc[0]
    if pd.notna(num):
        fv = float(num)
        if np.isfinite(fv) and abs(fv - round(fv)) < 1e-8:
            return str(int(round(fv)))
        return str(fv)
    return str(u)


def plot_df_color_filter_mask(
    df: pd.DataFrame,
    color_col: str,
    color_filter: dict[str, Any] | None,
) -> np.ndarray:
    """Boolean mask aligned with ``df`` rows (True = include)."""
    if not color_filter:
        return np.ones(len(df), dtype=bool)
    kind = color_filter.get("kind")
    if kind == "threshold":
        level = float(color_filter["level"])
        use_max = bool(color_filter.get("use_max"))
        series = pd.to_numeric(df[color_col], errors="coerce")
        if use_max:
            mask = (series <= level) & series.notna()
        else:
            mask = (series >= level) & series.notna()
        return mask.to_numpy(dtype=bool)
    if kind == "range":
        r0 = float(color_filter["range_min"])
        r1 = float(color_filter["range_max"])
        if r0 > r1:
            r0, r1 = r1, r0
        series = pd.to_numeric(df[color_col], errors="coerce")
        inside = (series >= r0) & (series <= r1) & series.notna()
        if bool(color_filter.get("invert_range")):
            mask = (~inside) & series.notna()
        else:
            mask = inside
        return mask.to_numpy(dtype=bool)
    if kind == "discrete":
        keys_inc = frozenset(str(k) for k in color_filter.get("keys") or ())
        ser = df[color_col]
        keys_series = ser.map(lambda u, cc=color_col: covariate_row_filter_key(cc, u))
        return keys_series.isin(keys_inc).to_numpy(dtype=bool)
    return np.ones(len(df), dtype=bool)


def _discrete_legend_sort_tuple(key: str) -> tuple:
    if key == "__na__":
        return (2, "")
    try:
        return (0, float(key))
    except ValueError:
        return (1, str(key))


def _stable_discrete_covariate_hex_map(
    plot_df: pd.DataFrame,
    color_col: str,
    discrete_label_colors: dict[str, str] | None,
) -> dict[str, str]:
    """Filter-key → marker hex; palette order matches :func:`covariate_legend_context_payload`."""
    raw_full = plot_df[color_col]
    _, uniques_f = pd.factorize(raw_full, sort=True)
    pal = _dashboard_chimerax_colors(max(len(uniques_f), 1))

    def hex_for_key(fk: str, palette_idx: int) -> str:
        override = discrete_label_colors.get(fk) if discrete_label_colors else None
        if isinstance(override, str) and re.match(
            r"^#[0-9a-fA-F]{6}$", override.strip()
        ):
            return override.strip().lower()
        return pal[palette_idx % len(pal)]

    stable_hex_by_key: dict[str, str] = {}
    for idx, u in enumerate(uniques_f):
        if pd.isna(u):
            continue
        fk = covariate_row_filter_key(color_col, u)
        stable_hex_by_key[fk] = hex_for_key(fk, idx)

    if bool(raw_full.isna().any()):
        fk_na = "__na__"
        ov_na = discrete_label_colors.get(fk_na) if discrete_label_colors else None
        if isinstance(ov_na, str) and re.match(r"^#[0-9a-fA-F]{6}$", ov_na.strip()):
            stable_hex_by_key[fk_na] = ov_na.strip().lower()
        else:
            stable_hex_by_key[fk_na] = "#aab4bf"

    return stable_hex_by_key


def _scatter_discrete_marker_arrays(
    plot_df: pd.DataFrame,
    sub_df: pd.DataFrame,
    color_col: str,
    discrete_label_colors: dict[str, str] | None,
) -> tuple[list[str], list[str]]:
    """Per-row ChimeraX hex colours + covariate filter keys for Plotly ``marker.color`` / ``customdata``."""
    lookup = _stable_discrete_covariate_hex_map(
        plot_df, color_col, discrete_label_colors
    )
    ser = sub_df[color_col]
    n = len(sub_df)
    colors: list[str] = []
    keys: list[str] = []
    for i in range(n):
        v = ser.iloc[i]
        if pd.isna(v):
            fk = "__na__"
        else:
            fk = covariate_row_filter_key(color_col, v)
        colors.append(lookup.get(fk, "#aab4bf"))
        keys.append(fk)
    return colors, keys


def discrete_category_counts_by_filter_key(
    plot_df: pd.DataFrame, color_col: str
) -> dict[str, int]:
    """Covariate filter-key → particle count in the full ``plot_df`` analysis table."""
    if color_col not in plot_df.columns:
        return {}
    s = plot_df[color_col]
    if not _lower_color_series_is_discrete(s):
        return {}
    codes_arr, uniques = pd.factorize(s, sort=True)
    codes = np.asarray(codes_arr, dtype=np.int64)
    out: dict[str, int] = {}
    for idx, u in enumerate(uniques):
        if pd.isna(u):
            continue
        fk = covariate_row_filter_key(color_col, u)
        out[fk] = int(np.sum(codes == idx))
    if bool(s.isna().any()):
        out["__na__"] = int(s.isna().sum())
    return out


def covariate_legend_context_payload(
    exp: DashboardExperiment,
    column: str,
    *,
    plot_df: pd.DataFrame | None = None,
    violin_sample_max: int = 12_000,
) -> dict[str, Any]:
    """JSONable payload for building the colour histogram / discrete toggles."""
    df = exp.plot_df if plot_df is None else plot_df
    if column not in df.columns:
        raise ValueError("Unknown column.")
    s = df[column]
    if _lower_color_series_is_discrete(s):
        _codes_arr, uniques = pd.factorize(s, sort=True)
        lookup = _stable_discrete_covariate_hex_map(df, column, None)
        counts_map = discrete_category_counts_by_filter_key(df, column)
        cats: list[dict[str, Any]] = []
        for idx, u in enumerate(uniques):
            if pd.isna(u):
                continue
            fk = covariate_row_filter_key(column, u)
            cnt = counts_map.get(fk, 0)
            if cnt == 0:
                continue
            cats.append(
                {
                    "key": fk,
                    "label": _lower_legend_entry_label(column, u),
                    "count": cnt,
                    "color": lookup.get(fk, "#aab4bf"),
                }
            )
        if bool(s.isna().any()):
            cats.append(
                {
                    "key": "__na__",
                    "label": "(missing)",
                    "count": counts_map.get("__na__", int(s.isna().sum())),
                    "color": lookup.get("__na__", "#aab4bf"),
                }
            )
        cats.sort(key=lambda c: _discrete_legend_sort_tuple(str(c["key"])))
        n_pts = len(df)
        point_cap = max(n_pts, 1)
        _, chip_blend_alpha = _mpl_pair_grid_marker_s_alpha(
            n_pts,
            point_cap=point_cap,
            dense_s=2.5,
            dense_a=0.2,
            sparse_s=10.0,
            sparse_a=0.78,
        )
        return {
            "mode": "discrete",
            "categories": cats,
            "chip_blend_alpha": float(chip_blend_alpha),
            "chip_blend_bg": _DASHBOARD_CREAM,
        }
    vals = pd.to_numeric(s, errors="coerce").dropna().to_numpy(dtype=np.float64)
    if vals.size == 0:
        return {"mode": "continuous", "values": []}
    rng = np.random.default_rng(0)
    if vals.size > violin_sample_max:
        vals = rng.choice(vals, size=violin_sample_max, replace=False)
    return {"mode": "continuous", "values": vals.tolist()}
