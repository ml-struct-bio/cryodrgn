"""Load particle-indexed covariate numpy arrays from ``.pkl`` files into ``plot_df``."""

from __future__ import annotations

import os
import re
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd

from cryodrgn import utils

if TYPE_CHECKING:
    from cryodrgn.dashboard.data import DashboardExperiment

_COL_SANITIZE_RE = re.compile(r"[^\w]+")
CovariateKind = Literal["numeric", "labels"]


def _sanitize_column_stem(stem: str) -> str:
    s = _COL_SANITIZE_RE.sub("_", stem).strip("_")
    if not s:
        s = "covariate"
    if s[0].isdigit():
        s = f"cov_{s}"
    return s


def _unique_column_name(plot_df_columns: set[str], base: str) -> str:
    if base not in plot_df_columns:
        return base
    n = 2
    while f"{base}_{n}" in plot_df_columns:
        n += 1
    return f"{base}_{n}"


def _parse_particle_covariate_raw(
    raw: object, path: str
) -> tuple[np.ndarray, CovariateKind]:
    """Validate pickle contents and classify as numeric or discrete label values."""
    if isinstance(raw, dict):
        raise ValueError(f"Expected an array in {path!r}; got dict.")
    if isinstance(raw, np.ndarray):
        arr = raw
    elif isinstance(raw, (list, tuple)):
        try:
            arr = np.asarray(raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Expected an array in {path!r}; got {type(raw).__name__}."
            ) from exc
    else:
        raise ValueError(
            f"Expected a numpy array in {path!r}; got {type(raw).__name__}."
        )
    if arr.ndim not in (1, 2):
        raise ValueError(
            "Covariate array must be 1-D (one value per particle) or 2-D "
            f"(one row per particle); got shape {tuple(arr.shape)!r}."
        )
    if np.issubdtype(arr.dtype, np.number):
        return np.asarray(arr), "numeric"
    if (
        arr.dtype == object
        or pd.api.types.is_string_dtype(arr)
        or arr.dtype.kind in "OUSb"
    ):
        return np.asarray(arr, dtype=object), "labels"
    raise ValueError(
        f"Covariate values must be numeric or text labels; got dtype {arr.dtype!r}."
    )


def _label_series_from_column(values: np.ndarray) -> pd.Series:
    """Normalize a 1-D label column for discrete dashboard colouring."""

    def _one(value: object) -> object:
        if value is None:
            return pd.NA
        if isinstance(value, float) and np.isnan(value):
            return pd.NA
        if isinstance(value, bytes):
            return value.decode("utf-8", errors="replace")
        return str(value)

    flat = np.asarray(values, dtype=object).ravel()
    return pd.Series((_one(v) for v in flat), dtype="string")


def invalidate_covariate_column_caches(exp: DashboardExperiment) -> None:
    """Drop cached column lists after mutating ``plot_df``."""
    exp.__dict__.pop("numeric_columns", None)
    exp.__dict__.pop("color_covariate_columns", None)


def merge_covariate_pkl(
    exp: DashboardExperiment, path: str
) -> tuple[list[str], list[str]]:
    """Load ``path`` and append column(s) to ``exp.plot_df``.

    Returns ``(column_names, discrete_column_names)``. Non-numeric label arrays are
    always stored as discrete text covariates. Raises :class:`ValueError` when the
    file format does not match dashboard particle-index expectations.
    """
    abs_path = os.path.abspath(path)
    if not abs_path.lower().endswith(".pkl"):
        raise ValueError("Select a .pkl file.")
    if not os.path.isfile(abs_path):
        raise ValueError(f"File not found on server: {abs_path}")

    raw = utils.load_pkl(abs_path)
    arr, kind = _parse_particle_covariate_raw(raw, abs_path)
    n_particles = len(exp.plot_df)
    if arr.shape[0] != n_particles:
        raise ValueError(
            f"Covariate length {arr.shape[0]} does not match particle count "
            f"{n_particles} (plot_df row count)."
        )

    prefix = _sanitize_column_stem(os.path.splitext(os.path.basename(abs_path))[0])
    existing = set(exp.plot_df.columns)
    new_cols: list[str] = []
    discrete_cols: list[str] = []

    if arr.ndim == 1:
        col = _unique_column_name(existing, prefix)
        if kind == "numeric":
            exp.plot_df[col] = np.asarray(arr, dtype=np.float64).ravel()
        else:
            exp.plot_df[col] = _label_series_from_column(arr)
            discrete_cols.append(col)
            if col not in exp.user_covariate_columns:
                exp.user_covariate_columns.append(col)
        new_cols.append(col)
    else:
        for j in range(arr.shape[1]):
            base = f"{prefix}_{j}" if arr.shape[1] > 1 else prefix
            col = _unique_column_name(existing, base)
            existing.add(col)
            if kind == "numeric":
                exp.plot_df[col] = np.asarray(arr[:, j], dtype=np.float64)
            else:
                exp.plot_df[col] = _label_series_from_column(arr[:, j])
                discrete_cols.append(col)
                if col not in exp.user_covariate_columns:
                    exp.user_covariate_columns.append(col)
            new_cols.append(col)

    invalidate_covariate_column_caches(exp)
    return new_cols, discrete_cols
