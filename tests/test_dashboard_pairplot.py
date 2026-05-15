"""Pair plot interface tests."""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest

from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.plots import (
    _continuous_series_stats,
    mpl_cmap_for_palette,
    normalize_continuous_palette,
    pair_grid_figure_aspect_ratio,
    pair_grid_margin_fractions_for_js,
    pair_grid_png,
    pair_grid_skeleton_placeholder_layout,
)
from cryodrgn.dashboard.plots_color_covariate import covariate_row_filter_key


class TestDashboardPairPlot:
    """``pair_grid_png`` rendering via both the Flask route and direct import."""

    @pytest.mark.parametrize(
        "color_col,diagonal_emb,upper_style",
        [
            ("labels", "umap", "scatter"),
            ("znorm", "pc", "hex"),
            ("UMAP1", "pc", "scatter"),
        ],
    )
    def test_api_pairplot(
        self,
        flask_client,
        color_col: str,
        diagonal_emb: str,
        upper_style: str,
    ) -> None:
        r = flask_client.post(
            "/api/pairplot",
            json={
                "color_col": color_col,
                "diagonal_emb": diagonal_emb,
                "upper_style": upper_style,
            },
        )
        assert r.status_code == 200, r.get_json()
        js = r.get_json()
        # The route returns the PNG as base64 alongside per-cell geometry.
        assert js.get("png_b64"), "pairplot response missing png_b64"
        assert js.get("cells"), "pairplot response missing cells"
        assert len(js["cells"]) == 16  # zdim=4 -> 4x4 grid

    def test_api_pairplot_accepts_continuous_palette(self, flask_client) -> None:
        r = flask_client.post(
            "/api/pairplot",
            json={
                "color_col": "znorm",
                "diagonal_emb": "pc",
                "upper_style": "hex",
                "palette": "plasma",
            },
        )
        assert r.status_code == 200, r.get_json()
        assert r.get_json().get("png_b64")

    def test_api_pairplot_accepts_discrete_color_filter(
        self, flask_client, dashboard_experiment: DashboardExperiment
    ) -> None:
        """First discrete selection must change the PNG (server path + client lastPayload sync)."""
        s = dashboard_experiment.plot_df["labels"].dropna()
        if s.nunique() < 2:
            pytest.skip("fixture needs >=2 label values for discrete filter regression")
        u = s.value_counts().index[-1]
        fk = covariate_row_filter_key("labels", u)
        base = {"color_col": "labels", "diagonal_emb": "umap", "upper_style": "scatter"}
        r0 = flask_client.post("/api/pairplot", json=base)
        r1 = flask_client.post(
            "/api/pairplot",
            json={**base, "color_filter": {"kind": "discrete", "keys": [fk]}},
        )
        assert r0.status_code == 200, r0.get_json()
        assert r1.status_code == 200, r1.get_json()
        assert r0.get_json().get("png_b64") != r1.get_json().get("png_b64")

    def test_pair_grid_png_is_deterministic(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        """Same inputs should give byte-identical PNG output."""
        png_a, cells_a = pair_grid_png(
            dashboard_experiment,
            lower_color_col="labels",
            diagonal_emb="umap",
            upper_style="scatter",
        )
        png_b, cells_b = pair_grid_png(
            dashboard_experiment,
            lower_color_col="labels",
            diagonal_emb="umap",
            upper_style="scatter",
        )
        assert png_a == png_b, "pair_grid_png is non-deterministic"
        assert cells_a == cells_b
        # zdim == 4 -> 4x4 cells.
        assert len(cells_a) == dashboard_experiment.z.shape[1] ** 2


class TestNormalizeContinuousPalette:
    @pytest.mark.parametrize(
        "raw,expected",
        [
            (None, "Viridis"),
            ("", "Viridis"),
            ("viridis", "Viridis"),
            ("TURBO", "Turbo"),
            ("  Plasma  ", "Plasma"),
            ("not_a_palette", "Viridis"),
            (123, "Viridis"),
        ],
    )
    def test_cases(self, raw: object, expected: str) -> None:
        pal_raw = raw if (raw is None or isinstance(raw, str)) else str(raw)
        assert normalize_continuous_palette(pal_raw) == expected

    def test_mpl_cmap_for_palette(self) -> None:
        assert mpl_cmap_for_palette("Viridis") == "viridis"
        # Unknown falls back to viridis.
        assert mpl_cmap_for_palette("not_a_palette") == "viridis"


class TestContinuousSeriesStats:
    def test_constant_series_has_finite_span(self) -> None:
        s = pd.Series([2.0, 2.0, 2.0])
        _, cmin, cmax = _continuous_series_stats(s)
        assert cmin < cmax

    def test_all_nan_falls_back(self) -> None:
        s = pd.Series([np.nan, np.nan])
        vals, cmin, cmax = _continuous_series_stats(s)
        assert cmin == 0.0 and cmax == 1.0
        assert np.isnan(vals).all()

    def test_mixed_series_uses_extrema(self) -> None:
        s = pd.Series([1.0, np.nan, 3.0, 5.0])
        _, cmin, cmax = _continuous_series_stats(s)
        assert cmin == 1.0 and cmax == 5.0


class TestPairGridHexAndSkeleton:
    def test_hex_style_is_png_and_deterministic(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        png_a, cells_a = pair_grid_png(
            dashboard_experiment,
            lower_color_col="labels",
            diagonal_emb="pc",
            upper_style="hex",
        )
        png_b, cells_b = pair_grid_png(
            dashboard_experiment,
            lower_color_col="labels",
            diagonal_emb="pc",
            upper_style="hex",
        )
        assert png_a[:8] == b"\x89PNG\r\n\x1a\n"
        assert png_a == png_b
        assert cells_a == cells_b

    def test_placeholder_layout_shape(self) -> None:
        cells = pair_grid_skeleton_placeholder_layout(3)
        assert len(cells) == 9
        for c in cells:
            assert {"x0", "y0", "x1", "y1"} <= set(c.keys())

    def test_placeholder_zdim_zero_is_empty(self) -> None:
        assert pair_grid_skeleton_placeholder_layout(0) == []

    def test_margin_fractions_for_js(self) -> None:
        m = pair_grid_margin_fractions_for_js()
        for k in (
            "left_m",
            "top_m",
            "bottom_m",
            "right_axes",
            "edge_gap",
            "right_gap",
        ):
            assert k in m
            assert isinstance(m[k], float)

    def test_placeholder_cell_bboxes_match_png(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        zdim = int(dashboard_experiment.z.shape[1])
        ph = pair_grid_skeleton_placeholder_layout(zdim)
        _, cells = pair_grid_png(
            dashboard_experiment,
            lower_color_col="labels",
            diagonal_emb="pc",
            upper_style="scatter",
        )
        assert len(ph) == len(cells)
        for a, b in zip(ph, cells, strict=True):
            for k in ("x0", "y0", "x1", "y1"):
                assert abs(float(a[k]) - float(b[k])) < 1e-9

    def test_figure_aspect_ratio_uses_zdim(self) -> None:
        assert pair_grid_figure_aspect_ratio(0) == 1.0
        assert pair_grid_figure_aspect_ratio(3) > 0


class TestSavePairPlotPng:
    def test_writes_png(self, flask_client, tmp_path) -> None:
        r = flask_client.post(
            "/api/save_pairplot_png",
            json={
                "color_col": "labels",
                "diagonal_emb": "pc",
                "upper_style": "scatter",
                "filename": "pytest_pairplot.png",
            },
        )
        assert r.status_code == 200
        js = r.get_json()
        assert js["ok"] is True
        assert os.path.isfile(js["path"])
        with open(js["path"], "rb") as fh:
            assert fh.read(8) == b"\x89PNG\r\n\x1a\n"
