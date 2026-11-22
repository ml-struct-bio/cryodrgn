"""Pair plot interface tests."""

from __future__ import annotations

import json
import os
import re

import numpy as np
import pandas as pd
import pytest

from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.plots import (
    _continuous_series_stats,
    pair_grid_figure_aspect_ratio,
    pair_grid_margin_fractions_for_js,
    pair_grid_png,
    pair_grid_skeleton_placeholder_layout,
)
from cryodrgn.dashboard.plots_color_covariate import covariate_row_filter_key

pytestmark = pytest.mark.dashboard


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


class TestPairGridHelpers:
    @pytest.mark.parametrize(
        "values,cmin,cmax",
        [
            ([2.0, 2.0, 2.0], None, None),  # constant → finite span
            ([np.nan, np.nan], 0.0, 1.0),
            ([1.0, np.nan, 3.0, 5.0], 1.0, 5.0),
        ],
        ids=["constant", "all_nan", "mixed"],
    )
    def test_continuous_series_stats(
        self, values: list, cmin: float | None, cmax: float | None
    ) -> None:
        s = pd.Series(values)
        vals, lo, hi = _continuous_series_stats(s)
        if cmin is None:
            assert lo < hi
        else:
            assert lo == cmin and hi == cmax
        if all(isinstance(v, float) and np.isnan(v) for v in values):
            assert np.isnan(vals).all()

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


class TestPairplotBrowserSmoke:
    """Headless Chromium: pair-grid PNG refresh and upper-triangle style toggle."""

    pytestmark = pytest.mark.browser

    def test_pairplot_image_loads_and_hex_style_redraws(
        self, playwright_page, dashboard_live_url
    ) -> None:
        from tests.conftest import dashboard_smoke_pairplot

        out = dashboard_smoke_pairplot(playwright_page, dashboard_live_url)
        assert out["src_changed"]

    def test_color_legend_fills_the_screen_with_the_menu_top_margin(
        self, playwright_page, dashboard_live_url
    ) -> None:
        """The right-hand legend is as tall as the screen, inset by the left menu's top gap."""
        from tests.conftest import _dashboard_smoke_wait_pairplot_image

        page = playwright_page
        page.goto(f"{dashboard_live_url.rstrip('/')}/pairplot")
        _dashboard_smoke_wait_pairplot_image(page, timeout_ms=60_000)
        page.wait_for_function(
            """() => {
              var panel = document.getElementById('pair-color-legend-panel');
              return panel && !panel.hidden && panel.getBoundingClientRect().height > 100;
            }"""
        )
        box = page.evaluate(
            """() => {
              var card = document.getElementById('pair-controls-card');
              var panel = document.getElementById('pair-color-legend-panel');
              var cardTop = card.getBoundingClientRect().top;
              var panelBox = panel.getBoundingClientRect();
              return {
                cardTop: cardTop,
                panelTop: panelBox.top,
                bottomGap: window.innerHeight - panelBox.bottom,
                scrollY: window.scrollY || 0
              };
            }"""
        )
        assert box["scrollY"] == 0
        assert abs(box["panelTop"] - box["cardTop"]) <= 2
        assert abs(box["bottomGap"] - box["cardTop"]) <= 2
        visible = page.evaluate(
            """() => {
              var sw = document.getElementById('pair-color-discrete-switches');
              var panel = document.getElementById('pair-color-legend-panel');
              var nodes = sw.querySelectorAll('.cryo-cc-discrete-cell');
              var panelBox = panel.getBoundingClientRect();
              var hidden = 0;
              for (var i = 0; i < nodes.length; i++) {
                var r = nodes[i].getBoundingClientRect();
                var inside = r.height > 2
                  && r.top >= panelBox.top - 1
                  && r.bottom <= panelBox.bottom + 1
                  && r.left >= panelBox.left - 1
                  && r.right <= panelBox.right + 1;
                if (!inside) hidden++;
              }
              var fits = sw.scrollHeight <= window.innerHeight - panelBox.top - 12;
              return {n: nodes.length, hidden: hidden, fits: fits};
            }"""
        )
        assert visible["n"] > 0
        assert visible["fits"]
        assert visible["hidden"] == 0
        page.evaluate(
            """() => {
              var sw = document.getElementById('pair-color-discrete-switches');
              var sample = sw.querySelector('.cryo-cc-discrete-cell');
              for (var i = 0; i < 7; i++) {
                var cell = sample.cloneNode(true);
                cell.setAttribute('data-clone', '1');
                sw.appendChild(cell);
              }
              window.dispatchEvent(new Event('resize'));
            }"""
        )
        filled = page.evaluate(
            """() => {
              var sw = document.getElementById('pair-color-discrete-switches');
              var panel = document.getElementById('pair-color-legend-panel');
              var nodes = sw.querySelectorAll('.cryo-cc-discrete-cell');
              var panelBox = panel.getBoundingClientRect();
              var hidden = 0;
              for (var i = 0; i < nodes.length; i++) {
                var r = nodes[i].getBoundingClientRect();
                var inside = r.height > 2
                  && r.top >= panelBox.top - 1
                  && r.bottom <= panelBox.bottom + 1
                  && r.right <= panelBox.right + 1;
                if (!inside) hidden++;
              }
              return {
                n: nodes.length,
                hidden: hidden,
                scrolls: sw.scrollHeight > sw.clientHeight + 1
              };
            }"""
        )
        assert filled["n"] >= visible["n"]
        assert filled["hidden"] == 0, filled
        assert not filled["scrolls"], filled

    def test_continuous_legend_uses_the_screen_height(
        self, playwright_page, dashboard_live_url
    ) -> None:
        """A continuous covariate histogram fills the legend down to the screen."""
        from tests.conftest import _dashboard_smoke_wait_pairplot_image

        page = playwright_page
        page.goto(f"{dashboard_live_url.rstrip('/')}/pairplot")
        _dashboard_smoke_wait_pairplot_image(page, timeout_ms=60_000)
        page.click('input[name="color_cov"][value="znorm"]')
        page.wait_for_function(
            """() => {
              var wrap = document.getElementById('pair-color-continuous-wrap');
              var hist = document.getElementById('pair-color-hist');
              return wrap && !wrap.hidden && hist
                && hist.getBoundingClientRect().height > 160;
            }"""
        )
        box = page.evaluate(
            """() => {
              var card = document.getElementById('pair-controls-card');
              var row = document.querySelector('.cryo-dash-page--pairplot .cryo-dash-row');
              var panel = document.getElementById('pair-color-legend-panel');
              var hist = document.getElementById('pair-color-hist');
              var footer = panel.querySelector('.cryo-cc-legend__footer--continuous');
              var panelBox = panel.getBoundingClientRect();
              var histBox = hist.getBoundingClientRect();
              var footerTop = footer.getBoundingClientRect().top;
              return {
                rowTop: row.getBoundingClientRect().top,
                cardTop: card.getBoundingClientRect().top,
                panelTop: panelBox.top,
                bottomGap: window.innerHeight - panelBox.bottom,
                histGap: footerTop - histBox.bottom,
                histH: histBox.height
              };
            }"""
        )
        assert abs(box["panelTop"] - box["rowTop"]) <= 2
        assert box["cardTop"] - box["panelTop"] > 24
        assert box["bottomGap"] <= 16
        assert box["histH"] > 160
        assert box["histGap"] <= 24


class TestPairGridPointRaster:
    """Pair-grid points are composited into images instead of drawn as 64 scatters.

    That took a grid of 84k points from 10.6 s to 3.4 s (1.6 s on 4 threads). The
    images must look like Matplotlib's scatter: same marker size and opacity, drawn in
    the same order.
    """

    @staticmethod
    def _render(mode: str, x, y, rgb, size: float, alpha: float) -> np.ndarray:
        import matplotlib.pyplot as plt

        from cryodrgn.dashboard.plots_pair_grid import (
            _add_points_image,
            _composite_points,
            _panel_pixels,
        )

        fig = plt.figure(figsize=(2.5, 2.5), dpi=100)
        ax = fig.add_axes((0, 0, 1, 1))
        ax.set_xlim(-3, 3)
        ax.set_ylim(-3, 3)
        ax.set_axis_off()
        px = _panel_pixels(ax)
        if mode == "raster":
            img = _composite_points(px, x, y, ("rgb", rgb), alpha, size, 100)
            _add_points_image(ax, px, img)
        else:
            # Matplotlib with every point moved to the centre of its pixel, as the
            # raster places them.
            s = px["bw"] / 6.0
            x = (np.floor(px["bx0"] + (x + 3) * s) + 0.5 - px["bx0"]) / s - 3
            y = (np.floor(px["by0"] + (y + 3) * s) + 0.5 - px["by0"]) / s - 3
            ax.scatter(x, y, c=rgb, s=size, alpha=alpha, linewidths=0)
        fig.canvas.draw()
        image = np.asarray(fig.canvas.buffer_rgba())[..., :3].astype(int)
        plt.close(fig)
        return image

    @pytest.mark.parametrize(
        "points, size, alpha", [(3000, 6.0, 0.3), (20000, 2.5, 0.2)]
    )
    def test_composites_like_matplotlib(self, points: int, size: float, alpha: float):
        """Same marker area, rim and stacking as ``scatter``; only the snapping of
        points to pixel centres differs (2-4 levels on average against unsnapped
        Matplotlib, the same as unsnapped against snapped Matplotlib)."""
        rng = np.random.default_rng(3)
        x, y = rng.normal(size=points), rng.normal(size=points)
        rgb = rng.random((points, 3))
        mpl = self._render("matplotlib", x, y, rgb, size, alpha)
        raster = self._render("raster", x, y, rgb, size, alpha)
        diff = np.abs(mpl - raster)
        assert diff.mean() < 1.0
        assert np.percentile(diff, 99) <= 6

    def test_points_drawn_later_are_on_top(self) -> None:
        from cryodrgn.dashboard.plots_pair_grid import _composite_points

        px = {
            "x0": 0.0, "x1": 10.0, "y0": 0.0, "y1": 10.0,
            "bx0": 0.0, "by0": 0.0, "bw": 10.0, "bh": 10.0,
            "c0": 0, "c1": 10, "r0": 0, "r1": 10,
        }  # fmt: skip
        red, blue = [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]
        x = y = np.array([5.5, 5.5])
        size = 0.36  # at 72 dpi: 0.6 px across, 28% of one pixel
        img = _composite_points(px, x, y, ("rgb", np.array([red, blue])), 0.5, size, 72)
        top = _composite_points(px, x, y, ("rgb", np.array([blue, red])), 0.5, size, 72)
        pixel = img[4, 5]
        assert pixel[3] > 0
        assert pixel[2] > pixel[0], "blue was drawn last"
        assert top[4, 5][0] > top[4, 5][2], "red was drawn last"
        # The palette-index path composites exactly like per-point colours.
        lut = np.array([red, blue])
        indexed = _composite_points(
            px, x, y, ("lut", np.array([0, 1]), lut), 0.5, size, 72
        )
        assert np.allclose(indexed, img)

    def test_marker_footprint_covers_the_marker_area(self) -> None:
        from cryodrgn.dashboard.plots_pair_grid import _marker_footprint

        for diameter in (0.8, 2.36, 2.64, 6.3):
            _, _, frac = _marker_footprint(diameter)
            assert frac.sum() == pytest.approx(np.pi * diameter**2 / 4, rel=0.05)

    def test_threads_do_not_change_the_image(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        kw = dict(lower_color_col="labels", diagonal_emb="pc", upper_style="scatter")
        one, _ = pair_grid_png(dashboard_experiment, **kw, workers=1)
        four, _ = pair_grid_png(dashboard_experiment, **kw, workers=4)
        assert one == four

    def test_switching_back_and_saving_reuse_the_grid(
        self, flask_client, monkeypatch: pytest.MonkeyPatch, tmp_path
    ) -> None:
        from cryodrgn.dashboard import plots_pair_grid

        renders = []
        real = plots_pair_grid.pair_grid_png
        monkeypatch.setattr(
            plots_pair_grid,
            "pair_grid_png",
            lambda *a, **k: renders.append(k.get("upper_style")) or real(*a, **k),
        )
        monkeypatch.setattr(
            plots_pair_grid,
            "_PAIR_GRID_CACHE",
            type(plots_pair_grid._PAIR_GRID_CACHE)(),
        )
        body = {"color_col": "znorm", "diagonal_emb": "pc", "upper_style": "scatter"}
        first = flask_client.post("/api/pairplot", json=body).get_json()
        flask_client.post("/api/pairplot", json={**body, "upper_style": "hex"})
        again = flask_client.post("/api/pairplot", json=body).get_json()
        assert again["png_b64"] == first["png_b64"]
        saved = flask_client.post(
            "/api/save_pairplot_png", json={**body, "filename": "reused.png"}
        )
        assert saved.status_code == 200
        assert renders == ["scatter", "hex"], "the scatter grid was rendered once"


class TestPairplotStaleReplies:
    """A reply for settings the user has left behind must not reach the screen.

    Hex, then straight back to Scatter: the hex image used to land under the Scatter
    button, because going back to the image on screen sent no request.
    """

    pytestmark = pytest.mark.browser

    def test_going_back_to_the_image_on_screen_drops_the_request_in_flight(
        self, playwright_page, dashboard_live_url
    ) -> None:
        from tests.conftest import _dashboard_smoke_wait_pairplot_image

        page = playwright_page
        held = []

        def hold_hex(route) -> None:
            body = json.loads(route.request.post_data or "{}")
            if body.get("upper_style") == "hex":
                held.append(route)
            else:
                route.continue_()

        page.route("**/api/pairplot", hold_hex)
        try:
            page.goto(f"{dashboard_live_url.rstrip('/')}/pairplot")
            src0 = _dashboard_smoke_wait_pairplot_image(page, timeout_ms=60_000)
            page.click('input[name="upper_style"][value="hex"]')
            page.wait_for_timeout(300)
            assert len(held) == 1, "the hex request must be in flight"
            page.click('input[name="upper_style"][value="scatter"]')
            page.wait_for_timeout(300)
            try:
                held[0].continue_()  # the page has aborted it; let it through anyway
            except Exception:
                pass
            page.wait_for_timeout(2500)
            state = page.evaluate(
                """() => ({
                  src: document.getElementById('pairplot').src,
                  busy: document.getElementById('pairplot-viewport')
                    .getAttribute('aria-busy')
                })"""
            )
            assert state["src"] == src0, "the hex image replaced the scatter grid"
            assert state["busy"] == "false"
        finally:
            page.unroute("**/api/pairplot")


class TestPairplotColourChoiceConsistency:
    """Every entry in the pair-grid colour menu must render and describe itself.

    ``pairplot_page`` builds ``color_choices`` by filtering ``color_covariate_columns``,
    while ``/api/pairplot`` and the legend endpoint validate colours independently.
    The existing coverage spot-checks three columns, which cannot detect a menu entry
    that no longer renders.

    The menu is read back from the rendered radio group rather than recomputed, so it
    reflects whatever covariates the running server actually has loaded.
    """

    @staticmethod
    def _menu(client) -> tuple[list[str], list[str]]:
        body = client.get("/pairplot").get_data(as_text=True)
        choices = [
            c
            for c in re.findall(r'name="color_cov" value="([^"]+)"', body)
            if c != "none"
        ]
        assert choices, "pair grid offered no colour choices"
        discrete = re.search(r"var discreteColorCols = (\[.*?\]);", body, re.S)
        return choices, json.loads(discrete.group(1)) if discrete else []

    def test_every_colour_choice_renders_a_grid(self, flask_client) -> None:
        choices, _ = self._menu(flask_client)
        for column in choices:
            r = flask_client.post("/api/pairplot", json={"color_col": column})
            assert r.status_code == 200, (
                f"menu offers colour {column!r} but the grid will not render it: "
                f"{r.get_data(as_text=True)[:200]}"
            )
            assert r.get_json().get("png_b64")

    def test_menu_and_legend_agree_on_which_choices_are_discrete(
        self, flask_client
    ) -> None:
        """The page picks swatches from one list; the legend answers from another."""
        choices, discrete = self._menu(flask_client)
        for column in choices:
            payload = flask_client.post(
                "/api/covariate_legend_context", json={"column": column}
            ).get_json()
            expected = "discrete" if column in discrete else "continuous"
            assert payload.get("mode") == expected, (
                f"pair grid treats {column!r} as {expected} but the legend endpoint "
                f"reports {payload.get('mode')!r}"
            )

    def test_latent_axes_are_kept_out_of_the_colour_menu(self, flask_client) -> None:
        """``z*`` are the grid's own axes, so offering them as colours is meaningless."""
        choices, _ = self._menu(flask_client)
        assert not [c for c in choices if re.fullmatch(r"z\d+", c)]
