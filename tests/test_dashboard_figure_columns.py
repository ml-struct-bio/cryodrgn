"""Compact "columns" wire format for the scatter endpoints.

The contract under test: for any request, the figure decoded from
``wire=columns`` equals the figure the JSON response carries — coordinates to
float32 precision, everything else exactly (colours, ``customdata`` rows, hover
strings, layout) — and a follow-up request that names the ``base_key`` it holds
gets the colour-dependent columns only.
"""

from __future__ import annotations

import json
import random
import struct
from typing import Any

import numpy as np
import pytest

from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.figure_columns import (
    MAGIC,
    MIMETYPE,
    G5Text,
    NullableFloats,
    Rows,
    decode_figure_columns,
    encode_figure_columns,
)
from cryodrgn.dashboard.palette_config import _DASHBOARD_PALETTE_TO_MPL
from cryodrgn.dashboard.plots_color_covariate import (
    _sample_colorscale_hex,
    _sample_colorscale_hex_reference,
)
from cryodrgn.dashboard.plots_scatter import scatter3d_z_json, scatter_json
from tests.conftest import decode_plotly_figure

pytestmark = pytest.mark.dashboard


def _plain(value: Any) -> Any:
    """Arrays → lists, recursively, so two figures compare with ``==``."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {k: _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _assert_same_figure(json_fig: dict, col_fig: dict) -> None:
    a = _plain(decode_plotly_figure(json_fig))
    b = _plain(col_fig)
    ta, tb = a["data"][0], b["data"][0]
    assert set(ta) == set(tb)
    for axis in ("x", "y", "z"):
        if axis in ta:
            np.testing.assert_allclose(
                np.asarray(tb.pop(axis), dtype=float),
                np.asarray(ta.pop(axis), dtype=float),
                rtol=1e-6,
                atol=1e-9,
            )
    # NaN never equals itself: compare numeric per-point colours separately.
    ca, cb = ta["marker"].get("color"), tb["marker"].get("color")
    if isinstance(ca, list) and ca and isinstance(ca[0], float):
        np.testing.assert_array_equal(np.asarray(cb, float), np.asarray(ca, float))
        ta["marker"].pop("color")
        tb["marker"].pop("color")
    assert ta == tb
    assert a["layout"] == b["layout"]


def _color_columns(exp: DashboardExperiment) -> list[str | None]:
    return [None, *exp.color_covariate_columns]


class TestScatterColumnsMatchFigureJson:
    def test_explorer_scatter_every_colour_column(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        exp = dashboard_experiment
        xcol, ycol = exp.numeric_columns[0], exp.numeric_columns[1]
        for color in _color_columns(exp):
            js = json.loads(scatter_json(exp, xcol, ycol, color, max_points=None))
            body = scatter_json(exp, xcol, ycol, color, max_points=None, wire="columns")
            assert isinstance(body, bytes) and body.startswith(MAGIC)
            fig, _info = decode_figure_columns(body)
            _assert_same_figure(js, fig)

    def test_explorer_scatter_subsample_preselect_and_overrides(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        exp = dashboard_experiment
        n = len(exp.plot_df)
        kw: dict[str, Any] = dict(
            max_points=max(2, n // 2),
            preselect_plot_df_rows=[0, 1, n - 1],
            marker_size=1.8,
            marker_opacity=0.5,
            continuous_palette="Plasma",
            discrete_label_colors={"1": "#ff0000"},
        )
        for color in ("labels", "znorm"):
            js = json.loads(scatter_json(exp, "z0", "z1", color, **kw))
            fig, _ = decode_figure_columns(
                scatter_json(exp, "z0", "z1", color, wire="columns", **kw)
            )
            _assert_same_figure(js, fig)

    def test_latent3d_scatter_colours_and_filters(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        exp = dashboard_experiment
        if int(exp.z.shape[1]) < 3:
            pytest.skip("needs three latent dimensions")
        first_key = str(int(sorted(exp.plot_df["labels"].unique())[0]) + 1)
        variants: list[dict[str, Any]] = [
            dict(color=None),
            dict(color="labels"),
            dict(color="znorm", continuous_palette="Turbo"),
            dict(color="labels", no_subsample=True),
            dict(
                color="labels", color_filter={"kind": "discrete", "keys": [first_key]}
            ),
            dict(
                color="znorm",
                color_filter={
                    "kind": "threshold",
                    "level": float(exp.plot_df["znorm"].median()),
                    "use_max": False,
                },
            ),
        ]
        for v in variants:
            color = v.pop("color")
            js = json.loads(scatter3d_z_json(exp, "z0", "z1", "z2", color, **v))
            fig, _ = decode_figure_columns(
                scatter3d_z_json(exp, "z0", "z1", "z2", color, wire="columns", **v)
            )
            _assert_same_figure(js, fig)


class TestBaseColumnReuse:
    def test_colour_change_omits_coordinates(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        exp = dashboard_experiment
        xcol, ycol = exp.numeric_columns[0], exp.numeric_columns[1]
        first = scatter_json(exp, xcol, ycol, None, max_points=None, wire="columns")
        _fig0, info0 = decode_figure_columns(first)
        assert not info0["base_omitted"]
        assert {"x", "y", "cd.index", "cd.row"} <= set(info0["base"])

        second = scatter_json(
            exp,
            xcol,
            ycol,
            "labels",
            max_points=None,
            wire="columns",
            have=info0["base_key"],
        )
        full = scatter_json(exp, xcol, ycol, "labels", max_points=None, wire="columns")
        assert len(second) < len(full)
        fig2, info2 = decode_figure_columns(second, cached=info0["base"])
        assert info2["base_omitted"] and info2["base_key"] == info0["base_key"]
        js = json.loads(scatter_json(exp, xcol, ycol, "labels", max_points=None))
        _assert_same_figure(js, fig2)

    def test_stale_key_gets_a_complete_response(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        exp = dashboard_experiment
        cols = exp.numeric_columns
        _f, info = decode_figure_columns(
            scatter_json(exp, cols[0], cols[1], None, max_points=None, wire="columns")
        )
        other_axes = scatter_json(
            exp,
            cols[1],
            cols[0],
            None,
            max_points=None,
            wire="columns",
            have=info["base_key"],
        )
        fig, info_b = decode_figure_columns(other_axes)  # decodes without a cache
        assert not info_b["base_omitted"] and info_b["base_key"] != info["base_key"]
        assert len(fig["data"][0]["x"]) == len(exp.plot_df)

    def test_missing_cache_is_an_error_not_a_wrong_figure(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        exp = dashboard_experiment
        cols = exp.numeric_columns
        _f, info = decode_figure_columns(
            scatter_json(exp, cols[0], cols[1], None, max_points=None, wire="columns")
        )
        body = scatter_json(
            exp,
            cols[0],
            cols[1],
            None,
            max_points=None,
            wire="columns",
            have=info["base_key"],
        )
        with pytest.raises(KeyError):
            decode_figure_columns(body)


class TestScatterColumnsApi:
    def test_json_stays_the_default(self, flask_client) -> None:
        r = flask_client.get("/api/scatter?x=z0&y=z1&color=labels")
        assert r.status_code == 200 and r.mimetype == "application/json"
        assert "data" in r.get_json()

    def test_wire_columns_get_and_post(self, flask_client) -> None:
        r = flask_client.get("/api/scatter?x=z0&y=z1&color=labels&wire=columns")
        assert r.status_code == 200 and r.mimetype == MIMETYPE
        fig, info = decode_figure_columns(r.data)
        js = flask_client.get("/api/scatter?x=z0&y=z1&color=labels").get_json()
        _assert_same_figure(js, fig)

        r2 = flask_client.get(
            f"/api/scatter?x=z0&y=z1&color=znorm&wire=columns&have={info['base_key']}"
        )
        fig2, info2 = decode_figure_columns(r2.data, cached=info["base"])
        assert info2["base_omitted"]
        _assert_same_figure(
            flask_client.get("/api/scatter?x=z0&y=z1&color=znorm").get_json(), fig2
        )

        body = {"x": "z0", "y": "z1", "z": "z2", "color": "labels"}
        r3 = flask_client.post("/api/scatter3d_z", json={**body, "wire": "columns"})
        if r3.status_code == 400:
            pytest.skip("fixture has fewer than three latent dimensions")
        assert r3.mimetype == MIMETYPE
        fig3, _ = decode_figure_columns(r3.data)
        _assert_same_figure(
            flask_client.post("/api/scatter3d_z", json=body).get_json(), fig3
        )

    def test_errors_stay_json(self, flask_client) -> None:
        r = flask_client.get("/api/scatter?x=nope&y=z1&wire=columns")
        assert r.status_code == 400 and r.get_json()["error"]


class TestEncoderColumnKinds:
    def _roundtrip(self, trace: dict, n: int, **kw: Any) -> tuple[dict, dict, bytes]:
        body = encode_figure_columns({"title": None}, trace, n_points=n, key="k", **kw)
        fig, info = decode_figure_columns(body)
        return fig["data"][0], info, body

    def test_kinds_and_sharing(self) -> None:
        n = 500
        rng = np.random.default_rng(0)
        vals = rng.normal(size=n)
        vals[::50] = np.nan
        cats = np.asarray([f"c{i % 7}" for i in range(n)], dtype=object)
        trace = dict(
            type="scattergl",
            x=rng.normal(size=n),
            y=rng.normal(size=n).astype(np.float32),
            text=[""] * n,
            marker=dict(color=vals, size=[2.0 if i % 3 else 0.0 for i in range(n)]),
            customdata=Rows(
                [
                    np.arange(n),
                    rng.integers(0, 10**6, n),
                    NullableFloats(vals),
                    cats,
                    G5Text("v: ", "v: —", vals),
                ],
                ["row", "index", "value", "fk", "hover"],
                base=("row", "index"),
            ),
        )
        tr, info, body = self._roundtrip(trace, n)
        header = json.loads(
            body[len(MAGIC) + 4 :][: struct.unpack_from("<I", body, len(MAGIC))[0]]
        )
        kinds = [c["kind"] for c in header["cols"]]
        assert {"arange", "num", "dict", "const", "g5", "rows"} <= set(kinds)
        # marker.color, the customdata value and the hover source are one buffer
        offsets = [
            c["offset"]
            for c in header["cols"]
            if c["kind"] == "num" and c["dtype"] == "f8"
        ]
        assert len(set(offsets)) == 1
        assert tr["x"].dtype == np.float32 and tr["text"] == [""] * n
        assert tr["marker"]["size"] == trace["marker"]["size"]
        row = tr["customdata"][50]
        assert row[0] == 50 and row[2] is None and row[3] == "c1" and row[4] == "v: —"
        assert tr["customdata"][1][4] == "v: " + "%.5g" % vals[1]
        assert set(info["base"]) == {"x", "y", "cd.row", "cd.index"}

    def test_empty_figure(self) -> None:
        tr, _info, _ = self._roundtrip(
            dict(
                type="scattergl",
                x=np.zeros(0),
                y=np.zeros(0),
                customdata=Rows([np.zeros(0, int)], ["row"]),
            ),
            0,
        )
        assert len(tr["x"]) == 0 and tr["customdata"] == []

    def test_g5_text_matches_c_formatting(self) -> None:
        rng = random.Random(3)
        vals = [
            0.0,
            12344.5,
            12345.5,
            99999.5,
            1e-5,
            9.99995e-5,
            123456.0,
            float("nan"),
            float("inf"),
        ]
        vals += [rng.uniform(-1, 1) * 10 ** rng.randint(-9, 9) for _ in range(2000)]
        out = G5Text("p=", "none", np.asarray(vals)).to_list()
        for v, s in zip(vals, out):
            assert s == ("p=" + "%.5g" % v if np.isfinite(v) else "none")


class TestVectorizedColorscale:
    @pytest.mark.parametrize("palette", sorted(_DASHBOARD_PALETTE_TO_MPL))
    def test_identical_to_plotly_sample_colorscale(self, palette: str) -> None:
        rng = np.random.default_rng(7)
        t = np.concatenate(
            [
                rng.random(4000),
                np.linspace(0, 1, 511),
                np.arange(256) / 255.0,
                [0.0, 1.0],
            ]
        )
        assert _sample_colorscale_hex(palette, t) == _sample_colorscale_hex_reference(
            palette, t
        )
