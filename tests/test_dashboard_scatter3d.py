"""Dashboard 3-D scatterplot (Plotly ``scatter3d``) API, plot JSON, and client contracts.

Covers ``scatter3d_z_json``, ``/api/scatter3d_z*``, latent-3D preview/GIF routes, and
``plotly_scatter3d_scene.js`` / ``latent_3d.html`` orbit-preservation behaviour.
"""

from __future__ import annotations

import base64
import json
from io import BytesIO
from typing import Any

import numpy as np
import pandas as pd
import pytest
from PIL import Image
from plotly.colors import sample_colorscale

from cryodrgn.dashboard import plots_scatter as plots_scatter_mod
from cryodrgn.dashboard.column_names import (
    VOL_LANDSCAPE_3D_PLOT_DF_ROW,
    VOL_LANDSCAPE_IS_SKETCH_CENTROID,
    VOL_LANDSCAPE_NEAREST_SKETCH_VOL,
)
from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.plots_color_covariate import (
    covariate_row_filter_key,
    numeric_array_to_plotly_hex,
    plotly_color_to_hex,
    _stable_discrete_covariate_hex_map,
)
from cryodrgn.dashboard.plots_figure_utils import (
    _subsample_preserving_sketch_centroids,
)
from cryodrgn.dashboard.plots_scatter import (
    _SCATTER3D_FILTERED_MARKER_SIZE_ARRAY_BOOST,
    _scatter3d_apply_filter_visibility_sizes,
    _scatter3d_glyph_count_for_filter,
    _scatter3d_marker_sizes_with_legend_filter,
    scatter3d_discrete_level_png_bytes,
    scatter3d_z_json,
)
from cryodrgn.dashboard.plot_gif_utils import png_base64_frames_to_gif_bytes
from tests.conftest import (
    js_function_body,
    png_b64_rgb,
    read_dashboard_static_js,
    read_latent_3d_html,
)

LATENT_Z_AXES = frozenset({"z0", "z1", "z2"})
DEFAULT_POINT_CAP = 120_000


def _scatter3d_figure(
    exp: DashboardExperiment,
    x: str,
    y: str,
    z: str,
    color: str | None,
    **kwargs: Any,
) -> dict[str, Any]:
    return json.loads(scatter3d_z_json(exp, x, y, z, color, **kwargs))


def _expected_dashboard_scatter3d_glyph(
    n_glyph: int,
    *,
    point_cap: int = DEFAULT_POINT_CAP,
    volume_landscape_3d_style: bool = False,
) -> tuple[float, float]:
    from cryodrgn.dashboard.plots_figure_utils import _scatter3d_marker_size_opacity
    from cryodrgn.dashboard.plots_scatter import _dashboard_scatter3d_glyph_visual_scale

    msize, mopacity = _scatter3d_marker_size_opacity(n_glyph, point_cap=point_cap)
    if volume_landscape_3d_style:
        msize *= (1.0 - 0.31) * (1.0 - 0.13) * (1.0 - 0.13) * 0.8
    msize, mopacity = _dashboard_scatter3d_glyph_visual_scale(msize, mopacity)
    if volume_landscape_3d_style:
        msize *= 1.13 * 1.728 * 1.11
        mopacity = float(max(0.0, min(1.0, mopacity * 0.81)))
    return msize, mopacity


def _trace_visible_glyph(trace: dict[str, Any]) -> tuple[float, float, int, int]:
    n_total = len(trace["x"])
    op = float(trace["marker"]["opacity"])
    sz = trace["marker"]["size"]
    if isinstance(sz, list):
        vis = [float(s) for s in sz if float(s) > 0.0]
        assert vis, "expected at least one visible point"
        assert all(s == vis[0] for s in vis), "visible points share one glyph size"
        return vis[0], op, len(vis), n_total
    return float(sz), op, n_total, n_total


def _assert_filter_preserves_plotly_xyz(unfiltered: dict, filtered: dict) -> None:
    u = unfiltered["data"][0]
    f = filtered["data"][0]
    assert len(f["x"]) == len(u["x"])
    assert f["x"] == u["x"]
    assert f["y"] == u["y"]
    assert f["z"] == u["z"]
    for axis_name in ("xaxis", "yaxis", "zaxis"):
        scene = filtered["layout"]["scene"]
        assert "range" not in (scene.get(axis_name) or {})
    sizes = f["marker"]["size"]
    assert isinstance(sizes, list)
    assert len(sizes) == len(f["x"])
    assert any(float(s) == 0.0 for s in sizes)
    assert any(float(s) > 0.0 for s in sizes)
    colors = f["marker"]["color"]
    assert isinstance(colors, list)
    assert not any(str(c).strip().lower().startswith("rgba") for c in colors[:50])
    op = f["marker"]["opacity"]
    assert isinstance(op, (int, float))
    vis_sizes = [float(s) for s in sizes if float(s) > 0.0]
    u_sz = u["marker"]["size"]
    u_base = float(u_sz) if not isinstance(u_sz, list) else float(max(u_sz))
    if vis_sizes and len(vis_sizes) == len(sizes):
        assert max(vis_sizes) == pytest.approx(u_base)
    elif vis_sizes and len(vis_sizes) < len(sizes):
        assert max(vis_sizes) > u_base


@pytest.fixture(scope="module")
def plotly_scatter3d_scene_js() -> str:
    return read_dashboard_static_js("plotly_scatter3d_scene.js")


class TestScatter3dPlotGifUtils:
    """``plot_gif_utils.png_base64_frames_to_gif_bytes`` (latent-3D GIF frame assembly)."""

    def test_png_base64_frames_to_gif_bytes_round_trip(self) -> None:
        a = png_b64_rgb(rgb=(200, 30, 40))
        b = png_b64_rgb(rgb=(30, 180, 60))
        gif = png_base64_frames_to_gif_bytes([a, b], durations_ms=[50, 120])
        assert gif[:6] in (b"GIF87a", b"GIF89a")
        im = Image.open(BytesIO(gif))
        assert im.n_frames >= 2
        im.seek(1)

    def test_png_base64_frames_requires_two_frames(self) -> None:
        one = png_b64_rgb()
        with pytest.raises(ValueError, match="At least two"):
            png_base64_frames_to_gif_bytes([one])

    def test_png_base64_accepts_data_url_prefix(self) -> None:
        raw = png_b64_rgb()
        framed = "data:image/png;base64," + raw
        gif = png_base64_frames_to_gif_bytes([framed, raw], durations_ms=40)
        assert len(gif) > 32


class TestScatter3dApiEndpoints:
    """``/api/scatter3d_z*`` and latent-3D GIF / preview routes."""

    def test_api_scatter3d_z(self, flask_client) -> None:
        r = flask_client.get("/api/scatter3d_z?x=z0&y=z1&z=z2&color=znorm")
        assert r.status_code == 200
        assert r.get_json()["data"]

    def test_api_scatter3d_z_landscape_full_without_outputs_is_400(
        self, flask_client
    ) -> None:
        r = flask_client.get(
            "/api/scatter3d_z_landscape_full?x=z0&y=z1&z=z2&color=none"
        )
        assert r.status_code == 400
        err = r.get_json().get("error", "")
        assert (
            "vol_pca" in err.lower()
            or "landscape" in err.lower()
            or "three" in err.lower()
        )

    def test_api_latent3d_landscape_full_discrete_gif_without_outputs_is_400(
        self, flask_client
    ) -> None:
        r = flask_client.post(
            "/api/latent3d_landscape_full_discrete_gif",
            json={
                "x": "z0",
                "y": "z1",
                "z": "z2",
                "color": "labels",
                "discrete_keys": ["0", "1"],
            },
        )
        assert r.status_code == 400
        err = (r.get_json() or {}).get("error", "")
        assert "landscape" in err.lower() or "vol_pca" in err.lower()

    def test_api_latent3d_preview_png(self, flask_client) -> None:
        r = flask_client.get("/api/latent3d_preview.png?x=z0&y=z1&z=z2&color=znorm")
        assert r.status_code == 200
        assert r.data[:8] == b"\x89PNG\r\n\x1a\n"

    def test_api_latent3d_plot_gif_from_png_frames(self, flask_client) -> None:
        a = png_b64_rgb(rgb=(10, 20, 30))
        b = png_b64_rgb(rgb=(200, 180, 40))
        r = flask_client.post(
            "/api/latent3d_plot_gif_from_png_frames",
            json={"frames": [a, b], "durations_ms": [40, 80]},
        )
        assert r.status_code == 200, r.get_data(as_text=True)
        js = r.get_json()
        assert "gif_b64" in js
        raw = base64.standard_b64decode(js["gif_b64"])
        assert raw[:6] in (b"GIF87a", b"GIF89a")

    def test_api_latent3d_discrete_gif_requires_discrete_keys(
        self, flask_client
    ) -> None:
        r = flask_client.post(
            "/api/latent3d_discrete_gif",
            json={"x": "z0", "y": "z1", "z": "z2", "color": "labels"},
        )
        assert r.status_code == 400
        err = (r.get_json() or {}).get("error", "")
        assert "discrete_keys" in err.lower()


class TestScatter3dApiErrorPaths:
    def test_scatter3d_bad_color(self, flask_client) -> None:
        r = flask_client.get("/api/scatter3d_z?x=z0&y=z1&z=z2&color=does_not_exist")
        assert r.status_code == 400

    def test_latent3d_non_numeric_elev(self, flask_client) -> None:
        r = flask_client.get("/api/latent3d_preview.png?x=z0&y=z1&z=z2&elev=nope")
        assert r.status_code == 400


class TestScatter3dPlotJson:
    """``scatter3d_z_json`` figure contracts (server-side Plotly)."""

    def test_scatter3d_volume_landscape_style_scales_marker_size(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        """Regression: 3D volume-landscape UI keeps glyphs smaller than the default n-curve (relative factor unchanged after dashboard glyph scale)."""

        base = _scatter3d_figure(
            dashboard_experiment,
            "z0",
            "z1",
            "z2",
            None,
            volume_landscape_3d_style=False,
        )
        vol = _scatter3d_figure(
            dashboard_experiment,
            "z0",
            "z1",
            "z2",
            None,
            volume_landscape_3d_style=True,
        )
        sz_b = base["data"][0]["marker"]["size"]
        sz_v = vol["data"][0]["marker"]["size"]
        assert isinstance(sz_b, (int, float))
        assert isinstance(sz_v, (int, float))
        factor = (1.0 - 0.31) * (1.0 - 0.13) ** 2 * 0.8 * 1.13 * 1.728 * 1.11
        assert float(sz_v) == pytest.approx(float(sz_b) * factor, rel=1e-6, abs=1e-9)
        assert float(sz_v) < float(sz_b)

    def test_numeric_array_to_plotly_hex_batched_matches_per_point(
        self,
    ) -> None:
        """Batched ``sample_colorscale`` must match legacy per-point mapping."""
        vals = np.array([0.0, 0.5, 1.0, np.nan, -0.2, 1.2], dtype=np.float64)
        batched = numeric_array_to_plotly_hex(vals, "Viridis", vmin=0.0, vmax=1.0)
        legacy: list[str] = []
        for i in range(len(vals)):
            v = vals[i]
            if not np.isfinite(v):
                legacy.append("#9ca3af")
                continue
            t = max(0.0, min(1.0, float(v)))
            legacy.append(plotly_color_to_hex(sample_colorscale("Viridis", [t])[0]))
        assert batched == legacy
        assert all(str(c).startswith("#") for c in batched)

    def test_scatter3d_continuous_covariate_uses_per_point_hex_marker_colors(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        """Continuous scatter3d markers must be per-point hex so letter overlays match WebGL."""

        e = dashboard_experiment
        allow = frozenset({"z0", "z1", "z2"})
        fig = _scatter3d_figure(
            e,
            "z0",
            "z1",
            "z2",
            "z0",
            plot_df=e.plot_df.iloc[:80],
            xyz_axes_allowed=allow,
            volume_landscape_3d_style=True,
            continuous_palette="Viridis",
        )
        mk = fig["data"][0]["marker"]
        colors = mk["color"]
        assert isinstance(colors, list)
        assert len(colors) == 80
        assert all(
            isinstance(c, str) and (c.startswith("#") or c.startswith("rgb"))
            for c in colors
        )
        assert "colorscale" not in mk
        meta = (fig.get("layout") or {}).get("meta") or {}
        assert meta.get("cdrgn_color_mode") == "continuous"
        assert meta.get("cdrgn_continuous_palette") == "Viridis"

    def test_scatter3d_vol_landscape_overlay_trace_uses_markers_plus_text(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        """3D volume-landscape GIF UI: trace starts as markers+text so selection restyle skips mode flip."""
        e = dashboard_experiment
        sub = e.plot_df.iloc[:40].copy()
        sub[VOL_LANDSCAPE_NEAREST_SKETCH_VOL] = (
            np.arange(len(sub), dtype=np.int64) % 3 + 1
        )
        fig = _scatter3d_figure(
            e,
            "z0",
            "z1",
            "z2",
            "none",
            plot_df=sub,
            xyz_axes_allowed=LATENT_Z_AXES,
            volume_landscape_3d_style=True,
        )
        tr = fig["data"][0]
        assert tr["mode"] == "markers+text"
        assert tr["text"] == [""] * len(sub)
        assert tr.get("textposition") == "middle center"
        tf = tr.get("textfont") or {}
        assert tf.get("size") == pytest.approx(
            36 * 0.8 * 0.75 * 1.5 * 0.8 * 0.8 * 0.8 * 0.8
        )
        assert tf.get("color") == "#1a1a1a"
        meta = (fig.get("layout") or {}).get("meta") or {}
        assert meta.get("cdrgn_landscape_vol_animation") is True
        assert meta.get("cdrgn_landscape_sketch_centroid_cd") is not True

    def test_scatter3d_vol_landscape_sketch_centroid_customdata(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        e = dashboard_experiment
        sub = e.plot_df.iloc[:40].copy()
        sub[VOL_LANDSCAPE_3D_PLOT_DF_ROW] = np.arange(40, dtype=np.int64)
        sub[VOL_LANDSCAPE_NEAREST_SKETCH_VOL] = (
            np.arange(len(sub), dtype=np.int64) % 3 + 1
        )
        sub[VOL_LANDSCAPE_IS_SKETCH_CENTROID] = (
            (np.arange(len(sub), dtype=np.int64) % 5 == 0)
        ).astype(np.int64)
        allow = LATENT_Z_AXES
        fig = _scatter3d_figure(
            e,
            "z0",
            "z1",
            "z2",
            "none",
            plot_df=sub,
            xyz_axes_allowed=allow,
            volume_landscape_3d_style=True,
        )
        meta = (fig.get("layout") or {}).get("meta") or {}
        assert meta.get("cdrgn_landscape_sketch_centroid_cd") is True
        row0 = fig["data"][0]["customdata"][0]
        assert int(row0[-1]) in (1, 2, 3)
        assert int(row0[-2]) in (0, 1)

    def test_subsample_preserving_sketch_centroids_keeps_centroid_rows(self) -> None:
        n = 200
        cent = np.zeros(n, dtype=np.int64)
        cent[[5, 42]] = 1
        df = pd.DataFrame({VOL_LANDSCAPE_IS_SKETCH_CENTROID: cent})
        sub, idx = _subsample_preserving_sketch_centroids(df, 20, seed=1)
        assert len(sub) == 20
        assert int(sub[VOL_LANDSCAPE_IS_SKETCH_CENTROID].sum()) == 2
        assert 5 in idx and 42 in idx

    def test_scatter3d_without_vol_overlay_trace_stays_markers_only(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        e = dashboard_experiment
        sub = e.plot_df.iloc[:40].copy()
        fig = _scatter3d_figure(
            e,
            "z0",
            "z1",
            "z2",
            "none",
            plot_df=sub,
            xyz_axes_allowed=LATENT_Z_AXES,
            volume_landscape_3d_style=True,
        )
        tr = fig["data"][0]
        assert tr["mode"] == "markers"
        assert "text" not in tr
        meta = (fig.get("layout") or {}).get("meta") or {}
        assert meta.get("cdrgn_landscape_vol_animation") is not True


class TestScatter3dLegendFilterAxisStability:
    """All interactive Plotly 3-D views use ``scatter3d_z_json`` (fixed subsample + size mask)."""

    def test_scatter3d_filter_helpers_exported(self) -> None:
        assert hasattr(plots_scatter_mod, "_scatter3d_subsample_full_table")
        assert hasattr(plots_scatter_mod, "_scatter3d_filter_visibility_on_subsample")
        assert hasattr(plots_scatter_mod, "_scatter3d_apply_filter_visibility_sizes")
        assert hasattr(plots_scatter_mod, "_scatter3d_glyph_count_for_filter")
        assert hasattr(plots_scatter_mod, "_scatter3d_marker_sizes_with_legend_filter")

    def test_scatter3d_apply_filter_visibility_sizes_keeps_scalar_when_all_visible(
        self,
    ) -> None:
        import numpy as np

        marker = {"size": 2.5, "opacity": 0.4}
        all_vis = np.ones(8, dtype=bool)
        _scatter3d_apply_filter_visibility_sizes(marker, all_vis, 2.5)
        assert marker["size"] == 2.5
        part = all_vis.copy()
        part[0] = False
        _scatter3d_apply_filter_visibility_sizes(marker, part, 2.5, 4.0)
        assert isinstance(marker["size"], list)
        assert marker["size"][0] == 0.0
        assert marker["size"][1:] == [4.0] * 7
        assert marker.get("sizemode") == "diameter"
        assert marker.get("sizeref") == 1.0

    def test_scatter3d_filtered_visible_size_scales_from_subsample_baseline(
        self,
    ) -> None:
        import numpy as np

        n_sub = 10_000
        part = np.ones(n_sub, dtype=bool)
        part[:1000] = False
        base, _, vis = _scatter3d_marker_sizes_with_legend_filter(
            n_sub,
            part,
            point_cap=120_000,
            volume_landscape_3d_style=False,
        )
        assert vis > base
        assert vis >= base * _SCATTER3D_FILTERED_MARKER_SIZE_ARRAY_BOOST
        all_on = np.ones(n_sub, dtype=bool)
        base2, _, vis2 = _scatter3d_marker_sizes_with_legend_filter(
            n_sub,
            all_on,
            point_cap=120_000,
            volume_landscape_3d_style=False,
        )
        assert vis2 == base2

    def test_scatter3d_ignores_discrete_filter_on_continuous_column(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        """Stale discrete keys after a covariate switch must not empty the plot."""

        e = dashboard_experiment
        fk = covariate_row_filter_key("labels", e.plot_df["labels"].dropna().iloc[0])
        stale = {"kind": "discrete", "keys": [fk]}
        unfiltered = _scatter3d_figure(e, "z0", "z1", "z2", "znorm", color_filter=None)
        mismatched = _scatter3d_figure(e, "z0", "z1", "z2", "znorm", color_filter=stale)
        u = unfiltered["data"][0]
        m = mismatched["data"][0]
        assert m["x"] == u["x"]
        assert m["marker"]["size"] == u["marker"]["size"]
        assert m["marker"]["opacity"] == u["marker"]["opacity"]

    def test_scatter3d_glyph_count_for_filter_curve(self) -> None:
        import numpy as np

        assert _scatter3d_glyph_count_for_filter(1000, None) == 1000
        all_vis = np.ones(1000, dtype=bool)
        assert _scatter3d_glyph_count_for_filter(1000, all_vis) == 1000
        part = all_vis.copy()
        part[::2] = False
        assert _scatter3d_glyph_count_for_filter(1000, part) == 500

    def test_scatter3d_full_subsample_filter_matches_no_colour_glyph(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        """All legend toggles on (every subsample row visible) → same glyph size as colour=none."""

        e = dashboard_experiment
        sub = e.plot_df.iloc[:500].copy()
        allow = frozenset({"z0", "z1", "z2"})
        base_kw = dict(plot_df=sub, xyz_axes_allowed=allow)
        no_colour = _scatter3d_figure(
            e, "z0", "z1", "z2", "none", color_filter=None, **base_kw
        )
        lookup = _stable_discrete_covariate_hex_map(e.plot_df, "labels", None)
        all_keys = sorted(lookup.keys())
        all_on = _scatter3d_figure(
            e,
            "z0",
            "z1",
            "z2",
            "labels",
            color_filter={"kind": "discrete", "keys": all_keys},
            **base_kw,
        )
        u_sz, u_op, _, _ = _trace_visible_glyph(no_colour["data"][0])
        f_tr = all_on["data"][0]
        f_sz, f_op, n_vis, n_total = _trace_visible_glyph(f_tr)
        assert n_vis == n_total
        assert f_sz == pytest.approx(u_sz)
        assert f_op == pytest.approx(u_op)
        assert not isinstance(f_tr["marker"]["size"], list)

    def test_scatter3d_unfiltered_glyph_uses_full_subsample_count(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:

        e = dashboard_experiment
        point_cap = 120_000
        fig = _scatter3d_figure(e, "z0", "z1", "z2", "znorm", color_filter=None)
        tr = fig["data"][0]
        vis_sz, vis_op, n_vis, n_total = _trace_visible_glyph(tr)
        assert n_vis == n_total
        assert not isinstance(tr["marker"]["size"], list)
        exp_sz, exp_op = _expected_dashboard_scatter3d_glyph(
            n_total, point_cap=point_cap
        )
        assert vis_sz == pytest.approx(exp_sz)
        assert vis_op == pytest.approx(exp_op)

    def test_scatter3d_partial_filter_enlarges_glyph_size(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        """Fewer visible particles → larger glyphs; full subsample unchanged."""

        import pandas as pd

        e = dashboard_experiment
        point_cap = 120_000
        unfiltered = _scatter3d_figure(e, "z0", "z1", "z2", "znorm", color_filter=None)
        u_tr = unfiltered["data"][0]
        u_sz, u_op, _, n_sub = _trace_visible_glyph(u_tr)
        exp_sz, exp_op = _expected_dashboard_scatter3d_glyph(n_sub, point_cap=point_cap)
        assert u_sz == pytest.approx(exp_sz)
        assert u_op == pytest.approx(exp_op)

        znorm = pd.to_numeric(e.plot_df["znorm"], errors="coerce")
        range_cf = {
            "kind": "range",
            "range_min": float(znorm.quantile(0.25)),
            "range_max": float(znorm.quantile(0.75)),
        }
        ranged = _scatter3d_figure(e, "z0", "z1", "z2", "znorm", color_filter=range_cf)
        r_sz, r_op, n_vis_range, _ = _trace_visible_glyph(ranged["data"][0])
        assert n_vis_range < n_sub
        assert r_sz > u_sz
        assert r_op > u_op
        assert isinstance(ranged["data"][0]["marker"]["size"], list)

        labels = e.plot_df["labels"].dropna()
        fk = covariate_row_filter_key("labels", labels.iloc[0])
        disc_cf = {"kind": "discrete", "keys": [fk]}
        discrete = _scatter3d_figure(
            e, "z0", "z1", "z2", "labels", color_filter=disc_cf
        )
        d_sz, d_op, n_vis_disc, _ = _trace_visible_glyph(discrete["data"][0])
        assert n_vis_disc < n_vis_range
        assert d_sz > r_sz
        assert d_sz > u_sz

    def test_scatter3d_filter_glyph_scaling_volume_landscape_style(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:

        e = dashboard_experiment
        sub = e.plot_df.iloc[:200].copy()
        allow = frozenset({"z0", "z1", "z2"})
        point_cap = 120_000
        kwargs = dict(
            plot_df=sub,
            xyz_axes_allowed=allow,
            volume_landscape_3d_style=True,
        )
        unfiltered = _scatter3d_figure(
            e, "z0", "z1", "z2", "labels", color_filter=None, **kwargs
        )
        u_sz, u_op, _, n_sub = _trace_visible_glyph(unfiltered["data"][0])
        fk = covariate_row_filter_key("labels", sub["labels"].dropna().iloc[0])
        filtered = _scatter3d_figure(
            e,
            "z0",
            "z1",
            "z2",
            "labels",
            color_filter={"kind": "discrete", "keys": [fk]},
            **kwargs,
        )
        f_sz, f_op, n_vis, _ = _trace_visible_glyph(filtered["data"][0])
        exp_u_sz, exp_u_op = _expected_dashboard_scatter3d_glyph(
            n_sub,
            point_cap=point_cap,
            volume_landscape_3d_style=True,
        )
        assert n_vis < n_sub
        assert u_sz == pytest.approx(exp_u_sz)
        assert u_op == pytest.approx(exp_u_op)
        assert f_sz > u_sz
        assert f_op > u_op

    def test_latent_3d_scatter_filter_preserves_xyz(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:

        import pandas as pd

        e = dashboard_experiment
        znorm = pd.to_numeric(e.plot_df["znorm"], errors="coerce")
        cf = {
            "kind": "range",
            "range_min": float(znorm.quantile(0.25)),
            "range_max": float(znorm.quantile(0.75)),
        }
        unfiltered = _scatter3d_figure(e, "z0", "z1", "z2", "znorm", color_filter=None)
        filtered = _scatter3d_figure(e, "z0", "z1", "z2", "znorm", color_filter=cf)
        _assert_filter_preserves_plotly_xyz(unfiltered, filtered)

        labels = e.plot_df["labels"].dropna()
        fk = covariate_row_filter_key("labels", labels.iloc[0])
        disc_cf = {"kind": "discrete", "keys": [fk]}
        disc = _scatter3d_figure(e, "z0", "z1", "z2", "labels", color_filter=disc_cf)
        _assert_filter_preserves_plotly_xyz(unfiltered, disc)

    def test_volume_landscape_3d_scatter_filter_preserves_xyz(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        """``/api/scatter3d_z_landscape_full`` → same ``scatter3d_z_json`` contract."""

        e = dashboard_experiment
        sub = e.plot_df.iloc[:120].copy()
        allow = frozenset({"z0", "z1", "z2"})
        unfiltered = _scatter3d_figure(
            e,
            "z0",
            "z1",
            "z2",
            "labels",
            plot_df=sub,
            xyz_axes_allowed=allow,
            volume_landscape_3d_style=True,
            color_filter=None,
        )
        fk = covariate_row_filter_key("labels", sub["labels"].dropna().iloc[0])
        filtered = _scatter3d_figure(
            e,
            "z0",
            "z1",
            "z2",
            "labels",
            plot_df=sub,
            xyz_axes_allowed=allow,
            volume_landscape_3d_style=True,
            color_filter={"kind": "discrete", "keys": [fk]},
        )
        _assert_filter_preserves_plotly_xyz(unfiltered, filtered)

    def test_no_subsample_discrete_gif_frame_preserves_xyz(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        """Discrete GIF frame capture (``no_subsample``) uses opacity, not row drop."""

        e = dashboard_experiment
        sub = e.plot_df.iloc[:80]
        fk = covariate_row_filter_key("labels", sub["labels"].dropna().iloc[0])
        unfiltered = _scatter3d_figure(
            e,
            "z0",
            "z1",
            "z2",
            "labels",
            plot_df=sub,
            no_subsample=True,
            color_filter=None,
        )
        filtered = _scatter3d_figure(
            e,
            "z0",
            "z1",
            "z2",
            "labels",
            plot_df=sub,
            no_subsample=True,
            color_filter={"kind": "discrete", "keys": [fk]},
        )
        assert len(unfiltered["data"][0]["x"]) == len(sub)
        _assert_filter_preserves_plotly_xyz(unfiltered, filtered)
        point_cap = max(len(sub), 1)
        u_sz, _, _, n_all = _trace_visible_glyph(unfiltered["data"][0])
        f_sz, _, n_vis, _ = _trace_visible_glyph(filtered["data"][0])
        exp_u, _ = _expected_dashboard_scatter3d_glyph(n_all, point_cap=point_cap)
        assert n_vis < n_all
        assert u_sz == pytest.approx(exp_u)
        assert f_sz > u_sz

    def test_discrete_level_matplotlib_export_subsets_rows_by_design(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        """Static GIF frames are not interactive legend toggles — subset is intentional."""
        e = dashboard_experiment
        sub = e.plot_df.iloc[:200]
        fk = covariate_row_filter_key("labels", sub["labels"].dropna().iloc[0])
        png = scatter3d_discrete_level_png_bytes(
            sub,
            "z0",
            "z1",
            "z2",
            "labels",
            fk,
            discrete_label_colors=None,
            scene_axis_titles=("z0", "z1", "z2"),
            exp=e,
            xyz_axes_allowed=None,
        )
        assert isinstance(png, bytes) and len(png) > 100


class TestScatter3dLatent3dTemplateContracts:
    """Template / vol-PCA JS contracts for the shared ``latent_3d.html`` shell."""

    def test_latent_3d_vol_anim_wires_preview_loading_callback(self) -> None:
        """ChimeraX runs use ``setPreviewAnimLoading`` on the preview stack; main ``setRendering`` is scatter reload only."""
        text = read_latent_3d_html()
        assert "setPreviewAnimLoading:" in text
        assert "latent3dVolAnimApi.afterPlotRedraw" in text
        assert "l3dva-preview-anim-overlay" in text
        assert "function setRendering(on" in text
        assert "setRendering(true" in text
        assert "setLatent3dPreviewAnimLoading" in text
        assert "CryoLatent3dVolLandscapeAnim.boot" in text
        vol_text = read_dashboard_static_js("latent3d_landscape_vol_animations.js")
        assert "reapplySelectionHighlightIfNeeded" in vol_text
        assert "estimatedChimeraxViewMatrixText" in vol_text
        assert "viewRotationsAreActive" in vol_text
        assert "syncViewMatrixField" in vol_text
        assert "applyViewMatrixFromField" in vol_text
        assert "view-matrix-input" in text
        assert 'name="l3dva-gif-mode"' in text
        assert 'value="disabled"' in text
        assert "animationsEnabled" in vol_text
        assert "volSelectionBlockedForAdd" in vol_text
        assert "ChimeraX view matrix: rendering" not in vol_text
        assert "referenceScatter3dBaseMarkerSize" in vol_text

    def test_landscape_volpca_selection_sizes_respect_trace_marker(self) -> None:
        vol_text = read_dashboard_static_js("landscape_volpca.js")
        assert "referenceScatter3dBaseMarkerSize" in vol_text
        assert "cdrgnVolSelectionOverlay" in vol_text
        assert "volIdToPointIndex" in vol_text
        assert "updateVolSelectionOverlay" in vol_text
        assert "volMontageLabel" in vol_text
        assert "syncStableSelectionLabels" in vol_text
        assert "suppressPlotlySelectedUntil" in vol_text
        assert "beginClickToggleSelection" in vol_text
        assert "plotlySelectedSuppressed" in vol_text
        assert "plotlySelectedFromLassoOrBox" in vol_text
        assert "volFromClickEvent" in vol_text
        assert "volLassoBoxGestureActive" in vol_text
        assert "plotly_selecting" in vol_text
        assert "baseTracePointsFromSelectedEvent" in vol_text
        assert "traceForPlotlyPoint" in vol_text
        assert "applyBaseTraceOpacityDimming" in vol_text
        assert "VOLSKETCH_LASSO_DEBOUNCE_MS" in vol_text
        assert "beginVolsketchLassoDimmingGesture" in vol_text
        start = vol_text.find('gd.on("plotly_selected", function(ev)')
        assert start != -1
        mid = vol_text.find("volsketchLassoSelectTimer = setTimeout(function()", start)
        assert mid != -1
        end = vol_text.find("}, VOLSKETCH_LASSO_DEBOUNCE_MS);", mid)
        debounce_block = vol_text[mid:end]
        assert "handleVolSketchPlotlySelected" in debounce_block
        assert "selectedpoints: [indices.slice()]" in vol_text
        assert "assignVolMontageLabelsSorted" in vol_text
        assert "resetVolMontageLabels" in vol_text

    def test_latent_3d_loads_scatter3d_scene_module_before_vol_anim(self) -> None:
        text = read_latent_3d_html()
        assert "plotly_scatter3d_scene.js" in text
        assert "lastLatent3dPersistedSceneSnap" in text
        assert "latent3dSceneSnapHasCamera" in text
        assert "inset: 0" in text
        assert "background: transparent" in text
        p3 = text.index("plotly_scatter3d_scene.js")
        assert p3 < text.index("(function() {")
        js_text = read_dashboard_static_js("plotly_scatter3d_scene.js")
        assert "CryoPlotlyScatter3dScene" in js_text
        assert "applySnapshotToFigureLayout" in js_text
        assert "sceneRelayoutPatchCameraOnly" in js_text
        # Colour-covariate changes await legend ``refresh`` before vol-anim tail + scene restore.
        assert "pRefresh" in text and "cRefresh.then" in text


class TestPlotlyScatter3dSceneCameraPreserve:
    """Static contracts for ``plotly_scatter3d_scene.js`` orbit preservation."""

    def test_scene_relayout_patch_camera_only_is_camera_only(
        self, plotly_scatter3d_scene_js: str
    ) -> None:
        body = js_function_body(
            plotly_scatter3d_scene_js,
            "function sceneRelayoutPatchCameraOnly",
            "function restore(",
        )
        assert '"scene.camera"' in body
        assert "scene.xaxis.range" not in body
        assert "scene.yaxis.range" not in body

    def test_marker_restyle_preserving_camera_uses_plotly_update(
        self, plotly_scatter3d_scene_js: str
    ) -> None:
        body = js_function_body(
            plotly_scatter3d_scene_js,
            "function traceMarkerRestylePreservingCamera",
            "function traceCoordsRestyleFromFigure",
        )
        assert "Plotly.update(gd, upd, layoutPatch, [0])" in body
        assert "snapHasAxisRanges(pinSnap)" in body
        assert "sceneRelayoutPatchFromSnap(pinSnap)" in body
        assert "sceneRelayoutPatchCameraOnly(pinSnap)" in body

    def test_relayout_preserve_view_patch_uses_pin_ranges_not_server_layout(
        self, plotly_scatter3d_scene_js: str
    ) -> None:
        body = js_function_body(
            plotly_scatter3d_scene_js,
            "function relayoutPreserveViewPatch",
            "function applySnapshotToFigureLayout",
        )
        assert "sceneRelayoutPatchFromSnap(pin)" in body
        assert "sceneAxisRangeRelayoutPatchFromLayout" not in body

    def test_resolve_camera_prefers_webgl_get_camera(
        self, plotly_scatter3d_scene_js: str
    ) -> None:
        body = js_function_body(
            plotly_scatter3d_scene_js,
            "function resolveCameraFromGd",
            "function snapshot(",
        )
        assert "scene._scene.getCamera" in body

    def test_apply_snapshot_camera_only_skips_axis_ranges(
        self, plotly_scatter3d_scene_js: str
    ) -> None:
        body = js_function_body(
            plotly_scatter3d_scene_js,
            "function applySnapshotCameraOnlyToFigureLayout",
            "function layoutPatchWithoutScene",
        )
        assert "fig.layout.scene.camera" in body
        assert "xaxis" not in body or "applyRange" not in body


class TestColorCovariateLegendFilterColumnBinding:
    """Legend filters must not leak across covariate column changes."""

    def test_get_filter_for_api_discrete_requires_matching_column(self) -> None:
        text = read_dashboard_static_js("color_covariate_legend.js")
        body = js_function_body(
            text,
            "CryoColorCovariateLegend.prototype.getFilterForApi",
            "CryoColorCovariateLegend.prototype.clearThreshold",
        )
        assert "_discreteCol !== col" in body
        refresh = text.split('if (data.mode === "discrete")', 1)[1]
        assert "self._discreteCol = col" in refresh


class TestScatter3dLatent3dCameraSnapBack:
    """Regression guards for latent-3D viewing-angle snap-back fixes."""

    def test_set_rendering_does_not_restore_camera_before_overlay(self) -> None:
        text = read_latent_3d_html()
        fn = text.split("function setRendering(on, opts)", 1)[1]
        body = fn.split("function selectedLatentPalette", 1)[0]
        assert "Do not relayout the camera before showing the badge" in body
        assert "restoreCameraOnly" not in body

    def test_same_axes_reload_uses_nonblocking_rendering_overlay(self) -> None:
        text = read_latent_3d_html()
        fn = text.split("function loadPlot()", 1)[1]
        body = fn.split("function buildLatent3dPostPayload", 1)[0]
        assert "useNonblockingOverlay" in body
        assert "latent3dSameAxesAsLast" in body
        assert "setRendering(true, { nonblocking: useNonblockingOverlay })" in body

    def test_load_plot_captures_live_scene_pin_before_fetch(self) -> None:
        text = read_latent_3d_html()
        fn = text.split("function loadPlot()", 1)[1]
        body = fn.split("fetch(scatter3dUrl", 1)[0]
        assert "latent3dCaptureScenePinPreferLive()" in body
        assert (
            "P3S.snapshot(gd)"
            in text.split("function latent3dCaptureScenePinPreferLive", 1)[1].split(
                "function latent3dCameraWatchdog", 1
            )[0]
        )

    def test_colour_covariate_tail_awaits_legend_refresh_before_overlay_off(
        self,
    ) -> None:
        text = read_latent_3d_html()
        assert "knocks the scatter3d camera back" in text
        tail = text.split("function runLatent3dPostResizeTail", 1)[1]
        assert "cRefresh.then(function()" in tail
        idx_refresh = tail.index("cRefresh.then")
        idx_overlay_off = tail.index("setRendering(false)")
        assert idx_refresh < idx_overlay_off

    def test_same_axes_marker_restyle_preserves_camera(self) -> None:
        text = read_latent_3d_html()
        fn = text.split("function latent3dApplyFigurePreservingScene", 1)[1]
        body = fn.split("function latent3dStabilizeSceneAfterOverlay", 1)[0]
        assert "traceMarkerRestylePreservingCamera(fig, gd, viewPin)" in body
        assert "prepareFigureSceneForPinnedViewRedraw(fig, gd, pin)" in body
        assert "relayoutPreserveViewPatch(viewPin)" in body

    def test_same_axes_filter_restyle_does_not_apply_server_axis_ranges(self) -> None:
        text = read_latent_3d_html()
        fn = text.split("function latent3dApplyFigurePreservingScene", 1)[1]
        body = fn.split("function latent3dStabilizeSceneAfterOverlay", 1)[0]
        assert "relayoutCameraAndAxisPatch(pin, fig.layout, true)" not in body

    def test_legend_colour_filter_preserves_user_3d_axis_limits_regression(
        self, plotly_scatter3d_scene_js: str
    ) -> None:
        """Client-side freeze on ``latent_3d.html`` (server contract: ``TestScatter3dLegendFilterAxisStability``)."""
        html = read_latent_3d_html()
        scene_js = plotly_scatter3d_scene_js
        apply_body = html.split("function latent3dApplyFigurePreservingScene", 1)[
            1
        ].split("function latent3dStabilizeSceneAfterOverlay", 1)[0]
        assert "prepareFigureSceneForPinnedViewRedraw(fig, gd, pin)" in apply_body
        assert "relayoutPreserveViewPatch(viewPin)" in apply_body
        assert "relayoutCameraAndAxisPatch(pin, fig.layout, true)" not in apply_body
        marker_fn = scene_js.split("function traceMarkerRestylePreservingCamera", 1)[1]
        assert (
            "snapHasAxisRanges(pinSnap)"
            in marker_fn.split("function traceCoordsRestyleFromFigure", 1)[0]
        )
        prep_fn = scene_js.split("function prepareFigureSceneForPinnedViewRedraw", 1)[1]
        prep_body = prep_fn.split("function restorePinnedView", 1)[0]
        assert "stripSceneAxisRangesFromFigureLayout(fig)" in prep_body
        assert "effectiveViewSnapForRedraw(gd, pin)" in prep_body
        enforce_fn = scene_js.split("function scheduleEnforceCamera", 1)[1]
        enforce_body = enforce_fn.split("function markerRestyleUpdateFromFigure", 1)[0]
        assert "restorePinnedView(gd, pin)" in enforce_body

    def test_rendering_overlay_is_not_stacked_over_webgl_plot(self) -> None:
        text = read_latent_3d_html()
        assert "latent3d-rendering-overlay" in text
        assert "cryo-plot-rendering-overlay--nonblocking" in text
        # Full-bleed veil must not use backdrop-filter over the WebGL canvas.
        head = text.split("{% block content %}", 1)[0]
        assert "backdrop-filter" not in head or "No backdrop-filter" in head
        assert "#latent3d-rendering-overlay" in head
        assert "inset: 0" in head

    def test_afterplot_does_not_restore_stale_load_pin_after_user_orbit(self) -> None:
        text = read_latent_3d_html()
        fn = text.split("function latent3dOnPlotlyAfterplot", 1)[1]
        body = fn.split("function latent3dPersistSceneAfterUserOrbit", 1)[0]
        assert "User orbited during scene hold" in body
        assert "cameraEyeDiffers(loadPin, liveAfter)" in body
        persist = text.split("function latent3dPersistSceneAfterUserOrbit", 1)[1]
        persist = persist.split("function latent3dOnPlotlyAfterplot", 1)[0]
        assert "gd._cryoLatent3dLoadScenePin = holdPin" in persist

    def test_volanim_discrete_legend_does_not_expand_from_plot_width(self) -> None:
        """Discrete mode must not set ``--cryo-discrete-legend-w`` from plot stack (squeezes the 3D view)."""
        text = read_latent_3d_html()
        fn = text.split("function syncLatent3dDiscreteLegendTargetWidth", 1)[1]
        body = fn.split("function syncLatent3dLegendMiddleWidth", 1)[0]
        assert "latent3dIsVolanimLayout()" in body
        assert 'removeProperty("--cryo-discrete-legend-w")' in body
        volanim = text.split(
            ".cryo-dash-row--latent3d-volanim .latent3d-plot-legend-band", 1
        )[1]
        volanim = volanim.split("{% endif %}", 1)[0]
        assert "grid-template-columns: minmax(0, 1fr) minmax(0, 9.25rem)" in volanim
        assert "cryo-cc-discrete-switches--pair-2col" in volanim


class TestLatent3dVolAnimScatter3dCameraPreserve:
    """Vol-landscape animation UI must not reset the main scatter3d orbit on relayout."""

    def test_montage_relayout_only_patches_camera_when_annotations_present(
        self,
    ) -> None:
        text = read_dashboard_static_js("latent3d_landscape_vol_animations.js")
        fn = text.split("function relayoutMontageAnnotations", 1)[1]
        body = fn.split("function refreshSelectionHighlight", 1)[0]
        assert "captureVolAnimScenePinLiveOnly() || scenePin" in body
        assert "if (anns.length && camPin && camPin.camera" in body
        assert "Patching scene.camera when only clearing" in text

    def test_capture_vol_anim_scene_pin_prefers_live_snapshot(self) -> None:
        text = read_dashboard_static_js("latent3d_landscape_vol_animations.js")
        fn = text.split("function captureVolAnimScenePin", 1)[1]
        body = fn.split("function scheduleRestoreVolAnimScenePin", 1)[0]
        assert "P3S.snapshot(gd)" in body
        assert "P3S.getPinnedSnap(gd)" in body

    def test_selection_highlight_extends_scene_hold_when_load_pin_set(self) -> None:
        text = read_dashboard_static_js("latent3d_landscape_vol_animations.js")
        fn = text.split("function refreshSelectionHighlight", 1)[1]
        body = fn.split("function reapplySelectionHighlightIfNeeded", 1)[0]
        assert "gd._cryoLatent3dLoadScenePin" in body
        assert "gd._cryoLatent3dSceneHoldUntil = Date.now()" in body

    def test_vol_anim_restore_uses_camera_only_api(self) -> None:
        text = read_dashboard_static_js("latent3d_landscape_vol_animations.js")
        assert "P3S.restoreCameraOnly" in text
        assert "P3S.scheduleEnforceCamera" in text

    def test_cycle_segment_colors_prefer_server_backgrounds(self) -> None:
        text = read_dashboard_static_js("latent3d_landscape_vol_animations.js")
        fn = text.split("function refreshPreviewOverlayLetterColors", 1)[1]
        body = fn.split("function boldAnchorPointForVol", 1)[0]
        assert "segment_covariate_backgrounds" in body
        assert "fromServer" in body

    def test_refresh_selection_highlight_prefers_live_scene_pin(self) -> None:
        """GIF / montage tail must not restore double-click / load pin over the live orbit."""
        text = read_dashboard_static_js("latent3d_landscape_vol_animations.js")
        fn = text.split("function refreshSelectionHighlight", 1)[1]
        body = fn.split("function toggleVol", 1)[0]
        assert "_cryoLatent3dLoadScenePin || captureVolAnimScenePin" not in body
        assert "var scenePin = captureVolAnimScenePinLiveOnly()" in body
        assert "captureVolAnimScenePinLiveOnly() || scenePin" in text
        assert "volAnimSelOverlayDepth >= 1" in body
        assert "mergeSnapsPreferPointerdown(volAnimPointerdownSnap, pin)" in text
        cap = text.split("function captureVolAnimScenePin()", 1)[1]
        cap = cap.split("function scheduleRestoreVolAnimScenePin", 1)[0]
        assert "Date.now() - volAnimPointerdownAt < 1200" in cap

    def test_montage_letter_colors_use_discrete_and_continuous_plot_modes(self) -> None:
        text = read_dashboard_static_js("latent3d_landscape_vol_animations.js")
        assert "function discreteColorLegendMap" in text
        assert "function continuousCssFromValue" in text
        assert "montageLetterFontColor(trace, gd, i)" in text
        assert "montagePlotLetterStrokeForFill" in text
        assert "rebuildRotateFrameBadgeBackgrounds" in text
        assert "applyBadgeLetterColor(badge, letterBg)" in text
        assert "#ffffff" in text or "#ffffff" in text
        assert "applyBadgeLetterColor" in text
        assert "sketchCentroidPointIndexForVol" in text
        assert "VOL_MONTAGE_PLOT_LETTER_PX_SECONDARY" in text
        assert "montagePlotLetterFontSizePx" in text
        assert "VOL_MONTAGE_PLOT_LETTER_PX * 0.8" in text
        assert "resolvePlotPalette()" in text
        assert "refreshPreviewOverlayLetterColors" in text
