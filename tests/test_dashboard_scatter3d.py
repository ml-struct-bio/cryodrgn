"""Dashboard 3-D scatterplot (Plotly ``scatter3d``) API, plot JSON, and client contracts.

Covers ``scatter3d_z_json``, ``/api/scatter3d_z*``, latent-3D preview/GIF routes, and
``plotly_scatter3d_scene.js`` / ``latent_3d.html`` orbit-preservation behaviour.
"""

from __future__ import annotations

import base64
import json
import re
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
    _scatter3d_marker_sizes_with_legend_filter,
    scatter3d_discrete_level_png_bytes,
    scatter3d_z_json,
)
from cryodrgn.dashboard.plot_gif_utils import png_base64_frames_to_gif_bytes
from tests.conftest import (
    DASHBOARD_ANALYZE_EPOCH,
    decode_plotly_figure,
    decode_plotly_value,
    plotly_trace_array,
    png_b64_rgb,
)

pytestmark = pytest.mark.dashboard

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
    return decode_plotly_figure(
        json.loads(scatter3d_z_json(exp, x, y, z, color, **kwargs))
    )


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
    n_total = len(plotly_trace_array(trace, "x"))
    op = float(trace["marker"]["opacity"])
    sz = decode_plotly_value(trace["marker"]["size"])
    if isinstance(sz, list):
        vis = [float(s) for s in sz if float(s) > 0.0]
        assert vis, "expected at least one visible point"
        assert all(s == vis[0] for s in vis), "visible points share one glyph size"
        return vis[0], op, len(vis), n_total
    return float(sz), op, n_total, n_total


def _assert_filter_preserves_plotly_xyz(unfiltered: dict, filtered: dict) -> None:
    u = unfiltered["data"][0]
    f = filtered["data"][0]
    ux = plotly_trace_array(u, "x")
    fx = plotly_trace_array(f, "x")
    assert len(fx) == len(ux)
    assert fx == ux
    assert plotly_trace_array(f, "y") == plotly_trace_array(u, "y")
    assert plotly_trace_array(f, "z") == plotly_trace_array(u, "z")
    for axis_name in ("xaxis", "yaxis", "zaxis"):
        scene = filtered["layout"]["scene"]
        assert "range" not in (scene.get(axis_name) or {})
    sizes = decode_plotly_value(f["marker"]["size"])
    assert isinstance(sizes, list)
    assert len(sizes) == len(fx)
    assert any(float(s) == 0.0 for s in sizes)
    assert any(float(s) > 0.0 for s in sizes)
    colors = decode_plotly_value(f["marker"]["color"])
    assert isinstance(colors, list)
    assert not any(str(c).strip().lower().startswith("rgba") for c in colors[:50])
    op = f["marker"]["opacity"]
    assert isinstance(op, (int, float))
    vis_sizes = [float(s) for s in sizes if float(s) > 0.0]
    u_sz = decode_plotly_value(u["marker"]["size"])
    u_base = float(u_sz) if not isinstance(u_sz, list) else float(max(u_sz))
    if vis_sizes and len(vis_sizes) == len(sizes):
        assert max(vis_sizes) == pytest.approx(u_base)
    elif vis_sizes and len(vis_sizes) < len(sizes):
        assert max(vis_sizes) > u_base


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

    def test_api_scatter3d_z_landscape_full_with_mock_outputs(
        self,
        flask_client_landscape_full,
        dashboard_workdir_with_landscape_full: str,
    ) -> None:
        from cryodrgn.dashboard.landscape_full_3d import landscape_full_3d_ready

        assert landscape_full_3d_ready(
            dashboard_workdir_with_landscape_full, DASHBOARD_ANALYZE_EPOCH
        )
        r = flask_client_landscape_full.get(
            "/api/scatter3d_z_landscape_full"
            "?x=landscape_vol_PC1&y=landscape_vol_PC2&z=landscape_vol_PC3&color=none"
        )
        assert r.status_code == 200, r.get_data(as_text=True)[:500]
        assert r.get_json()["data"]

    def test_landscape_full_scatter_covariate_roundtrip(
        self, flask_client_landscape_full
    ) -> None:
        """Colour by vol-PC uses display labels from ``covariate_labels``."""
        from cryodrgn.dashboard import covariate_labels

        r = flask_client_landscape_full.get(
            "/api/scatter3d_z_landscape_full"
            "?x=landscape_vol_PC1&y=landscape_vol_PC2&z=landscape_vol_PC3"
            "&color=landscape_vol_PC1"
        )
        assert r.status_code == 200, r.get_data(as_text=True)[:500]
        fig = r.get_json()
        assert fig["data"]
        layout = fig.get("layout") or {}
        # Axis / colourbar titles should use the pretty Vol PC label when available.
        pretty = covariate_labels.covariate_display_name("landscape_vol_PC1")
        text_blob = json.dumps(layout)
        assert (
            "Vol PC" in text_blob
            or pretty in text_blob
            or "landscape_vol_PC1" in text_blob
        )
        marker = fig["data"][0].get("marker") or {}
        assert marker.get("color") is not None or marker.get("colorscale") is not None

    def test_landscape_full_3d_page_renders_with_mock_outputs(
        self, flask_client_landscape_full
    ) -> None:
        r = flask_client_landscape_full.get("/landscape-full-3d")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert 'id="latent3d"' in body
        assert "/vendor/plotly.min.js" in body
        assert "analyze_landscape_full" not in body.lower()

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

    def test_scatter3d_vol_landscape_centroid_cd_before_nearest_with_color(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        """Colour covariates must not displace the centroid flag from length-2.

        ``latent3d_landscape_vol_animations.js`` circles montage letters using
        ``customdata[length-2] === 1``. Packing the colour column after the
        centroid flag made circling depend on colour values.
        """
        e = dashboard_experiment
        sub = e.plot_df.iloc[:40].copy()
        sub[VOL_LANDSCAPE_3D_PLOT_DF_ROW] = np.arange(40, dtype=np.int64)
        sub[VOL_LANDSCAPE_NEAREST_SKETCH_VOL] = (
            np.arange(len(sub), dtype=np.int64) % 3 + 1
        )
        cent = (np.arange(len(sub), dtype=np.int64) % 5 == 0).astype(np.int64)
        sub[VOL_LANDSCAPE_IS_SKETCH_CENTROID] = cent
        # Distinct from 0/1 so a misplaced colour column cannot look like a flag.
        sub["umap1"] = np.linspace(10.0, 20.0, len(sub), dtype=np.float64)
        allow = frozenset(list(LATENT_Z_AXES) + ["umap1"])
        fig = _scatter3d_figure(
            e,
            "z0",
            "z1",
            "z2",
            "umap1",
            plot_df=sub,
            xyz_axes_allowed=allow,
            volume_landscape_3d_style=True,
        )
        cd = np.asarray(fig["data"][0]["customdata"])
        assert cd.shape[1] == 4
        # [plot_row, colour, centroid_flag, nearest_vol] — colour must not sit at length-2.
        assert np.all((cd[:, 1] >= 10.0) & (cd[:, 1] <= 20.0))
        assert set(cd[:, -2].astype(np.int64).tolist()) <= {0, 1}
        assert int(cd[:, -2].sum()) == int(cent.sum())
        assert set(cd[:, -1].astype(np.int64).tolist()) <= {1, 2, 3}
        assert not np.allclose(cd[:, -2], cd[:, 1])

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
        assert plotly_trace_array(m, "x") == plotly_trace_array(u, "x")
        assert m["marker"]["size"] == u["marker"]["size"]
        assert m["marker"]["opacity"] == u["marker"]["opacity"]

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
        assert len(plotly_trace_array(unfiltered["data"][0], "x")) == len(sub)
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


class TestScatter3dBrowserSmoke:
    """Headless Chromium: latent-3D and landscape-full-3D Plotly flows."""

    pytestmark = pytest.mark.browser

    def test_scatter3d_and_discrete_legend(
        self, playwright_page, dashboard_live_url
    ) -> None:
        from tests.conftest import dashboard_smoke_latent_3d

        out = dashboard_smoke_latent_3d(playwright_page, dashboard_live_url)
        assert out["trace0type"] == "scatter3d"
        assert out["discrete_legend_toggles"] >= 1

    def test_camera_stable_on_discrete_covariate_change(
        self, playwright_page, dashboard_live_url
    ) -> None:
        from tests.conftest import dashboard_smoke_latent_3d_camera_on_covariate_change

        out = dashboard_smoke_latent_3d_camera_on_covariate_change(
            playwright_page, dashboard_live_url
        )
        assert out["camera_stable"]

    @pytest.mark.parametrize(
        "smoke_name,fixture_url,key,min_val",
        [
            (
                "dashboard_smoke_landscape_full_3d",
                "dashboard_landscape_full_live_url",
                "count",
                1,
            ),
            (
                "dashboard_smoke_landscape_full_3d_clear_selection",
                "dashboard_landscape_full_live_url",
                "annotations_cleared",
                True,
            ),
        ],
        ids=["random_selection", "clear_selection"],
    )
    def test_landscape_full_3d_smokes(
        self,
        smoke_name: str,
        fixture_url: str,
        key: str,
        min_val,
        request,
        playwright_page,
    ) -> None:
        import tests.conftest as cf

        smoke_fn = getattr(cf, smoke_name)
        url = request.getfixturevalue(fixture_url)
        out = smoke_fn(playwright_page, url)
        assert out is not None
        if isinstance(min_val, bool):
            assert out[key] is min_val
        else:
            assert out[key] >= min_val


@pytest.fixture(scope="class")
def oriented_page(playwright_isolated_browser, dashboard_live_url):
    """One orbited latent-3D page shared by ``TestLatent3dCameraPreservation``."""
    from tests.conftest import (
        dashboard_open_latent_3d,
        dashboard_orbit_scene_camera,
        dashboard_scene_camera,
    )

    context = playwright_isolated_browser.new_context(
        viewport={"width": 1400, "height": 900}
    )
    page = context.new_page()
    dashboard_open_latent_3d(page, dashboard_live_url)
    default_camera = dashboard_scene_camera(page, "latent3d")
    orbited = dashboard_orbit_scene_camera(page, "latent3d")
    yield page, default_camera, orbited
    context.close()


class TestLatent3dCameraPreservation:
    """A user's orbit must survive every control that redraws the 3D scene.

    Losing the camera on redraw has been fixed repeatedly: a selection restyled the
    trace ``mode`` and reset orbit and zoom, axis ranges read back as ``undefined``,
    and Plotly's own reset fired before the pose could be captured. Comparing cameras
    from the default pose cannot detect a reset *back* to that default, so the scene
    is dragged somewhere distinctive first and every check also asserts the pose has
    not fallen back.

    One page serves the whole class: each 3D scene holds WebGL resources that are not
    released between contexts, and after a few of them a drag no longer orbits at all
    under SwiftShader.
    """

    pytestmark = pytest.mark.browser

    @staticmethod
    def _settled_camera(page):
        from tests.conftest import (
            DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
            _dashboard_smoke_wait_latent3d_overlay_hidden,
            dashboard_scene_camera,
        )

        _dashboard_smoke_wait_latent3d_overlay_hidden(
            page, timeout_ms=DASHBOARD_BROWSER_FAST_TIMEOUT_MS
        )
        return dashboard_scene_camera(page, "latent3d")

    def _assert_camera_held(self, oriented_page, before, after):
        """The pose must be unchanged, and must not have fallen back to the default."""
        from tests.conftest import dashboard_cameras_match

        _, default_camera, _ = oriented_page
        assert dashboard_cameras_match(before, after), f"{before!r} -> {after!r}"
        assert not dashboard_cameras_match(
            default_camera, after
        ), "camera reset to the default pose"

    def _to_discrete_colour(self, page):
        from tests.conftest import (
            DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
            _dashboard_smoke_set_select_value,
        )

        _dashboard_smoke_set_select_value(page, "sc", "labels")
        page.wait_for_function(
            """() => {
              var s = document.getElementById('latent3d-color-discrete-switches');
              return s && s.querySelectorAll('button, label, input').length > 0;
            }""",
            timeout=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        return self._settled_camera(page)

    def test_dragging_the_scene_orbits_the_camera(self, oriented_page):
        """Guards the rest of the class: a no-op orbit would make it vacuous."""
        from tests.conftest import dashboard_cameras_match

        _, default_camera, orbited = oriented_page
        assert orbited["eye"] and orbited["up"]
        assert not dashboard_cameras_match(default_camera, orbited)

    def test_orbit_survives_a_discrete_covariate_change(self, oriented_page):
        page, _, _ = oriented_page
        before = self._settled_camera(page)
        after = self._to_discrete_colour(page)
        self._assert_camera_held(oriented_page, before, after)

    def test_orbit_survives_toggling_a_discrete_level(self, oriented_page):
        page, _, _ = oriented_page
        before = self._to_discrete_colour(page)

        page.locator("#latent3d-color-discrete-switches input").first.click()
        page.wait_for_timeout(600)
        self._assert_camera_held(oriented_page, before, self._settled_camera(page))

    def test_orbit_survives_inverting_the_discrete_selection(self, oriented_page):
        page, _, _ = oriented_page
        before = self._to_discrete_colour(page)

        invert = page.locator("#latent3d-btn-discrete-invert")
        if invert.count() == 0 or not invert.is_visible():
            pytest.skip("discrete invert control not offered on this dataset")
        invert.click()
        page.wait_for_timeout(600)
        self._assert_camera_held(oriented_page, before, self._settled_camera(page))

    def test_orbit_survives_clearing_the_colour_column(self, oriented_page):
        from tests.conftest import _dashboard_smoke_set_select_value

        page, _, _ = oriented_page
        before = self._to_discrete_colour(page)

        _dashboard_smoke_set_select_value(page, "sc", "none")
        page.wait_for_timeout(600)
        self._assert_camera_held(oriented_page, before, self._settled_camera(page))

    def test_scene_still_holds_points_after_the_redraws(self, oriented_page):
        """A preserved camera pointing at an empty scene would still be a regression."""
        from tests.conftest import (
            DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
            _dashboard_smoke_wait_plot_ready,
        )

        page, _, _ = oriented_page
        self._to_discrete_colour(page)
        info = _dashboard_smoke_wait_plot_ready(
            page, "latent3d", timeout_ms=DASHBOARD_BROWSER_FAST_TIMEOUT_MS
        )
        assert info["n"] > 0
        assert info["type"] == "scatter3d"


class TestLandscapeFull3dColumnConsistency:
    """The 3D volume-landscape page must not offer columns its API rejects.

    The page injects ``axis_cols`` and ``numeric_cols`` into the same
    ``latent_3d.html`` menus that ``/api/scatter3d_z_landscape_full`` validates
    separately. This is the sibling of the landscape vol-PCA colour menu, where the
    two lists had drifted apart and every mixed-case column was refused.
    """

    @staticmethod
    def _page_columns(client):
        """The axis and colour lists the page actually hands to its menus."""
        body = client.get("/landscape-full-3d").get_data(as_text=True)
        axes = re.search(r"var axisCols = (\[.*?\]);", body, re.S)
        colours = re.search(r"var cols = (\[.*?\]);", body, re.S)
        assert axes and colours, "landscape-full-3D page did not embed its menus"
        return json.loads(axes.group(1)), json.loads(colours.group(1))

    @staticmethod
    def _scatter(client, axes, color):
        return client.get(
            "/api/scatter3d_z_landscape_full",
            query_string={"x": axes[0], "y": axes[1], "z": axes[2], "color": color},
        )

    def test_every_advertised_axis_can_be_plotted(
        self, flask_client_landscape_full
    ) -> None:
        axes, _ = self._page_columns(flask_client_landscape_full)
        assert len(axes) >= 3, "the 3D landscape page needs three axes"

        for axis in axes:
            others = [a for a in axes if a != axis][:2]
            r = self._scatter(flask_client_landscape_full, [axis, *others], "none")
            assert r.status_code == 200, (
                f"page offers axis {axis!r} but the scatter rejects it: "
                f"{r.get_data(as_text=True)[:200]}"
            )

    def test_every_advertised_colour_column_is_accepted(
        self, flask_client_landscape_full
    ) -> None:
        axes, colours = self._page_columns(flask_client_landscape_full)
        assert colours, "the 3D landscape page offered no colour columns"

        for column in ["none", *colours]:
            r = self._scatter(flask_client_landscape_full, axes, column)
            assert r.status_code == 200, (
                f"page offers colour {column!r} but the scatter rejects it: "
                f"{r.get_data(as_text=True)[:200]}"
            )

    def test_advertised_colours_actually_colour_the_points(
        self, flask_client_landscape_full
    ) -> None:
        """Accepting a colour and then ignoring it is the quieter half of this bug."""
        from tests.conftest import decode_plotly_figure

        axes, colours = self._page_columns(flask_client_landscape_full)
        for column in colours:
            fig = decode_plotly_figure(
                self._scatter(flask_client_landscape_full, axes, column).get_json()
            )
            marker = fig["data"][0].get("marker") or {}
            assert isinstance(
                marker.get("color"), list
            ), f"colour {column!r} was accepted but left a single flat marker colour"
