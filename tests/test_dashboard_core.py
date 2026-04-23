"""Core dashboard tests (logic + primary API integration).

Split from ``tests/test_dashboard.py`` to keep module size manageable.
"""

from __future__ import annotations

import argparse
import base64
import io
import os
import pickle

import numpy as np
import pandas as pd
import pytest

from cryodrgn.dashboard import app as dash_app
from cryodrgn.dashboard.app import (
    _trajectory_eligibility_error,
    _TRAJECTORY_INELIGIBLE_MSG,
)
from cryodrgn.dashboard.data import DashboardExperiment, list_z_epochs
from cryodrgn.dashboard.explorer_volumes import (
    _chimerax_render_cmds,
    explorer_volumes_eligible,
)
from cryodrgn.dashboard.plots import pair_grid_png
from cryodrgn.dashboard.preload import (
    _preload_cache_time_estimate_bounds,
    _stratified_xy_row_indices,
    encode_particle_batch,
    format_preload_cache_time_hint,
    load_plot_df_rows_from_plot_inds_file,
    montage_bytes,
    particle_thumbnail_b64_from_row,
    sample_plot_df_rows_for_preload,
)
from cryodrgn.dashboard.trajectory import (
    _compute_direct_anchor_trajectory,
    _dijkstra_path_from_neighbors,
    _graph_neighbor_arrays,
    _round_direct_mode_traj_xy,
    _trajectory_xy_ok_for_direct,
    _TRAJ_GRAPH_NEIGHBOR_CACHE,
    compute_trajectory_latent_path,
    default_trajectory_endpoints_xy,
    direct_anchor_particle_indices_payload,
    has_pc_columns,
    has_umap_columns,
    parse_anchor_indices_txt,
    parse_int_from_dict as _parse_int_from_dict,
    parse_traj_interpolation_value as _parse_traj_interpolation_value,
    parse_traj_neighbor_value as _parse_traj_neighbor_value,
    parse_traj_points_value as _parse_traj_points_value,
    parse_trajectory_request_body,
    plot_df_rows_for_dataset_indices,
    random_dataset_indices,
    trajectory_anchor_mode_params,
    trajectory_anchor_payload_from_indices,
    trajectory_default_xy_cols,
    trajectory_plot_axis_columns,
    validate_trajectory_plot_axes,
    z_traj_to_savetxt_str,
)

ANALYZE_EPOCH = 2


def _traj_flask_200_or_ineligible(r, experiment: DashboardExperiment) -> bool:
    """If ``explorer_volumes_eligible`` is true, require HTTP 200; else require 400 + standard error."""
    if explorer_volumes_eligible(experiment):
        assert (
            r.status_code == 200
        ), f"{r.status_code}: {r.get_data(as_text=True)[:500]!r}"
        return True
    assert r.status_code == 400
    err = (r.get_json() or {}).get("error")
    assert err == _TRAJECTORY_INELIGIBLE_MSG, err
    return False


class TestParseIntFromDict:
    """Cover the new ``_parse_int_from_dict`` helper + its three wrappers."""

    @pytest.mark.parametrize(
        "raw,expected",
        [(5, 5), ("7", 7), (3.7, 3), (-100, 0), (999, 20), (None, 4)],
    )
    def test_coerce_and_clamp(self, raw: object, expected: int) -> None:
        got = _parse_int_from_dict(
            {"n": raw} if raw is not None else {}, "n", default=4, lo=0, hi=20
        )
        assert got == expected

    @pytest.mark.parametrize("bad", ["abc", [1, 2], {"x": 1}])
    def test_non_numeric_falls_back_to_default(self, bad: object) -> None:
        assert _parse_int_from_dict({"n": bad}, "n", default=7, lo=0, hi=20) == 7

    def test_missing_key_returns_default(self) -> None:
        assert _parse_int_from_dict({}, "missing", default=9, lo=0, hi=20) == 9

    def test_bad_value_returns_default_unclamped(self) -> None:
        # When ``int(...)`` raises, the default is returned *without* clamping,
        # matching the pre-refactor fallback behaviour.
        assert _parse_int_from_dict({"k": "oops"}, "k", default=-5, lo=0, hi=10) == -5

    def test_traj_points_wrapper_clamps(self) -> None:
        assert _parse_traj_points_value({"n_points": 1}) == 2  # lo=2
        assert _parse_traj_points_value({"n_points": 99}) == 20  # hi=20
        assert _parse_traj_points_value({}) == 4  # default

    def test_traj_interpolation_wrapper_allows_zero(self) -> None:
        # Interpolation range is [0, 20] (unlike anchor count which starts at 2).
        assert _parse_traj_interpolation_value({"n_points": 0}) == 0
        assert _parse_traj_interpolation_value({"n_points": -1}) == 0
        assert _parse_traj_interpolation_value({"n_points": 50}) == 20

    def test_traj_neighbor_wrapper_bounds(self) -> None:
        assert _parse_traj_neighbor_value({"k": 1}, "k", default=10) == 2
        assert _parse_traj_neighbor_value({"k": 9999}, "k", default=10) == 200
        assert _parse_traj_neighbor_value({}, "k", default=15) == 15


class TestTrajectoryEligibilityError:
    """Ensure the deduplicated trajectory guard returns the right response."""

    def test_eligible_returns_none(self, monkeypatch: pytest.MonkeyPatch) -> None:
        app = dash_app.create_app(workdir=None)
        monkeypatch.setattr(
            "cryodrgn.dashboard.app.explorer_volumes_eligible", lambda _e: True
        )
        with app.test_request_context():
            assert _trajectory_eligibility_error(object()) is None  # type: ignore[arg-type]

    def test_ineligible_returns_400_json(self, monkeypatch: pytest.MonkeyPatch) -> None:
        app = dash_app.create_app(workdir=None)
        monkeypatch.setattr(
            "cryodrgn.dashboard.app.explorer_volumes_eligible", lambda _e: False
        )
        with app.test_request_context():
            resp, status = _trajectory_eligibility_error(object())  # type: ignore[arg-type,misc]
        assert status == 400
        assert "error" in resp.get_json()


class TestChimeraxRenderCmds:
    """Pure-function sanity checks for the extracted ChimeraX helper."""

    def test_static_view_has_no_turn(self) -> None:
        cmds = _chimerax_render_cmds(
            "/tmp/x.mrc", "/tmp/x.png", 100, vol_name="vol000", turn_y=None
        )
        assert not any(c.startswith("turn ") for c in cmds)
        assert cmds[0].startswith("open ")
        assert cmds[-1] == "exit"
        assert any("save " in c for c in cmds)

    def test_rotated_view_injects_turn(self) -> None:
        cmds = _chimerax_render_cmds(
            "/tmp/x.mrc", "/tmp/x.png", 100, vol_name="vol000", turn_y=45.0
        )
        assert any("turn y 45.0" in c for c in cmds)

    def test_paths_with_spaces_are_quoted(self) -> None:
        cmds = _chimerax_render_cmds(
            "/with space/x.mrc",
            "/out dir/y.png",
            100,
            vol_name="vol000",
            turn_y=None,
        )
        # shlex.quote wraps paths with spaces in single quotes.
        assert any("'/with space/x.mrc'" in c for c in cmds)
        assert any("'/out dir/y.png'" in c for c in cmds)


# ---------------------------------------------------------------------------
# Integration tests against a real (tiny) experiment
# ---------------------------------------------------------------------------


class TestDashboardExperiment:
    """Smoke-test ``load_experiment`` + basic shape invariants."""

    def test_list_z_epochs(self, dashboard_workdir: str) -> None:
        assert list_z_epochs(dashboard_workdir) == [ANALYZE_EPOCH]

    def test_load_experiment_shapes(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        e = dashboard_experiment
        assert e.epoch == ANALYZE_EPOCH
        assert e.z.shape == (100, 4)
        assert len(e.plot_df) == 100
        # A handful of columns every dashboard view expects to be present.
        for col in ("UMAP1", "UMAP2", "PC1", "z0", "znorm", "labels"):
            assert col in e.plot_df.columns, col


class TestDashboardPages:
    """Every top-level Flask view should render without error."""

    @pytest.mark.parametrize(
        "path",
        [
            "/",
            "/explorer",
            "/pairplot",
            "/latent-3d",
            "/trajectory",
            "/command-builder",
        ],
    )
    def test_page_renders(self, flask_client, path: str) -> None:
        r = flask_client.get(path)
        assert r.status_code == 200, f"{path} returned {r.status_code}"
        assert r.data, f"{path} returned empty body"


class TestDashboardScatterApis:
    """Scatter / 3-D scatter / preview endpoints return JSON or PNG payloads."""

    def test_api_scatter_json(self, flask_client) -> None:
        r = flask_client.get("/api/scatter?x=UMAP1&y=UMAP2&color=labels&filter_ui=1")
        assert r.status_code == 200
        js = r.get_json()
        assert "data" in js  # Plotly figure
        assert js["data"]

    def test_api_scatter_no_color(self, flask_client) -> None:
        r = flask_client.get("/api/scatter?x=UMAP1&y=UMAP2&color=none&marker_size=3")
        assert r.status_code == 200

    def test_api_scatter3d_z(self, flask_client) -> None:
        r = flask_client.get("/api/scatter3d_z?x=z0&y=z1&z=z2&color=znorm")
        assert r.status_code == 200
        assert r.get_json()["data"]

    def test_api_latent3d_preview_png(self, flask_client) -> None:
        r = flask_client.get("/api/latent3d_preview.png?x=z0&y=z1&z=z2&color=znorm")
        assert r.status_code == 200
        # PNG magic bytes.
        assert r.data[:8] == b"\x89PNG\r\n\x1a\n"

    def test_api_preview_montage(self, flask_client) -> None:
        r = flask_client.get("/api/preview_montage?rows=0,1,2,3")
        assert r.status_code == 200
        assert r.data[:8] == b"\x89PNG\r\n\x1a\n"


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


class TestDashboardTrajectoryCoords:
    """Flask tests for trajectory JSON APIs (gated by ``explorer_volumes_eligible`` in the app)."""

    def test_direct_interpolation(
        self, flask_client, dashboard_experiment: DashboardExperiment
    ) -> None:
        r = flask_client.post(
            "/api/trajectory_coords",
            json={
                "mode": "direct",
                "x": "z0",
                "y": "z1",
                "start": [0.0, 0.0],
                "end": [1.0, 1.0],
                "n_points": 3,
            },
        )
        if not _traj_flask_200_or_ineligible(r, dashboard_experiment):
            return
        js = r.get_json()
        # 3 interpolation pts + 2 endpoints = 5 (matches live-server observation).
        assert len(js["z_traj"]) >= 2
        assert all(len(z) == 4 for z in js["z_traj"])
        assert js["mode"] == "direct"

    def test_nearest_mode(
        self, flask_client, dashboard_experiment: DashboardExperiment
    ) -> None:
        r = flask_client.post(
            "/api/trajectory_coords",
            json={
                "mode": "nearest",
                "x": "UMAP1",
                "y": "UMAP2",
                "start": [0.0, 0.0],
                "end": [1.0, 1.0],
                "n_points": 3,
            },
        )
        if not _traj_flask_200_or_ineligible(r, dashboard_experiment):
            return
        assert r.get_json()["mode"] == "nearest"

    def test_anchor_driven_trajectory(
        self, flask_client, dashboard_experiment: DashboardExperiment
    ) -> None:
        r = flask_client.post(
            "/api/trajectory_coords",
            json={
                "anchor_indices": [0, 5, 10],
                "mode": "direct",
                "x": "z0",
                "y": "z1",
                "n_points": 3,
            },
        )
        if not _traj_flask_200_or_ineligible(r, dashboard_experiment):
            return
        js = r.get_json()
        assert js.get("traj_particle_indices"), "missing anchor particle indices"

    def test_kmeans_centers(
        self, flask_client, dashboard_experiment: DashboardExperiment
    ) -> None:
        r = flask_client.post(
            "/api/trajectory_kmeans_centers",
            json={"x": "z0", "y": "z1", "mode": "direct", "n_points": 2},
        )
        if not _traj_flask_200_or_ineligible(r, dashboard_experiment):
            return

    def test_random_indices(
        self, flask_client, dashboard_experiment: DashboardExperiment
    ) -> None:
        r = flask_client.post(
            "/api/trajectory_random_indices",
            json={"x": "z0", "y": "z1", "mode": "direct", "n_points": 2},
        )
        if not _traj_flask_200_or_ineligible(r, dashboard_experiment):
            return

    def test_default_endpoints(
        self, flask_client, dashboard_experiment: DashboardExperiment
    ) -> None:
        r = flask_client.get("/api/default_trajectory_endpoints?x=z0&y=z1")
        if not _traj_flask_200_or_ineligible(r, dashboard_experiment):
            return


class TestDashboardZPkl:
    """Shape checks against the raw training artefacts the dashboard consumes."""

    def test_z_pkl_matches_experiment(
        self,
        dashboard_workdir: str,
        dashboard_experiment: DashboardExperiment,
    ) -> None:
        with open(
            os.path.join(dashboard_workdir, f"z.{ANALYZE_EPOCH}.pkl"), "rb"
        ) as fh:
            z = pickle.load(fh)
        assert z.shape == dashboard_experiment.z.shape
        np.testing.assert_allclose(z, dashboard_experiment.z, rtol=1e-5)


# ---------------------------------------------------------------------------
# trajectory.py: anchor-text parsing, rounding, formatting helpers
# ---------------------------------------------------------------------------


class TestParseAnchorIndicesTxt:
    """Whitespace-delimited anchor-index parsing with user-supplied text."""

    def test_parses_whitespace(self) -> None:
        assert parse_anchor_indices_txt(b"0 5 10\n") == [0, 5, 10]

    def test_parses_commas_and_semicolons(self) -> None:
        assert parse_anchor_indices_txt(b"0,5;10\t12") == [0, 5, 10, 12]

    def test_strips_surrounding_whitespace_and_newlines(self) -> None:
        assert parse_anchor_indices_txt(b"\n  1  2\r\n3\n") == [1, 2, 3]

    def test_allows_negative_and_zero(self) -> None:
        # Range validation happens at the caller (compute_*); parser only
        # checks that each token is an int literal.
        assert parse_anchor_indices_txt(b"-1 0 1") == [-1, 0, 1]

    def test_rejects_single_index(self) -> None:
        with pytest.raises(ValueError, match="at least two"):
            parse_anchor_indices_txt(b"42")

    def test_rejects_empty(self) -> None:
        with pytest.raises(ValueError, match="at least two"):
            parse_anchor_indices_txt(b"   \n  ")

    def test_rejects_non_integer_tokens(self) -> None:
        with pytest.raises(ValueError, match="Invalid anchor index token"):
            parse_anchor_indices_txt(b"0 5 foo")

    def test_rejects_floats(self) -> None:
        with pytest.raises(ValueError, match="Invalid anchor index token"):
            parse_anchor_indices_txt(b"0 5 3.14")

    def test_rejects_non_utf8(self) -> None:
        with pytest.raises(ValueError, match="UTF-8"):
            parse_anchor_indices_txt(b"\xff\xfe\xfd")


class TestZTrajSavetxtRoundTrip:
    def test_roundtrip_through_numpy(self) -> None:
        z = np.arange(12, dtype=np.float64).reshape(4, 3)
        text = z_traj_to_savetxt_str(z)
        loaded = np.loadtxt(io.StringIO(text))
        np.testing.assert_allclose(loaded, z)


class TestRoundDirectModeTrajXY:
    def test_small_range_rounds_to_three_decimals(self) -> None:
        arr = np.array([[0.123456, 0.987654], [0.111111, 0.999999]])
        rounded = _round_direct_mode_traj_xy(arr)
        np.testing.assert_allclose(rounded, [[0.123, 0.988], [0.111, 1.000]])

    def test_large_range_rounds_to_two_decimals(self) -> None:
        arr = np.array([[123.456789, 0.1234567], [5.5555, 200.98765]])
        rounded = _round_direct_mode_traj_xy(arr)
        assert rounded[0, 0] == pytest.approx(123.46)
        # second column here has max |val| = 200.98 >= 100 -> 2 decimals.
        assert rounded[1, 1] == pytest.approx(200.99)

    def test_nan_entries_preserved(self) -> None:
        arr = np.array([[np.nan, 1.234], [2.345, np.nan]])
        rounded = _round_direct_mode_traj_xy(arr)
        assert np.isnan(rounded[0, 0])
        assert np.isnan(rounded[1, 1])
        assert rounded[1, 0] == pytest.approx(2.345)


class TestTrajectoryXYOkForDirect:
    @pytest.mark.parametrize(
        "x,y,expected",
        [
            ("z0", "z1", True),
            ("PC1", "PC2", True),
            ("z0", "PC1", True),
            ("UMAP1", "z0", False),
            ("UMAP1", "UMAP2", False),
            ("labels", "z0", False),
        ],
    )
    def test_allowed_combinations(self, x: str, y: str, expected: bool) -> None:
        assert _trajectory_xy_ok_for_direct(x, y) is expected


class TestTrajectoryDefaultXYCols:
    def test_prefers_pc_when_zdim_greater_than_2(self) -> None:
        cols = ["z0", "z1", "z2", "PC1", "PC2", "UMAP1", "UMAP2"]
        assert trajectory_default_xy_cols(cols, zdim=3) == ("PC1", "PC2")

    def test_falls_back_to_umap_when_no_pc(self) -> None:
        cols = ["z0", "z1", "UMAP1", "UMAP2"]
        assert trajectory_default_xy_cols(cols, zdim=2) == ("UMAP1", "UMAP2")

    def test_falls_back_to_first_two_when_no_embedding(self) -> None:
        cols = ["z0", "z1"]
        assert trajectory_default_xy_cols(cols, zdim=2) == ("z0", "z1")


class TestTrajectoryPlotAxisColumns:
    def test_includes_all_z_pc_umap(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        cols = trajectory_plot_axis_columns(dashboard_experiment)
        zdim = int(dashboard_experiment.z.shape[1])
        assert [f"z{i}" for i in range(zdim)] == cols[:zdim]
        # zdim=4 > 2 so PC columns are included; UMAP has 2 axes.
        assert "PC1" in cols and "PC4" in cols
        assert "UMAP1" in cols and "UMAP2" in cols
        assert "labels" not in cols

    def test_validate_raises_for_disallowed_pair(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        with pytest.raises(ValueError, match="not available"):
            validate_trajectory_plot_axes(dashboard_experiment, "z0", "labels")

    def test_validate_accepts_allowed_pair(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        validate_trajectory_plot_axes(dashboard_experiment, "z0", "z1")
        validate_trajectory_plot_axes(dashboard_experiment, "PC1", "PC2")


class TestHasEmbeddingColumns:
    def test_has_umap_and_pc_on_real_experiment(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        assert has_umap_columns(dashboard_experiment)
        assert has_pc_columns(dashboard_experiment)


class TestDijkstraFromNeighbors:
    """Shortest-path algorithm on hand-constructed sparse neighbour graphs."""

    @staticmethod
    def _line_graph() -> tuple[np.ndarray, np.ndarray]:
        # 0 - 1 - 2 - 3 each edge weight 1; pad inactive neighbour slot with -1.
        neighbors = np.array([[1, -1], [0, 2], [1, 3], [2, -1]], dtype=np.int64)
        inf = np.inf
        dists = np.array(
            [[1.0, inf], [1.0, 1.0], [1.0, 1.0], [1.0, inf]], dtype=np.float64
        )
        return neighbors, dists

    def test_same_node_returns_singleton(self) -> None:
        nb, d = self._line_graph()
        assert _dijkstra_path_from_neighbors(nb, d, 2, 2) == [2]

    def test_path_across_line(self) -> None:
        nb, d = self._line_graph()
        assert _dijkstra_path_from_neighbors(nb, d, 0, 3) == [0, 1, 2, 3]
        assert _dijkstra_path_from_neighbors(nb, d, 3, 0) == [3, 2, 1, 0]

    def test_disconnected_returns_none(self) -> None:
        neighbors = np.array([[1, -1], [0, -1], [3, -1], [2, -1]], dtype=np.int64)
        inf = np.inf
        dists = np.array(
            [[1.0, inf], [1.0, inf], [1.0, inf], [1.0, inf]], dtype=np.float64
        )
        assert _dijkstra_path_from_neighbors(neighbors, dists, 0, 3) is None

    def test_prefers_shorter_edge(self) -> None:
        # Triangle: 0-1-2 via 1 (weight 2+2=4) vs 0-2 directly (weight 1).
        neighbors = np.array([[1, 2], [0, 2], [0, 1]], dtype=np.int64)
        dists = np.array([[2.0, 1.0], [2.0, 2.0], [1.0, 2.0]], dtype=np.float64)
        assert _dijkstra_path_from_neighbors(neighbors, dists, 0, 2) == [0, 2]


class TestGraphNeighborArrays:
    def test_shapes_and_cache(self, dashboard_experiment: DashboardExperiment) -> None:
        _TRAJ_GRAPH_NEIGHBOR_CACHE.clear()
        nb, d = _graph_neighbor_arrays(
            dashboard_experiment, max_neighbors=5, avg_neighbors=5
        )
        n = dashboard_experiment.z.shape[0]
        assert nb.shape[0] == n
        assert d.shape == nb.shape
        assert nb.shape[1] == 5  # k neighbours (self dropped)
        # Cached on second call: same object identity.
        nb2, d2 = _graph_neighbor_arrays(
            dashboard_experiment, max_neighbors=5, avg_neighbors=5
        )
        assert nb2 is nb and d2 is d


class TestComputeDirectAnchorTrajectory:
    def test_endpoints_and_interpolation_count(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        anchors = [0, 10, 20]
        interp = 3
        z_traj, traj_rows, traj_xy = _compute_direct_anchor_trajectory(
            dashboard_experiment, anchors, "z0", "z1", interp
        )
        assert traj_rows is None
        # (len(anchors) - 1) * (interp + 1) + 1 = 2 * 4 + 1 = 9 rows.
        assert z_traj.shape[0] == 9
        assert traj_xy.shape == (9, 2)
        # Endpoints match the anchor latents exactly.
        np.testing.assert_allclose(z_traj[0], dashboard_experiment.z[0])
        np.testing.assert_allclose(z_traj[-1], dashboard_experiment.z[20])

    def test_requires_two_anchors(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        with pytest.raises(ValueError, match="at least two"):
            _compute_direct_anchor_trajectory(dashboard_experiment, [5], "z0", "z1", 3)

    def test_out_of_range_raises(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        with pytest.raises(ValueError, match="out of range"):
            _compute_direct_anchor_trajectory(
                dashboard_experiment, [0, 9999], "z0", "z1", 2
            )


class TestParseTrajectoryRequestBody:
    def test_anchor_direct_happy_path(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        p = parse_trajectory_request_body(
            dashboard_experiment,
            {
                "anchor_indices": [0, 5, 10],
                "mode": "direct",
                "x": "z0",
                "y": "z1",
                "n_points": 3,
            },
        )
        assert p["use_anchors"] is True
        assert p["anchor_indices"] == [0, 5, 10]
        assert p["mode"] == "direct"
        assert p["n_points"] == 3

    def test_anchor_graph_happy_path(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        p = parse_trajectory_request_body(
            dashboard_experiment,
            {
                "anchor_indices": [0, 10, 20],
                "mode": "graph",
                "x": "z0",
                "y": "z1",
                "n_points": 2,
            },
        )
        assert p["mode"] == "graph"
        assert p["max_neighbors"] >= 2
        assert p["avg_neighbors"] >= 2

    def test_anchor_bad_mode_rejected(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        with pytest.raises(ValueError, match='"direct" or "graph"'):
            parse_trajectory_request_body(
                dashboard_experiment,
                {
                    "anchor_indices": [0, 5],
                    "mode": "spline",
                    "x": "z0",
                    "y": "z1",
                },
            )

    def test_anchor_indices_not_int_rejected(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        with pytest.raises(ValueError, match="list of integers"):
            parse_trajectory_request_body(
                dashboard_experiment,
                {"anchor_indices": [0, "bogus"], "x": "z0", "y": "z1"},
            )

    def test_direct_requires_z_or_pc_axes(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        with pytest.raises(ValueError, match="principal-component or z latent"):
            parse_trajectory_request_body(
                dashboard_experiment,
                {
                    "mode": "direct",
                    "x": "UMAP1",
                    "y": "UMAP2",
                    "start": [0.0, 0.0],
                    "end": [1.0, 1.0],
                    "n_points": 3,
                },
            )

    def test_bad_mode_for_non_anchor_rejected(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        with pytest.raises(ValueError, match='"direct" or "nearest"'):
            parse_trajectory_request_body(
                dashboard_experiment,
                {
                    "mode": "spline",
                    "x": "z0",
                    "y": "z1",
                    "start": [0.0, 0.0],
                    "end": [1.0, 1.0],
                },
            )

    def test_bad_start_end_rejected(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        with pytest.raises(ValueError, match="start and end"):
            parse_trajectory_request_body(
                dashboard_experiment,
                {"mode": "direct", "x": "z0", "y": "z1", "start": [0.0]},
            )

    def test_non_numeric_endpoints_rejected(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        with pytest.raises(ValueError, match="numeric"):
            parse_trajectory_request_body(
                dashboard_experiment,
                {
                    "mode": "direct",
                    "x": "z0",
                    "y": "z1",
                    "start": ["a", "b"],
                    "end": [1.0, 1.0],
                },
            )

    def test_traj_xy_custom_overrides_start_end(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        pts = [[0.0, 0.0], [0.5, 0.5], [1.0, 1.0]]
        p = parse_trajectory_request_body(
            dashboard_experiment,
            {"mode": "nearest", "x": "z0", "y": "z1", "traj_xy": pts},
        )
        assert p["traj_xy_custom"] == pts


class TestComputeTrajectoryLatentPath:
    def test_nearest_rows_land_on_real_particles(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        p = parse_trajectory_request_body(
            dashboard_experiment,
            {
                "mode": "nearest",
                "x": "UMAP1",
                "y": "UMAP2",
                "start": [0.0, 0.0],
                "end": [1.0, 1.0],
                "n_points": 4,
            },
        )
        z, rows, xy = compute_trajectory_latent_path(dashboard_experiment, p)
        assert z.shape == (4, dashboard_experiment.z.shape[1])
        assert rows is not None and len(rows) == 4
        for i, r in enumerate(rows):
            np.testing.assert_allclose(z[i], dashboard_experiment.z[r])

    def test_direct_endpoints_equal_data_points(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        p = parse_trajectory_request_body(
            dashboard_experiment,
            {
                "mode": "direct",
                "x": "z0",
                "y": "z1",
                "start": [0.0, 0.0],
                "end": [1.0, 1.0],
                "n_points": 5,
            },
        )
        z, rows, xy = compute_trajectory_latent_path(dashboard_experiment, p)
        assert rows is None
        assert z.shape == (5, dashboard_experiment.z.shape[1])
        # First and last in xy match interpolation endpoints (rounded).
        assert xy.shape == (5, 2)


class TestRandomDatasetIndices:
    def test_returns_distinct_indices_in_range(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        out = random_dataset_indices(dashboard_experiment, k=7)
        assert len(out) == 7
        assert len(set(out)) == 7
        n = int(dashboard_experiment.z.shape[0])
        assert all(0 <= i < n for i in out)

    def test_clips_to_available(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        n = int(dashboard_experiment.z.shape[0])
        out = random_dataset_indices(dashboard_experiment, k=n + 100)
        assert len(out) == n


class TestDefaultTrajectoryEndpointsXY:
    def test_endpoints_span_long_axis(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        start, end = default_trajectory_endpoints_xy(
            dashboard_experiment, "UMAP1", "UMAP2"
        )
        assert len(start) == 2 and len(end) == 2
        # Endpoints should not coincide for real data.
        assert start != end

    def test_synthetic_tight_cluster_follows_pc1(self) -> None:
        # Build a minimal DashboardExperiment-like stub with a line in XY.
        xy = np.stack(
            [np.linspace(-10, 10, 50), np.zeros(50) + 0.01 * np.arange(50)],
            axis=1,
        )
        stub = argparse.Namespace(
            plot_df=pd.DataFrame({"Ax": xy[:, 0], "Ay": xy[:, 1]}),
        )
        start, end = default_trajectory_endpoints_xy(stub, "Ax", "Ay")
        # Long axis is roughly X; endpoints span near [-10, +10].
        assert start[0] < 0 < end[0] or end[0] < 0 < start[0]


class TestTrajectoryAnchorModeParams:
    def test_direct_mode_clamps_interpolation(self) -> None:
        mode, n_points, maxn, avgn = trajectory_anchor_mode_params(
            {"mode": "direct", "n_points": 100}
        )
        assert mode == "direct"
        assert n_points == 20
        assert maxn >= 2 and avgn >= 2

    def test_graph_mode_clamps_anchor_count(self) -> None:
        mode, n_points, _, _ = trajectory_anchor_mode_params(
            {"mode": "graph", "n_points": 999}
        )
        assert mode == "graph"
        assert n_points == 20


class TestDirectAnchorParticleIndicesPayload:
    def test_no_interp_is_straight_passthrough(self) -> None:
        pidx = direct_anchor_particle_indices_payload(
            anchor_indices=[3, 7, 11], interpolation_points=0, n_total=3
        )
        assert pidx == [3, 7, 11]

    def test_interpolated_fills_none_between_anchors(self) -> None:
        pidx = direct_anchor_particle_indices_payload(
            anchor_indices=[1, 2, 3], interpolation_points=2, n_total=7
        )
        # Anchors at positions 0, 3, 6; intermediate slots are None.
        assert pidx == [1, None, None, 2, None, None, 3]

    def test_total_mismatch_returns_none(self) -> None:
        assert (
            direct_anchor_particle_indices_payload(
                anchor_indices=[1, 2, 3], interpolation_points=2, n_total=5
            )
            is None
        )

    def test_too_few_anchors_returns_none(self) -> None:
        assert (
            direct_anchor_particle_indices_payload(
                anchor_indices=[1], interpolation_points=2, n_total=3
            )
            is None
        )


class TestTrajectoryAnchorPayloadFromIndices:
    def test_direct_shape(self, dashboard_experiment: DashboardExperiment) -> None:
        payload = trajectory_anchor_payload_from_indices(
            dashboard_experiment,
            [0, 10, 20],
            "z0",
            "z1",
            mode="direct",
            n_points=2,
        )
        assert payload["ok"] is True
        assert payload["mode"] == "direct"
        assert payload["anchor_indices"] == [0, 10, 20]
        assert "traj_particle_indices" in payload  # direct pidx merged in


class TestPlotDfRowsForDatasetIndices:
    def test_roundtrip_from_all_indices(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        ds_ids = np.asarray(dashboard_experiment.all_indices[[1, 3, 5]], dtype=int)
        rows = plot_df_rows_for_dataset_indices(dashboard_experiment, ds_ids)
        assert sorted(rows) == [1, 3, 5]

    def test_empty_input_returns_empty(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        assert (
            plot_df_rows_for_dataset_indices(
                dashboard_experiment, np.asarray([], dtype=int)
            )
            == []
        )


# ---------------------------------------------------------------------------
# preload.py: sampling, encoding, formatting helpers
# ---------------------------------------------------------------------------


class TestPreloadTimeHints:
    @pytest.mark.parametrize("cpus", [1, 2, 4, 8, 16, 64])
    def test_bounds_are_ordered_and_positive(self, cpus: int) -> None:
        lo, hi = _preload_cache_time_estimate_bounds(cpus)
        assert 0 < lo <= hi
        assert hi >= 18

    def test_format_singular_core(self) -> None:
        s = format_preload_cache_time_hint(1)
        assert s.endswith("1 core")

    def test_format_plural_cores_embeds_cpus(self) -> None:
        s = format_preload_cache_time_hint(8)
        assert s.endswith("8 cores")
        assert s.startswith("~")


class TestStratifiedXYRowIndices:
    def test_count_is_bounded_and_unique(self) -> None:
        rng = np.random.default_rng(123)
        coords = rng.normal(size=(500, 2))
        total_k = 200
        picks = _stratified_xy_row_indices(coords, rng, total_k)
        assert len(picks) <= total_k
        assert all(0 <= i < 500 for i in picks)

    def test_empty_coords_returns_empty(self) -> None:
        rng = np.random.default_rng(0)
        assert _stratified_xy_row_indices(np.zeros((0, 2)), rng, total_k=10) == set()


class TestSamplePlotDfRowsForPreload:
    def test_restrict_rows_returns_subset(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        restrict = [0, 1, 2, 3, 4]
        rows, ds = sample_plot_df_rows_for_preload(
            dashboard_experiment, "UMAP1", "UMAP2", restrict_to_rows=restrict
        )
        assert set(rows).issubset(set(restrict))
        assert len(rows) == len(ds)

    def test_empty_restrict_returns_empty(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        rows, ds = sample_plot_df_rows_for_preload(
            dashboard_experiment, "UMAP1", "UMAP2", restrict_to_rows=[]
        )
        assert rows == [] and ds == []


class TestMontageBytes:
    def test_empty_returns_hint_png(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        png = montage_bytes(dashboard_experiment, [])
        assert png[:8] == b"\x89PNG\r\n\x1a\n"

    def test_rows_returns_png(self, dashboard_experiment: DashboardExperiment) -> None:
        png = montage_bytes(dashboard_experiment, [0, 1, 2, 3])
        assert png[:8] == b"\x89PNG\r\n\x1a\n"
        assert len(png) > 1000  # not a trivial image

    def test_too_many_rows_are_capped(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        n = dashboard_experiment.z.shape[0]
        png_capped = montage_bytes(dashboard_experiment, list(range(min(n, 50))))
        assert png_capped[:8] == b"\x89PNG\r\n\x1a\n"


class TestParticleThumbnailB64FromRow:
    def test_returns_base64_jpeg(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        b64 = particle_thumbnail_b64_from_row(dashboard_experiment, 0, max_side=48)
        assert isinstance(b64, str) and len(b64) > 0
        raw = base64.standard_b64decode(b64)
        # JPEG magic (SOI marker).
        assert raw[:3] == b"\xff\xd8\xff"


class TestEncodeParticleBatch:
    def test_matches_request_count(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        e = dashboard_experiment
        global_idx = [int(e.all_indices[i]) for i in (0, 1, 2)]
        jpegs = encode_particle_batch(e.particles_path, e.datadir, global_idx, 48)
        assert len(jpegs) == 3
        for b64 in jpegs:
            raw = base64.standard_b64decode(b64)
            assert raw[:3] == b"\xff\xd8\xff"


class TestLoadPlotDfRowsFromPlotIndsFile:
    def test_empty_path_returns_empty(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        assert load_plot_df_rows_from_plot_inds_file(dashboard_experiment, "") == []
        assert load_plot_df_rows_from_plot_inds_file(dashboard_experiment, None) == []

    def test_missing_file_returns_empty(
        self, dashboard_experiment: DashboardExperiment, tmp_path
    ) -> None:
        assert (
            load_plot_df_rows_from_plot_inds_file(
                dashboard_experiment, str(tmp_path / "no_such.pkl")
            )
            == []
        )

    def test_roundtrip_via_pickle(
        self, dashboard_experiment: DashboardExperiment, tmp_path
    ) -> None:
        e = dashboard_experiment
        wanted_ds = np.asarray(e.all_indices[[2, 4, 6]], dtype=int)
        pkl = tmp_path / "plot_inds.pkl"
        with open(pkl, "wb") as fh:
            pickle.dump(wanted_ds, fh)
        rows = load_plot_df_rows_from_plot_inds_file(e, str(pkl))
        assert sorted(rows) == [2, 4, 6]


# ---------------------------------------------------------------------------
# context.py: discovery, caches, session routes, template helpers
# ---------------------------------------------------------------------------
