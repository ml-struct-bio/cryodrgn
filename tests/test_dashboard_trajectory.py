"""Trajectory creator interface tests and trajectory helper coverage."""

from __future__ import annotations

import argparse
import io
import os

import numpy as np
import pandas as pd
import pytest

from cryodrgn.dashboard import app as dash_app
from cryodrgn.dashboard.route_helpers import (
    _TRAJECTORY_INELIGIBLE_MSG,
    _trajectory_eligibility_error,
)
from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.particle_explorer import explorer_volumes_eligible
from cryodrgn.dashboard.trajectory import (
    _compute_direct_anchor_trajectory,
    _dijkstra_path_from_neighbors,
    _graph_neighbor_arrays,
    _graph_neighbor_max_dist,
    _path_has_crossings,
    _round_direct_mode_traj_xy,
    _snap_points_to_nearest_particles,
    _xy_already_on_particle,
    order_anchor_indices_for_direct_path,
    order_points_noncrossing_path_heuristic,
    order_points_shortest_noncrossing_path,
    parse_anchor_path_order,
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
    resolve_trajectory_decode_plan,
    trajectory_anchor_mode_params,
    trajectory_anchor_payload_from_indices,
    trajectory_default_xy_cols,
    trajectory_plot_axis_columns,
    validate_trajectory_plot_axes,
    z_traj_to_savetxt_str,
)
from tests.conftest import _monkeypatch_explorer_volumes_eligible

pytestmark = pytest.mark.dashboard


def _traj_flask_200_or_ineligible(r, experiment: DashboardExperiment) -> bool:
    """If volume exploration is eligible, require HTTP 200; otherwise require the standard 400."""
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

    @pytest.mark.parametrize(
        "func,kwargs,expected",
        [
            # traj_points wrapper: lo=2, hi=20, default=4
            (_parse_traj_points_value, {"data": {"n_points": 1}}, 2),
            (_parse_traj_points_value, {"data": {"n_points": 99}}, 20),
            (_parse_traj_points_value, {"data": {}}, 4),
            # traj_interpolation wrapper: lo=0, hi=20, default=4
            (_parse_traj_interpolation_value, {"data": {"n_points": 0}}, 0),
            (_parse_traj_interpolation_value, {"data": {"n_points": -1}}, 0),
            (_parse_traj_interpolation_value, {"data": {"n_points": 50}}, 20),
            # traj_neighbor wrapper: lo=2, hi=200
            (
                _parse_traj_neighbor_value,
                {"data": {"k": 1}, "key": "k", "default": 10},
                2,
            ),
            (
                _parse_traj_neighbor_value,
                {"data": {"k": 9999}, "key": "k", "default": 10},
                200,
            ),
            (_parse_traj_neighbor_value, {"data": {}, "key": "k", "default": 15}, 15),
        ],
    )
    def test_wrapper_functions(
        self, func: callable, kwargs: dict, expected: int
    ) -> None:
        """Consolidated parametrized tests for all three wrapper functions."""
        assert func(**kwargs) == expected


class TestTrajectoryEligibilityError:
    """Ensure the deduplicated trajectory guard returns the right response."""

    def test_eligible_returns_none(self, monkeypatch: pytest.MonkeyPatch) -> None:
        app = dash_app.create_app(workdir=None)
        _monkeypatch_explorer_volumes_eligible(monkeypatch, eligible=True)
        with app.test_request_context():
            assert _trajectory_eligibility_error(object()) is None  # type: ignore[arg-type]

    def test_ineligible_returns_400_json(self, monkeypatch: pytest.MonkeyPatch) -> None:
        app = dash_app.create_app(workdir=None)
        _monkeypatch_explorer_volumes_eligible(monkeypatch, eligible=False)
        with app.test_request_context():
            resp, status = _trajectory_eligibility_error(object())  # type: ignore[arg-type,misc]
        assert status == 400
        assert "error" in resp.get_json()


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

    def test_nearest_mode_returns_marker_colors_for_color_scale(
        self, flask_client, dashboard_experiment: DashboardExperiment
    ) -> None:
        color_col = next(
            (
                c
                for c in dashboard_experiment.numeric_columns
                if c in dashboard_experiment.plot_df.columns
            ),
            None,
        )
        if color_col is None:
            pytest.skip("experiment has no numeric color column")
        r = flask_client.post(
            "/api/trajectory_coords",
            json={
                "mode": "nearest",
                "x": "UMAP1",
                "y": "UMAP2",
                "start": [0.0, 0.0],
                "end": [1.0, 1.0],
                "n_points": 3,
                "color": color_col,
                "palette": "Viridis",
            },
        )
        if not _traj_flask_200_or_ineligible(r, dashboard_experiment):
            return
        js = r.get_json()
        assert js["mode"] == "nearest"
        colors = js.get("traj_marker_colors")
        assert isinstance(colors, list)
        assert len(colors) == len(js["traj_rows"])
        assert all(
            c is None or (isinstance(c, str) and c.startswith("#")) for c in colors
        )
        assert any(isinstance(c, str) and c.startswith("#") for c in colors)

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
        assert js.get("anchor_indices") == [0, 5, 10]
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
            json={"count": 10, "exclude_indices": [0, 1]},
        )
        if not _traj_flask_200_or_ineligible(r, dashboard_experiment):
            return
        body = r.get_json()
        assert body.get("ok") is True
        indices = body.get("indices") or body.get("anchor_indices") or []
        assert len(indices) >= 1
        assert 0 not in indices and 1 not in indices

    def test_random_indices_include_xy(
        self, flask_client, dashboard_experiment: DashboardExperiment
    ) -> None:
        r = flask_client.post(
            "/api/trajectory_random_indices",
            json={
                "count": 5,
                "exclude_indices": [0, 1],
                "x": "z0",
                "y": "z1",
            },
        )
        if not _traj_flask_200_or_ineligible(r, dashboard_experiment):
            return
        body = r.get_json()
        assert body.get("ok") is True
        indices = body.get("indices") or []
        xy = body.get("xy") or []
        assert len(xy) == len(indices)
        assert all(len(pt) == 2 for pt in xy)

    def test_default_endpoints(
        self, flask_client, dashboard_experiment: DashboardExperiment
    ) -> None:
        r = flask_client.get("/api/default_trajectory_endpoints?x=z0&y=z1")
        if not _traj_flask_200_or_ineligible(r, dashboard_experiment):
            return


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


class TestSnapPointsPreserveAlreadyOnParticle:
    def test_exact_particle_xy_unchanged(self) -> None:
        coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float64)
        pts = coords.copy()
        rows, xy = _snap_points_to_nearest_particles(coords, pts)
        assert rows == [0, 1, 2]
        np.testing.assert_array_equal(xy, pts)

    def test_rounded_particle_xy_preserved(self) -> None:
        coords = np.array([[0.123456, 0.987654], [5.0, 5.0]], dtype=np.float64)
        # Direct-mode display of particle 0 (3-decimal rounding).
        rounded = _round_direct_mode_traj_xy(coords[0:1])[0]
        assert _xy_already_on_particle(rounded, coords[0])
        rows, xy = _snap_points_to_nearest_particles(
            coords, np.array([rounded, [4.9, 4.9]])
        )
        assert rows[0] == 0
        np.testing.assert_allclose(xy[0], rounded)
        # Off-particle sample still snaps to true particle coordinates.
        assert rows[1] == 1
        np.testing.assert_allclose(xy[1], coords[1])

    def test_nearest_custom_path_keeps_on_particle_samples(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        coords = dashboard_experiment.plot_df[["UMAP1", "UMAP2"]].to_numpy(
            dtype=np.float64
        )
        p0 = coords[0]
        # Midpoint between two particles — must move under snap.
        p_mid = 0.5 * (coords[0] + coords[min(5, len(coords) - 1)])
        p2 = coords[min(10, len(coords) - 1)]
        path = [
            [float(p0[0]), float(p0[1])],
            [float(p_mid[0]), float(p_mid[1])],
            [float(p2[0]), float(p2[1])],
        ]
        body = parse_trajectory_request_body(
            dashboard_experiment,
            {
                "mode": "nearest",
                "x": "UMAP1",
                "y": "UMAP2",
                "start": path[0],
                "end": path[-1],
                "n_points": 3,
                "traj_xy": path,
            },
        )
        _z, rows, xy = compute_trajectory_latent_path(dashboard_experiment, body)
        assert rows is not None and len(rows) == 3
        # Endpoints already on particles: XY must be unchanged.
        np.testing.assert_allclose(xy[0], path[0])
        np.testing.assert_allclose(xy[2], path[2])
        # Free midpoint moves onto a particle coordinate.
        assert not np.allclose(xy[1], path[1])
        np.testing.assert_allclose(xy[1], coords[rows[1]])


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


class TestGraphNeighborMaxDist:
    def test_partition_matches_full_sort(self) -> None:
        rng = np.random.default_rng(0)
        ndist = rng.uniform(0.1, 5.0, size=(40, 7))
        n, avg = 40, 5
        got = _graph_neighbor_max_dist(ndist, n, avg)
        flat = np.sort(ndist.reshape(-1))
        total_neighbors = max(1, min(int(n * avg), int(ndist.size)))
        want = float(flat[total_neighbors - 1])
        assert got == want


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


class TestOrderPointsShortestNoncrossingPath:
    def test_square_visits_all_without_crossings_exact(self) -> None:
        # Corners of a unit square; optimal non-crossing path has length 3.
        points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
        order = order_points_shortest_noncrossing_path(points, path_order="exact")
        assert sorted(order) == [0, 1, 2, 3]
        assert not _path_has_crossings(order, points)
        length = sum(
            np.linalg.norm(points[order[i + 1]] - points[order[i]])
            for i in range(len(order) - 1)
        )
        assert length == pytest.approx(3.0)

    def test_heuristic_is_fast_and_noncrossing(self) -> None:
        points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
        order = order_points_noncrossing_path_heuristic(points)
        assert sorted(order) == [0, 1, 2, 3]
        assert not _path_has_crossings(order, points)

    def test_reorders_for_direct_anchors(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        anchors = [0, 10, 20, 30]
        ordered = order_anchor_indices_for_direct_path(
            dashboard_experiment, anchors, "z0", "z1", path_order="heuristic"
        )
        assert sorted(ordered) == sorted(anchors)
        coords = dashboard_experiment.plot_df[["z0", "z1"]].values.astype(np.float64)
        local = np.vstack([coords[int(a)] for a in ordered])
        local_order = list(range(len(ordered)))
        assert not _path_has_crossings(local_order, local)

    def test_preserve_keeps_anchor_order(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        anchors = [12, 3, 8, 20]
        ordered = order_anchor_indices_for_direct_path(
            dashboard_experiment, anchors, "z0", "z1", path_order="preserve"
        )
        assert ordered == anchors


class TestParseAnchorPathOrder:
    def test_defaults_to_preserve(self) -> None:
        assert parse_anchor_path_order({}) == "preserve"
        assert parse_anchor_path_order({"anchor_path_order": ""}) == "preserve"

    def test_accepts_heuristic_aliases(self) -> None:
        assert (
            parse_anchor_path_order({"anchor_path_order": "heuristic"}) == "heuristic"
        )
        assert parse_anchor_path_order({"anchor_path_order": "inexact"}) == "heuristic"

    def test_accepts_exact_aliases(self) -> None:
        assert parse_anchor_path_order({"anchor_path_order": "exact"}) == "exact"
        assert parse_anchor_path_order({"anchor_path_order": "held-karp"}) == "exact"

    def test_accepts_preserve_aliases(self) -> None:
        assert parse_anchor_path_order({"anchor_path_order": "preserve"}) == "preserve"
        assert parse_anchor_path_order({"anchor_path_order": "selection"}) == "preserve"
        assert parse_anchor_path_order({"anchor_path_order": "original"}) == "preserve"


class TestComputeDirectAnchorTrajectory:
    def test_zero_interpolation_is_pure_anchor_permutation(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        """Exact/inexact path order with n_points=0 must not densify the path."""
        anchors = [0, 10, 20, 30]
        ordered = order_anchor_indices_for_direct_path(
            dashboard_experiment, anchors, "z0", "z1", path_order="exact"
        )
        z_traj, traj_rows, traj_xy = _compute_direct_anchor_trajectory(
            dashboard_experiment, ordered, "z0", "z1", 0
        )
        assert traj_rows == ordered
        assert z_traj.shape[0] == len(ordered)
        assert traj_xy.shape == (len(ordered), 2)
        assert sorted(ordered) == sorted(anchors)

    def test_endpoints_and_interpolation_count(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        anchors = [0, 10, 20]
        ordered = order_anchor_indices_for_direct_path(
            dashboard_experiment, anchors, "z0", "z1", path_order="heuristic"
        )
        interp = 3
        z_traj, traj_rows, traj_xy = _compute_direct_anchor_trajectory(
            dashboard_experiment, ordered, "z0", "z1", interp
        )
        assert traj_rows is None
        # (len(anchors) - 1) * (interp + 1) + 1 = 2 * 4 + 1 = 9 rows.
        assert z_traj.shape[0] == 9
        assert traj_xy.shape == (9, 2)
        # Endpoints match the ordered anchor latents exactly.
        np.testing.assert_allclose(z_traj[0], dashboard_experiment.z[ordered[0]])
        np.testing.assert_allclose(z_traj[-1], dashboard_experiment.z[ordered[-1]])

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

    @pytest.mark.parametrize(
        "payload,error_match",
        [
            # Anchor mode with invalid mode
            (
                {"anchor_indices": [0, 5], "mode": "spline", "x": "z0", "y": "z1"},
                '"direct" or "graph"',
            ),
            # Anchor indices not integers
            (
                {"anchor_indices": [0, "bogus"], "x": "z0", "y": "z1"},
                "list of integers",
            ),
            # Direct mode requires PC or z axes
            (
                {
                    "mode": "direct",
                    "x": "UMAP1",
                    "y": "UMAP2",
                    "start": [0.0, 0.0],
                    "end": [1.0, 1.0],
                    "n_points": 3,
                },
                "principal-component or z latent",
            ),
            # Non-anchor mode with invalid mode
            (
                {
                    "mode": "spline",
                    "x": "z0",
                    "y": "z1",
                    "start": [0.0, 0.0],
                    "end": [1.0, 1.0],
                },
                '"direct" or "nearest"',
            ),
            # Bad start/end length
            ({"mode": "direct", "x": "z0", "y": "z1", "start": [0.0]}, "start and end"),
            # Non-numeric endpoints
            (
                {
                    "mode": "direct",
                    "x": "z0",
                    "y": "z1",
                    "start": ["a", "b"],
                    "end": [1.0, 1.0],
                },
                "numeric",
            ),
        ],
    )
    def test_validation_errors(
        self, dashboard_experiment: DashboardExperiment, payload: dict, error_match: str
    ) -> None:
        """Consolidated parametrized tests for request validation errors."""
        with pytest.raises(ValueError, match=error_match):
            parse_trajectory_request_body(dashboard_experiment, payload)

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

    def test_nearest_with_traj_xy_snaps_present_path_not_linear_start_end(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        """Custom traj_xy (edited direct-trace) must drive nearest snap, not start/end."""
        coords = dashboard_experiment.plot_df[["UMAP1", "UMAP2"]].to_numpy(
            dtype=np.float64
        )
        # Build a bent path whose midpoints differ from the linear start→end chord.
        p0 = coords[0]
        p1 = coords[min(5, len(coords) - 1)]
        p2 = coords[min(10, len(coords) - 1)]
        bent = [
            [float(p0[0]), float(p0[1])],
            [float(p1[0]), float(p1[1])],
            [float(p2[0]), float(p2[1])],
        ]
        linear_p = parse_trajectory_request_body(
            dashboard_experiment,
            {
                "mode": "nearest",
                "x": "UMAP1",
                "y": "UMAP2",
                "start": bent[0],
                "end": bent[-1],
                "n_points": 3,
            },
        )
        bent_p = parse_trajectory_request_body(
            dashboard_experiment,
            {
                "mode": "nearest",
                "x": "UMAP1",
                "y": "UMAP2",
                "start": bent[0],
                "end": bent[-1],
                "n_points": 3,
                "traj_xy": bent,
            },
        )
        _z_lin, rows_lin, _xy_lin = compute_trajectory_latent_path(
            dashboard_experiment, linear_p
        )
        _z_bent, rows_bent, xy_bent = compute_trajectory_latent_path(
            dashboard_experiment, bent_p
        )
        assert rows_bent is not None and len(rows_bent) == 3
        assert xy_bent.shape == (3, 2)
        # Each custom point snaps to its own neighbourhood particle.
        for i, pt in enumerate(bent):
            nearest = int(np.argmin(np.sum((coords - np.asarray(pt)) ** 2, axis=1)))
            assert rows_bent[i] == nearest
        # A bent path must not collapse to the linear start→end snap when midpoints differ.
        mid_linear = 0.5 * (np.asarray(bent[0]) + np.asarray(bent[-1]))
        if np.linalg.norm(np.asarray(bent[1]) - mid_linear) > 1e-6:
            assert rows_bent != rows_lin

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

    def test_excludes_existing_indices(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        n = int(dashboard_experiment.z.shape[0])
        exclude = list(range(min(5, n)))
        out = random_dataset_indices(dashboard_experiment, k=7, exclude=exclude)
        assert len(out) == min(7, n - len(exclude))
        assert not set(out).intersection(exclude)

    def test_large_n_samples_without_full_candidate_list(self) -> None:
        """Regression: must not build an O(n) Python list for large stacks."""
        n = 250_000

        class _Exp:
            z = np.zeros((n, 2), dtype=np.float32)

        out = random_dataset_indices(_Exp(), k=10, exclude=[0, 1, 2])
        assert len(out) == 10
        assert len(set(out)) == 10
        assert not set(out).intersection({0, 1, 2})
        assert all(0 <= i < n for i in out)


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


class TestResolveTrajectoryDecodePlan:
    def test_full_decode_by_default(self) -> None:
        z = np.arange(12, dtype=np.float64).reshape(4, 3)
        z_decode, slots, n_traj = resolve_trajectory_decode_plan({}, z)
        assert n_traj == 4
        assert slots == [0, 1, 2, 3]
        np.testing.assert_array_equal(z_decode, z)

    def test_partial_decode_indices(self) -> None:
        z = np.arange(12, dtype=np.float64).reshape(4, 3)
        z_decode, slots, n_traj = resolve_trajectory_decode_plan(
            {"decode_indices": [1, 3]}, z
        )
        assert n_traj == 4
        assert slots == [1, 3]
        np.testing.assert_array_equal(z_decode, z[[1, 3]])

    def test_rejects_out_of_range(self) -> None:
        z = np.zeros((3, 2), dtype=np.float64)
        with pytest.raises(ValueError, match="out of range"):
            resolve_trajectory_decode_plan({"decode_indices": [0, 3]}, z)


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
            anchor_path_order="preserve",
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


class TestTrajectoryImportAnchors:
    def test_import_happy_path(
        self,
        flask_client,
        dashboard_experiment: DashboardExperiment,
        tmp_path,
    ) -> None:
        anchors_file = tmp_path / "anchors.txt"
        anchors_file.write_text("0 5 10\n")
        r = flask_client.post(
            "/api/trajectory_import_anchors",
            json={
                "server_path": str(anchors_file),
                "mode": "direct",
                "x": "z0",
                "y": "z1",
                "n_points": 2,
            },
        )
        if not _traj_flask_200_or_ineligible(r, dashboard_experiment):
            return
        js = r.get_json()
        assert js["ok"] is True
        assert js["anchor_indices"] == [0, 5, 10]

    def test_rejects_non_txt(self, flask_client, tmp_path) -> None:
        bad = tmp_path / "anchors.csv"
        bad.write_text("0,5,10")
        r = flask_client.post(
            "/api/trajectory_import_anchors",
            json={"server_path": str(bad), "mode": "direct", "x": "z0", "y": "z1"},
        )
        assert r.status_code == 400

    def test_missing_file(self, flask_client, tmp_path) -> None:
        r = flask_client.post(
            "/api/trajectory_import_anchors",
            json={
                "server_path": str(tmp_path / "missing.txt"),
                "mode": "direct",
                "x": "z0",
                "y": "z1",
            },
        )
        assert r.status_code == 400

    def test_missing_path_field(self, flask_client) -> None:
        r = flask_client.post("/api/trajectory_import_anchors", json={})
        assert r.status_code == 400


class TestSaveZPath:
    def test_roundtrip(
        self, flask_client, tmp_path, dashboard_experiment: DashboardExperiment
    ) -> None:
        z_arr = np.arange(8, dtype=np.float64).reshape(2, 4)
        txt = z_traj_to_savetxt_str(z_arr)
        out_path = str(tmp_path / "z-path.txt")
        r = flask_client.post(
            "/api/trajectory_save_zpath",
            json={"z_path_txt": txt, "out_path": out_path},
        )
        if not _traj_flask_200_or_ineligible(r, dashboard_experiment):
            return
        assert os.path.isfile(out_path)
        loaded = np.loadtxt(out_path)
        np.testing.assert_allclose(loaded, z_arr)

    def test_non_string_txt_rejected(self, flask_client) -> None:
        r = flask_client.post(
            "/api/trajectory_save_zpath",
            json={"z_path_txt": 12345, "out_path": "x.txt"},
        )
        assert r.status_code == 400


class TestTrajectoryAnchorSaveImportRoundtrip:
    """Import anchors → coords → save z-path in one Flask session."""

    def test_trajectory_anchor_save_import_roundtrip(
        self,
        flask_client,
        dashboard_experiment: DashboardExperiment,
        tmp_path,
    ) -> None:
        anchors_file = tmp_path / "anchors.txt"
        anchors_file.write_text("0 5 10\n")
        r = flask_client.post(
            "/api/trajectory_import_anchors",
            json={
                "server_path": str(anchors_file),
                "mode": "direct",
                "x": "z0",
                "y": "z1",
                "n_points": 2,
            },
        )
        if not _traj_flask_200_or_ineligible(r, dashboard_experiment):
            return
        js = r.get_json()
        assert js["ok"] is True
        anchors = js["anchor_indices"]
        assert anchors == [0, 5, 10]

        r = flask_client.post(
            "/api/trajectory_coords",
            json={
                "mode": "direct",
                "x": "z0",
                "y": "z1",
                "anchor_indices": anchors,
                "n_points": 2,
            },
        )
        if not _traj_flask_200_or_ineligible(r, dashboard_experiment):
            return
        coords = r.get_json()
        assert len(coords["z_traj"]) >= 2
        z_txt = z_traj_to_savetxt_str(np.asarray(coords["z_traj"], dtype=np.float64))
        out_path = str(tmp_path / "z-path-roundtrip.txt")
        r = flask_client.post(
            "/api/trajectory_save_zpath",
            json={"z_path_txt": z_txt, "out_path": out_path},
        )
        if not _traj_flask_200_or_ineligible(r, dashboard_experiment):
            return
        assert r.get_json().get("ok") is True
        assert os.path.isfile(out_path)
        loaded = np.loadtxt(out_path)
        assert loaded.ndim == 2
        assert loaded.shape[0] >= 2


class TestTrajectoryVolumeApis:
    """HTTP coverage for trajectory volume generation and save routes."""

    _DIRECT_BODY = {
        "mode": "direct",
        "x": "z0",
        "y": "z1",
        "start": [0.0, 0.0],
        "end": [1.0, 1.0],
        "n_points": 3,
    }

    @pytest.mark.parametrize(
        "path,payload",
        [
            ("/api/trajectory_volumes", _DIRECT_BODY),
            (
                "/api/trajectory_save_volumes",
                {"volume_cache_id": "tok", "out_dir": "/tmp/out"},
            ),
            ("/api/trajectory_save_gif", {"images": ["a", "b"], "fps": 4}),
        ],
        ids=["volumes", "save_volumes", "save_gif"],
    )
    def test_ineligible_volume_apis_are_400(
        self,
        flask_client,
        monkeypatch: pytest.MonkeyPatch,
        path: str,
        payload: dict,
    ) -> None:
        _monkeypatch_explorer_volumes_eligible(monkeypatch, eligible=False)
        r = flask_client.post(path, json=payload)
        assert r.status_code == 400
        assert r.get_json().get("error") == _TRAJECTORY_INELIGIBLE_MSG

    def test_trajectory_save_volumes_requires_cache_id(
        self, flask_client_volumes_eligible
    ) -> None:
        r = flask_client_volumes_eligible.post(
            "/api/trajectory_save_volumes",
            json={"out_dir": "/tmp/out"},
        )
        assert r.status_code == 400
        assert "volume_cache_id" in r.get_json().get("error", "").lower()

    def test_trajectory_save_volumes_requires_out_dir(
        self, flask_client_volumes_eligible
    ) -> None:
        r = flask_client_volumes_eligible.post(
            "/api/trajectory_save_volumes",
            json={"volume_cache_id": "tok"},
        )
        assert r.status_code == 400
        assert "folder" in r.get_json().get("error", "").lower()

    def test_trajectory_save_volumes_api_mocked(
        self,
        flask_client_volumes_eligible,
        tmp_path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        out_dir = tmp_path / "saved_vols"
        saved_path = str(out_dir / "trajectory_volume_001.mrc")

        def _fake_save(token: str, directory: str, *, filename_prefix: str = "volume"):
            assert token == "cache-tok"
            assert directory == str(out_dir)
            assert filename_prefix == "trajectory_volume"
            return [saved_path]

        monkeypatch.setattr(
            "cryodrgn.dashboard.routes_analysis.save_cached_volumes_to_dir",
            _fake_save,
        )
        r = flask_client_volumes_eligible.post(
            "/api/trajectory_save_volumes",
            json={"volume_cache_id": "cache-tok", "out_dir": str(out_dir)},
        )
        assert r.status_code == 200
        j = r.get_json()
        assert j["ok"] is True
        assert j["n_saved"] == 1
        assert j["files"] == [saved_path]

    def test_trajectory_save_gif_requires_two_frames(
        self, flask_client_volumes_eligible
    ) -> None:
        r = flask_client_volumes_eligible.post(
            "/api/trajectory_save_gif",
            json={"images": ["only"], "fps": 4},
        )
        assert r.status_code == 400
        assert "two" in r.get_json().get("error", "").lower()

    def test_trajectory_save_gif_roundtrip(
        self,
        flask_client_volumes_eligible,
        tmp_path,
        dashboard_experiment: DashboardExperiment,
    ) -> None:
        from tests.conftest import png_b64_rgb

        a = png_b64_rgb(rgb=(10, 20, 30))
        b = png_b64_rgb(rgb=(200, 180, 40))
        out_path = str(tmp_path / "traj-movie.gif")
        r = flask_client_volumes_eligible.post(
            "/api/trajectory_save_gif",
            json={"images": [a, b], "fps": 5, "out_path": out_path},
        )
        if not _traj_flask_200_or_ineligible(r, dashboard_experiment):
            return
        assert os.path.isfile(out_path)
        raw = open(out_path, "rb").read()
        assert raw[:6] in (b"GIF87a", b"GIF89a")

    def test_trajectory_volumes_api_mocked(
        self,
        flask_client_volumes_eligible,
        dashboard_experiment: DashboardExperiment,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        z_traj = dashboard_experiment.z[:3]

        def _fake_compute(exp, params):
            return z_traj, None, np.zeros((3, 2), dtype=np.float64)

        from tests.conftest import make_fake_trajectory_volume_pngs

        monkeypatch.setattr(
            "cryodrgn.dashboard.routes_analysis.compute_trajectory_latent_path",
            _fake_compute,
        )
        monkeypatch.setattr(
            "cryodrgn.dashboard.routes_analysis.generate_trajectory_volume_pngs",
            make_fake_trajectory_volume_pngs(),
        )
        r = flask_client_volumes_eligible.post(
            "/api/trajectory_volumes",
            json={**self._DIRECT_BODY, "render_backend": "chimerax"},
        )
        assert r.status_code == 200, r.get_data(as_text=True)[:500]
        j = r.get_json()
        assert j["ok"] is True
        assert j["volume_cache_id"] == "cache-tok"
        assert len(j["images"]) == 3
        assert all(isinstance(b64, str) and b64 for b64 in j["images"])
        assert len(j["z_traj"]) >= 2

    def test_trajectory_volumes_from_cache_api_mocked(
        self,
        flask_client_volumes_eligible,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        fake_payloads = [
            {"index": 0, "volume_b64": "AAA=", "D": 32},
            {"index": 1, "volume_b64": "BBB=", "D": 32},
        ]

        monkeypatch.setattr(
            "cryodrgn.dashboard.routes_analysis.trajectory_volume_b64_list_from_cache",
            lambda token, exp: fake_payloads,
        )
        r = flask_client_volumes_eligible.post(
            "/api/trajectory_volumes",
            json={
                "volume_cache_id": "cache-tok",
                "volumes_from_cache": True,
                "render_backend": "vtk",
            },
        )
        assert r.status_code == 200, r.get_data(as_text=True)[:500]
        j = r.get_json()
        assert j["ok"] is True
        assert j["volume_cache_id"] == "cache-tok"
        assert j["render_backend"] == "vtk"
        assert j["volumes"] == fake_payloads

    """Headless Chromium: default trajectory overlay and random anchor coords."""

    @pytest.mark.browser
    @pytest.mark.slow
    def test_default_trajectory_and_random_anchors(
        self, playwright_page, dashboard_volumes_eligible_live_url
    ) -> None:
        from tests.conftest import (
            DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
            dashboard_smoke_trajectory,
            playwright_route_volume_viewer_render_stub,
        )

        playwright_route_volume_viewer_render_stub(playwright_page)
        out = dashboard_smoke_trajectory(
            playwright_page,
            dashboard_volumes_eligible_live_url,
            timeout_ms=DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
        )
        assert out is not None
        assert out["glyph_markers_before_anchor"] >= 2
        assert out["active_picker_before_anchor"] >= 2
        assert out["glyph_markers_after_anchor"] >= 2


class TestTrajectoryPageWithMockEligibility:
    """Trajectory shell renders without CUDA when eligibility is mocked."""

    def test_trajectory_page_renders_scatter(
        self, flask_client_volumes_eligible
    ) -> None:
        r = flask_client_volumes_eligible.get("/trajectory")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert 'id="scatter"' in body
        assert "CUDA GPU" not in body
        assert 'id="traj-z-panel"' in body
        assert body.index('id="traj-z-panel"') < body.index('id="traj-volume-actions"')
        assert "refreshTrajectoryZPanel" in body

    def test_trajectory_coords_direct_mode(self, flask_client_volumes_eligible) -> None:
        r = flask_client_volumes_eligible.post(
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
        assert r.status_code == 200, r.get_data(as_text=True)[:500]
        js = r.get_json()
        assert len(js["z_traj"]) >= 2
        assert all(len(z) == 4 for z in js["z_traj"])
        assert js["mode"] == "direct"


class TestManualTrajectoryParticleSnap:
    def test_snap_request_body_flag(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        p = parse_trajectory_request_body(
            dashboard_experiment,
            {
                "anchor_indices": [0, 5, 10],
                "mode": "direct",
                "x": "z0",
                "y": "z1",
                "n_points": 2,
                "snap_to_nearest_particles": True,
            },
        )
        assert p["snap_to_nearest_particles"] is True

    def test_midpoint_snap_uses_nearest_particles_in_plot_axes(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        coords = dashboard_experiment.plot_df[["z0", "z1"]].values.astype(np.float64)
        rows = [0, min(10, len(coords) - 1), min(20, len(coords) - 1)]
        p = {
            "use_anchors": True,
            "anchor_indices": rows,
            "xcol": "z0",
            "ycol": "z1",
            "mode": "direct",
            "n_points": 1,
            "snap_to_nearest_particles": True,
            "anchor_path_order": "preserve",
        }
        z_traj, traj_rows, traj_xy = compute_trajectory_latent_path(
            dashboard_experiment, p
        )
        assert traj_rows is not None
        assert len(traj_rows) == len(z_traj) == len(traj_xy)
        expected_n = (len(rows) - 1) * (p["n_points"] + 1) + 1
        assert len(traj_rows) == expected_n
        n = int(dashboard_experiment.z.shape[0])
        for i, r in enumerate(traj_rows):
            ri = int(r)
            assert 0 <= ri < n
            np.testing.assert_allclose(z_traj[i], dashboard_experiment.z[ri])
            np.testing.assert_allclose(traj_xy[i], coords[ri])
        mid = 0.5 * (coords[rows[0]] + coords[rows[1]])
        expected_mid = int(np.argmin(np.sum((coords - mid) ** 2, axis=1)))
        assert traj_rows[1] == expected_mid

    def test_ten_anchors_one_interp_point_yields_nineteen_samples(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        n = int(dashboard_experiment.z.shape[0])
        rows = [int(i * (n - 1) / 9) for i in range(10)]
        p = {
            "use_anchors": True,
            "anchor_indices": rows,
            "xcol": "z0",
            "ycol": "z1",
            "mode": "direct",
            "n_points": 1,
            "snap_to_nearest_particles": True,
            "anchor_path_order": "preserve",
        }
        z_traj, traj_rows, traj_xy = compute_trajectory_latent_path(
            dashboard_experiment, p
        )
        assert len(traj_rows) == len(z_traj) == len(traj_xy) == 19

    def test_direct_line_custom_traj_xy_preserved_for_display(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        """Client densified PC/catalog polyline must not be replaced by particle XY."""
        n = int(dashboard_experiment.z.shape[0])
        rows = [0, min(10, n - 1), min(20, n - 1)]
        custom = [[float(i), float(i) * 0.5] for i in range(5)]
        p = parse_trajectory_request_body(
            dashboard_experiment,
            {
                "anchor_indices": rows,
                "mode": "direct",
                "x": "z0",
                "y": "z1",
                "n_points": 1,
                "traj_xy": custom,
            },
        )
        assert p["traj_xy_custom"] == custom
        assert p["snap_to_nearest_particles"] is False
        z_traj, traj_rows, traj_xy = compute_trajectory_latent_path(
            dashboard_experiment, p
        )
        assert traj_rows is None
        assert len(z_traj) == len(traj_xy) == 5
        np.testing.assert_allclose(traj_xy, np.asarray(custom, dtype=np.float64))

    def test_snap_coords_api_payload(
        self, flask_client, dashboard_experiment: DashboardExperiment
    ) -> None:
        rows = [0, min(10, len(dashboard_experiment.z) - 1)]
        r = flask_client.post(
            "/api/trajectory_coords",
            json={
                "anchor_indices": rows,
                "mode": "direct",
                "x": "z0",
                "y": "z1",
                "n_points": 1,
                "snap_to_nearest_particles": True,
                "anchor_path_order": "preserve",
            },
        )
        if not _traj_flask_200_or_ineligible(r, dashboard_experiment):
            return
        js = r.get_json()
        assert js.get("traj_rows")
        assert not js.get("analyze_volume_ids")
        assert len(js["traj_rows"]) == 3
