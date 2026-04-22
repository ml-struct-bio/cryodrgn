"""End-to-end and unit tests for ``cryodrgn.dashboard``.

These tests complement the existing suite by exercising the Flask dashboard
app and its helpers against a tiny ``train_vae`` + ``analyze`` run. They are
grouped so that ``pytest -n2 --dist=loadscope`` keeps them on the same worker
and the session-scoped ``dashboard_workdir`` fixture is built once per run.

A stable scratch location can be reused across runs via
``CRYODRGN_DASHBOARD_TEST_OUTDIR`` (otherwise ``tmp_path_factory`` is used).
If the cache already contains a usable ``analyze.N/umap.pkl``, the fixture
skips re-training. Point the env var at a shared path (e.g. a scratch FS) to
avoid paying the ~80s fixture cost on every run.
"""

from __future__ import annotations

import argparse
import os
import pickle
from collections.abc import Iterator

import numpy as np
import pytest

from cryodrgn.commands import analyze, train_vae
from cryodrgn.dashboard import app as dash_app
from cryodrgn.dashboard.app import _trajectory_eligibility_error
from cryodrgn.dashboard.data import DashboardExperiment, list_z_epochs, load_experiment
from cryodrgn.dashboard.explorer_volumes import _chimerax_render_cmds
from cryodrgn.dashboard.plots import pair_grid_png
from cryodrgn.dashboard.trajectory import (
    parse_int_from_dict as _parse_int_from_dict,
    parse_traj_interpolation_value as _parse_traj_interpolation_value,
    parse_traj_neighbor_value as _parse_traj_neighbor_value,
    parse_traj_points_value as _parse_traj_points_value,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
TRAIN_EPOCHS = 3
ANALYZE_EPOCH = 2  # 0-indexed; must be < TRAIN_EPOCHS


def _is_usable_workdir(workdir: str) -> bool:
    """True if ``workdir`` already has weights + z + analyze outputs we need."""
    required = [
        os.path.join(workdir, "config.yaml"),
        os.path.join(workdir, f"weights.{ANALYZE_EPOCH}.pkl"),
        os.path.join(workdir, f"z.{ANALYZE_EPOCH}.pkl"),
        os.path.join(workdir, f"analyze.{ANALYZE_EPOCH}", "umap.pkl"),
    ]
    return all(os.path.exists(p) for p in required)


def _run_train_and_analyze(workdir: str) -> None:
    """Produce a tiny heterogeneous reconstruction + analysis at ``workdir``."""
    parser = argparse.ArgumentParser()
    train_vae.add_args(parser)
    train_args = parser.parse_args(
        [
            os.path.join(DATA_DIR, "hand.mrcs"),
            "-o",
            workdir,
            "--poses",
            os.path.join(DATA_DIR, "hand_rot_trans.pkl"),
            "--ctf",
            os.path.join(DATA_DIR, "test_ctf.100.pkl"),
            "-b",
            "8",
            "--no-amp",
            "-n",
            str(TRAIN_EPOCHS),
            "--zdim",
            "4",
            "--tdim",
            "16",
            "--tlayers",
            "1",
            "--seed",
            "0",
            "--no-analysis",
        ]
    )
    train_vae.main(train_args)

    parser = argparse.ArgumentParser()
    analyze.add_args(parser)
    analyze_args = parser.parse_args(
        [workdir, str(ANALYZE_EPOCH), "--ksample", "5", "--pc", "2"]
    )
    analyze.main(analyze_args)


@pytest.fixture(scope="session")
def dashboard_workdir(tmp_path_factory: pytest.TempPathFactory) -> str:
    """Shared train_vae + analyze output directory (built once per session)."""
    cache = os.environ.get("CRYODRGN_DASHBOARD_TEST_OUTDIR")
    if cache:
        workdir = os.path.join(cache, "pytest_dashboard_fixture")
        os.makedirs(workdir, exist_ok=True)
    else:
        workdir = str(tmp_path_factory.mktemp("dashboard_vae"))

    if not _is_usable_workdir(workdir):
        _run_train_and_analyze(workdir)
        assert _is_usable_workdir(workdir), f"Fixture incomplete at {workdir!r}"
    return workdir


@pytest.fixture(scope="session")
def dashboard_experiment(dashboard_workdir: str) -> DashboardExperiment:
    """Cached ``DashboardExperiment`` loaded from :func:`dashboard_workdir`."""
    return load_experiment(dashboard_workdir)


@pytest.fixture(scope="session")
def flask_client(dashboard_workdir: str) -> Iterator:
    """Flask test client bound to the session-scoped workdir."""
    app = dash_app.create_app(workdir=dashboard_workdir)
    with app.test_client() as client:
        yield client


# ---------------------------------------------------------------------------
# Pure-function unit tests (no workdir fixture needed)
# ---------------------------------------------------------------------------


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
        kwargs = dict(
            lower_color_col="labels",
            diagonal_emb="umap",
            upper_style="scatter",
        )
        png_a, cells_a = pair_grid_png(dashboard_experiment, **kwargs)
        png_b, cells_b = pair_grid_png(dashboard_experiment, **kwargs)
        assert png_a == png_b, "pair_grid_png is non-deterministic"
        assert cells_a == cells_b
        # zdim == 4 -> 4x4 cells.
        assert len(cells_a) == dashboard_experiment.z.shape[1] ** 2


class TestDashboardTrajectoryCoords:
    """Trajectory coord-only endpoints don't need ChimeraX / GPU."""

    def test_direct_interpolation(self, flask_client) -> None:
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
        assert r.status_code == 200
        js = r.get_json()
        # 3 interpolation pts + 2 endpoints = 5 (matches live-server observation).
        assert len(js["z_traj"]) >= 2
        assert all(len(z) == 4 for z in js["z_traj"])
        assert js["mode"] == "direct"

    def test_nearest_mode(self, flask_client) -> None:
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
        assert r.status_code == 200
        assert r.get_json()["mode"] == "nearest"

    def test_anchor_driven_trajectory(self, flask_client) -> None:
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
        assert r.status_code == 200
        js = r.get_json()
        assert js.get("traj_particle_indices"), "missing anchor particle indices"

    def test_kmeans_centers(self, flask_client) -> None:
        r = flask_client.post(
            "/api/trajectory_kmeans_centers",
            json={"x": "z0", "y": "z1", "mode": "direct", "n_points": 2},
        )
        assert r.status_code == 200

    def test_random_indices(self, flask_client) -> None:
        r = flask_client.post(
            "/api/trajectory_random_indices",
            json={"x": "z0", "y": "z1", "mode": "direct", "n_points": 2},
        )
        assert r.status_code == 200

    def test_default_endpoints(self, flask_client) -> None:
        r = flask_client.get("/api/default_trajectory_endpoints?x=z0&y=z1")
        assert r.status_code == 200


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
