"""Tests for the volume sketched landscape (vol PCA) dashboard UI and APIs.

Uses the same session-scoped ``dashboard_workdir`` fixture as
``test_dashboard_core`` / ``conftest`` (epoch **2** analyze outputs).

Parallel runs: safe under ``pytest -n2 --dist=loadscope`` (and similar). The
shared training output is only read while copying; the writable tree lives
under ``tmp_path_factory`` (per-worker basetemp). ``conftest`` already
serializes building ``dashboard_workdir`` with a file lock across xdist
workers. Landscape-specific fixtures are **module**-scoped so loadscope can
keep dependent tests on one worker and we only **copytree** once per module
per worker; HTTP clients stay **function**-scoped so Flask sessions do not
leak between tests.

``tmp_path_factory.mktemp()`` creates an **empty** directory; ``shutil.copytree``
requires the destination leaf not to exist, so copies always land in
``<mktemp>/cryo_out``. Optional ``landscape.*`` trees are stripped from copies
when tests assume analyze_landscape was never run (golden fixtures may include
them).
"""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import numpy as np
import pytest

from cryodrgn import utils
from cryodrgn.dashboard import app as dash_app
from cryodrgn.dashboard.data import DashboardExperiment, load_experiment
from cryodrgn.dashboard.landscape_volpca import (
    kmeans_sorted_vol_indices,
    landscape_analysis_ready,
    landscape_dir_for_epoch,
    landscape_volpca_scatter_json,
    list_landscape_epochs,
    meta_for_api,
)

# Must match ``DASHBOARD_ANALYZE_EPOCH`` in ``tests/conftest.py``.
_ANALYZE_EPOCH = 2
_LANDSCAPE_K = 3


def _strip_landscape_dirs(workdir: str) -> None:
    """Remove ``landscape.N`` trees so copies match a run without analyze_landscape."""
    for name in list(os.listdir(workdir)):
        if name.startswith("landscape.") and os.path.isdir(os.path.join(workdir, name)):
            shutil.rmtree(os.path.join(workdir, name), ignore_errors=True)


def _copy_dashboard_workdir(
    dashboard_workdir: str,
    tmp_path_factory: pytest.TempPathFactory,
    parent_name: str,
) -> str:
    """Copy into ``parent_name/cryo_out`` (``mktemp`` only creates the parent — *not* the leaf)."""
    parent = tmp_path_factory.mktemp(parent_name)
    dst = Path(parent) / "cryo_out"
    shutil.copytree(dashboard_workdir, dst)
    return str(dst)


def _write_minimal_landscape(workdir: str, epoch: int = _ANALYZE_EPOCH) -> str:
    """Create ``landscape.{epoch}/vol_pca_K.pkl`` + ``kmeansK/`` sketch files.

    Returns the absolute ``landscape_dir`` path.
    """
    land = os.path.join(workdir, f"landscape.{epoch}")
    km = os.path.join(land, f"kmeans{_LANDSCAPE_K}")
    os.makedirs(km, exist_ok=True)
    pc = np.array(
        [[0.0, 0.1], [1.0, -0.5], [0.2, 0.3]],
        dtype=np.float64,
    )
    utils.save_pkl(pc, os.path.join(land, f"vol_pca_{_LANDSCAPE_K}.pkl"))
    for i in range(1, _LANDSCAPE_K + 1):
        p = os.path.join(km, f"vol_{i:03d}.mrc")
        with open(p, "wb"):
            pass
    centers_path = os.path.join(km, "centers_ind.txt")
    with open(centers_path, "w", encoding="utf-8") as fh:
        for row in range(_LANDSCAPE_K):
            fh.write(f"{row}\n")
    umap_full = np.array(
        [[0.0, 1.0], [2.0, 3.0], [4.0, 5.0], [6.0, 7.0], [8.0, 9.0]],
        dtype=np.float64,
    )
    utils.save_pkl(umap_full, os.path.join(land, "umap.pkl"))
    return land


@pytest.fixture(scope="module")
def dashboard_workdir_with_landscape(
    dashboard_workdir: str,
    tmp_path_factory: pytest.TempPathFactory,
) -> str:
    """Copy the shared dashboard output once per module (per xdist worker)."""
    d = _copy_dashboard_workdir(
        dashboard_workdir, tmp_path_factory, "workdir_landscape"
    )
    _strip_landscape_dirs(d)
    _write_minimal_landscape(d)
    return d


@pytest.fixture(scope="module")
def dashboard_workdir_plain_copy(
    dashboard_workdir: str,
    tmp_path_factory: pytest.TempPathFactory,
) -> str:
    """Like the session workdir but never includes ``landscape.*`` (golden tree may)."""
    d = _copy_dashboard_workdir(
        dashboard_workdir, tmp_path_factory, "wd_plain_no_landscape"
    )
    _strip_landscape_dirs(d)
    return d


@pytest.fixture(scope="module")
def experiment_landscape(
    dashboard_workdir_with_landscape: str,
) -> DashboardExperiment:
    return load_experiment(dashboard_workdir_with_landscape)


@pytest.fixture(scope="module")
def experiment_plain_no_landscape(
    dashboard_workdir_plain_copy: str,
) -> DashboardExperiment:
    return load_experiment(dashboard_workdir_plain_copy)


@pytest.fixture
def flask_client_landscape(dashboard_workdir_with_landscape: str):
    """New Flask client per test (avoid shared cookie/session state)."""
    app = dash_app.create_app(workdir=dashboard_workdir_with_landscape)
    with app.test_client() as client:
        yield client


@pytest.fixture
def flask_client_no_landscape(dashboard_workdir_plain_copy: str):
    """App bound to a tree with no ``landscape.N`` (not the session golden path)."""
    app = dash_app.create_app(workdir=dashboard_workdir_plain_copy)
    with app.test_client() as client:
        yield client


class TestLandscapeVolpcaKmeansHelpers:
    """``vol_mean.mrc`` must not become a fake vol 000 index."""

    def test_kmeans_sorted_vol_indices_excludes_vol_mean(self, tmp_path) -> None:
        d = tmp_path / "kmeans100"
        d.mkdir()
        (d / "vol_mean.mrc").write_bytes(b"x")
        (d / "vol_001.mrc").write_bytes(b"x")
        (d / "vol_100.mrc").write_bytes(b"x")
        assert kmeans_sorted_vol_indices(str(d)) == [1, 100]


class TestLandscapeVolpcaFilesystemHelpers:
    def test_list_landscape_epochs_empty(self, tmp_path) -> None:
        assert list_landscape_epochs(str(tmp_path)) == []

    def test_list_landscape_epochs_finds_folder(
        self, dashboard_workdir_with_landscape: str
    ) -> None:
        assert list_landscape_epochs(dashboard_workdir_with_landscape) == [
            _ANALYZE_EPOCH
        ]

    def test_landscape_analysis_ready(
        self, dashboard_workdir_with_landscape: str
    ) -> None:
        assert landscape_analysis_ready(
            dashboard_workdir_with_landscape,
            _ANALYZE_EPOCH,
        )

    def test_landscape_dir_for_epoch(
        self, dashboard_workdir_with_landscape: str
    ) -> None:
        d = landscape_dir_for_epoch(dashboard_workdir_with_landscape, _ANALYZE_EPOCH)
        assert d.endswith(f"landscape.{_ANALYZE_EPOCH}")
        assert os.path.isdir(d)


class TestLandscapeVolpcaMetaAndScatter:
    def test_meta_ok(self, experiment_landscape: DashboardExperiment) -> None:
        m = meta_for_api(experiment_landscape)
        assert m["ok"] is True
        assert m["landscape_epoch"] == _ANALYZE_EPOCH
        assert m["n_volumes"] == _LANDSCAPE_K
        assert m["n_pc"] == 2
        assert m["n_umap"] == 2
        assert m["has_state_color"] is False
        assert m["kmeans_k"] == _LANDSCAPE_K
        opts = m["color_options"]
        assert any(o["value"] == "none" for o in opts)
        assert any(o["value"] == "state" for o in opts)

    def test_meta_no_landscape_dirs(
        self, experiment_plain_no_landscape: DashboardExperiment
    ) -> None:
        m = meta_for_api(experiment_plain_no_landscape)
        assert m["ok"] is False
        assert "No landscape" in (m.get("error") or "")

    def test_scatter_json_marker_size_matches_sketch_constant(
        self, experiment_landscape: DashboardExperiment
    ) -> None:
        """Regression: vol-PCA scatter JSON round-trips and markers match ``_VOLSKETCH_SCATTER_MARKER``."""
        from cryodrgn.dashboard import landscape_volpca as lv

        land = landscape_dir_for_epoch(
            experiment_landscape.workdir,
            _ANALYZE_EPOCH,
        )
        js = landscape_volpca_scatter_json(
            land,
            experiment_landscape,
            axis_x="pc:0",
            axis_y="pc:1",
            color_mode="none",
            continuous_palette=None,
        )
        fig = json.loads(js)
        trace = fig["data"][0]
        assert trace["type"] == "scattergl"
        assert "ids" in trace
        assert len(trace["ids"]) == _LANDSCAPE_K
        sz = trace["marker"]["size"]
        assert float(sz) == pytest.approx(
            float(lv._VOLSKETCH_SCATTER_MARKER), rel=0, abs=1e-9
        )
        assert float(sz) == pytest.approx(9.0 * (1.0 - 0.13) * 1.3, rel=0, abs=1e-9)

    def test_scatter_rejects_bad_axes(
        self, experiment_landscape: DashboardExperiment
    ) -> None:
        land = landscape_dir_for_epoch(
            experiment_landscape.workdir,
            _ANALYZE_EPOCH,
        )
        with pytest.raises(ValueError, match="distinct"):
            landscape_volpca_scatter_json(
                land,
                experiment_landscape,
                axis_x="pc:0",
                axis_y="pc:0",
                color_mode="none",
            )

    def test_scatter_umap_axes(self, experiment_landscape: DashboardExperiment) -> None:
        land = landscape_dir_for_epoch(
            experiment_landscape.workdir,
            _ANALYZE_EPOCH,
        )
        js = landscape_volpca_scatter_json(
            land,
            experiment_landscape,
            axis_x="umap:0",
            axis_y="umap:1",
            color_mode="none",
        )
        fig = json.loads(js)
        trace = fig["data"][0]
        assert trace["x"] == [0.0, 2.0, 4.0]
        assert trace["y"] == [1.0, 3.0, 5.0]


class TestLandscapeVolpcaFlaskRoutes:
    def test_meta_route(self, flask_client_landscape) -> None:
        r = flask_client_landscape.get("/api/landscape_volpca/meta")
        assert r.status_code == 200
        m = r.get_json()
        assert m["ok"] is True
        assert m["n_volumes"] == _LANDSCAPE_K

    def test_scatter_route(self, flask_client_landscape) -> None:
        r = flask_client_landscape.get(
            "/api/landscape_volpca/scatter",
            query_string={"axis_x": "pc:0", "axis_y": "pc:1", "color": "none"},
        )
        assert r.status_code == 200
        fig = json.loads(r.get_data(as_text=True))
        assert fig["data"][0]["type"] == "scattergl"

    def test_scatter_route_legacy_pc_params(self, flask_client_landscape) -> None:
        r = flask_client_landscape.get(
            "/api/landscape_volpca/scatter",
            query_string={"pc_x": 0, "pc_y": 1, "color": "none"},
        )
        assert r.status_code == 200
        fig = json.loads(r.get_data(as_text=True))
        assert fig["data"][0]["type"] == "scattergl"

    def test_scatter_route_bad_color(self, flask_client_landscape) -> None:
        r = flask_client_landscape.get(
            "/api/landscape_volpca/scatter",
            query_string={"pc_x": 0, "pc_y": 1, "color": "__not_a_column__"},
        )
        assert r.status_code == 400
        assert "error" in (r.get_json() or {})

    def test_page_renders_with_landscape(self, flask_client_landscape) -> None:
        r = flask_client_landscape.get("/landscape-volpca")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        # ``volsketch-grid`` also appears in always-loaded CSS; require the live grid div.
        assert '<div class="volsketch-grid">' in body
        assert 'id="volsketch"' in body

    def test_page_no_landscape_message(self, flask_client_no_landscape) -> None:
        r = flask_client_no_landscape.get("/landscape-volpca")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "analyze_landscape" in body
        assert "No landscape.N folders" in body
        assert '<div class="volsketch-grid">' not in body
        assert 'id="volsketch"' not in body

    def test_index_links_landscape_when_ready(self, flask_client_landscape) -> None:
        r = flask_client_landscape.get("/")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "/landscape-volpca" in body

    def test_scatter_requires_both_axis_params(self, flask_client_landscape) -> None:
        r = flask_client_landscape.get(
            "/api/landscape_volpca/scatter",
            query_string={"axis_x": "pc:0", "color": "none"},
        )
        assert r.status_code == 400
        assert "axis_y" in (r.get_json() or {}).get("error", "")

    def test_scatter_state_color_unavailable_is_400(
        self, flask_client_landscape
    ) -> None:
        r = flask_client_landscape.get(
            "/api/landscape_volpca/scatter",
            query_string={"axis_x": "pc:0", "axis_y": "pc:1", "color": "state"},
        )
        assert r.status_code == 400
        assert "state" in (r.get_json() or {}).get("error", "").lower()

    def test_generate_animations_requires_volumes(self, flask_client_landscape) -> None:
        r = flask_client_landscape.post(
            "/api/landscape_volpca/generate_animations",
            json={},
        )
        assert r.status_code == 400

    def test_generate_animations_rejects_bad_view_rotations(
        self, flask_client_landscape
    ) -> None:
        r = flask_client_landscape.post(
            "/api/landscape_volpca/generate_animations",
            json={"vol_indices": [1], "view_rotations": [1, 2, 3]},
        )
        assert r.status_code == 400

    def test_save_animations_requires_token(self, flask_client_landscape) -> None:
        r = flask_client_landscape.post(
            "/api/landscape_volpca/save_animations",
            json={},
        )
        assert r.status_code == 400

    def test_save_animations_rejects_bad_out_dir(self, flask_client_landscape) -> None:
        r = flask_client_landscape.post(
            "/api/landscape_volpca/save_animations",
            json={"token": "nope", "out_dir": 123},
        )
        assert r.status_code == 400
