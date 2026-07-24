"""Core dashboard smoke, context, session, app-shell, template contracts, and small helpers."""

from __future__ import annotations

import os
import pickle
import re
import tempfile
import numpy as np
import pytest

import yaml

from cryodrgn.dashboard import app as dash_app
from cryodrgn.dashboard import column_names
from cryodrgn.dashboard import covariate_labels
from cryodrgn.dashboard import palette_config
from cryodrgn.dashboard.context import (
    EXP_CACHE,
    EXPERIMENT_STORE,
    PRELOAD_CACHE,
    _abbrev_middle_token,
    _argv_four_command_lines,
    _cmd_argv_for_nav_display,
    _config_has_cryodrgn_cmd,
    _run_log_cryodrgn_version,
    _workdir_options,
    abbrev_middle,
    nav_interface_title,
    active_workdir,
    clear_experiment_caches,
    discover_cryodrgn_workdirs,
    epochs_for_workdir,
    resolve_epoch,
)
from cryodrgn.dashboard.data import DashboardExperiment, list_z_epochs, load_experiment
from cryodrgn.dashboard.trajectory import _TRAJ_GRAPH_NEIGHBOR_CACHE
from tests.conftest import (
    _DASHBOARD_DEFAULT_TEST_CACHE,
    _DASHBOARD_FIXTURE_SUBDIR,
    _dashboard_is_usable_workdir,
    _dashboard_resolve_fixture_workdir,
    decode_plotly_figure,
    decode_plotly_value,
    read_dashboard_template,
    torch_cuda_reports_available_but_broken,
)

pytestmark = pytest.mark.dashboard

ANALYZE_EPOCH = 2


class TestDashboardFixtureCache:
    def test_session_workdir_is_usable(self, dashboard_workdir: str) -> None:
        assert _dashboard_is_usable_workdir(dashboard_workdir)

    def test_complete_default_cache_skips_tmp_build(
        self, tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delenv("CRYODRGN_DASHBOARD_TEST_OUTDIR", raising=False)
        default = os.path.join(_DASHBOARD_DEFAULT_TEST_CACHE, _DASHBOARD_FIXTURE_SUBDIR)
        workdir, shared = _dashboard_resolve_fixture_workdir(tmp_path_factory)
        if _dashboard_is_usable_workdir(default):
            assert workdir == default
            assert shared is True
        else:
            assert shared is False
            assert workdir != default

    def test_torch_cuda_broken_gpu_probe_returns_bool(self) -> None:
        assert isinstance(torch_cuda_reports_available_but_broken(), bool)


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

    def test_load_experiment_without_umap(
        self, dashboard_workdir_plain_copy: str
    ) -> None:
        umap_path = os.path.join(
            dashboard_workdir_plain_copy,
            f"analyze.{ANALYZE_EPOCH}",
            "umap.pkl",
        )
        os.remove(umap_path)
        e = load_experiment(dashboard_workdir_plain_copy)
        assert e.umap is None
        assert "UMAP1" not in e.plot_df.columns
        assert "z0" in e.plot_df.columns


class TestLegacyRedirects:
    """Old bookmark URLs should still land on the right shell."""

    def test_abinit_builder_redirects_to_command_builder(self, flask_client) -> None:
        r = flask_client.get("/abinit-builder", follow_redirects=False)
        assert r.status_code == 302
        assert r.headers.get("Location", "").endswith("/command-builder")

    def test_filter_redirects_to_explorer(self, flask_client) -> None:
        r = flask_client.get("/filter", follow_redirects=False)
        assert r.status_code == 302
        assert r.headers.get("Location", "").endswith("/explorer")


class TestDashboardPages:
    """Every top-level Flask view should render without error."""

    @pytest.mark.parametrize(
        "path",
        [
            "/",
            "/explorer",
            "/pairplot",
            "/latent-3d",
            "/landscape-full-3d",
            "/trajectory",
            "/command-builder",
            "/landscape-volpca",
        ],
    )
    def test_page_renders(self, flask_client, path: str) -> None:
        r = flask_client.get(path)
        assert r.status_code == 200, f"{path} returned {r.status_code}"
        assert r.data, f"{path} returned empty body"

    def test_volume_viewer_url_redirects_to_trajectory(self, flask_client) -> None:
        r = flask_client.get("/volume-viewer", follow_redirects=False)
        assert r.status_code in (301, 302, 303, 307, 308)
        assert "/trajectory" in (r.headers.get("Location") or "")


class TestDashboardScatterApis:
    """2-D scatter and shared preview endpoints (3-D scatter: ``test_dashboard_scatter3d``)."""

    def test_api_scatter_json(self, flask_client) -> None:
        r = flask_client.get("/api/scatter?x=UMAP1&y=UMAP2&color=labels&filter_ui=1")
        assert r.status_code == 200
        js = r.get_json()
        assert "data" in js  # Plotly figure
        assert js["data"]

    def test_api_scatter_no_color(self, flask_client) -> None:
        r = flask_client.get("/api/scatter?x=UMAP1&y=UMAP2&color=none&marker_size=3")
        assert r.status_code == 200

    def test_api_scatter_discrete_label_colors(self, flask_client) -> None:
        import json
        from urllib.parse import quote

        overrides = json.dumps({"1": "#ff0000"})
        r = flask_client.get(
            "/api/scatter?x=UMAP1&y=UMAP2&color=labels&filter_ui=1"
            + "&discrete_label_colors="
            + quote(overrides)
        )
        assert r.status_code == 200
        js = r.get_json()
        colors = js["data"][0]["marker"]["color"]
        assert isinstance(colors, list)
        assert "#ff0000" in colors
        meta = js.get("layout", {}).get("meta", {})
        assert meta.get("cdrgn_color_legend", {}).get("type") == "discrete"
        assert meta["cdrgn_color_legend"]["items"]

    def test_api_scatter_honors_explicit_max_points(self, flask_client) -> None:
        r = flask_client.get("/api/scatter?x=UMAP1&y=UMAP2&color=none&max_points=2")
        assert r.status_code == 200
        js = decode_plotly_figure(r.get_json())
        assert len(decode_plotly_value(js["data"][0]["customdata"])) == 2

    def test_api_covariate_legend_context_landscape_scope_without_outputs_is_400(
        self, flask_client
    ) -> None:
        r = flask_client.post(
            "/api/covariate_legend_context",
            json={"column": "z0", "scope": "landscape_full_sampled"},
        )
        assert r.status_code == 400

    def test_api_preview_montage(self, flask_client) -> None:
        r = flask_client.get("/api/preview_montage?rows=0,1,2,3")
        assert r.status_code == 200
        assert r.data[:8] == b"\x89PNG\r\n\x1a\n"

    def test_api_scatter_preselect_rows_highlights_points(self, flask_client) -> None:
        r = flask_client.get(
            "/api/scatter?x=UMAP1&y=UMAP2&color=labels&preselect_rows=3,7,11"
        )
        assert r.status_code == 200
        js = r.get_json()
        meta = (js.get("layout") or {}).get("meta") or {}
        pre_trace_idx = meta.get("cdrgn_preselected")
        assert pre_trace_idx is not None and len(pre_trace_idx) == 3
        cd = decode_plotly_value(js["data"][0]["customdata"])
        rows_highlighted = {int(cd[int(i)][1]) for i in pre_trace_idx}
        assert rows_highlighted == {3, 7, 11}

    def test_api_scatter_invalid_preselect_rows_is_400(self, flask_client) -> None:
        r = flask_client.get(
            "/api/scatter?x=UMAP1&y=UMAP2&color=labels&preselect_rows=0,foo"
        )
        assert r.status_code == 400
        assert "preselect" in (r.get_json().get("error") or "").lower()

    def test_api_scatter_explorer_flag_uses_explorer_cap_path(
        self, flask_client
    ) -> None:
        r = flask_client.get(
            "/api/scatter?x=UMAP1&y=UMAP2&color=none&explorer_scatter=1"
        )
        assert r.status_code == 200
        js = decode_plotly_figure(r.get_json())
        assert len(decode_plotly_value(js["data"][0]["customdata"])) == 100

    @pytest.mark.parametrize(
        "extra_query",
        [
            "explorer_scatter=1",
            "filter_ui=1",
            "",
        ],
    )
    def test_api_scatter_main_trace_is_scattergl(
        self, flask_client, extra_query: str
    ) -> None:
        """Browser scatter must use Plotly Scattergl — SVG Scatter does not scale to explorer caps."""
        base = "x=UMAP1&y=UMAP2&color=labels"
        q = f"{base}&{extra_query}" if extra_query else base
        r = flask_client.get(f"/api/scatter?{q}")
        assert r.status_code == 200, r.get_data(as_text=True)
        js = r.get_json()
        assert js["data"], js
        assert js["data"][0].get("type") == "scattergl"

    def test_api_scatter_full_returns_entire_df(self, flask_client) -> None:
        r = flask_client.get("/api/scatter?x=UMAP1&y=UMAP2&color=none&full=1")
        assert r.status_code == 200
        js = decode_plotly_figure(r.get_json())
        assert len(decode_plotly_value(js["data"][0]["customdata"])) == 100


class TestDashboardScatterCapHelpers:
    """Caps shared by the CLI ``--filter-max`` flag and explorer scatter."""

    @pytest.mark.parametrize(
        "env_value,expected_cap,expected_from_env,func",
        [
            ("120000", 120_000, True, dash_app._particle_explorer_scatter_max_points),
            (
                "not_an_int",
                200_000,
                False,
                dash_app._particle_explorer_scatter_max_points,
            ),
            ("nope", 500_000, False, dash_app._filter_ui_scatter_max_points),
            ("800000", 800_000, True, dash_app._particle_explorer_scatter_max_points),
            (
                "30000000",
                2_000_000,
                True,
                dash_app._particle_explorer_scatter_max_points,
            ),  # clamped to max
            (
                "10",
                50_000,
                True,
                dash_app._particle_explorer_scatter_max_points,
            ),  # clamped to min
        ],
    )
    def test_scatter_cap_env_variations(
        self,
        monkeypatch: pytest.MonkeyPatch,
        env_value: str,
        expected_cap: int,
        expected_from_env: bool,
        func: callable,
    ) -> None:
        """Consolidated test for scatter cap helpers with various env values."""
        monkeypatch.setenv("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", env_value)
        assert func() == expected_cap
        monkeypatch.delenv("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", raising=False)

    def test_scatter_cap_base_helper_directly(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Test the base _scatter_cap helper that powers all cap functions."""
        from cryodrgn.dashboard.route_helpers import _scatter_cap

        # Valid env value
        monkeypatch.setenv("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", "75000")
        cap, from_env = _scatter_cap(100_000)
        assert cap == 75_000
        assert from_env is True

        # Invalid env value falls back to default
        monkeypatch.setenv("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", "invalid")
        cap, from_env = _scatter_cap(100_000)
        assert cap == 100_000
        assert from_env is False

        # Empty env uses default
        monkeypatch.delenv("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", raising=False)
        cap, from_env = _scatter_cap(100_000)
        assert cap == 100_000
        assert from_env is False

    def test_explorer_scatter_cap_from_env_indicator(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Test the boolean indicator for env override presence."""
        monkeypatch.setenv("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", "12345")
        assert dash_app._particle_explorer_scatter_cap_from_env() is True
        monkeypatch.delenv("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", raising=False)

        monkeypatch.setenv("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", "bad")
        assert dash_app._particle_explorer_scatter_cap_from_env() is False
        monkeypatch.delenv("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", raising=False)


class TestApiTryHelper:
    """Tests for the _api_try helper that standardizes API error handling."""

    @pytest.fixture
    def api_try(self):
        """Fixture to run _api_try tests within a Flask app context."""
        from cryodrgn.dashboard.route_helpers import _api_try

        app = dash_app.create_app(workdir=None)
        with app.app_context():
            yield _api_try

    def test_api_try_success_returns_result(self, api_try) -> None:
        def success_fn():
            return {"data": "success"}

        result = api_try(success_fn, "test error")
        assert result == {"data": "success"}

    def test_api_try_valueerror_returns_400(self, api_try) -> None:
        def raise_valueerror():
            raise ValueError("invalid input")

        response, status = api_try(raise_valueerror, "test error")
        assert status == 400
        assert "error" in response.get_json()
        assert "invalid input" in response.get_json()["error"]

    def test_api_try_runtimeerror_returns_500(self, api_try) -> None:
        def raise_runtimeerror():
            raise RuntimeError("server error")

        response, status = api_try(raise_runtimeerror, "test error")
        assert status == 500
        assert "error" in response.get_json()
        assert "server error" in response.get_json()["error"]

    def test_api_try_exception_returns_500_with_logging(self, api_try) -> None:
        def raise_exception():
            raise Exception("unexpected error")

        # Create a mock logger
        class MockLogger:
            def __init__(self):
                self.messages = []

            def exception(self, msg):
                self.messages.append(msg)

        mock_logger = MockLogger()
        response, status = api_try(
            raise_exception, "custom error message", logger=mock_logger
        )

        assert status == 500
        assert "error" in response.get_json()
        assert "unexpected error" in response.get_json()["error"]
        assert "custom error message" in mock_logger.messages

    def test_api_try_no_logger_doesnt_fail(self, api_try) -> None:
        def raise_exception():
            raise Exception("error without logger")

        # Should not raise even with no logger
        response, status = api_try(raise_exception, "test message")
        assert status == 500
        assert "error" in response.get_json()


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


class TestListZEpochs:
    """Epoch discovery from z.N.pkl + analyze.N directory pairs."""

    @pytest.mark.parametrize(
        "setup_func,expected",
        [
            # Missing workdir returns empty
            (lambda p: None, []),
            # z.N.pkl without matching analyze.N is skipped
            (
                lambda p: [
                    (p / "z.1.pkl").write_bytes(b""),
                    (p / "z.2.pkl").write_bytes(b""),
                    (p / "analyze.2").mkdir(),
                ],
                [2],
            ),
            # Multiple epochs are sorted
            (
                lambda p: [
                    [
                        (p / f"z.{ep}.pkl").write_bytes(b""),
                        (p / f"analyze.{ep}").mkdir(),
                    ]
                    for ep in (5, 1, 3)
                ],
                [1, 3, 5],
            ),
        ],
    )
    def test_epoch_discovery(self, tmp_path, setup_func, expected) -> None:
        """Consolidated tests for epoch listing with various directory states."""
        if setup_func is not None:
            setup_func(tmp_path)
        assert list_z_epochs(str(tmp_path)) == expected


class TestDashboardExperimentExtras:
    def test_can_preview_particles_is_true(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        assert dashboard_experiment.can_preview_particles is True

    def test_numeric_columns_exclude_index(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        cols = dashboard_experiment.numeric_columns
        assert "index" not in cols
        assert "UMAP1" in cols and "PC1" in cols


# ---------------------------------------------------------------------------
# context.py/app.py: context, session routes, and app-shell behaviour
# ---------------------------------------------------------------------------


class TestConfigHasCryodrgnCmd:
    @pytest.mark.parametrize(
        "cfg,expected",
        [
            ({"cmd": "cryodrgn train_vae ..."}, True),
            ({"cmd": ["cryodrgn", "train_vae", "..."]}, True),
            ({"cmd": "relion_refine"}, False),
            ({}, False),
            ({"cmd": None}, False),
            ("not a dict", False),
        ],
    )
    def test_detects_cryodrgn(self, cfg: object, expected: bool) -> None:
        assert _config_has_cryodrgn_cmd(cfg) is expected


class TestRunLogCryodrgnVersion:
    """``run.log`` parsing for the dashboard ``trained with cryoDRGN …`` stamp.

    The first line is often ``sys.argv`` (``.../cryodrgn train_vae ...``). The parser
    must not treat ``train_vae`` as a version string.
    """

    def test_prefers_version_line_after_argv(self, tmp_path) -> None:
        log = tmp_path / "run.log"
        log.write_text(
            "/opt/conda/bin/cryodrgn train_vae particles.txt -o out\n"
            "cryoDRGN 4.2.2.dev3+g7c39d6d22.d20260426\n",
            encoding="utf-8",
        )
        full, short, title = _run_log_cryodrgn_version(str(tmp_path))
        assert full == "4.2.2.dev3+g7c39d6d22.d20260426"
        assert short == "4.2.2.dev3"
        assert title is not None
        assert "4.2.2.dev3" in title
        assert str(log.resolve()) in title or str(log) in title

    def test_skips_subcommand_without_leading_digit(self, tmp_path) -> None:
        (tmp_path / "run.log").write_text(
            "prefix cryodrgn train_vae suffix\n",
            encoding="utf-8",
        )
        assert _run_log_cryodrgn_version(str(tmp_path)) == (None, None, None)

    def test_matches_logged_line_with_level_prefix(self, tmp_path) -> None:
        (tmp_path / "run.log").write_text(
            "INFO:__main__:cryoDRGN 3.1.0\n",
            encoding="utf-8",
        )
        full, short, title = _run_log_cryodrgn_version(str(tmp_path))
        assert full == "3.1.0"
        assert short == "3.1.0"
        assert title and "3.1.0" in title

    def test_missing_run_log_returns_none(self, tmp_path) -> None:
        assert _run_log_cryodrgn_version(str(tmp_path)) == (None, None, None)

    def test_empty_workdir_returns_none(self) -> None:
        assert _run_log_cryodrgn_version("") == (None, None, None)


class TestDiscoverCryodrgnWorkdirs:
    def test_only_returns_cryodrgn_workdirs(self, tmp_path) -> None:
        good = tmp_path / "run_vae"
        good.mkdir()
        (good / "config.yaml").write_text(
            yaml.safe_dump({"cmd": ["cryodrgn", "train_vae"]})
        )

        other_cmd = tmp_path / "relion_run"
        other_cmd.mkdir()
        (other_cmd / "config.yaml").write_text(
            yaml.safe_dump({"cmd": ["relion_refine"]})
        )

        no_cfg = tmp_path / "bare"
        no_cfg.mkdir()

        bad_yaml = tmp_path / "broken"
        bad_yaml.mkdir()
        (bad_yaml / "config.yaml").write_text(": : :\nnot valid")

        file_not_dir = tmp_path / "notdir.yaml"
        file_not_dir.write_text("cmd: cryodrgn train_vae")

        found = discover_cryodrgn_workdirs(str(tmp_path))
        assert [os.path.basename(p) for p in found] == ["run_vae"]

    def test_nonexistent_cwd_returns_empty(self, tmp_path) -> None:
        assert discover_cryodrgn_workdirs(str(tmp_path / "no_such")) == []


class TestWorkdirOptions:
    def test_relative_labels(self, tmp_path) -> None:
        (tmp_path / "a").mkdir()
        (tmp_path / "b").mkdir()
        abs_a = str((tmp_path / "a").resolve())
        abs_b = str((tmp_path / "b").resolve())
        opts = _workdir_options([abs_a, abs_b], str(tmp_path.resolve()))
        assert opts[0]["value"] == abs_a
        assert opts[0]["label"] == "a"
        assert opts[1]["label"] == "b"


class TestEpochsForWorkdir:
    def test_returns_analyzed_epochs_sorted(self, dashboard_workdir: str) -> None:
        # Re-fetch via helper; should reuse the cached list.
        epochs = epochs_for_workdir(dashboard_workdir)
        assert epochs == sorted(epochs)
        assert ANALYZE_EPOCH in epochs


class TestContextDisplayHelpers:
    @pytest.mark.parametrize(
        "text,maxlen,expected",
        [
            ("hello", 30, "hello"),
            (None, 30, ""),
            ("abcdef", 3, "abc"),
        ],
    )
    def test_abbrev_middle_simple_cases(
        self, text: str | None, maxlen: int, expected: str
    ) -> None:
        assert abbrev_middle(text, maxlen=maxlen) == expected

    def test_abbrev_middle_long_uses_middle_ellipsis(self) -> None:
        out = abbrev_middle("a" * 20 + "b" * 20, maxlen=20)
        assert len(out) == 20
        assert "\u2026" in out
        assert out.startswith("a")
        assert out.endswith("b")

    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("3D Visualizer", "3D visualizer"),
            ("PARTICLE explorer", "particle explorer"),
            ("3-D latent space", "3D latent space"),
            ("", ""),
            (None, ""),
        ],
    )
    def test_nav_interface_title(self, raw: str | None, expected: str) -> None:
        assert nav_interface_title(raw) == expected

    @pytest.mark.parametrize(
        "argv,expected",
        [
            (
                ["/usr/bin/python", "-m", "cryodrgn", "train_vae", "-o", "out"],
                ["cryodrgn", "train_vae", "-o", "out"],
            ),
            (
                ["/opt/envs/cdrgn/bin/cryodrgn", "train_vae", "-o", "out"],
                ["cryodrgn", "train_vae", "-o", "out"],
            ),
            (
                ["/usr/bin/python", "/opt/envs/cdrgn/bin/cryodrgn", "analyze", "10"],
                ["cryodrgn", "analyze", "10"],
            ),
            ([], []),
            (["bash", "foo.sh"], ["bash", "foo.sh"]),
        ],
    )
    def test_cmd_argv_for_nav_display(
        self, argv: list[str], expected: list[str]
    ) -> None:
        assert _cmd_argv_for_nav_display(argv) == expected

    @pytest.mark.parametrize(
        "argv,expected",
        [
            ([], []),
            (["cryodrgn"], ["cryodrgn"]),
        ],
    )
    def test_argv_four_command_lines_simple_cases(
        self, argv: list[str], expected: list[str]
    ) -> None:
        assert _argv_four_command_lines(argv) == expected

    def test_argv_four_command_lines_groups_head_and_abbreviates(self) -> None:
        out = _argv_four_command_lines(
            ["cryodrgn", "train_vae", "particles.mrcs", "-o", "out", "--zdim", "8"]
        )
        assert out[0] == "cryodrgn train_vae"
        assert 1 <= len(out) - 1 <= 3

        long_out = _argv_four_command_lines(
            ["cryodrgn", "train_vae", "/a/" + "x" * 300]
        )
        assert any("\u2026" in line for line in long_out)

    def test_abbrev_middle_token_keeps_short_and_abbreviates_long(self) -> None:
        assert _abbrev_middle_token("short") == "short"
        out = _abbrev_middle_token("/" + "x" * 200, maxlen=50)
        assert len(out) == 50
        assert "\u2026" in out


class TestExperimentStore:
    def test_get_experiment_caches_by_key(
        self, dashboard_workdir: str, dashboard_experiment: DashboardExperiment
    ) -> None:
        EXPERIMENT_STORE.clear_all()
        a = EXPERIMENT_STORE.get_experiment(
            dashboard_workdir, dashboard_experiment.epoch, -1
        )
        b = EXPERIMENT_STORE.get_experiment(
            dashboard_workdir, dashboard_experiment.epoch, -1
        )
        assert a is b
        assert len(EXPERIMENT_STORE.experiments) == 1

    def test_clear_preloads_for_experiment_is_scoped(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        EXPERIMENT_STORE.clear_all()
        exp = dashboard_experiment
        ep = int(exp.epoch)
        km = int(exp.kmeans_folder_id)
        PRELOAD_CACHE[(ep, km, "z0", "z1", None)] = ([1], ["img"], 0.1)
        PRELOAD_CACHE[(ep + 1, km, "z0", "z1", None)] = ([2], ["img2"], 0.2)
        n = EXPERIMENT_STORE.clear_preloads_for_experiment(exp)
        assert n == 1
        assert len(PRELOAD_CACHE) == 1
        assert PRELOAD_CACHE[(ep + 1, km, "z0", "z1", None)][0] == [2]


class TestClearExperimentCaches:
    def test_clears_all_three_caches(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        EXP_CACHE[("x", 0, -1)] = dashboard_experiment
        PRELOAD_CACHE[(0, -1, "z0", "z1", None)] = ([1], ["img"], 0.1)
        from scipy.sparse import csr_matrix

        _TRAJ_GRAPH_NEIGHBOR_CACHE[("wd", 0, 2, 2)] = (
            np.zeros((1, 1), dtype=np.int64),
            np.zeros((1, 1), dtype=np.float64),
            csr_matrix((1, 1)),
        )
        clear_experiment_caches()
        assert not EXP_CACHE
        assert not PRELOAD_CACHE
        assert not _TRAJ_GRAPH_NEIGHBOR_CACHE


class TestDiscoveryRootApp:
    def test_create_app_discovery_root_lists_child_runs(self, tmp_path) -> None:
        run = tmp_path / "run1"
        (run / "analyze.0").mkdir(parents=True)
        (run / "z.0.pkl").write_bytes(b"")
        app = dash_app.create_app(workdir=None, discovery_root=str(tmp_path))
        discovered = app.config["DASHBOARD_DISCOVERED_WORKDIRS"]
        assert os.path.abspath(str(run)) in discovered
        assert app.config["DASHBOARD_DISCOVERY_CWD"] == os.path.abspath(str(tmp_path))

    def test_create_app_discovery_root_empty_raises(self, tmp_path) -> None:
        empty = tmp_path / "empty"
        empty.mkdir()
        with pytest.raises(ValueError, match="No cryoDRGN analyzed"):
            dash_app.create_app(workdir=None, discovery_root=str(empty))


class TestActiveWorkdirAndResolveEpoch:
    def test_active_workdir_returns_dashboard_workdir(
        self, dashboard_workdir: str
    ) -> None:
        app = dash_app.create_app(workdir=dashboard_workdir)
        with app.test_request_context():
            assert active_workdir(app) == dashboard_workdir
            assert resolve_epoch(app) == ANALYZE_EPOCH

    def test_no_workdir_returns_none(self) -> None:
        app = dash_app.create_app(workdir=None)
        with app.test_request_context():
            assert active_workdir(app) is None


class TestFlaskErrorPaths:
    def test_scatter_bad_axis(self, flask_client) -> None:
        r = flask_client.get("/api/scatter?x=labels&y=notacolumn")
        assert r.status_code == 400
        assert "error" in r.get_json()

    def test_scatter_bad_color(self, flask_client) -> None:
        r = flask_client.get("/api/scatter?x=UMAP1&y=UMAP2&color=bogus")
        assert r.status_code == 400

    def test_preview_montage_non_integer_rows(self, flask_client) -> None:
        r = flask_client.get("/api/preview_montage?rows=0,foo,2")
        assert r.status_code == 400

    def test_pairplot_missing_color_col(self, flask_client) -> None:
        r = flask_client.post("/api/pairplot", json={"diagonal_emb": "pc"})
        assert r.status_code == 400

    def test_pairplot_bogus_diagonal(self, flask_client) -> None:
        r = flask_client.post(
            "/api/pairplot",
            json={
                "color_col": "labels",
                "diagonal_emb": "bogus",
                "upper_style": "scatter",
            },
        )
        assert r.status_code == 400

    def test_pairplot_bogus_upper(self, flask_client) -> None:
        r = flask_client.post(
            "/api/pairplot",
            json={
                "color_col": "labels",
                "diagonal_emb": "umap",
                "upper_style": "stack",
            },
        )
        assert r.status_code == 400

    def test_pairplot_z_as_color_rejected(self, flask_client) -> None:
        r = flask_client.post(
            "/api/pairplot",
            json={"color_col": "z0", "diagonal_emb": "umap", "upper_style": "hex"},
        )
        assert r.status_code == 400

    def test_save_selection_empty_rows(self, flask_client) -> None:
        r = flask_client.post("/api/save_selection", json={"rows": []})
        assert r.status_code == 400

    def test_save_selection_out_of_range(self, flask_client) -> None:
        r = flask_client.post("/api/save_selection", json={"rows": [0, 999999]})
        assert r.status_code == 400


class TestListServerFiles:
    def test_lists_workdir_contents(self, flask_client) -> None:
        r = flask_client.get("/api/list_server_files")
        assert r.status_code == 200
        js = r.get_json()
        assert js["ok"] is True
        assert "entries" in js

    def test_bad_dir_returns_400(self, flask_client, tmp_path) -> None:
        r = flask_client.get(
            "/api/list_server_files?dir=" + str(tmp_path / "does_not_exist")
        )
        assert r.status_code == 400


class TestSetEpochEndpoint:
    def test_post_same_epoch_succeeds(self, flask_client) -> None:
        r = flask_client.post("/api/set_epoch", json={"epoch": ANALYZE_EPOCH})
        assert r.status_code == 200
        assert r.get_json() == {"ok": True, "epoch": ANALYZE_EPOCH}

    def test_invalid_epoch_is_400(self, flask_client) -> None:
        r = flask_client.post("/api/set_epoch", json={"epoch": 99999})
        assert r.status_code == 400

    def test_non_integer_epoch_is_400(self, flask_client) -> None:
        r = flask_client.post("/api/set_epoch", json={"epoch": "abc"})
        assert r.status_code == 400

    def test_missing_epoch_is_400(self, flask_client) -> None:
        r = flask_client.post("/api/set_epoch", json={})
        assert r.status_code == 400


class TestSetWorkdirEndpoint:
    def test_invalid_workdir_rejected_by_default_app(self, flask_client) -> None:
        r = flask_client.post(
            "/api/set_workdir", json={"workdir": "/nope/does/not/exist"}
        )
        assert r.status_code == 400

    def test_clear_workdir_rejected_in_bound_mode(self, flask_client) -> None:
        # Server was started with a workdir, so clearing is not allowed.
        r = flask_client.post("/api/set_workdir", json={"workdir": ""})
        assert r.status_code == 400

    def test_valid_workdir_switch_in_builder_only_mode(
        self, dashboard_workdir: str
    ) -> None:
        app = dash_app.create_app(workdir=None)
        # Pretend the session discovered our fixture workdir.
        app.config["DASHBOARD_DISCOVERED_WORKDIRS"] = [dashboard_workdir]
        with app.test_client() as client:
            r = client.post("/api/set_workdir", json={"workdir": dashboard_workdir})
            assert r.status_code == 200, r.get_json()
            js = r.get_json()
            assert js["ok"] is True
            assert js["workdir"] == dashboard_workdir

    def test_clear_in_builder_only_mode_succeeds(self, dashboard_workdir: str) -> None:
        app = dash_app.create_app(workdir=None)
        app.config["DASHBOARD_DISCOVERED_WORKDIRS"] = [dashboard_workdir]
        with app.test_client() as client:
            r = client.post("/api/set_workdir", json={"workdir": ""})
            assert r.status_code == 200


class TestChimeraXPathApi:
    def test_set_chimerax_path_rejected_when_boot_env_set(
        self, dashboard_workdir: str, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("CHIMERAX_PATH", "/opt/chimerax/bin/ChimeraX")
        app = dash_app.create_app(workdir=dashboard_workdir)
        with app.test_client() as c:
            r = c.post("/api/set_chimerax_path", json={"path": "/tmp/should-not-apply"})
            assert r.status_code == 400

    def test_set_chimerax_path_accepts_executable(
        self, dashboard_workdir: str, monkeypatch: pytest.MonkeyPatch, tmp_path
    ) -> None:
        monkeypatch.delenv("CHIMERAX_PATH", raising=False)
        fake = tmp_path / "fake_chimerax"
        fake.write_text("#!/bin/sh\nexit 0\n")
        fake.chmod(0o755)
        app = dash_app.create_app(workdir=dashboard_workdir)
        try:
            with app.test_client() as c:
                r = c.post("/api/set_chimerax_path", json={"path": str(fake)})
                assert r.status_code == 200
                j = r.get_json()
                assert j is not None
                assert j.get("ok") is True
                resolved = os.path.realpath(str(fake))
                assert j.get("path") == resolved
                assert os.environ.get("CHIMERAX_PATH") == resolved
        finally:
            os.environ.pop("CHIMERAX_PATH", None)

    def test_set_chimerax_path_rejects_non_executable(
        self, dashboard_workdir: str, monkeypatch: pytest.MonkeyPatch, tmp_path
    ) -> None:
        monkeypatch.delenv("CHIMERAX_PATH", raising=False)
        fake = tmp_path / "not_run"
        fake.write_text("not a script")
        fake.chmod(0o644)
        app = dash_app.create_app(workdir=dashboard_workdir)
        with app.test_client() as c:
            r = c.post("/api/set_chimerax_path", json={"path": str(fake)})
            assert r.status_code == 400


class TestRoutesTableIntegrity:
    def test_every_entry_is_callable(self) -> None:
        for rule, view_func, methods in dash_app._ROUTES:
            assert callable(view_func), rule
            # Methods can be a single string or a tuple of strings
            if isinstance(methods, str):
                assert methods, rule
            else:
                assert isinstance(methods, tuple) and methods, rule

    def test_create_app_registers_every_route(self, dashboard_workdir: str) -> None:
        app = dash_app.create_app(workdir=dashboard_workdir)
        rules = {r.rule for r in app.url_map.iter_rules()}
        for rule, _func, _methods in dash_app._ROUTES:
            assert rule in rules, f"missing rule: {rule}"


class TestIndexTemplateNavLinks:
    def test_index_landscape_is_inactive_without_analyze_landscape(
        self, flask_client_no_landscape
    ) -> None:
        r = flask_client_no_landscape.get("/")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "/landscape-volpca" not in body
        assert "Volume sketched landscape explorer" in body
        assert "analyze_landscape" in body


class TestCommandBuilderOnlyMode:
    def test_builder_only_index_renders(self, dashboard_workdir: str) -> None:
        app = dash_app.create_app(workdir=None)
        with app.test_client() as client:
            r = client.get("/")
            assert r.status_code == 200
            r2 = client.get("/command-builder")
            assert r2.status_code == 200

    def test_builder_only_explorer_redirects_home(self) -> None:
        app = dash_app.create_app(workdir=None)
        with app.test_client() as client:
            r = client.get("/explorer", follow_redirects=False)
            assert r.status_code == 302
            assert r.headers["Location"].endswith("/")


class TestDiscreteCovariateLegendContracts:
    """Guards for discrete covariate toggle legends (CryoColorCovariateLegend hosts).

    The pair-plot generator once advanced ``lastPayload`` before a new ``blob:`` PNG was
    actually shown (manual ``el.onload()`` right after ``src`` while ``complete`` still
    described the old bitmap), so the first discrete toggle looked like a no-op.
    """

    _HOST_TEMPLATES = (
        "pair_grid.html",
        "latent_3d.html",
        "_particle_explorer_montage_load_save.html",
    )

    def _read_template(self, rel: str) -> str:
        return read_dashboard_template(rel)

    @pytest.mark.parametrize("rel", _HOST_TEMPLATES)
    def test_discrete_legend_hosts_set_notify_on_refresh_false(self, rel: str) -> None:
        """``refresh()`` must not fire ``onFilterChange``; hosts redraw via ``afterLayout`` / UI."""
        text = self._read_template(rel)
        marker = "new CryoColorCovariateLegend({"
        starts = [m.start() for m in re.finditer(re.escape(marker), text)]
        assert (
            len(starts) == 1
        ), f"{rel}: expected exactly one {marker!r}, got {len(starts)}"
        pos = starts[0]
        window = text[pos : pos + 30_000]
        assert "notifyOnRefresh: false" in window, (
            f"{rel}: CryoColorCovariateLegend options must set notifyOnRefresh: false "
            "within ~30k chars of the constructor (avoids duplicate redraws vs lastPayload)."
        )

    def test_pairplot_no_synthetic_img_onload_after_blob_src(self) -> None:
        text = self._read_template("pair_grid.html")
        bad = re.search(
            r"el\.src\s*=\s*nextUrl\s*;[\s\S]{0,400}?"
            r"el\.complete\s*&&\s*el\.naturalWidth[\s\S]{0,120}?"
            r"el\.onload\s*\(\s*\)",
            text,
        )
        assert bad is None, (
            "pair_grid.html: do not synthesize img onload from complete/naturalWidth right "
            "after assigning a new blob src — lastPayload can desync from the visible bitmap."
        )


class TestChimeraxAnimation:
    """Unified ChimeraX rotation GIF helpers (``cryodrgn.dashboard.chimerax_animation``)."""

    def test_explorer_rotation_gif_frames_default(self) -> None:
        from cryodrgn.dashboard.chimerax_animation import explorer_rotation_gif_frames

        assert explorer_rotation_gif_frames(8) == 16
        assert explorer_rotation_gif_frames(1) == 4
        assert explorer_rotation_gif_frames(32) == 64
        assert explorer_rotation_gif_frames(100) == 64

    def test_rotation_turn_y_increments_constant_speed(self) -> None:
        from cryodrgn.dashboard.chimerax_animation import rotation_turn_y_increments

        assert rotation_turn_y_increments(4) == [0.0, 90.0, 90.0, 90.0]
        assert rotation_turn_y_increments(8) == [0.0] + [45.0] * 7
        inc = rotation_turn_y_increments(16)
        assert inc[0] == 0.0
        assert inc[1:] == pytest.approx([22.5] * 15)

    def test_parallel_jobs_caps_at_eight(self) -> None:
        from cryodrgn.dashboard.chimerax_animation import (
            landscape_rotate_parallel_jobs,
            parallel_jobs,
        )

        assert parallel_jobs(16, 12) == 8
        assert parallel_jobs(4, 12) == 4
        assert landscape_rotate_parallel_jobs(16, 12) == 8

    def test_landscape_rotate_render_strategy_hybrid(self) -> None:
        from cryodrgn.dashboard.chimerax_animation import (
            landscape_rotate_render_strategy,
        )

        assert landscape_rotate_render_strategy(1, 16, 8) == "serial_session"
        assert landscape_rotate_render_strategy(3, 16, 8) == "parallel_session"
        assert landscape_rotate_render_strategy(2, 8, 4) == "serial_session"

    def test_landscape_rotate_parallel_jobs_capped_at_eight(self) -> None:
        from cryodrgn.dashboard.chimerax_animation import (
            landscape_rotate_parallel_jobs,
        )

        assert landscape_rotate_parallel_jobs(16, 5) == 5
        assert landscape_rotate_parallel_jobs(16, 12) == 8
        assert landscape_rotate_parallel_jobs(4, 12) == 4

    def test_chimerax_animation_meta(self) -> None:
        from cryodrgn.dashboard.chimerax_animation import chimerax_animation_meta

        meta = chimerax_animation_meta(8)
        assert meta["chimerax_cpus_effective"] == 8
        assert meta["chimerax_parallel_volume_cap"] == 8

    def test_rotation_session_cmds_single_exit_incremental_turns(self) -> None:
        from cryodrgn.dashboard.chimerax_animation import (
            chimerax_rotation_session_cmds,
            rotation_turn_y_increments,
        )

        turns = rotation_turn_y_increments(4)
        frames = [(f"/tmp/f{i}.png", float(turns[i]), None) for i in range(4)]
        cmds, matrix_log = chimerax_rotation_session_cmds(
            "/tmp/x.mrc",
            frames,
            100,
            volume_color="#336699",
        )
        assert cmds.count("exit") == 1
        assert sum(1 for c in cmds if c.startswith("open ")) == 1
        assert sum(1 for c in cmds if c.startswith("save ")) == 4
        assert turns[1:] == [90.0, 90.0, 90.0]
        assert any("turn y 90" in c for c in cmds)
        assert "surface dust all size 10" in cmds
        assert "lighting soft" in cmds
        assert matrix_log is None

    def test_render_rotating_gif_dispatches_to_session(self, monkeypatch) -> None:
        from cryodrgn.dashboard.chimerax_animation import render_rotating_gif

        called: list[str] = []

        def _fake_session(*_a, **_k):
            called.append("session")
            return None

        monkeypatch.setattr(
            "cryodrgn.dashboard.chimerax_animation.render_rotating_gif_single_session",
            _fake_session,
        )
        render_rotating_gif("/tmp/x.mrc", "/tmp/out.gif", gif_frames=8, ncpus=4)
        assert called == ["session"]

    def test_use_chimerax_xvfb_defaults_and_env(self, monkeypatch) -> None:
        from cryodrgn.dashboard import chimerax_animation as cx

        monkeypatch.delenv("CRYODRGN_CHIMERAX_XVFB", raising=False)
        monkeypatch.delenv("DISPLAY", raising=False)
        assert cx.use_chimerax_xvfb() is True
        monkeypatch.setenv("DISPLAY", ":0")
        assert cx.use_chimerax_xvfb() is False
        monkeypatch.setenv("CRYODRGN_CHIMERAX_XVFB", "1")
        assert cx.use_chimerax_xvfb() is True
        monkeypatch.setenv("CRYODRGN_CHIMERAX_XVFB", "0")
        assert cx.use_chimerax_xvfb() is False

    def test_run_chimerax_cmds_uses_xvfb_when_headless(self, monkeypatch) -> None:
        from cryodrgn.dashboard import chimerax_animation as cx

        monkeypatch.setenv("CHIMERAX_PATH", "/opt/ChimeraX")
        monkeypatch.delenv("CRYODRGN_CHIMERAX_XVFB", raising=False)
        monkeypatch.delenv("DISPLAY", raising=False)
        monkeypatch.setattr(cx, "_ensure_shared_xvfb_display", lambda: ":91")
        seen: list[dict] = []

        class _Proc:
            returncode = 0
            stdout = "ok"
            stderr = ""

        def _fake_run(cmd, **kwargs):
            seen.append({"cmd": cmd, "env": kwargs.get("env")})
            return _Proc()

        monkeypatch.setattr(cx.subprocess, "run", _fake_run)
        out, err = cx.run_chimerax_cmds(["open /tmp/a.mrc", "exit"])
        assert out == "ok" and err == ""
        assert len(seen) == 1
        assert "--offscreen" not in seen[0]["cmd"]
        assert "/opt/ChimeraX" in seen[0]["cmd"]
        assert seen[0]["env"] is not None
        assert seen[0]["env"].get("DISPLAY") == ":91"

    def test_run_chimerax_cmds_retries_xvfb_after_opengl_failure(
        self, monkeypatch
    ) -> None:
        from cryodrgn.dashboard import chimerax_animation as cx

        monkeypatch.setenv("CHIMERAX_PATH", "/opt/ChimeraX")
        monkeypatch.delenv("CRYODRGN_CHIMERAX_XVFB", raising=False)
        monkeypatch.setenv("DISPLAY", ":99")
        monkeypatch.setattr(cx, "_ensure_shared_xvfb_display", lambda: ":91")
        seen: list[dict] = []

        class _Proc:
            def __init__(self, code: int, err: str):
                self.returncode = code
                self.stdout = ""
                self.stderr = err

        def _fake_run(cmd, **kwargs):
            seen.append({"cmd": cmd, "env": kwargs.get("env")})
            if "--offscreen" in cmd:
                return _Proc(
                    0,
                    "LimitationError: Unable to save images because "
                    "OpenGL rendering is not available",
                )
            return _Proc(0, "")

        monkeypatch.setattr(cx.subprocess, "run", _fake_run)
        cx.run_chimerax_cmds(["save /tmp/x.png", "exit"])
        assert len(seen) == 2
        assert "--offscreen" in seen[0]["cmd"]
        assert "--offscreen" not in seen[1]["cmd"]
        assert seen[1]["env"] is not None
        assert seen[1]["env"].get("DISPLAY") == ":91"


class TestDashboardModules:
    """Imports, pure helpers, and landscape-full helpers (folded from ``test_dashboard_modules``)."""

    def test_dashboard_import_order_avoids_plots_route_helpers_cycle(self) -> None:
        """``route_helpers`` must not depend on ``plots`` (historical circular import)."""
        import cryodrgn.dashboard.route_helpers as rh  # noqa: PLC0415

        assert (
            rh.normalize_continuous_palette
            is palette_config.normalize_continuous_palette
        )
        import cryodrgn.dashboard.plots as plots  # noqa: PLC0415

        assert plots.scatter_json is not None
        import cryodrgn.dashboard.landscape_full_3d as lf3  # noqa: PLC0415

        assert (
            lf3.VOL_LANDSCAPE_3D_PLOT_DF_ROW
            == column_names.VOL_LANDSCAPE_3D_PLOT_DF_ROW
        )

    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("landscape_vol_PC12", "Vol PC12"),
            ("landscape_vol_cluster", "Vol cluster"),
            ("landscape_vol_UMAP1", "Vol UMAP1"),
            ("landscape_vol_umap2", "Vol UMAP2"),
        ],
    )
    def test_covariate_display_names(self, raw: str, expected: str) -> None:
        assert covariate_labels.covariate_display_name(raw) == expected

    def test_landscape_vol_pc_labels_and_display_map_include_variance(self) -> None:
        evr = np.array([0.453, 0.12], dtype=np.float64)
        assert covariate_labels.landscape_vol_pc_pretty_label(1, evr) == (
            "Vol PC1 (45.3%)"
        )
        assert covariate_labels.landscape_vol_pc_pretty_label(2, evr) == (
            "Vol PC2 (12.0%)"
        )
        assert covariate_labels.landscape_vol_pc_pretty_label(3, evr) == "Vol PC3"

        display = covariate_labels.covariate_display_map(
            ["landscape_vol_PC1", "z0"],
            vol_pc_explained_variance_ratio=np.array([0.1], dtype=np.float64),
        )
        assert display["landscape_vol_PC1"] == "Vol PC1 (10.0%)"
        assert display["z0"] == "z0"

    @pytest.mark.parametrize(
        "raw,expected,cmap",
        [
            (None, "Viridis", "viridis"),
            ("", "Viridis", "viridis"),
            ("viridis", "Viridis", "viridis"),
            ("PLASMA", "Plasma", "plasma"),
            ("TURBO", "Turbo", "turbo"),
            ("  Plasma  ", "Plasma", "plasma"),
            ("not_a_real_palette", "Viridis", "viridis"),
        ],
    )
    def test_continuous_palette_normalization(
        self, raw: str | None, expected: str, cmap: str
    ) -> None:
        assert palette_config.normalize_continuous_palette(raw) == expected
        assert palette_config.mpl_cmap_for_palette(expected) == cmap

    def test_vol_landscape_plot_df_row_constant_matches_column_names(self) -> None:
        assert column_names.VOL_LANDSCAPE_3D_PLOT_DF_ROW == "_dashboard_plot_df_row"

    def test_merge_covariate_pkl_1d(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        from cryodrgn.dashboard.covariate_pkl import merge_covariate_pkl

        e = dashboard_experiment
        n = len(e.plot_df)
        path = os.path.join(e.workdir, "test_covariate_1d.pkl")
        with open(path, "wb") as fh:
            pickle.dump(np.arange(n, dtype=np.float64), fh)
        cols, discrete = merge_covariate_pkl(e, path)
        assert cols == ["test_covariate_1d"]
        assert discrete == []
        assert "test_covariate_1d" in e.numeric_columns
        assert e.plot_df["test_covariate_1d"].shape == (n,)

    def test_merge_covariate_pkl_wrong_length(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        from cryodrgn.dashboard.covariate_pkl import merge_covariate_pkl

        e = dashboard_experiment
        path = os.path.join(e.workdir, "bad_len.pkl")
        with open(path, "wb") as fh:
            pickle.dump(np.zeros(3, dtype=np.float64), fh)
        with pytest.raises(ValueError, match="does not match particle count"):
            merge_covariate_pkl(e, path)

    def test_merge_covariate_pkl_not_array(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        from cryodrgn.dashboard.covariate_pkl import merge_covariate_pkl

        e = dashboard_experiment
        path = os.path.join(e.workdir, "not_array.pkl")
        with open(path, "wb") as fh:
            pickle.dump({"a": 1}, fh)
        with pytest.raises(ValueError, match="Expected an array"):
            merge_covariate_pkl(e, path)

    def test_merge_covariate_pkl_string_labels(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        from cryodrgn.dashboard.covariate_pkl import merge_covariate_pkl
        from cryodrgn.dashboard.plots_color_covariate import (
            _lower_color_series_is_discrete,
        )

        e = dashboard_experiment
        n = len(e.plot_df)
        labels = np.array(
            ["open", "closed"] * (n // 2) + ["open"] * (n % 2), dtype=object
        )
        path = os.path.join(e.workdir, "state_labels.pkl")
        with open(path, "wb") as fh:
            pickle.dump(labels, fh)
        cols, discrete = merge_covariate_pkl(e, path)
        assert cols == ["state_labels"]
        assert discrete == ["state_labels"]
        assert "state_labels" in e.color_covariate_columns
        assert "state_labels" not in e.numeric_columns
        assert _lower_color_series_is_discrete(e.plot_df["state_labels"])

    def test_merge_covariate_pkl_string_labels_2d(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        from cryodrgn.dashboard.covariate_pkl import merge_covariate_pkl

        e = dashboard_experiment
        n = len(e.plot_df)
        arr = np.column_stack(
            [
                np.array(["a", "b"] * (n // 2) + ["a"] * (n % 2), dtype=object),
                np.array(["x", "y"] * (n // 2) + ["x"] * (n % 2), dtype=object),
            ]
        )
        path = os.path.join(e.workdir, "multi_labels.pkl")
        with open(path, "wb") as fh:
            pickle.dump(arr, fh)
        cols, discrete = merge_covariate_pkl(e, path)
        assert cols == ["multi_labels_0", "multi_labels_1"]
        assert discrete == cols
        assert all(c in e.user_covariate_columns for c in cols)

    def test_api_load_covariate_pkl_string_labels(
        self, flask_client, dashboard_experiment: DashboardExperiment
    ) -> None:
        e = dashboard_experiment
        n = len(e.plot_df)
        path = os.path.join(e.workdir, "api_labels.pkl")
        with open(path, "wb") as fh:
            pickle.dump(
                np.array(["type_a", "type_b"] * (n // 2) + ["type_a"] * (n % 2)), fh
            )
        r = flask_client.post("/api/load_covariate_pkl", json={"path": path})
        assert r.status_code == 200
        j = r.get_json()
        assert j["ok"] is True
        assert j["discrete_columns"] == ["api_labels"]
        assert "api_labels" in j["color_cols"]
        r_scatter = flask_client.get(
            "/api/scatter?x=UMAP1&y=UMAP2&color=api_labels&filter_ui=1"
        )
        assert r_scatter.status_code == 200

    def test_api_load_covariate_pkl(
        self, flask_client, dashboard_experiment: DashboardExperiment
    ) -> None:
        e = dashboard_experiment
        n = len(e.plot_df)
        path = os.path.join(e.workdir, "api_cov.pkl")
        with open(path, "wb") as fh:
            pickle.dump(np.linspace(0, 1, n), fh)
        r = flask_client.post("/api/load_covariate_pkl", json={"path": path})
        assert r.status_code == 200
        j = r.get_json()
        assert j["ok"] is True
        assert j["primary_column"] == "api_cov"
        assert "api_cov" in j["numeric_cols"]

    def test_api_load_covariate_pkl_outside_workdir(
        self, flask_client, dashboard_experiment: DashboardExperiment
    ) -> None:
        e = dashboard_experiment
        n = len(e.plot_df)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "outside_cov.pkl")
            with open(path, "wb") as fh:
                pickle.dump(np.linspace(0, 1, n), fh)
            assert not os.path.abspath(path).startswith(
                os.path.abspath(e.workdir) + os.sep
            )
            r = flask_client.post("/api/load_covariate_pkl", json={"path": path})
            assert r.status_code == 200
            j = r.get_json()
            assert j["ok"] is True
            assert j["primary_column"] == "outside_cov"
            assert "outside_cov" in j["numeric_cols"]

    def test_api_load_covariate_pkl_bad_format(
        self, flask_client, dashboard_experiment: DashboardExperiment
    ) -> None:
        e = dashboard_experiment
        path = os.path.join(e.workdir, "api_bad.pkl")
        with open(path, "wb") as fh:
            pickle.dump([1, 2, 3], fh)
        r = flask_client.post("/api/load_covariate_pkl", json={"path": path})
        assert r.status_code == 400
        assert "error" in r.get_json()

    def test_api_list_server_files_lists_pkl(
        self, flask_client, dashboard_experiment: DashboardExperiment
    ) -> None:
        e = dashboard_experiment
        path = os.path.join(e.workdir, "listed_cov.pkl")
        with open(path, "wb") as fh:
            pickle.dump(np.zeros(len(e.plot_df)), fh)
        r = flask_client.get("/api/list_server_files?kinds=pkl&dir=" + e.workdir)
        assert r.status_code == 200
        names = [
            ent["name"] for ent in r.get_json()["entries"] if ent["type"] == "file"
        ]
        assert "listed_cov.pkl" in names


class TestDashboardLandscapeHelpers:
    def test_landscape_full_ready_false_without_outputs(self) -> None:
        from cryodrgn.dashboard import landscape_full_3d  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as tmp:
            assert not landscape_full_3d.landscape_full_ready(tmp, 0)
            assert not landscape_full_3d.landscape_full_3d_ready(tmp, 0)

    def test_landscape_full_numeric_covariates_skip_animation_helper_columns(
        self,
    ) -> None:
        import pandas as pd

        from cryodrgn.dashboard import landscape_full_3d as lf3
        from cryodrgn.dashboard.column_names import (
            VOL_LANDSCAPE_3D_PLOT_DF_ROW,
            VOL_LANDSCAPE_IS_SKETCH_CENTROID,
            VOL_LANDSCAPE_NEAREST_SKETCH_VOL,
        )

        df = pd.DataFrame(
            {
                "index": [0, 1],
                "landscape_vol_PC1": [0.0, 1.0],
                VOL_LANDSCAPE_3D_PLOT_DF_ROW: [10, 11],
                VOL_LANDSCAPE_NEAREST_SKETCH_VOL: [2, 3],
            }
        )
        cols = lf3.landscape_full_sampled_numeric_covariates(df)
        assert "landscape_vol_PC1" in cols
        assert VOL_LANDSCAPE_NEAREST_SKETCH_VOL not in cols
        assert VOL_LANDSCAPE_IS_SKETCH_CENTROID not in cols
        assert VOL_LANDSCAPE_3D_PLOT_DF_ROW not in cols


class TestDecodePlotlyTypedArrays:
    """``decode_plotly_value`` matches Plotly 6 binary trace blobs in tests."""

    def test_decode_plotly_value_accepts_string_and_list_shape(self) -> None:
        import base64

        arr = np.array([[1, 2], [3, 4]], dtype=np.int32)
        b64 = base64.standard_b64encode(arr.tobytes()).decode("ascii")
        blob_str = {"dtype": "i4", "bdata": b64, "shape": "2,2"}
        blob_list = {"dtype": "i4", "bdata": b64, "shape": [2, 2]}
        assert decode_plotly_value(blob_str) == [[1, 2], [3, 4]]
        assert decode_plotly_value(blob_list) == [[1, 2], [3, 4]]


class TestBundledPlotlyJs:
    """Plotly.js is served from the installed Python package (not a CDN)."""

    def test_vendor_route_returns_js(self) -> None:
        from cryodrgn.dashboard.bundled_plotly import bundled_plotly_js_version

        app = dash_app.create_app(workdir=None)
        with app.test_client() as client:
            r = client.get("/vendor/plotly.min.js")
        assert r.status_code == 200
        assert "javascript" in (r.content_type or "")
        assert len(r.data) > 1_000_000
        assert b"Plotly" in r.data
        assert r.headers.get("X-Plotly-Version") == bundled_plotly_js_version()

    def test_command_builder_page_references_vendor_plotly(self) -> None:
        app = dash_app.create_app(workdir=None)
        with app.test_client() as client:
            r = client.get("/command-builder")
        assert r.status_code == 200
        html = r.data.decode("utf-8", errors="replace")
        assert "cdn.plot.ly" not in html
        assert "/vendor/plotly.min.js" in html

    @pytest.mark.parametrize("path", ["/explorer", "/latent-3d"])
    def test_dashboard_pages_reference_vendor_plotly(
        self, flask_client, path: str
    ) -> None:
        r = flask_client.get(path)
        assert r.status_code == 200
        html = r.data.decode("utf-8", errors="replace")
        assert "cdn.plot.ly" not in html
        assert "/vendor/plotly.min.js" in html


class TestDashboardSessionIntegrationFlows:
    """Multi-step Flask flows spanning session, covariates, and discovery."""

    def test_session_workdir_epoch_roundtrip(self, dashboard_workdir: str) -> None:
        from cryodrgn.dashboard.context import EXPERIMENT_STORE

        app = dash_app.create_app(workdir=None)
        app.config["DASHBOARD_DISCOVERED_WORKDIRS"] = [dashboard_workdir]
        with app.test_client() as client:
            r = client.post("/api/set_workdir", json={"workdir": dashboard_workdir})
            assert r.status_code == 200, r.get_json()
            assert r.get_json()["ok"] is True

            r = client.post("/api/set_epoch", json={"epoch": ANALYZE_EPOCH})
            assert r.status_code == 200
            assert r.get_json() == {"ok": True, "epoch": ANALYZE_EPOCH}

            EXPERIMENT_STORE.clear_all()
            r = client.get("/api/scatter?x=UMAP1&y=UMAP2&color=none")
            assert r.status_code == 200
            assert r.get_json()["data"]
            assert EXPERIMENT_STORE.experiments

    def test_discovery_root_workdir_switching(
        self, dashboard_workdir: str, tmp_path
    ) -> None:
        import shutil

        parent = tmp_path / "experiments"
        parent.mkdir()
        run_a = parent / "run_a"
        run_b = parent / "run_b"
        shutil.copytree(dashboard_workdir, run_a)
        shutil.copytree(dashboard_workdir, run_b)

        app = dash_app.create_app(workdir=None, discovery_root=str(parent))
        discovered = app.config["DASHBOARD_DISCOVERED_WORKDIRS"]
        assert os.path.abspath(str(run_a)) in discovered
        assert os.path.abspath(str(run_b)) in discovered

        with app.test_client() as client:
            r = client.post(
                "/api/set_workdir", json={"workdir": os.path.abspath(str(run_a))}
            )
            assert r.status_code == 200, r.get_json()
            r_scatter = client.get("/api/scatter?x=z0&y=z1&color=none")
            assert r_scatter.status_code == 200
            n_a = len(r_scatter.get_json()["data"][0].get("x") or [])

            r = client.post(
                "/api/set_workdir", json={"workdir": os.path.abspath(str(run_b))}
            )
            assert r.status_code == 200, r.get_json()
            r_scatter_b = client.get("/api/scatter?x=z0&y=z1&color=none")
            assert r_scatter_b.status_code == 200
            n_b = len(r_scatter_b.get_json()["data"][0].get("x") or [])
            assert n_a == n_b
            assert n_b > 0

    def test_covariate_filter_save_selection_flow(
        self, flask_client, dashboard_experiment: DashboardExperiment, tmp_path
    ) -> None:
        e = dashboard_experiment
        n = len(e.plot_df)
        path = os.path.join(e.workdir, "flow_cov.pkl")
        with open(path, "wb") as fh:
            pickle.dump(np.linspace(0, 1, n), fh)

        r = flask_client.post("/api/load_covariate_pkl", json={"path": path})
        assert r.status_code == 200, r.get_json()
        assert r.get_json()["primary_column"] == "flow_cov"

        r = flask_client.post(
            "/api/covariate_threshold_rows",
            json={"column": "flow_cov", "level": 0.5, "use_max": False},
        )
        assert r.status_code == 200, r.get_json()
        rows = r.get_json()["rows"]
        assert isinstance(rows, list) and rows

        dest = tmp_path / "sel_out"
        dest.mkdir()
        r = flask_client.post(
            "/api/save_selection",
            json={
                "rows": rows[: min(5, len(rows))],
                "basename": "flow_sel",
                "sel_dir": str(dest),
            },
        )
        assert r.status_code == 200, r.get_json()
        saved = dest / "flow_sel.pkl"
        assert saved.is_file()
        with open(saved, "rb") as fh:
            loaded = pickle.load(fh)
        assert loaded.size >= 1


class TestDashboardIndexBrowserSmoke:
    """Landing page navigation cards."""

    pytestmark = pytest.mark.browser

    def test_landing_cards_nav_and_trajectory_gating(
        self, playwright_page, dashboard_live_url
    ) -> None:
        from tests.conftest import dashboard_smoke_index

        out = dashboard_smoke_index(playwright_page, dashboard_live_url)
        assert out["landing_links"] >= 4
        assert out["mentions_volume_landscapes"]
        if out["has_trajectory"]:
            pytest.skip("trajectory card active on this runner (GPU + weights present)")
        assert out["trajectory_ineligible_note"]

    @pytest.mark.parametrize(
        "fixture_url,expect_key",
        [
            ("dashboard_volumes_eligible_live_url", "has_trajectory"),
            ("dashboard_landscape_volpca_live_url", "has_landscape_volpca"),
        ],
        ids=["trajectory_eligible", "landscape_volpca"],
    )
    def test_landing_feature_cards_when_outputs_present(
        self, fixture_url: str, expect_key: str, request, playwright_page
    ) -> None:
        from tests.conftest import dashboard_smoke_index

        url = request.getfixturevalue(fixture_url)
        out = dashboard_smoke_index(playwright_page, url)
        assert out[expect_key]
