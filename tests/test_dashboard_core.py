"""Core dashboard smoke, context, session, app-shell, template contracts, and small helpers."""

from __future__ import annotations

import base64
import os
import pickle
import re
import tempfile
from io import BytesIO
from pathlib import Path

import numpy as np
import pytest
import yaml
from PIL import Image

from cryodrgn.dashboard import app as dash_app
from cryodrgn.dashboard import column_names
from cryodrgn.dashboard import covariate_labels
from cryodrgn.dashboard import palette_config
from cryodrgn.dashboard.context import (
    EXP_CACHE,
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
from cryodrgn.dashboard.data import DashboardExperiment, list_z_epochs
from cryodrgn.dashboard.particle_explorer import explorer_volumes_eligible
from cryodrgn.dashboard.plot_gif_utils import png_base64_frames_to_gif_bytes
from cryodrgn.dashboard.trajectory import _TRAJ_GRAPH_NEIGHBOR_CACHE

ANALYZE_EPOCH = 2


def _png_b64_rgb(
    w: int = 8, h: int = 8, rgb: tuple[int, int, int] = (200, 30, 40)
) -> str:
    """Solid-colour PNG as standard base64 (GIF helpers + latent-3D GIF API tests)."""
    buf = BytesIO()
    Image.new("RGB", (w, h), rgb).save(buf, format="PNG")
    return base64.standard_b64encode(buf.getvalue()).decode("ascii")


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

    def test_api_scatter_honors_explicit_max_points(self, flask_client) -> None:
        r = flask_client.get("/api/scatter?x=UMAP1&y=UMAP2&color=none&max_points=2")
        assert r.status_code == 200
        js = r.get_json()
        assert len(js["data"][0]["customdata"]) == 2

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

    def test_api_covariate_legend_context_landscape_scope_without_outputs_is_400(
        self, flask_client
    ) -> None:
        r = flask_client.post(
            "/api/covariate_legend_context",
            json={"column": "z0", "scope": "landscape_full_sampled"},
        )
        assert r.status_code == 400

    def test_api_latent3d_preview_png(self, flask_client) -> None:
        r = flask_client.get("/api/latent3d_preview.png?x=z0&y=z1&z=z2&color=znorm")
        assert r.status_code == 200
        # PNG magic bytes.
        assert r.data[:8] == b"\x89PNG\r\n\x1a\n"

    def test_api_latent3d_plot_gif_from_png_frames(self, flask_client) -> None:
        a = _png_b64_rgb(rgb=(10, 20, 30))
        b = _png_b64_rgb(rgb=(200, 180, 40))
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
        cd = js["data"][0]["customdata"]
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
        js = r.get_json()
        assert len(js["data"][0]["customdata"]) == 100

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
        js = r.get_json()
        assert len(js["data"][0]["customdata"]) == 100


class TestPlotGifUtils:
    """``plot_gif_utils.png_base64_frames_to_gif_bytes`` (PNG-frame GIF assembly)."""

    def test_png_base64_frames_to_gif_bytes_round_trip(self) -> None:
        a = _png_b64_rgb(rgb=(200, 30, 40))
        b = _png_b64_rgb(rgb=(30, 180, 60))
        gif = png_base64_frames_to_gif_bytes([a, b], durations_ms=[50, 120])
        assert gif[:6] in (b"GIF87a", b"GIF89a")
        im = Image.open(BytesIO(gif))
        assert im.n_frames >= 2
        im.seek(1)

    def test_png_base64_frames_requires_two_frames(self) -> None:
        one = _png_b64_rgb()
        with pytest.raises(ValueError, match="At least two"):
            png_base64_frames_to_gif_bytes([one])

    def test_png_base64_accepts_data_url_prefix(self) -> None:
        raw = _png_b64_rgb()
        framed = "data:image/png;base64," + raw
        gif = png_base64_frames_to_gif_bytes([framed, raw], durations_ms=40)
        assert len(gif) > 32


class TestDashboardScatterCapHelpers:
    """Caps shared by the CLI ``--filter-max`` flag and explorer scatter."""

    def test_explorer_scatter_cap_reads_filter_max_env(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", "120000")
        assert dash_app._particle_explorer_scatter_max_points() == 120_000
        assert dash_app._particle_explorer_scatter_cap_from_env() is True
        monkeypatch.delenv("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", raising=False)

    def test_explorer_scatter_cap_invalid_env_falls_back_to_default(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", "not_an_int")
        assert dash_app._particle_explorer_scatter_max_points() == 200_000
        assert dash_app._particle_explorer_scatter_cap_from_env() is False
        monkeypatch.delenv("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", raising=False)

    def test_filter_ui_cap_invalid_env_uses_half_million_default(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", "nope")
        assert dash_app._filter_ui_scatter_max_points() == 500_000
        monkeypatch.delenv("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", raising=False)


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
    def test_missing_workdir_returns_empty(self, tmp_path) -> None:
        assert list_z_epochs(str(tmp_path / "missing")) == []

    def test_requires_matching_analyze_dir(self, tmp_path) -> None:
        # ``z.N.pkl`` without ``analyze.N/`` should be skipped.
        (tmp_path / "z.1.pkl").write_bytes(b"")
        (tmp_path / "z.2.pkl").write_bytes(b"")
        (tmp_path / "analyze.2").mkdir()
        assert list_z_epochs(str(tmp_path)) == [2]

    def test_multi_epoch_sorted(self, tmp_path) -> None:
        for ep in (5, 1, 3):
            (tmp_path / f"z.{ep}.pkl").write_bytes(b"")
            (tmp_path / f"analyze.{ep}").mkdir()
        assert list_z_epochs(str(tmp_path)) == [1, 3, 5]


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


class TestAbbrevMiddle:
    def test_short_unchanged(self) -> None:
        assert abbrev_middle("hello", 30) == "hello"

    def test_long_uses_middle_ellipsis(self) -> None:
        s = "a" * 20 + "b" * 20
        out = abbrev_middle(s, maxlen=20)
        assert len(out) == 20
        assert "\u2026" in out
        assert out.startswith("a")
        assert out.endswith("b")

    def test_none_returns_empty(self) -> None:
        assert abbrev_middle(None) == ""

    def test_small_maxlen_truncates_plainly(self) -> None:
        assert abbrev_middle("abcdef", maxlen=3) == "abc"


class TestNavInterfaceTitle:
    def test_lowercase_and_3d_token(self) -> None:
        assert nav_interface_title("3D Visualizer") == "3D visualizer"
        assert nav_interface_title("PARTICLE explorer") == "particle explorer"

    def test_hyphenated_three_d(self) -> None:
        assert nav_interface_title("3-D latent space") == "3D latent space"

    def test_empty(self) -> None:
        assert nav_interface_title("") == ""
        assert nav_interface_title(None) == ""


class TestAbbrevMiddleToken:
    def test_short_unchanged(self) -> None:
        assert _abbrev_middle_token("short") == "short"

    def test_long_has_ellipsis(self) -> None:
        s = "/" + "x" * 200
        out = _abbrev_middle_token(s, maxlen=50)
        assert len(out) == 50
        assert "\u2026" in out


class TestCmdArgvForNavDisplay:
    def test_python_m_cryodrgn(self) -> None:
        assert _cmd_argv_for_nav_display(
            ["/usr/bin/python", "-m", "cryodrgn", "train_vae", "-o", "out"]
        ) == ["cryodrgn", "train_vae", "-o", "out"]

    def test_entrypoint_first(self) -> None:
        assert _cmd_argv_for_nav_display(
            ["/opt/envs/cdrgn/bin/cryodrgn", "train_vae", "-o", "out"]
        ) == ["cryodrgn", "train_vae", "-o", "out"]

    def test_python_wrapper_second(self) -> None:
        assert _cmd_argv_for_nav_display(
            ["/usr/bin/python", "/opt/envs/cdrgn/bin/cryodrgn", "analyze", "10"]
        ) == ["cryodrgn", "analyze", "10"]

    def test_empty_returns_empty(self) -> None:
        assert _cmd_argv_for_nav_display([]) == []

    def test_unrecognised_is_unchanged(self) -> None:
        assert _cmd_argv_for_nav_display(["bash", "foo.sh"]) == ["bash", "foo.sh"]


class TestArgvFourCommandLines:
    def test_empty(self) -> None:
        assert _argv_four_command_lines([]) == []

    def test_single_token(self) -> None:
        assert _argv_four_command_lines(["cryodrgn"]) == ["cryodrgn"]

    def test_head_is_two_tokens(self) -> None:
        out = _argv_four_command_lines(
            ["cryodrgn", "train_vae", "particles.mrcs", "-o", "out", "--zdim", "8"]
        )
        assert out[0] == "cryodrgn train_vae"
        # Remaining lines are non-empty and at most 3 more lines.
        assert 1 <= len(out) - 1 <= 3

    def test_long_token_is_abbreviated(self) -> None:
        long_path = "/a/" + "x" * 300
        out = _argv_four_command_lines(["cryodrgn", "train_vae", long_path])
        # Abbreviated somewhere in the rendered argv.
        assert any("\u2026" in line for line in out)


class TestClearExperimentCaches:
    def test_clears_all_three_caches(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        EXP_CACHE[("x", 0, -1)] = dashboard_experiment
        PRELOAD_CACHE[(0, -1, "z0", "z1", None)] = ([1], ["img"], 0.1)
        _TRAJ_GRAPH_NEIGHBOR_CACHE[("wd", 0, 2, 2)] = (
            np.zeros((1, 1), dtype=np.int64),
            np.zeros((1, 1), dtype=np.float64),
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

    def test_scatter3d_bad_color(self, flask_client) -> None:
        r = flask_client.get("/api/scatter3d_z?x=z0&y=z1&z=z2&color=does_not_exist")
        assert r.status_code == 400

    def test_latent3d_non_numeric_elev(self, flask_client) -> None:
        r = flask_client.get("/api/latent3d_preview.png?x=z0&y=z1&z=z2&elev=nope")
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


class TestRoutesTableIntegrity:
    def test_every_entry_is_callable(self) -> None:
        for rule, view_func, methods in dash_app._ROUTES:
            assert callable(view_func), rule
            assert isinstance(methods, tuple) and methods

    def test_create_app_registers_every_route(self, dashboard_workdir: str) -> None:
        app = dash_app.create_app(workdir=dashboard_workdir)
        rules = {r.rule for r in app.url_map.iter_rules()}
        for rule, _func, _methods in dash_app._ROUTES:
            assert rule in rules, f"missing rule: {rule}"


class TestIndexTemplateNavLinks:
    def test_index_has_expected_nav_links(
        self, flask_client, dashboard_experiment: DashboardExperiment
    ) -> None:
        r = flask_client.get("/")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        for link in (
            "/explorer",
            "/pairplot",
            "/latent-3d",
            "/command-builder",
        ):
            assert link in body, f"nav link {link!r} missing from /"
        assert "3D volume landscapes" in body
        if explorer_volumes_eligible(dashboard_experiment):
            assert "/trajectory" in body, "nav link '/trajectory' missing from /"
        else:
            assert "CUDA-enabled machine" in body
            assert "Trajectory creator" in body
            assert "/trajectory" not in body

    def test_index_landscape_is_inactive_without_analyze_landscape(
        self, flask_client
    ) -> None:
        r = flask_client.get("/")
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

    _REPO_ROOT = Path(__file__).resolve().parents[1]
    _DASH_TEMPLATES = _REPO_ROOT / "cryodrgn" / "dashboard" / "templates"
    _HOST_TEMPLATES = (
        "pair_grid.html",
        "latent_3d.html",
        "_particle_explorer_montagejs.html",
    )

    def _read_template(self, rel: str) -> str:
        p = self._DASH_TEMPLATES / rel
        assert p.is_file(), f"missing template {p}"
        return p.read_text(encoding="utf-8")

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

    @pytest.mark.parametrize("rel", _HOST_TEMPLATES)
    def test_discrete_legend_hosts_wire_discrete_dom(self, rel: str) -> None:
        text = self._read_template(rel)
        assert "discreteSwitches:" in text
        assert "discreteWrap:" in text

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

    def test_color_covariate_legend_fits_discrete_switch_column_widths(self) -> None:
        js = (
            self._REPO_ROOT
            / "cryodrgn"
            / "dashboard"
            / "static"
            / "js"
            / "color_covariate_legend.js"
        ).read_text(encoding="utf-8")
        assert "fitDiscreteSwitchColumnWidths" in js
        assert "0.79" in js
        assert "--cryo-discrete-col-w" in js
        assert "--cryo-discrete-cell-max-w" in js
        assert 'removeProperty("--cryo-discrete-cell-min-h")' in js
        assert "auto-fill" in js or "scrollWidth" in js
        css = (
            self._REPO_ROOT
            / "cryodrgn"
            / "dashboard"
            / "static"
            / "css"
            / "cryo_cc_legend_palette_menu.css"
        ).read_text(encoding="utf-8")
        assert "80svh" in css
        assert "overflow-y: auto" in css
        assert "justify-content: center" in css
        assert "--cryo-discrete-grid-gap-y: 0.24rem" in css
        assert "auto-fill" in css
        assert "--cryo-discrete-col-w" in css
        assert "cryo-cc-discrete-panel--collapsed" in css
        assert ":not(.cryo-cc-discrete-toggle-heading-row)" in css
        assert (
            "#color-discrete-wrap.cryo-cc-discrete-wrap--show.cryo-cc-discrete-panel--collapsed"
            in css
        )
        assert "#pair-color-discrete-wrap.cryo-cc-discrete-wrap--show" in css
        assert "cryo-cc-discrete-toggle-title" in css
        assert "setDiscretePanelCollapsed" in js
        assert "_ensureDiscreteCollapseToggle" in js
        assert "_upgradeDiscreteCollapseTitleButton" in js
        assert "_discretePanelCollapsed" in js
        pe_css = (
            self._REPO_ROOT
            / "cryodrgn"
            / "dashboard"
            / "static"
            / "css"
            / "particle_explorer.css"
        ).read_text(encoding="utf-8")
        montage = (
            self._REPO_ROOT
            / "cryodrgn"
            / "dashboard"
            / "templates"
            / "_particle_explorer_montage.html"
        ).read_text(encoding="utf-8")
        montagejs = (
            self._REPO_ROOT
            / "cryodrgn"
            / "dashboard"
            / "templates"
            / "_particle_explorer_montagejs.html"
        ).read_text(encoding="utf-8")
        assert "cryo-explorer-discrete-hides-scatter-palette-radios" in pe_css
        assert (
            "#scatter-palette-radios.cryo-palette-select--options-pane {\n"
            "    display: none !important;"
        ) in pe_css
        assert "cryo-cc-discrete-toggle-title" in montage
        assert "Toggle selection" in montage
        assert "K-means legend" not in montagejs

    def test_color_covariate_legend_refresh_respects_notify_on_refresh(self) -> None:
        js = (
            self._REPO_ROOT
            / "cryodrgn"
            / "dashboard"
            / "static"
            / "js"
            / "color_covariate_legend.js"
        ).read_text(encoding="utf-8")
        assert re.search(
            r"if\s*\(\s*self\.notifyOnRefresh\s*&&\s*!\s*suppressNotify\s*\)\s*self\._notify\s*\(",
            js,
        ), "finishLegendLayout must keep notifyOnRefresh + suppressNotify guard"

    def test_pairplot_draw_still_assigns_blob_src(self) -> None:
        text = self._read_template("pair_grid.html")
        assert "el.src = nextUrl" in text
        assert "el.onload = function()" in text


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

    def test_vol_landscape_plot_df_row_constant_matches_column_names(self) -> None:
        assert column_names.VOL_LANDSCAPE_3D_PLOT_DF_ROW == "_dashboard_plot_df_row"

    def test_covariate_display_name_landscape_vol_pc(self) -> None:
        assert (
            covariate_labels.covariate_display_name("landscape_vol_PC12") == "Vol PC12"
        )

    def test_covariate_display_name_landscape_vol_cluster(self) -> None:
        assert (
            covariate_labels.covariate_display_name("landscape_vol_cluster")
            == "Vol cluster"
        )

    def test_landscape_vol_pc_pretty_label_with_variance(self) -> None:
        evr = np.array([0.453, 0.12], dtype=np.float64)
        assert (
            covariate_labels.landscape_vol_pc_pretty_label(1, evr) == "Vol PC1 (45.3%)"
        )
        assert (
            covariate_labels.landscape_vol_pc_pretty_label(2, evr) == "Vol PC2 (12.0%)"
        )
        assert covariate_labels.landscape_vol_pc_pretty_label(3, evr) == "Vol PC3"

    def test_covariate_display_map_vol_pc_variance(self) -> None:
        evr = np.array([0.1], dtype=np.float64)
        m = covariate_labels.covariate_display_map(
            ["landscape_vol_PC1", "z0"],
            vol_pc_explained_variance_ratio=evr,
        )
        assert m["landscape_vol_PC1"] == "Vol PC1 (10.0%)"
        assert m["z0"] == "z0"

    @pytest.mark.parametrize(
        "raw,expected",
        [
            (None, "Viridis"),
            ("", "Viridis"),
            ("viridis", "Viridis"),
            ("PLASMA", "Plasma"),
            ("not_a_real_palette", "Viridis"),
        ],
    )
    def test_normalize_continuous_palette(self, raw: str | None, expected: str) -> None:
        assert palette_config.normalize_continuous_palette(raw) == expected

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

    def test_scatter3d_volume_landscape_style_scales_marker_size(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        """Regression: 3D volume-landscape UI keeps glyphs smaller than the default n-curve (relative factor unchanged after dashboard glyph scale)."""
        import json

        from cryodrgn.dashboard.plots_scatter import scatter3d_z_json

        base = json.loads(
            scatter3d_z_json(
                dashboard_experiment,
                "z0",
                "z1",
                "z2",
                None,
                volume_landscape_3d_style=False,
            )
        )
        vol = json.loads(
            scatter3d_z_json(
                dashboard_experiment,
                "z0",
                "z1",
                "z2",
                None,
                volume_landscape_3d_style=True,
            )
        )
        sz_b = base["data"][0]["marker"]["size"]
        sz_v = vol["data"][0]["marker"]["size"]
        assert isinstance(sz_b, (int, float))
        assert isinstance(sz_v, (int, float))
        factor = (1.0 - 0.31) * (1.0 - 0.13) ** 2 * 0.8 * 1.13 * 1.728 * 1.11
        assert float(sz_v) == pytest.approx(float(sz_b) * factor, rel=1e-6, abs=1e-9)
        assert float(sz_v) < float(sz_b)

    def test_scatter3d_continuous_covariate_uses_per_point_hex_marker_colors(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        """Continuous scatter3d markers must be per-point hex so letter overlays match WebGL."""
        import json

        from cryodrgn.dashboard.plots_scatter import scatter3d_z_json

        e = dashboard_experiment
        allow = frozenset({"z0", "z1", "z2"})
        fig = json.loads(
            scatter3d_z_json(
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
        import json

        import numpy as np

        from cryodrgn.dashboard.column_names import VOL_LANDSCAPE_NEAREST_SKETCH_VOL
        from cryodrgn.dashboard.plots_scatter import scatter3d_z_json

        e = dashboard_experiment
        sub = e.plot_df.iloc[:40].copy()
        sub[VOL_LANDSCAPE_NEAREST_SKETCH_VOL] = (
            np.arange(len(sub), dtype=np.int64) % 3 + 1
        )
        allow = frozenset({"z0", "z1", "z2"})
        js = scatter3d_z_json(
            e,
            "z0",
            "z1",
            "z2",
            "none",
            plot_df=sub,
            xyz_axes_allowed=allow,
            volume_landscape_3d_style=True,
        )
        fig = json.loads(js)
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
        import json

        import numpy as np

        from cryodrgn.dashboard.column_names import (
            VOL_LANDSCAPE_3D_PLOT_DF_ROW,
            VOL_LANDSCAPE_IS_SKETCH_CENTROID,
            VOL_LANDSCAPE_NEAREST_SKETCH_VOL,
        )
        from cryodrgn.dashboard.plots_scatter import scatter3d_z_json

        e = dashboard_experiment
        sub = e.plot_df.iloc[:40].copy()
        sub[VOL_LANDSCAPE_3D_PLOT_DF_ROW] = np.arange(40, dtype=np.int64)
        sub[VOL_LANDSCAPE_NEAREST_SKETCH_VOL] = (
            np.arange(len(sub), dtype=np.int64) % 3 + 1
        )
        sub[VOL_LANDSCAPE_IS_SKETCH_CENTROID] = (
            (np.arange(len(sub), dtype=np.int64) % 5 == 0)
        ).astype(np.int64)
        allow = frozenset({"z0", "z1", "z2"})
        js = scatter3d_z_json(
            e,
            "z0",
            "z1",
            "z2",
            "none",
            plot_df=sub,
            xyz_axes_allowed=allow,
            volume_landscape_3d_style=True,
        )
        fig = json.loads(js)
        meta = (fig.get("layout") or {}).get("meta") or {}
        assert meta.get("cdrgn_landscape_sketch_centroid_cd") is True
        row0 = fig["data"][0]["customdata"][0]
        assert int(row0[-1]) in (1, 2, 3)
        assert int(row0[-2]) in (0, 1)

    def test_subsample_preserving_sketch_centroids_keeps_centroid_rows(self) -> None:
        import numpy as np
        import pandas as pd

        from cryodrgn.dashboard.column_names import VOL_LANDSCAPE_IS_SKETCH_CENTROID
        from cryodrgn.dashboard.plots_figure_utils import (
            _subsample_preserving_sketch_centroids,
        )

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
        import json

        from cryodrgn.dashboard.plots_scatter import scatter3d_z_json

        e = dashboard_experiment
        sub = e.plot_df.iloc[:40].copy()
        allow = frozenset({"z0", "z1", "z2"})
        js = scatter3d_z_json(
            e,
            "z0",
            "z1",
            "z2",
            "none",
            plot_df=sub,
            xyz_axes_allowed=allow,
            volume_landscape_3d_style=True,
        )
        fig = json.loads(js)
        tr = fig["data"][0]
        assert tr["mode"] == "markers"
        assert "text" not in tr
        meta = (fig.get("layout") or {}).get("meta") or {}
        assert meta.get("cdrgn_landscape_vol_animation") is not True

    def test_latent_3d_vol_anim_wires_preview_loading_callback(self) -> None:
        """ChimeraX runs use ``setPreviewAnimLoading`` on the preview stack; main ``setRendering`` is scatter reload only."""
        root = Path(__file__).resolve().parents[1]
        p = root / "cryodrgn" / "dashboard" / "templates" / "latent_3d.html"
        text = p.read_text(encoding="utf-8")
        assert "setPreviewAnimLoading:" in text
        assert "latent3dVolAnimApi.afterPlotRedraw" in text
        assert "l3dva-preview-anim-overlay" in text
        assert "function setRendering(on" in text
        assert "setRendering(true" in text
        assert "setLatent3dPreviewAnimLoading" in text
        assert "CryoLatent3dVolLandscapeAnim.boot" in text
        vol_js = (
            root
            / "cryodrgn"
            / "dashboard"
            / "static"
            / "js"
            / "latent3d_landscape_vol_animations.js"
        )
        vol_text = vol_js.read_text(encoding="utf-8")
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
        root = Path(__file__).resolve().parents[1]
        js = root / "cryodrgn" / "dashboard" / "static" / "js" / "landscape_volpca.js"
        assert "referenceScatter3dBaseMarkerSize" in js.read_text(encoding="utf-8")

    def test_latent_3d_loads_scatter3d_scene_module_before_vol_anim(self) -> None:
        root = Path(__file__).resolve().parents[1]
        text = (
            root / "cryodrgn" / "dashboard" / "templates" / "latent_3d.html"
        ).read_text(encoding="utf-8")
        assert "plotly_scatter3d_scene.js" in text
        assert "lastLatent3dPersistedSceneSnap" in text
        assert "latent3dSceneSnapHasCamera" in text
        assert "inset: 0" in text
        assert "background: transparent" in text
        p3 = text.index("plotly_scatter3d_scene.js")
        assert p3 < text.index("(function() {")
        js = (
            root
            / "cryodrgn"
            / "dashboard"
            / "static"
            / "js"
            / "plotly_scatter3d_scene.js"
        )
        assert js.is_file()
        js_text = js.read_text(encoding="utf-8")
        assert "CryoPlotlyScatter3dScene" in js_text
        assert "applySnapshotToFigureLayout" in js_text
        assert "sceneRelayoutPatchCameraOnly" in js_text
        # Colour-covariate changes await legend ``refresh`` before vol-anim tail + scene restore.
        assert "pRefresh" in text and "cRefresh.then" in text


class TestPlotlyScatter3dSceneCameraPreserve:
    """Static contracts for ``plotly_scatter3d_scene.js`` orbit preservation."""

    @staticmethod
    def _scene_js() -> str:
        root = Path(__file__).resolve().parents[1]
        return (
            root
            / "cryodrgn"
            / "dashboard"
            / "static"
            / "js"
            / "plotly_scatter3d_scene.js"
        ).read_text(encoding="utf-8")

    def test_scene_relayout_patch_camera_only_is_camera_only(self) -> None:
        text = self._scene_js()
        fn = text.split("function sceneRelayoutPatchCameraOnly", 1)[1]
        body = fn.split("function restore(", 1)[0]
        assert '"scene.camera"' in body
        assert "scene.xaxis.range" not in body
        assert "scene.yaxis.range" not in body

    def test_marker_restyle_preserving_camera_uses_plotly_update(self) -> None:
        text = self._scene_js()
        fn = text.split("function traceMarkerRestylePreservingCamera", 1)[1]
        body = fn.split("function traceCoordsRestyleFromFigure", 1)[0]
        assert "Plotly.update(gd, upd, layoutPatch, [0])" in body
        assert "sceneRelayoutPatchCameraOnly(pinSnap)" in body

    def test_resolve_camera_prefers_webgl_get_camera(self) -> None:
        text = self._scene_js()
        fn = text.split("function resolveCameraFromGd", 1)[1]
        body = fn.split("function snapshot(", 1)[0]
        assert "scene._scene.getCamera" in body

    def test_apply_snapshot_camera_only_skips_axis_ranges(self) -> None:
        text = self._scene_js()
        fn = text.split("function applySnapshotCameraOnlyToFigureLayout", 1)[1]
        body = fn.split("function layoutPatchWithoutScene", 1)[0]
        assert "fig.layout.scene.camera" in body
        assert "xaxis" not in body or "applyRange" not in body


class TestLatent3dScatter3dCameraSnapBack:
    """Regression guards for latent-3D viewing-angle snap-back fixes."""

    @staticmethod
    def _latent_html() -> str:
        root = Path(__file__).resolve().parents[1]
        return (
            root / "cryodrgn" / "dashboard" / "templates" / "latent_3d.html"
        ).read_text(encoding="utf-8")

    def test_set_rendering_does_not_restore_camera_before_overlay(self) -> None:
        text = self._latent_html()
        fn = text.split("function setRendering(on, opts)", 1)[1]
        body = fn.split("function selectedLatentPalette", 1)[0]
        assert "Do not relayout the camera before showing the badge" in body
        assert "restoreCameraOnly" not in body

    def test_same_axes_reload_uses_nonblocking_rendering_overlay(self) -> None:
        text = self._latent_html()
        fn = text.split("function loadPlot()", 1)[1]
        body = fn.split("function buildLatent3dPostPayload", 1)[0]
        assert "useNonblockingOverlay" in body
        assert "latent3dSameAxesAsLast" in body
        assert "setRendering(true, { nonblocking: useNonblockingOverlay })" in body

    def test_load_plot_captures_live_scene_pin_before_fetch(self) -> None:
        text = self._latent_html()
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
        text = self._latent_html()
        assert "knocks the scatter3d camera back" in text
        tail = text.split("function runLatent3dPostResizeTail", 1)[1]
        assert "cRefresh.then(function()" in tail
        idx_refresh = tail.index("cRefresh.then")
        idx_overlay_off = tail.index("setRendering(false)")
        assert idx_refresh < idx_overlay_off

    def test_same_axes_marker_restyle_preserves_camera(self) -> None:
        text = self._latent_html()
        fn = text.split("function latent3dApplyFigurePreservingScene", 1)[1]
        body = fn.split("function latent3dStabilizeSceneAfterOverlay", 1)[0]
        assert "traceMarkerRestylePreservingCamera(fig, gd, pin)" in body
        assert "applySnapshotCameraOnlyToFigureLayout(fig, pin)" in body

    def test_rendering_overlay_is_not_stacked_over_webgl_plot(self) -> None:
        text = self._latent_html()
        assert "latent3d-rendering-overlay" in text
        assert "cryo-plot-rendering-overlay--nonblocking" in text
        # Full-bleed veil must not use backdrop-filter over the WebGL canvas.
        head = text.split("{% block content %}", 1)[0]
        assert "backdrop-filter" not in head or "No backdrop-filter" in head
        assert "#latent3d-rendering-overlay" in head
        assert "inset: 0" in head

    def test_afterplot_does_not_restore_stale_load_pin_after_user_orbit(self) -> None:
        text = self._latent_html()
        fn = text.split("function latent3dOnPlotlyAfterplot", 1)[1]
        body = fn.split("function latent3dPersistSceneAfterUserOrbit", 1)[0]
        assert "User orbited during scene hold" in body
        assert "cameraEyeDiffers(loadPin, liveAfter)" in body
        persist = text.split("function latent3dPersistSceneAfterUserOrbit", 1)[1]
        persist = persist.split("function latent3dOnPlotlyAfterplot", 1)[0]
        assert "gd._cryoLatent3dLoadScenePin = holdPin" in persist

    def test_volanim_discrete_legend_does_not_expand_from_plot_width(self) -> None:
        """Discrete mode must not set ``--cryo-discrete-legend-w`` from plot stack (squeezes the 3D view)."""
        text = self._latent_html()
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
