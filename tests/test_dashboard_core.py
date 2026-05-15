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
        assert VOL_LANDSCAPE_3D_PLOT_DF_ROW not in cols
