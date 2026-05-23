"""Command builder interface and dashboard CLI tests."""

from __future__ import annotations

import argparse
import logging
import os

import pytest

from cryodrgn.commands import dashboard as dash_cli, train_vae
from cryodrgn.dashboard.command_builder_cli_help import (
    help_map_from_command_py,
    load_cli_help_maps,
)
from cryodrgn.dashboard.command_builder_data import (
    COMMAND_BUILDER_REQUIRED_FIELD_TITLES,
    COMMAND_BUILDER_SCHEMA,
    default_outdir_for_command,
)
from cryodrgn.dashboard.context import command_builder_template_kwargs
from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.github_pages import (
    build_github_pages_site,
    render_command_builder_html,
)


class TestDefaultOutdirForCommand:
    def test_basename(self) -> None:
        assert default_outdir_for_command("train_vae") == "001_train_vae"

    def test_with_workdir(self, tmp_path) -> None:
        wd = str(tmp_path / "run")
        assert default_outdir_for_command("abinit", wd) == os.path.join(
            wd, "001_abinit"
        )


class TestCommandBuilderTemplateKwargs:
    def test_no_experiment_uses_defaults(self) -> None:
        kw = command_builder_template_kwargs(None)
        assert kw["default_zdim"] == 8
        assert kw["default_outdir_abinit"] == "001_abinit"
        assert kw["default_outdir_train_vae"] == "001_train_vae"
        assert kw["default_poses"] == ""
        assert "command_builder_schema" in kw
        assert "command_builder_required_field_titles" in kw

    def test_with_experiment_uses_config(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        kw = command_builder_template_kwargs(dashboard_experiment)
        assert kw["default_zdim"] == int(
            dashboard_experiment.train_configs["model_args"]["zdim"]
        )
        assert kw["default_outdir_abinit"].endswith("001_abinit")
        assert kw["default_outdir_train_vae"].endswith("001_train_vae")
        assert kw["default_outdir_train_nn"].endswith("001_train_nn")
        assert kw["default_outdir_train_dec"].endswith("001_train_dec")


class TestCommandBuilderSchemaIntegrity:
    """Guard against drift between dashboard schema and real cryoDRGN CLIs."""

    def test_schema_covers_all_four_commands(self) -> None:
        assert set(COMMAND_BUILDER_SCHEMA.keys()) == {
            "abinit",
            "train_vae",
            "train_nn",
            "train_dec",
        }

    def test_groups_have_descriptions(self) -> None:
        for cmd, groups in COMMAND_BUILDER_SCHEMA.items():
            for g in groups:
                assert g.get("description"), f"{cmd}:{g['title']} missing description"

    @pytest.mark.parametrize("cmd", ["abinit", "train_vae", "train_nn", "train_dec"])
    def test_arg_ids_are_unique(self, cmd: str) -> None:
        ids = [a["id"] for g in COMMAND_BUILDER_SCHEMA[cmd] for a in g["args"]]
        assert len(ids) == len(set(ids)), f"duplicate ids in {cmd!r}"

    @pytest.mark.parametrize("cmd", ["abinit", "train_vae", "train_nn", "train_dec"])
    def test_every_cli_flag_has_help_entry(self, cmd: str) -> None:
        """Every ``--flag`` listed in the schema must exist in the CLI help map."""
        help_map = load_cli_help_maps().get(cmd, {})
        missing: list[str] = []
        for group in COMMAND_BUILDER_SCHEMA[cmd]:
            for arg in group["args"]:
                if arg.get("w") == "no_amp":
                    continue
                flags = arg.get("cli") or []
                # At least one of the declared CLI tokens should exist in help.
                if flags and not any(c in help_map for c in flags):
                    missing.append(f"{arg['id']}:{flags}")
        assert not missing, f"{cmd}: schema flags missing from argparse help: {missing}"


class TestRequiredFieldTitles:
    def test_contains_expected_keys(self) -> None:
        expected = {
            "ab_particles",
            "ab_out",
            "ab_zdim",
            "vae_particles",
            "vae_out",
            "vae_poses",
            "vae_zdim",
            "nn_particles",
            "nn_out",
            "nn_poses",
            "dec_particles",
            "dec_out",
            "dec_poses",
            "dec_zdim",
        }
        assert expected <= set(COMMAND_BUILDER_REQUIRED_FIELD_TITLES.keys())

    def test_all_values_nonempty(self) -> None:
        for k, v in COMMAND_BUILDER_REQUIRED_FIELD_TITLES.items():
            assert isinstance(v, str) and v.strip(), f"empty help for {k}"


class TestHelpMapFromCommandPy:
    def test_train_vae_has_common_flags(self) -> None:
        from pathlib import Path

        p = Path(train_vae.__file__)
        m = help_map_from_command_py(p)
        assert "-o" in m or "--outdir" in m
        assert "--zdim" in m
        assert "--poses" in m


class TestDashboardCLI:
    @staticmethod
    def _parse(argv: list[str]) -> argparse.Namespace:
        parser = argparse.ArgumentParser()
        dash_cli.add_args(parser)
        return parser.parse_args(argv)

    def test_parses_minimal_args(self) -> None:
        ns = self._parse([])
        assert ns.outdir is None
        assert ns.epoch == -1
        assert ns.kmeans == -1
        assert ns.port == 5050
        assert ns.cpus == 4
        assert ns.verbose == 0
        assert ns.command_builder is False
        assert ns.particle_selection is False

    def test_parses_with_outdir(self) -> None:
        ns = self._parse(["/tmp/some_run", "--epoch", "3", "--kmeans", "5"])
        assert ns.outdir.endswith("some_run")
        assert ns.epoch == 3
        assert ns.kmeans == 5

    def test_view_flag_aliases(self) -> None:
        ns = self._parse(["--filter"])
        assert ns.particle_selection is True

    def test_verbose_count_levels(self) -> None:
        assert self._parse(["-v"]).verbose == 1
        assert self._parse(["-vv"]).verbose == 2
        assert self._parse(["-vvv"]).verbose == 3

    def test_view_flags_are_mutually_exclusive(self) -> None:
        parser = argparse.ArgumentParser()
        dash_cli.add_args(parser)
        with pytest.raises(SystemExit):
            parser.parse_args(["--pair-grid", "--command-builder"])

    def test_filter_max_points_sets_env(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", raising=False)
        monkeypatch.setattr(dash_cli, "run_server", lambda **kw: None)
        ns = self._parse(
            ["--filter-max", "123000", "--no-browser", "--command-builder"]
        )
        dash_cli.main(ns)
        assert os.environ.get("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS") == "123000"

    def test_builder_only_with_experiment_view_raises(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(dash_cli, "run_server", lambda **kw: None)
        ns = self._parse(["--particle-selection", "--no-browser"])
        with pytest.raises(ValueError, match="need an output directory"):
            dash_cli.main(ns)

    def test_builder_only_with_command_builder_ok(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        called = {}

        def fake_run_server(**kwargs) -> None:
            called.update(kwargs)

        monkeypatch.setattr(dash_cli, "run_server", fake_run_server)
        ns = self._parse(["--command-builder", "--no-browser"])
        dash_cli.main(ns)
        assert called["workdir"] is None
        assert called["port"] == 5050

    def test_main_invokes_run_server_with_outdir(
        self,
        monkeypatch: pytest.MonkeyPatch,
        dashboard_workdir: str,
    ) -> None:
        called = {}

        def fake_run_server(**kwargs) -> None:
            called.update(kwargs)

        monkeypatch.setattr(dash_cli, "run_server", fake_run_server)
        ns = self._parse([dashboard_workdir, "--no-browser"])
        dash_cli.main(ns)
        assert called["workdir"] == dashboard_workdir
        assert called["epoch"] == -1
        assert called.get("discovery_root") is None

    def test_parent_outdir_sets_discovery_mode(
        self, tmp_path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        run = tmp_path / "run1"
        (run / "analyze.0").mkdir(parents=True)
        (run / "z.0.pkl").write_bytes(b"")
        called: dict = {}

        def fake_run_server(**kwargs) -> None:
            called.update(kwargs)

        monkeypatch.setattr(dash_cli, "run_server", fake_run_server)
        ns = self._parse([str(tmp_path), "--no-browser"])
        dash_cli.main(ns)
        assert called["workdir"] is None
        assert called["discovery_root"] == os.path.abspath(str(tmp_path))

    def test_single_child_outdir_passes_workdir_not_discovery(
        self, tmp_path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        run = tmp_path / "run1"
        (run / "analyze.0").mkdir(parents=True)
        (run / "z.0.pkl").write_bytes(b"")
        called: dict = {}

        def fake_run_server(**kwargs) -> None:
            called.update(kwargs)

        monkeypatch.setattr(dash_cli, "run_server", fake_run_server)
        ns = self._parse([str(run), "--no-browser"])
        dash_cli.main(ns)
        assert called["workdir"] == os.path.abspath(str(run))
        assert called.get("discovery_root") is None

    def test_outdir_not_a_directory_raises(self, tmp_path) -> None:
        f = tmp_path / "notadir"
        f.write_text("x")
        ns = self._parse([str(f), "--no-browser"])
        with pytest.raises(ValueError, match="Not a directory"):
            dash_cli.main(ns)

    def test_main_configures_logging_from_verbose(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        called = {}

        def fake_cfg(verbosity: int) -> None:
            called["verbosity"] = verbosity

        monkeypatch.setattr(dash_cli, "_configure_dashboard_logging", fake_cfg)
        monkeypatch.setattr(dash_cli, "run_server", lambda **kw: None)
        ns = self._parse(["--command-builder", "--no-browser", "-vv"])
        dash_cli.main(ns)
        assert called["verbosity"] == 2

    def test_default_logging_suppresses_werkzeug_internal(self) -> None:
        dash_cli._configure_dashboard_logging(0)
        assert (
            logging.getLogger("werkzeug._internal").getEffectiveLevel()
            >= logging.WARNING
        )


class TestGithubPagesCommandBuilder:
    def test_render_includes_schema_and_builder_ui(self) -> None:
        html = render_command_builder_html(
            base_path="/cryodrgn/",
            repo_url="https://github.com/ml-struct-bio/cryodrgn",
        )
        assert '<base href="/cryodrgn/"/>' in html
        assert 'class="cmd-builder-github-pages"' in html
        assert "github-pages-note" in html
        assert "var CMD_SCHEMA" in html
        assert "cryodrgn abinit" in html
        assert 'class="cmd-group-card"' in html
        assert "cmd-group-card-trigger" in html
        assert "cmd-group-card-desc" in html
        assert "cmd-builder-required-bar" in html
        assert "cmd-builder-primary-pair" in html
        assert "cmd-builder-primary-region--dataset" in html
        assert "cmd-builder-primary-region--training" in html
        assert "initGithubPagesGroupCards" in html
        assert "nav-reconstruction-cmd" in html
        assert html.count('id="cmd-type"') == 1
        assert "No experiment loaded" not in html
        assert '<span class="nav-page-title">' not in html
        assert 'class="cmd-builder-program"' not in html
        assert "cmd-outdir-stepper" in html
        assert 'value="001_abinit"' in html
        assert "plot.ly" not in html.lower()
        assert 'fetch("/api/set_workdir"' not in html

    def test_build_writes_index_and_nojekyll(self, tmp_path) -> None:
        out = build_github_pages_site(
            tmp_path / "site",
            base_path="/cryodrgn/",
            repo_url="https://github.com/ml-struct-bio/cryodrgn",
        )
        assert (out / "index.html").is_file()
        assert (out / ".nojekyll").is_file()
        body = (out / "index.html").read_text(encoding="utf-8")
        assert "cmd-form" in body
