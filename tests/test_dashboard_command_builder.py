"""Command builder interface and dashboard CLI tests."""

from __future__ import annotations

import argparse
import logging
import os
import pytest
from pathlib import Path

from cryodrgn.commands import dashboard as dash_cli, train_vae
from cryodrgn.dashboard.command_builder_cli_help import (
    defaults_map_from_command_py,
    display_name_from_arg_help,
    help_map_from_command_py,
    jinja_arg_display_name,
    load_cli_help_maps,
    load_command_module_docstrings,
    resolve_help_default_placeholders,
    resolved_help_for_flag,
)
from cryodrgn.dashboard.command_builder_data import (
    COMMAND_BUILDER_BATCH_SIZE_ARG_IDS,
    COMMAND_BUILDER_COMMAND_KEYS,
    COMMAND_BUILDER_MANUSCRIPT_LABELS,
    COMMAND_BUILDER_MANUSCRIPT_URLS,
    COMMAND_BUILDER_REQUIRED_FIELD_TITLES,
    COMMAND_BUILDER_SCHEMA,
    arg_is_batch_size_denominated,
    arg_is_epoch_denominated,
    arg_is_num_epochs,
    arg_show_display_name,
    batch_size_arg_ids_for_command,
    default_outdir_for_command,
)
from cryodrgn.dashboard.context import command_builder_template_kwargs
from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.command_builder_page import _github_repo_release_url


class TestCommandModuleDocstrings:
    def test_loads_summaries_for_all_commands(self) -> None:
        docs = load_command_module_docstrings()
        assert set(docs.keys()) == {
            "abinit",
            "abinit_het_old",
            "abinit_homo_old",
            "train_vae",
            "train_nn",
            "train_dec",
            "backproject_voxel",
            "analyze",
            "analyze_landscape",
            "analyze_landscape_full",
        }
        assert "cryoDRGN-AI" in docs["abinit"]
        assert "VAE" in docs["train_vae"]
        assert "neural net" in docs["train_nn"]
        assert docs["train_dec"] == "Train an autodecoder"
        assert "backprojection" in docs["backproject_voxel"].lower()
        assert "latent space" in docs["analyze"].lower()
        assert "comparing volumes" in docs["analyze_landscape"].lower()
        assert "latent space" in docs["analyze_landscape_full"].lower()


class TestArgIsNumEpochs:
    def test_num_epochs_flags(self) -> None:
        assert arg_is_num_epochs({"cli": ["-n", "--num-epochs"], "w": "number"})

    def test_other_epoch_flags(self) -> None:
        assert not arg_is_num_epochs({"cli": ["--epochs-sgd"], "w": "number"})


class TestArgIsEpochDenominated:
    def test_epoch_flags(self) -> None:
        assert not arg_is_epoch_denominated(
            {"cli": ["-n", "--num-epochs"], "w": "number"}
        )
        assert arg_is_epoch_denominated({"cli": ["--epochs-sgd"], "w": "number"})
        assert arg_is_epoch_denominated({"cli": ["--l-ramp-epochs"], "w": "number"})
        assert arg_is_epoch_denominated({"cli": ["--pose-only-phase"], "w": "number"})

    def test_non_epoch_flags(self) -> None:
        assert not arg_is_epoch_denominated({"cli": ["--pretrain"], "w": "number"})
        assert not arg_is_epoch_denominated(
            {"cli": ["--n-imgs-pretrain"], "w": "number"}
        )
        assert not arg_is_epoch_denominated({"cli": ["--batch-size"], "w": "number"})


class TestArgIsBatchSizeDenominated:
    def test_batch_size_flags(self) -> None:
        assert arg_is_batch_size_denominated(
            {"cli": ["-b", "--batch-size"], "w": "number"}
        )
        assert arg_is_batch_size_denominated(
            {"cli": ["--batch-size-hps"], "w": "number"}
        )
        assert arg_is_batch_size_denominated(
            {"cli": ["--test-batch-size"], "w": "number"}
        )

    def test_non_batch_size_flags(self) -> None:
        assert not arg_is_batch_size_denominated(
            {"cli": ["-n", "--num-epochs"], "w": "number"}
        )
        assert not arg_is_batch_size_denominated({"cli": ["--pretrain"], "w": "number"})


class TestBatchSizeArgIdsForCommand:
    def test_abinit_lists_all_batch_size_fields(self) -> None:
        ids = batch_size_arg_ids_for_command("abinit")
        assert ids == [
            "ab_batch_size_hps",
            "ab_batch_size_known_poses",
            "ab_batch_size_sgd",
        ]

    def test_analyze_commands_have_no_batch_size_fields(self) -> None:
        assert batch_size_arg_ids_for_command("analyze") == []
        assert batch_size_arg_ids_for_command("analyze_landscape") == []

    def test_landscape_full_includes_train_and_test_batch_sizes(self) -> None:
        assert batch_size_arg_ids_for_command("analyze_landscape_full") == [
            "alfull_batch_size",
            "alfull_test_batch_size",
        ]

    def test_precomputed_map_covers_all_commands(self) -> None:
        assert set(COMMAND_BUILDER_BATCH_SIZE_ARG_IDS.keys()) == set(
            COMMAND_BUILDER_COMMAND_KEYS
        )


class TestArgShowDisplayName:
    def test_hides_when_unit_suffix_present(self) -> None:
        assert not arg_show_display_name(
            {"cli": ["-n", "--num-epochs"], "w": "number", "help": "x"}
        )
        assert not arg_show_display_name(
            {"cli": ["--epochs-sgd"], "w": "number", "help": "x"}
        )
        assert not arg_show_display_name(
            {"cli": ["--pose-only-phase"], "w": "number", "help": "x"}
        )
        assert not arg_show_display_name(
            {"cli": ["--batch-size-hps"], "w": "number", "help": "x"}
        )

    def test_keeps_short_batch_flag(self) -> None:
        assert arg_show_display_name(
            {"cli": ["-b", "--batch-size"], "w": "number", "help": "x"}
        )

    def test_shows_without_unit_suffix(self) -> None:
        assert arg_show_display_name({"cli": ["--wd"], "w": "text", "help": "x"})


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
        assert kw["default_workdir"] == ""
        assert kw["default_epoch"] == ""

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
        assert kw["default_workdir"] == dashboard_experiment.workdir
        assert kw["default_epoch"] == str(dashboard_experiment.epoch)


class TestCommandBuilderSchemaIntegrity:
    """Guard against drift between dashboard schema and real cryoDRGN CLIs."""

    def test_schema_covers_all_commands(self) -> None:
        assert set(COMMAND_BUILDER_SCHEMA.keys()) == set(COMMAND_BUILDER_COMMAND_KEYS)

    def test_groups_have_descriptions(self) -> None:
        for cmd, groups in COMMAND_BUILDER_SCHEMA.items():
            for g in groups:
                assert g.get("description"), f"{cmd}:{g['title']} missing description"

    @pytest.mark.parametrize("cmd", list(COMMAND_BUILDER_COMMAND_KEYS))
    def test_arg_ids_are_unique(self, cmd: str) -> None:
        ids = [a["id"] for g in COMMAND_BUILDER_SCHEMA[cmd] for a in g["args"]]
        assert len(ids) == len(set(ids)), f"duplicate ids in {cmd!r}"

    @pytest.mark.parametrize("cmd", list(COMMAND_BUILDER_COMMAND_KEYS))
    def test_every_cli_flag_has_help_entry(self, cmd: str) -> None:
        """Every ``--flag`` listed in the schema must exist in the CLI help map."""
        help_map = load_cli_help_maps().get(cmd, {})
        missing: list[str] = []
        for group in COMMAND_BUILDER_SCHEMA[cmd]:
            for arg in group["args"]:
                if arg.get("w") == "no_amp":
                    continue
                flags = arg.get("cli") or []
                # At least one of the declared CLI tokens should exist in help,
                # or the schema carries an explicit ``help`` string (e.g. --ctf-alg).
                if (
                    flags
                    and not any(c in help_map for c in flags)
                    and not arg.get("help")
                ):
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
            "bpv_particles",
            "bpv_out",
            "bpv_poses",
            "dec_particles",
            "dec_out",
            "dec_poses",
            "dec_zdim",
            "ana_workdir",
            "ana_epoch",
            "alsc_workdir",
            "alfull_workdir",
        }
        assert expected <= set(COMMAND_BUILDER_REQUIRED_FIELD_TITLES.keys())

    def test_all_values_nonempty(self) -> None:
        for k, v in COMMAND_BUILDER_REQUIRED_FIELD_TITLES.items():
            assert isinstance(v, str) and v.strip(), f"empty help for {k}"


class TestHelpMapFromCommandPy:
    def test_train_vae_has_common_flags(self) -> None:
        p = Path(train_vae.__file__)
        m = help_map_from_command_py(p)
        assert "-o" in m or "--outdir" in m
        assert "--zdim" in m
        assert "--poses" in m


class TestResolveHelpDefaults:
    def test_abinit_max_threads_default_in_tooltip(self) -> None:
        from pathlib import Path

        import cryodrgn.commands.abinit as abinit_mod

        dm = defaults_map_from_command_py(Path(abinit_mod.__file__))
        hm = help_map_from_command_py(Path(abinit_mod.__file__))
        assert dm["--max-threads"] == 16
        t = resolved_help_for_flag(hm, dm, "--max-threads")
        assert t is not None
        assert "%(default)s" not in t
        assert "(default: 16)" in t

    def test_resolve_placeholder_substitution(self) -> None:
        assert (
            resolve_help_default_placeholders("Batch (default: %(default)s)", 64)
            == "Batch (default: 64)"
        )


class TestDisplayNameFromArgHelp:
    def test_common_training_flags(self) -> None:
        assert (
            display_name_from_arg_help(
                "Weight decay for the optimizer (default: %(default)s)"
            )
            == "weight decay"
        )
        assert (
            display_name_from_arg_help(
                "Number of total epochs to train for (default: %(default)s)"
            )
            == "total epochs"
        )
        assert display_name_from_arg_help("Flag for lazy data loading.") == (
            "lazy data loading"
        )
        assert (
            display_name_from_arg_help(
                "Learning rate for the pose table optimizer (default: %(default)s)"
            )
            == "learning rate (pose table)"
        )

    def test_jinja_filter_accepts_arg_dict(self) -> None:
        wd_arg = next(
            a
            for g in COMMAND_BUILDER_SCHEMA["abinit"]
            for a in g["args"]
            if a.get("cli") == ["--wd"]
        )
        assert jinja_arg_display_name(wd_arg) == "weight decay"

    def test_ind_override(self) -> None:
        ind_arg = next(
            a
            for g in COMMAND_BUILDER_SCHEMA["abinit"]
            for a in g["args"]
            if a.get("cli") == ["--ind"]
        )
        assert jinja_arg_display_name(ind_arg) == "filtering .pkl"

    def test_max_threads_suppressed(self) -> None:
        mt_arg = next(
            a
            for g in COMMAND_BUILDER_SCHEMA["abinit"]
            for a in g["args"]
            if a.get("cli") == ["--max-threads"]
        )
        assert not arg_show_display_name(mt_arg)
        assert jinja_arg_display_name(mt_arg) == "threads"


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


class TestCommandBuilderManuscriptLabels:
    def test_reconstruction_commands_mapped(self) -> None:
        assert COMMAND_BUILDER_MANUSCRIPT_LABELS["train_nn"] == "cryoDRGN1\nmanuscript"
        assert COMMAND_BUILDER_MANUSCRIPT_LABELS["train_vae"] == "cryoDRGN1\nmanuscript"
        assert (
            COMMAND_BUILDER_MANUSCRIPT_LABELS["abinit_het_old"]
            == "cryoDRGN2\nmanuscript"
        )
        assert COMMAND_BUILDER_MANUSCRIPT_LABELS["abinit"] == "cryoDRGN-AI\nmanuscript"
        assert (
            COMMAND_BUILDER_MANUSCRIPT_LABELS["train_dec"] == "cryoDRGN-AI\nmanuscript"
        )

    def test_analyze_commands_not_mapped(self) -> None:
        for key in ("analyze", "analyze_landscape", "backproject_voxel"):
            assert key not in COMMAND_BUILDER_MANUSCRIPT_LABELS


class TestCommandBuilderManuscriptUrls:
    def test_reconstruction_commands_mapped(self) -> None:
        assert (
            COMMAND_BUILDER_MANUSCRIPT_URLS["abinit"]
            == "https://www.nature.com/articles/s41592-025-02720-4"
        )
        assert COMMAND_BUILDER_MANUSCRIPT_URLS["train_dec"].endswith(
            "s41592-025-02720-4"
        )
        assert COMMAND_BUILDER_MANUSCRIPT_URLS["train_vae"].endswith(
            "s41592-020-01049-4"
        )
        assert (
            "ICCV_2021_paper.pdf" in COMMAND_BUILDER_MANUSCRIPT_URLS["abinit_het_old"]
        )

    def test_analyze_commands_not_mapped(self) -> None:
        for key in ("analyze", "analyze_landscape", "backproject_voxel"):
            assert key not in COMMAND_BUILDER_MANUSCRIPT_URLS


class TestGithubRepoReleaseUrl:
    def test_pep440_prerelease_gets_hyphen(self) -> None:
        url = _github_repo_release_url(
            "https://github.com/ml-struct-bio/cryodrgn",
            "4.3.0a8",
        )
        assert url.endswith("/tree/4.3.0-a8")

    def test_tag_with_hyphen_unchanged(self) -> None:
        url = _github_repo_release_url(
            "https://github.com/ml-struct-bio/cryodrgn",
            "4.3.0-a8",
        )
        assert url.endswith("/tree/4.3.0-a8")

    def test_rc_prerelease(self) -> None:
        url = _github_repo_release_url(
            "https://github.com/ml-struct-bio/cryodrgn",
            "4.3.0rc1",
        )
        assert url.endswith("/tree/4.3.0-rc1")


class TestCommandBuilderStaticSite:
    """GitHub Pages bundle omits Plotly (command builder does not use charts)."""

    def test_rendered_html_has_no_plotly_script(self) -> None:
        from cryodrgn.dashboard.command_builder_page import render_command_builder_html

        html = render_command_builder_html()
        assert "cdn.plot.ly" not in html
        assert "/vendor/plotly.min.js" not in html
        assert "plotly.min.js" not in html


class TestCommandBuilderBrowserSmoke:
    """Headless Chromium: command line preview updates when command type changes."""

    def test_cmd_type_switch_updates_preview(
        self, playwright_page, dashboard_live_url
    ) -> None:
        from tests.conftest import dashboard_smoke_command_builder

        out = dashboard_smoke_command_builder(playwright_page, dashboard_live_url)
        assert out["initial_has_abinit"]
        assert out["switched_to_train_vae"]
