"""Command builder interface and dashboard CLI tests."""

from __future__ import annotations

import argparse
import logging
import os

import pytest

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
    COMMAND_BUILDER_COMMAND_KEYS,
    COMMAND_BUILDER_MANUSCRIPT_LABELS,
    COMMAND_BUILDER_MANUSCRIPT_URLS,
    COMMAND_BUILDER_REQUIRED_FIELD_TITLES,
    COMMAND_BUILDER_SCHEMA,
    arg_is_batch_size_denominated,
    arg_is_epoch_denominated,
    arg_is_num_epochs,
    arg_show_display_name,
    default_outdir_for_command,
)
from cryodrgn.dashboard.context import command_builder_template_kwargs
from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.command_builder_page import (
    _github_repo_release_url,
    build_command_builder_page_site,
    render_command_builder_html,
)


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
        from pathlib import Path

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


class TestCommandBuilderPage:
    def test_render_includes_schema_and_builder_ui(self) -> None:
        html = render_command_builder_html(
            base_path="/cryodrgn/",
            repo_url="https://github.com/ml-struct-bio/cryodrgn",
        )
        assert '<base href="/cryodrgn/"/>' in html
        assert 'class="cmd-builder-github-pages"' in html
        assert "github-pages-note" not in html
        assert "nav-github-release-link" in html
        assert 'id="nav-manuscript-link"' in html
        assert "nav-manuscript-link" in html
        assert "CMD_MANUSCRIPT_URLS" in html
        assert "CMD_MANUSCRIPT_LABELS" in html
        assert 'id="nav-manuscript-label"' in html
        assert "cryoDRGN-AI" in html
        assert "s41592-025-02720-4" in html
        assert "s41592-020-01049-4" in html
        assert "ICCV_2021_paper.pdf" in html
        assert "syncManuscriptLink" in html
        assert '<span class="nav-brand-home-arrow"' not in html
        assert (
            "color: #fff"
            in html.split("body.cmd-builder-github-pages .cmd-output pre")[1][:500]
        )
        assert (
            '"Cascadia Code"'
            in html.split("body.cmd-builder-github-pages .cmd-output pre")[1][:500]
        )
        assert "github.com/ml-struct-bio/cryodrgn/tree/" in html
        assert "var CMD_SCHEMA" in html
        assert "cryodrgn abinit" in html
        assert 'class="cmd-group-card"' in html
        assert "cmd-group-card-trigger" in html
        assert "cmd-group-card-desc" in html
        assert "cmd-builder-required-bar" in html
        assert "cmd-builder-primary-region-content" in html
        req_heading_css = html.split(".cmd-builder-required-heading")[1][:200]
        assert "calc(var(--gp-nav-cmd-label-font) * 1.25)" in req_heading_css
        assert "font-variant: normal" in req_heading_css
        assert "font-weight: 700" in req_heading_css
        assert "Dataset<br>loading" in html
        assert "dataset-split > .cmd-builder-primary-col--checks" in html
        checks_rule = html.split("dataset-split > .cmd-builder-primary-col--checks", 1)[
            1
        ]
        assert "align-items: center" in checks_rule[:250]
        assert "Training<br>parameters" in html
        assert "border-radius: 0" in html.split(".cmd-builder-required-cell")[1][:400]
        assert "cmd-builder-primary-pair" in html
        assert "cmd-builder-primary-region--dataset" in html
        assert "cmd-builder-dataset-core" in html
        assert "cmd-builder-dataset-rest-fields" in html
        ab_panel = html.split('id="panel-abinit"', 1)[1].split(
            'id="panel-train_vae"', 1
        )[0]
        core_pos = ab_panel.index('<div class="cmd-builder-dataset-core">')
        assert ab_panel.index('id="ab_ind"', core_pos) < ab_panel.index(
            'id="ab_datadir"', core_pos
        )
        rest_pos = ab_panel.index(
            '<div class="cmd-builder-dataset-rest-fields', core_pos
        )
        rest_chunk = ab_panel[rest_pos : rest_pos + 2000]
        assert 'id="ab_ind"' not in rest_chunk
        assert 'id="ab_datadir"' not in rest_chunk
        assert "display: contents" in html
        train_panel = html.split('id="panel-abinit"', 1)[1].split(
            "cmd-builder-primary-region--training", 1
        )[1]
        assert "cmd-arg-unit--num-epochs" in train_panel
        assert '<span class="cmd-arg-unit-line">total</span>' in train_panel
        assert 'cmd-arg-unit--scaled" aria-hidden="true">epochs</span>' in train_panel
        assert "cmd-builder-primary-region--training" in html
        assert "initCommandBuilderPageGroupCards" in html
        assert "nav-reconstruction-cmd" in html
        assert html.count('id="cmd-type"') == 1
        assert "No experiment loaded" not in html
        assert '<span class="nav-page-title">' not in html
        assert 'class="cmd-builder-program"' not in html
        assert "cmd-outdir-stepper" in html
        assert "Output folder" in html
        assert 'value="001_abinit"' in html
        assert 'value="001_abinit/"' not in html
        assert "cmd-outdir-input-suffix" in html
        assert "nav-cmd-doc" in html
        assert "nav-cmd-doc-text" in html
        assert "CMD_COMMAND_DOCS" in html
        assert "heterogeneous reconstruction" in html
        assert "abinit_het_old" in html
        assert "abinit_homo_old" in html
        assert 'value="001_abinit_het_old"' in html
        assert "backproject_voxel" in html
        assert 'value="001_backproject_voxel"' in html
        assert "analyze_landscape" in html
        assert 'id="ana_workdir"' in html
        assert "cryodrgn analyze" in html
        assert "optgroup" in html
        assert "plot.ly" not in html.lower()
        assert 'fetch("/api/set_workdir"' not in html
        assert 'id="cmd-colorize"' in html
        assert "Colorize" in html
        assert "cmd-colorize-rainbow" in html
        assert "renderColoredCommand" in html
        assert 'id="cmd-copy"' in html
        assert "Copy" in html
        assert "btn-copy-icon" in html
        assert "btn-copy-label-line" in html
        assert 'stroke="#000"' in html
        assert "max-height: 10dvh" in html
        assert "--gp-font-scale-root: calc(1.221 * 0.9)" in html
        assert "--gp-nav-cmd-label-font:" in html
        assert "--gp-arg-region-title-font:" in html
        assert "CMD_WRAP_MIN_CHARS = 80" in html
        assert "CMD_LINE_MAX_CHARS = 100" in html
        assert "layoutCommandLineGroups" in html
        assert "CMD_LINE_HANG_SPACES" in html
        assert "cmd-line-hang" in html
        assert "--cmd-line-hang-indent: 2ch" in html
        assert "updateCommandDisplay" in html
        assert "Advanced parameters" in html
        assert "cmd-builder-advanced-region" in html
        assert "gap: calc(0.45rem * 3 * 1.13)" in html
        assert "min-width: calc(4.75rem * 1.3)" in html
        assert "max-height: 17dvh" in html
        assert "cmd-builder-cmd-dock" in html
        assert "cmd-arg-unit--scaled" in html
        assert "cmd-arg-unit--num-epochs" in html
        assert 'id="ab_n"' in html
        ab_n_chunk = html[html.index('id="ab_n"') : html.index('id="ab_n"') + 500]
        assert "cmd-arg-unit--num-epochs" in ab_n_chunk
        assert '<span class="cmd-arg-unit-line">total</span>' in ab_n_chunk
        assert "cmd-arg-unit--scaled" in html.split('id="ab_epochs_sgd"', 1)[1][:400]
        pose_only_chunk = html.split('id="ab_pose_only_phase"', 1)[1][:450]
        assert (
            'cmd-arg-unit--scaled" aria-hidden="true">epochs</span>' in pose_only_chunk
        )
        ab_n_label = html.split('for="ab_n"', 1)[1][:280]
        assert "cmd-arg-display-name" not in ab_n_label
        assert "cmd-arg-display-name" in html
        assert "font-style: italic" in html.split(".cmd-arg-display-name")[1][:120]
        assert 'class="cmd-arg-display-name">weight decay</span>' in html
        assert (
            '<code>--wd</code><span class="cmd-arg-display-name">weight decay</span>'
            in html
        )
        vae_b_chunk = html.split('for="vae_b"', 1)[1][:500]
        assert 'cmd-arg-unit--scaled" aria-hidden="true">images</span>' in vae_b_chunk
        assert "cmd-arg-display-name" in vae_b_chunk
        assert "command-builder-manuscript.png" in html
        assert "nav-manuscript-icon" in html
        assert "nav-header-icons" in html
        assert "nav-header-icon-label" in html
        assert ">code</span>" in html
        assert (
            COMMAND_BUILDER_MANUSCRIPT_LABELS["abinit"]
            in html.split('id="nav-manuscript-label"', 1)[1][:120]
        )
        assert "command builder" in html
        gh_svg_css = html.split(
            "body.cmd-builder-github-pages .nav-github-release-link svg {"
        )[1][:120]
        assert "var(--gp-nav-github-icon-scale)" in gh_svg_css
        assert "--gp-nav-github-icon-scale: 0.576" in html
        assert "nav-header-icon-item--code" in html
        assert "nav-header-icon-item--paper" in html
        assert "flex-direction: row" in html.split(".nav-header-icon-item {")[1][:120]
        paper_item = html.split('id="nav-manuscript-icon-item"', 1)[1][:500]
        assert "nav-manuscript-link" in paper_item.split("nav-header-icon-label", 1)[0]
        assert 'href="https://cryodrgn.cs.princeton.edu/"' in html
        assert 'class="nav-brand-link"' in html
        gp_logo_hover = html.split(
            "body.cmd-builder-github-pages a.nav-brand-link:hover .nav-logo", 1
        )[1][:200]
        assert "filter: none" in gp_logo_hover
        assert "transform: none" in gp_logo_hover
        assert ".nav a.nav-brand-link:hover .nav-logo" in html
        assert "brightness(1.14) drop-shadow" in html
        assert (
            "body.cmd-builder-github-pages a.nav-brand-link:hover .nav-brand-line2-main"
            in html
        )
        assert (
            "background: transparent"
            in html.split(
                "body.cmd-builder-github-pages a.nav-brand-link:hover .nav-brand-line2-main"
            )[1][:400]
        )
        assert (
            "gap: 0.02rem"
            in html.split("body.cmd-builder-github-pages header.nav .nav-brand-text")[
                1
            ][:80]
        )
        assert "gap: 0.02rem" in html.split(".nav-header-icon-item--code")[1][:80]
        nav_icons_css = html.split(".nav-header-icons {", 1)[1][:320]
        assert "grid-column: 2 / 5" in nav_icons_css
        assert "margin-inline-start: 35%" in nav_icons_css
        assert "translateX(-50%)" in nav_icons_css
        ms_link_css = html.split(
            "body.cmd-builder-github-pages .nav-manuscript-link {", 1
        )[1][:360]
        assert "overflow: visible" in ms_link_css
        assert "gp-nav-manuscript-icon-zoom" not in html
        assert "gp-nav-manuscript-icon-width-scale" not in html
        assert "--gp-nav-bar-logo-size: min(3.88rem, 12vw)" in html
        assert "--gp-nav-header-icon-size: calc(2.142rem * 1.2)" in html
        assert "* 0.8 * var(--gp-nav-font-scale)" in html
        assert "var(--gp-nav-cmd-label-font) * 0.9" in html
        paper_item_css = html.split(
            "body.cmd-builder-github-pages .nav-header-icon-item--paper {", 1
        )[1][:80]
        assert "gap: 0.35rem" in paper_item_css
        gh_hover = html.split(
            "body.cmd-builder-github-pages .nav-github-release-link:hover,", 1
        )[1][:420]
        assert "background: transparent" in gh_hover
        icon_hover_fx = html.split(
            "body.cmd-builder-github-pages .nav-header-icon-item--code:hover .nav-github-release-link svg,",
            1,
        )[1][:700]
        assert "filter: brightness(1.14)" in icon_hover_fx
        assert "transform: scale(1.07)" in icon_hover_fx
        ms_img_css = html.rsplit(
            "body.cmd-builder-github-pages .nav-manuscript-link img.nav-manuscript-icon {",
            1,
        )[1][:200]
        assert "max-height: var(--gp-nav-bar-logo-size)" in ms_img_css
        cmd_dock_css = html.split(
            "body.cmd-builder-github-pages .cmd-builder-cmd-dock {"
        )[1][:520]
        assert "border-radius: 8px" in cmd_dock_css
        assert "rgba(255, 228, 232, 0.92)" in cmd_dock_css
        assert 'id="nav-manuscript-icon-item"' in html
        assert "grid-template-columns: auto 1fr auto 1fr auto" in html
        assert "filtering .pkl" in html
        assert "%(default)s" not in html
        assert "(default: 16)" in html.split('id="ab_max_threads"', 1)[1][:400]
        zdim_chunk = html.split('for="ab_zdim"', 1)[1][:500]
        assert "conformation latent space" in zdim_chunk

    def test_local_base_path_does_not_break_static_assets(self) -> None:
        html = render_command_builder_html(
            base_path="/",
            repo_url="https://github.com/ml-struct-bio/cryodrgn",
        )
        head = html.split("</head>", 1)[0]
        assert '<base href="/"/>' in head or '<base href="/">' in head
        assert "github.com" not in head.split("<base ", 1)[1].split(">", 1)[0]
        assert 'src="static/img/command-builder-manuscript.png"' in html
        assert 'href="https://cryodrgn.cs.princeton.edu/"' in html

    def test_build_writes_index_and_nojekyll(self, tmp_path) -> None:
        out = build_command_builder_page_site(
            tmp_path / "site",
            base_path="/cryodrgn/",
            repo_url="https://github.com/ml-struct-bio/cryodrgn",
        )
        assert (out / "index.html").is_file()
        assert (out / ".nojekyll").is_file()
        body = (out / "index.html").read_text(encoding="utf-8")
        assert "cmd-form" in body
