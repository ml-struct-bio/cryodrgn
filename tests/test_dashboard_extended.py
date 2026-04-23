"""Extended dashboard tests (context, CLI, and extra API paths).

Split from ``tests/test_dashboard.py`` to keep module size manageable.
"""

from __future__ import annotations

import argparse
import logging
import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import yaml

from cryodrgn.commands import dashboard as dash_cli, train_vae
from cryodrgn.dashboard import app as dash_app
from cryodrgn.dashboard.app import _TRAJECTORY_INELIGIBLE_MSG
from cryodrgn.dashboard.command_builder_cli_help import (
    help_map_from_command_py,
    load_cli_help_maps,
)
from cryodrgn.dashboard.command_builder_data import (
    COMMAND_BUILDER_REQUIRED_FIELD_TITLES,
    COMMAND_BUILDER_SCHEMA,
)
from cryodrgn.dashboard.context import (
    EXP_CACHE,
    PRELOAD_CACHE,
    _abbrev_middle_token,
    _argv_four_command_lines,
    _cmd_argv_for_nav_display,
    _config_has_cryodrgn_cmd,
    _workdir_options,
    abbrev_middle,
    active_workdir,
    clear_experiment_caches,
    command_builder_template_kwargs,
    discover_cryodrgn_workdirs,
    epochs_for_workdir,
    resolve_epoch,
)
from cryodrgn.dashboard.data import DashboardExperiment, list_z_epochs
from cryodrgn.dashboard.explorer_volumes import (
    _config_yaml_path,
    _is_drgnai_config,
    _mpl_retrim_png,
    _sorted_vol_mrc_paths,
    explorer_volumes_eligible,
)
from cryodrgn.dashboard.plots import (
    _continuous_series_stats,
    mpl_cmap_for_palette,
    normalize_continuous_palette,
    pair_grid_png,
    pair_grid_skeleton_placeholder_layout,
)
from cryodrgn.dashboard.trajectory import (
    _TRAJ_GRAPH_NEIGHBOR_CACHE,
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


class TestCommandBuilderTemplateKwargs:
    def test_no_experiment_uses_defaults(self) -> None:
        kw = command_builder_template_kwargs(None)
        assert kw["default_zdim"] == 8
        assert kw["default_outdir_abinit"] == "abinit_run"
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
        assert kw["default_outdir_abinit"].endswith("abinit_run")
        assert kw["default_outdir_train"].endswith("train_next")


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


# ---------------------------------------------------------------------------
# plots.py: palettes and continuous-series stats
# ---------------------------------------------------------------------------


class TestNormalizeContinuousPalette:
    @pytest.mark.parametrize(
        "raw,expected",
        [
            (None, "Viridis"),
            ("", "Viridis"),
            ("viridis", "Viridis"),
            ("TURBO", "Turbo"),
            ("  Plasma  ", "Plasma"),
            ("not_a_palette", "Viridis"),
            (123, "Viridis"),
        ],
    )
    def test_cases(self, raw: object, expected: str) -> None:
        pal_raw = raw if (raw is None or isinstance(raw, str)) else str(raw)
        assert normalize_continuous_palette(pal_raw) == expected

    def test_mpl_cmap_for_palette(self) -> None:
        assert mpl_cmap_for_palette("Viridis") == "viridis"
        # Unknown falls back to viridis.
        assert mpl_cmap_for_palette("not_a_palette") == "viridis"


class TestContinuousSeriesStats:
    def test_constant_series_has_finite_span(self) -> None:
        s = pd.Series([2.0, 2.0, 2.0])
        _, cmin, cmax = _continuous_series_stats(s)
        assert cmin < cmax

    def test_all_nan_falls_back(self) -> None:
        s = pd.Series([np.nan, np.nan])
        vals, cmin, cmax = _continuous_series_stats(s)
        assert cmin == 0.0 and cmax == 1.0
        assert np.isnan(vals).all()

    def test_mixed_series_uses_extrema(self) -> None:
        s = pd.Series([1.0, np.nan, 3.0, 5.0])
        _, cmin, cmax = _continuous_series_stats(s)
        assert cmin == 1.0 and cmax == 5.0


class TestPairGridHexAndSkeleton:
    def test_hex_style_is_png_and_deterministic(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        png_a, cells_a = pair_grid_png(
            dashboard_experiment,
            lower_color_col="labels",
            diagonal_emb="pc",
            upper_style="hex",
        )
        png_b, cells_b = pair_grid_png(
            dashboard_experiment,
            lower_color_col="labels",
            diagonal_emb="pc",
            upper_style="hex",
        )
        assert png_a[:8] == b"\x89PNG\r\n\x1a\n"
        assert png_a == png_b
        assert cells_a == cells_b

    def test_placeholder_layout_shape(self) -> None:
        cells = pair_grid_skeleton_placeholder_layout(3)
        assert len(cells) == 9
        for c in cells:
            assert {"x0", "y0", "x1", "y1"} <= set(c.keys())

    def test_placeholder_zdim_zero_is_empty(self) -> None:
        assert pair_grid_skeleton_placeholder_layout(0) == []


# ---------------------------------------------------------------------------
# explorer_volumes.py: helpers that don't require ChimeraX / CUDA
# ---------------------------------------------------------------------------


class TestIsDrgnaiConfig:
    def test_recognises_drgnai(self) -> None:
        assert _is_drgnai_config({"data_norm_mean": 0.0, "data_norm_std": 1.0})

    def test_classic_config_is_false(self) -> None:
        assert not _is_drgnai_config({"cmd": ["cryodrgn", "train_vae"]})


class TestConfigYamlPath:
    def test_returns_yaml_when_present(self, dashboard_workdir: str) -> None:
        assert _config_yaml_path(dashboard_workdir).endswith("config.yaml")

    def test_falls_back_to_pkl(self, tmp_path) -> None:
        (tmp_path / "config.pkl").write_bytes(b"\x00")
        assert _config_yaml_path(str(tmp_path)).endswith("config.pkl")

    def test_missing_raises(self, tmp_path) -> None:
        with pytest.raises(FileNotFoundError):
            _config_yaml_path(str(tmp_path))


class TestSortedVolMrcPaths:
    def test_sorts_by_index_and_caps_count(self, tmp_path) -> None:
        for i in [10, 2, 5, 7]:
            (tmp_path / f"vol_{i:03d}.mrc").write_bytes(b"\x00")
        out = _sorted_vol_mrc_paths(str(tmp_path), 3)
        assert [os.path.basename(p) for p in out] == [
            "vol_002.mrc",
            "vol_005.mrc",
            "vol_007.mrc",
        ]

    def test_insufficient_volumes_raises(self, tmp_path) -> None:
        (tmp_path / "vol_001.mrc").write_bytes(b"\x00")
        with pytest.raises(RuntimeError, match="Expected 5 volumes"):
            _sorted_vol_mrc_paths(str(tmp_path), 5)


class TestExplorerVolumesEligible:
    def test_false_without_weights(
        self, dashboard_experiment: DashboardExperiment, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(
            "cryodrgn.dashboard.explorer_volumes.torch_cuda_available", lambda: True
        )
        weight_path = os.path.join(
            dashboard_experiment.workdir, f"weights.{dashboard_experiment.epoch}.pkl"
        )
        orig_isfile = os.path.isfile

        def fake_isfile(path: str) -> bool:
            if os.path.abspath(path) == os.path.abspath(weight_path):
                return False
            return orig_isfile(path)

        monkeypatch.setattr("os.path.isfile", fake_isfile)
        assert not explorer_volumes_eligible(dashboard_experiment)

    def test_true_with_weights_and_cuda(
        self,
        dashboard_experiment: DashboardExperiment,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path,
    ) -> None:
        monkeypatch.setattr(
            "cryodrgn.dashboard.explorer_volumes.torch_cuda_available", lambda: True
        )
        weight = os.path.join(
            dashboard_experiment.workdir, f"weights.{dashboard_experiment.epoch}.pkl"
        )
        created = False
        if not os.path.exists(weight):
            with open(weight, "wb") as fh:
                fh.write(b"stub")
            created = True
        try:
            assert explorer_volumes_eligible(dashboard_experiment)
        finally:
            if created:
                os.remove(weight)

    def test_false_without_cuda(
        self, dashboard_experiment: DashboardExperiment, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(
            "cryodrgn.dashboard.explorer_volumes.torch_cuda_available", lambda: False
        )
        assert not explorer_volumes_eligible(dashboard_experiment)


class TestMplRetrimPng:
    def test_rewrites_png_in_place(self, tmp_path) -> None:
        fig, ax = plt.subplots(figsize=(3, 3))
        ax.plot([0, 1], [0, 1])
        p = tmp_path / "x.png"
        fig.savefig(p)
        plt.close(fig)
        before = p.read_bytes()
        assert before[:8] == b"\x89PNG\r\n\x1a\n"
        _mpl_retrim_png(str(p), dpi=80)
        after = p.read_bytes()
        assert after[:8] == b"\x89PNG\r\n\x1a\n"
        assert after != before


# ---------------------------------------------------------------------------
# data.py: list_z_epochs edge cases beyond the smoke check
# ---------------------------------------------------------------------------


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
# command_builder_data.py: schema / help drift guardrails
# ---------------------------------------------------------------------------


class TestCommandBuilderSchemaIntegrity:
    """Guard against drift between dashboard schema and real cryoDRGN CLIs."""

    def test_schema_covers_all_four_commands(self) -> None:
        assert set(COMMAND_BUILDER_SCHEMA.keys()) == {
            "abinit",
            "train_vae",
            "train_nn",
            "train_dec",
        }

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


# ---------------------------------------------------------------------------
# app.py: Flask error paths, session switching, and route table integrity
# ---------------------------------------------------------------------------


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


class TestSaveSelectionRoundTrip:
    def test_pkl_roundtrip(self, flask_client, tmp_path) -> None:
        r = flask_client.post(
            "/api/save_selection",
            json={
                "rows": [0, 1, 2, 3],
                "basename": "pytest_sel",
                "sel_dir": str(tmp_path),
            },
        )
        assert r.status_code == 200, r.get_json()
        js = r.get_json()
        assert js["ok"] is True
        assert os.path.isfile(js["path"])
        with open(js["path"], "rb") as fh:
            saved = pickle.load(fh)
        assert saved.size == 4


class TestSavePairPlotPng:
    def test_writes_png(self, flask_client, tmp_path) -> None:
        r = flask_client.post(
            "/api/save_pairplot_png",
            json={
                "color_col": "labels",
                "diagonal_emb": "pc",
                "upper_style": "scatter",
                "filename": "pytest_pairplot.png",
            },
        )
        assert r.status_code == 200
        js = r.get_json()
        assert js["ok"] is True
        assert os.path.isfile(js["path"])
        with open(js["path"], "rb") as fh:
            assert fh.read(8) == b"\x89PNG\r\n\x1a\n"


class TestApiPreloadImages:
    def test_get_small_selection(self, flask_client) -> None:
        r = flask_client.get(
            "/api/preload_images?x=UMAP1&y=UMAP2&selected_rows=0,1,2,3"
        )
        assert r.status_code == 200
        js = r.get_json()
        assert len(js["rows"]) == len(js["images"])

    def test_post_body(self, flask_client) -> None:
        r = flask_client.post(
            "/api/preload_images",
            json={"x": "UMAP1", "y": "UMAP2", "selected_rows": [0, 1, 2, 3]},
        )
        assert r.status_code == 200
        js = r.get_json()
        assert len(js["rows"]) == len(js["images"])

    def test_bad_selected_rows_rejected(self, flask_client) -> None:
        r = flask_client.post(
            "/api/preload_images",
            json={"selected_rows": ["bogus"]},
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
        if explorer_volumes_eligible(dashboard_experiment):
            assert "/trajectory" in body, "nav link '/trajectory' missing from /"
        else:
            assert "CUDA-enabled machine" in body
            assert "Trajectory creator" in body
            assert "/trajectory" not in body


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


# ---------------------------------------------------------------------------
# commands/dashboard.py: CLI argument parsing + ``main`` guard
# ---------------------------------------------------------------------------


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
        monkeypatch.setattr(dash_app, "run_server", lambda **kw: None)
        ns = self._parse(
            ["--filter-max", "123000", "--no-browser", "--command-builder"]
        )
        dash_cli.main(ns)
        assert os.environ.get("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS") == "123000"

    def test_builder_only_with_experiment_view_raises(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(dash_app, "run_server", lambda **kw: None)
        ns = self._parse(["--particle-selection", "--no-browser"])
        with pytest.raises(ValueError, match="need an output directory"):
            dash_cli.main(ns)

    def test_builder_only_with_command_builder_ok(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        called = {}

        def fake_run_server(**kwargs) -> None:
            called.update(kwargs)

        monkeypatch.setattr(dash_app, "run_server", fake_run_server)
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

        monkeypatch.setattr(dash_app, "run_server", fake_run_server)
        ns = self._parse([dashboard_workdir, "--no-browser"])
        dash_cli.main(ns)
        assert called["workdir"] == dashboard_workdir
        assert called["epoch"] == -1

    def test_main_configures_logging_from_verbose(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        called = {}

        def fake_cfg(verbosity: int) -> None:
            called["verbosity"] = verbosity

        monkeypatch.setattr(dash_cli, "_configure_dashboard_logging", fake_cfg)
        monkeypatch.setattr(dash_app, "run_server", lambda **kw: None)
        ns = self._parse(["--command-builder", "--no-browser", "-vv"])
        dash_cli.main(ns)
        assert called["verbosity"] == 2

    def test_default_logging_suppresses_werkzeug_internal(self) -> None:
        dash_cli._configure_dashboard_logging(0)
        assert (
            logging.getLogger("werkzeug._internal").getEffectiveLevel()
            >= logging.WARNING
        )
