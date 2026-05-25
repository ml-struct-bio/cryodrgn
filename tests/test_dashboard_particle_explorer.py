"""Particle explorer interface, image preload, and particle preview tests."""

from __future__ import annotations

import base64
import os
import pickle
import re

import matplotlib.pyplot as plt
import numpy as np
import pytest

from cryodrgn.dashboard.context import PRELOAD_CACHE
from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.particle_explorer import (
    _chimerax_render_cmds,
    chimerax_view_matrix_camera_arg,
    format_chimerax_view_matrix_display,
    _config_yaml_path,
    _is_drgnai_config,
    _mpl_retrim_png,
    _sorted_vol_mrc_paths,
    explorer_volumes_eligible,
    montage_cell_label,
)
from cryodrgn.dashboard.preload import (
    DEFAULT_PRELOAD_IMAGE_LIMIT,
    _hybrid_random_knn_spaced_local_indices,
    _preload_cache_time_estimate_bounds,
    encode_particle_batch,
    explorer_cache_size_power10_step,
    explorer_initial_preload_image_limit,
    format_preload_cache_time_hint,
    load_plot_df_rows_from_plot_inds_file,
    montage_bytes,
    particle_thumbnail_b64_from_row,
    sample_plot_df_rows_for_preload,
)


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
        assert any("view #1 orient" in c for c in cmds)

    def test_rotated_view_injects_turn(self) -> None:
        cmds = _chimerax_render_cmds(
            "/tmp/x.mrc", "/tmp/x.png", 100, vol_name="vol000", turn_y=45.0
        )
        assert any("turn y 45.0" in c for c in cmds)
        assert any(c.strip() == "volume center #1" for c in cmds)

    def test_chimerax_volume_color_spec_normalizes_hex(self) -> None:
        from cryodrgn.dashboard.particle_explorer import chimerax_volume_color_spec

        assert chimerax_volume_color_spec("#ff5500") == "#ff5500"
        assert chimerax_volume_color_spec("ff5500") == "#ff5500"
        assert chimerax_volume_color_spec("cornflowerblue") == "cornflowerblue"

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

    def test_volume_color_override(self) -> None:
        cmds = _chimerax_render_cmds(
            "/tmp/x.mrc",
            "/tmp/x.png",
            100,
            vol_name="vol000",
            turn_y=None,
            volume_color="#ff5500",
        )
        assert any("volume color #ff5500" in c for c in cmds)

    def test_view_matrix_camera_arg_parses_display_text(self) -> None:
        text = "camera 1,0,0,0, 0,1,0,0, 0,0,1,0"
        arg = chimerax_view_matrix_camera_arg(text)
        assert arg == "1,0,0,0,0,1,0,0,0,0,1,0"
        assert format_chimerax_view_matrix_display(arg) == "camera " + arg

    def test_view_matrix_camera_arg_rejects_short_input(self) -> None:
        with pytest.raises(ValueError, match="12 numbers"):
            chimerax_view_matrix_camera_arg("camera 1,2,3")

    def test_view_matrix_camera_overrides_turns_in_cmds(self) -> None:
        cmds = _chimerax_render_cmds(
            "/tmp/x.mrc",
            "/tmp/x.png",
            100,
            vol_name="vol000",
            turn_y=None,
            view_turns=[("x", 90.0)],
            view_matrix_camera="1,0,0,0,0,1,0,0,0,0,1,0",
        )
        assert any("view matrix camera 1,0,0,0,0,1,0,0,0,0,1,0" in c for c in cmds)
        assert not any(c.startswith("turn x") for c in cmds)


class TestPreloadTimeHints:
    @pytest.mark.parametrize("cpus", [1, 2, 4, 8, 16, 64])
    def test_bounds_are_ordered_and_positive(self, cpus: int) -> None:
        lo, hi = _preload_cache_time_estimate_bounds(cpus)
        assert 0 < lo <= hi
        assert hi >= 18

    def test_format_singular_core(self) -> None:
        s = format_preload_cache_time_hint(1)
        assert s.endswith("1 core")

    def test_format_plural_cores_embeds_cpus(self) -> None:
        s = format_preload_cache_time_hint(8)
        assert s.endswith("8 cores")
        assert s.startswith("~")


class TestHybridRandomKnnSpacedPreload:
    def test_count_is_bounded_and_unique(self) -> None:
        rng = np.random.default_rng(123)
        coords = rng.normal(size=(500, 2))
        total_k = 200
        picks = _hybrid_random_knn_spaced_local_indices(coords, rng, total_k)
        assert len(picks) <= total_k
        assert len(picks) == len(set(picks))
        assert all(0 <= i < 500 for i in picks)

    def test_empty_coords_returns_empty(self) -> None:
        rng = np.random.default_rng(0)
        assert (
            _hybrid_random_knn_spaced_local_indices(np.zeros((0, 2)), rng, total_k=10)
            == set()
        )

    def test_grid_spacing_takes_one_per_bin_first(
        self,
    ) -> None:
        from cryodrgn.dashboard.preload import _grid_spaced_outlier_pick

        coords = np.array(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [0.0, 1.0],
                [1.0, 1.0],
                [0.5, 0.5],
            ],
            dtype=np.float64,
        )
        order = np.array([4, 0, 1, 2, 3])
        picks = _grid_spaced_outlier_pick(
            coords, order, want=4, exclude=set(), n_bins=2
        )
        assert picks[0] == 4
        assert len(picks) == 4
        assert len(set(picks)) == 4


class TestExplorerCacheSizePower10Step:
    def test_largest_power_of_ten_at_most_five_percent(self) -> None:
        assert explorer_cache_size_power10_step(0) == 1
        assert explorer_cache_size_power10_step(100) == 1
        assert explorer_cache_size_power10_step(1000) == 10
        assert explorer_cache_size_power10_step(10_000) == 100
        assert explorer_cache_size_power10_step(50_000) == 1000
        assert explorer_cache_size_power10_step(100_000) == 1000


class TestExplorerInitialPreloadLimit:
    def test_zero_points(self) -> None:
        assert explorer_initial_preload_image_limit(0) == 1

    def test_default_matches_step_capped_by_plotted_count(self) -> None:
        assert explorer_initial_preload_image_limit(100) == 1
        assert explorer_initial_preload_image_limit(1000) == 10
        assert explorer_initial_preload_image_limit(49_999) == 1000
        assert explorer_initial_preload_image_limit(50_000) == 1000
        assert explorer_initial_preload_image_limit(50_001) == 1000
        assert explorer_initial_preload_image_limit(100_000) == 1000

    def test_small_plotted_count(self) -> None:
        assert explorer_initial_preload_image_limit(5) == 1


class TestSamplePlotDfRowsForPreload:
    def test_default_limit_is_1000(self) -> None:
        assert DEFAULT_PRELOAD_IMAGE_LIMIT == 1000

    def test_restrict_rows_returns_subset(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        restrict = [0, 1, 2, 3, 4]
        rows, ds = sample_plot_df_rows_for_preload(
            dashboard_experiment, "UMAP1", "UMAP2", restrict_to_rows=restrict
        )
        assert set(rows).issubset(set(restrict))
        assert len(rows) == len(ds)

    def test_empty_restrict_returns_empty(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        rows, ds = sample_plot_df_rows_for_preload(
            dashboard_experiment, "UMAP1", "UMAP2", restrict_to_rows=[]
        )
        assert rows == [] and ds == []

    def test_max_images_limits_sample_size(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        rows, ds = sample_plot_df_rows_for_preload(
            dashboard_experiment, "UMAP1", "UMAP2", max_images=2
        )
        assert len(rows) <= 2
        assert len(ds) == len(rows)

    def test_exclude_rows_are_not_resampled(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        rows, _ = sample_plot_df_rows_for_preload(
            dashboard_experiment,
            "UMAP1",
            "UMAP2",
            restrict_to_rows=[0, 1, 2, 3, 4],
            exclude_rows=[0, 1, 2],
        )
        assert rows
        assert set(rows).isdisjoint({0, 1, 2})


class TestMontageBytes:
    def test_empty_returns_hint_png(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        png = montage_bytes(dashboard_experiment, [])
        assert png[:8] == b"\x89PNG\r\n\x1a\n"

    def test_rows_returns_png(self, dashboard_experiment: DashboardExperiment) -> None:
        png = montage_bytes(dashboard_experiment, [0, 1, 2, 3])
        assert png[:8] == b"\x89PNG\r\n\x1a\n"
        assert len(png) > 1000  # not a trivial image

    def test_too_many_rows_are_capped(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        n = dashboard_experiment.z.shape[0]
        png_capped = montage_bytes(dashboard_experiment, list(range(min(n, 50))))
        assert png_capped[:8] == b"\x89PNG\r\n\x1a\n"


class TestParticleThumbnailB64FromRow:
    def test_returns_base64_jpeg(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        b64 = particle_thumbnail_b64_from_row(dashboard_experiment, 0, max_side=48)
        assert isinstance(b64, str) and len(b64) > 0
        raw = base64.standard_b64decode(b64)
        # JPEG magic (SOI marker).
        assert raw[:3] == b"\xff\xd8\xff"


class TestEncodeParticleBatch:
    def test_matches_request_count(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        e = dashboard_experiment
        global_idx = [int(e.all_indices[i]) for i in (0, 1, 2)]
        jpegs = encode_particle_batch(e.particles_path, e.datadir, global_idx, 48)
        assert len(jpegs) == 3
        for b64 in jpegs:
            raw = base64.standard_b64decode(b64)
            assert raw[:3] == b"\xff\xd8\xff"


class TestLoadPlotDfRowsFromPlotIndsFile:
    def test_empty_path_returns_empty(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        assert load_plot_df_rows_from_plot_inds_file(dashboard_experiment, "") == []
        assert load_plot_df_rows_from_plot_inds_file(dashboard_experiment, None) == []

    def test_missing_file_returns_empty(
        self, dashboard_experiment: DashboardExperiment, tmp_path
    ) -> None:
        assert (
            load_plot_df_rows_from_plot_inds_file(
                dashboard_experiment, str(tmp_path / "no_such.pkl")
            )
            == []
        )

    def test_roundtrip_via_pickle(
        self, dashboard_experiment: DashboardExperiment, tmp_path
    ) -> None:
        e = dashboard_experiment
        wanted_ds = np.asarray(e.all_indices[[2, 4, 6]], dtype=int)
        pkl = tmp_path / "plot_inds.pkl"
        with open(pkl, "wb") as fh:
            pickle.dump(wanted_ds, fh)
        rows = load_plot_df_rows_from_plot_inds_file(e, str(pkl))
        assert sorted(rows) == [2, 4, 6]


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
            "cryodrgn.dashboard.particle_explorer.torch_cuda_available", lambda: True
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
            "cryodrgn.dashboard.particle_explorer.torch_cuda_available", lambda: True
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
            "cryodrgn.dashboard.particle_explorer.torch_cuda_available", lambda: False
        )
        assert not explorer_volumes_eligible(dashboard_experiment)


class TestMontageCellLabel:
    def test_matches_particle_explorer_linear_scheme(self) -> None:
        assert montage_cell_label(0) == "A"
        assert montage_cell_label(22) == "Z"
        assert montage_cell_label(23) == "AA"
        assert montage_cell_label(24) == "AB"


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


class TestApiCovariateThresholdRows:
    def test_threshold_api_returns_rows(self, flask_client) -> None:
        r = flask_client.post(
            "/api/covariate_threshold_rows",
            json={"column": "UMAP1", "level": 1e30, "use_max": False},
        )
        assert r.status_code == 200
        js = r.get_json()
        assert isinstance(js["rows"], list)
        assert js["n"] == len(js["rows"])

    def test_use_max_selects_points_below_level(self, flask_client) -> None:
        r = flask_client.post(
            "/api/covariate_threshold_rows",
            json={"column": "UMAP1", "level": 1e30, "use_max": True},
        )
        assert r.status_code == 200
        js = r.get_json()
        assert js["n"] == 100

    def test_invalid_column_rejected(self, flask_client) -> None:
        r = flask_client.post(
            "/api/covariate_threshold_rows",
            json={"column": "not_a_column", "level": 0.0, "use_max": False},
        )
        assert r.status_code == 400

    def test_index_column_rejected(self, flask_client) -> None:
        r = flask_client.post(
            "/api/covariate_threshold_rows",
            json={"column": "index", "level": 0.0, "use_max": False},
        )
        assert r.status_code == 400

    def test_invalid_level_rejected(self, flask_client) -> None:
        r = flask_client.post(
            "/api/covariate_threshold_rows",
            json={"column": "UMAP1", "level": "nope", "use_max": False},
        )
        assert r.status_code == 400

    def test_range_selects_inside_interval(self, flask_client) -> None:
        r = flask_client.post(
            "/api/covariate_threshold_rows",
            json={
                "column": "UMAP1",
                "range_min": -1e30,
                "range_max": 1e30,
                "invert_range": False,
            },
        )
        assert r.status_code == 200
        js = r.get_json()
        assert js["n"] == len(js["rows"])
        assert js["n"] == 100

    def test_range_invert_outside_full_interval_is_empty(self, flask_client) -> None:
        r = flask_client.post(
            "/api/covariate_threshold_rows",
            json={
                "column": "UMAP1",
                "range_min": -1e30,
                "range_max": 1e30,
                "invert_range": True,
            },
        )
        assert r.status_code == 200
        js = r.get_json()
        assert js["n"] == 0

    def test_range_invalid_bounds_rejected(self, flask_client) -> None:
        r = flask_client.post(
            "/api/covariate_threshold_rows",
            json={
                "column": "UMAP1",
                "range_min": "x",
                "range_max": 1.0,
                "invert_range": False,
            },
        )
        assert r.status_code == 400


class TestApiExplorerVolumeMedia:
    def test_ineligible_returns_400(self, flask_client, monkeypatch) -> None:
        monkeypatch.setattr(
            "cryodrgn.dashboard.routes_explorer.explorer_volumes_eligible",
            lambda _e: False,
        )
        r = flask_client.post(
            "/api/explorer_volume_media",
            json={"rows": [0], "mode": "static"},
        )
        assert r.status_code == 400

    def test_empty_rows_returns_400(self, flask_client, monkeypatch) -> None:
        monkeypatch.setattr(
            "cryodrgn.dashboard.routes_explorer.explorer_volumes_eligible",
            lambda _e: True,
        )
        r = flask_client.post(
            "/api/explorer_volume_media",
            json={"rows": [], "mode": "static"},
        )
        assert r.status_code == 400

    def test_bad_mode_returns_400(self, flask_client, monkeypatch) -> None:
        monkeypatch.setattr(
            "cryodrgn.dashboard.routes_explorer.explorer_volumes_eligible",
            lambda _e: True,
        )
        r = flask_client.post(
            "/api/explorer_volume_media",
            json={"rows": [0, 1], "mode": "movie"},
        )
        assert r.status_code == 400
        assert "mode" in (r.get_json() or {}).get("error", "").lower()


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

    def test_post_requested_cache_size(self, flask_client) -> None:
        r = flask_client.post(
            "/api/preload_images",
            json={"x": "UMAP1", "y": "UMAP2", "cache_size": 2},
        )
        assert r.status_code == 200
        js = r.get_json()
        assert len(js["rows"]) <= 2
        assert len(js["rows"]) == len(js["images"])

    def test_larger_cache_size_appends_to_existing_cache(self, flask_client) -> None:
        PRELOAD_CACHE.clear()
        r1 = flask_client.post(
            "/api/preload_images",
            json={"x": "UMAP1", "y": "UMAP2", "cache_size": 2},
        )
        assert r1.status_code == 200
        first_rows = r1.get_json()["rows"]
        assert len(first_rows) <= 2
        assert len(PRELOAD_CACHE) == 1

        r2 = flask_client.post(
            "/api/preload_images",
            json={"x": "UMAP1", "y": "UMAP2", "cache_size": 4},
        )
        assert r2.status_code == 200
        second_rows = r2.get_json()["rows"]
        assert second_rows[: len(first_rows)] == first_rows
        key = next(iter(PRELOAD_CACHE))
        assert PRELOAD_CACHE[key][0] == second_rows

    def test_initial_rows_separate_cache_from_full_pool(self, flask_client) -> None:
        PRELOAD_CACHE.clear()
        r1 = flask_client.post(
            "/api/preload_images",
            json={
                "x": "UMAP1",
                "y": "UMAP2",
                "cache_size": 2,
                "initial_rows": [0, 1],
            },
        )
        assert r1.status_code == 200
        first_rows = r1.get_json()["rows"]
        assert set(first_rows).issubset({0, 1})

        r2 = flask_client.post(
            "/api/preload_images",
            json={"x": "UMAP1", "y": "UMAP2", "cache_size": 4},
        )
        assert r2.status_code == 200
        assert len(PRELOAD_CACHE) == 2

    def test_initial_rows_extends_same_restriction_cache(self, flask_client) -> None:
        PRELOAD_CACHE.clear()
        r1 = flask_client.post(
            "/api/preload_images",
            json={
                "x": "UMAP1",
                "y": "UMAP2",
                "cache_size": 2,
                "initial_rows": [0, 1],
            },
        )
        assert r1.status_code == 200
        first_rows = r1.get_json()["rows"]
        assert len(PRELOAD_CACHE) == 1

        r2 = flask_client.post(
            "/api/preload_images",
            json={
                "x": "UMAP1",
                "y": "UMAP2",
                "cache_size": 4,
                "initial_rows": [0, 1],
            },
        )
        assert r2.status_code == 200
        second_rows = r2.get_json()["rows"]
        assert second_rows[: len(first_rows)] == first_rows
        assert len(PRELOAD_CACHE) == 1

    def test_bad_selected_rows_rejected(self, flask_client) -> None:
        r = flask_client.post(
            "/api/preload_images",
            json={"selected_rows": ["bogus"]},
        )
        assert r.status_code == 400

    def test_bad_cache_size_rejected(self, flask_client) -> None:
        r = flask_client.post(
            "/api/preload_images",
            json={"x": "UMAP1", "y": "UMAP2", "cache_size": 0},
        )
        assert r.status_code == 400

    def test_single_response_caps_thumbnail_batch_size(self, flask_client) -> None:
        from cryodrgn.dashboard.preload import MAX_PRELOAD_IMAGES_PER_HTTP_RESPONSE

        PRELOAD_CACHE.clear()
        r = flask_client.post(
            "/api/preload_images",
            json={"x": "UMAP1", "y": "UMAP2", "cache_size": 5000},
        )
        assert r.status_code == 200
        js = r.get_json()
        assert len(js["rows"]) <= MAX_PRELOAD_IMAGES_PER_HTTP_RESPONSE
        assert len(js["rows"]) == len(js["images"])

    def test_restrict_to_scatter_limits_pool(
        self,
        flask_client,
        dashboard_experiment: DashboardExperiment,
    ) -> None:
        from cryodrgn.dashboard.plots import plot_df_row_indices_for_explorer_scatter

        PRELOAD_CACHE.clear()
        cap = 12
        r = flask_client.post(
            "/api/preload_images",
            json={
                "x": "UMAP1",
                "y": "UMAP2",
                "cache_size": 50,
                "restrict_to_scatter_plot": True,
                "scatter_max_points": cap,
            },
        )
        assert r.status_code == 200
        js = r.get_json()
        pool = set(
            plot_df_row_indices_for_explorer_scatter(dashboard_experiment.plot_df, cap)
        )
        assert pool
        assert set(js["rows"]).issubset(pool)

    def test_selected_rows_only_off_scatter_subsample_returns_no_images(
        self,
        flask_client,
        dashboard_experiment: DashboardExperiment,
    ) -> None:
        """Selection-only preload intersects with scatter pool; off-pool rows encode nothing.

        Matches chunked ``Add selection images`` when some picks are not in the
        downsampled scatter — the UI should plateau instead of erroring.
        """
        from cryodrgn.dashboard.plots import plot_df_row_indices_for_explorer_scatter

        PRELOAD_CACHE.clear()
        cap = 10
        pool = set(
            plot_df_row_indices_for_explorer_scatter(dashboard_experiment.plot_df, cap)
        )
        assert len(pool) == cap
        n = len(dashboard_experiment.plot_df)
        off_scatter = next(i for i in range(n) if i not in pool)

        r = flask_client.post(
            "/api/preload_images",
            json={
                "x": "UMAP1",
                "y": "UMAP2",
                "cache_size": 5,
                "restrict_to_scatter_plot": True,
                "scatter_max_points": cap,
                "selected_rows": [off_scatter],
            },
        )
        assert r.status_code == 200
        js = r.get_json()
        assert js["rows"] == []
        assert js["images"] == []


class TestParticleExplorerTemplateRegressions:
    def test_particle_page_exposes_cache_progress_and_image_grid_shell(
        self, flask_client
    ) -> None:
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert 'id="image-cache-progress"' in body
        assert 'role="progressbar"' in body
        assert 'id="image-grid-menu-toggle"' in body
        assert "Image grid" in body
        assert re.search(
            r"Load\s+particle\s+images\s+into\s+cache\s+to\s+display\s+in\s+grid",
            body,
        )
        assert 'id="cryo-explorer-image-cache-fieldset"' in body
        assert 'id="cryo-explorer-cache-load-body"' in body
        assert 'id="cryo-explorer-cache-load-progress-wrap"' in body
        assert "Loading particle images into cache" in body
        assert 'id="btn-expand-cache"' in body
        assert "Build new" in body
        assert 'id="btn-clear-image-cache"' in body
        assert 'id="montage-cache-size-label-text"' in body
        assert 'id="btn-cache-selection-uncached"' in body
        assert "Add 0 selection images to cache" in body
        assert 'id="color-discrete-switches"' in body
        assert "js/cryo_cc_legend_primitives.js" in body
        assert "js/color_covariate_legend.js" in body
        p_idx = body.find("js/cryo_cc_legend_primitives.js")
        c_idx = body.find("js/color_covariate_legend.js")
        assert p_idx != -1 and c_idx != -1 and p_idx < c_idx
        assert "new CryoColorCovariateLegend" in body
        assert "function showGridHighlightsEnabled()" in body

    def test_explorer_legend_static_assets_served_by_flask(self, flask_client) -> None:
        """Wheel/sdist must ship nested ``static/js/*`` (``static/*`` alone omits them)."""
        r_prim = flask_client.get("/static/js/cryo_cc_legend_primitives.js")
        assert r_prim.status_code == 200, r_prim.status_code
        assert b"CryoCcLegendPrimitives" in r_prim.data
        r_leg = flask_client.get("/static/js/color_covariate_legend.js")
        assert r_leg.status_code == 200, r_leg.status_code
        assert b"CryoColorCovariateLegend" in r_leg.data

    def test_lasso_box_selection_union_and_deselect_preserves(
        self, flask_client
    ) -> None:
        """Selection UX regression checks for disjoint lasso/box drags.

        We can't run the browser JS here, so we assert the template contains
        the key behaviours:
        - lasso/box selections accumulate (union) across multiple drags
        - Plotly deselect events do not wipe the stored selection
          (only the explicit "Clear selection" button should do that)
        - lasso/box accumulation only happens when the previous selection mode
          was already geometric/``lasso`` (range/toggle selections overwrite)
        """
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)

        # Accumulation logic is guarded by "prevMode === 'lasso'" so lasso clears
        # range/toggle selection instead of unioning with it.
        assert 'prevMode === "lasso"' in body
        assert "var mergedSet = new Set(baseRows);" in body
        assert "Re-apply accumulated union after region geometry is committed" in body

        # Range/toggle selections overwrite the current selection (they don't
        # accumulate with lasso drags).
        assert "saveSelectionRows = Array.from(rowSet);" in body
        assert 'statusPrefix === "Color threshold"' in body
        assert 'statusPrefix === "Discrete"' in body

        start = body.find('gd.on("plotly_deselect", function()')
        assert start != -1, "missing plotly_deselect handler"
        end = body.find('gd.addEventListener("mouseleave"', start)
        assert end != -1 and end > start, "could not bound plotly_deselect block"
        deselect_block = body[start:end]

        assert (
            'Keep the accumulated lasso/box selection unless the explicit "Clear selection"'
            in deselect_block
        )
        # Previous behaviour wiped selection state here; we intentionally keep it.
        assert "saveSelectionRows = []" not in deselect_block

        # Explicit Clear selection button should still wipe the stored selection.
        clear_start = body.find("function clearExplorerSelection()")
        assert clear_start != -1, "missing clearExplorerSelection()"
        clear_end = body.find(
            "function rowsMatchingColorThresholdFromTrace", clear_start
        )
        assert (
            clear_end != -1 and clear_end > clear_start
        ), "could not bound clear block"
        clear_block = body[clear_start:clear_end]
        assert (
            "applyRowsSelection([], undefined, LASSO_SELECTION_DEBOUNCE_MS + 120);"
            in clear_block
        )
        assert "clearScatterGeometricSelection();" in clear_block
        assert (
            "syncCommittedScatterRegionOverlays({ clearSelections: true });"
            in clear_block
        )
        assert "pickedClearMontage" in clear_block
        assert "nRegClear" in body
        assert (
            "syncMontageResampleFromCacheButton();\n"
            "    updateParticleSelFieldset();" in clear_block
        )

        # Scattergl can leave lasso/box paint until dragmode is nudged and transient shapes are stripped.
        assert "pulseScatterglDragmodeToFlushSelectionPaint" in body
        assert "scheduleScatterglSelectionPaintFlushAfterClear" in body
        assert "scatterExplorerShapesPurgeForGeometryClear" in body
        assert "shapeMatchesCryoCommittedScatterRegion" in body
        assert "plotlyRelayoutShapesHardReplaceThenPatch" in body
        assert "setMontageCellSelectionBorderStyle" in body
        assert "2px 3px 3px 3px" in body

    def test_full_cache_load_suppresses_montage_and_plot_highlight_updates(
        self, flask_client
    ) -> None:
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "suppressMontageUpdate: true" in body
        assert "suppressPlotGridHighlights = true" in body
        assert "IMAGE_CACHE_HTTP_CHUNK_MAX" in body
        assert "preloadFetchErrorMessage" in body

    def test_scatter_double_click_replaces_montage_slot_a_without_selection(
        self, flask_client
    ) -> None:
        """Double-click assigns the particle to montage A and may grow the cache."""
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "handleScatterPointDoubleClick" in body
        assert "replaceMontageSlotAt" in body
        assert "paintMontageCellAtIndex" in body
        assert "assignScatterRowToMontageSlotA" in body
        assert "scatterBlockSingleClickUntil" in body
        assert "montageResampleSuppressed" in body
        assert "armMontageResampleSuppress" in body
        lasso_snap = body.find("function applyLassoSelectionFromSnapshot()")
        assert lasso_snap != -1
        lasso_block = body[lasso_snap : lasso_snap + 1200]
        assert "montageResampleSuppressed()" in lasso_block
        assert "suppressMontageRefresh: true" in body
        assert "updateMontageOrdered" in body
        assert "suppressSelectionAfterScatterDoubleClick" in body
        click_start = body.find('gd.on("plotly_click", function(ev)')
        assert click_start != -1
        click_block = body[click_start : click_start + 2200]
        assert "handleScatterPointDoubleClick(pt)" in click_block
        assert "updateMontage(nbs)" in click_block

    def test_queue_highlight_restyle_merges_pending_xy_update_with_styling_patch(
        self, flask_client
    ) -> None:
        """Regression test for montage cache resample drift.

        `updateMontage()` enqueues a full grid-highlight restyle payload (x/y + styling).
        A later call like `refreshGridHighlightMarkerStylesFromLastRows()` enqueues a
        styling-only patch; if the pending payload is overwritten, the highlighted points
        drift because the x/y update is dropped.
        """
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "function queueHighlightRestyle(restyleData)" in body
        assert "highlightRestyleRaf != null && pendingHighlightRestyle" in body
        assert "pendingHighlightRestyle[k] = restyleData[k];" in body

    def test_grid_letter_highlights_constant_opacity_and_covariate_outline(
        self, flask_client
    ) -> None:
        """Grid-letter overlay uses fixed marker opacity; selection is fill vs hollow only."""
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "cryo-grid-highlight-marker-policy" in body
        assert "GRID_HIGHLIGHT_MARKER_OPACITY" in body
        assert "appendGridHighlightMarkerRestyle" in body
        assert "GRID_HIGHLIGHT_CLEAR_FILL" in body
        assert "refreshGridHighlightMarkerStylesFromLastRows" in body
        assert "treatAllGridPointsAsSelected" in body
        assert "cryo-grid-highlight-no-dim" in body
        assert 'type: "scatter"' in body
        assert "multiGeom" in body
        assert 'restyleData["textfont.color"]' in body

    def test_committed_scatter_shapes_use_between_layer_for_grid_letters(
        self, flask_client
    ) -> None:
        """Committed lasso/box fills use Plotly layer "between" so letter markers draw on top."""
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert 'layer: "between"' in body
        assert "CDRGN_COMMIT_REGION_LINE_WIDTH" in body

    def test_multi_region_lasso_overlay_chips_solo_and_colour_wheel(
        self, flask_client
    ) -> None:
        """Disjoint regions use HTML overlay chips with solo ``1`` + colour wheel (not Plotly text)."""
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert 'id="scatter-region-chips-overlay"' in body
        assert "soloCommittedScatterRegion" in body
        assert "__cdrgnScatterRegion:" in body
        assert "scatterRegionChipPointerDown" in body
        assert "cryo-explorer-scatter-region-chip__actions" in body
        assert "scheduleScatterRegionLabelChipsSync" in body
        assert "removeCommittedScatterRegion" in body
        assert "cryo-explorer-scatter-region-chip__remove" in body
        assert "cryo-explorer-scatter-region-chip__save" in body
        assert "saveCommittedScatterRegion" in body
        assert "openRegionSelectionFileBrowser" in body
        assert "pendingRegionSaveRows" in body

    def test_multi_region_selection_pie_coloured_slices(self, flask_client) -> None:
        """Multi-region pie: per-region slices; total % stays black."""
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "buildMultiRegionSelPieBackground" in body
        assert "applySelectionPieVisual" in body
        assert "cryo-explorer-sel-pie--multi-region" in body
        assert "sel-pie-union-outline" not in body
        assert "updateSelPieUnionOutline" not in body

    def test_dashboard_save_buttons_use_floppy_disk_icon(self, flask_client) -> None:
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "cryo_dashboard_icons.js" in body
        assert "CryoDashboardIcons" in body

    def test_save_floppy_icon_is_outline_reference_style(self) -> None:
        from pathlib import Path

        js = Path("cryodrgn/dashboard/static/js/cryo_dashboard_icons.js").read_text(
            encoding="utf-8"
        )
        assert 'fill=\\"none\\"' in js
        assert 'stroke=\\"currentColor\\"' in js
        assert "M3.75 6.25" in js
        assert 'x=\\"2.6\\"' not in js

    def test_multi_region_overlap_preserves_existing_regions_on_commit(
        self, flask_client
    ) -> None:
        """Overlapping lassos: unchanged regions keep prior rows; no polygon clipping."""
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "rowArr = prevSnap[ci].rows.slice()" in body
        assert "polygon_clipping.umd.js" not in body
        assert "overlayRaw" not in body
        assert "clipScatterRawExcludingNewerRegions" not in body
        assert "Overlapping lassos: each region keeps" in body

    def test_selection_save_file_browser_is_modal_popup(self, flask_client) -> None:
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert 'id="sel-file-browser-panel"' in body
        assert "cryo-explorer-save-modal" in body
        assert 'role="dialog"' in body
        assert "cryo-explorer-save-modal__backdrop" in body
        assert "syncSelectionSaveModalTitle" in body
        assert "cryo-explorer-save-modal-open" in body

    def test_plotly_selected_commits_regions_before_montage_pool_refresh(
        self, flask_client
    ) -> None:
        """Image cache + multi-lasso: rebuild ``committedScatterRegions`` before montage refresh.

        ``updateMontage`` / ``appendGridHighlightMarkerRestyle`` read ``committedScatterRegions``;
        ``applyLassoSelectionFromSnapshot`` must run only after
        ``syncCommittedScatterRegionOverlays``, then grid-letter colours catch up via
        ``refreshGridHighlightMarkerStylesFromLastRows``.
        """
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        start = body.find('gd.on("plotly_selected", function(ev)')
        assert start != -1, "missing plotly_selected handler"
        mid = body.find("lassoSelectionDebounceTimer = setTimeout(function()", start)
        assert mid != -1, "missing debounced plotly_selected block"
        end = body.find("}, LASSO_SELECTION_DEBOUNCE_MS);", mid)
        assert end != -1 and end > mid, "could not bound plotly_selected debounce block"
        debounce_block = body[mid:end]
        sync_marker = "syncCommittedScatterRegionOverlays({ clearSelections: true });"
        snap_marker = "applyLassoSelectionFromSnapshot();"
        refresh_marker = "refreshGridHighlightMarkerStylesFromLastRows();"
        pos_sync = debounce_block.find(sync_marker)
        pos_snap = debounce_block.find(snap_marker)
        pos_refresh = debounce_block.find(refresh_marker)
        assert pos_sync != -1, debounce_block[:400]
        assert pos_snap != -1
        assert pos_refresh != -1
        assert pos_sync < pos_snap < pos_refresh, (
            "expected order: syncCommittedScatterRegionOverlays → "
            "applyLassoSelectionFromSnapshot → refreshGridHighlightMarkerStylesFromLastRows"
        )
        assert "selectionsSnapshotForCommit" in body
        assert "dedupeConsecutiveEqualScatterShapes" in body
        pos_comment = debounce_block.find(
            "Re-apply accumulated union after region geometry is committed"
        )
        pos_restyle = debounce_block.find("applyScatterSelectionHighlight(selectedTi);")
        assert pos_comment != -1
        assert (
            pos_restyle != -1
        ), "expected applyScatterSelectionHighlight in debounce block"
        pos_selection_apply = pos_restyle
        assert pos_sync < pos_comment < pos_selection_apply, (
            "region overlays must sync before selection highlight apply "
            "(restyle/overlay can clear layout.selections)"
        )

    def test_multi_region_lasso_combo_uses_selectedpoints_dimming_in_debounce(
        self, flask_client
    ) -> None:
        """Scattergl selection uses compact ``selectedpoints`` restyle (trace map for row→index)."""
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "applyScatterSelectionHighlight" in body
        assert "rowToTraceIndexMap" in body
        assert "rowsUnionFromCommittedScatterRegions" in body
        assert "mergePersistedScatterSelectionShapes" in body
        start = body.find('gd.on("plotly_selected", function(ev)')
        assert start != -1
        mid = body.find("lassoSelectionDebounceTimer = setTimeout(function()", start)
        end = body.find("}, LASSO_SELECTION_DEBOUNCE_MS);", mid)
        debounce_block = body[mid:end]
        assert "applyScatterSelectionHighlight(selectedTi);" in debounce_block
        assert "selectedpoints" in body
        assert "cdrgnSelectionOverlay" not in body

    def test_multi_region_row_membership_recomputed_from_geometry(
        self, flask_client
    ) -> None:
        """Each region's ``rows`` must follow lasso geometry, not only the last ``rowsSnap``."""
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "recomputeCommittedScatterRegionRowsFromGeometry" in body
        assert "rowsUnionFromCommittedScatterRegions" in body
        assert "scatterRawShapeContainsDataXY" in body
        assert "pointInPolygonXY" in body
        assert "traceIndexForPlotDfRow" in body

    def test_multi_region_montage_and_grid_use_scatter_region_line_colour(
        self, flask_client
    ) -> None:
        """Montage + grid-letter styling for no-covariate multi-region tracks ``scatterRegionPlotStyle``."""
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "function selectionRegionMontageStyles(regionIdx)" in body
        assert "scatterRegionPlotStyle(regionIdx).line" in body
        assert "discreteLabelMontageStyles(lineHex)" in body
        assert "borderNoCov = scatterRegionPlotStyle(ridxM).line;" in body
        assert "scatterRegionPlotStyle(rIdxgeom).line" in body
        assert '(multiGeom ? "#94a3b8" : ACCENT)' in body

    def test_scatter_region_overlay_chips_compact_vertical_css(
        self, flask_client
    ) -> None:
        """Region count chips stay short while preserving count font + icon metrics."""
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert (
            ".cryo-dash-page--particle-explorer "
            ".cryo-explorer-scatter-region-chip.cryo-cc-discrete-cell--plastic {"
        ) in body
        assert "inset 0 0 0 1px rgba(255, 255, 255, 0.28)" in body
        assert ".cryo-explorer-scatter-region-chip__row" in body
        assert (
            "cdrgn: chip row — count label vertically centered with action icons"
            in body
        )
        assert "padding: 0.14rem 0.22rem 0.18rem 0.22rem" in body
        assert (
            ".cryo-explorer-scatter-region-chip .cryo-cc-discrete-switch-label" in body
        )
        assert "max-height: 0" in body
        assert "font-size: 0.7rem" in body
        assert (
            ".cryo-explorer-scatter-region-chip .cryo-cc-discrete-switch-count" in body
        )
        assert "font-size: 0.6rem" in body
        assert "line-height: 1" in body
        assert "transform: scale(1.524)" in body
        assert "width: 0.56rem" in body
        assert "height: 0.56rem" in body

    def test_montage_cards_use_top_meta_band_and_tight_image_margins(
        self, flask_client
    ) -> None:
        """Letter + covariate/idx top-aligned in ``cryo-montage-meta``; tight cell padding and zero gap to image."""
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "cryo-montage-meta" in body
        assert "cryo-montage-meta-right" in body
        assert "cryo-montage-cell--light" in body
        assert 'className = "cryo-montage-meta"' in body
        assert 'className = "cryo-montage-meta-right"' in body
        assert (
            "  .cryo-montage-img-wrap {\n"
            "    position: relative;\n"
            "    aspect-ratio: 1;\n"
            "    flex: 1 1 auto;\n"
            "    margin: 0;\n"
            "    min-height: 0;\n"
            "  }"
        ) in body
        assert (
            "  .cryo-montage-cell.cryo-montage-cell--light {\n"
            "    background: var(--paper, #faf8f4);\n"
            "    padding: 1px;\n"
            "    gap: 0;\n"
            "  }"
        ) in body
        assert "cryo-montage-footer" not in body
        assert (
            "  .cryo-montage-meta {\n"
            "    display: flex;\n"
            "    flex-direction: row;\n"
            "    align-items: flex-start;\n"
        ) in body
        assert "#clear-explorer-selection:disabled" in body
        assert 'lbl.style.display = "inline-flex"' in body
        assert 'lbl.style.alignSelf = "flex-start"' in body
        assert "scaleRem(metaSize, 1.1)" in body
        assert "continuousMontageStylesFromT" in body
        assert "letterFontRem * 0.11" in body
        assert "1.22 / labStr.length" in body
        assert 'meta.style.alignItems = "flex-start"' in body
        assert 'meta.style.gap = "0"' in body
        assert "lbl: lbl" in body


class TestPreloadDeltaResponses:
    def test_delta_response_returns_only_new_images_and_total_cached(
        self, flask_client
    ) -> None:
        PRELOAD_CACHE.clear()
        first = flask_client.post(
            "/api/preload_images",
            json={"x": "UMAP1", "y": "UMAP2", "cache_size": 2},
        )
        assert first.status_code == 200, first.get_json()
        first_js = first.get_json()
        assert len(first_js["rows"]) <= 2

        second = flask_client.post(
            "/api/preload_images",
            json={
                "x": "UMAP1",
                "y": "UMAP2",
                "cache_size": 4,
                "response_mode": "delta",
            },
        )
        assert second.status_code == 200, second.get_json()
        second_js = second.get_json()
        assert second_js["total_cached"] >= len(first_js["rows"])
        assert set(second_js["rows"]).isdisjoint(set(first_js["rows"]))
        assert len(second_js["rows"]) == len(second_js["images"])
        assert "batch_elapsed" in second_js
        assert isinstance(second_js["batch_elapsed"], int | float)
        assert second_js["batch_elapsed"] >= 0

    def test_delta_response_empty_when_target_already_cached(
        self, flask_client
    ) -> None:
        PRELOAD_CACHE.clear()
        first = flask_client.post(
            "/api/preload_images",
            json={"x": "UMAP1", "y": "UMAP2", "cache_size": 3},
        )
        assert first.status_code == 200, first.get_json()
        cached_count = first.get_json()["total_cached"]

        second = flask_client.post(
            "/api/preload_images",
            json={
                "x": "UMAP1",
                "y": "UMAP2",
                "cache_size": 3,
                "response_mode": "delta",
            },
        )
        assert second.status_code == 200, second.get_json()
        assert second.get_json()["rows"] == []
        assert second.get_json()["images"] == []
        assert second.get_json()["total_cached"] == cached_count
        assert second.get_json().get("batch_elapsed") == 0.0

    def test_delta_response_empty_when_cache_size_below_server_rows(
        self, flask_client
    ) -> None:
        """Smaller ``cache_size`` than the stored PRELOAD_CACHE row count → empty delta.

        The particle explorer must treat this as a plateau (no client-side error), not
        assume every chunk grows ``cachedImageCount()``.
        """
        PRELOAD_CACHE.clear()
        first = flask_client.post(
            "/api/preload_images",
            json={"x": "UMAP1", "y": "UMAP2", "cache_size": 6},
        )
        assert first.status_code == 200, first.get_json()
        cached_count = first.get_json()["total_cached"]
        assert cached_count >= 1

        second = flask_client.post(
            "/api/preload_images",
            json={
                "x": "UMAP1",
                "y": "UMAP2",
                "cache_size": 2,
                "response_mode": "delta",
            },
        )
        assert second.status_code == 200, second.get_json()
        js = second.get_json()
        assert js["rows"] == []
        assert js["images"] == []
        assert js["total_cached"] == cached_count

    def test_post_invalidate_cache_clears_epoch_preloads(self, flask_client) -> None:
        PRELOAD_CACHE.clear()
        first = flask_client.post(
            "/api/preload_images",
            json={"x": "UMAP1", "y": "UMAP2", "cache_size": 6},
        )
        assert first.status_code == 200, first.get_json()
        assert len(PRELOAD_CACHE) >= 1

        inv = flask_client.post(
            "/api/preload_images",
            json={"invalidate_cache": True},
        )
        assert inv.status_code == 200, inv.get_json()
        assert inv.get_json().get("ok") is True
        assert PRELOAD_CACHE == {}

        again = flask_client.post(
            "/api/preload_images",
            json={"x": "UMAP1", "y": "UMAP2", "cache_size": 3},
        )
        assert again.status_code == 200, again.get_json()
        assert again.get_json()["total_cached"] <= 3

    def test_delta_nothing_new_when_extending_with_same_cap_and_full_cache(
        self, flask_client
    ) -> None:
        PRELOAD_CACHE.clear()
        first = flask_client.post(
            "/api/preload_images",
            json={"x": "UMAP1", "y": "UMAP2", "cache_size": 4},
        )
        assert first.status_code == 200, first.get_json()
        n0 = first.get_json()["total_cached"]

        second = flask_client.post(
            "/api/preload_images",
            json={
                "x": "UMAP1",
                "y": "UMAP2",
                "cache_size": 8,
                "response_mode": "delta",
            },
        )
        assert second.status_code == 200, second.get_json()
        n1 = second.get_json()["total_cached"]
        assert n1 >= n0

        third = flask_client.post(
            "/api/preload_images",
            json={
                "x": "UMAP1",
                "y": "UMAP2",
                "cache_size": 8,
                "response_mode": "delta",
            },
        )
        assert third.status_code == 200, third.get_json()
        js = third.get_json()
        assert js["rows"] == []
        assert js["images"] == []
        assert js["total_cached"] == n1

    @pytest.mark.parametrize("bad_cache_size", [True, False, 2.5, "2.5"])
    def test_non_integer_cache_sizes_are_rejected(
        self, flask_client, bad_cache_size: object
    ) -> None:
        r = flask_client.post(
            "/api/preload_images",
            json={"x": "UMAP1", "y": "UMAP2", "cache_size": bad_cache_size},
        )
        assert r.status_code == 400
