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


class TestApiExplorerVolumeMedia:
    def test_ineligible_returns_400(self, flask_client, monkeypatch) -> None:
        monkeypatch.setattr(
            "cryodrgn.dashboard.app.explorer_volumes_eligible", lambda _e: False
        )
        r = flask_client.post(
            "/api/explorer_volume_media",
            json={"rows": [0], "mode": "static"},
        )
        assert r.status_code == 400

    def test_empty_rows_returns_400(self, flask_client, monkeypatch) -> None:
        monkeypatch.setattr(
            "cryodrgn.dashboard.app.explorer_volumes_eligible", lambda _e: True
        )
        r = flask_client.post(
            "/api/explorer_volume_media",
            json={"rows": [], "mode": "static"},
        )
        assert r.status_code == 400

    def test_bad_mode_returns_400(self, flask_client, monkeypatch) -> None:
        monkeypatch.setattr(
            "cryodrgn.dashboard.app.explorer_volumes_eligible", lambda _e: True
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
        assert re.search(r"Build\s+image\s+cache\s+to\s+populate\s+grid", body)
        assert "Build cache" in body
        assert 'id="montage-cache-size-label-text"' in body
        assert 'id="btn-cache-selection-uncached"' in body
        assert "Add 0 selection images to cache" in body
        assert 'id="color-discrete-switches"' in body
        assert "function showGridHighlightsEnabled()" in body

    def test_full_cache_load_suppresses_montage_and_plot_highlight_updates(
        self, flask_client
    ) -> None:
        r = flask_client.get("/explorer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "suppressMontageUpdate: true" in body
        assert "suppressPlotGridHighlights = true" in body


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

    @pytest.mark.parametrize("bad_cache_size", [True, False, 2.5, "2.5"])
    def test_non_integer_cache_sizes_are_rejected(
        self, flask_client, bad_cache_size: object
    ) -> None:
        r = flask_client.post(
            "/api/preload_images",
            json={"x": "UMAP1", "y": "UMAP2", "cache_size": bad_cache_size},
        )
        assert r.status_code == 400
