"""Particle explorer interface, image preload, and particle preview tests."""

from __future__ import annotations

import base64
import os
import pickle
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from cryodrgn.dashboard.context import (
    PRELOAD_CACHE,
    EXPERIMENT_STORE,
    clear_experiment_caches,
)
from cryodrgn.dashboard.data import DashboardExperiment
from cryodrgn.dashboard.particle_explorer import (
    _chimerax_render_cmds,
    _coerce_z_decode_matrix,
    _config_yaml_path,
    _decode_z_values_to_vol_paths,
    _is_drgnai_config,
    _mpl_retrim_png,
    _register_vol_mrc_cache,
    _sorted_vol_mrc_paths,
    _VOL_CACHE_LOCK,
    _VOL_MRC_CACHE,
    chimerax_view_matrix_camera_arg,
    explorer_volumes_eligible,
    format_chimerax_view_matrix_display,
    generate_montage_volume_pngs,
    generate_trajectory_volume_pngs,
    montage_cell_label,
    save_cached_volumes_to_dir,
    torch_cuda_available,
    volume_cell_gif_from_cache,
)
from tests.conftest import (
    _monkeypatch_explorer_volumes_eligible,
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

pytestmark = pytest.mark.dashboard


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
        assert any("volume #1 sdLevel 2" in c for c in cmds)
        assert "surface dust all size 10" in cmds
        assert "lighting soft" in cmds

    def test_volume_level_override(self) -> None:
        cmds = _chimerax_render_cmds(
            "/tmp/x.mrc",
            "/tmp/x.png",
            100,
            vol_name="vol000",
            turn_y=None,
            volume_level=0.42,
        )
        assert any("volume #1 level 0.42" in c for c in cmds)
        assert not any("sdLevel" in c for c in cmds)

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
        assert not any("view #1 orient" in c for c in cmds)

    def test_vector_axis_turn_emits_single_turn_cmd(self) -> None:
        cmds = _chimerax_render_cmds(
            "/tmp/x.mrc",
            "/tmp/x.png",
            100,
            vol_name="vol000",
            turn_y=None,
            view_turns=[("0.4,0.8,0.4472", 27.0)],
        )
        assert any("view #1 orient" in c for c in cmds)
        turn_cmds = [c for c in cmds if c.startswith("turn ")]
        assert len(turn_cmds) == 1, turn_cmds
        # Axis is renormalised to a unit vector and kept as a comma-separated triple.
        assert turn_cmds[0].startswith("turn ")
        axis = turn_cmds[0].split()[1]
        comps = [float(v) for v in axis.split(",")]
        assert len(comps) == 3
        assert abs(sum(c * c for c in comps) - 1.0) < 1e-4


class TestPreloadTimeHints:
    @pytest.mark.parametrize("cpus", [1, 64])
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

    def test_rejects_negative_index(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            montage_cell_label(-1)


class TestTorchCudaAvailable:
    def test_import_error_returns_false(self, monkeypatch: pytest.MonkeyPatch) -> None:
        import builtins

        real_import = builtins.__import__

        def _fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "torch":
                raise ImportError("no torch")
            return real_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", _fake_import)
        assert torch_cuda_available() is False


class TestExplorerVolumesEligibleExtra:
    def test_false_when_tilt_enc_mode(
        self, dashboard_experiment: DashboardExperiment, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from dataclasses import replace

        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer.torch_cuda_available", lambda: True
        )
        exp = replace(dashboard_experiment, enc_mode="tilt")
        assert not explorer_volumes_eligible(exp)


@pytest.fixture(autouse=True)
def _clear_particle_explorer_vol_cache() -> None:
    from cryodrgn.dashboard.particle_explorer import _vol_cache_evict_unlocked

    with _VOL_CACHE_LOCK:
        for tok in list(_VOL_MRC_CACHE):
            _vol_cache_evict_unlocked(tok)
    yield
    with _VOL_CACHE_LOCK:
        for tok in list(_VOL_MRC_CACHE):
            _vol_cache_evict_unlocked(tok)


class TestVolumeCaches:
    """Direct coverage for ``cryodrgn.dashboard.volume_caches`` (no Flask)."""

    def test_mrc_cache_ttl_evicts_and_removes_dir(self, tmp_path: Path) -> None:
        from cryodrgn.dashboard.volume_caches import VolumeMrcCache

        cache = VolumeMrcCache(max_entries=8, ttl_s=0.01)
        mrc_dir = tmp_path / "vols"
        mrc_dir.mkdir()
        vol = mrc_dir / "vol_001.mrc"
        vol.write_bytes(b"\x00")
        token = cache.register(str(mrc_dir), [str(vol)], (0,))
        assert cache.get_meta(token) is not None
        import time

        time.sleep(0.02)
        assert cache.get_meta(token) is None
        assert not mrc_dir.exists()

    def test_mrc_cache_max_entries_evicts_oldest_dir(self, tmp_path: Path) -> None:
        from cryodrgn.dashboard.volume_caches import VolumeMrcCache

        cache = VolumeMrcCache(max_entries=2, ttl_s=3600.0)
        tokens: list[str] = []
        dirs: list[Path] = []
        for i in range(3):
            d = tmp_path / f"mrc_{i}"
            d.mkdir()
            p = d / "vol_001.mrc"
            p.write_bytes(b"\x00")
            dirs.append(d)
            tokens.append(cache.register(str(d), [str(p)], (i,)))
        assert len(cache.entries) == 2
        assert tokens[0] not in cache.entries
        assert not dirs[0].exists()
        assert dirs[1].exists() and dirs[2].exists()

    def test_mrc_cache_evict_unlocked_removes_dir(self, tmp_path: Path) -> None:
        from cryodrgn.dashboard.volume_caches import VolumeMrcCache

        cache = VolumeMrcCache()
        mrc_dir = tmp_path / "gone"
        mrc_dir.mkdir()
        vol = mrc_dir / "vol_001.mrc"
        vol.write_bytes(b"\x00")
        token = cache.register(str(mrc_dir), [str(vol)], (0,))
        with cache.lock:
            cache.evict_unlocked(token)
        assert token not in cache.entries
        assert not mrc_dir.exists()

    def test_job_store_progress_snapshot_and_ttl(self) -> None:
        from cryodrgn.dashboard.volume_caches import VolumeJobStore

        store = VolumeJobStore(progress_ttl_s=0.01, partial_ttl_s=0.01)
        store.progress_register("tok", total=10, workers=2, phase="decode")
        store.progress_set_done("tok", 4)
        snap = store.progress_snapshot("tok")
        assert snap is not None
        assert snap["done"] == 4
        assert snap["total"] == 10
        assert snap["n_gpus"] == 2
        import time

        time.sleep(0.02)
        assert store.progress_snapshot("tok") is None

    def test_volume_generation_token_reuse_via_register(
        self,
        dashboard_experiment: DashboardExperiment,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Two montage generations register distinct tokens; epoch switch clears experiments."""

        def _fake_decode(
            exp: DashboardExperiment, z_values: np.ndarray, mrc_dir: str
        ) -> list[str]:
            os.makedirs(mrc_dir, exist_ok=True)
            paths: list[str] = []
            for i in range(len(z_values)):
                p = os.path.join(mrc_dir, f"vol_{i + 1:03d}.mrc")
                Path(p).write_bytes(b"\x00")
                paths.append(p)
            return paths

        def _fake_pngs(tasks, chimerax_cpus=1):
            from PIL import Image

            paths = []
            for _idx, _mrc, png, _dpi in tasks:
                Image.new("RGB", (4, 4), color=(9, 8, 7)).save(png)
                paths.append(png)
            return paths

        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer._decode_z_values_to_vol_paths",
            _fake_decode,
        )
        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer.render_static_pngs_parallel",
            _fake_pngs,
        )
        _pngs1, tok1 = generate_montage_volume_pngs(
            dashboard_experiment, [0], chimerax_cpus=1
        )
        _pngs2, tok2 = generate_montage_volume_pngs(
            dashboard_experiment, [1], chimerax_cpus=1
        )
        assert tok1 != tok2
        assert tok1 in _VOL_MRC_CACHE and tok2 in _VOL_MRC_CACHE

        EXPERIMENT_STORE.clear_all()
        EXPERIMENT_STORE.get_experiment(
            dashboard_experiment.workdir, dashboard_experiment.epoch, -1
        )
        assert EXPERIMENT_STORE.experiments
        clear_experiment_caches()
        assert not EXPERIMENT_STORE.experiments
        # Volume MRC tokens are process-local and intentionally survive epoch cache clears.
        assert tok1 in _VOL_MRC_CACHE


class TestParticleExplorerVolumeCache:
    def test_generate_montage_and_gif_from_cache(
        self,
        dashboard_experiment: DashboardExperiment,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        def _fake_decode(
            exp: DashboardExperiment, z_values: np.ndarray, mrc_dir: str
        ) -> list[str]:
            os.makedirs(mrc_dir, exist_ok=True)
            paths: list[str] = []
            for i in range(len(z_values)):
                p = os.path.join(mrc_dir, f"vol_{i + 1:03d}.mrc")
                Path(p).write_bytes(b"\x00")
                paths.append(p)
            return paths

        def _fake_pngs(tasks, chimerax_cpus=1):
            from PIL import Image

            paths = []
            for _idx, _mrc, png, _dpi in tasks:
                Image.new("RGB", (6, 6), color=(1, 2, 3)).save(png)
                paths.append(png)
            return paths

        def _fake_gif(mrc_path: str, out_gif: str, **kwargs) -> None:
            Path(out_gif).write_bytes(b"GIF89a")

        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer._decode_z_values_to_vol_paths",
            _fake_decode,
        )
        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer.render_static_pngs_parallel",
            _fake_pngs,
        )
        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer.render_rotating_gif",
            _fake_gif,
        )
        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer.cx.render_rotating_gif",
            _fake_gif,
        )
        pngs, token = generate_montage_volume_pngs(
            dashboard_experiment, [0, 1], chimerax_cpus=2
        )
        assert len(pngs) == 2
        assert pngs[0][:8] == b"\x89PNG\r\n\x1a\n"
        gif = volume_cell_gif_from_cache(
            token, 0, rows_expected=(0, 1), gif_frames=6, chimerax_cpus=2
        )
        assert gif.startswith(b"GIF89a")

    def test_save_cached_volumes_to_dir(
        self,
        dashboard_experiment: DashboardExperiment,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        mrc_dir = tmp_path / "mrcs"
        mrc_dir.mkdir()
        vols = [mrc_dir / "vol_001.mrc", mrc_dir / "vol_002.mrc"]
        for p in vols:
            p.write_bytes(b"\x00")
        token = _register_vol_mrc_cache(str(mrc_dir), [str(p) for p in vols], (0, 1))
        out = tmp_path / "saved"
        saved = save_cached_volumes_to_dir(token, str(out), filename_prefix="cell")
        assert len(saved) == 2
        assert all(os.path.isfile(p) for p in saved)
        assert saved[0].endswith("cell_001.mrc")

    def test_volume_cache_validation_errors(self, tmp_path: Path) -> None:
        mrc_dir = tmp_path / "one"
        mrc_dir.mkdir()
        vol = mrc_dir / "vol_001.mrc"
        vol.write_bytes(b"\x00")
        token = _register_vol_mrc_cache(str(mrc_dir), [str(vol)], (0,))
        with pytest.raises(ValueError, match="Unknown or expired"):
            volume_cell_gif_from_cache("bad-token", 0, rows_expected=(0,))
        with pytest.raises(ValueError, match="do not match"):
            volume_cell_gif_from_cache(token, 0, rows_expected=(1,))
        with pytest.raises(ValueError, match="out of range"):
            volume_cell_gif_from_cache(token, 3, rows_expected=(0,))
        with pytest.raises(ValueError, match="Choose an output folder"):
            save_cached_volumes_to_dir(token, "")

    def test_cache_expiry_evicts_entry(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        mrc_dir = tmp_path / "mrcs"
        mrc_dir.mkdir()
        vol = mrc_dir / "vol_001.mrc"
        vol.write_bytes(b"\x00")
        token = _register_vol_mrc_cache(str(mrc_dir), [str(vol)], (0,))
        with _VOL_CACHE_LOCK:
            _VOL_MRC_CACHE[token]["t0"] = 0.0
        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer.time.monotonic",
            lambda: 8000.0,
        )
        with pytest.raises(ValueError, match="expired"):
            volume_cell_gif_from_cache(token, 0, rows_expected=(0,))

    def test_cache_prunes_when_max_entries_exceeded(self, tmp_path: Path) -> None:
        from cryodrgn.dashboard.particle_explorer import _VOL_CACHE_MAX_ENTRIES

        tokens: list[str] = []
        for i in range(_VOL_CACHE_MAX_ENTRIES + 1):
            d = tmp_path / f"mrc_{i}"
            d.mkdir()
            p = d / "vol_001.mrc"
            p.write_bytes(b"\x00")
            tokens.append(_register_vol_mrc_cache(str(d), [str(p)], (i,)))
        assert len(_VOL_MRC_CACHE) <= _VOL_CACHE_MAX_ENTRIES
        assert tokens[0] not in _VOL_MRC_CACHE


class TestParticleExplorerVolumeGeneration:
    def test_generate_trajectory_volume_pngs_mocked(
        self,
        dashboard_experiment: DashboardExperiment,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        zdim = dashboard_experiment.z.shape[1]
        z_traj = dashboard_experiment.z[:3]

        def _fake_decode(
            exp: DashboardExperiment,
            z_values: np.ndarray,
            mrc_dir: str,
            *,
            progress_token: str | None = None,
        ) -> list[str]:
            os.makedirs(mrc_dir, exist_ok=True)
            paths: list[str] = []
            for i in range(len(z_values)):
                p = os.path.join(mrc_dir, f"vol_{i + 1:03d}.mrc")
                Path(p).write_bytes(b"\x00")
                paths.append(p)
            return paths

        def _fake_cycle(views, chimerax_cpus=1):
            from PIL import Image

            paths = []
            for view in views:
                Image.new("RGB", (4, 4)).save(view.out_png)
                paths.append(view.out_png)
            return (
                paths,
                "camera matrix" if views and views[0].report_view_matrix else None,
            )

        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer._decode_z_values_to_vol_paths",
            _fake_decode,
        )
        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer.render_landscape_cycle_static_views",
            _fake_cycle,
        )
        monkeypatch.setattr(
            "cryodrgn.dashboard.chimerax_animation.resolve_chimerax_volume_level",
            lambda _path, level=None: 0.5 if level is None else float(level),
        )
        pngs, token = generate_trajectory_volume_pngs(
            dashboard_experiment, z_traj, chimerax_cpus=2, pipeline=False
        )
        assert len(pngs) == 3
        assert token
        pngs2, _vm = generate_trajectory_volume_pngs(
            dashboard_experiment,
            z_traj,
            chimerax_cpus=4,
            view_turns=[("y", 15.0)],
            pipeline=False,
        )
        assert len(pngs2) == 3
        with pytest.raises(ValueError, match="must be"):
            generate_trajectory_volume_pngs(
                dashboard_experiment, np.zeros((2, zdim + 1))
            )

    def test_pipelined_decode_and_chimerax_pngs_completes(
        self,
        dashboard_experiment: DashboardExperiment,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        from cryodrgn.dashboard.particle_explorer import (
            _pipelined_decode_and_chimerax_pngs,
            volume_job_partial_register,
            volume_job_partial_snapshot,
            volume_job_progress_register_pipeline,
            volume_job_progress_snapshot,
        )

        z_traj = dashboard_experiment.z[:4]
        n = len(z_traj)
        token = "test-pipeline-token"

        def _fake_parallel(
            exp: DashboardExperiment,
            z_values: np.ndarray,
            mrc_dir: str,
            n_gpus: int,
        ) -> None:
            os.makedirs(mrc_dir, exist_ok=True)
            for i in range(len(z_values)):
                p = os.path.join(mrc_dir, f"vol_{i + 1:03d}.mrc")
                Path(p).write_bytes(b"\x00" * 1024)

        def _fake_render(mrc_path: str, out_png: str, **kwargs: object) -> str | None:
            from PIL import Image

            Image.new("RGB", (4, 4)).save(out_png)
            return "camera matrix" if kwargs.get("report_view_matrix") else None

        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer._decode_z_values_parallel_impl",
            _fake_parallel,
        )
        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer.cuda_gpu_count_for_decode",
            lambda: 2,
        )
        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer.render_static_png",
            _fake_render,
        )
        monkeypatch.setattr(
            "cryodrgn.dashboard.chimerax_animation.resolve_chimerax_volume_level",
            lambda _path, level=None: 0.5,
        )

        mrc_dir = tempfile.mkdtemp()
        png_dir = tempfile.mkdtemp()
        volume_job_progress_register_pipeline(token, n, n_gpus=2, n_cpus=4)
        volume_job_partial_register(token, n)

        pngs, vm, vol_files = _pipelined_decode_and_chimerax_pngs(
            dashboard_experiment,
            z_traj,
            mrc_dir,
            png_dir,
            chimerax_cpus=4,
            progress_token=token,
        )
        assert len(pngs) == n
        assert len(vol_files) == n
        snap = volume_job_progress_snapshot(token)
        assert snap is not None
        assert snap["decode_done"] == n
        assert snap["render_done"] == n
        partial = volume_job_partial_snapshot(token)
        assert partial is not None
        assert partial["complete"] is True
        assert len(partial["images"]) == n

    def test_rerender_chimerax_pngs_from_volume_cache(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        from cryodrgn.dashboard.particle_explorer import (
            _register_vol_mrc_cache,
            rerender_chimerax_pngs_from_volume_cache,
        )

        calls: list[int] = []
        mrc_a = tmp_path / "a.mrc"
        mrc_b = tmp_path / "b.mrc"
        mrc_a.write_bytes(b"\x00")
        mrc_b.write_bytes(b"\x00")

        def _fake_cycle(views, chimerax_cpus=1):
            calls.append(int(chimerax_cpus))
            from PIL import Image

            paths = []
            for view in views:
                Image.new("RGB", (4, 4)).save(view.out_png)
                paths.append(view.out_png)
            return paths, None

        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer.render_landscape_cycle_static_views",
            _fake_cycle,
        )
        monkeypatch.setattr(
            "cryodrgn.dashboard.chimerax_animation.resolve_chimerax_volume_level",
            lambda _path, level=None: 0.5 if level is None else float(level),
        )
        token = _register_vol_mrc_cache(str(tmp_path), [str(mrc_a), str(mrc_b)], ())
        blobs, _ = rerender_chimerax_pngs_from_volume_cache(
            token, chimerax_cpus=6, view_turns=[("x", 30.0)]
        )
        assert len(blobs) == 2
        assert calls == [6]

    def test_generate_montage_rejects_invalid_rows(
        self, dashboard_experiment: DashboardExperiment
    ) -> None:
        with pytest.raises(ValueError, match="Invalid plot row"):
            generate_montage_volume_pngs(dashboard_experiment, [-1])
        with pytest.raises(ValueError, match="Invalid plot row"):
            generate_montage_volume_pngs(dashboard_experiment, [])

    def test_decode_z_values_to_vol_paths_classic_branch(
        self,
        dashboard_experiment: DashboardExperiment,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        z = dashboard_experiment.z[:2]
        mrc_dir = tmp_path / "decode"
        called: list[str] = []

        def _fake_worker(
            exp: DashboardExperiment,
            z_values: np.ndarray,
            out_dir: str,
            *,
            device: int = 0,
            vol_start_index: int = 1,
        ) -> None:
            called.append(out_dir)
            os.makedirs(out_dir, exist_ok=True)
            for i in range(len(z_values)):
                Path(out_dir, f"vol_{vol_start_index + i:03d}.mrc").write_bytes(b"\x00")

        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer.cuda_gpu_count_for_decode",
            lambda: 1,
        )
        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer._decode_z_values_worker",
            _fake_worker,
        )
        paths = _decode_z_values_to_vol_paths(dashboard_experiment, z, str(mrc_dir))
        assert called == [str(mrc_dir)]
        assert len(paths) == 2

    def test_coerce_z_decode_matrix_single_row(self) -> None:
        z = np.arange(4, dtype=np.float64)
        out = _coerce_z_decode_matrix(z, 4)
        assert out.shape == (1, 4)

    def test_coerce_z_decode_matrix_two_rows(self) -> None:
        z = np.arange(8, dtype=np.float64).reshape(2, 4)
        out = _coerce_z_decode_matrix(z, 4)
        np.testing.assert_array_equal(out, z)

    def test_decode_classic_writes_unique_zfile(
        self,
        dashboard_experiment: DashboardExperiment,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        z = dashboard_experiment.z[:1]
        zfiles: list[str] = []

        def _fake_gen_volumes(weights, cfg, zfile, outdir, **kwargs):
            zfiles.append(zfile)
            os.makedirs(outdir, exist_ok=True)
            Path(outdir, "vol_001.mrc").write_bytes(b"\x00")

        monkeypatch.setattr(
            "cryodrgn.analysis.gen_volumes",
            _fake_gen_volumes,
        )
        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer._is_drgnai_config",
            lambda _cfg: False,
        )
        from cryodrgn.dashboard.particle_explorer import _decode_z_values_classic

        _decode_z_values_classic(
            dashboard_experiment, z, str(tmp_path), vol_start_index=3
        )
        assert len(zfiles) == 1
        assert zfiles[0].endswith("z_values_0003.txt")

    def test_sorted_vol_mrc_paths_skips_subdirs(self, tmp_path: Path) -> None:
        (tmp_path / "vol_001.mrc").write_bytes(b"\x00")
        (tmp_path / "vol_002.mrc").write_bytes(b"\x00")
        (tmp_path / "subdir").mkdir()
        out = _sorted_vol_mrc_paths(str(tmp_path), 2)
        assert len(out) == 2

    def test_generate_montage_cleanup_on_decode_failure(
        self,
        dashboard_experiment: DashboardExperiment,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        def _boom(*_a, **_k):
            raise RuntimeError("decode failed")

        monkeypatch.setattr(
            "cryodrgn.dashboard.particle_explorer._decode_z_values_to_vol_paths",
            _boom,
        )
        with pytest.raises(RuntimeError, match="decode failed"):
            generate_montage_volume_pngs(dashboard_experiment, [0])


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

    def test_save_selection_custom_path_roundtrip(
        self, flask_client, dashboard_experiment: DashboardExperiment, tmp_path
    ) -> None:
        """``sel_dir`` + custom basename writes both main and inverse pickles."""
        dest = tmp_path / "custom_out"
        dest.mkdir()
        rows = [0, 2, 4]
        r = flask_client.post(
            "/api/save_selection",
            json={
                "rows": rows,
                "basename": "region_a",
                "sel_dir": str(dest),
                "save_inverse": True,
            },
        )
        assert r.status_code == 200, r.get_json()
        js = r.get_json()
        main_path = dest / "region_a.pkl"
        inv_path = dest / "region_a_inverse.pkl"
        assert js["path"] == str(main_path)
        assert js["inverse_path"] == str(inv_path)
        assert main_path.is_file()
        assert inv_path.is_file()
        expected = np.asarray(dashboard_experiment.all_indices[rows], dtype=int)
        all_ds = np.asarray(dashboard_experiment.all_indices, dtype=int)
        with open(main_path, "rb") as fh:
            selected = pickle.load(fh)
        with open(inv_path, "rb") as fh:
            inverse = pickle.load(fh)
        assert np.array_equal(selected, expected)
        assert np.array_equal(inverse, np.setdiff1d(all_ds, expected))


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
        _monkeypatch_explorer_volumes_eligible(monkeypatch, eligible=False)
        r = flask_client.post(
            "/api/explorer_volume_media",
            json={"rows": [0], "mode": "static"},
        )
        assert r.status_code == 400

    def test_empty_rows_returns_400(self, flask_client, monkeypatch) -> None:
        _monkeypatch_explorer_volumes_eligible(monkeypatch, eligible=True)
        r = flask_client.post(
            "/api/explorer_volume_media",
            json={"rows": [], "mode": "static"},
        )
        assert r.status_code == 400

    def test_bad_mode_returns_400(self, flask_client, monkeypatch) -> None:
        _monkeypatch_explorer_volumes_eligible(monkeypatch, eligible=True)
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


class TestParticleExplorerBrowserSmoke:
    """Headless Chromium: cache/montage, colour legend, panels, and save modal."""

    pytestmark = pytest.mark.browser

    @pytest.mark.parametrize(
        "smoke_name,checks",
        [
            (
                "dashboard_smoke_particle_explorer",
                ("cached_images>=1", "grid_images>=1", "scatter_letters.count>=1"),
            ),
            (
                "dashboard_smoke_particle_explorer_color_covariate",
                ("discrete_toggles>=1", "points_after_color==scatter_points"),
            ),
            (
                "dashboard_smoke_particle_explorer_cache_expand",
                ("expanded_cached>initial_cached",),
            ),
        ],
        ids=["cache_montage", "color_covariate", "cache_expand"],
    )
    def test_plotly_smokes(
        self,
        smoke_name: str,
        checks: tuple[str, ...],
        playwright_page,
        dashboard_live_url,
    ) -> None:
        import tests.conftest as cf

        out = getattr(cf, smoke_name)(playwright_page, dashboard_live_url)
        for check in checks:
            if check == "points_after_color==scatter_points":
                assert out["points_after_color"] == out["scatter_points"]
            elif check == "expanded_cached>initial_cached":
                assert out["expanded_cached"] > out["initial_cached"]
            elif check.endswith(">=1"):
                key = check[: -len(">=1")]
                if "." in key:
                    cur = out
                    for part in key.split("."):
                        cur = cur[part]
                    assert cur >= 1
                else:
                    assert out[key] >= 1
            else:
                raise AssertionError(f"unknown check {check!r}")

    def test_selection_save_modal_opens(
        self, playwright_page, dashboard_live_url
    ) -> None:
        from tests.conftest import (
            DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
            _dashboard_smoke_wait_plot_ready,
        )

        base = dashboard_live_url.rstrip("/")
        playwright_page.goto(
            f"{base}/explorer",
            wait_until="domcontentloaded",
            timeout=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        _dashboard_smoke_wait_plot_ready(
            playwright_page, "scatter", timeout_ms=DASHBOARD_BROWSER_FAST_TIMEOUT_MS
        )
        playwright_page.evaluate(
            """() => {
              var fs = document.getElementById('particle-sel-fieldset');
              if (fs) fs.disabled = false;
            }"""
        )
        playwright_page.click("#save-indices-custom")
        playwright_page.wait_for_function(
            """() => {
              var panel = document.getElementById('sel-file-browser-panel');
              return panel && !panel.hidden;
            }""",
            timeout=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        assert playwright_page.evaluate(
            """() => document.body.classList.contains('cryo-explorer-save-modal-open')"""
        )

    def test_volume_explorer_panels_attached(
        self, playwright_page, dashboard_volumes_eligible_live_url
    ) -> None:
        from tests.conftest import dashboard_smoke_particle_explorer_panels

        out = dashboard_smoke_particle_explorer_panels(
            playwright_page, dashboard_volumes_eligible_live_url
        )
        assert out["panels"] >= 12
        assert out["volumes_panel"]

    def test_volumes_panel_omitted_when_ineligible(
        self, playwright_page, dashboard_volumes_ineligible_live_url
    ) -> None:
        from tests.conftest import dashboard_smoke_particle_explorer_no_volumes_panel

        out = dashboard_smoke_particle_explorer_no_volumes_panel(
            playwright_page, dashboard_volumes_ineligible_live_url
        )
        assert out["volumes_panel"] is False
