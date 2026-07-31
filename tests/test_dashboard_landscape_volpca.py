"""Tests for the volume sketched landscape (vol PCA) dashboard UI and APIs.

Uses the same session-scoped ``dashboard_workdir`` fixture as
``test_dashboard_core`` / ``conftest`` (epoch **2** analyze outputs).

Parallel runs: safe under ``pytest -n2 --dist=loadscope`` (and similar). The
shared training output is only read while copying; the writable tree lives
under ``tmp_path_factory`` (per-worker basetemp). ``conftest`` already
serializes building ``dashboard_workdir`` with a file lock across xdist
workers. Landscape-specific fixtures are **module**-scoped so loadscope can
keep dependent tests on one worker and we only **copytree** once per module
per worker; HTTP clients stay **function**-scoped so Flask sessions do not
leak between tests.

``tmp_path_factory.mktemp()`` creates an **empty** directory; ``shutil.copytree``
requires the destination leaf not to exist, so copies always land in
``<mktemp>/cryo_out``. Optional ``landscape.*`` trees are stripped from copies
when tests assume analyze_landscape was never run (golden fixtures may include
them).
"""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import numpy as np
import pytest

from cryodrgn import utils
from cryodrgn.dashboard.data import DashboardExperiment, load_experiment
from cryodrgn.dashboard.landscape_volpca import (
    _cycle_gif_from_png_paths,
    _gif_first_frame_to_png,
    _landscape_cycle_gif_timing,
    _lookup_rotate_keyframe_pngs,
    _sketch_continuous_covariate_per_volume_values,
    animation_payload_b64,
    clear_landscape_animation_cache,
    generate_landscape_volume_animations,
    kmeans_sorted_vol_indices,
    landscape_analysis_ready,
    landscape_color_options,
    landscape_dir_for_epoch,
    landscape_vol_state_hex_by_vol_index,
    landscape_volpca_scatter_figure,
    landscape_volpca_scatter_json,
    list_landscape_epochs,
    load_landscape_sketch_umap_coords,
    load_pca_explained_variance,
    load_sketch_centroid_plot_df_rows,
    load_sketch_state_labels,
    load_vol_pca_matrix,
    meta_for_api,
    normalize_landscape_view_rotations,
    parse_landscape_axis_spec,
    resolve_kmeans_sketch_bundle,
    save_landscape_animations,
    sketch_plot_color_covariate_variable_label,
    sketch_vol_color_covariate_overlay_text,
    sketch_vol_marker_hex_by_vol_index,
    sketch_vol_rotate_frames_covariate_overlay,
    vol_mrc_path,
)
from tests.conftest import (
    DASHBOARD_ANALYZE_EPOCH as _ANALYZE_EPOCH,
    decode_plotly_figure,
    plotly_trace_array,
)

pytestmark = pytest.mark.dashboard

_LANDSCAPE_K = 3


class _FakeVolPcaObj:
    """Stand-in for sklearn PCA with ``explained_variance_ratio_`` (pickle-safe)."""

    explained_variance_ratio_ = np.array([0.453, 0.12], dtype=np.float64)


def _write_landscape_bundle(
    root: Path,
    *,
    epoch: int = _ANALYZE_EPOCH,
    k: int = _LANDSCAPE_K,
    with_states: bool = False,
    with_evr: bool = False,
    bad_centers: bool = False,
) -> tuple[str, str]:
    """Minimal ``landscape.{epoch}`` tree under ``root`` for unit tests."""
    land = root / f"landscape.{epoch}"
    if land.is_dir():
        shutil.rmtree(land)
    km = land / f"kmeans{k}"
    km.mkdir(parents=True)
    pc = np.array(
        [[0.0, 0.1], [1.0, -0.5], [0.2, 0.3]],
        dtype=np.float64,
    )
    utils.save_pkl(pc, land / f"vol_pca_{k}.pkl")
    for i in range(1, k + 1):
        (km / f"vol_{i:03d}.mrc").write_bytes(b"\x00")
    centers_path = km / "centers_ind.txt"
    if bad_centers:
        centers_path.write_text("0\n1\n", encoding="utf-8")
    else:
        with centers_path.open("w", encoding="utf-8") as fh:
            for row in range(k):
                fh.write(f"{row}\n")
    umap_full = np.array(
        [[0.0, 1.0], [2.0, 3.0], [4.0, 5.0], [6.0, 7.0], [8.0, 9.0]],
        dtype=np.float64,
    )
    utils.save_pkl(umap_full, land / "umap.pkl")
    if with_states:
        sk = land / "sketch_clustering_0.5"
        sk.mkdir()
        utils.save_pkl(np.array([0, 1, 2], dtype=np.int64), sk / "state_labels.pkl")
    if with_evr:
        utils.save_pkl(_FakeVolPcaObj(), land / "vol_pca_obj.pkl")
    return str(land), str(km)


@pytest.fixture(autouse=True)
def _reset_landscape_animation_cache() -> None:
    yield
    clear_landscape_animation_cache()


@pytest.fixture
def landscape_bundle_tmp(
    tmp_path: Path, dashboard_workdir: str
) -> tuple[DashboardExperiment, str, str]:
    """Dashboard experiment plus a writable landscape tree on disk."""
    work = tmp_path / "cryo_out"
    shutil.copytree(dashboard_workdir, work)
    for child in work.iterdir():
        if child.is_dir() and child.name.startswith("landscape."):
            shutil.rmtree(child)
    land, km = _write_landscape_bundle(work)
    exp = load_experiment(str(work))
    return exp, land, km


@pytest.fixture(scope="module")
def experiment_landscape(
    dashboard_workdir_with_landscape_volpca: str,
) -> DashboardExperiment:
    return load_experiment(dashboard_workdir_with_landscape_volpca)


@pytest.fixture(scope="module")
def experiment_plain_no_landscape(
    dashboard_workdir_plain_copy: str,
) -> DashboardExperiment:
    return load_experiment(dashboard_workdir_plain_copy)


class TestLandscapeVolpcaKmeansHelpers:
    """``vol_mean.mrc`` must not become a fake vol 000 index."""

    def test_kmeans_sorted_vol_indices_excludes_vol_mean(self, tmp_path) -> None:
        d = tmp_path / "kmeans100"
        d.mkdir()
        (d / "vol_mean.mrc").write_bytes(b"x")
        (d / "vol_001.mrc").write_bytes(b"x")
        (d / "vol_100.mrc").write_bytes(b"x")
        assert kmeans_sorted_vol_indices(str(d)) == [1, 100]


class TestLandscapeVolpcaFilesystemHelpers:
    def test_list_landscape_epochs_empty(self, tmp_path) -> None:
        assert list_landscape_epochs(str(tmp_path)) == []

    def test_list_landscape_epochs_finds_folder(
        self, dashboard_workdir_with_landscape_volpca: str
    ) -> None:
        assert list_landscape_epochs(dashboard_workdir_with_landscape_volpca) == [
            _ANALYZE_EPOCH
        ]

    def test_landscape_analysis_ready(
        self, dashboard_workdir_with_landscape_volpca: str
    ) -> None:
        assert landscape_analysis_ready(
            dashboard_workdir_with_landscape_volpca,
            _ANALYZE_EPOCH,
        )

    def test_landscape_dir_for_epoch(
        self, dashboard_workdir_with_landscape_volpca: str
    ) -> None:
        d = landscape_dir_for_epoch(
            dashboard_workdir_with_landscape_volpca, _ANALYZE_EPOCH
        )
        assert d.endswith(f"landscape.{_ANALYZE_EPOCH}")
        assert os.path.isdir(d)


class TestLandscapeVolpcaMetaAndScatter:
    def test_meta_ok(self, experiment_landscape: DashboardExperiment) -> None:
        m = meta_for_api(experiment_landscape)
        assert m["ok"] is True
        assert m["landscape_epoch"] == _ANALYZE_EPOCH
        assert m["n_volumes"] == _LANDSCAPE_K
        assert m["n_pc"] == 2
        assert m["n_umap"] == 2
        assert m["has_state_color"] is False
        assert m["kmeans_k"] == _LANDSCAPE_K
        opts = m["color_options"]
        assert any(o["value"] == "none" for o in opts)
        assert any(o["value"] == "state" for o in opts)

    def test_meta_no_landscape_dirs(
        self, experiment_plain_no_landscape: DashboardExperiment
    ) -> None:
        m = meta_for_api(experiment_plain_no_landscape)
        assert m["ok"] is False
        assert "No landscape" in (m.get("error") or "")

    def test_scatter_json_marker_size_matches_sketch_constant(
        self, experiment_landscape: DashboardExperiment
    ) -> None:
        """Regression: vol-PCA scatter JSON round-trips and markers match ``_VOLSKETCH_SCATTER_MARKER``."""
        from cryodrgn.dashboard import landscape_volpca as lv

        land = landscape_dir_for_epoch(
            experiment_landscape.workdir,
            _ANALYZE_EPOCH,
        )
        js = landscape_volpca_scatter_json(
            land,
            experiment_landscape,
            axis_x="pc:0",
            axis_y="pc:1",
            color_mode="none",
            continuous_palette=None,
        )
        fig = json.loads(js)
        trace = fig["data"][0]
        assert trace["type"] == "scattergl"
        assert "ids" in trace
        assert len(trace["ids"]) == _LANDSCAPE_K
        sz = trace["marker"]["size"]
        assert float(sz) == pytest.approx(
            float(lv._VOLSKETCH_SCATTER_MARKER), rel=0, abs=1e-9
        )
        assert float(sz) == pytest.approx(9.0 * (1.0 - 0.13) * 1.3, rel=0, abs=1e-9)

    def test_scatter_rejects_bad_axes(
        self, experiment_landscape: DashboardExperiment
    ) -> None:
        land = landscape_dir_for_epoch(
            experiment_landscape.workdir,
            _ANALYZE_EPOCH,
        )
        with pytest.raises(ValueError, match="distinct"):
            landscape_volpca_scatter_json(
                land,
                experiment_landscape,
                axis_x="pc:0",
                axis_y="pc:0",
                color_mode="none",
            )

    def test_scatter_umap_axes(self, experiment_landscape: DashboardExperiment) -> None:
        land = landscape_dir_for_epoch(
            experiment_landscape.workdir,
            _ANALYZE_EPOCH,
        )
        js = landscape_volpca_scatter_json(
            land,
            experiment_landscape,
            axis_x="umap:0",
            axis_y="umap:1",
            color_mode="none",
        )
        fig = decode_plotly_figure(json.loads(js))
        trace = fig["data"][0]
        assert plotly_trace_array(trace, "x") == [0.0, 2.0, 4.0]
        assert plotly_trace_array(trace, "y") == [1.0, 3.0, 5.0]


class TestSketchVolContinuousColors:
    """ChimeraX cycle colours must differ per sketch volume for continuous covariates."""

    def test_landscape_vol_pc_values_differ_per_volume(
        self, experiment_landscape: DashboardExperiment
    ) -> None:
        land = landscape_dir_for_epoch(
            experiment_landscape.workdir,
            _ANALYZE_EPOCH,
        )
        cont = _sketch_continuous_covariate_per_volume_values(
            land,
            experiment_landscape,
            "landscape_vol_PC1",
            [1, 2, 3],
        )
        assert cont is not None
        vals, cmin, cmax = cont
        assert len(vals) == 3
        assert len({float(v) for v in vals if np.isfinite(v)}) == 3
        assert cmax > cmin

    def test_sketch_vol_marker_hex_distinct_for_vol_pc(
        self, experiment_landscape: DashboardExperiment
    ) -> None:
        land = landscape_dir_for_epoch(
            experiment_landscape.workdir,
            _ANALYZE_EPOCH,
        )
        hexes = sketch_vol_marker_hex_by_vol_index(
            land,
            experiment_landscape,
            color_mode="landscape_vol_PC1",
        )
        assert len(hexes) == _LANDSCAPE_K
        assert len(set(hexes.values())) == _LANDSCAPE_K

    def test_landscape_vol_pc_column_name_case_insensitive(
        self, experiment_landscape: DashboardExperiment
    ) -> None:
        land = landscape_dir_for_epoch(
            experiment_landscape.workdir,
            _ANALYZE_EPOCH,
        )
        lower = sketch_vol_marker_hex_by_vol_index(
            land,
            experiment_landscape,
            color_mode="landscape_vol_pc1",
        )
        upper = sketch_vol_marker_hex_by_vol_index(
            land,
            experiment_landscape,
            color_mode="landscape_vol_PC1",
        )
        assert lower == upper


class TestLandscapeVolpcaPureHelpers:
    def test_normalize_landscape_view_rotations_empty(self) -> None:
        assert normalize_landscape_view_rotations(None) == []
        assert normalize_landscape_view_rotations({}) == []

    def test_normalize_landscape_view_rotations_axes(self) -> None:
        out = normalize_landscape_view_rotations({"x": 90.0, "y": 0, "z": -45})
        assert out == [("x", 90.0), ("z", -45.0)]

    def test_normalize_landscape_view_rotations_rejects_bad_type(self) -> None:
        with pytest.raises(ValueError, match="object"):
            normalize_landscape_view_rotations([1, 2, 3])

    def test_normalize_landscape_view_rotations_rejects_non_finite(self) -> None:
        with pytest.raises(ValueError, match="finite"):
            normalize_landscape_view_rotations({"x": float("nan")})

    def test_list_landscape_epochs_missing_workdir(self, tmp_path: Path) -> None:
        assert list_landscape_epochs(str(tmp_path / "nope")) == []

    def test_resolve_kmeans_sketch_bundle_missing(self, tmp_path: Path) -> None:
        d = tmp_path / "landscape.1"
        d.mkdir()
        with pytest.raises(FileNotFoundError):
            resolve_kmeans_sketch_bundle(str(d))

    def test_kmeans_sorted_vol_indices_empty_dir(self, tmp_path: Path) -> None:
        d = tmp_path / "kmeans"
        d.mkdir()
        assert kmeans_sorted_vol_indices(str(d)) == []

    def test_load_vol_pca_matrix_bad_ndim(self, tmp_path: Path) -> None:
        land, _ = _write_landscape_bundle(tmp_path)
        utils.save_pkl(np.array([1.0, 2.0]), Path(land) / "vol_pca_3.pkl")
        with pytest.raises(ValueError, match="2-D"):
            load_vol_pca_matrix(land, 3)

    def test_load_pca_explained_variance_missing(self, tmp_path: Path) -> None:
        land, _ = _write_landscape_bundle(tmp_path)
        assert load_pca_explained_variance(land) is None

    def test_load_pca_explained_variance_present(self, tmp_path: Path) -> None:
        land, _ = _write_landscape_bundle(tmp_path, with_evr=True)
        evr = load_pca_explained_variance(land)
        assert evr is not None
        assert len(evr) == 2

    def test_load_sketch_state_labels_missing(self, tmp_path: Path) -> None:
        land, _ = _write_landscape_bundle(tmp_path)
        assert load_sketch_state_labels(land) is None

    def test_load_sketch_state_labels_present(self, tmp_path: Path) -> None:
        land, _ = _write_landscape_bundle(tmp_path, with_states=True)
        states = load_sketch_state_labels(land)
        assert states is not None
        assert list(states) == [0, 1, 2]

    def test_landscape_vol_state_hex_without_states(self, tmp_path: Path) -> None:
        land, km = _write_landscape_bundle(tmp_path)
        assert landscape_vol_state_hex_by_vol_index(land, km, 3) is None

    def test_landscape_vol_state_hex_with_states(self, tmp_path: Path) -> None:
        land, km = _write_landscape_bundle(tmp_path, with_states=True)
        m = landscape_vol_state_hex_by_vol_index(land, km, 3)
        assert m is not None
        assert len(m) == 3

    def test_vol_mrc_path_missing(self, tmp_path: Path) -> None:
        land, km = _write_landscape_bundle(tmp_path)
        with pytest.raises(FileNotFoundError):
            vol_mrc_path(km, 99)

    def test_sketch_plot_color_covariate_variable_label(self) -> None:
        assert sketch_plot_color_covariate_variable_label("none") is None
        assert (
            sketch_plot_color_covariate_variable_label("state") == "Agglomerative state"
        )
        assert sketch_plot_color_covariate_variable_label("PC1") == "latent-space PC1"

    def test_landscape_color_options_includes_numeric(
        self, experiment_landscape: DashboardExperiment
    ) -> None:
        opts = landscape_color_options(experiment_landscape)
        values = {o["value"] for o in opts}
        assert "none" in values
        assert "state" in values
        assert any(v not in ("none", "state") for v in values)
        labels_option = next(o for o in opts if o["value"] == "labels")
        assert labels_option["label"] == "k-means labels"

    def test_load_sketch_centroid_plot_df_rows_errors(self, tmp_path: Path) -> None:
        land, km = _write_landscape_bundle(tmp_path, bad_centers=True)
        with pytest.raises(ValueError, match="centers_ind"):
            load_sketch_centroid_plot_df_rows(km, 3)

    def test_load_landscape_sketch_umap_coords_none_without_umap(
        self, tmp_path: Path
    ) -> None:
        land, km = _write_landscape_bundle(tmp_path)
        os.remove(os.path.join(land, "umap.pkl"))
        assert load_landscape_sketch_umap_coords(land, km, 3) is None

    def test_parse_landscape_axis_spec(self) -> None:
        assert parse_landscape_axis_spec("pc:0") == ("pc", 0)
        assert parse_landscape_axis_spec("UMAP:1") == ("umap", 1)
        with pytest.raises(ValueError, match="Invalid axis"):
            parse_landscape_axis_spec("bogus")
        with pytest.raises(ValueError, match="kind"):
            parse_landscape_axis_spec("foo:0")

    def test_landscape_cycle_gif_timing(self) -> None:
        fpv, ms = _landscape_cycle_gif_timing(4)
        assert fpv == 4
        assert ms >= 30
        fpv_hi, _ = _landscape_cycle_gif_timing(100)
        assert fpv_hi == 30


class TestLandscapeVolpcaScatterColorModes:
    def test_scatter_pc_axes_with_explained_variance_title(
        self, landscape_bundle_tmp: tuple[DashboardExperiment, str, str]
    ) -> None:
        exp, land, _km = landscape_bundle_tmp
        utils.save_pkl(_FakeVolPcaObj(), os.path.join(land, "vol_pca_obj.pkl"))
        fig = landscape_volpca_scatter_figure(
            land,
            exp,
            axis_x=("pc", 0),
            axis_y=("pc", 1),
            color_mode="none",
        )
        assert "45.3%" in fig.layout.xaxis.title.text

    def test_scatter_labels_color(
        self, experiment_landscape: DashboardExperiment
    ) -> None:
        land = landscape_dir_for_epoch(experiment_landscape.workdir, _ANALYZE_EPOCH)
        js = landscape_volpca_scatter_json(
            land,
            experiment_landscape,
            axis_x="pc:0",
            axis_y="pc:1",
            color_mode="labels",
        )
        trace = json.loads(js)["data"][0]
        assert trace["marker"]["color"]

    def test_scatter_unknown_color_falls_back_to_gray(
        self, experiment_landscape: DashboardExperiment
    ) -> None:
        land = landscape_dir_for_epoch(experiment_landscape.workdir, _ANALYZE_EPOCH)
        js = landscape_volpca_scatter_json(
            land,
            experiment_landscape,
            axis_x="pc:0",
            axis_y="pc:1",
            color_mode="not_a_real_column",
        )
        trace = json.loads(js)["data"][0]
        assert trace["marker"]["color"] == "#4a5568"

    def test_scatter_state_color_with_states(
        self, landscape_bundle_tmp: tuple[DashboardExperiment, str, str]
    ) -> None:
        exp, land, _km = landscape_bundle_tmp
        sk = Path(land) / "sketch_clustering_0.5"
        sk.mkdir()
        utils.save_pkl(np.array([0, 1, 2], dtype=np.int64), sk / "state_labels.pkl")
        js = landscape_volpca_scatter_json(
            land,
            exp,
            axis_x="pc:0",
            axis_y="pc:1",
            color_mode="state",
        )
        trace = json.loads(js)["data"][0]
        assert len(trace["marker"]["color"]) == _LANDSCAPE_K

    def test_scatter_umap_axis_out_of_range(
        self, experiment_landscape: DashboardExperiment
    ) -> None:
        land = landscape_dir_for_epoch(experiment_landscape.workdir, _ANALYZE_EPOCH)
        with pytest.raises(ValueError, match="UMAP axis"):
            landscape_volpca_scatter_figure(
                land,
                experiment_landscape,
                axis_x=("umap", 99),
                axis_y=("umap", 0),
                color_mode="none",
            )


class TestSketchVolOverlayHelpers:
    def test_sketch_vol_marker_hex_labels(
        self, experiment_landscape: DashboardExperiment
    ) -> None:
        land = landscape_dir_for_epoch(experiment_landscape.workdir, _ANALYZE_EPOCH)
        hexes = sketch_vol_marker_hex_by_vol_index(
            land,
            experiment_landscape,
            color_mode="labels",
        )
        assert len(hexes) == _LANDSCAPE_K

    def test_sketch_vol_marker_hex_none_returns_empty(
        self, experiment_landscape: DashboardExperiment
    ) -> None:
        land = landscape_dir_for_epoch(experiment_landscape.workdir, _ANALYZE_EPOCH)
        assert (
            sketch_vol_marker_hex_by_vol_index(
                land, experiment_landscape, color_mode="none"
            )
            == {}
        )

    def test_sketch_vol_color_overlay_text_pc(
        self, experiment_landscape: DashboardExperiment
    ) -> None:
        land = landscape_dir_for_epoch(experiment_landscape.workdir, _ANALYZE_EPOCH)
        _k, km = resolve_kmeans_sketch_bundle(land)
        txt = sketch_vol_color_covariate_overlay_text(
            land,
            km,
            _LANDSCAPE_K,
            experiment_landscape,
            "landscape_vol_PC1",
            1,
        )
        assert txt is not None
        assert txt != "—"

    def test_sketch_vol_color_overlay_unknown_vol(self, experiment_landscape) -> None:
        land = landscape_dir_for_epoch(experiment_landscape.workdir, _ANALYZE_EPOCH)
        _k, km = resolve_kmeans_sketch_bundle(land)
        assert (
            sketch_vol_color_covariate_overlay_text(
                land, km, _LANDSCAPE_K, experiment_landscape, "none", 999
            )
            is None
        )

    def test_sketch_vol_rotate_frames_overlay(
        self, experiment_landscape: DashboardExperiment
    ) -> None:
        from dataclasses import replace

        from cryodrgn.dashboard.landscape_full_3d import (
            VOL_LANDSCAPE_NEAREST_SKETCH_VOL,
        )

        land = landscape_dir_for_epoch(experiment_landscape.workdir, _ANALYZE_EPOCH)
        _k, km = resolve_kmeans_sketch_bundle(land)
        df = experiment_landscape.plot_df.copy()
        df[VOL_LANDSCAPE_NEAREST_SKETCH_VOL] = 1
        exp = replace(experiment_landscape, plot_df=df)
        rot = sketch_vol_rotate_frames_covariate_overlay(
            exp,
            land,
            km,
            _LANDSCAPE_K,
            "UMAP1",
            1,
            8,
            "Viridis",
        )
        assert rot is not None
        assert rot["style"] == "rotate_frames"
        assert len(rot["frame_covariate_texts"]) == 8

    def test_continuous_covariate_from_plot_df(
        self, experiment_landscape: DashboardExperiment
    ) -> None:
        land = landscape_dir_for_epoch(experiment_landscape.workdir, _ANALYZE_EPOCH)
        col = next(
            c
            for c in experiment_landscape.numeric_columns
            if c.lower() not in ("labels",)
        )
        cont = _sketch_continuous_covariate_per_volume_values(
            land, experiment_landscape, col, [1, 2, 3]
        )
        assert cont is not None


class TestLandscapeVolpcaAnimationPipeline:
    @staticmethod
    def _fake_batch_rotate(specs, gif_frames, chimerax_cpus):
        from PIL import Image

        results = []
        for spec in specs:
            img = Image.new("RGB", (12, 12), color=(200, 10, 10))
            img.save(spec["out_gif"], save_all=True, duration=50, loop=0)
            vm = "1,0,0,0,0,1,0,0,0,0,1,0" if spec.get("report_view_matrix") else None
            results.append((spec["vol_key"], vm))
        return results

    @staticmethod
    def _fake_cycle_render(views, chimerax_cpus):
        from PIL import Image

        paths = []
        vm = None
        for v in views:
            img = Image.new("RGB", (12, 12), color=(10, 200, 10))
            img.save(v.out_png)
            paths.append(v.out_png)
            if v.report_view_matrix:
                vm = "matrix snapshot"
        return paths, vm

    def test_generate_rotate_each_animations(
        self,
        landscape_bundle_tmp: tuple[DashboardExperiment, str, str],
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        exp, land, _km = landscape_bundle_tmp
        monkeypatch.setattr(
            "cryodrgn.dashboard.landscape_volpca.batch_landscape_rotate_gifs",
            self._fake_batch_rotate,
        )
        token, files, rendered = generate_landscape_volume_animations(
            land,
            [1, 2],
            exp=exp,
            mode="rotate_each",
            gif_frames=6,
            plot_color_mode="landscape_vol_PC1",
        )
        assert rendered == [1, 2]
        assert len(files) == 2
        assert all(f["kind"] == "rotate" for f in files)
        payload = animation_payload_b64(token)
        assert len(payload) == 2
        assert all("gif_b64" in p for p in payload)

    def test_generate_cycle_animations_multi_volume(
        self,
        landscape_bundle_tmp: tuple[DashboardExperiment, str, str],
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        exp, land, _km = landscape_bundle_tmp
        monkeypatch.setattr(
            "cryodrgn.dashboard.landscape_volpca.render_landscape_cycle_static_views",
            self._fake_cycle_render,
        )
        token, files, rendered = generate_landscape_volume_animations(
            land,
            [1, 2, 3],
            exp=exp,
            mode="cycle",
            plot_color_mode="landscape_vol_PC2",
        )
        assert rendered == [1, 2, 3]
        assert len(files) == 1
        assert files[0]["kind"] == "cycle"
        assert files[0]["vol"] is None
        saved = save_landscape_animations(
            token, None, exp=exp, landscape_epoch=_ANALYZE_EPOCH
        )
        assert any(p.endswith("cycle.gif") for p in saved)

    def test_generate_cycle_single_volume(
        self,
        landscape_bundle_tmp: tuple[DashboardExperiment, str, str],
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        exp, land, _km = landscape_bundle_tmp
        monkeypatch.setattr(
            "cryodrgn.dashboard.landscape_volpca.render_landscape_cycle_static_views",
            self._fake_cycle_render,
        )
        _token, files, _ = generate_landscape_volume_animations(
            land,
            [2],
            exp=exp,
            mode="cycle",
        )
        assert len(files) == 1
        assert files[0]["vol"] == 2

    def test_generate_rejects_mode_both(
        self, landscape_bundle_tmp: tuple[DashboardExperiment, str, str]
    ) -> None:
        exp, land, _km = landscape_bundle_tmp
        with pytest.raises(ValueError, match="both"):
            generate_landscape_volume_animations(land, [1], exp=exp, mode="both")

    def test_lookup_rotate_keyframes_reuse(
        self,
        landscape_bundle_tmp: tuple[DashboardExperiment, str, str],
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ) -> None:
        exp, land, _km = landscape_bundle_tmp
        monkeypatch.setattr(
            "cryodrgn.dashboard.landscape_volpca.batch_landscape_rotate_gifs",
            self._fake_batch_rotate,
        )
        token1, _, _ = generate_landscape_volume_animations(
            land,
            [1],
            exp=exp,
            mode="rotate_each",
            plot_color_mode="none",
        )
        token2, files2, _ = generate_landscape_volume_animations(
            land,
            [1],
            exp=exp,
            mode="cycle",
            reuse_rotate_keyframes_token=token1,
            plot_color_mode="none",
        )
        assert len(files2) == 1
        reuse = _lookup_rotate_keyframe_pngs(
            token1, land, "none", "none", None, [], None
        )
        assert reuse is not None
        assert reuse[0]

    def test_cycle_gif_from_png_paths(self, tmp_path: Path) -> None:
        from PIL import Image

        pngs = []
        for i in range(2):
            p = tmp_path / f"f{i}.png"
            Image.new("RGB", (8, 8), color=(i * 80, 0, 0)).save(p)
            pngs.append(str(p))
        out = tmp_path / "out.gif"
        _cycle_gif_from_png_paths(pngs, str(out), cycle_frames_per_vol=3)
        assert out.is_file()
        _gif_first_frame_to_png(str(out), str(tmp_path / "frame0.png"))
        assert (tmp_path / "frame0.png").is_file()

    def test_meta_wrong_epoch(self, dashboard_workdir_plain_copy: str) -> None:
        exp = load_experiment(dashboard_workdir_plain_copy)
        m = meta_for_api(exp)
        assert m["ok"] is False

    def test_animation_payload_unknown_token(self) -> None:
        with pytest.raises(ValueError, match="Unknown or expired"):
            animation_payload_b64("not-a-real-token")

    def test_save_animations_unknown_token(
        self, experiment_landscape: DashboardExperiment
    ) -> None:
        with pytest.raises(ValueError, match="Unknown or expired"):
            save_landscape_animations(
                "nope",
                None,
                exp=experiment_landscape,
                landscape_epoch=_ANALYZE_EPOCH,
            )

    def test_generate_empty_volume_list_raises(
        self, landscape_bundle_tmp: tuple[DashboardExperiment, str, str]
    ) -> None:
        exp, land, _km = landscape_bundle_tmp
        with pytest.raises(ValueError, match="at least one"):
            generate_landscape_volume_animations(land, [], exp=exp, mode="cycle")


class TestLandscapeVolpcaFlaskRoutes:
    def test_meta_route(self, flask_client_landscape) -> None:
        r = flask_client_landscape.get("/api/landscape_volpca/meta")
        assert r.status_code == 200
        m = r.get_json()
        assert m["ok"] is True
        assert m["n_volumes"] == _LANDSCAPE_K

    def test_scatter_route(self, flask_client_landscape) -> None:
        r = flask_client_landscape.get(
            "/api/landscape_volpca/scatter",
            query_string={"axis_x": "pc:0", "axis_y": "pc:1", "color": "none"},
        )
        assert r.status_code == 200
        fig = json.loads(r.get_data(as_text=True))
        assert fig["data"][0]["type"] == "scattergl"

    def test_scatter_route_legacy_pc_params(self, flask_client_landscape) -> None:
        r = flask_client_landscape.get(
            "/api/landscape_volpca/scatter",
            query_string={"pc_x": 0, "pc_y": 1, "color": "none"},
        )
        assert r.status_code == 200
        fig = json.loads(r.get_data(as_text=True))
        assert fig["data"][0]["type"] == "scattergl"

    def test_scatter_route_bad_color(self, flask_client_landscape) -> None:
        r = flask_client_landscape.get(
            "/api/landscape_volpca/scatter",
            query_string={"pc_x": 0, "pc_y": 1, "color": "__not_a_column__"},
        )
        assert r.status_code == 400
        assert "error" in (r.get_json() or {})

    def test_page_renders_with_landscape(self, flask_client_landscape) -> None:
        r = flask_client_landscape.get("/landscape-volpca")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        # ``volsketch-grid`` also appears in always-loaded CSS; require the live grid div.
        assert '<div class="volsketch-grid">' in body
        assert 'id="volsketch"' in body

    def test_page_no_landscape_message(self, flask_client_no_landscape) -> None:
        r = flask_client_no_landscape.get("/landscape-volpca")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "analyze_landscape" in body
        assert "No landscape.N folders" in body
        assert '<div class="volsketch-grid">' not in body
        assert 'id="volsketch"' not in body

    # ``state`` colouring is covered end to end by
    # ``TestLandscapeVolpcaWorkflow.test_state_colour_is_offered_and_refused_by_the_same_flag``,
    # which checks both branches of ``has_state_color`` rather than just the refusal.
    @pytest.mark.parametrize(
        "method,path,payload,expected_error",
        [
            (
                "GET",
                "/api/landscape_volpca/scatter",
                {"axis_x": "pc:0", "color": "none"},
                "axis_y",
            ),
            (
                "POST",
                "/api/landscape_volpca/generate_animations",
                {},
                None,
            ),
            (
                "POST",
                "/api/landscape_volpca/generate_animations",
                {"vol_indices": [1], "view_rotations": [1, 2, 3]},
                None,
            ),
            (
                "POST",
                "/api/landscape_volpca/generate_animations",
                {"vol_indices": [1], "view_matrix": "camera 1,2"},
                "12 numbers",
            ),
            ("POST", "/api/landscape_volpca/save_animations", {}, None),
            (
                "POST",
                "/api/landscape_volpca/save_animations",
                {"token": "nope", "out_dir": 123},
                None,
            ),
        ],
        ids=[
            "scatter_missing_axis_y",
            "animations_without_volumes",
            "animations_bad_view_rotations",
            "animations_bad_view_matrix",
            "save_without_token",
            "save_bad_out_dir",
        ],
    )
    def test_malformed_requests_are_refused(
        self, flask_client_landscape, method, path, payload, expected_error
    ) -> None:
        if method == "GET":
            r = flask_client_landscape.get(path, query_string=payload)
        else:
            r = flask_client_landscape.post(path, json=payload)
        assert r.status_code == 400, r.get_data(as_text=True)[:200]
        if expected_error:
            assert expected_error in (r.get_json() or {}).get("error", "")

    def test_generate_animations_cycle_happy_path(
        self,
        flask_client_landscape,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        def _fake_cycle(views, chimerax_cpus):
            from PIL import Image

            paths = []
            for v in views:
                Image.new("RGB", (8, 8)).save(v.out_png)
                paths.append(v.out_png)
            return paths, None

        monkeypatch.setattr(
            "cryodrgn.dashboard.landscape_volpca.render_landscape_cycle_static_views",
            _fake_cycle,
        )
        r = flask_client_landscape.post(
            "/api/landscape_volpca/generate_animations",
            json={"vol_indices": [1, 2], "mode": "cycle", "plot_color_mode": "none"},
        )
        assert r.status_code == 200, r.get_json()
        js = r.get_json()
        assert js.get("token")
        assert js.get("items")


@pytest.fixture(scope="module")
def dashboard_landscape_live_url(dashboard_landscape_volpca_live_url: str):
    """Alias kept for this module's browser smoke tests."""
    return dashboard_landscape_volpca_live_url


class TestLandscapeVolpcaPlotlyBrowserSmoke:
    """Headless Chromium: vol PCA overlay, clear, and axis reload."""

    pytestmark = pytest.mark.browser

    @pytest.mark.parametrize(
        "smoke_name,key,expect",
        [
            ("dashboard_smoke_landscape_volpca", "overlayTexts", 1),
            (
                "dashboard_smoke_landscape_volpca_clear_selection",
                "overlay_cleared",
                True,
            ),
            ("dashboard_smoke_landscape_volpca_axis_reload", "points", 1),
        ],
        ids=["overlay_letters", "clear_selection", "axis_reload"],
    )
    def test_volpca_browser_smokes(
        self,
        smoke_name: str,
        key: str,
        expect,
        playwright_page,
        dashboard_landscape_live_url,
    ) -> None:
        import tests.conftest as cf

        out = getattr(cf, smoke_name)(playwright_page, dashboard_landscape_live_url)
        if isinstance(expect, bool):
            assert out[key] is expect
        else:
            assert out[key] >= expect


# ---------------------------------------------------------------------------
# Landscape metadata / endpoint consistency workflow
# ---------------------------------------------------------------------------


class TestLandscapeVolpcaWorkflow:
    """Everything ``meta`` advertises must be usable by the endpoints that follow it.

    The landscape page builds its axis and colour menus purely from ``meta``, so an
    option listed there but rejected later is a dead control in the UI.
    """

    @staticmethod
    def _meta(client) -> dict:
        r = client.get("/api/landscape_volpca/meta")
        assert r.status_code == 200, r.get_data(as_text=True)[:400]
        return r.get_json()

    @staticmethod
    def _scatter(client, **params):
        return client.get("/api/landscape_volpca/scatter", query_string=params)

    def test_every_axis_meta_advertises_is_accepted_by_the_scatter_endpoint(
        self, flask_client_landscape
    ):
        meta = self._meta(flask_client_landscape)
        axes = [f"pc:{i}" for i in range(int(meta.get("n_pc") or 0))]
        axes += [f"umap:{i}" for i in range(int(meta.get("n_umap") or 0))]
        assert len(axes) >= 2, "landscape fixture advertised too few axes"

        for axis in axes:
            partner = next(other for other in axes if other != axis)
            assert (
                self._scatter(
                    flask_client_landscape, axis_x=axis, axis_y=partner
                ).status_code
                == 200
            ), f"meta offers axis {axis!r} but scatter rejects it"

    def test_every_colour_option_meta_advertises_is_accepted(
        self, flask_client_landscape
    ):
        meta = self._meta(flask_client_landscape)
        options = [opt.get("value") for opt in (meta.get("color_options") or []) if opt]
        assert options, "landscape fixture advertised no colour options"
        if not meta.get("has_state_color"):
            # The page disables this entry from the same flag; see the test below.
            options = [value for value in options if value != "state"]

        for value in options:
            r = self._scatter(
                flask_client_landscape, axis_x="pc:0", axis_y="pc:1", color=value
            )
            assert r.status_code == 200, (
                f"meta offers colour {value!r} but scatter rejects it: "
                f"{r.get_data(as_text=True)[:200]}"
            )
            if value == "none":
                continue
            marker = decode_plotly_figure(r.get_json())["data"][0].get("marker") or {}
            assert isinstance(marker.get("color"), list), (
                f"colour {value!r} was accepted but produced a single flat marker "
                "colour, so the choice had no effect on the plot"
            )

    def test_state_colour_is_offered_and_refused_by_the_same_flag(
        self, flask_client_landscape
    ):
        """``color_options`` always lists agglomerative state; ``has_state_color``
        decides whether it works, and the page disables the entry from that flag. If
        the two ever disagree the menu grows an option that only produces an error."""
        meta = self._meta(flask_client_landscape)
        values = [opt.get("value") for opt in (meta.get("color_options") or [])]
        assert "state" in values

        r = self._scatter(
            flask_client_landscape, axis_x="pc:0", axis_y="pc:1", color="state"
        )
        if meta.get("has_state_color"):
            assert r.status_code == 200
        else:
            assert r.status_code == 400
            assert "state" in (r.get_json().get("error") or "").lower()

    def test_scatter_plots_exactly_the_advertised_number_of_volumes(
        self, flask_client_landscape
    ):
        meta = self._meta(flask_client_landscape)
        n_volumes = int(meta.get("n_volumes") or 0)
        if n_volumes <= 0:
            pytest.skip("landscape fixture reports no volumes")

        fig = decode_plotly_figure(
            self._scatter(
                flask_client_landscape, axis_x="pc:0", axis_y="pc:1"
            ).get_json()
        )
        plotted = sum(len(trace.get("x") or []) for trace in fig["data"])
        assert (
            plotted == n_volumes
        ), f"meta advertises {n_volumes} volumes but the scatter plots {plotted}"

    def test_an_axis_beyond_the_advertised_range_is_refused(
        self, flask_client_landscape
    ):
        meta = self._meta(flask_client_landscape)
        beyond = f"pc:{int(meta.get('n_pc') or 0) + 5}"
        r = self._scatter(flask_client_landscape, axis_x=beyond, axis_y="pc:0")
        assert r.status_code == 400, "an out-of-range component should not render"
