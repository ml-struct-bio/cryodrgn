"""Tests for trajectory-integrated volume display and analyze-volume APIs."""

from __future__ import annotations

import base64
from pathlib import Path

import numpy as np
import pytest

from cryodrgn.dashboard.volume_slice_viewer import (
    _png_is_mostly_blank,
    analyze_volume_by_id_payload,
    analyze_volume_markers_payload,
    analyze_volumes_batch_payload,
    analyze_volumes_catalog_payload,
    apply_reconstruction_window,
    default_analyze_volume_id,
    discover_analyze_volume_catalog,
    discover_analyze_volume_markers,
    downsample_volume_box_average,
    PLOT3D_TARGET_D,
    reconstruction_window_params,
    spherical_window_mask_3d,
    volume_array_b64,
    vtk_transfer_volume_payload,
)


_VTK_BUNDLE = (
    Path(__file__).resolve().parents[1]
    / "cryodrgn/dashboard/static/js/volume_raycast_vtk.bundle.js"
)

# Matches ``CryoVolume3dUtils.PLOT3D_TARGET_D`` / client box-average downsample.
_PLOT3D_TARGET_D = PLOT3D_TARGET_D


def _downsample_volume_box_average(
    vol: np.ndarray, target_d: int = _PLOT3D_TARGET_D
) -> np.ndarray:
    """Test helper mirroring the dashboard downsample export."""
    return downsample_volume_box_average(vol, target_d)


class TestVolumeSliceViewerPure:
    def test_vtk_bundle_shipped_with_package(self) -> None:
        """3D raycast bundle is vendored in static/ (no separate npm install)."""
        assert (
            _VTK_BUNDLE.is_file()
        ), f"Missing {_VTK_BUNDLE.name} — run scripts/vendor_volume_raycast_bundle.sh"
        assert _VTK_BUNDLE.stat().st_size > 100_000
        body = _VTK_BUNDLE.read_text(encoding="utf-8", errors="ignore")
        assert "CryoVolumeRaycastView" in body

    def test_volume_array_b64_roundtrip(self) -> None:
        vol = np.arange(8, dtype=np.float32).reshape(2, 2, 2)
        b64 = volume_array_b64(vol)
        raw = base64.standard_b64decode(b64)
        back = np.frombuffer(raw, dtype=np.float32).reshape(2, 2, 2)
        np.testing.assert_array_equal(back, vol)

    def test_png_is_mostly_blank_detects_uniform_white(self, tmp_path: Path) -> None:
        from PIL import Image

        path = tmp_path / "white.png"
        Image.new("RGB", (64, 64), color=(255, 255, 255)).save(path)
        assert _png_is_mostly_blank(str(path)) is True

    def test_png_is_mostly_blank_accepts_structured_image(self, tmp_path: Path) -> None:
        from PIL import Image

        path = tmp_path / "pattern.png"
        arr = np.zeros((64, 64, 3), dtype=np.uint8)
        arr[:32, :, 0] = 200
        Image.fromarray(arr).save(path)
        assert _png_is_mostly_blank(str(path)) is False

    def test_downsample_volume_box_average(self) -> None:
        vol = np.ones((8, 8, 8), dtype=np.float32) * 3.0
        ds = downsample_volume_box_average(vol, target_d=128)
        assert ds.shape == (128, 128, 128)
        np.testing.assert_allclose(ds, 3.0)
        same = downsample_volume_box_average(ds, target_d=128)
        np.testing.assert_array_equal(same, ds)

    def test_vtk_transfer_volume_payload_downsamples_large_cubes(self) -> None:
        vol = np.arange(16**3, dtype=np.float32).reshape(16, 16, 16)
        payload = vtk_transfer_volume_payload(vol, target_d=8)
        assert payload["D"] == 8
        assert payload["source_D"] == 16
        assert payload["downsample"] == "box_average"
        raw = base64.standard_b64decode(payload["volume_b64"])
        back = np.frombuffer(raw, dtype=np.float32).reshape(8, 8, 8)
        expected = downsample_volume_box_average(vol, target_d=8)
        np.testing.assert_allclose(back, expected)

    def test_analyze_volumes_batch_returns_target_d(self, dashboard_experiment) -> None:
        catalog = discover_analyze_volume_catalog(dashboard_experiment)
        vol_ids = [catalog[0]["id"]]
        batch = analyze_volumes_batch_payload(
            dashboard_experiment, vol_ids, n_cpus=1, target_d=128
        )
        assert batch["target_d"] == 128
        entry = batch["volumes"][vol_ids[0]]
        assert entry["D"] <= 128
        assert "source_D" in entry

    def test_spherical_window_mask_3d_soft_edge(self) -> None:
        mask = spherical_window_mask_3d(D=8, in_rad=0.5, out_rad=1.0)
        assert mask.shape == (8, 8, 8)
        assert mask.max() == pytest.approx(1.0)
        assert mask.min() == pytest.approx(0.0)
        center = mask[4, 4, 4]
        corner = mask[0, 0, 0]
        assert center > corner

    def test_apply_reconstruction_window_zeros_outside_sphere(
        self, dashboard_experiment
    ) -> None:
        enabled, in_rad, out_rad = reconstruction_window_params(dashboard_experiment)
        assert enabled is True
        assert in_rad == pytest.approx(0.85)
        assert out_rad == pytest.approx(0.99)
        vol = np.ones((16, 16, 16), dtype=np.float32)
        masked = apply_reconstruction_window(vol, dashboard_experiment)
        mask = spherical_window_mask_3d(D=16, in_rad=in_rad, out_rad=out_rad)
        np.testing.assert_allclose(masked, mask)
        assert masked[0, 0, 0] == pytest.approx(0.0)
        assert masked[8, 8, 8] == pytest.approx(1.0)


class TestVolumeSliceViewerRoutes:
    def test_discover_analyze_catalog(self, dashboard_experiment) -> None:
        catalog = discover_analyze_volume_catalog(dashboard_experiment)
        kinds = {e["kind"] for e in catalog}
        assert "kmeans" in kinds
        assert "pc" in kinds

    def test_analyze_volumes_payload(self, dashboard_experiment) -> None:
        payload = analyze_volumes_catalog_payload(dashboard_experiment)
        assert payload["ok"] is True
        assert payload["catalog"]
        assert payload["markers"]
        assert len(payload["markers"]) == len(payload["catalog"])
        assert "volumes" not in payload
        km = [
            e for e in payload["catalog"] if e.get("kind") == "kmeans" and "znorm" in e
        ]
        assert km
        assert payload["default_vol_id"] == min(km, key=lambda e: e["znorm"])["id"]

    def test_default_analyze_volume_id_lowest_znorm(self, dashboard_experiment) -> None:
        payload = analyze_volumes_catalog_payload(
            dashboard_experiment, include_markers=False
        )
        km = [
            e for e in payload["catalog"] if e.get("kind") == "kmeans" and "znorm" in e
        ]
        assert km
        assert (
            default_analyze_volume_id(payload["catalog"])
            == min(km, key=lambda e: e["znorm"])["id"]
        )

    def test_analyze_volumes_payload_fast_catalog(self, dashboard_experiment) -> None:
        payload = analyze_volumes_catalog_payload(
            dashboard_experiment, include_markers=False
        )
        assert payload["ok"] is True
        assert payload["catalog"]
        assert payload["markers"] == []

    def test_analyze_volume_markers_payload(self, dashboard_experiment) -> None:
        payload = analyze_volume_markers_payload(dashboard_experiment)
        catalog = discover_analyze_volume_catalog(dashboard_experiment)
        assert payload["ok"] is True
        assert len(payload["markers"]) == len(catalog)

    def test_analyze_volume_markers(self, dashboard_experiment) -> None:
        catalog = discover_analyze_volume_catalog(dashboard_experiment)
        markers = discover_analyze_volume_markers(dashboard_experiment, catalog)
        assert len(markers) == len(catalog)
        kinds = {m["kind"] for m in markers}
        assert "kmeans" in kinds
        assert "pc" in kinds
        km = [m for m in markers if m["kind"] == "kmeans"]
        assert all(m["plot_row"] is not None for m in km)
        assert all(m["label"].startswith("K") for m in km)
        assert all(int(m["label"][1:]) >= 1 for m in km)

    def test_kmeans_volume_ids_for_anchor_indices(self, dashboard_experiment) -> None:
        from cryodrgn.dashboard.volume_slice_viewer import (
            kmeans_volume_ids_for_anchor_indices,
        )

        markers = discover_analyze_volume_markers(dashboard_experiment)
        km = [m for m in markers if m["kind"] == "kmeans"]
        assert len(km) >= 2
        anchor_rows = [int(m["plot_row"]) for m in km[:3]]
        vol_ids = kmeans_volume_ids_for_anchor_indices(
            dashboard_experiment, anchor_rows
        )
        assert vol_ids is not None
        assert len(vol_ids) == len(anchor_rows)
        assert all(v.startswith("kmeans:") for v in vol_ids)

    def test_analyze_volume_by_id(self, dashboard_experiment) -> None:
        catalog = discover_analyze_volume_catalog(dashboard_experiment)
        assert catalog
        vol_id = catalog[0]["id"]
        payload = analyze_volume_by_id_payload(dashboard_experiment, vol_id)
        assert payload["ok"] is True
        assert payload["id"] == vol_id
        assert payload["volume_b64"]
        assert payload["D"] > 0

    def test_analyze_volumes_batch_matches_single(self, dashboard_experiment) -> None:
        catalog = discover_analyze_volume_catalog(dashboard_experiment)
        vol_ids = [e["id"] for e in catalog[:3]]
        batch = analyze_volumes_batch_payload(dashboard_experiment, vol_ids, n_cpus=2)
        assert batch["ok"] is True
        for vol_id in vol_ids:
            single = analyze_volume_by_id_payload(dashboard_experiment, vol_id)
            assert batch["volumes"][vol_id]["volume_b64"] == single["volume_b64"]

    def test_analyze_volumes_api(self, flask_client) -> None:
        r = flask_client.get("/api/volume_viewer/analyze_volumes")
        assert r.status_code == 200
        j = r.get_json()
        assert j["ok"] is True
        assert j["catalog"]
        assert j["markers"]
        assert len(j["markers"]) == len(j["catalog"])
        assert j.get("volumes") in (None, {})
        km = [e for e in j["catalog"] if e.get("kind") == "kmeans" and "znorm" in e]
        assert km
        assert j["default_vol_id"] == min(km, key=lambda e: e["znorm"])["id"]

    def test_analyze_volumes_api_fast_catalog(self, flask_client) -> None:
        r = flask_client.get("/api/volume_viewer/analyze_volumes?include_markers=0")
        assert r.status_code == 200
        j = r.get_json()
        assert j["ok"] is True
        assert j["catalog"]
        assert j["markers"] == []
        km = [e for e in j["catalog"] if e.get("kind") == "kmeans" and "znorm" in e]
        assert km
        assert j["default_vol_id"] == min(km, key=lambda e: e["znorm"])["id"]

    def test_analyze_markers_api(self, flask_client) -> None:
        cat = flask_client.get(
            "/api/volume_viewer/analyze_volumes?include_markers=0"
        ).get_json()
        r = flask_client.get("/api/volume_viewer/analyze_markers")
        assert r.status_code == 200
        j = r.get_json()
        assert j["ok"] is True
        assert len(j["markers"]) == len(cat["catalog"])

    def test_analyze_volume_api(self, flask_client) -> None:
        cat = flask_client.get("/api/volume_viewer/analyze_volumes").get_json()
        assert cat["catalog"]
        vol_id = cat["catalog"][0]["id"]
        r = flask_client.get(f"/api/volume_viewer/analyze_volume?id={vol_id}")
        assert r.status_code == 200
        j = r.get_json()
        assert j["ok"] is True
        assert j["id"] == vol_id
        assert j["volume_b64"]

    def test_analyze_volumes_batch_api(self, flask_client) -> None:
        cat = flask_client.get(
            "/api/volume_viewer/analyze_volumes?include_markers=0"
        ).get_json()
        vol_ids = [e["id"] for e in cat["catalog"][:2]]
        r = flask_client.post(
            "/api/volume_viewer/analyze_volumes_batch",
            json={"ids": vol_ids},
        )
        assert r.status_code == 200
        j = r.get_json()
        assert j["ok"] is True
        for vol_id in vol_ids:
            assert vol_id in j["volumes"]
            assert j["volumes"][vol_id]["volume_b64"]

    def test_analyze_volumes_batch_api_rejects_empty_ids(self, flask_client) -> None:
        r = flask_client.post(
            "/api/volume_viewer/analyze_volumes_batch",
            json={"ids": []},
        )
        assert r.status_code == 400
        assert "ids" in r.get_json().get("error", "").lower()

    def test_analyze_volumes_batch_api_rejects_missing_ids(self, flask_client) -> None:
        r = flask_client.post(
            "/api/volume_viewer/analyze_volumes_batch",
            json={},
        )
        assert r.status_code == 400
        assert "ids" in r.get_json().get("error", "").lower()

    def test_analyze_volumes_chimerax_batch_api_rejects_empty_ids(
        self, flask_client
    ) -> None:
        r = flask_client.post(
            "/api/volume_viewer/analyze_volumes_chimerax_batch",
            json={"ids": []},
        )
        assert r.status_code == 400
        assert "ids" in r.get_json().get("error", "").lower()

    def test_volume_viewer_url_redirects_to_trajectory(self, flask_client) -> None:
        r = flask_client.get("/volume-viewer", follow_redirects=False)
        assert r.status_code in (301, 302, 303, 307, 308)
        assert "/trajectory" in (r.headers.get("Location") or "")

    def test_volume_viewer_url_redirects_when_volumes_eligible(
        self, flask_client_volumes_eligible
    ) -> None:
        r = flask_client_volumes_eligible.get("/volume-viewer", follow_redirects=False)
        assert r.status_code in (301, 302, 303, 307, 308)
        assert "/trajectory" in (r.headers.get("Location") or "")

    def test_trajectory_page_integrates_volume_display(
        self, flask_client_volumes_eligible
    ) -> None:
        r = flask_client_volumes_eligible.get("/trajectory")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "Trajectory creator" in body
        assert "vslice-canvas" in body
        assert "vslice-volume-picker-rows" in body
        assert "btn-vslice-pan-up" in body
        assert "vslice-vtk-container" in body
        assert "vslice-view-mode-3d" not in body
        assert "traj-vol-backend-vtk" in body
        assert "traj-vol-backend-chimerax" in body
        assert "traj-vol-backend-slice" in body
        assert "traj-vol-backend-hint" in body
        assert (
            'id="traj-vol-backend-chimerax"' in body
            and "checked" in body.split("traj-vol-backend-chimerax")[1].split(">")[0]
        )
        assert (
            'id="traj-vol-backend-chimerax"' in body
            and "disabled" in body.split("traj-vol-backend-chimerax")[1].split(">")[0]
        )
        assert "btn-traj-vol-dock-below" in body
        assert "traj-vol-volume-nav" in body
        assert "traj-vol-volume-slider" in body
        assert "traj-vol-expanded-host" in body
        assert "traj-vol-popout-modal" in body
        assert "traj-mode-manual" in body
        assert 'id="traj-mode-manual" value="manual" checked' in body.replace(
            " ", ""
        ) or (
            'id="traj-mode-manual"' in body
            and "checked" in body.split("traj-mode-manual")[1].split(">")[0]
        )
        assert "Choosing trajectory waypoints" in body
        assert "using particles and particle sets" in body
        assert "Tracing trajectory path directly" in body
        assert "traj-scatter-interp-nearest" in body
        assert "traj-scatter-interp-direct" in body
        assert "btn-traj-manual-volume" in body
        assert "btn-rerender-volumes" in body
        assert "deferManualWaypointVolumeRerender" in body
        assert "finishManualVolumeDeselection" in body
        assert "plotRowsForManualVolIds" in body
        assert "syncManualVolumeViewerForWaypointSelection" in body
        assert "manualSnappedDecodePathActive" in body
        assert "preserveManualSnappedDecodeVolumeCatalog" in body
        assert "reverseManualInterpolatedVolumeCatalogPreservingAnchors" in body
        assert "syncManualAnchorIndicesToVolIdOrder" in body
        assert "compactManualAnchorVolumesFromSnapshot" in body
        assert "manualInterpolatedVolumeReversePending" in body
        assert "reverseDecodedTrajectoryVolumesInMemory" in body
        assert "Render volumes" in body
        assert "btn-traj-manual-graph" in body
        assert "manual-snap-n-points" in body
        assert "manual-graph-n-points" in body
        assert "Add points along" in body
        assert "direct line" in body
        assert "Add points using" in body
        assert "graph traversal" in body
        assert "Add trajectory points" in body
        assert "through interpolation" in body
        assert "traj-manual-interp-stack" in body
        assert "cryo-traj-manual-mode-fields" in body
        assert "traj-manual-interp-graph" not in body
        assert "lastDisplayedVolumeFocusIndex" in body
        assert "rememberDisplayedVolumeFocusIndex" in body
        assert "restoreRememberedVolumeFocusIndex" in body
        assert "snapVolumeFocusToActiveTick" in body
        assert "trajectoryAddPointsBusy" in body
        assert "snapFocusToReadyTick" in body
        assert "trajectoryPathPointCountForHighlight" in body
        assert "btn-traj-graph" not in body
        assert "traj-mode-alt" not in body
        assert "volume_slice_canvas.js" in body
        assert "trajectory_volume_display.js" in body
        assert "stableBackendChrome: true" in body
        assert "cryo-traj-vol-backend-panel--inactive" in body
        assert "cryo-traj-vol-backend-panel--reserved" in body
        assert 'id="vslice-chimerax-view-controls"' in body
        assert "vslice-iso-controls" in body
        assert "traj-chimerax-view-rotate-row" in body
        assert "volume_3d_utils.js" in body
        assert "/api/volume_viewer/analyze_volumes" in body

    def test_decode_api_requires_row(self, flask_client_volumes_eligible) -> None:
        r = flask_client_volumes_eligible.post(
            "/api/volume_viewer/decode",
            json={},
        )
        assert r.status_code == 400
        assert "row" in r.get_json().get("error", "").lower()

    def test_decode_api_success(
        self,
        flask_client_volumes_eligible,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        fake = np.zeros((8, 8, 8), dtype=np.float32)
        monkeypatch.setattr(
            "cryodrgn.dashboard.volume_slice_viewer._decode_volume_array",
            lambda _exp, row: fake,
        )
        r = flask_client_volumes_eligible.post(
            "/api/volume_viewer/decode",
            json={"row": 0},
        )
        assert r.status_code == 200
        j = r.get_json()
        assert j["ok"] is True
        assert j["row"] == 0
        assert j["D"] == 8
        assert j["volume_b64"]
        assert j["volume_dtype"] == "float32"
        raw = base64.standard_b64decode(j["volume_b64"])
        back = np.frombuffer(raw, dtype=np.float32).reshape(8, 8, 8)
        np.testing.assert_array_equal(back, fake)

    def test_landing_has_no_standalone_volume_viewer_card(
        self, flask_client_volumes_eligible
    ) -> None:
        r = flask_client_volumes_eligible.get("/")
        body = r.get_data(as_text=True)
        assert "Trajectory creator" in body
        assert "/trajectory" in body
        assert 'href="/volume-viewer"' not in body
        assert ">Volume viewer</h2>" not in body


class TestTrajectoryVolumeBrowserSmoke:
    """Headless Chromium: manual k-means/PC picker and VTK chrome on trajectory page."""

    pytestmark = [pytest.mark.browser, pytest.mark.slow]

    @pytest.fixture(autouse=True)
    def _stub_volume_viewer_apis(self, playwright_page):
        from tests.conftest import playwright_route_volume_viewer_render_stub

        playwright_route_volume_viewer_render_stub(playwright_page)

    def test_catalog_and_canvas_load(
        self, playwright_page, dashboard_volumes_eligible_live_url
    ) -> None:
        from tests.conftest import dashboard_smoke_volume_viewer

        out = dashboard_smoke_volume_viewer(
            playwright_page, dashboard_volumes_eligible_live_url
        )
        assert out is not None
        assert out["volume_picker_buttons"] >= 1
        assert out["canvas_present"] is True
        assert (out.get("scatter_points") or 0) > 0

    def test_manual_mode_trajectory_overlay_on_load(
        self, playwright_page, dashboard_volumes_eligible_live_url
    ) -> None:
        from tests.conftest import (
            DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
            _dashboard_smoke_volume_viewer_ready,
        )

        ready = _dashboard_smoke_volume_viewer_ready(
            playwright_page,
            dashboard_volumes_eligible_live_url,
            timeout_ms=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        assert ready is not None
        playwright_page.wait_for_function(
            """() => {
              var overlay = document.getElementById('traj-glyph-overlay');
              return !!(overlay && overlay.querySelector('.cryo-traj-glyph-path'));
            }""",
            timeout=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        state = playwright_page.evaluate(
            """() => {
              var overlay = document.getElementById('traj-glyph-overlay');
              var gd = document.getElementById('scatter');
              var lineTraces = 0;
              if (gd && gd.data) {
                for (var i = 1; i < gd.data.length; i++) {
                  if (gd.data[i] && gd.data[i].mode && gd.data[i].mode.indexOf('lines') >= 0) {
                    lineTraces++;
                  }
                }
              }
              return {
                svgPath: !!(overlay && overlay.querySelector('.cryo-traj-glyph-path')),
                glyphMarkers: overlay
                  ? overlay.querySelectorAll('.cryo-traj-glyph-marker').length
                  : 0,
                plotlyLineTraces: lineTraces,
                activePickerBtns: document.querySelectorAll(
                  '.cryo-vslice-vol-btn--active'
                ).length,
                anchorPickHidden: (function() {
                  var el = document.getElementById('traj-anchor-pick');
                  return !!(el && el.hidden);
                })(),
                manualPickerVisible: (function() {
                  var el = document.getElementById('traj-manual-picker');
                  return !!(el && !el.hidden);
                })()
              };
            }"""
        )
        assert state["svgPath"] is True
        assert state["activePickerBtns"] >= 2
        assert state["glyphMarkers"] >= 2
        assert state["plotlyLineTraces"] >= 1
        assert state["manualPickerVisible"] is True
        gen_btn = playwright_page.evaluate(
            """() => {
              var btn = document.getElementById('btn-generate-volumes');
              var saveBtn = document.getElementById('btn-save-volumes');
              var rerenderBtn = document.getElementById('btn-rerender-volumes');
              return {
                present: !!btn,
                hidden: btn ? btn.hidden : true,
                disabled: btn ? btn.disabled : true,
                savePresent: !!saveBtn,
                saveHidden: saveBtn ? saveBtn.hidden : true,
                saveDisabled: saveBtn ? saveBtn.disabled : true,
                rerenderPresent: !!rerenderBtn,
                rerenderHidden: rerenderBtn ? rerenderBtn.hidden : true,
                rerenderDisabled: rerenderBtn ? rerenderBtn.disabled : true
              };
            }"""
        )
        assert gen_btn["present"] is True
        assert gen_btn["hidden"] is False
        assert gen_btn["disabled"] is True
        assert gen_btn["savePresent"] is True
        assert gen_btn["saveHidden"] is False
        assert gen_btn["saveDisabled"] is True
        assert gen_btn["rerenderPresent"] is True
        assert gen_btn["rerenderHidden"] is False
        assert gen_btn["rerenderDisabled"] is False

    def test_manual_mode_chimerax_default_with_analyze_volumes(
        self, playwright_page, dashboard_volumes_eligible_live_url
    ) -> None:
        from tests.conftest import (
            DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
            _dashboard_smoke_volume_viewer_ready,
            dashboard_smoke_rerender_manual_volumes,
        )

        ready = _dashboard_smoke_volume_viewer_ready(
            playwright_page,
            dashboard_volumes_eligible_live_url,
            timeout_ms=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        assert ready is not None
        dashboard_smoke_rerender_manual_volumes(
            playwright_page,
            timeout_ms=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        playwright_page.wait_for_function(
            """() => {
              var cx = document.getElementById('traj-vol-backend-chimerax');
              var preview = document.getElementById('vslice-chimerax-preview');
              return !!(cx && cx.checked && !cx.disabled
                && preview && !preview.hidden && preview.src);
            }""",
            timeout=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )

    def test_3d_mode_chrome_and_bundle(
        self, playwright_page, dashboard_volumes_eligible_live_url
    ) -> None:
        from tests.conftest import (
            DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
            _dashboard_smoke_volume_viewer_ready,
            dashboard_smoke_rerender_manual_volumes,
        )

        ready = _dashboard_smoke_volume_viewer_ready(
            playwright_page,
            dashboard_volumes_eligible_live_url,
            timeout_ms=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        assert ready is not None
        dashboard_smoke_rerender_manual_volumes(
            playwright_page,
            timeout_ms=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        playwright_page.wait_for_function(
            """() => {
              var vtkBackend = document.getElementById('traj-vol-backend-vtk');
              return !!(vtkBackend && !vtkBackend.disabled);
            }""",
            timeout=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        playwright_page.click("label[for='traj-vol-backend-vtk']")
        playwright_page.wait_for_function(
            """() => {
              var vtkBackend = document.getElementById('traj-vol-backend-vtk');
              return !!(vtkBackend && vtkBackend.checked && !vtkBackend.disabled);
            }""",
            timeout=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        chrome = playwright_page.evaluate(
            """() => {
              var vtkBackend = document.getElementById('traj-vol-backend-vtk');
              var sliceBackend = document.getElementById('traj-vol-backend-slice');
              var vtk = document.getElementById('vslice-vtk-container');
              var iso = document.getElementById('vslice-iso-controls');
              var hint = document.getElementById('traj-vol-backend-hint');
              return {
                vtk_backend_checked: !!(vtkBackend && vtkBackend.checked),
                vtk_backend_disabled: !!(vtkBackend && vtkBackend.disabled),
                vtk_container_present: !!vtk,
                iso_present: !!iso,
                hint_visible: !!(hint && !hint.hidden),
                vtk_container_hidden: vtk ? vtk.hidden : true
              };
            }"""
        )
        assert chrome["vtk_backend_checked"] is True
        assert chrome["vtk_backend_disabled"] is False
        assert chrome["hint_visible"] is False
        assert chrome["vtk_container_present"] is True
        assert chrome["iso_present"] is True
        assert chrome["vtk_container_hidden"] is False
        slice_controls = playwright_page.evaluate(
            """() => {
              var row = document.getElementById('vslice-slice-controls-row');
              return {
                present: !!row,
                hidden: row ? row.hidden : true
              };
            }"""
        )
        assert slice_controls["present"] is True
        assert slice_controls["hidden"] is True
        pad = playwright_page.evaluate(
            """() => {
              var col = document.getElementById('vslice-pad-column');
              return {
                pad_present: !!col,
                pad_visible: !!(col && !col.hidden),
                pan_up: !!document.getElementById('btn-vslice-pan-up')
              };
            }"""
        )
        assert pad["pad_present"] is True
        assert pad["pad_visible"] is True
        assert pad["pan_up"] is True

    def test_manual_volume_picker_toggles_selection(
        self, playwright_page, dashboard_volumes_eligible_live_url
    ) -> None:
        from tests.conftest import (
            DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
            dashboard_smoke_volume_viewer_picker_switch,
        )

        out = dashboard_smoke_volume_viewer_picker_switch(
            playwright_page,
            dashboard_volumes_eligible_live_url,
            timeout_ms=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        assert out is not None
        if not out.get("switched"):
            pytest.skip(out.get("reason", "not enough volumes to switch"))
        assert out["from_label"]
        assert out["to_label"]

    def test_manual_mode_chimerax_gallery_renders_images(
        self, playwright_page, dashboard_volumes_eligible_live_url
    ) -> None:
        from tests.conftest import (
            DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
            _dashboard_smoke_volume_viewer_ready,
            dashboard_smoke_rerender_manual_volumes,
            fulfill_volume_viewer_render_route,
        )

        tiny_png = (
            "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8"
            "z8BQDwAEhQGAhKmMIQAAAABJRU5ErkJggg=="
        )

        def _route_handler(route):
            if (
                route.request.method == "POST"
                and "analyze_volumes_chimerax_batch" in route.request.url
            ):
                body = (
                    '{"ok": true, "images": ["' + tiny_png + '", "' + tiny_png + '"]}'
                )
                route.fulfill(
                    status=200,
                    content_type="application/json",
                    body=body,
                )
            elif not fulfill_volume_viewer_render_route(route):
                route.continue_()

        playwright_page.route("**/api/volume_viewer/**", _route_handler)
        ready = _dashboard_smoke_volume_viewer_ready(
            playwright_page,
            dashboard_volumes_eligible_live_url,
            timeout_ms=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        assert ready is not None
        dashboard_smoke_rerender_manual_volumes(
            playwright_page,
            timeout_ms=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        playwright_page.wait_for_function(
            """() => {
              var cx = document.getElementById('traj-vol-backend-chimerax');
              return !!(cx && !cx.disabled);
            }""",
            timeout=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        playwright_page.click("label[for='traj-vol-backend-chimerax']")
        playwright_page.wait_for_function(
            """() => {
              var col = document.getElementById('traj-vol-column');
              var preview = document.getElementById('vslice-chimerax-preview');
              if (!preview || preview.hidden || !preview.src) return false;
              if (col && !col.hidden) return false;
              return preview.src.indexOf("blob:") === 0
                || preview.src.indexOf("data:image/png;base64,") === 0;
            }""",
            timeout=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        layout = playwright_page.evaluate(
            """() => {
              function disp(id) {
                var el = document.getElementById(id);
                return el ? window.getComputedStyle(el).display : null;
              }
              return {
                vtk: disp('vslice-vtk-container'),
                canvas: disp('vslice-canvas'),
              };
            }"""
        )
        assert layout["vtk"] == "none", f"vtk container not hidden: {layout}"
        assert layout["canvas"] == "none", f"slice canvas not hidden: {layout}"
        playwright_page.click("#btn-traj-vol-dock-below")
        playwright_page.wait_for_function(
            """() => {
              var modal = document.getElementById('traj-vol-popout-modal');
              var col = document.getElementById('traj-vol-column');
              if (!modal || modal.hidden || !col || col.hidden) return false;
              var imgs = col.querySelectorAll('img.cryo-traj-vol-main');
              if (imgs.length < 2) return false;
              for (var i = 0; i < imgs.length; i++) {
                var src = imgs[i].src || "";
                if (src.indexOf("blob:") !== 0
                  && src.indexOf("data:image/png;base64,") !== 0) {
                  return false;
                }
              }
              return true;
            }""",
            timeout=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )

    def test_manual_mode_chimerax_rotate_controls_send_view_turns(
        self, playwright_page, dashboard_volumes_eligible_live_url
    ) -> None:
        import json
        import time

        from tests.conftest import (
            DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
            _dashboard_smoke_volume_viewer_ready,
            dashboard_smoke_rerender_manual_volumes,
            fulfill_volume_viewer_render_route,
        )

        tiny_png = (
            "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8"
            "z8BQDwAEhQGAhKmMIQAAAABJRU5ErkJggg=="
        )
        captured: list[dict] = []

        def _route_handler(route):
            if (
                route.request.method == "POST"
                and "analyze_volumes_chimerax_batch" in route.request.url
            ):
                captured.append(json.loads(route.request.post_data or "{}"))
                route.fulfill(
                    status=200,
                    content_type="application/json",
                    body=(
                        '{"ok": true, "images": ["' + tiny_png + '"], '
                        '"view_matrix": "camera 1,0,0,0,0,1,0,0,0,0,1,0"}'
                    ),
                )
            elif not fulfill_volume_viewer_render_route(route):
                route.continue_()

        playwright_page.route("**/api/volume_viewer/**", _route_handler)
        ready = _dashboard_smoke_volume_viewer_ready(
            playwright_page,
            dashboard_volumes_eligible_live_url,
            timeout_ms=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        assert ready is not None
        dashboard_smoke_rerender_manual_volumes(
            playwright_page,
            timeout_ms=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        playwright_page.click("label[for='traj-vol-backend-chimerax']")
        playwright_page.wait_for_function(
            """() => {
              var preview = document.getElementById('vslice-chimerax-preview');
              var rotateBtn = document.querySelector('[data-traj-chimerax-view-axis="y"]');
              return !!(preview && !preview.hidden && preview.src
                && rotateBtn && !rotateBtn.disabled);
            }""",
            timeout=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        assert captured, "expected initial ChimeraX batch request"
        initial = captured[0]
        assert not initial.get(
            "view_matrix"
        ), f"initial switch should not sync VTK view matrix: {initial!r}"
        assert not initial.get(
            "view_turns"
        ), f"initial switch should use default ChimeraX view: {initial!r}"
        playwright_page.click('[data-traj-chimerax-view-axis="y"]')
        deadline = time.time() + 60.0
        while time.time() < deadline and len(captured) < 2:
            playwright_page.wait_for_timeout(250)
        assert len(captured) >= 2, "expected ChimeraX re-render after Rotate Y"
        rotate_body = captured[-1]
        turns = rotate_body.get("view_turns") or []
        by_axis = {
            str(t.get("axis")).lower(): float(t.get("degrees", 0)) for t in turns
        }
        assert (
            abs(by_axis.get("y", 0) - 180) < 1.0
        ), f"unexpected y turn after rotate: {turns}"

    def test_manual_mode_vtk_expand_then_chimerax_shows_aside_preview(
        self, playwright_page, dashboard_volumes_eligible_live_url
    ) -> None:
        from tests.conftest import (
            DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
            _dashboard_smoke_volume_viewer_ready,
            dashboard_smoke_rerender_manual_volumes,
            fulfill_volume_viewer_render_route,
        )

        tiny_png = (
            "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8"
            "z8BQDwAEhQGAhKmMIQAAAABJRU5ErkJggg=="
        )

        def _route_handler(route):
            if (
                route.request.method == "POST"
                and "analyze_volumes_chimerax_batch" in route.request.url
            ):
                route.fulfill(
                    status=200,
                    content_type="application/json",
                    body='{"ok": true, "images": ["' + tiny_png + '"]}',
                )
            elif not fulfill_volume_viewer_render_route(route):
                route.continue_()

        playwright_page.route("**/api/volume_viewer/**", _route_handler)
        ready = _dashboard_smoke_volume_viewer_ready(
            playwright_page,
            dashboard_volumes_eligible_live_url,
            timeout_ms=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        assert ready is not None
        dashboard_smoke_rerender_manual_volumes(
            playwright_page,
            timeout_ms=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        playwright_page.wait_for_function(
            """() => {
              var vtk = document.getElementById('traj-vol-backend-vtk');
              var vtkHost = document.getElementById('vslice-vtk-container');
              return !!(vtk && vtk.checked && !vtk.disabled
                && vtkHost && !vtkHost.hidden);
            }""",
            timeout=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )
        dock_btn = playwright_page.query_selector("#btn-traj-vol-dock-below")
        if dock_btn and not dock_btn.is_hidden():
            playwright_page.click("#btn-traj-vol-dock-below")
            playwright_page.wait_for_function(
                """() => {
                  var modal = document.getElementById('traj-vol-popout-modal');
                  var host = document.getElementById('traj-vol-expanded-host');
                  var row = document.getElementById('vslice-display-row');
                  return !!(modal && !modal.hidden && host && !host.hidden
                    && row && host.contains(row));
                }""",
                timeout=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
            )
        playwright_page.click("label[for='traj-vol-backend-chimerax']")
        playwright_page.wait_for_function(
            """() => {
              var shell = document.getElementById('traj-vol-aside-shell');
              var row = document.getElementById('vslice-display-row');
              var expanded = document.getElementById('traj-vol-expanded-host');
              var preview = document.getElementById('vslice-chimerax-preview');
              if (!shell || shell.hidden || !row || row.hidden) return false;
              if (!shell.contains(row)) return false;
              if (expanded && expanded.contains(row)) return false;
              if (!preview || preview.hidden || !preview.src) return false;
              return preview.src.indexOf("blob:") === 0
                || preview.src.indexOf("data:image/png;base64,") === 0;
            }""",
            timeout=DASHBOARD_BROWSER_FAST_TIMEOUT_MS,
        )

    def test_manual_mode_chimerax_vtk_roundtrip_restores_vtk_preview(
        self, playwright_page, dashboard_volumes_eligible_live_url
    ) -> None:
        from tests.conftest import (
            DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
            _dashboard_smoke_volume_viewer_ready,
            dashboard_smoke_rerender_manual_volumes,
        )

        ready = _dashboard_smoke_volume_viewer_ready(
            playwright_page,
            dashboard_volumes_eligible_live_url,
            timeout_ms=DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
        )
        assert ready is not None
        dashboard_smoke_rerender_manual_volumes(
            playwright_page,
            timeout_ms=DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
        )
        playwright_page.wait_for_function(
            """() => {
              var vtk = document.getElementById('vslice-vtk-container');
              return !!(vtk && !vtk.hidden);
            }""",
            timeout=DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
        )
        playwright_page.click("label[for='traj-vol-backend-chimerax']")
        playwright_page.wait_for_function(
            """() => {
              var preview = document.getElementById('vslice-chimerax-preview');
              return !!(preview && !preview.hidden && preview.src);
            }""",
            timeout=DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
        )
        playwright_page.click("label[for='traj-vol-backend-vtk']")
        playwright_page.wait_for_function(
            """() => {
              var vtk = document.getElementById('vslice-vtk-container');
              var canvas = document.getElementById('vslice-canvas');
              var preview = document.getElementById('vslice-chimerax-preview');
              if (!vtk || vtk.hidden) return false;
              if (window.getComputedStyle(vtk).display === 'none') return false;
              if (preview && !preview.hidden) return false;
              if (canvas && !canvas.hidden) return false;
              var row = document.getElementById('vslice-display-row');
              return !!(row && !row.hidden);
            }""",
            timeout=DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
        )
