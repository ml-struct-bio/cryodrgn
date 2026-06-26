"""Tests for the dashboard volume viewer."""

from __future__ import annotations

import base64
from pathlib import Path

import numpy as np
import pytest

from cryodrgn.dashboard.volume_slice_viewer import (
    analyze_volume_by_id_payload,
    analyze_volume_markers_payload,
    analyze_volumes_batch_payload,
    analyze_volumes_catalog_payload,
    apply_reconstruction_window,
    default_analyze_volume_id,
    discover_analyze_volume_catalog,
    discover_analyze_volume_markers,
    reconstruction_window_params,
    spherical_window_mask_3d,
    volume_array_b64,
)


_VTK_BUNDLE = (
    Path(__file__).resolve().parents[1]
    / "cryodrgn/dashboard/static/js/volume_raycast_vtk.bundle.js"
)

# Matches ``CryoVolume3dUtils.PLOT3D_TARGET_D`` / client box-average downsample.
_PLOT3D_TARGET_D = 128


def _downsample_volume_box_average(
    vol: np.ndarray, target_d: int = _PLOT3D_TARGET_D
) -> np.ndarray:
    """Box-average resample to ``target_d``³ (mirrors ``volume_3d_utils.js``)."""
    vol = np.asarray(vol, dtype=np.float32)
    d = int(vol.shape[0])
    target_d = int(target_d)
    if d == target_d:
        return vol
    scale = d / target_d
    out = np.zeros((target_d, target_d, target_d), dtype=np.float32)
    for iz in range(target_d):
        z0 = int(np.floor(iz * scale))
        z1 = int(min(d, np.ceil((iz + 1) * scale)))
        for iy in range(target_d):
            y0 = int(np.floor(iy * scale))
            y1 = int(min(d, np.ceil((iy + 1) * scale)))
            for ix in range(target_d):
                x0 = int(np.floor(ix * scale))
                x1 = int(min(d, np.ceil((ix + 1) * scale)))
                block = vol[x0:x1, y0:y1, z0:z1]
                out[ix, iy, iz] = float(block.mean()) if block.size else 0.0
    return out


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

    def test_downsample_volume_box_average(self) -> None:
        vol = np.ones((8, 8, 8), dtype=np.float32) * 3.0
        ds = _downsample_volume_box_average(vol, target_d=128)
        assert ds.shape == (128, 128, 128)
        np.testing.assert_allclose(ds, 3.0)
        same = _downsample_volume_box_average(ds, target_d=128)
        np.testing.assert_array_equal(same, ds)

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

    def test_page_requires_volumes_eligible(self, flask_client) -> None:
        r = flask_client.get("/volume-viewer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "Volume viewer" in body or "CUDA GPU" in body

    def test_page_renders_when_volumes_eligible(
        self, flask_client_volumes_eligible
    ) -> None:
        r = flask_client_volumes_eligible.get("/volume-viewer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "Volume viewer" in body
        assert "vslice-canvas" in body
        assert "vslice-volume-picker" in body
        assert "vslice-volume-picker-rows" in body
        assert "btn-vslice-pan-up" in body
        assert "btn-vslice-zoom-in" in body
        assert "vslice-rotation-lock" in body
        assert "vslice-slice-contrast" in body
        assert "vslice-marker-glyphs-overlay" in body
        assert "vslice-status-slot" in body
        assert "vslice-controls-dock" in body
        assert "volume_slice_canvas.js" in body
        assert "volume_3d_utils.js" in body
        assert "volume_raycast_vtk.bundle.js" in body
        assert "vslice-vtk-container" in body
        assert "vslice-view-mode-3d" in body
        assert "vslice-iso-level" in body
        assert "Isosurface level" in body
        assert "vslice-volume-nav" in body
        assert "vslice-volume-nav-row" in body
        assert "btn-vslice-vol-prev" in body
        assert "btn-vslice-vol-next" in body
        assert "vslice-fix-grid-size" in body
        assert "cryo-vslice-view-mode-switch" in body
        assert "/api/volume_viewer/analyze_volumes" in body
        assert "/api/volume_viewer/analyze_markers" in body
        assert "/api/volume_viewer/analyze_volume" in body
        assert "/api/volume_viewer/decode" in body
        assert "Double-click" in body
        assert "Create indices file" not in body
        assert "cryo-explorer-save-btn" not in body
        assert "btn-vslice-decode-selection" not in body

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

    def test_landing_lists_volume_viewer(self, flask_client_volumes_eligible) -> None:
        r = flask_client_volumes_eligible.get("/")
        body = r.get_data(as_text=True)
        assert "Volume viewer" in body
        assert "/volume-viewer" in body


class TestVolumeSliceViewerBrowserSmoke:
    """Headless Chromium: analyze catalog, slice canvas, and 3D raycast."""

    def test_catalog_and_canvas_load(
        self, playwright_page, dashboard_volumes_eligible_live_url
    ) -> None:
        from tests.conftest import dashboard_smoke_volume_viewer

        out = dashboard_smoke_volume_viewer(
            playwright_page, dashboard_volumes_eligible_live_url
        )
        assert out is not None
        assert out["volume_picker_buttons"] >= 1
        assert out["canvas_visible"] is True
        assert (out.get("scatter_points") or 0) > 0

    def test_3d_mode_chrome_and_bundle(
        self, playwright_page, dashboard_volumes_eligible_live_url
    ) -> None:
        from tests.conftest import dashboard_smoke_volume_viewer_3d_chrome

        out = dashboard_smoke_volume_viewer_3d_chrome(
            playwright_page, dashboard_volumes_eligible_live_url
        )
        assert out is not None
        assert out["mode3d"] is True
        assert out["vtk_visible"] is True
        assert out["slice_canvas_hidden"] is True
        assert out["iso_visible"] is True
        assert out["slice_controls_hidden"] is True
        assert out["bundle_loaded"] is True
        assert "camera" in out["reset_label"].lower()

    def test_3d_raycast_render_when_webgl_available(
        self, playwright_page, dashboard_volumes_eligible_live_url
    ) -> None:
        from tests.conftest import dashboard_smoke_volume_viewer_3d_render

        out = dashboard_smoke_volume_viewer_3d_render(
            playwright_page, dashboard_volumes_eligible_live_url
        )
        assert out is not None
        if out.get("webgl_unavailable"):
            pytest.skip("WebGL unavailable for vtk.js after scattergl (headless CI)")
        assert out["vtk_canvas"] is True
        assert out["iso_level"] is not None
        assert 0 <= out["iso_level"] <= 100

    def test_3d_mode_roundtrip_to_slices(
        self, playwright_page, dashboard_volumes_eligible_live_url
    ) -> None:
        from tests.conftest import dashboard_smoke_volume_viewer_3d_roundtrip

        out = dashboard_smoke_volume_viewer_3d_roundtrip(
            playwright_page, dashboard_volumes_eligible_live_url
        )
        assert out is not None
        assert out["back_to_2d"] is True
        assert "slice" in out["reset_button_label_2d"].lower()

    def test_volume_picker_switches_selection(
        self, playwright_page, dashboard_volumes_eligible_live_url
    ) -> None:
        from tests.conftest import dashboard_smoke_volume_viewer_picker_switch

        out = dashboard_smoke_volume_viewer_picker_switch(
            playwright_page, dashboard_volumes_eligible_live_url
        )
        assert out is not None
        if not out.get("switched"):
            pytest.skip(out.get("reason", "not enough volumes to switch"))
        assert out["from_label"]
        assert out["to_label"]
        assert out["from_label"] != out["to_label"]

    def test_volume_nav_visible_in_3d_multi_select(
        self, playwright_page, dashboard_volumes_eligible_live_url
    ) -> None:
        from tests.conftest import dashboard_smoke_volume_viewer_nav_3d

        out = dashboard_smoke_volume_viewer_nav_3d(
            playwright_page, dashboard_volumes_eligible_live_url
        )
        assert out is not None
        if out.get("reason") == "fewer_than_two_volumes":
            pytest.skip("fewer than two volumes in catalog")
        assert out["nav_active"] is True
        assert out["nav_height"] > 0
        assert out["nav_label"] and out["nav_label"] != "—"
