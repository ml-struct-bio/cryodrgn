"""Tests for the dashboard volume slice viewer."""

from __future__ import annotations

import base64
import numpy as np
from cryodrgn.dashboard.volume_slice_viewer import (
    analyze_volume_by_id_payload,
    analyze_volume_markers_payload,
    analyze_volumes_batch_payload,
    analyze_volumes_catalog_payload,
    discover_analyze_volume_catalog,
    discover_analyze_volume_markers,
    volume_array_b64,
)


class TestVolumeSliceViewerPure:
    def test_volume_array_b64_roundtrip(self) -> None:
        vol = np.arange(8, dtype=np.float32).reshape(2, 2, 2)
        b64 = volume_array_b64(vol)
        raw = base64.standard_b64decode(b64)
        back = np.frombuffer(raw, dtype=np.float32).reshape(2, 2, 2)
        np.testing.assert_array_equal(back, vol)


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
        r = flask_client.get("/api/volume_slice_viewer/analyze_volumes")
        assert r.status_code == 200
        j = r.get_json()
        assert j["ok"] is True
        assert j["catalog"]
        assert j["markers"]
        assert len(j["markers"]) == len(j["catalog"])
        assert j.get("volumes") in (None, {})

    def test_analyze_volumes_api_fast_catalog(self, flask_client) -> None:
        r = flask_client.get(
            "/api/volume_slice_viewer/analyze_volumes?include_markers=0"
        )
        assert r.status_code == 200
        j = r.get_json()
        assert j["ok"] is True
        assert j["catalog"]
        assert j["markers"] == []

    def test_analyze_markers_api(self, flask_client) -> None:
        cat = flask_client.get(
            "/api/volume_slice_viewer/analyze_volumes?include_markers=0"
        ).get_json()
        r = flask_client.get("/api/volume_slice_viewer/analyze_markers")
        assert r.status_code == 200
        j = r.get_json()
        assert j["ok"] is True
        assert len(j["markers"]) == len(cat["catalog"])

    def test_analyze_volume_api(self, flask_client) -> None:
        cat = flask_client.get("/api/volume_slice_viewer/analyze_volumes").get_json()
        assert cat["catalog"]
        vol_id = cat["catalog"][0]["id"]
        r = flask_client.get(f"/api/volume_slice_viewer/analyze_volume?id={vol_id}")
        assert r.status_code == 200
        j = r.get_json()
        assert j["ok"] is True
        assert j["id"] == vol_id
        assert j["volume_b64"]

    def test_page_requires_volumes_eligible(self, flask_client) -> None:
        r = flask_client.get("/volume-slice-viewer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "Volume slice viewer" in body or "CUDA GPU" in body

    def test_page_renders_when_volumes_eligible(
        self, flask_client_volumes_eligible
    ) -> None:
        r = flask_client_volumes_eligible.get("/volume-slice-viewer")
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert "Volume slice viewer" in body
        assert "vslice-canvas" in body
        assert "vslice-volume-picker" in body
        assert "vslice-volume-picker-rows" in body
        assert "btn-vslice-pan-up" in body
        assert "btn-vslice-zoom-in" in body
        assert "vslice-slice-contrast" in body
        assert "vslice-marker-glyphs-overlay" in body
        assert "vslice-status-slot" in body
        assert "vslice-controls-dock" in body
        assert "volume_slice_canvas.js" in body
        assert "/api/volume_slice_viewer/analyze_volumes" in body
        assert "/api/volume_slice_viewer/analyze_markers" in body
        assert "/api/volume_slice_viewer/analyze_volume" in body
        assert "/api/volume_slice_viewer/decode" in body
        assert "Double-click" in body
        assert "Create indices file" not in body
        assert "cryo-explorer-save-btn" not in body
        assert "btn-vslice-decode-selection" not in body

    def test_decode_api_requires_row(self, flask_client_volumes_eligible) -> None:
        r = flask_client_volumes_eligible.post(
            "/api/volume_slice_viewer/decode",
            json={},
        )
        assert r.status_code == 400
        assert "row" in r.get_json().get("error", "").lower()

    def test_landing_lists_volume_slice_viewer(
        self, flask_client_volumes_eligible
    ) -> None:
        r = flask_client_volumes_eligible.get("/")
        body = r.get_data(as_text=True)
        assert "Volume slice viewer" in body
        assert "/volume-slice-viewer" in body
