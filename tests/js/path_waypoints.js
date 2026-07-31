/**
 * Waypoint selection and visit order: one anchor per selected volume, marker
 * coordinates, and tours that interleave catalogue and custom rows.
 *
 * Run by ``TestTrajectoryScenarioScripts`` in ``tests/test_dashboard_trajectory.py``.
 */
"use strict";

const {
  VolumeState, VolumeStateUtils, Session, Mutations, Path, Display,
  DirectTraceUiState, assert, elStub
} = require("./_harness.js");

// Waypoint rebuild keeps one anchor per selected volume, even when two volumes
// snap to the same scatter particle.
const duplicateStore = {
  trajectoryMode: "manual",
  manualSelectedVolIds: ["v1", "v2"],
  anchorIndicesActive: null,
  anchorTrajXY: null,
  trajPlotRows: null,
  manualInterpTrajectoryCount: 0,
  manualVolumeSnapActive: false,
  manualParticleSnapPathActive: false,
  graphTraversalActive: false,
  manualMarkersByVolId: {
    v1: { plot_row: 7 },
    v2: { plot_row: 7 }
  },
  manualActiveCustomPlotRows: [],
  manualCustomPlotRows: [],
  anchorKmeansVolumeIds: [],
  currentPath: null
};
const duplicatePath = Path.ensureCurrentPath(duplicateStore, {
  rowXYForPlotRow: row => [Number(row), Number(row) + 0.5]
});
const duplicateResult = duplicatePath.rebuildFromSelection({ force: true });
assert(duplicateResult.rows.join(",") === "7,7", "duplicate selected rows preserved");
assert(duplicateStore.anchorTrajXY.length === 2, "duplicate rows keep duplicate XY samples");
assert(duplicatePath.sampleCount() === 2, "duplicate rows count as two samples");

duplicateStore.manualActiveCustomPlotRows = [7, 8];
const customResult = duplicatePath.rebuildFromSelection({ force: true });
assert(customResult.rows.join(",") === "7,7,8", "custom row appends only when distinct");
assert(duplicateStore.anchorTrajXY.length === 3, "distinct custom row adds XY sample");

// plot_row alone is enough for densify APIs — do not require marker.xy / rowXY.
const plotRowOnlyStore = {
  trajectoryMode: "manual",
  manualSelectedVolIds: ["v1", "v2", "v3"],
  anchorIndicesActive: null,
  anchorTrajXY: null,
  trajPlotRows: null,
  manualInterpTrajectoryCount: 0,
  manualVolumeSnapActive: false,
  manualParticleSnapPathActive: false,
  graphTraversalActive: false,
  manualMarkersByVolId: {
    v1: { plot_row: 1 },
    v2: { plot_row: 2 },
    v3: { plot_row: 3 }
  },
  manualActiveCustomPlotRows: [],
  manualCustomPlotRows: [],
  anchorKmeansVolumeIds: [],
  currentPath: null
};
const plotRowOnlyPath = Path.ensureCurrentPath(plotRowOnlyStore, {});
const plotRowOnlySel = plotRowOnlyPath.anchorRowsFromSelection({ includeCustom: false });
assert(plotRowOnlySel.rows.join(",") === "1,2,3", "plot_row-only markers resolve without xy");
assert(plotRowOnlySel.xy === null, "xy omitted when coordinates are unavailable");
plotRowOnlyPath.armInterpolation("volume");
assert(plotRowOnlyStore.manualVolumeSnapActive === true, "direct-line arm sets volume snap");
assert(plotRowOnlyStore.manualParticleSnapPathActive === false, "direct-line does not arm particle snap");
assert(plotRowOnlyStore.graphTraversalActive === false, "direct-line arm clears graph");

// Latent-only markers (xy, no plot_row) must not invent negative dataset indices.
const latentOnlyStore = {
  trajectoryMode: "manual",
  manualSelectedVolIds: ["pc:a", "pc:b"],
  manualMarkersByVolId: {
    "pc:a": { xy: [0, 0] },
    "pc:b": { xy: [1, 1] }
  },
  manualActiveCustomPlotRows: [],
  manualCustomPlotRows: [],
  currentPath: null
};
const latentOnlyPath = Path.ensureCurrentPath(latentOnlyStore, {});
const latentOnlySel = latentOnlyPath.anchorRowsFromSelection({ includeCustom: false });
assert(latentOnlySel.rows.length === 0, "xy-only markers are not fake plot rows");

// Visit-order waypoints must include Other / random indices with catalog vols.
const visitStore = {
  trajectoryMode: "manual",
  manualSelectedVolIds: ["pc2:0", "pc2:1", "pc2:2"],
  manualMarkersByVolId: {
    "pc2:0": { plot_row: 10 },
    "pc2:1": { plot_row: 20 },
    "pc2:2": { plot_row: 30 }
  },
  manualActiveCustomPlotRows: [101, 102, 103],
  manualCustomPlotRows: [101, 102, 103],
  currentPath: null
};
const visitPath = Path.ensureCurrentPath(visitStore, {
  rowXYForPlotRow: row => [Number(row), Number(row)]
});
const visitSel = visitPath.anchorRowsFromSelection({ includeCustom: true });
assert(
  visitSel.rows.join(",") === "10,20,30,101,102,103",
  "visit-order selection includes PC + random indices"
);
const catalogOnly = visitPath.anchorRowsFromSelection({ includeCustom: false });
assert(catalogOnly.rows.join(",") === "10,20,30", "catalog-only excludes random indices");
visitPath.rebuildFromSelection({ force: true });
assert(
  visitStore.anchorIndicesActive.join(",") === "10,20,30,101,102,103",
  "rebuildFromSelection keeps random indices on the path"
);

// Exact / Inexact must keep an interleaved catalog + Other tour (not regroup
// to catalog-then-custom).
const visitOrderRows = [101, 10, 102, 20, 103, 30];
visitStore.manualSelectedVolIds = ["pc2:0", "pc2:1", "pc2:2"];
visitStore.manualActiveCustomPlotRows = [101, 102, 103];
const visitOrderSession = new Session({
  store: visitStore,
  pathHooks: { rowXYForPlotRow: row => [Number(row), Number(row)] },
  hooks: {
    reorderSelectionFromRows: rows => {
      const volsByRow = { 10: ["pc2:0"], 20: ["pc2:1"], 30: ["pc2:2"] };
      const nextVols = [];
      const nextCustom = [];
      rows.forEach(row => {
        row = Number(row);
        if (volsByRow[row] && volsByRow[row].length) nextVols.push(volsByRow[row].shift());
        if ([101, 102, 103].includes(row)) nextCustom.push(row);
      });
      visitStore.manualSelectedVolIds = nextVols;
      visitStore.manualActiveCustomPlotRows = nextCustom;
    }
  },
  volumes: new VolumeState()
});
const visitOrderResult = visitOrderSession.setVisitOrderPath("exact", visitOrderRows);
assert(visitOrderResult.ok, "visit-order mutation ok");
assert(
  visitStore.anchorIndicesActive.join(",") === visitOrderRows.join(","),
  "visit-order keeps interleaved PC + random tour"
);
const afterVisitSel = visitOrderSession.path().anchorRowsFromSelection({ includeCustom: true });
assert(
  afterVisitSel.rows.join(",") === visitOrderRows.join(","),
  "selection helper preserves interleaved visit order"
);
visitOrderSession.rebuildFromSelection({ force: true });
assert(
  visitStore.anchorIndicesActive.join(",") === visitOrderRows.join(","),
  "rebuild after visit-order does not regroup catalog-then-custom"
);

// Visit-order must never drop decode/render identity for Other (non-catalog)
// slots. Session.syncVolumesToPath's fallback (no full per-slot ``ids`` layout
// available — e.g. a session bootstrapped without slotIdsForPath /
// volumePathLayout hooks) used to fall back to volumes.alignToPath with only
// the compact catalog id list, silently collapsing an already-aligned mixed
// catalog+Other path down to catalog-only anchors and orphaning every Other
// slot's decoded/rendered media at its old position.
const visitDecodeStore = {
  trajectoryMode: "manual",
  manualSelectedVolIds: ["pc3:0", "pc3:1", "pc3:2"],
  manualMarkersByVolId: {
    "pc3:0": { plot_row: 40 },
    "pc3:1": { plot_row: 50 },
    "pc3:2": { plot_row: 60 }
  },
  manualActiveCustomPlotRows: [201, 202, 203],
  manualCustomPlotRows: [201, 202, 203],
  currentPath: null
};
const visitDecodeVolumes = new VolumeState();
const visitDecodeInitialIds = ["pc3:0", "custom:201", "pc3:1", "custom:202", "pc3:2", "custom:203"];
visitDecodeVolumes.replaceSlots(
  visitDecodeInitialIds,
  visitDecodeInitialIds.map(id => ({ volume_b64: "V-" + id })),
  visitDecodeInitialIds.map(id => "img-" + id)
);
const visitDecodeSession = new Session({
  store: visitDecodeStore,
  pathHooks: { rowXYForPlotRow: row => [Number(row), Number(row)] },
  hooks: {
    // Deliberately no slotIdsForPath / volumePathLayout hooks here — this is
    // the exact configuration that exposes the fallback in syncVolumesToPath.
    reorderSelectionFromRows: rows => {
      const volsByRow = { 40: ["pc3:0"], 50: ["pc3:1"], 60: ["pc3:2"] };
      const nextVols = [];
      const nextCustom = [];
      rows.forEach(row => {
        row = Number(row);
        if (volsByRow[row] && volsByRow[row].length) nextVols.push(volsByRow[row].shift());
        if ([201, 202, 203].includes(row)) nextCustom.push(row);
      });
      visitDecodeStore.manualSelectedVolIds = nextVols;
      visitDecodeStore.manualActiveCustomPlotRows = nextCustom;
    }
  },
  volumes: visitDecodeVolumes
});
assert(visitDecodeVolumes.decodedCount() === 6 && visitDecodeVolumes.renderedCount() === 6,
  "visit-order regression fixture starts fully decoded and rendered");
const visitDecodeNewOrder = [201, 40, 202, 50, 203, 60];
const visitDecodeResult = visitDecodeSession.setVisitOrderPath("exact", visitDecodeNewOrder);
assert(visitDecodeResult.ok, "visit-order dispatch ok without layout hooks");
assert(visitDecodeVolumes.decodedCount() === 6,
  "visit-order without layout hooks must not drop decode readiness for Other slots");
assert(visitDecodeVolumes.renderedCount() === 6,
  "visit-order without layout hooks must not drop render readiness for Other slots");
var visitDecodeDebts = visitDecodeVolumes.decodeRenderDebts();
assert(visitDecodeDebts.decode === 0 && visitDecodeDebts.render === 0,
  "visit-order alone must never reintroduce Decode/Render debt on an already-ready path");

// Shared plot_row catalog volumes must keep distinct marker XY when rebuilding
// (Reset after visit-order must not hide point 3 under point 4).
const dupStore = {
  trajectoryMode: "manual",
  manualSelectedVolIds: ["pc1:0", "pc1:1", "pc1:2", "pc1:3"],
  manualMarkersByVolId: {
    "pc1:0": { plot_row: 10, xy: [0, 0] },
    "pc1:1": { plot_row: 20, xy: [1, 0] },
    // Two PC samples share a nearest-particle row but keep distinct latent XY.
    "pc1:2": { plot_row: 30, xy: [2, 0] },
    "pc1:3": { plot_row: 30, xy: [3, 0] }
  },
  manualActiveCustomPlotRows: [],
  manualCustomPlotRows: [],
  anchorIndicesActive: [30, 10, 30, 20],
  currentPath: null
};
const dupPath = Path.ensureCurrentPath(dupStore, {
  rowXYForPlotRow: row => [Number(row), Number(row)]
});
const dupSel = dupPath.anchorRowsFromSelection({ includeCustom: true });
assert(dupSel.rows.join(",") === "30,10,30,20", "prefer active order with duplicate plot rows");
assert(dupSel.xy && dupSel.xy.length === 4, "duplicate plot_row selection has xy");
assert(
  dupSel.xy[0][0] === 2 && dupSel.xy[2][0] === 3,
  "duplicate plot_row slots keep distinct catalog marker XY"
);
dupPath.rebuildFromSelection({ force: true, preferActiveOrder: false });
assert(
  dupStore.anchorIndicesActive.join(",") === "10,20,30,30",
  "preferActiveOrder:false restores catalog volume order"
);
assert(
  dupStore.anchorTrajXY[2][0] === 2 && dupStore.anchorTrajXY[3][0] === 3,
  "catalog rebuild keeps distinct XY for shared plot_row volumes"
);

console.log("path_waypoints: ok");
