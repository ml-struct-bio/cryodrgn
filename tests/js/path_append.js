/**
 * Appending random waypoints to an existing path, for waypoint mode, direct
 * tracing, and editable polylines.
 *
 * Run by ``TestTrajectoryScenarioScripts`` in ``tests/test_dashboard_trajectory.py``.
 */
"use strict";

const {
  VolumeState, Session, Path, assert
} = require("./_harness.js");

// Add-random must append waypoints, not splice into a stale densified polyline.
const appendIds = Array.from({ length: 10 }, (_, i) => "pc1:" + i);
const appendStore = {
  trajectoryMode: "manual",
  manualSelectedVolIds: appendIds,
  manualMarkersByVolId: appendIds.reduce((acc, id, i) => {
    acc[id] = { plot_row: i };
    return acc;
  }, {}),
  manualActiveCustomPlotRows: [],
  manualCustomPlotRows: [],
  anchorIndicesActive: appendIds.map((_, i) => i),
  anchorTrajXY: Array.from({ length: 19 }, (_, i) => [i, i]),
  trajPlotRows: Array.from({ length: 19 }, (_, i) => i),
  manualInterpTrajectoryCount: 19,
  manualVolumeSnapActive: true,
  manualParticleSnapPathActive: false,
  graphTraversalActive: false,
  currentPath: null
};
const appendSession = new Session({
  store: appendStore,
  pathHooks: { rowXYForPlotRow: row => [Number(row), Number(row)] },
  hooks: {
    appendPlotRows: rows => {
      rows.forEach(row => {
        row = Number(row);
        if (!Number.isFinite(row)) return;
        if (!appendStore.manualCustomPlotRows.includes(row)) {
          appendStore.manualCustomPlotRows.push(row);
        }
        if (!appendStore.manualActiveCustomPlotRows.includes(row)) {
          appendStore.manualActiveCustomPlotRows.push(row);
        }
      });
      return true;
    }
  },
  volumes: new VolumeState()
});
const appendRows = Array.from({ length: 10 }, (_, i) => 100 + i);
const appendResult = appendSession.appendPlotRows(appendRows);
assert(appendResult.ok, "append random waypoints ok");
assert(appendStore.manualActiveCustomPlotRows.length === 10, "ten random rows active");
assert(appendStore.anchorIndicesActive.length === 20, "10 catalog + 10 random waypoints");
assert(appendStore.trajPlotRows === null, "append rebuild clears stale densified samples");
assert(appendSession.path().sampleCount() === 20, "append sample count is waypoint count");

// Path-only anchors (live path without catalog/custom selection) must survive append.
const pathOnlyStore = {
  trajectoryMode: "manual",
  manualSelectedVolIds: [],
  manualMarkersByVolId: {},
  manualActiveCustomPlotRows: [],
  manualCustomPlotRows: [],
  anchorIndicesActive: [1, 5, 10, 20, 30],
  anchorTrajXY: [[1, 1], [5, 5], [10, 10], [20, 20], [30, 30]],
  trajPlotRows: null,
  manualInterpTrajectoryCount: 0,
  manualVolumeSnapActive: false,
  manualParticleSnapPathActive: false,
  graphTraversalActive: false,
  currentPath: null
};
const pathOnlySession = new Session({
  store: pathOnlyStore,
  pathHooks: { rowXYForPlotRow: row => [Number(row), Number(row)] },
  hooks: {
    appendPlotRows: rows => {
      rows.forEach(row => {
        row = Number(row);
        if (!Number.isFinite(row)) return;
        if (!pathOnlyStore.manualCustomPlotRows.includes(row)) {
          pathOnlyStore.manualCustomPlotRows.push(row);
        }
        if (!pathOnlyStore.manualActiveCustomPlotRows.includes(row)) {
          pathOnlyStore.manualActiveCustomPlotRows.push(row);
        }
      });
      return true;
    }
  },
  volumes: new VolumeState()
});
const pathOnlyAppend = Array.from({ length: 10 }, (_, i) => 100 + i);
// Simulate ensureExistingWaypointsTrackedBeforeRandomAppend (template layer).
(pathOnlyStore.anchorIndicesActive || []).forEach(row => {
  row = Number(row);
  if (!Number.isFinite(row)) return;
  if (!pathOnlyStore.manualCustomPlotRows.includes(row)) {
    pathOnlyStore.manualCustomPlotRows.push(row);
  }
  if (!pathOnlyStore.manualActiveCustomPlotRows.includes(row)) {
    pathOnlyStore.manualActiveCustomPlotRows.push(row);
  }
});
const pathOnlyResult = pathOnlySession.appendPlotRows(pathOnlyAppend);
assert(pathOnlyResult.ok, "path-only append ok");
assert(
  pathOnlyStore.anchorIndicesActive.length === 15,
  "path-only append keeps prior anchors + random"
);

// Direct-trace Add-random must use the same Session append + rebuild path as waypoints.
const directAppendIds = Array.from({ length: 10 }, (_, i) => "pc1:" + i);
const directAppendStore = {
  trajectoryMode: "direct",
  manualSelectedVolIds: directAppendIds,
  manualMarkersByVolId: directAppendIds.reduce((acc, id, i) => {
    acc[id] = { plot_row: i };
    return acc;
  }, {}),
  manualActiveCustomPlotRows: [],
  manualCustomPlotRows: [],
  anchorIndicesActive: directAppendIds.map((_, i) => i),
  anchorTrajXY: Array.from({ length: 19 }, (_, i) => [i, i]),
  trajPlotRows: Array.from({ length: 19 }, (_, i) => i),
  editableTrajXY: Array.from({ length: 19 }, (_, i) => [i, i]),
  manualInterpTrajectoryCount: 19,
  startXY: [0, 0],
  endXY: [9, 9],
  directEndpointVolumeIds: [],
  currentPath: null
};
const directAppendSession = new Session({
  store: directAppendStore,
  pathHooks: { rowXYForPlotRow: row => [Number(row), Number(row)] },
  hooks: {
    appendPlotRows: rows => {
      rows.forEach(row => {
        row = Number(row);
        if (!Number.isFinite(row)) return;
        if (!directAppendStore.manualCustomPlotRows.includes(row)) {
          directAppendStore.manualCustomPlotRows.push(row);
        }
        if (!directAppendStore.manualActiveCustomPlotRows.includes(row)) {
          directAppendStore.manualActiveCustomPlotRows.push(row);
        }
      });
      return true;
    }
  },
  volumes: new VolumeState()
});
const directAppendRows = Array.from({ length: 10 }, (_, i) => 200 + i);
const directAppendResult = directAppendSession.appendPlotRows(directAppendRows);
assert(directAppendResult.ok, "direct-trace append random waypoints ok");
assert(directAppendStore.manualActiveCustomPlotRows.length === 10, "direct-trace ten random rows active");
assert(directAppendStore.anchorIndicesActive.length === 29, "direct-trace append extends anchor rows");
assert(directAppendStore.trajPlotRows && directAppendStore.trajPlotRows.length === 29,
  "direct-trace append extends densified samples");
assert(directAppendStore.editableTrajXY && directAppendStore.editableTrajXY.length === 29,
  "direct-trace append extends editable polyline");
assert(directAppendSession.path().kind === "direct", "direct-trace append keeps direct path kind");
assert(directAppendSession.path().sampleCount() === 29, "direct-trace append sample count grows");

// Direct-trace editable polyline without pre-materialized plot rows must append.
const directEditableOnlyStore = {
  trajectoryMode: "direct",
  manualSelectedVolIds: [],
  manualMarkersByVolId: {},
  manualActiveCustomPlotRows: [],
  manualCustomPlotRows: [],
  anchorIndicesActive: null,
  anchorTrajXY: null,
  trajPlotRows: null,
  editableTrajXY: Array.from({ length: 10 }, (_, i) => [i, i * 0.5]),
  manualInterpTrajectoryCount: 10,
  startXY: [0, 0],
  endXY: [9, 4.5],
  directEndpointVolumeIds: [],
  currentPath: null
};
const directEditableOnlySession = new Session({
  store: directEditableOnlyStore,
  pathHooks: {
    rowXYForPlotRow: row => [Number(row), Number(row) * 0.5],
    nearestPlotRowForXY: x => Math.round(Number(x)),
    plotRowsForPathXY: xy => xy.map((_, i) => i)
  },
  hooks: {
    appendPlotRows: rows => {
      rows.forEach(row => {
        row = Number(row);
        if (!Number.isFinite(row)) return;
        if (!directEditableOnlyStore.manualCustomPlotRows.includes(row)) {
          directEditableOnlyStore.manualCustomPlotRows.push(row);
        }
        if (!directEditableOnlyStore.manualActiveCustomPlotRows.includes(row)) {
          directEditableOnlyStore.manualActiveCustomPlotRows.push(row);
        }
      });
      return true;
    }
  },
  volumes: new VolumeState()
});
const editableOnlyAppend = Array.from({ length: 5 }, (_, i) => 500 + i);
const editableOnlyResult = directEditableOnlySession.appendPlotRows(editableOnlyAppend);
assert(editableOnlyResult.ok, "direct editable-only append ok");
assert(
  directEditableOnlyStore.editableTrajXY
    && directEditableOnlyStore.editableTrajXY.length === 15,
  "direct editable-only append extends latent polyline"
);
assert(
  directEditableOnlyStore.trajPlotRows
    && directEditableOnlyStore.trajPlotRows.length === 15,
  "direct editable-only append materializes plot rows"
);

// Direct-trace with no live geometry still rebuilds from catalog + custom selection.
const directRebuildStore = {
  trajectoryMode: "direct",
  manualSelectedVolIds: directAppendIds,
  manualMarkersByVolId: directAppendIds.reduce((acc, id, i) => {
    acc[id] = { plot_row: i };
    return acc;
  }, {}),
  manualActiveCustomPlotRows: [],
  manualCustomPlotRows: [],
  anchorIndicesActive: null,
  anchorTrajXY: null,
  trajPlotRows: null,
  editableTrajXY: null,
  manualInterpTrajectoryCount: 0,
  startXY: null,
  endXY: null,
  directEndpointVolumeIds: [],
  currentPath: null
};
const directRebuildSession = new Session({
  store: directRebuildStore,
  pathHooks: { rowXYForPlotRow: row => [Number(row), Number(row)] },
  hooks: {
    appendPlotRows: rows => {
      rows.forEach(row => {
        row = Number(row);
        if (!Number.isFinite(row)) return;
        if (!directRebuildStore.manualCustomPlotRows.includes(row)) {
          directRebuildStore.manualCustomPlotRows.push(row);
        }
        if (!directRebuildStore.manualActiveCustomPlotRows.includes(row)) {
          directRebuildStore.manualActiveCustomPlotRows.push(row);
        }
      });
      return true;
    }
  },
  volumes: new VolumeState()
});
const directRebuildRows = Array.from({ length: 10 }, (_, i) => 300 + i);
const directRebuildResult = directRebuildSession.appendPlotRows(directRebuildRows);
assert(directRebuildResult.ok, "direct-trace rebuild append ok");
assert(directRebuildStore.anchorIndicesActive.length === 20,
  "direct-trace rebuild uses catalog + random waypoints");
assert(directRebuildSession.path().sampleCount() === 20,
  "direct-trace rebuild sample count is waypoint count");

console.log("path_append: ok");
