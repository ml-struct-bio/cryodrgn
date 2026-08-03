/**
 * Session-level path mutations: reverse, render batching, and the rebuild that
 * keeps volume readiness aligned through a mutation rather than after it.
 *
 * Run by ``TestTrajectoryScenarioScripts`` in ``tests/test_dashboard_trajectory.py``.
 */
"use strict";

const {
  VolumeState, Session, Mutations, Path, assert
} = require("./_harness.js");

// Session + reverse mutation with waypoint store
const store = {
  trajectoryMode: "manual",
  manualSelectedVolIds: ["a", "b", "c"],
  anchorIndicesActive: [1, 2, 3],
  anchorTrajXY: [[0, 0], [1, 1], [2, 2]],
  trajPlotRows: [1, 2, 3],
  manualInterpTrajectoryCount: 0,
  manualVolumeSnapActive: false,
  manualParticleSnapPathActive: false,
  graphTraversalActive: false,
  manualMarkersByVolId: {},
  manualActiveCustomPlotRows: [],
  manualCustomPlotRows: [],
  anchorKmeansVolumeIds: [],
  currentPath: null
};
const session = new Session({
  store,
  pathHooks: {},
  volumes: new VolumeState()
});
session.volumes().replaceSlots(
  ["1", "2", "3"],
  [{ volume_b64: "A" }, { volume_b64: "B" }, null],
  ["ia", "ib", null]
);
const path = Path.ensureCurrentPath(store, {});
assert(path && path.kind === "waypoint", "waypoint path");
const decodeBefore = session.volumes().missingDecodeIndices().length;
const renderBefore = session.volumes().missingRenderIndices().length;
const result = session.dispatch(Mutations.reverse(), { silent: true });
assert(result.ok, "reverse ok");
assert(store.manualSelectedVolIds.join(",") === "a,b,c",
  "catalog selection is not reversed (avoids remapping slot ids)");
assert(store.anchorIndicesActive.join(",") === "3,2,1", "path anchors reversed");
assert(store.trajPlotRows.join(",") === "3,2,1", "path sample rows reversed");
assert(session.volumes().ids().join(",") === "3,2,1",
  "volume view is a pure reverse of prior ids");
assert(session.volumes().missingDecodeIndices().length === decodeBefore,
  "decode debt unchanged by reverse");
assert(session.volumes().missingRenderIndices().length === renderBefore,
  "render debt unchanged by reverse");
assert(session.renderDebtCount() === session.volumes().renderDebtCount(),
  "Session.renderDebtCount falls back to VolumeState");
assert(session.inactiveVolumeTickCount() === session.renderDebtCount(),
  "Session inactiveVolumeTickCount ≡ renderDebtCount");
assert(session.decodeDebtCount() === session.volumes().decodeDebtCount(),
  "Session.decodeDebtCount falls back to VolumeState via inactive ∩ undecoded");
var sessionDebts = session.decodeRenderDebts();
assert(sessionDebts.decode === session.decodeDebtCount()
    && sessionDebts.render === session.renderDebtCount()
    && sessionDebts.decode <= sessionDebts.render,
  "Session.decodeRenderDebts decode ≤ render by construction");
session.hooks.renderDebtCount = function() { return 7; };
assert(session.renderDebtCount() === 7,
  "Session.renderDebtCount prefers display hook (slider invariant)");
delete session.hooks.renderDebtCount;
session.hooks.inactiveTickIndices = function() { return [0, 2]; };
session.hooks.isSlotDecoded = function(i) { return i === 2; };
var hookedDebts = session.decodeRenderDebts();
assert(hookedDebts.render === 2 && hookedDebts.decode === 1
    && hookedDebts.decodeIndices.join(",") === "0",
  "Session hooks: decode is inactive ticks minus decoded slots");
assert(hookedDebts.decode <= hookedDebts.render,
  "hooked decodeRenderDebts cannot exceed inactive ticks");
delete session.hooks.inactiveTickIndices;
delete session.hooks.isSlotDecoded;
// renderVolumes batch indices ≡ inactive ticks (progress / button one number).
session.hooks.inactiveTickIndices = function() { return [1, 3, 4]; };
var capturedRenderOpts = null;
session._pipeline = {
  render: function(opts) { capturedRenderOpts = opts || {}; return Promise.resolve({ ok: true }); }
};
session.renderVolumes({});
assert(capturedRenderOpts && capturedRenderOpts.indices
    && capturedRenderOpts.indices.join(",") === "1,3,4",
  "Session.renderVolumes uses inactiveVolumeTickIndices as the batch");
session.renderVolumes({ forceAll: true });
assert(capturedRenderOpts && capturedRenderOpts.forceAll === true
    && capturedRenderOpts.indices == null,
  "forceAll leaves indices unset for pipeline expansion");
delete session.hooks.inactiveTickIndices;
session._pipeline = null;

// DirectTracePath delegates endpoint linearization to the shared helper.
const line = Path.linearXYFromEndpoints([0, 0], [2, 4], 3);
assert(JSON.stringify(line) === JSON.stringify([[0, 0], [1, 2], [2, 4]]), "linear helper");
const directStore = {
  trajectoryMode: "direct",
  anchorIndicesActive: [1, 2],
  anchorTrajXY: [[9, 9], [10, 10]],
  editableTrajXY: null,
  trajPlotRows: [1, 2],
  manualInterpTrajectoryCount: 0,
  startXY: [0, 0],
  endXY: [2, 4],
  directEndpointVolumeIds: [],
  currentPath: null
};
const directPath = Path.ensureCurrentPath(directStore, {});
assert(directPath.buildLinearPath(3), "direct build linear path");
assert(directStore.anchorIndicesActive === null, "direct build clears anchors");
assert(directStore.editableTrajXY.length === 3, "direct build writes editable samples");
assert(directStore.manualInterpTrajectoryCount === 3, "direct build tracks sample count");

// Session rebuild aligns volume readiness through the mutation, not a template fallback.
const rebuildStore = {
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
    v1: { plot_row: 10 },
    v2: { plot_row: 20 }
  },
  manualActiveCustomPlotRows: [],
  manualCustomPlotRows: [],
  anchorKmeansVolumeIds: [],
  currentPath: null
};
const rebuildVolumes = new VolumeState();
rebuildVolumes.replaceSlots(["v2", "v1"], [{ volume_b64: "B" }, { volume_b64: "A" }], []);
const rebuildSession = new Session({
  store: rebuildStore,
  pathHooks: { rowXYForPlotRow: row => [Number(row), Number(row) + 1] },
  hooks: { slotIdsForPath: () => rebuildStore.manualSelectedVolIds.slice() },
  volumes: rebuildVolumes
});
const rebuild = rebuildSession.rebuildFromSelection({ force: true });
assert(rebuild.ok && rebuild.rows.join(",") === "10,20", "session rebuild rows");
assert(rebuildVolumes.slotIdAt(0) === "v1", "session rebuild aligns first id");
assert(rebuildVolumes.volumes()[0].volume_b64 === "A", "session rebuild preserves first volume");

// Readiness follows slot identity through the reverse, rather than staying put.
assert(session.volumes().isDecoded(2) && session.volumes().isRendered(2),
  "former first slot readiness follows id to the end");
assert(!session.volumes().isDecoded(0) && !session.volumes().isRendered(0),
  "former last slot stays inactive at the start");

console.log("session_mutations: ok");
