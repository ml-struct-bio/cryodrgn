/**
 * Keeping volume slots aligned with a mixed catalogue/custom path across
 * densify, deselection, append, restore, and reverse.
 *
 * Run by ``TestTrajectoryScenarioScripts`` in ``tests/test_dashboard_trajectory.py``.
 */
"use strict";

const {
  VolumeState, VolumeStateUtils, Session, Mutations, Path, Display,
  DirectTraceUiState, assert, elStub
} = require("./_harness.js");

// Mixed PC + Other waypoints: syncVolumesToPath must keep 20 full-path slots,
// not invent a densified 10-anchor layout that leaves half the ticks inactive.
const mixedIds = Array.from({ length: 10 }, (_, i) => "pc1:" + i)
  .concat(Array.from({ length: 10 }, (_, i) => "custom:" + (100 + i)));
const mixedVolumes = new VolumeState();
mixedVolumes.replaceSlots(
  mixedIds.slice(0, 10),
  mixedIds.slice(0, 10).map((_, i) => ({ volume_b64: "V" + i })),
  []
);
const mixedSession = new Session({
  store: {
    trajectoryMode: "manual",
    manualSelectedVolIds: mixedIds.slice(0, 10),
    manualMarkersByVolId: mixedIds.slice(0, 10).reduce((acc, id, i) => {
      acc[id] = { plot_row: i };
      return acc;
    }, {}),
    manualActiveCustomPlotRows: Array.from({ length: 10 }, (_, i) => 100 + i),
    manualCustomPlotRows: Array.from({ length: 10 }, (_, i) => 100 + i),
    anchorIndicesActive: Array.from({ length: 20 }, (_, i) => (i < 10 ? i : 100 + (i - 10))),
    currentPath: null
  },
  pathHooks: { rowXYForPlotRow: row => [Number(row), Number(row)] },
  hooks: {
    volumePathLayout: () => ({
      pathN: 20,
      nAnchors: 20,
      nInterp: 0,
      ids: mixedIds.slice(),
      compactIds: mixedIds.slice()
    }),
    slotIdsForPath: () => mixedIds.slice()
  },
  volumes: mixedVolumes
});
const mixedSync = mixedSession.syncVolumesToPath({ pathN: 20 });
assert(mixedSync.ok && mixedVolumes.slotCount() === 20, "mixed path sync keeps 20 slots");
assert(mixedVolumes.slotIdAt(0) === "pc1:0", "mixed path keeps catalog slot 0");
assert(mixedVolumes.slotIdAt(10) === "custom:100", "mixed path keeps custom slot 10");
assert(mixedVolumes.isDecoded(0), "catalog decode preserved after mixed sync");
assert(!mixedVolumes.isDecoded(10), "custom slot starts undecoded");

const mixedPath = mixedSession.path();
const mixedRows = Array.from({ length: 20 }, (_, i) => (i < 10 ? i : 100 + (i - 10)));
const inflatedXY = Array.from({ length: 39 }, (_, i) => [i, -i]);
mixedPath.applyCoordsPayload({
  anchor_indices: mixedRows.slice(),
  traj_xy: inflatedXY,
  z_traj: inflatedXY
}, { preserveWaypointGeometry: true });
assert(mixedSession.sampleCount() === 20, "undensified decode apply keeps 20 waypoint samples");
assert(!mixedSession.path().hasInteriorSamples, "undensified decode apply does not create interiors");

mixedPath.applyCoordsPayload({
  anchor_indices: mixedRows.slice(),
  traj_xy: inflatedXY,
  z_traj: inflatedXY
});
assert(mixedSession.path().sampleCount() === 39, "densify-armed decode apply may grow samples");
assert(mixedSession.path().hasInteriorSamples, "densified decode apply records interiors");

// After Decode/Render of an interleaved PC+Other path, deselecting Other must
// rematch ChimeraX frames by stable slot id (not by index).
const rematchIds = [];
const rematchVols = [];
const rematchImgs = [];
for (let i = 0; i < 10; i++) {
  rematchIds.push("pc1:" + i);
  rematchVols.push({ volume_b64: "PC" + i });
  rematchImgs.push("img-pc-" + i);
  rematchIds.push("custom:" + (100 + i));
  rematchVols.push({ volume_b64: "CU" + i });
  rematchImgs.push("img-cu-" + i);
}
const rematchState = new VolumeState();
rematchState.replaceSlots(rematchIds, rematchVols, rematchImgs);
rematchState.setGenerated(true, "cache-rematch");
const catalogOnlyIds = Array.from({ length: 10 }, (_, i) => "pc1:" + i);
rematchState.alignToIds(catalogOnlyIds);
assert(rematchState.slotCount() === 10, "deselect Other shrinks to catalog slots");
assert(rematchState.slotIdAt(0) === "pc1:0", "catalog slot 0 kept after Other deselection");
assert(rematchState.images()[0] === "img-pc-0", "catalog image 0 rematched by id");
assert(rematchState.images()[3] === "img-pc-3", "catalog image 3 rematched by id");
assert(rematchState.volumes()[3].volume_b64 === "PC3", "catalog volume 3 rematched by id");
assert(!rematchState.images().some((img) => String(img || "").indexOf("img-cu-") === 0),
  "custom images dropped after Other deselection");

// Visit-order catalog-only rematch: selection order differs from path order.
const visitOrderIds = ["pc1:2", "pc1:0", "pc1:1"];
const visitOrderState = new VolumeState();
visitOrderState.replaceSlots(
  ["pc1:0", "pc1:1", "pc1:2"],
  [{ volume_b64: "A" }, { volume_b64: "B" }, { volume_b64: "C" }],
  ["imgA", "imgB", "imgC"]
);
visitOrderState.alignToIds(visitOrderIds);
assert(visitOrderState.images().join(",") === "imgC,imgA,imgB",
  "visit-order realign rematches images by id");
assert(visitOrderState.volumes().map(v => v.volume_b64).join(",") === "C,A,B",
  "visit-order realign rematches volumes by id");

// Empty incoming write must not be modeled as wiping renders: alignToIds keeps
// frames when shrinking an interleaved path to catalog-only ids.
const preserveState = new VolumeState();
preserveState.replaceSlots(rematchIds, rematchVols, rematchImgs);
preserveState.setGenerated(true, "cache-preserve");
const beforeRendered = preserveState.renderedCount();
assert(beforeRendered === 20, "interleaved path starts fully rendered");
preserveState.alignToIds(catalogOnlyIds);
assert(preserveState.renderedCount() === 10, "catalog-only shrink keeps all PC renders");
assert(preserveState.missingRenderIndices().length === 0, "no render debt after Other deselection");

// Path-slot authority after Other deselection: the retained PC rows must follow
// the live path order, not manualSelectedVolIds / picker order.
const afterOtherStore = {
  trajectoryMode: "manual",
  manualSelectedVolIds: ["pc1:0", "pc1:1", "pc1:2"],
  manualMarkersByVolId: {
    "pc1:0": { plot_row: 0 },
    "pc1:1": { plot_row: 1 },
    "pc1:2": { plot_row: 2 }
  },
  manualActiveCustomPlotRows: [],
  manualCustomPlotRows: [100, 101],
  // This is the filtered path after an Exact / Inexact interleaving such as
  // [100, 2, 101, 0, 1] followed by deselecting Other.
  anchorIndicesActive: [2, 0, 1],
  anchorTrajXY: [[2, 2], [0, 0], [1, 1]],
  trajPlotRows: null,
  currentPath: null
};
const afterOtherSession = new Session({
  store: afterOtherStore,
  pathHooks: { rowXYForPlotRow: row => [Number(row), Number(row)] },
  volumes: new VolumeState()
});
assert(
  afterOtherSession.waypointRows().join(",") === "2,0,1",
  "session waypointRows prefer filtered live path order after Other deselection"
);
const afterOtherIds = afterOtherSession.waypointRows().map(row => "pc1:" + row);
const afterOtherVolumes = new VolumeState();
afterOtherVolumes.replaceSlots(
  ["custom:100", "pc1:2", "custom:101", "pc1:0", "pc1:1"],
  [
    { volume_b64: "CU100" },
    { volume_b64: "PC2" },
    { volume_b64: "CU101" },
    { volume_b64: "PC0" },
    { volume_b64: "PC1" }
  ],
  ["imgCU100", "imgPC2", "imgCU101", "imgPC0", "imgPC1"]
);
afterOtherVolumes.alignToIds(afterOtherIds);
assert(afterOtherVolumes.ids().join(",") === "pc1:2,pc1:0,pc1:1",
  "slider ids follow path-derived PC order after Other deselection");
assert(afterOtherVolumes.images().join(",") === "imgPC2,imgPC0,imgPC1",
  "PC images remain attached to path-derived ids after Other deselection");
assert(afterOtherVolumes.missingRenderIndices().length === 0,
  "path-derived catalog-only slider has no false render debt");

// Densify / interpolate / snap must commit new samples as active Other waypoints
// (same path model as Add-random): path length === waypoint count, no interior layer.
const densifyPromoteStore = {
  trajectoryMode: "manual",
  manualSelectedVolIds: Array.from({ length: 3 }, (_, i) => "pc1:" + i),
  manualMarkersByVolId: Object.fromEntries(
    Array.from({ length: 3 }, (_, i) => ["pc1:" + i, { plot_row: i }])
  ),
  manualActiveCustomPlotRows: [],
  manualCustomPlotRows: [],
  manualVolumeSnapActive: true,
  manualParticleSnapPathActive: false,
  graphTraversalActive: false,
  manualInterpTrajectoryCount: 0,
  lastLatentTrajectoryPointCount: 0,
  // Simulated densify result: 3 anchors + 1 interior per segment → 5 samples.
  // Interiors are particle rows 100, 101 (to be promoted to Other).
  anchorIndicesActive: [0, 1, 2],
  trajPlotRows: [0, 100, 1, 101, 2],
  anchorTrajXY: [[0, 0], [100, 100], [1, 1], [101, 101], [2, 2]],
  currentPath: null
};
const densifyPromoteSession = new Session({
  store: densifyPromoteStore,
  pathHooks: { rowXYForPlotRow: row => [Number(row), Number(row)] },
  volumes: new VolumeState()
});
const densifyPromotePath = densifyPromoteSession.path();
assert(densifyPromotePath.hasInteriorSamples, "pre-promote densify has sample interiors");
// Promote: activate interiors as Other and install full ordered waypoints.
densifyPromoteStore.manualActiveCustomPlotRows = [100, 101];
densifyPromoteStore.manualCustomPlotRows = [100, 101];
densifyPromoteStore.manualVolumeSnapActive = false;
densifyPromoteStore.manualInterpTrajectoryCount = 0;
densifyPromotePath.setAnchorRows([0, 100, 1, 101, 2]);
densifyPromotePath.setAnchorPathXY([
  [0, 0], [100, 100], [1, 1], [101, 101], [2, 2]
]);
densifyPromotePath.setSamplePlotRows(null);
densifyPromotePath.setInterpolatedCount(0);
assert(!densifyPromotePath.hasInteriorSamples,
  "after densify→Other promote, path length equals waypoints (no interior layer)");
assert(densifyPromoteSession.waypointCount() === 5,
  "densify interiors become waypoints like Add-random");
assert(densifyPromoteSession.sampleCount() === 5,
  "sampleCount matches waypoint count after densify promote");
assert(
  densifyPromoteStore.manualActiveCustomPlotRows.join(",") === "100,101",
  "densify interiors are active Other rows"
);
const densifyIds = densifyPromoteSession.waypointRows().map(row =>
  row < 100 ? ("pc1:" + row) : ("custom:" + row)
);
const densifyVolumes = new VolumeState();
densifyVolumes.alignToIds(densifyIds);
assert(densifyVolumes.slotCount() === 5, "densify-as-Other uses one slot per waypoint");
assert(densifyVolumes.ids().join(",") === "pc1:0,custom:100,pc1:1,custom:101,pc1:2",
  "densify-as-Other slot ids interleave catalog and custom");

// Add-random onto a loaded PC1 catalog: rematch keeps PC decoded; only Other
// slots owe decode (Decode 10, not Decode 20).
const appendDecodeStoreIds = Array.from({ length: 10 }, (_, i) => "pc1:" + i);
const appendDecodeVols = appendDecodeStoreIds.map((_, i) => ({ volume_b64: "PC" + i }));
const appendDecodeImgs = appendDecodeStoreIds.map((_, i) => "img-pc-" + i);
const appendDecodeState = new VolumeState();
appendDecodeState.replaceSlots(appendDecodeStoreIds, appendDecodeVols, appendDecodeImgs);
const appendMixedIds = [];
for (let i = 0; i < 10; i++) {
  appendMixedIds.push("pc1:" + i);
  appendMixedIds.push("custom:" + (200 + i));
}
appendDecodeState.alignToIds(appendMixedIds);
assert(appendDecodeState.slotCount() === 20, "append expands to PC + Other slots");
assert(appendDecodeState.decodedCount() === 10, "PC catalog volumes stay decoded after append");
assert(appendDecodeState.missingDecodeIndices().length === 10,
  "only Other / random slots owe decode after append");
assert(
  appendDecodeState.missingDecodeIndices().every((idx) => idx % 2 === 1),
  "missing decode indices are the Other slots"
);

// Path grow (PC1×10 + 10 Other) must rematch catalog ChimeraX by durable id, not XY.
const growPcIds = Array.from({ length: 10 }, (_, i) => "pc1:" + i);
const growPcImgs = growPcIds.map((_, i) => "img-pc-" + i);
const growState = new VolumeState();
growState.replaceSlots(
  growPcIds,
  growPcIds.map((_, i) => ({ volume_b64: "PC" + i })),
  growPcImgs
);
const grownMixedIds = growPcIds.concat(
  Array.from({ length: 10 }, (_, i) => "custom:" + (300 + i))
);
growState.alignToIds(grownMixedIds);
assert(growState.slotCount() === 20, "alignToIds grows to PC + Other slots");
assert(growState.renderedCount() === 10, "catalog ChimeraX survives path grow");
for (let gpi = 0; gpi < 10; gpi++) {
  assert(growState.isRendered(gpi), "PC slot keeps ChimeraX after grow");
}

// Add-random → clear Other must restore undensified catalog waypoints (PC1 × 10),
// not a densified (nAnchors-1)*(nInterp+1)+1 path such as 19.
const restoreStore = {
  trajectoryMode: "manual",
  manualSelectedVolIds: Array.from({ length: 10 }, (_, i) => "pc1:" + i),
  manualMarkersByVolId: Object.fromEntries(
    Array.from({ length: 10 }, (_, i) => ["pc1:" + i, { plot_row: i }])
  ),
  manualActiveCustomPlotRows: Array.from({ length: 10 }, (_, i) => 1000 + i),
  manualCustomPlotRows: Array.from({ length: 10 }, (_, i) => 1000 + i),
  manualVolumeSnapActive: false,
  manualParticleSnapPathActive: false,
  graphTraversalActive: false,
  manualInterpTrajectoryCount: 0,
  lastLatentTrajectoryPointCount: 0,
  anchorIndicesActive: null,
  anchorTrajXY: null,
  trajPlotRows: null,
  currentPath: null
};
const restoreSession = new Session({
  store: restoreStore,
  pathHooks: {
    rowXYForPlotRow: row => [Number(row), Number(row)]
  },
  volumes: new VolumeState()
});
const restorePath = restoreSession.path();
assert(restorePath && restorePath.kind === "waypoint", "restore path is waypoint mode");
const rebuiltMixed = restorePath.rebuildFromSelection({ force: true });
assert(rebuiltMixed.ok && rebuiltMixed.rows.length === 20,
  "PC1 + 10 random rebuilds to 20 undensified waypoints");
assert(!restorePath.hasInteriorSamples, "mixed append path is not densified");
assert(restoreSession.waypointCount() === 20, "session reports 20 mixed waypoints");

// Simulate deselecting Other: clear active customs, rebuild from catalog selection.
restoreStore.manualActiveCustomPlotRows = [];
const rebuiltCatalog = restorePath.rebuildFromSelection({ force: true, preferActiveOrder: false });
assert(rebuiltCatalog.ok && rebuiltCatalog.rows.length === 10,
  "clearing Other restores 10 catalog waypoints");
assert(rebuiltCatalog.rows.join(",") === "0,1,2,3,4,5,6,7,8,9",
  "restored catalog order matches PC1 selection");
assert(!restorePath.hasInteriorSamples, "restored catalog path has no interiors");
assert(restorePath.getInterpolatedCount() === 0, "restored catalog clears interpolated count");
assert(restoreSession.waypointCount() === 10, "session reports compact PC1 after clearing Other");
assert(restoreSession.sampleCount() === 10, "sampleCount stays compact after clearing Other");
// Path length equals anchors — undensified bookkeeping must not claim densify.
restoreStore.manualInterpTrajectoryCount = 0;
assert(
  !(restoreStore.manualInterpTrajectoryCount > 10
    || (restoreStore.anchorTrajXY && restoreStore.anchorTrajXY.length > 10)),
  "no densified geometry remains after clearing Other"
);

// Reverse must keep decode/render readiness with slot identity (not index).
const reverseReady = new VolumeState();
reverseReady.replaceSlots(
  ["pc1:0", "custom:1", "pc1:2", "custom:3"],
  [
    { volume_b64: "A" },
    null,
    { volume_b64: "C" },
    null
  ],
  ["imgA", null, "imgC", null]
);
assert(reverseReady.missingDecodeIndices().length === 2, "two Other slots owe decode");
assert(reverseReady.missingRenderIndices().length === 2, "two Other slots owe render");
reverseReady.reverse();
assert(reverseReady.ids().join(",") === "custom:3,pc1:2,custom:1,pc1:0",
  "reverse permutes ids with the path");
assert(reverseReady.isDecoded(1) && reverseReady.isRendered(1),
  "PC media stays ready after reverse (by identity)");
assert(reverseReady.isDecoded(3) && reverseReady.isRendered(3),
  "other PC media stays ready after reverse");
assert(reverseReady.missingDecodeIndices().length === 2,
  "decode debt unchanged by reverse");
assert(reverseReady.missingRenderIndices().length === 2,
  "render debt unchanged by reverse");

// Reverse must keep null-id interior media with its geometry (not flip XY alone).
const reverseGeom = new VolumeState();
const revPath = [[0, 0], [1, 0], [2, 0], [3, 0]];
reverseGeom.replaceSlots(
  ["pc1:0", null, null, "pc1:9"],
  [
    { volume_b64: "E0", decoded: true },
    { volume_b64: "I1", decoded: true },
    { volume_b64: "I2", decoded: true },
    { volume_b64: "E1", decoded: true }
  ],
  ["p0", "p1", "p2", "p3"],
  revPath
);
const i1Before = reverseGeom.volumes()[1].volume_b64;
const i1T = Number(reverseGeom.volumes()[1].path_t);
reverseGeom.reverse();
assert(reverseGeom.slotIdAt(0) === "pc1:9" && reverseGeom.slotIdAt(3) === "pc1:0",
  "reverse swaps endpoint catalog ids");
assert(reverseGeom.volumes()[2].volume_b64 === i1Before,
  "interior that was at index 1 moves to index 2 under reverse");
assert(Number(reverseGeom.volumes()[2].path_t).toFixed(3) === "0.667",
  "path_t is restamped to the new slot index after reverse");
assert(Number(reverseGeom.volumes()[1].path_t).toFixed(3) === "0.333",
  "other interior path_t restamped");
assert(
  reverseGeom.slotXyAt(0)[0] === 3 && reverseGeom.volumes()[0].volume_b64 === "E1",
  "endpoint media stays paired with reversed XY"
);
assert(i1T.toFixed(3) === "0.333", "precondition: interior started at t=1/3");

// Compact catalog reverse: extract anchor media from a pre-reverse snapshot,
// reverse the compact arrays, and re-seed onto the densified layout.
const compactRev = new VolumeState();
const compactNAnchors = 3;
const compactNInterp = 1;
const compactPathN = (compactNAnchors - 1) * (compactNInterp + 1) + 1;
const snapVols = new Array(compactPathN).fill(null);
const snapImgs = new Array(compactPathN).fill(null);
const snapIds = new Array(compactPathN).fill(null);
snapVols[0] = { volume_b64: "A", decoded: true };
snapVols[2] = { volume_b64: "B", decoded: true };
snapVols[4] = { volume_b64: "C", decoded: true };
snapIds[0] = "pc1:0";
snapIds[2] = "pc1:5";
snapIds[4] = "pc1:9";
compactRev.alignToPath(null, {
  pathN: compactPathN,
  nAnchors: compactNAnchors,
  nInterp: compactNInterp
});
compactRev.reverseCompactCatalog({
  nAnchors: compactNAnchors,
  nInterp: compactNInterp,
  pathN: compactPathN,
  snapshot: { volumes: snapVols, images: snapImgs, ids: snapIds }
});
assert(compactRev.slotIdAt(0) === "pc1:9",
  "compact reverse anchors first slot to former last id");
assert(compactRev.slotIdAt(4) === "pc1:0",
  "compact reverse anchors last slot to former first id");
assert(compactRev.volumes()[0].volume_b64 === "C",
  "first anchor gets reversed endpoint media");
assert(compactRev.volumes()[4].volume_b64 === "A",
  "last anchor gets reversed endpoint media");

// Original order: PC1 / PC2 / kmeans subsets sort by native index; Other ignored.
const catalogMeta = {
  "pc1:2": { kind: "pc", pc: 1, sample_index: 2 },
  "pc1:0": { kind: "pc", pc: 1, sample_index: 0 },
  "pc1:1": { kind: "pc", pc: 1, sample_index: 1 },
  "pc2:1": { kind: "pc", pc: 2, sample_index: 1 },
  "pc2:0": { kind: "pc", pc: 2, sample_index: 0 },
  "kmeans:2": { kind: "kmeans", cluster_label: 2 },
  "kmeans:0": { kind: "kmeans", cluster_label: 0 },
  "kmeans:1": { kind: "kmeans", cluster_label: 1 }
};
assert(
  Path.sortCatalogVolumeIdsByOriginalIndex(
    ["pc1:2", "pc1:0", "pc1:1"],
    catalogMeta
  ).join(",") === "pc1:0,pc1:1,pc1:2",
  "PC1 original order is sample_index ascending"
);
assert(
  Path.sortCatalogVolumeIdsByOriginalIndex(
    ["kmeans:2", "kmeans:0", "kmeans:1"],
    catalogMeta
  ).join(",") === "kmeans:0,kmeans:1,kmeans:2",
  "kmeans original order is cluster_label ascending"
);
assert(
  Path.sortCatalogVolumeIdsByOriginalIndex(
    ["pc2:1", "pc1:2", "kmeans:1", "pc1:0", "pc2:0", "kmeans:0"],
    catalogMeta
  ).join(",") === "kmeans:0,kmeans:1,pc1:0,pc1:2,pc2:0,pc2:1",
  "mixed catalog sorts by subset then native index"
);

const interleavedSlots = [
  { kind: "custom", row: 101 },
  { kind: "catalog", volId: "pc1:2" },
  { kind: "custom", row: 102 },
  { kind: "catalog", volId: "pc1:0" },
  { kind: "catalog", volId: "pc1:1" },
  { kind: "custom", row: 103 }
];
const interleavedOrig = Path.reorderSlotsByOriginalCatalogIndex(
  interleavedSlots,
  catalogMeta
);
assert(
  interleavedOrig.volIds.join(",") === "pc1:0,pc1:1,pc1:2",
  "Original order rewrites PC1 slots by sample_index"
);
assert(
  interleavedOrig.customRows.join(",") === "101,102,103",
  "Original order leaves Other slots in place"
);
assert(
  interleavedOrig.slots.map(s => s.kind === "custom" ? ("c:" + s.row) : s.volId)
    .join(",") === "c:101,pc1:0,c:102,pc1:1,pc1:2,c:103",
  "Other positions are ignored when ordering the PC1 subset"
);

const kmSlots = [
  { kind: "catalog", volId: "kmeans:2" },
  { kind: "custom", row: 50 },
  { kind: "catalog", volId: "kmeans:0" },
  { kind: "catalog", volId: "kmeans:1" }
];
const kmOrig = Path.reorderSlotsByOriginalCatalogIndex(kmSlots, catalogMeta);
assert(
  kmOrig.slots.map(s => s.kind === "custom" ? ("c:" + s.row) : s.volId)
    .join(",") === "kmeans:0,c:50,kmeans:1,kmeans:2",
  "kmeans subset is index-ordered while Other stays put"
);

// Direct-trace densify: two endpoint catalog ids on a 3-point path (Decode 1).
const DirectTraceCtor = Path.DirectTracePath;
assert(DirectTraceCtor, "DirectTracePath is exported");
const directTraceStore = {
  trajectoryMode: "direct",
  directEndpointVolumeIds: ["pc1:0", "pc1:9"],
  editableTrajXY: [[0, 0], [0.5, 0.5], [1, 1]],
  anchorIndicesActive: null,
  anchorTrajXY: null,
  trajPlotRows: null,
  startXY: [0, 0],
  endXY: [1, 1],
  manualInterpTrajectoryCount: 3,
  manualSelectedVolIds: [],
  manualTrajectoryVolumeIds: [],
  anchorKmeansVolumeIds: [],
  manualMarkersByVolId: {},
  manualCustomPlotRows: [],
  manualActiveCustomPlotRows: [],
  currentPath: null
};
const directTracePath = new DirectTraceCtor(directTraceStore, {});
assert(directTracePath.isDirectTraceMode(), "direct path reports direct-trace mode");
assert(
  directTracePath.getVolumeSlotIds(3).join(",") === "pc1:0,,pc1:9",
  "direct densify places endpoints at 0 and N-1"
);
const insertRes = directTracePath.addCoordinate([0.25, 0.25], { nearestSegment: true });
assert(insertRes && insertRes.ok && insertRes.action === "insert", "insert adds an interior point");
assert(
  directTracePath.getEndpointVolumeIds().join(",") === "pc1:0,pc1:9",
  "insert preserves endpoint catalog ids"
);
assert(directTracePath.sampleCount() === 4, "insert grows path to 4 samples");
// Restore 3-point geometry for densify seed checks below.
directTraceStore.editableTrajXY = [[0, 0], [0.5, 0.5], [1, 1]];
directTraceStore.manualInterpTrajectoryCount = 3;
directTracePath.setEndpointVolumeIds(["pc1:0", "pc1:9"]);
assert(
  directTracePath.getSelectedVolumeIds().join(",") === "pc1:0,pc1:9",
  "direct compact ids are the two endpoints"
);
const directTraceVolumes = new VolumeState();
directTraceVolumes.seedCompactCatalog(
  directTracePath.getSelectedVolumeIds(),
  [],
  [],
  { pathN: 3, nAnchors: 2, nInterp: 1 }
);
assert(directTraceVolumes.slotCount() === 3, "direct densify yields pathN slots");
assert(
  String(directTraceVolumes.ids()[0]) === "pc1:0" && String(directTraceVolumes.ids()[2]) === "pc1:9",
  "seedCompactCatalog places endpoint ids on densified direct path"
);
assert(
  directTraceVolumes.ids()[1] == null,
  "direct interior slot has null id (needs decode)"
);
directTraceVolumes.setDecoded(0, { decoded: true, catalog_id: "pc1:0" });
directTraceVolumes.setDecoded(2, { decoded: true, catalog_id: "pc1:9" });
assert(
  directTraceVolumes.missingDecodeIndices().join(",") === "1",
  "catalog endpoints do not inflate Decode debt on a 3-point direct path"
);
// Freely dragged endpoints keep stale until Decode; seed must not wipe that.
directTraceVolumes.markStale(0);
directTraceVolumes.seedCompactCatalog(
  ["pc1:0", "pc1:9"],
  [{ decoded: true, catalog_id: "pc1:0" }, { decoded: true, catalog_id: "pc1:9" }],
  [{ image_b64: "old-end-0" }, { image_b64: "old-end-2" }],
  { pathN: 3, nAnchors: 2, nInterp: 1 }
);
assert(directTraceVolumes.isStale(0), "seedCompactCatalog preserves stale endpoint");
assert(
  !directTraceVolumes.images()[0],
  "seedCompactCatalog does not reseed media onto a stale endpoint"
);
assert(
  String(directTraceVolumes.ids()[2]) === "pc1:9",
  "seedCompactCatalog still refreshes non-stale endpoint"
);
assert(
  !!directTraceVolumes.images()[2],
  "seedCompactCatalog reseeds media onto a fresh endpoint"
);
// Decode with volume_b64 must not clear stale / allow catalog image reseed.
directTraceVolumes.markStale(2);
directTraceVolumes.setRendered(2, null);
directTraceVolumes.setDecoded(2, { volume_b64: "vtk-end-2", catalog_id: "pc1:9" });
assert(directTraceVolumes.isStale(2), "setDecoded keeps stale until ChimeraX");
directTraceVolumes.seedCompactCatalog(
  ["pc1:0", "pc1:9"],
  [{ decoded: true, catalog_id: "pc1:0" }, { volume_b64: "vtk-end-2", catalog_id: "pc1:9" }],
  [{ image_b64: "old-end-0" }, { image_b64: "stale-catalog-end" }],
  { pathN: 3, nAnchors: 2, nInterp: 1 }
);
assert(directTraceVolumes.isStale(2), "seed skips stale decoded-without-image endpoint");
assert(
  !directTraceVolumes.images()[2],
  "seed does not reseed catalog PNG onto decoded-without-image endpoint"
);
directTraceVolumes.setRendered(2, "fresh-cx-end");
assert(!directTraceVolumes.isStale(2), "setRendered clears stale after ChimeraX");
assert(
  directTraceVolumes.images()[2] === "fresh-cx-end",
  "ChimeraX frame lands on previously stale endpoint"
);
// alignToPath must not drop stale on catalog-id endpoints (pc1:*).
directTraceVolumes.markStale(0);
directTraceVolumes.setRendered(0, null);
directTraceVolumes.seedCompactCatalog(
  ["pc1:0", "pc1:9"],
  [{ decoded: true, catalog_id: "pc1:0" }, { decoded: true, catalog_id: "pc1:9" }],
  [{ image_b64: "resurrect-0" }, { image_b64: "resurrect-2" }],
  { pathN: 3, nAnchors: 2, nInterp: 1 }
);
assert(directTraceVolumes.isStale(0), "seed/align preserves stale on catalog-id endpoint");
assert(
  !directTraceVolumes.images()[0],
  "seed does not resurrect catalog PNG onto stale catalog-id endpoint"
);

// DirectTraceUiState: controls alignment + reverse invalidation maps.
const DirectUi = DirectTraceUiState;
assert(DirectUi, "CryoDirectTraceUiState is exported");
const ui = new DirectUi();
ui.markStale(1);
ui.markDetached(2);
ui.markDragInvalidated(1);
ui.reverseInvalidationMaps(4);
assert(ui.isStale(2), "stale index 1 reverses to 2 on pathN=4");
assert(ui.isDetached(1), "detached index 2 reverses to 1 on pathN=4");
assert(ui.isDragInvalidated(2), "drag-invalidated index 1 reverses to 2");
const uiClear = new DirectUi();
uiClear.markStale(1);
uiClear.markDetached(2);
uiClear.markDragInvalidated(1);
uiClear.clearInvalidationAtIndices([1]);
assert(uiClear.isDetached(2), "clearInvalidationAtIndices preserves detached endpoint labels");
assert(!uiClear.isStale(1), "clearInvalidationAtIndices clears stale");
assert(!uiClear.isDragInvalidated(1), "clearInvalidationAtIndices clears drag-invalidated");
const uiRemap = new DirectUi();
uiRemap.markDetached(0);
uiRemap.markDetached(3);
uiRemap.markStale(3);
uiRemap.remapFlagsBySlotMap({ "0": 0, "1": 4, "2": 7, "3": 11 });
assert(uiRemap.isDetached(0) && uiRemap.isDetached(11),
  "densify rematch map moves detached endpoints to new indices");
assert(!uiRemap.isDetached(3), "orphan detached at old end index is dropped");
assert(uiRemap.isStale(11) && !uiRemap.isStale(3),
  "densify rematch map moves stale with media");
let controlN = 3;
const aligned = ui.alignControlsToLivePath({
  editableTrajXY: [[0, 0], [0.3, 0.3], [0.6, 0.6], [1, 1]],
  controlNPoints: controlN,
  setControlNPoints: function (n) { controlN = n; }
});
assert(aligned.pathN === 4 && controlN === 4 && aligned.controlUpdated,
  "insert-grown path bumps Trajectory-controls nPoints");
let controlKeep = 8;
const noShrink = ui.alignControlsToLivePath({
  editableTrajXY: [[0, 0], [0.5, 0.5], [1, 1]],
  controlNPoints: controlKeep,
  setControlNPoints: function (n) { controlKeep = n; }
});
assert(controlKeep === 8 && noShrink.pathN === 8 && !noShrink.controlUpdated,
  "nPoints control is not pulled down to a shorter editable path");
assert(
  ui.rebuildPointCountForDirectMode({
    editableTrajXY: [[0, 0], [0.5, 0.5], [1, 1], [1.5, 1.5]],
    controlNPoints: 3
  }) === 4,
  "→direct rebuild uses max(control, editable) length"
);
assert(ui.pathDisplayLengthMismatch(5, 3), "detects path vs slider length mismatch");
assert(!ui.pathDisplayLengthMismatch(4, 4), "matching lengths are not a mismatch");

// Nearest snap: stamp must not rebind decode XY; off-path indices drive tick inactivity.
const nearestSnap = new VolumeState();
const freePath = [[0, 0], [0.5, 0.25], [1, 1]];
nearestSnap.replaceSlots(
  ["pc1:0", null, "pc1:9"],
  [
    { volume_b64: "E0", decoded: true },
    { volume_b64: "FREE", decoded: true },
    { volume_b64: "E1", decoded: true }
  ],
  ["p0", "p1", "p2"],
  freePath
);
assert(
  nearestSnap.decodedXyAt(1)[0] === 0.5 && nearestSnap.decodedXyAt(1)[1] === 0.25,
  "decodedXyAt reads vol.traj_xy for free interior"
);
const snappedPath = [[0, 0], [0.7, 0.7], [1, 1]];
nearestSnap.stampPathXy(snappedPath);
assert(
  nearestSnap.decodedXyAt(1)[0] === 0.5 && nearestSnap.decodedXyAt(1)[1] === 0.25,
  "stampPathXy keeps decode traj_xy when live path snaps away"
);
assert(
  nearestSnap.slotXyAt(1)[0] === 0.7 && nearestSnap.slotXyAt(1)[1] === 0.7,
  "stampPathXy still updates live slot XY"
);
assert(
  nearestSnap.indicesOffDecodedPath(snappedPath).join(",") === "1",
  "only the snapped-away interior is off the decoded path"
);
assert(
  nearestSnap.indicesOffDecodedPath(freePath).length === 0,
  "original free path still matches all decoded geometries"
);
// Already-on-particle: full-precision particle XY vs rounded decode must not
// deactivate the tick (matches server direct-mode rounding half-ulp, per-axis).
assert(
  nearestSnap.indicesOffDecodedPath([
    [0, 0],
    [0.5004, 0.2504],
    [1, 1]
  ]).length === 0,
  "sub-rounding drift on both axes does not mark off-path"
);
// Nearest snap: invalidate by movement, not decode-vs-particle XY mismatch.
assert(
  nearestSnap.indicesMovedBetweenPaths(freePath, freePath).length === 0,
  "unmoved path keeps all decoded ticks"
);
assert(
  nearestSnap.indicesMovedBetweenPaths(freePath, [
    [0, 0],
    [0.5003, 0.2503],
    [1, 1]
  ]).length === 0,
  "already-on-particle rounding drift is not treated as movement"
);
assert(
  nearestSnap.indicesMovedBetweenPaths(freePath, snappedPath).join(",") === "1",
  "only the sample that snapped away is marked moved"
);
nearestSnap.setDecoded(1, null);
nearestSnap.setRendered(1, null);
nearestSnap.markStale(1);
assert(!nearestSnap.isDecoded(1) && nearestSnap.isStale(1),
  "cleared off-path slot is not decoded (inactive tick)");
assert(
  nearestSnap.indicesOffDecodedPath(snappedPath).length === 0,
  "no media left off-path after invalidate"
);

// Choose-waypoints: endpoint displacement vs default PC1 vol 1 / vol 10 pair.
function selftestPathEndpointIdsFromRows(rows, selectedIds, markersByVolId) {
  if (!rows || rows.length < 2) return [];
  const volsByRow = {};
  (selectedIds || []).forEach((id) => {
    const marker = markersByVolId[id];
    if (!marker || marker.plot_row == null) return;
    const pr = Number(marker.plot_row);
    if (!Number.isFinite(pr)) return;
    if (!volsByRow[pr]) volsByRow[pr] = [];
    volsByRow[pr].push(String(id));
  });
  const firstList = volsByRow[rows[0]];
  const lastList = volsByRow[rows[rows.length - 1]];
  const firstId = firstList && firstList.length ? firstList[0] : null;
  const lastId = lastList && lastList.length ? lastList[lastList.length - 1] : null;
  return firstId && lastId ? [firstId, lastId] : [];
}
function selftestEndpointsDisplaced(endpoints, defaults) {
  return endpoints.length >= 2 && defaults.length >= 2
    && (endpoints[0] !== defaults[0] || endpoints[1] !== defaults[1]);
}
const pc10Rows = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
const pc10Ids = Array.from({ length: 10 }, (_, i) => "pc1:" + i);
const pc10Markers = Object.fromEntries(
  pc10Ids.map((id, i) => [id, { plot_row: i }])
);
const defaultPc1Endpoints = ["pc1:0", "pc1:9"];
const intactEndpoints = selftestPathEndpointIdsFromRows(pc10Rows, pc10Ids, pc10Markers);
assert(
  !selftestEndpointsDisplaced(intactEndpoints, defaultPc1Endpoints),
  "default PC1×10 path endpoints match vol 1 / vol 10"
);
const reversedPc10Rows = pc10Rows.slice().reverse();
const reversedEndpoints = selftestPathEndpointIdsFromRows(
  reversedPc10Rows,
  pc10Ids,
  pc10Markers
);
assert(
  selftestEndpointsDisplaced(reversedEndpoints, defaultPc1Endpoints),
  "reversed PC1×10 path endpoints are displaced from defaults"
);
assert(
  reversedEndpoints.join(",") === "pc1:9,pc1:0",
  "reversed path reports swapped endpoint catalog ids"
);

// Direct-trace particle-set selection: rebuildFromSelection places one point
// per catalog volume at marker XY (choose-waypoints parity) — never densify.
const directPcStore = {
  trajectoryMode: "direct",
  startXY: [0, 0],
  endXY: [1, 1],
  editableTrajXY: [[0, 0], [0.33, 0.33], [0.66, 0.66], [1, 1]],
  directEndpointVolumeIds: ["pc1:0", "pc1:9"],
  manualSelectedVolIds: Array.from({ length: 10 }, (_, i) => "pc1:" + i),
  manualMarkersByVolId: Object.fromEntries(
    Array.from({ length: 10 }, (_, i) => [
      "pc1:" + i,
      { plot_row: 100 + i, xy: [i * 0.1, i * 0.2], vol_id: "pc1:" + i }
    ])
  ),
  manualActiveCustomPlotRows: [],
  manualCustomPlotRows: [],
  anchorIndicesActive: null,
  anchorTrajXY: null,
  trajPlotRows: null,
  manualInterpTrajectoryCount: 4,
  lastLatentTrajectoryPointCount: 4,
  currentPath: null
};
const directPcPath = new DirectTraceCtor(directPcStore, {});
assert(directPcPath.isDirectTraceMode(), "particle-set direct path stays direct");
const directPcRebuild = directPcPath.rebuildFromSelection({ force: true });
assert(directPcRebuild.ok, "direct-trace particle-set rebuild ok");
assert(directPcRebuild.rows && directPcRebuild.rows.length === 10,
  "direct-trace particle-set path has one point per PC1 volume");
assert(
  directPcRebuild.xy
    && directPcRebuild.xy.length === 10
    && Math.abs(directPcRebuild.xy[0][0] - 0) < 1e-9
    && Math.abs(directPcRebuild.xy[9][0] - 0.9) < 1e-9
    && Math.abs(directPcRebuild.xy[9][1] - 1.8) < 1e-9,
  "direct-trace particle-set points use catalog marker XY"
);
assert(!directPcStore.editableTrajXY, "particle-set rebuild clears editable polyline");
assert(directPcStore.manualInterpTrajectoryCount === 0,
  "particle-set rebuild clears densify count");

console.log("session_path_sync: ok");
