/**
 * Snapshot / restore Plotly ``scatter3d`` scene camera and axis ranges (orbit + zoom box).
 * Used when ``Plotly.react`` replaces trace data (e.g. colour covariate) without changing axes.
 */
(function(global) {
  "use strict";

  function snapshot(gd) {
    try {
      var scene = gd._fullLayout && gd._fullLayout.scene;
      if (!scene) return null;
      var snap = {};
      var camSrc = null;
      if (scene.camera && typeof scene.camera === "object") {
        camSrc = scene.camera;
      } else {
        try {
          var ls = gd.layout && gd.layout.scene && gd.layout.scene.camera;
          if (ls && typeof ls === "object") camSrc = ls;
        } catch (eCam) { /* ignore */ }
      }
      if (camSrc) {
        snap.camera = JSON.parse(JSON.stringify(camSrc));
      }
      function copyRange(axis) {
        if (!axis) return null;
        if (axis.range && axis.range.length >= 2) {
          return [axis.range[0], axis.range[1]];
        }
        try {
          var rl = axis._rl;
          if (rl && rl.length >= 2) return [rl[0], rl[1]];
        } catch (e) { /* ignore */ }
        return null;
      }
      snap.xrange = copyRange(scene.xaxis);
      snap.yrange = copyRange(scene.yaxis);
      snap.zrange = copyRange(scene.zaxis);
      return snap;
    } catch (e) {
      return null;
    }
  }

  function sceneRelayoutPatchFromSnap(snap) {
    if (!snap) return {};
    var patch = {};
    if (snap.camera) patch["scene.camera"] = snap.camera;
    if (snap.xrange) {
      patch["scene.xaxis.range"] = snap.xrange;
      patch["scene.xaxis.autorange"] = false;
    }
    if (snap.yrange) {
      patch["scene.yaxis.range"] = snap.yrange;
      patch["scene.yaxis.autorange"] = false;
    }
    if (snap.zrange) {
      patch["scene.zaxis.range"] = snap.zrange;
      patch["scene.zaxis.autorange"] = false;
    }
    return patch;
  }

  /**
   * Layout patch for ``Plotly.update`` alongside trace edits: **camera only**.
   * Applying ``scene.*axis.range`` without ``scene.camera`` in the same update can reset the
   * scatter3d orbit (often noticeable on the first volume selection after orbiting).
   */
  function sceneRelayoutPatchCameraOnly(snap) {
    if (!snap || !snap.camera || typeof snap.camera !== "object") return {};
    return { "scene.camera": snap.camera };
  }

  function restore(gd, snap) {
    if (!snap || typeof Plotly === "undefined" || !Plotly.relayout) {
      return Promise.resolve();
    }
    var patch = sceneRelayoutPatchFromSnap(snap);
    if (!Object.keys(patch).length) return Promise.resolve();
    return Plotly.relayout(gd, patch).catch(function() {});
  }

  function scheduleRestore(gd, snap, onDone) {
    requestAnimationFrame(function() {
      var p = restore(gd, snap);
      var chain = p && typeof p.then === "function" ? p.catch(function() {}) : Promise.resolve();
      chain.then(function() {
        requestAnimationFrame(function() {
          var p2 = restore(gd, snap);
          var chain2 = p2 && typeof p2.then === "function" ? p2.catch(function() {}) : Promise.resolve();
          chain2.then(function() {
            if (typeof onDone === "function") {
              try {
                onDone();
              } catch (eDone) { /* ignore */ }
            }
          });
        });
      });
    });
  }

  /**
   * Patch ``fig.layout.scene`` on an incoming figure JSON so ``Plotly.react`` applies the
   * user's orbit/zoom instead of the server defaults (colour-only redraws).
   */
  function applySnapshotToFigureLayout(fig, snap) {
    if (!fig || !fig.layout || !snap) return;
    try {
      fig.layout.scene = fig.layout.scene || {};
      var scn = fig.layout.scene;
      if (snap.camera && typeof snap.camera === "object") {
        scn.camera = JSON.parse(JSON.stringify(snap.camera));
      }
      function applyRange(axisName, rangeArr) {
        if (!rangeArr || rangeArr.length < 2) return;
        scn[axisName] = scn[axisName] || {};
        scn[axisName].range = [rangeArr[0], rangeArr[1]];
        scn[axisName].autorange = false;
      }
      applyRange("xaxis", snap.xrange);
      applyRange("yaxis", snap.yrange);
      applyRange("zaxis", snap.zrange);
    } catch (e) { /* ignore */ }
  }

  global.CryoPlotlyScatter3dScene = {
    snapshot: snapshot,
    applySnapshotToFigureLayout: applySnapshotToFigureLayout,
    sceneRelayoutPatchFromSnap: sceneRelayoutPatchFromSnap,
    sceneRelayoutPatchCameraOnly: sceneRelayoutPatchCameraOnly,
    restore: restore,
    scheduleRestore: scheduleRestore,
  };
})(typeof window !== "undefined" ? window : globalThis);
