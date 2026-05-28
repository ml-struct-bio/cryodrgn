/**
 * Snapshot / restore Plotly ``scatter3d`` scene camera and axis ranges (orbit + zoom box).
 * Used when ``Plotly.react`` replaces trace data (e.g. colour covariate) without changing axes.
 */
(function(global) {
  "use strict";

  function snapHasCamera(snap) {
    if (!snap || !snap.camera || typeof snap.camera !== "object") return false;
    var eye = snap.camera.eye;
    return !!(eye && typeof eye === "object");
  }

  function snapHasAxisRanges(snap) {
    return !!(
      snap
      && (
        (snap.xrange && snap.xrange.length >= 2)
        || (snap.yrange && snap.yrange.length >= 2)
        || (snap.zrange && snap.zrange.length >= 2)
      )
    );
  }

  function getPinnedSnap(gd) {
    if (!gd) return null;
    var pin = gd._cryoLatent3dScenePin;
    return snapHasCamera(pin) ? pin : null;
  }

  function plotlyDefaultSceneCamera() {
    return {
      up: { x: 0, y: 0, z: 1 },
      center: { x: 0, y: 0, z: 0 },
      eye: { x: 1.25, y: 1.25, z: 1.25 },
    };
  }

  function normalizeSceneCamera(partial) {
    var d = plotlyDefaultSceneCamera();
    if (!partial || typeof partial !== "object") {
      return JSON.parse(JSON.stringify(d));
    }
    var up = partial.up && typeof partial.up === "object"
      ? { x: +partial.up.x, y: +partial.up.y, z: +partial.up.z }
      : d.up;
    var center = partial.center && typeof partial.center === "object"
      ? { x: +partial.center.x, y: +partial.center.y, z: +partial.center.z }
      : d.center;
    var eye = partial.eye && typeof partial.eye === "object"
      ? { x: +partial.eye.x, y: +partial.eye.y, z: +partial.eye.z }
      : d.eye;
    return { up: up, center: center, eye: eye };
  }

  function cloneSnap(snap) {
    if (!snapHasCamera(snap)) return null;
    return {
      camera: JSON.parse(JSON.stringify(normalizeSceneCamera(snap.camera))),
      xrange: snap.xrange ? snap.xrange.slice() : snap.xrange,
      yrange: snap.yrange ? snap.yrange.slice() : snap.yrange,
      zrange: snap.zrange ? snap.zrange.slice() : snap.zrange,
    };
  }

  function setPinnedSnap(gd, snap) {
    if (!gd || !snapHasCamera(snap)) return;
    gd._cryoLatent3dScenePin = cloneSnap(snap);
  }

  function beginRestore(gd) {
    if (!gd) return;
    gd._cryoSceneRestoreDepth = (gd._cryoSceneRestoreDepth || 0) + 1;
  }

  function endRestore(gd) {
    if (!gd || !gd._cryoSceneRestoreDepth) return;
    gd._cryoSceneRestoreDepth--;
    if (gd._cryoSceneRestoreDepth < 0) gd._cryoSceneRestoreDepth = 0;
  }

  function isRestoring(gd) {
    return !!(gd && gd._cryoSceneRestoreDepth > 0);
  }

  function cameraEyeDiffers(a, b, eps) {
    eps = eps == null ? 1e-5 : eps;
    if (!a || !b || !a.camera || !b.camera) return true;
    var ae = a.camera.eye;
    var be = b.camera.eye;
    if (!ae || !be) return true;
    return (
      Math.abs((ae.x || 0) - (be.x || 0)) > eps
      || Math.abs((ae.y || 0) - (be.y || 0)) > eps
      || Math.abs((ae.z || 0) - (be.z || 0)) > eps
    );
  }

  /**
   * Active orbit lives on the WebGL scene during drag; ``_fullLayout.scene.camera`` can lag
   * until mouseup, so prefer ``scene._scene.getCamera()`` when available.
   */
  function resolveCameraFromGd(gd) {
    try {
      var scene = gd._fullLayout && gd._fullLayout.scene;
      if (!scene) return null;
      if (scene._scene && typeof scene._scene.getCamera === "function") {
        return JSON.parse(JSON.stringify(normalizeSceneCamera(scene._scene.getCamera())));
      }
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
        return JSON.parse(JSON.stringify(normalizeSceneCamera(camSrc)));
      }
    } catch (e) { /* ignore */ }
    return null;
  }

  function snapshot(gd) {
    try {
      var scene = gd._fullLayout && gd._fullLayout.scene;
      if (!scene) return null;
      var snap = {};
      var cam = resolveCameraFromGd(gd);
      if (cam) {
        snap.camera = cam;
      }
      var layoutScene = null;
      try {
        layoutScene = gd.layout && gd.layout.scene;
      } catch (eLay) { /* ignore */ }
      function copyRange(axis, axisName) {
        if (axis) {
          if (axis.range && axis.range.length >= 2) {
            return [axis.range[0], axis.range[1]];
          }
          try {
            var rl = axis._rl;
            if (rl && rl.length >= 2) return [rl[0], rl[1]];
          } catch (e) { /* ignore */ }
        }
        try {
          var layAx = layoutScene && layoutScene[axisName];
          if (layAx && layAx.range && layAx.range.length >= 2) {
            return [layAx.range[0], layAx.range[1]];
          }
        } catch (e2) { /* ignore */ }
        return null;
      }
      snap.xrange = copyRange(scene.xaxis, "xaxis");
      snap.yrange = copyRange(scene.yaxis, "yaxis");
      snap.zrange = copyRange(scene.zaxis, "zaxis");
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

  /** Preserve orbit only — avoids pinning stale ``scene.*axis.range`` after trace restyle / colour redraw. */
  function restoreCameraOnly(gd, snap) {
    if (!snap || typeof Plotly === "undefined" || !Plotly.relayout) {
      return Promise.resolve();
    }
    var patch = sceneRelayoutPatchCameraOnly(snap);
    if (!Object.keys(patch).length) return Promise.resolve();
    beginRestore(gd);
    var p = Plotly.relayout(gd, patch).catch(function() {});
    if (p && typeof p.then === "function") {
      return p.then(
        function() {
          endRestore(gd);
        },
        function() {
          endRestore(gd);
        }
      );
    }
    endRestore(gd);
    return Promise.resolve();
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

  function scheduleRestoreCameraOnly(gd, snap, onDone) {
    requestAnimationFrame(function() {
      var p = restoreCameraOnly(gd, snap);
      var chain = p && typeof p.then === "function" ? p.catch(function() {}) : Promise.resolve();
      chain.then(function() {
        requestAnimationFrame(function() {
          var p2 = restoreCameraOnly(gd, snap);
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
   * Re-apply orbit after late layout (aside resize, annotation relayout, ``Plotly.Plots.resize``).
   * Triple restore spans two animation frames plus a timer tick so it runs after Plotly settles.
   */
  function scheduleEnforceCamera(gd, snap, onDone) {
    var pin = cloneSnap(snap);
    if (!pin) {
      if (typeof onDone === "function") {
        try {
          onDone();
        } catch (eDone) { /* ignore */ }
      }
      return;
    }
    function finish() {
      var live = snapshot(gd);
      if (snapHasCamera(live)) {
        setPinnedSnap(gd, live);
      } else if (snapHasCamera(pin)) {
        setPinnedSnap(gd, pin);
      }
      if (typeof onDone === "function") {
        try {
          onDone();
        } catch (eDone2) { /* ignore */ }
      }
    }
    function chainRestore(ms) {
      return new Promise(function(res) {
        setTimeout(function() {
          var p = restorePinnedView(gd, pin);
          var c = p && typeof p.then === "function" ? p.catch(function() {}) : Promise.resolve();
          c.then(res);
        }, ms);
      });
    }
    requestAnimationFrame(function() {
      var p0 = restorePinnedView(gd, pin);
      var c0 = p0 && typeof p0.then === "function" ? p0.catch(function() {}) : Promise.resolve();
      c0.then(function() {
        requestAnimationFrame(function() {
          var p1 = restorePinnedView(gd, pin);
          var c1 = p1 && typeof p1.then === "function" ? p1.catch(function() {}) : Promise.resolve();
          c1.then(function() {
            chainRestore(0)
              .then(function() { return chainRestore(50); })
              .then(function() { return chainRestore(150); })
              .then(finish);
          });
        });
      });
    });
  }

  function markerRestyleUpdateFromFigure(fig) {
    var tr = fig.data[0];
    var upd = {
      marker: [tr.marker],
      customdata: [tr.customdata],
    };
    if (tr.text !== undefined) upd.text = [tr.text];
    if (tr.textfont !== undefined) upd.textfont = [tr.textfont];
    if (tr.textposition !== undefined) upd.textposition = [tr.textposition];
    if (tr.hovertemplate !== undefined) upd.hovertemplate = [tr.hovertemplate];
    if (tr.mode !== undefined) upd.mode = [tr.mode];
    return upd;
  }

  function traceCanMarkerRestyle(fig, gd) {
    if (!fig || !fig.data || !fig.data[0] || !gd || !gd.data || !gd.data[0]) {
      return false;
    }
    var tr = fig.data[0];
    var live = gd.data[0];
    return !!(tr.x && live.x && tr.x.length === live.x.length);
  }

  function traceRestyleFromFigure(fig, gd) {
    if (!traceCanMarkerRestyle(fig, gd)) {
      return Promise.resolve(false);
    }
    if (typeof Plotly === "undefined" || !Plotly.restyle) {
      return Promise.resolve(false);
    }
    var upd = markerRestyleUpdateFromFigure(fig);
    var p = Plotly.restyle(gd, upd, [0]);
    return p && typeof p.then === "function" ? p.then(function() { return true; }) : Promise.resolve(true);
  }

  /**
   * Marker/colour restyle plus camera in one ``Plotly.update`` — scatter3d resets orbit if the
   * camera is restored only in a follow-up relayout after ``Plotly.restyle``.
   */
  function traceMarkerRestylePreservingCamera(fig, gd, pin) {
    if (!traceCanMarkerRestyle(fig, gd)) {
      return Promise.resolve(false);
    }
    if (typeof Plotly === "undefined") {
      return Promise.resolve(false);
    }
    var upd = markerRestyleUpdateFromFigure(fig);
    var pinSnap = snapHasCamera(pin) ? cloneSnap(pin) : null;
    if (pinSnap && typeof Plotly.update === "function") {
      var viewPatch = snapHasAxisRanges(pinSnap)
        ? sceneRelayoutPatchFromSnap(pinSnap)
        : sceneRelayoutPatchCameraOnly(pinSnap);
      var layoutPatch = Object.assign(viewPatch, disableSceneSelectionRelayoutPatch());
      beginRestore(gd);
      var pUpd = Plotly.update(gd, upd, layoutPatch, [0]);
      var chain = pUpd && typeof pUpd.then === "function" ? pUpd.catch(function() {}) : Promise.resolve();
      return chain.then(
        function() {
          endRestore(gd);
          return true;
        },
        function() {
          endRestore(gd);
          return true;
        }
      );
    }
    if (!Plotly.restyle) {
      return Promise.resolve(false);
    }
    var p = Plotly.restyle(gd, upd, [0]);
    return p && typeof p.then === "function" ? p.then(function() { return true; }) : Promise.resolve(true);
  }

  function traceCoordsRestyleFromFigure(fig) {
    if (!fig || !fig.data || !fig.data[0]) return {};
    return traceArrayUpdateFromFigure(fig);
  }

  function relayoutCameraAndAxisPatch(pin, layout, includeAxisRanges) {
    var patch = sceneRelayoutPatchCameraOnly(pin);
    if (includeAxisRanges) {
      if (snapHasAxisRanges(pin)) {
        var pinAxis = sceneRelayoutPatchFromSnap(pin);
        Object.keys(pinAxis).forEach(function(k) {
          if (k !== "scene.camera") patch[k] = pinAxis[k];
        });
      } else if (layout) {
        Object.assign(patch, sceneAxisRangeRelayoutPatchFromLayout(layout));
      }
    }
    Object.assign(patch, disableSceneSelectionRelayoutPatch());
    return patch;
  }

  /** Camera plus pinned axis limits — never server layout ranges on legend filter redraws. */
  function relayoutPreserveViewPatch(pin) {
    var patch = sceneRelayoutPatchFromSnap(pin);
    Object.assign(patch, disableSceneSelectionRelayoutPatch());
    return patch;
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

  function layoutPatchWithoutScene(layout) {
    var patch = {};
    if (!layout || typeof layout !== "object") return patch;
    Object.keys(layout).forEach(function(k) {
      if (k === "scene") return;
      patch[k] = layout[k];
    });
    return patch;
  }

  /** Axis limits from server layout without replacing ``scene.camera``. */
  function sceneAxisRangeRelayoutPatchFromLayout(layout) {
    var patch = {};
    var scn = layout && layout.scene;
    if (!scn || typeof scn !== "object") return patch;
    ["xaxis", "yaxis", "zaxis"].forEach(function(ax) {
      var axis = scn[ax];
      if (!axis || !axis.range || axis.range.length < 2) return;
      patch["scene." + ax + ".range"] = [axis.range[0], axis.range[1]];
      patch["scene." + ax + ".autorange"] = false;
    });
    return patch;
  }

  function traceArrayUpdateFromFigure(fig) {
    if (!fig || !fig.data || !fig.data[0]) return {};
    var tr = fig.data[0];
    var upd = {
      x: [tr.x],
      y: [tr.y],
      z: [tr.z],
      mode: [tr.mode],
      marker: [tr.marker],
      customdata: [tr.customdata],
    };
    if (tr.text !== undefined) upd.text = [tr.text];
    if (tr.textfont !== undefined) upd.textfont = [tr.textfont];
    if (tr.textposition !== undefined) upd.textposition = [tr.textposition];
    if (tr.hovertemplate !== undefined) upd.hovertemplate = [tr.hovertemplate];
    if (tr.ids !== undefined) upd.ids = [tr.ids];
    return upd;
  }

  function disableSceneSelectionRelayoutPatch() {
    return {
      "scene.dragmode": "orbit",
      selections: [],
    };
  }

  function stripSceneAxisRangesFromFigureLayout(fig) {
    if (!fig || !fig.layout || !fig.layout.scene) return;
    var scn = fig.layout.scene;
    ["xaxis", "yaxis", "zaxis"].forEach(function(axName) {
      var ax = scn[axName];
      if (!ax || typeof ax !== "object") return;
      delete ax.range;
      delete ax.autorange;
    });
  }

  /**
   * Merge ``pin`` with live axis limits from ``gd`` when the pin only has a camera
   * (common before discrete legend toggles that replace the trace via ``Plotly.react``).
   */
  function effectiveViewSnapForRedraw(gd, pin) {
    var out = cloneSnap(pin);
    var live = null;
    try {
      live = snapshot(gd);
    } catch (e) { /* ignore */ }
    if (!live) return out;
    if (!out) return cloneSnap(live);
    if (!out.camera && live.camera) {
      out.camera = JSON.parse(JSON.stringify(live.camera));
    }
    if (!snapHasAxisRanges(out)) {
      if (live.xrange) out.xrange = live.xrange.slice();
      if (live.yrange) out.yrange = live.yrange.slice();
      if (live.zrange) out.zrange = live.zrange.slice();
    }
    return out;
  }

  /**
   * Before ``Plotly.react`` on same-axes legend filters: drop server ``scene.*.range`` and
   * apply the user's current view (orbit + zoom box).
   */
  function prepareFigureSceneForPinnedViewRedraw(fig, gd, pin) {
    if (!fig || !fig.layout || !pin || !gd) return null;
    var view = effectiveViewSnapForRedraw(gd, pin);
    if (!view || !snapHasCamera(view)) return null;
    stripSceneAxisRangesFromFigureLayout(fig);
    applySnapshotToFigureLayout(fig, view);
    return view;
  }

  function restorePinnedView(gd, snap) {
    return snapHasAxisRanges(snap) ? restore(gd, snap) : restoreCameraOnly(gd, snap);
  }

  /** Orbit only for colour-covariate / marker restyle redraws (server axis limits stay authoritative). */
  function applySnapshotCameraOnlyToFigureLayout(fig, snap) {
    if (!fig || !fig.layout || !snap || !snap.camera || typeof snap.camera !== "object") return;
    try {
      fig.layout.scene = fig.layout.scene || {};
      fig.layout.scene.camera = JSON.parse(JSON.stringify(snap.camera));
    } catch (e) { /* ignore */ }
  }

  global.CryoPlotlyScatter3dScene = {
    snapshot: snapshot,
    snapHasCamera: snapHasCamera,
    snapHasAxisRanges: snapHasAxisRanges,
    getPinnedSnap: getPinnedSnap,
    setPinnedSnap: setPinnedSnap,
    isRestoring: isRestoring,
    cameraEyeDiffers: cameraEyeDiffers,
    applySnapshotToFigureLayout: applySnapshotToFigureLayout,
    applySnapshotCameraOnlyToFigureLayout: applySnapshotCameraOnlyToFigureLayout,
    stripSceneAxisRangesFromFigureLayout: stripSceneAxisRangesFromFigureLayout,
    effectiveViewSnapForRedraw: effectiveViewSnapForRedraw,
    prepareFigureSceneForPinnedViewRedraw: prepareFigureSceneForPinnedViewRedraw,
    restorePinnedView: restorePinnedView,
    sceneRelayoutPatchFromSnap: sceneRelayoutPatchFromSnap,
    sceneRelayoutPatchCameraOnly: sceneRelayoutPatchCameraOnly,
    restore: restore,
    restoreCameraOnly: restoreCameraOnly,
    scheduleRestore: scheduleRestore,
    scheduleRestoreCameraOnly: scheduleRestoreCameraOnly,
    scheduleEnforceCamera: scheduleEnforceCamera,
    cloneSnap: cloneSnap,
    layoutPatchWithoutScene: layoutPatchWithoutScene,
    sceneAxisRangeRelayoutPatchFromLayout: sceneAxisRangeRelayoutPatchFromLayout,
    traceArrayUpdateFromFigure: traceArrayUpdateFromFigure,
    traceCoordsRestyleFromFigure: traceCoordsRestyleFromFigure,
    traceRestyleFromFigure: traceRestyleFromFigure,
    traceMarkerRestylePreservingCamera: traceMarkerRestylePreservingCamera,
    normalizeSceneCamera: normalizeSceneCamera,
    plotlyDefaultSceneCamera: plotlyDefaultSceneCamera,
    resolveCameraFromGd: resolveCameraFromGd,
    relayoutCameraAndAxisPatch: relayoutCameraAndAxisPatch,
    relayoutPreserveViewPatch: relayoutPreserveViewPatch,
    disableSceneSelectionRelayoutPatch: disableSceneSelectionRelayoutPatch,
  };
})(typeof window !== "undefined" ? window : globalThis);
