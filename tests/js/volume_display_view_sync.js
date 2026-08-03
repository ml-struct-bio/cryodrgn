/**
 * The viewing angle shared between the VTK raycaster and ChimeraX renders,
 * including the snapshot frozen for the duration of one render batch.
 *
 * Run by ``TestTrajectoryScenarioScripts`` in ``tests/test_dashboard_trajectory.py``.
 */
"use strict";

const {
  Display, assert
} = require("./_harness.js");

// Display freezes one ChimeraX view snapshot for a render batch.
if (Display) {
  const display = new Display({
    asideHostEl: { appendChild: function() {} },
    captureChimeraxViewSnapshot: function() {
      return {
        view_matrix: "",
        view_turns: [{ axis: "y", degrees: 45 }],
        viewKey: "t:y:45.000",
        isoKey: "iso:auto",
        iso_level: null
      };
    }
  });
  const live = display.chimeraxViewSnapshot();
  assert(live.view_turns.length === 1 && live.view_turns[0].degrees === 45,
    "Display.chimeraxViewSnapshot reads live capture hook");
  display.beginChimeraxViewBatch({
    view_matrix: "",
    view_turns: [{ axis: "x", degrees: 10 }],
    viewKey: "t:x:10.000",
    isoKey: "iso:auto"
  });
  // Mutate the "live" hook — batch must stay frozen.
  display.captureChimeraxViewSnapshot = function() {
    return {
      view_turns: [{ axis: "z", degrees: 90 }],
      viewKey: "t:z:90.000",
      isoKey: "iso:auto"
    };
  };
  const frozen = display.chimeraxViewSnapshot();
  assert(frozen.viewKey === "t:x:10.000"
      && frozen.view_turns[0].axis === "x",
    "Display batch snapshot ignores live rotation changes");
  const payload = display.applyChimeraxViewToPayload({});
  assert(payload.view_turns && payload.view_turns[0].degrees === 10
      && !payload.view_matrix,
    "applyChimeraxViewToPayload uses frozen batch turns");
  display.endChimeraxViewBatch();
  assert(display.chimeraxViewBatchSnapshot() == null,
    "endChimeraxViewBatch clears freeze");

  // VTK ↔ ChimeraX backend switches share the on-screen viewing angle.
  var syncedToCx = null;
  var syncDisplay = new Display({
    asideHostEl: { appendChild: function() {} },
    onSyncVtkViewToChimerax: function(snap) { syncedToCx = snap; }
  });
  syncDisplay.backend = "vtk";
  syncDisplay.volumes = [{ volume_b64: "AAAA", D: 2 }];
  syncDisplay.raycastView = {
    volume: {},
    getChimeraxViewMatrixCamera: function() {
      return "1,0,0,0,0,1,0,0,0,0,1,0";
    },
    getChimeraxViewTurns: function() {
      return [{ axis: "y", degrees: 30 }];
    }
  };
  syncDisplay.setBackend("chimerax");
  assert(syncedToCx && syncedToCx.view_turns
      && syncedToCx.view_turns.length === 1
      && String(syncedToCx.view_turns[0].axis).indexOf(",") < 0
      && syncedToCx.view_turns[0].axis === "y"
      && Math.abs(syncedToCx.view_turns[0].degrees - 30) < 1e-6
      && !syncedToCx.view_matrix,
    "VTK→ChimeraX sync pushes live view turns to page hook");
  assert(syncDisplay._vtkChimeraxSyncTurns.length === 1
      && Math.abs(syncDisplay._vtkChimeraxSyncTurns[0].degrees - 30) < 1e-6,
    "VTK→ChimeraX sync stores turns on display");

  syncDisplay.setChimeraxRenderedViewMatrix(
    "camera 0,0,1,0,0,1,0,0,-1,0,0,0"
  );
  syncDisplay.chimeraxImages = ["png"];
  syncDisplay.backend = "chimerax";
  syncDisplay.setBackend("vtk");
  assert(syncDisplay._pendingApplyViewTurnsToVtk
      && syncDisplay._pendingApplyViewTurnsToVtk.some(function (t) {
        return t.axis === "y" && Math.abs(Math.abs(t.degrees) - 90) < 1.0;
      }),
    "ChimeraX→VTK sync converts rendered matrix rotation to turns");
  assert(syncDisplay._pendingApplyViewMatrixToVtk !== true,
    "ChimeraX→VTK sync must not apply translated ChimeraX absolute matrix");

  // Identity ChimeraX matrix (translation only) must not yank the VTK camera,
  // and must not re-apply a stale pre-ChimeraX VTK matrix on repeated switches.
  var identityCx = new Display({
    asideHostEl: { appendChild: function() {} }
  });
  identityCx.backend = "chimerax";
  identityCx.chimeraxImages = ["png"];
  identityCx._vtkViewMatrixBeforeChimerax =
    "0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000";
  identityCx.setChimeraxRenderedViewMatrix(
    "camera 1,0,0,178.6,0,1,0,152.72,0,0,1,781.56"
  );
  identityCx.setBackend("vtk");
  assert(!identityCx._pendingApplyViewMatrixToVtk
      && !(identityCx._pendingApplyViewTurnsToVtk
        && identityCx._pendingApplyViewTurnsToVtk.length),
    "identity ChimeraX matrix leaves VTK default framing (no stale matrix)");

  // VTK→ChimeraX sync prefers orient-relative turns (not zero-T matrices).
  var syncedToCx = null;
  var vtkToCx = new Display({
    asideHostEl: { appendChild: function() {} },
    onSyncVtkViewToChimerax: function(snap) { syncedToCx = snap; }
  });
  vtkToCx.backend = "vtk";
  vtkToCx.volumes = [{ volume_b64: "AAAA", D: 2 }];
  vtkToCx.raycastView = {
    volume: {},
    getChimeraxViewMatrixCamera: function() {
      return "0.707,0,0.707,0,0,1,0,0,-0.707,0,0.707,0";
    },
    getChimeraxViewTurns: function() {
      return [{ axis: "y", degrees: 45 }];
    }
  };
  vtkToCx.setBackend("chimerax");
  assert(syncedToCx
      && syncedToCx.view_turns
      && syncedToCx.view_turns.length === 1
      && Math.abs(syncedToCx.view_turns[0].degrees - 45) < 1e-6
      && !syncedToCx.view_matrix,
    "VTK→ChimeraX sync sends turns (not zero-T matrix)");

  var turnsOnly = new Display({
    asideHostEl: { appendChild: function() {} }
  });
  turnsOnly.backend = "chimerax";
  turnsOnly.chimeraxImages = ["png"];
  turnsOnly.setChimeraxRenderedViewTurns([{ axis: "x", degrees: 15 }]);
  turnsOnly.setBackend("vtk");
  assert(turnsOnly._pendingApplyViewTurnsToVtk
      && turnsOnly._pendingApplyViewTurnsToVtk[0].axis === "x"
      && Math.abs(turnsOnly._pendingApplyViewTurnsToVtk[0].degrees - 15) < 1e-6,
    "ChimeraX→VTK sync falls back to rendered view turns");
}

console.log("volume_display_view_sync: ok");
