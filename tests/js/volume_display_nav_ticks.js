/**
 * Slider and keyboard navigation over volume ticks: what counts as a ready
 * stopping point, and where an unready target snaps to.
 *
 * Run by ``TestTrajectoryScenarioScripts`` in ``tests/test_dashboard_trajectory.py``.
 */
"use strict";

const {
  VolumeState, VolumeStateUtils, Session, Mutations, Path, Display,
  DirectTraceUiState, assert, elStub
} = require("./_harness.js");

// Nav/tick readiness must be one predicate (_slotTickReadyAt): a slot the page
// hook reports ready (Session-backed) is a valid slider/keyboard stopping
// point even before its bytes are mirrored into this display's local arrays.
// Regression for: ticks render active but the slider skips past them.
if (Display) {
  var navReadyHook = { 2: true, 4: true };
  var navSlider = elStub({ value: "0" });
  navSlider.setAttribute = function (name, val) { this["_" + name] = val; };
  var navPrevBtn = elStub({ disabled: false });
  var navNextBtn = elStub({ disabled: false });
  var navDisplay = new Display({
    asideHostEl: { appendChild: function () {} },
    volumeSliderEl: navSlider,
    volumeNavEl: {
      classList: { toggle: function () {} },
      setAttribute: function () {},
      hidden: false
    },
    btnVolPrev: navPrevBtn,
    btnVolNext: navNextBtn,
    slotReadyAt: function (i) { return !!navReadyHook[i]; }
  });
  navDisplay.backend = "chimerax";
  navDisplay.expectedVolumeCount = 5;
  // No local bytes anywhere — the hook alone marks slots 2 and 4 ready.
  navDisplay.chimeraxImages = [null, null, null, null, null];

  assert(navDisplay._slotTickReadyAt(2) === true,
    "hook-ready slot reports tick-ready without local bytes");
  assert(navDisplay.inactiveTickIndices().join(",") === "0,1,3",
    "hook-ready slots 2 and 4 are excluded from the inactive-tick set");

  assert(navDisplay._snapFocusIndexToReady(2) === 2,
    "slider snap stops exactly on a tick the chrome renders active");
  assert(navDisplay._snapFocusIndexToReady(3) === 2 || navDisplay._snapFocusIndexToReady(3) === 4,
    "an unready target snaps to its nearest hook-ready neighbour");
  assert(navDisplay._snapFocusIndexToReady(0) === 2,
    "snapping searches outward in both directions for the nearest ready tick");

  navDisplay.setFocusIndex(2);
  assert(navDisplay.getFocusIndex() === 2,
    "setFocusIndex lands exactly on the hook-ready tick, not a neighbour");

  navDisplay._syncVolumeNavChrome();
  assert(navSlider.value === "2",
    "nav chrome parks the slider on the hook-ready focus with no local bytes");
  assert(navPrevBtn.disabled === false && navNextBtn.disabled === false,
    "Prev/Next stay enabled while more than one hook-ready tick exists");

  navDisplay._stepReadyFocus(1);
  assert(navDisplay.chimeraxFocusIndex === 4,
    "stepping forward from a hook-ready tick lands on the next one, skipping unready slot 3");

  // No hook at all: falls back to the raw per-backend check (pre-existing
  // behaviour), so a display with real local bytes still navigates correctly.
  var rawDisplay = new Display({ asideHostEl: { appendChild: function () {} } });
  rawDisplay.backend = "chimerax";
  rawDisplay.expectedVolumeCount = 3;
  rawDisplay.chimeraxImages = [null, "img-b64", null];
  assert(rawDisplay._snapFocusIndexToReady(0) === 1,
    "without a hook, snap falls back to raw chimerax image presence");
}

console.log("volume_display_nav_ticks: ok");
