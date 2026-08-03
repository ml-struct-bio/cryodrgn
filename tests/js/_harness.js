/**
 * Shared setup for the trajectory scenario scripts in this directory.
 *
 * Loads the dashboard's IIFE modules into one sandbox and hands back their
 * constructors, so each scenario file can open with a single ``require`` and get
 * straight to the behaviour it is describing.
 *
 * The scenario files are run by ``TestTrajectoryScenarioScripts`` in
 * ``tests/test_dashboard_trajectory.py``, which evaluates them in a real browser
 * realm. They also run under Node directly where it is available:
 *   node tests/js/<name>.js
 */
"use strict";

const fs = require("fs");
const pathModule = require("path");
const vm = require("vm");

const dir = pathModule.join(
  __dirname, "..", "..", "cryodrgn", "dashboard", "static", "js"
);
const files = [
  "trajectory_volume_state.js",
  "trajectory_volume_display.js",
  "trajectory_path.js",
  "trajectory_path_mutations.js",
  "trajectory_volume_pipeline.js",
  "trajectory_session.js",
  "trajectory_direct_trace_ui.js"
];

const sandbox = { console, window: {} };
sandbox.window = sandbox;
vm.createContext(sandbox);

for (const f of files) {
  const src = fs.readFileSync(pathModule.join(dir, f), "utf8");
  vm.runInContext(src, sandbox, { filename: f });
}

function assert(cond, msg) {
  if (!cond) throw new Error(msg || "assertion failed");
}

/** Plain object with the DOM listener surface TrajectoryVolumeDisplay wires up. */
function elStub(props) {
  const el = props || {};
  if (typeof el.addEventListener !== "function") el.addEventListener = function () {};
  if (typeof el.removeEventListener !== "function") el.removeEventListener = function () {};
  return el;
}

module.exports = {
  sandbox: sandbox,
  VolumeState: sandbox.CryoTrajectoryVolumeState,
  VolumeStateUtils: sandbox.CryoTrajectoryVolumeStateUtils,
  Session: sandbox.CryoTrajectorySession,
  Mutations: sandbox.CryoTrajectoryPathMutations,
  Path: sandbox.CryoTrajectoryPath,
  Display: sandbox.CryoTrajectoryVolumeDisplay,
  DirectTraceUiState: sandbox.CryoDirectTraceUiState,
  assert: assert,
  elStub: elStub
};
