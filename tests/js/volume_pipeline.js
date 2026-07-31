/**
 * Decode / render job orchestration: which slots a job asks for, and what happens
 * to a response that arrives after the path has moved on.
 *
 * The pipeline's fixtures are asynchronous fetch callbacks, so this stays a
 * scenario script. Assertions run inside promises, so the file exports the chain
 * and the harness awaits it.
 *
 * Run by ``TestTrajectoryScenarioScripts`` in ``tests/test_dashboard_trajectory.py``.
 */
"use strict";

const { VolumeState, assert } = require("./_harness.js");

const Pipeline = require("./_harness.js").sandbox.CryoTrajectoryVolumePipeline;
assert(Pipeline, "CryoTrajectoryVolumePipeline missing from the module set");

/** Four slots: 0 and 3 decoded and rendered, 1 and 2 owing both. */
function debtState() {
  const st = new VolumeState();
  st.replaceSlots(
    ["a", "b", "c", "d"],
    [{ volume_b64: "A" }, null, null, { volume_b64: "D" }],
    ["imgA", null, null, "imgD"]
  );
  return st;
}

function pipelineWith(hooks, volumes) {
  // ``arguments.length`` rather than a default, so a scenario can pass null on
  // purpose to exercise the no-volume-state guard.
  return new Pipeline({
    volumes: arguments.length > 1 ? volumes : debtState(),
    hooks: hooks || {}
  });
}

const checks = [];

// A pipeline with nothing to talk to reports why, rather than throwing.
checks.push(
  pipelineWith({}, null).decode({}).then(function (res) {
    assert(res.ok === false && res.reason === "no-volumes",
      "decode without a volume state reports no-volumes");
  })
);
checks.push(
  pipelineWith({}).decode({}).then(function (res) {
    assert(res.ok === false && res.reason === "no-fetch",
      "decode without a fetch hook reports no-fetch");
  })
);

// Default batch is the decode debt, not every slot.
let decodeAsked = null;
checks.push(
  pipelineWith({
    fetchDecode: function (req) {
      decodeAsked = req;
      return Promise.resolve({ volumes: [] });
    }
  }).decode({}).then(function () {
    assert(decodeAsked.indices.join(",") === "1,2",
      "decode asks for the debt slots only");
    assert(decodeAsked.ids.join(",") === "b,c",
      "decode sends the slot ids alongside the indices");
    assert(decodeAsked.forceAll === false, "forceAll defaults off");
  })
);

// forceAll widens the batch to the whole path.
let forcedAsked = null;
checks.push(
  pipelineWith({
    fetchDecode: function (req) {
      forcedAsked = req;
      return Promise.resolve({ volumes: [] });
    }
  }).decode({ forceAll: true }).then(function () {
    assert(forcedAsked.indices.join(",") === "0,1,2,3",
      "forceAll decodes every slot");
    assert(forcedAsked.forceAll === true, "forceAll is passed to the fetch hook");
  })
);

// An explicit index list wins over the computed debt.
let explicitAsked = null;
checks.push(
  pipelineWith({
    fetchDecode: function (req) {
      explicitAsked = req;
      return Promise.resolve({ volumes: [] });
    }
  }).decode({ indices: [3] }).then(function () {
    assert(explicitAsked.indices.join(",") === "3",
      "an explicit index list is used as given");
  })
);

// A response from a superseded generation must not be applied: a slower first
// decode landing after a newer one would otherwise overwrite fresher volumes.
const staleState = debtState();
let releaseFirst = null;
const firstDecode = new Promise(function (resolve) { releaseFirst = resolve; });
const stalePipeline = pipelineWith(
  {
    fetchDecode: function (req) {
      return req.stallMe ? firstDecode : Promise.resolve({ volumes: [] });
    }
  },
  staleState
);
const slow = stalePipeline.decode({ indices: [1], stallMe: true });
checks.push(
  stalePipeline.decode({ indices: [2] }).then(function (fresh) {
    assert(fresh.ok !== false || fresh.reason !== "stale",
      "the newest decode is the live one");
    releaseFirst({ volumes: [{ index: 1, volume_b64: "LATE" }] });
    return slow;
  }).then(function (res) {
    assert(res.ok === false && res.reason === "stale",
      "a superseded decode reports itself stale");
    assert(!staleState.isDecoded(1),
      "a stale decode must not write its volume into the state");
  })
);

// Render mirrors decode: the default batch is the inactive ChimeraX ticks, and
// forceAll widens it to the whole path. Each scenario captures into its own
// variable, since the checks run concurrently.
let renderDefaultAsked = null;
checks.push(
  pipelineWith({
    fetchRender: function (req) {
      renderDefaultAsked = req;
      return Promise.resolve({ images: [] });
    }
  }).render({}).then(function () {
    assert(renderDefaultAsked.indices.join(",") === "1,2",
      "render asks for the inactive ticks only");
    assert(renderDefaultAsked.ids.join(",") === "b,c",
      "render sends the slot ids alongside the indices");
  })
);

let renderForcedAsked = null;
checks.push(
  pipelineWith({
    fetchRender: function (req) {
      renderForcedAsked = req;
      return Promise.resolve({ images: [] });
    }
  }).render({ forceAll: true }).then(function () {
    assert(renderForcedAsked.forceAll === true, "render forwards forceAll");
    assert(renderForcedAsked.indices.join(",") === "0,1,2,3",
      "forceAll re-renders every slot, not just the inactive ticks");
  })
);

// An empty index list is treated as "no preference", not as "render nothing".
let renderEmptyAsked = null;
checks.push(
  pipelineWith({
    fetchRender: function (req) {
      renderEmptyAsked = req;
      return Promise.resolve({ images: [] });
    }
  }).render({ indices: [] }).then(function () {
    assert(renderEmptyAsked.indices.join(",") === "1,2",
      "an empty index list falls back to the inactive ticks");
  })
);

// Busy flags bracket a job so the page can gate its controls.
const busyState = debtState();
let releaseBusy = null;
const busyFetch = new Promise(function (resolve) { releaseBusy = resolve; });
const busyPipeline = pipelineWith(
  { fetchDecode: function () { return busyFetch; } },
  busyState
);
assert(!busyPipeline.isBusy(), "an idle pipeline is not busy");
const busyRun = busyPipeline.decode({});
assert(busyPipeline.isBusy() && busyPipeline.isDecoding(),
  "a pipeline with a decode in flight reports busy");
releaseBusy({ volumes: [] });
checks.push(
  busyRun.then(function () {
    assert(!busyPipeline.isBusy(), "the busy flag clears once the job settles");
  })
);

module.exports = Promise.all(checks).then(function () {
  console.log("volume_pipeline: ok");
});
