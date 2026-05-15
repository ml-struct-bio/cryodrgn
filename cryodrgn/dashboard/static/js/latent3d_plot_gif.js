/**
 * GIF export for the 3D Plotly scatter (volume landscape quick actions).
 */
(function(global) {
  "use strict";

  var ORBIT_FRAMES = 48;
  var ORBIT_FRAME_MS = 144;
  var DISCRETE_FRAME_MS = 960;
  var TOIMAGE_SCALE = 2;

  function vecSub(a, b) {
    return { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z };
  }
  function vecAdd(a, b) {
    return { x: a.x + b.x, y: a.y + b.y, z: a.z + b.z };
  }
  function vecScale(v, s) {
    return { x: v.x * s, y: v.y * s, z: v.z * s };
  }
  function vecLen(v) {
    return Math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
  }
  function vecNorm(v) {
    var L = vecLen(v) || 1;
    return { x: v.x / L, y: v.y / L, z: v.z / L };
  }
  function vecCross(a, b) {
    return {
      x: a.y * b.z - a.z * b.y,
      y: a.z * b.x - a.x * b.z,
      z: a.x * b.y - a.y * b.x
    };
  }
  function vecDot(a, b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }

  /** Rodrigues rotation of vector v around unit axis u. */
  function rotateAroundAxis(v, u, theta) {
    var cos = Math.cos(theta);
    var sin = Math.sin(theta);
    var uxv = vecCross(u, v);
    var udotv = vecDot(u, v);
    return vecAdd(vecAdd(vecScale(v, cos), vecScale(uxv, sin)), vecScale(u, udotv * (1 - cos)));
  }

  function plotlyDefaultSceneCamera() {
    return {
      up: { x: 0, y: 0, z: 1 },
      center: { x: 0, y: 0, z: 0 },
      eye: { x: 1.25, y: 1.25, z: 1.25 }
    };
  }

  function normalizeSceneCamera(partial) {
    var d = plotlyDefaultSceneCamera();
    if (!partial || typeof partial !== "object") {
      return JSON.parse(JSON.stringify(d));
    }
    var up = partial.up && typeof partial.up === "object"
      ? vecNorm(partial.up)
      : d.up;
    var center = partial.center && typeof partial.center === "object"
      ? { x: +partial.center.x, y: +partial.center.y, z: +partial.center.z }
      : d.center;
    var eye = partial.eye && typeof partial.eye === "object"
      ? { x: +partial.eye.x, y: +partial.eye.y, z: +partial.eye.z }
      : d.eye;
    return { up: up, center: center, eye: eye };
  }

  /**
   * Resolved camera: Plotly often keeps the active camera only on ``gd._fullLayout.scene``,
   * not on ``gd.layout.scene``, so reading layout alone fails on a fresh draw.
   */
  function resolveSceneCameraFromGd(gd) {
    var raw = null;
    try {
      var flCam = gd._fullLayout && gd._fullLayout.scene && gd._fullLayout.scene.camera;
      if (flCam && typeof flCam === "object") raw = flCam;
    } catch (e0) { /* ignore */ }
    if (!raw) {
      try {
        raw = gd.layout && gd.layout.scene && gd.layout.scene.camera;
      } catch (e1) { /* ignore */ }
    }
    try {
      return JSON.parse(JSON.stringify(normalizeSceneCamera(raw)));
    } catch (e2) {
      return JSON.parse(JSON.stringify(plotlyDefaultSceneCamera()));
    }
  }

  /**
   * Orbit eye around the scene center by rotating (eye - center) about the axis
   * that is the camera "up" projected into the plane perpendicular to the view ray
   * (vertical basis of the focal plane).
   */
  function orbitCameraAroundViewPlaneVertical(cam, theta) {
    var center = cam.center || { x: 0, y: 0, z: 0 };
    var eye = cam.eye || { x: 1.25, y: 1.25, z: 1.25 };
    var upRaw = cam.up || { x: 0, y: 0, z: 1 };
    var up0 = vecNorm(upRaw);
    var offset = vecSub(eye, center);
    var look = vecNorm(vecSub(center, eye));
    var axisProj = vecSub(up0, vecScale(look, vecDot(up0, look)));
    var axisLen = vecLen(axisProj);
    var axis;
    if (axisLen > 1e-6) {
      axis = vecScale(axisProj, 1 / axisLen);
    } else {
      var worldZ = { x: 0, y: 0, z: 1 };
      var ax1 = vecCross(look, worldZ);
      if (vecLen(ax1) < 1e-6) {
        ax1 = vecCross(look, { x: 0, y: 1, z: 0 });
      }
      ax1 = vecNorm(ax1);
      axis = vecNorm(vecCross(look, ax1));
    }
    var rotated = rotateAroundAxis(offset, axis, theta);
    return {
      center: center,
      up: up0,
      eye: vecAdd(center, rotated)
    };
  }

  function parseJsonResponse(r, text) {
    var j;
    try {
      j = JSON.parse(text);
    } catch (e) {
      throw new Error("API did not return JSON");
    }
    if (!r.ok) throw new Error(j.error || r.statusText);
    return j;
  }

  function postAssembleGif(url, frames, durationsMs) {
    var body = { frames: frames, durations_ms: durationsMs };
    return fetch(url, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body)
    }).then(function(r) {
      return r.text().then(function(text) {
        return parseJsonResponse(r, text);
      });
    });
  }

  function dataUrlToBareB64(dataUrl) {
    if (!dataUrl || typeof dataUrl !== "string") return "";
    var i = dataUrl.indexOf(",");
    return i >= 0 ? dataUrl.slice(i + 1) : dataUrl;
  }

  function capturePngDataUrl(gd) {
    var w = gd.offsetWidth || 640;
    var h = gd.offsetHeight || 480;
    return Plotly.toImage(gd, {
      format: "png",
      width: w,
      height: h,
      imageScale: TOIMAGE_SCALE
    });
  }

  function waitFrames(n) {
    return new Promise(function(resolve) {
      var c = 0;
      function tick() {
        c++;
        if (c >= n) resolve();
        else requestAnimationFrame(tick);
      }
      requestAnimationFrame(tick);
    });
  }

  function b64StandardToBlob(b64) {
    var bin = atob(b64);
    var n = bin.length;
    var arr = new Uint8Array(n);
    for (var i = 0; i < n; i++) arr[i] = bin.charCodeAt(i);
    return new Blob([arr], { type: "image/gif" });
  }

  function setBusy(btnOrbitLocal, btnDiscreteLocal, on) {
    if (btnOrbitLocal) btnOrbitLocal.disabled = !!on;
    if (btnDiscreteLocal) btnDiscreteLocal.disabled = !!on;
  }

  function clearBusy(btnOrbitLocal, btnDiscreteLocal) {
    if (btnOrbitLocal) btnOrbitLocal.disabled = false;
    if (btnDiscreteLocal) btnDiscreteLocal.disabled = false;
  }

  function openModal(modal, imgEl, gifB64) {
    var url = "data:image/gif;base64," + gifB64;
    if (imgEl.dataset.objectUrl) {
      try {
        URL.revokeObjectURL(imgEl.dataset.objectUrl);
      } catch (e0) { /* ignore */ }
      delete imgEl.dataset.objectUrl;
    }
    imgEl.src = url;
    modal.hidden = false;
    modal.setAttribute("aria-hidden", "false");
  }

  function closeModal(modal) {
    modal.hidden = true;
    modal.setAttribute("aria-hidden", "true");
  }

  function saveGifWithPicker(blob, suggestedName) {
    if (global.showSaveFilePicker && typeof global.showSaveFilePicker === "function") {
      return global.showSaveFilePicker({
        suggestedName: suggestedName,
        types: [
          {
            description: "GIF",
            accept: { "image/gif": [".gif"] }
          }
        ]
      })
        .then(function(handle) {
          return handle.createWritable();
        })
        .then(function(writable) {
          return writable.write(blob).then(function() {
            return writable.close();
          });
        });
    }
    var a = document.createElement("a");
    var u = URL.createObjectURL(blob);
    a.href = u;
    a.download = suggestedName;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(function() {
      try {
        URL.revokeObjectURL(u);
      } catch (e1) { /* ignore */ }
    }, 2500);
    return Promise.resolve();
  }

  function boot(opts) {
    var gd = opts.gd;
    var assembleUrl = opts.assembleUrl;
    var discreteGifServerUrl = opts.discreteGifServerUrl || null;
    var sx = opts.sx;
    var sy = opts.sy;
    var sz = opts.sz;
    var sc = opts.sc;
    var covariateDisplayMap = opts.covariateDisplayMap || {};
    var buildPayload = opts.buildPayload;
    var statusEl = opts.statusEl;

    var btnOrbit = document.getElementById("l3d-gif-orbit-360");
    var btnDiscrete = document.getElementById("l3d-gif-discrete-levels");
    var modal = document.getElementById("l3d-plot-gif-modal");
    var imgEl = document.getElementById("l3d-plot-gif-modal-preview");
    var btnSave = document.getElementById("l3d-plot-gif-modal-save");
    var discreteGifProgressEl = document.getElementById("l3d-discrete-gif-progress");

    var lastGifBlob = null;
    var lastSuggestedName = "cryodrgn_landscape_3d.gif";

    function axesOk() {
      var x = sx.value;
      var y = sy.value;
      var z = sz.value;
      return x !== y && x !== z && y !== z;
    }

    function flashStatus(msg) {
      if (!statusEl) return;
      var prev = statusEl.textContent;
      statusEl.textContent = msg;
      if (msg) {
        setTimeout(function() {
          if (statusEl.textContent === msg) statusEl.textContent = prev;
        }, 5200);
      }
    }

    function setDiscreteGifProgressVisible(visible) {
      if (!discreteGifProgressEl) return;
      discreteGifProgressEl.hidden = !visible;
      discreteGifProgressEl.setAttribute("aria-busy", visible ? "true" : "false");
    }

    function discreteCovariateDisplayName() {
      var col = sc.value;
      if (!col || col === "none") return "";
      if (Object.prototype.hasOwnProperty.call(covariateDisplayMap, col)) {
        return String(covariateDisplayMap[col]);
      }
      return String(col);
    }

    function syncDiscreteButton() {
      if (!btnDiscrete) return;
      var ok = false;
      var legendTitle = "";
      try {
        var m = gd.layout && gd.layout.meta;
        ok = !!(
          m
          && m.cdrgn_color_mode === "discrete"
          && sc.value
          && sc.value !== "none"
        );
        if (ok && m.cdrgn_color_legend && m.cdrgn_color_legend.title) {
          legendTitle = String(m.cdrgn_color_legend.title);
        }
      } catch (e) {
        ok = false;
      }
      btnDiscrete.disabled = !ok;
      var mid;
      if (!ok) {
        mid = "colour covariate";
      } else {
        mid = legendTitle || discreteCovariateDisplayName() || "colour covariate";
      }
      btnDiscrete.textContent = "Create\n" + mid + "\nanimation";
      var tip = ok
        ? (
          "Steps through each level of " + mid.replace(/"/g, "'")
          + " on the server (parallel Matplotlib); axis limits match the plot; view is the default 3D angle."
        )
        : "Choose a discrete colour covariate to enable this animation.";
      btnDiscrete.setAttribute("title", tip);
    }

    function runOrbitGif() {
      if (!axesOk()) {
        flashStatus("Choose three different axes before exporting a GIF.");
        return;
      }
      if (typeof Plotly === "undefined" || !Plotly.relayout || !Plotly.toImage) return;
      if (!gd.layout || !gd.layout.scene) {
        flashStatus("Wait until the 3D plot has finished loading, then try again.");
        return;
      }
      var cam0 = resolveSceneCameraFromGd(gd);
      setBusy(btnOrbit, btnDiscrete, true);
      var frames = [];
      var durs = [];
      var saved = cam0;
      var chain = Promise.resolve();
      for (var i = 0; i < ORBIT_FRAMES; i++) {
        (function(ii) {
          chain = chain.then(function() {
            var theta = (2 * Math.PI * ii) / ORBIT_FRAMES;
            var nextCam = orbitCameraAroundViewPlaneVertical(saved, theta);
            return Plotly.relayout(gd, { "scene.camera": nextCam });
          })
            .then(function() {
              return waitFrames(2);
            })
            .then(function() {
              return capturePngDataUrl(gd);
            })
            .then(function(dataUrl) {
              frames.push(dataUrlToBareB64(dataUrl));
              durs.push(ORBIT_FRAME_MS);
            });
        })(i);
      }
      chain = chain
        .then(function() {
          return Plotly.relayout(gd, { "scene.camera": saved });
        })
        .then(function() {
          return postAssembleGif(assembleUrl, frames, durs);
        })
        .then(function(res) {
          if (!res.gif_b64) throw new Error("Server did not return gif_b64.");
          lastGifBlob = b64StandardToBlob(res.gif_b64);
          lastSuggestedName = "cryodrgn_landscape_3d_orbit.gif";
          if (modal && imgEl) openModal(modal, imgEl, res.gif_b64);
        })
        .catch(function(err) {
          console.error(err);
          flashStatus(err.message || "Could not build orbit GIF.");
          return Plotly.relayout(gd, { "scene.camera": saved }).catch(function() {});
        })
        .then(function() {
          clearBusy(btnOrbit, btnDiscrete);
          syncDiscreteButton();
        });
    }

    function discreteKeysFromMeta() {
      try {
        var leg = gd.layout.meta.cdrgn_color_legend;
        if (!leg || leg.type !== "discrete" || !leg.items) return [];
        return leg.items.map(function(it) {
          return it.label;
        });
      } catch (e2) {
        return [];
      }
    }

    function runDiscreteLevelsGif() {
      if (!btnDiscrete || btnDiscrete.disabled) return;
      if (!axesOk()) {
        flashStatus("Choose three different axes before exporting a GIF.");
        return;
      }
      if (!discreteGifServerUrl) {
        flashStatus("Discrete-level GIF export is not available on this view.");
        return;
      }
      var keys = discreteKeysFromMeta();
      if (!keys.length) {
        flashStatus("No discrete colour levels to animate.");
        return;
      }
      setBusy(btnOrbit, btnDiscrete, true);
      setDiscreteGifProgressVisible(true);
      var payload = buildPayload({
        discrete_keys: keys,
        frame_duration_ms: DISCRETE_FRAME_MS
      });
      fetch(discreteGifServerUrl, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload)
      })
        .then(function(r) {
          return r.text().then(function(text) {
            return parseJsonResponse(r, text);
          });
        })
        .then(function(res) {
          if (!res.gif_b64) throw new Error("Server did not return gif_b64.");
          lastGifBlob = b64StandardToBlob(res.gif_b64);
          lastSuggestedName = "cryodrgn_landscape_3d_discrete_levels.gif";
          if (modal && imgEl) openModal(modal, imgEl, res.gif_b64);
        })
        .catch(function(err) {
          console.error(err);
          flashStatus(err.message || "Could not build discrete-level GIF.");
        })
        .then(function() {
          clearBusy(btnOrbit, btnDiscrete);
          syncDiscreteButton();
        })
        .finally(function() {
          setDiscreteGifProgressVisible(false);
        });
    }

    if (btnOrbit) btnOrbit.addEventListener("click", runOrbitGif);
    if (btnDiscrete) btnDiscrete.addEventListener("click", runDiscreteLevelsGif);

    if (modal) {
      modal.addEventListener("click", function(ev) {
        if (ev.target && ev.target.getAttribute("data-l3d-plot-gif-close")) {
          closeModal(modal);
        }
      });
    }
    if (btnSave) {
      btnSave.addEventListener("click", function() {
        if (!lastGifBlob) return;
        saveGifWithPicker(lastGifBlob, lastSuggestedName).catch(function(err) {
          console.error(err);
          flashStatus(err.message || "Save was cancelled or failed.");
        });
      });
    }

    global.addEventListener("keydown", function(ev) {
      if (ev.key === "Escape" && modal && !modal.hidden) {
        closeModal(modal);
      }
    });

    return { syncUi: syncDiscreteButton };
  }

  global.CryoLatent3dPlotGif = { boot: boot };
})(typeof window !== "undefined" ? window : this);
