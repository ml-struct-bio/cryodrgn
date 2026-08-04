/**
 * Shared helpers for landscape volPCA and latent-3D volume animations.
 *
 * Pure utilities live on ``CryoLandscapeVolShared`` directly. DOM-coupled
 * view / GIF chrome is built with ``createViewGifControls(ctx)``. Overlay
 * GIF clocks use ``attachOverlayClock(wrap, opts)``.
 */
(function (global) {
  "use strict";

  function postJson(url, body) {
    return fetch(url, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body || {}),
    }).then(function (r) {
      return r.json().then(function (j) {
        if (!r.ok) throw new Error(j.error || r.status);
        return j;
      });
    });
  }

  function plotlyTraceLength(trace, key) {
    var PLOTLY = global.CryoPlotlyArrays;
    if (!trace || !trace[key] || !PLOTLY) return 0;
    return PLOTLY.length(trace[key]);
  }

  function normalizeViewDegrees(degrees) {
    var d = Number(degrees);
    if (!isFinite(d)) return 0;
    d = d % 360;
    if (Math.abs(d) < 1e-9) return 0;
    return d;
  }

  function mat3Identity() {
    return [1, 0, 0, 0, 1, 0, 0, 0, 1];
  }

  function mat3Mul(a, b) {
    return [
      a[0] * b[0] + a[1] * b[3] + a[2] * b[6],
      a[0] * b[1] + a[1] * b[4] + a[2] * b[7],
      a[0] * b[2] + a[1] * b[5] + a[2] * b[8],
      a[3] * b[0] + a[4] * b[3] + a[5] * b[6],
      a[3] * b[1] + a[4] * b[4] + a[5] * b[7],
      a[3] * b[2] + a[4] * b[5] + a[5] * b[8],
      a[6] * b[0] + a[7] * b[3] + a[8] * b[6],
      a[6] * b[1] + a[7] * b[4] + a[8] * b[7],
      a[6] * b[2] + a[7] * b[5] + a[8] * b[8],
    ];
  }

  function mat3ForAxisTurn(axis, degrees) {
    var rad = (degrees * Math.PI) / 180;
    var c = Math.cos(rad);
    var s = Math.sin(rad);
    if (axis === "x") {
      return [1, 0, 0, 0, c, -s, 0, s, c];
    }
    if (axis === "y") {
      return [c, 0, s, 0, 1, 0, -s, 0, c];
    }
    return [c, -s, 0, s, c, 0, 0, 0, 1];
  }

  function shufflePick(arr, k) {
    var a = arr.slice();
    for (var i = a.length - 1; i > 0; i--) {
      var j = Math.floor(Math.random() * (i + 1));
      var t = a[i];
      a[i] = a[j];
      a[j] = t;
    }
    return a.slice(0, k);
  }

  /** A–Z omitting I, O, U — matches particle explorer / montage_cell_label. */
  var SAFE_LETTERS = (function () {
    var out = [];
    for (var c = 65; c <= 90; c++) {
      if (c !== 73 && c !== 79 && c !== 85) out.push(String.fromCharCode(c));
    }
    return out;
  })();

  function montageLabelAt(idx) {
    if (idx < SAFE_LETTERS.length) return SAFE_LETTERS[idx];
    var j = idx - SAFE_LETTERS.length;
    return SAFE_LETTERS[Math.floor(j / SAFE_LETTERS.length)]
      + SAFE_LETTERS[j % SAFE_LETTERS.length];
  }

  function validateViewMatrixText(text) {
    var raw = String(text || "").trim();
    if (!raw) {
      return { ok: false, msg: "Enter a ChimeraX view matrix (12 numbers)." };
    }
    var body = raw.toLowerCase().indexOf("camera") === 0 ? raw.slice(6).trim() : raw;
    // Allow integers, decimals (``.5`` or ``1.``), and scientific notation.
    var nums = body.match(/[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?/g);
    if (!nums || nums.length < 12) {
      return {
        ok: false,
        msg: "View matrix must contain 12 numbers (optionally prefixed with camera).",
      };
    }
    return { ok: true, text: raw };
  }

  /**
   * Attach a requestAnimationFrame overlay clock to a preview wrap.
   *
   * opts:
   *   frameDurationMs, totalFrames, applyFrameIdx(idx),
   *   indexFromFrame(frameIdx) -> idx,
   *   cancelProp (e.g. "_volsketchCancelOverlay")
   */
  function attachOverlayClock(wrap, opts) {
    opts = opts || {};
    if (!wrap) return null;
    var fd = Math.max(1, Math.floor(Number(opts.frameDurationMs)) || 1);
    var total = Math.max(1, Math.floor(Number(opts.totalFrames)) || 1);
    var applyFrameIdx = opts.applyFrameIdx;
    var indexFromFrame = opts.indexFromFrame || function (frameIdx) {
      return frameIdx;
    };
    var cancelProp = opts.cancelProp || "_cryoOverlayCancel";
    var rafId = 0;
    var bootRaf = 0;
    var startT = null;

    function tick(now) {
      if (!wrap.isConnected) {
        if (rafId) cancelAnimationFrame(rafId);
        rafId = 0;
        return;
      }
      if (startT === null) startT = now;
      var frameIdx = Math.floor((now - startT) / fd) % total;
      applyFrameIdx(indexFromFrame(frameIdx));
      rafId = requestAnimationFrame(tick);
    }

    function begin() {
      applyFrameIdx(indexFromFrame(0));
      startT = null;
      bootRaf = requestAnimationFrame(function () {
        bootRaf = 0;
        rafId = requestAnimationFrame(tick);
      });
    }

    function cancel() {
      if (bootRaf) cancelAnimationFrame(bootRaf);
      if (rafId) cancelAnimationFrame(rafId);
      bootRaf = 0;
      rafId = 0;
    }

    wrap[cancelProp] = cancel;

    return {
      begin: begin,
      cancel: cancel,
      startOnImage: function (img) {
        var dec = img && img.decode && img.decode();
        if (dec && typeof dec.then === "function") {
          dec.then(begin).catch(begin);
        } else if (img && img.complete) {
          begin();
        } else if (img) {
          img.addEventListener("load", begin, { once: true });
        } else {
          begin();
        }
      },
    };
  }

  /**
   * Factory for shared view-rotation / GIF-save chrome.
   *
   * ctx.state (mutable bag):
   *   viewRotations, lastChimeraxViewMatrix, chimeraxViewMatrixUnavailable,
   *   appliedViewMatrix, viewMatrixInputDirty, viewMatrixFieldFocused,
   *   landscapeAnimInFlight, lastAnimToken, lastBatchMode
   *
   * ctx.els: animateStatusEl, animateProgressEl, saveGifBtn, randomSelBtn,
   *   previewGrid, viewRotateRow, viewRotateAngleEl, viewMatrixInputEl,
   *   viewMatrixApplyBtn, viewRotateBtns (array)
   *
   * ctx.hooks:
   *   getGifMode(), getSelectedVols(), getMeta(),
   *   clearPreviewGrid(), cancelPendingGifWork(),
   *   scheduleAutoGif(reason), syncAnimOutputControls()
   */
  function createViewGifControls(ctx) {
    ctx = ctx || {};
    var state = ctx.state || {};
    var els = ctx.els || {};
    var hooks = ctx.hooks || {};

    function animationsEnabled() {
      return (typeof hooks.getGifMode === "function" ? hooks.getGifMode() : "disabled")
        !== "disabled";
    }

    function selectedSize() {
      var sel = typeof hooks.getSelectedVols === "function" ? hooks.getSelectedVols() : null;
      return sel && sel.size ? sel.size : 0;
    }

    function clearAppliedViewMatrix() {
      state.appliedViewMatrix = "";
    }

    function currentViewRotationPayload() {
      var vr = state.viewRotations || { x: 0, y: 0, z: 0 };
      return {
        x: normalizeViewDegrees(vr.x),
        y: normalizeViewDegrees(vr.y),
        z: normalizeViewDegrees(vr.z),
      };
    }

    function viewRotationsAreActive() {
      var r = currentViewRotationPayload();
      return Math.abs(r.x) > 1e-9 || Math.abs(r.y) > 1e-9 || Math.abs(r.z) > 1e-9;
    }

    function estimatedChimeraxViewMatrixText() {
      var vr = state.viewRotations || { x: 0, y: 0, z: 0 };
      var m = mat3Identity();
      ["x", "y", "z"].forEach(function (axis) {
        var deg = normalizeViewDegrees(vr[axis] || 0);
        if (Math.abs(deg) > 1e-9) {
          m = mat3Mul(m, mat3ForAxisTurn(axis, deg));
        }
      });
      var nums = [
        m[0], m[1], m[2], 0,
        m[3], m[4], m[5], 0,
        m[6], m[7], m[8], 0,
      ];
      return "camera " + nums.map(function (n) {
        return Number(n).toFixed(6);
      }).join(",");
    }

    function viewRotationSummary() {
      var r = currentViewRotationPayload();
      return "requested view turns X " + r.x.toFixed(1)
        + "°, Y " + r.y.toFixed(1)
        + "°, Z " + r.z.toFixed(1) + "°";
    }

    function currentViewMatrixDisplayText() {
      if (state.lastChimeraxViewMatrix) return state.lastChimeraxViewMatrix;
      if (viewRotationsAreActive()) return estimatedChimeraxViewMatrixText();
      return "";
    }

    function setAnimateStatus(msg, showProgress, isErr) {
      var animateStatusEl = els.animateStatusEl;
      var animateProgressEl = els.animateProgressEl;
      if (!animateStatusEl) return;
      animateStatusEl.textContent = msg || "";
      animateStatusEl.style.color = isErr ? "var(--error, #b42318)" : "";
      if (animateProgressEl) {
        animateProgressEl.hidden = !showProgress || !!isErr;
        if (showProgress && !isErr) {
          animateProgressEl.setAttribute("aria-busy", "true");
        } else {
          animateProgressEl.removeAttribute("aria-busy");
        }
      }
    }

    function showAnimBusySelectionMsg() {
      setAnimateStatus(
        "Wait for ChimeraX to finish the current animation preview before changing the selection.",
        false,
        true
      );
    }

    function syncViewMatrixField() {
      var viewMatrixInputEl = els.viewMatrixInputEl;
      if (!viewMatrixInputEl) return;
      if (state.viewMatrixFieldFocused || state.viewMatrixInputDirty) return;
      var text = currentViewMatrixDisplayText();
      if (state.chimeraxViewMatrixUnavailable && !text) {
        viewMatrixInputEl.placeholder = "ChimeraX view matrix not reported for this render.";
      } else if (!text) {
        viewMatrixInputEl.placeholder = (
          "camera n1,n2,... (12 numbers; available after animation loads)"
        );
      } else {
        viewMatrixInputEl.placeholder = "";
      }
      viewMatrixInputEl.value = text;
    }

    function syncViewRotationControls() {
      var viewRotateRow = els.viewRotateRow;
      var viewRotateAngleEl = els.viewRotateAngleEl;
      var viewRotateBtns = els.viewRotateBtns || [];
      var viewMatrixApplyBtn = els.viewMatrixApplyBtn;
      var viewMatrixInputEl = els.viewMatrixInputEl;
      var previewGrid = els.previewGrid;
      if (!animationsEnabled()) {
        var offHint = "Enable Cycle or Rotate animation output to render GIF previews.";
        if (viewRotateRow) {
          viewRotateRow.classList.add("is-disabled");
          viewRotateRow.title = offHint;
        }
        if (viewRotateAngleEl) {
          viewRotateAngleEl.disabled = true;
          viewRotateAngleEl.title = offHint;
        }
        viewRotateBtns.forEach(function (btn) {
          btn.disabled = true;
          btn.title = offHint;
        });
        if (viewMatrixApplyBtn) {
          viewMatrixApplyBtn.disabled = true;
          viewMatrixApplyBtn.title = offHint;
        }
        if (viewMatrixInputEl) {
          viewMatrixInputEl.disabled = true;
          viewMatrixInputEl.title = offHint;
        }
        return;
      }
      var hasLoadedPreview = !!(
        state.lastAnimToken
        && previewGrid
        && previewGrid.children.length > 0
        && !state.landscapeAnimInFlight
      );
      var hint = hasLoadedPreview
        ? "Rotate the loaded animation preview by the entered angle."
        : "Load an animation preview before rotating the view.";
      if (viewRotateRow) {
        viewRotateRow.classList.toggle("is-disabled", !hasLoadedPreview);
        viewRotateRow.title = hint;
      }
      if (viewRotateAngleEl) {
        viewRotateAngleEl.disabled = !hasLoadedPreview;
        viewRotateAngleEl.title = hint;
      }
      viewRotateBtns.forEach(function (btn) {
        btn.disabled = !hasLoadedPreview;
        btn.title = hint;
      });
      if (viewMatrixApplyBtn) {
        viewMatrixApplyBtn.disabled = !hasLoadedPreview;
        viewMatrixApplyBtn.title = hasLoadedPreview
          ? "Re-render previews using the view matrix in the field."
          : hint;
      }
      if (viewMatrixInputEl) {
        viewMatrixInputEl.disabled = !hasLoadedPreview;
        viewMatrixInputEl.title = hint;
      }
    }

    function syncSaveGifButton() {
      var saveGifBtn = els.saveGifBtn;
      var previewGrid = els.previewGrid;
      if (!saveGifBtn) return;
      if (!animationsEnabled()) {
        saveGifBtn.disabled = true;
        saveGifBtn.title = "Enable Cycle or Rotate animation output to render GIF previews.";
        if (typeof hooks.syncAnimOutputControls === "function") {
          hooks.syncAnimOutputControls();
        }
        return;
      }
      var hasPreviews = previewGrid && previewGrid.children.length > 0;
      var canSave = !!state.lastAnimToken && hasPreviews && !state.landscapeAnimInFlight;
      saveGifBtn.disabled = !canSave;
      if (canSave) {
        saveGifBtn.removeAttribute("title");
      } else if (state.landscapeAnimInFlight) {
        saveGifBtn.title = "Available after ChimeraX finishes rendering the current previews.";
      } else {
        saveGifBtn.title = "Select volumes and wait for GIF previews to finish rendering.";
      }
      if (typeof hooks.syncAnimOutputControls === "function") {
        hooks.syncAnimOutputControls();
      }
    }

    function clearAnimPreviewState() {
      if (typeof hooks.clearPreviewGrid === "function") {
        hooks.clearPreviewGrid();
      }
      state.lastAnimToken = null;
      state.lastBatchMode = null;
      state.lastChimeraxViewMatrix = "";
      state.chimeraxViewMatrixUnavailable = false;
      clearAppliedViewMatrix();
      state.viewMatrixInputDirty = false;
      syncViewMatrixField();
    }

    function finishSelectionWithoutAnimation(reason) {
      if (typeof hooks.cancelPendingGifWork === "function") {
        hooks.cancelPendingGifWork();
      }
      if (!selectedSize()) {
        clearAnimPreviewState();
        setAnimateStatus("", false, false);
      } else if (reason === "gif_mode") {
        clearAnimPreviewState();
        setAnimateStatus(
          "Animation disabled — volume selection still updates the plot.",
          false,
          false
        );
      }
      syncSaveGifButton();
      if (typeof hooks.syncAnimOutputControls === "function") {
        hooks.syncAnimOutputControls();
      }
    }

    function maxSelectable() {
      var META = typeof hooks.getMeta === "function" ? hooks.getMeta() : null;
      return META && META.n_volumes ? META.n_volumes : 10000;
    }

    function syncRandomSelButton() {
      var randomSelBtn = els.randomSelBtn;
      var META = typeof hooks.getMeta === "function" ? hooks.getMeta() : null;
      if (!randomSelBtn || !META || META.chimerax_cpus == null) return;
      var n = Math.max(1, Number(META.chimerax_cpus));
      randomSelBtn.textContent = "Choose " + n + " volume" + (n === 1 ? "" : "s") + " at random";
      randomSelBtn.disabled = false;
    }

    function applyViewMatrixFromField() {
      var viewMatrixInputEl = els.viewMatrixInputEl;
      if (!animationsEnabled()) {
        setAnimateStatus(
          "Enable Cycle or Rotate animation output to render GIF previews.",
          false,
          true
        );
        return;
      }
      if (!viewMatrixInputEl) return;
      var check = validateViewMatrixText(viewMatrixInputEl.value);
      if (!check.ok) {
        setAnimateStatus(check.msg, false, true);
        return;
      }
      state.appliedViewMatrix = check.text;
      state.viewMatrixInputDirty = false;
      syncViewMatrixField();
      setAnimateStatus("Applying custom ChimeraX view matrix.", false, false);
      if (selectedSize() > 0) {
        if (typeof hooks.scheduleAutoGif === "function") {
          hooks.scheduleAutoGif("view_matrix");
        }
      } else {
        syncSaveGifButton();
      }
    }

    function applyViewRotation(axis) {
      var viewRotateAngleEl = els.viewRotateAngleEl;
      if (!animationsEnabled()) {
        setAnimateStatus(
          "Enable Cycle or Rotate animation output to render GIF previews.",
          false,
          true
        );
        return;
      }
      var deg = viewRotateAngleEl ? Number(viewRotateAngleEl.value) : NaN;
      if (!isFinite(deg)) {
        setAnimateStatus("Enter a finite rotation angle in degrees.", false, true);
        return;
      }
      clearAppliedViewMatrix();
      state.viewMatrixInputDirty = false;
      if (!state.viewRotations) state.viewRotations = { x: 0, y: 0, z: 0 };
      state.viewRotations[axis] = normalizeViewDegrees(
        Number(state.viewRotations[axis] || 0) + deg
      );
      syncViewMatrixField();
      setAnimateStatus("Updated " + viewRotationSummary() + ".", false, false);
      if (selectedSize() > 0) {
        if (typeof hooks.scheduleAutoGif === "function") {
          hooks.scheduleAutoGif("view_rotation");
        }
      } else {
        syncSaveGifButton();
      }
    }

    return {
      animationsEnabled: animationsEnabled,
      clearAppliedViewMatrix: clearAppliedViewMatrix,
      currentViewRotationPayload: currentViewRotationPayload,
      viewRotationsAreActive: viewRotationsAreActive,
      estimatedChimeraxViewMatrixText: estimatedChimeraxViewMatrixText,
      viewRotationSummary: viewRotationSummary,
      currentViewMatrixDisplayText: currentViewMatrixDisplayText,
      setAnimateStatus: setAnimateStatus,
      showAnimBusySelectionMsg: showAnimBusySelectionMsg,
      syncViewMatrixField: syncViewMatrixField,
      syncViewRotationControls: syncViewRotationControls,
      syncSaveGifButton: syncSaveGifButton,
      clearAnimPreviewState: clearAnimPreviewState,
      finishSelectionWithoutAnimation: finishSelectionWithoutAnimation,
      maxSelectable: maxSelectable,
      syncRandomSelButton: syncRandomSelButton,
      applyViewMatrixFromField: applyViewMatrixFromField,
      applyViewRotation: applyViewRotation,
      validateViewMatrixText: validateViewMatrixText,
    };
  }

  global.CryoLandscapeVolShared = {
    postJson: postJson,
    plotlyTraceLength: plotlyTraceLength,
    normalizeViewDegrees: normalizeViewDegrees,
    mat3Identity: mat3Identity,
    mat3Mul: mat3Mul,
    mat3ForAxisTurn: mat3ForAxisTurn,
    shufflePick: shufflePick,
    SAFE_LETTERS: SAFE_LETTERS,
    montageLabelAt: montageLabelAt,
    validateViewMatrixText: validateViewMatrixText,
    attachOverlayClock: attachOverlayClock,
    createViewGifControls: createViewGifControls,
  };
})(typeof window !== "undefined" ? window : this);
