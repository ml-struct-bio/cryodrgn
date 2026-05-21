(function() {
  var gd = document.getElementById("volsketch");
  var plotStatusEl = document.getElementById("volsketch-plot-status");
  var animateStatusEl = document.getElementById("volsketch-animate-status");
  var animateProgressEl = document.getElementById("volsketch-animate-progress");
  var selSummaryEl = document.getElementById("volsketch-sel-summary");
  var selIndicesDetails = document.getElementById("volsketch-sel-indices-details");
  var selIndicesSummaryLabel = document.getElementById("volsketch-sel-indices-summary-label");
  var selIndicesList = document.getElementById("volsketch-sel-indices-list");
  var pcx = document.getElementById("volsketch-pcx");
  var pcy = document.getElementById("volsketch-pcy");
  var colorSel = document.getElementById("volsketch-color");
  var gifModeDescEl = document.getElementById("volsketch-mode-desc");
  var gifFrames = document.getElementById("volsketch-gif-frames");
  var gifFramesWrap = document.getElementById("volsketch-gif-frames-wrap");
  var previewGrid = document.getElementById("volsketch-preview-grid");
  var randomSelBtn = document.getElementById("volsketch-random-sel");
  var saveGifBtn = document.getElementById("volsketch-save-gif");
  var viewRotateRow = document.querySelector(".volsketch-view-rotate-row");
  var viewRotateAngleEl = document.getElementById("volsketch-view-rotate-angle");
  var viewMatrixInputEl = document.getElementById("volsketch-view-matrix-input");
  var viewMatrixApplyBtn = document.getElementById("volsketch-view-matrix-apply");
  var viewRotateBtns = Array.from(document.querySelectorAll("[data-volsketch-view-axis]"));

  function teardownVolsketchPreviews() {
    if (!previewGrid) return;
    previewGrid.querySelectorAll(".volsketch-anim-wrap").forEach(function(w) {
      if (w._volsketchCancelOverlay) {
        w._volsketchCancelOverlay();
        w._volsketchCancelOverlay = null;
      }
    });
  }

  function clearVolsketchPreviewGrid() {
    teardownVolsketchPreviews();
    if (previewGrid) previewGrid.innerHTML = "";
    syncViewRotationControls();
  }

  function volsketchHex6ToRgba(hex, a) {
    var m = /^#([0-9a-fA-F]{6})$/.exec(hex);
    if (!m) return null;
    var x = parseInt(m[1], 16);
    var r = (x >> 16) & 255;
    var g = (x >> 8) & 255;
    var b = x & 255;
    return "rgba(" + r + "," + g + "," + b + "," + a + ")";
  }

  function volsketchIsCircledMontageGlyphText(text) {
    var s = text == null ? "" : String(text);
    if (!s) return false;
    var code = s.charCodeAt(0);
    return code >= 0x24b6 && code <= 0x24cf;
  }

  /** Volume colour → circled glyph tint or round letter chip; covariate badges stay rectangular. */
  function volsketchApplyBadgeBackground(badge, color) {
    badge.style.color = "#1a1a1a";
    badge.style.backgroundColor = "";
    badge.style.borderColor = "";
    badge.classList.remove("volsketch-gif-overlay--on-tint");
    badge.classList.remove("volsketch-gif-overlay--circle-chip");
    if (typeof color !== "string" || !color.trim()) {
      return;
    }
    var c = color.trim();
    var rgba = null;
    if (/^#([0-9a-fA-F]{6})$/.test(c)) {
      rgba = volsketchHex6ToRgba(c, 0.88);
    } else {
      var rm = /^rgba?\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)/i.exec(c);
      if (rm) rgba = "rgba(" + rm[1] + "," + rm[2] + "," + rm[3] + ",0.88)";
    }
    if (!rgba) return;
    var label = badge.textContent || "";
    if (volsketchIsCircledMontageGlyphText(label)) {
      badge.style.color = rgba;
      return;
    }
    badge.classList.add("volsketch-gif-overlay--circle-chip");
    badge.style.backgroundColor = rgba;
    badge.style.borderColor = "rgba(27, 31, 36, 0.28)";
    badge.classList.add("volsketch-gif-overlay--on-tint");
  }

  function volsketchResetBadgeChrome(badge) {
    badge.style.backgroundColor = "";
    badge.style.borderColor = "";
    badge.style.color = "#1a1a1a";
    badge.classList.remove("volsketch-gif-overlay--on-tint");
    badge.classList.remove("volsketch-gif-overlay--circle-chip");
  }

  /** Keep pan/zoom across ``Plotly.restyle`` when selection styling updates (2-D scattergl). */
  function snapshotVolsketchAxesView(gd) {
    try {
      var fl = gd._fullLayout;
      if (!fl) return null;
      function copyRange(axis) {
        if (!axis || !axis.range || axis.range.length < 2) return null;
        return [axis.range[0], axis.range[1]];
      }
      return {
        xrange: copyRange(fl.xaxis),
        yrange: copyRange(fl.yaxis),
      };
    } catch (e) {
      return null;
    }
  }

  function restoreVolsketchAxesView(gd, snap) {
    if (!snap || typeof Plotly === "undefined" || !Plotly.relayout) {
      return Promise.resolve();
    }
    var patch = {};
    if (snap.xrange) patch["xaxis.range"] = snap.xrange;
    if (snap.yrange) patch["yaxis.range"] = snap.yrange;
    if (!Object.keys(patch).length) return Promise.resolve();
    return Plotly.relayout(gd, patch).catch(function() {});
  }

  /** In-dashboard labels only; GIF bytes on disk / save API stay unannotated. */
  function mountVolsketchPreviewOverlay(wrap, badge, covBadge, covVarEl, covValEl, img, spec) {
    wrap._volsketchCancelOverlay = null;
    function hideCov() {
      covVarEl.textContent = "";
      covValEl.textContent = "";
      covVarEl.hidden = true;
      covBadge.hidden = true;
      volsketchResetBadgeChrome(covBadge);
    }
    if (!spec || !spec.style) {
      badge.textContent = "";
      badge.hidden = true;
      volsketchResetBadgeChrome(badge);
      hideCov();
      return;
    }
    if (spec.style === "static") {
      badge.textContent = spec.text || "";
      badge.hidden = !badge.textContent;
      volsketchApplyBadgeBackground(badge, spec.badge_background);
      var ctext = spec.covariate_text != null ? String(spec.covariate_text).trim() : "";
      if (ctext) {
        var vLabStatic = (
          spec.covariate_variable_label != null
            ? String(spec.covariate_variable_label).trim()
            : ""
        );
        covVarEl.textContent = vLabStatic;
        covVarEl.hidden = !vLabStatic;
        covValEl.textContent = ctext;
        covBadge.hidden = false;
        volsketchApplyBadgeBackground(
          covBadge,
          spec.covariate_badge_background != null ? spec.covariate_badge_background : ""
        );
      } else {
        hideCov();
      }
      return;
    }
    if (spec.style === "rotate_frames") {
      var frameTexts = spec.frame_covariate_texts || [];
      var frameBgs = spec.frame_badge_backgrounds || [];
      var covVarRotate = (
        spec.covariate_variable_label != null
          ? String(spec.covariate_variable_label).trim()
          : ""
      );
      var fdRot = Math.max(1, parseInt(spec.frame_duration_ms, 10) || 100);
      var totalRot = Math.max(
        1,
        parseInt(spec.total_frames, 10) || frameTexts.length || 1
      );
      badge.textContent = spec.text || "";
      badge.hidden = !badge.textContent;
      volsketchApplyBadgeBackground(badge, spec.badge_background || "");
      if (!frameTexts.length) {
        hideCov();
        return;
      }
      var rafRot = 0;
      var bootRot = 0;
      var startRot = null;
      function applyRotateFrameIdx(fi) {
        fi = Math.max(0, Math.min(frameTexts.length - 1, fi | 0));
        var covt = frameTexts[fi] != null ? String(frameTexts[fi]).trim() : "";
        var bg = frameBgs[fi] != null ? frameBgs[fi] : "";
        covVarEl.textContent = covVarRotate;
        covVarEl.hidden = !covVarRotate;
        if (covt) {
          covValEl.textContent = covt;
          covBadge.hidden = false;
          volsketchApplyBadgeBackground(covBadge, bg);
        } else {
          hideCov();
        }
      }
      function tickRotate(now) {
        if (!wrap.isConnected) {
          if (rafRot) cancelAnimationFrame(rafRot);
          rafRot = 0;
          return;
        }
        if (startRot === null) startRot = now;
        var elapsed = now - startRot;
        var frameIdx = Math.floor(elapsed / fdRot) % totalRot;
        var covIdx = Math.min(
          frameTexts.length - 1,
          Math.floor((frameIdx * frameTexts.length) / totalRot)
        );
        applyRotateFrameIdx(covIdx);
        rafRot = requestAnimationFrame(tickRotate);
      }
      function beginRotateClock() {
        applyRotateFrameIdx(0);
        startRot = null;
        bootRot = requestAnimationFrame(function() {
          bootRot = 0;
          rafRot = requestAnimationFrame(tickRotate);
        });
      }
      wrap._volsketchCancelOverlay = function() {
        if (bootRot) cancelAnimationFrame(bootRot);
        if (rafRot) cancelAnimationFrame(rafRot);
        bootRot = 0;
        rafRot = 0;
      };
      var decRot = img.decode && img.decode();
      if (decRot && typeof decRot.then === "function") {
        decRot.then(beginRotateClock).catch(beginRotateClock);
      } else if (img.complete) {
        beginRotateClock();
      } else {
        img.addEventListener("load", beginRotateClock, { once: true });
      }
      return;
    }
    if (spec.style !== "cycle_segments") {
      badge.hidden = true;
      volsketchResetBadgeChrome(badge);
      hideCov();
      return;
    }
    var labels = spec.segment_labels || [];
    var segBgs = spec.segment_backgrounds || [];
    var segCovBgs = spec.segment_covariate_backgrounds || [];
    var segCov = spec.segment_covariate_texts || [];
    var covVarLabelCycle = (
      spec.covariate_variable_label != null
        ? String(spec.covariate_variable_label).trim()
        : ""
    );
    var fpv = Math.max(1, parseInt(spec.frames_per_segment, 10) || 8);
    var fd = Math.max(1, parseInt(spec.frame_duration_ms, 10) || 100);
    if (!labels.length) {
      badge.hidden = true;
      volsketchResetBadgeChrome(badge);
      hideCov();
      return;
    }
    badge.hidden = false;
    var totalF = Math.max(1, parseInt(spec.total_frames, 10) || labels.length * fpv);
    var rafId = 0;
    var bootRaf = 0;
    var startT = null;
    function applySegmentIdx(idx) {
      idx = Math.max(0, Math.min(labels.length - 1, idx | 0));
      badge.textContent = labels[idx] || "";
      var bg = segBgs[idx] != null ? segBgs[idx] : "";
      volsketchApplyBadgeBackground(badge, bg);
      var covt = "";
      if (segCov.length > idx && segCov[idx] != null) {
        covt = String(segCov[idx]).trim();
      }
      covVarEl.textContent = covVarLabelCycle;
      covVarEl.hidden = !covVarLabelCycle;
      if (covt) {
        covValEl.textContent = covt;
        covBadge.hidden = false;
        var covBg = segCovBgs.length > idx && segCovBgs[idx] != null ? segCovBgs[idx] : bg;
        volsketchApplyBadgeBackground(covBadge, covBg);
      } else {
        covVarEl.hidden = true;
        covVarEl.textContent = "";
        covValEl.textContent = "";
        covBadge.hidden = true;
        volsketchResetBadgeChrome(covBadge);
      }
    }
    function tick(now) {
      if (!wrap.isConnected) {
        if (rafId) cancelAnimationFrame(rafId);
        rafId = 0;
        return;
      }
      if (startT === null) startT = now;
      var elapsed = now - startT;
      var frameIdx = Math.floor(elapsed / fd) % totalF;
      var idx = Math.min(labels.length - 1, Math.floor(frameIdx / fpv));
      applySegmentIdx(idx);
      rafId = requestAnimationFrame(tick);
    }
    function beginClock() {
      applySegmentIdx(0);
      startT = null;
      bootRaf = requestAnimationFrame(function() {
        bootRaf = 0;
        rafId = requestAnimationFrame(tick);
      });
    }
    wrap._volsketchCancelOverlay = function() {
      if (bootRaf) cancelAnimationFrame(bootRaf);
      if (rafId) cancelAnimationFrame(rafId);
      bootRaf = 0;
      rafId = 0;
    };
    var dec = img.decode && img.decode();
    if (dec && typeof dec.then === "function") {
      dec.then(beginClock).catch(beginClock);
    } else if (img.complete) {
      beginClock();
    } else {
      img.addEventListener("load", beginClock, { once: true });
    }
  }

  var GIF_MODE_DESC = {
    cycle: (
      "Single GIF, fixed camera: steps through your selected volumes. "
      + "If you select more than fit, up to <strong>50</strong> are chosen at random for the cycle."
    ),
    rotate_each: (
      "One rotating GIF per selected sketch. "
      + "If you select more than fit, up to <strong>5</strong> volumes are chosen at random."
    ),
    disabled: (
      "No ChimeraX GIF previews — click or lasso volumes on the scatter plot to select them only."
    ),
  };

  function currentGifMode() {
    var el = document.querySelector('input[name="volsketch-gif-mode"]:checked');
    return el ? el.value : "disabled";
  }

  function animationsEnabled() {
    return currentGifMode() !== "disabled";
  }

  function cancelPendingGifWork() {
    if (debounceGifTimer) {
      clearTimeout(debounceGifTimer);
      debounceGifTimer = null;
    }
    gifReqGen++;
    landscapeAnimInFlight = false;
  }

  function clearAnimPreviewState() {
    clearVolsketchPreviewGrid();
    lastAnimToken = null;
    lastBatchMode = null;
    lastChimeraxViewMatrix = "";
    chimeraxViewMatrixUnavailable = false;
    clearAppliedViewMatrix();
    viewMatrixInputDirty = false;
    syncViewMatrixField();
  }

  function showAnimBusySelectionMsg() {
    setAnimateStatus(
      "Wait for ChimeraX to finish the current animation preview before changing the selection.",
      false,
      true
    );
  }

  function finishSelectionWithoutAnimation(reason) {
    cancelPendingGifWork();
    if (!selectedVols.size) {
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
    syncAnimOutputControls();
  }

  function syncAnimOutputControls() {
    syncGifModeDescription();
    syncGifFramesVisibility();
    if (randomSelBtn && META) {
      randomSelBtn.disabled = landscapeAnimInFlight;
      if (landscapeAnimInFlight) {
        randomSelBtn.title = "Wait for ChimeraX to finish the current animation preview.";
      } else {
        randomSelBtn.removeAttribute("title");
      }
    }
    syncViewRotationControls();
  }

  function syncGifModeDescription() {
    if (!gifModeDescEl) return;
    var m = currentGifMode();
    gifModeDescEl.innerHTML = GIF_MODE_DESC[m] || GIF_MODE_DESC.cycle;
  }

  function syncPreviewGridLayout() {
    if (!previewGrid) return;
    var cycle = currentGifMode() === "cycle";
    previewGrid.classList.toggle("volsketch-preview-grid--cycle", cycle);
  }

  function syncGifFramesVisibility() {
    if (!gifFramesWrap) return;
    var m = currentGifMode();
    var showRot = m === "rotate_each";
    gifFramesWrap.style.display = showRot ? "block" : "none";
    syncPreviewGridLayout();
  }
  syncAnimOutputControls();

  var META = null;
  var selectedVols = new Set();
  var lastAnimToken = null;
  var lastBatchMode = null;
  var clickTimer = null;
  var clickLastVol = null;
  var debounceGifTimer = null;
  var gifReqGen = 0;
  var debouncedGifReason = "selection";
  var viewRotations = { x: 0, y: 0, z: 0 };
  var lastChimeraxViewMatrix = "";
  var chimeraxViewMatrixUnavailable = false;
  var appliedViewMatrix = "";
  var viewMatrixInputDirty = false;
  var viewMatrixFieldFocused = false;
  /** True while a generate_animations request is in flight (disable save until rerender completes). */
  var landscapeAnimInFlight = false;

  function syncSaveGifButton() {
    if (!saveGifBtn) return;
    if (!animationsEnabled()) {
      saveGifBtn.disabled = true;
      saveGifBtn.title = "Enable Cycle or Rotate animation output to render GIF previews.";
      syncAnimOutputControls();
      return;
    }
    var hasPreviews = previewGrid && previewGrid.children.length > 0;
    var canSave = !!lastAnimToken && hasPreviews && !landscapeAnimInFlight;
    saveGifBtn.disabled = !canSave;
    if (canSave) {
      saveGifBtn.removeAttribute("title");
    } else if (landscapeAnimInFlight) {
      saveGifBtn.title = "Available after ChimeraX finishes rendering the current previews.";
    } else {
      saveGifBtn.title = "Select volumes and wait for GIF previews to finish rendering.";
    }
    syncAnimOutputControls();
  }

  function syncViewRotationControls() {
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
      viewRotateBtns.forEach(function(btn) {
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
      lastAnimToken
      && previewGrid
      && previewGrid.children.length > 0
      && !landscapeAnimInFlight
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
    viewRotateBtns.forEach(function(btn) {
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

  /** A–Z omitting I, O, U — matches particle explorer / montage_cell_label. */
  var SAFE_LETTERS = (function() {
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

  /** Vol index as in ``vol_NNN.mrc`` (``:03d`` in Python). */
  function volSketchFilenameIndex(v) {
    var n = typeof v === "number" ? v : parseInt(String(v), 10);
    if (isNaN(n)) return String(v);
    return String(n).padStart(3, "0");
  }

  function setPlotStatus(msg, isErr) {
    if (!plotStatusEl) return;
    plotStatusEl.textContent = msg || "";
    plotStatusEl.style.color = isErr ? "var(--error, #b42318)" : "";
  }

  function setAnimateStatus(msg, showProgress, isErr) {
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

  function maxSelectable() {
    return META && META.n_volumes ? META.n_volumes : 10000;
  }

  function syncRandomSelButton() {
    if (!randomSelBtn || !META || META.chimerax_cpus == null) return;
    var n = Math.max(1, Number(META.chimerax_cpus));
    randomSelBtn.textContent = "Choose " + n + " volume" + (n === 1 ? "" : "s") + " at random";
    randomSelBtn.disabled = false;
  }

  function normalizeViewDegrees(degrees) {
    var d = Number(degrees);
    if (!isFinite(d)) return 0;
    d = d % 360;
    if (Math.abs(d) < 1e-9) return 0;
    return d;
  }

  function viewRotationsAreActive() {
    var r = currentViewRotationPayload();
    return Math.abs(r.x) > 1e-9 || Math.abs(r.y) > 1e-9 || Math.abs(r.z) > 1e-9;
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

  function estimatedChimeraxViewMatrixText() {
    var m = mat3Identity();
    ["x", "y", "z"].forEach(function(axis) {
      var deg = normalizeViewDegrees(viewRotations[axis] || 0);
      if (Math.abs(deg) > 1e-9) {
        m = mat3Mul(m, mat3ForAxisTurn(axis, deg));
      }
    });
    var nums = [
      m[0], m[1], m[2], 0,
      m[3], m[4], m[5], 0,
      m[6], m[7], m[8], 0,
    ];
    return "camera " + nums.map(function(n) {
      return Number(n).toFixed(6);
    }).join(",");
  }

  function currentViewRotationPayload() {
    return {
      x: normalizeViewDegrees(viewRotations.x),
      y: normalizeViewDegrees(viewRotations.y),
      z: normalizeViewDegrees(viewRotations.z),
    };
  }

  function viewRotationSummary() {
    var r = currentViewRotationPayload();
    return "requested view turns X " + r.x.toFixed(1)
      + "°, Y " + r.y.toFixed(1)
      + "°, Z " + r.z.toFixed(1) + "°";
  }

  function currentViewMatrixDisplayText() {
    if (lastChimeraxViewMatrix) return lastChimeraxViewMatrix;
    if (viewRotationsAreActive()) return estimatedChimeraxViewMatrixText();
    return "";
  }

  function syncViewMatrixField() {
    if (!viewMatrixInputEl) return;
    if (viewMatrixFieldFocused || viewMatrixInputDirty) return;
    var text = currentViewMatrixDisplayText();
    if (chimeraxViewMatrixUnavailable && !text) {
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

  function validateViewMatrixText(text) {
    var raw = String(text || "").trim();
    if (!raw) {
      return { ok: false, msg: "Enter a ChimeraX view matrix (12 numbers)." };
    }
    var body = raw.toLowerCase().indexOf("camera") === 0 ? raw.slice(6).trim() : raw;
    var nums = body.match(/[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?/g);
    if (!nums || nums.length < 12) {
      return {
        ok: false,
        msg: "View matrix must contain 12 numbers (optionally prefixed with camera).",
      };
    }
    return { ok: true, text: raw };
  }

  function clearAppliedViewMatrix() {
    appliedViewMatrix = "";
  }

  function applyViewMatrixFromField() {
    if (!animationsEnabled()) {
      setAnimateStatus("Enable Cycle or Rotate animation output to render GIF previews.", false, true);
      return;
    }
    if (!viewMatrixInputEl) return;
    var check = validateViewMatrixText(viewMatrixInputEl.value);
    if (!check.ok) {
      setAnimateStatus(check.msg, false, true);
      return;
    }
    appliedViewMatrix = check.text;
    viewMatrixInputDirty = false;
    syncViewMatrixField();
    setAnimateStatus("Applying custom ChimeraX view matrix.", false, false);
    if (selectedVols.size > 0) {
      scheduleAutoGif("view_matrix");
    } else {
      syncSaveGifButton();
    }
  }

  function applyViewRotation(axis) {
    if (!animationsEnabled()) {
      setAnimateStatus("Enable Cycle or Rotate animation output to render GIF previews.", false, true);
      return;
    }
    var deg = viewRotateAngleEl ? Number(viewRotateAngleEl.value) : NaN;
    if (!isFinite(deg)) {
      setAnimateStatus("Enter a finite rotation angle in degrees.", false, true);
      return;
    }
    clearAppliedViewMatrix();
    viewMatrixInputDirty = false;
    viewRotations[axis] = normalizeViewDegrees(Number(viewRotations[axis] || 0) + deg);
    syncViewMatrixField();
    setAnimateStatus("Updated " + viewRotationSummary() + ".", false, false);
    if (selectedVols.size > 0) {
      scheduleAutoGif("view_rotation");
    } else {
      syncSaveGifButton();
    }
  }

  function customDataRowAsArray(row) {
    if (row == null) return [];
    if (Array.isArray(row)) return row;
    if (typeof row === "object" && row.length !== undefined) {
      try {
        return Array.from(row);
      } catch (e) {
        return [];
      }
    }
    return [row];
  }

  function volIdAtPointIndex(trace, i) {
    if (trace.ids && trace.ids.length === trace.x.length && trace.ids[i] != null && trace.ids[i] !== "") {
      var vid = parseInt(String(trace.ids[i]), 10);
      if (!isNaN(vid)) return vid;
    }
    var cd = trace.customdata;
    if (cd && cd.length === trace.x.length && cd[i] != null) {
      var row = customDataRowAsArray(cd[i]);
      if (row.length) return parseInt(String(row[0]), 10);
    }
    return NaN;
  }

  function allPlotVolIds() {
    if (!gd || !gd.data || !gd.data[0]) return [];
    var trace = gd.data[0];
    var n = trace.x ? trace.x.length : 0;
    if (!n) return [];
    var s = new Set();
    for (var i = 0; i < n; i++) {
      var v = volIdAtPointIndex(trace, i);
      if (!isNaN(v)) s.add(v);
    }
    return Array.from(s).sort(function(a, b) { return a - b; });
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

  function updateSelLabel() {
    var arr = Array.from(selectedVols).sort(function(a,b){ return a-b; });
    var nTotal = META && META.n_volumes != null ? Number(META.n_volumes) : NaN;
    if (selSummaryEl) {
      if (!isNaN(nTotal) && nTotal >= 0) {
        selSummaryEl.textContent = arr.length + "/" + nTotal + " sketched volumes selected";
      } else {
        selSummaryEl.textContent = arr.length
          ? (arr.length + " sketched volumes selected")
          : "No sketched volumes selected";
      }
    }
    if (selIndicesDetails && selIndicesList && selIndicesSummaryLabel) {
      selIndicesList.innerHTML = "";
      if (!arr.length) {
        selIndicesDetails.hidden = true;
        selIndicesDetails.open = false;
      } else {
        selIndicesDetails.hidden = false;
        selIndicesSummaryLabel.textContent = (
          "Selected volume indices (" + arr.length + ")"
        );
        var listParts = [];
        arr.forEach(function(v, i) {
          listParts.push(volSketchFilenameIndex(v) + " (" + montageLabelAt(i) + ")");
        });
        var listLabel = listParts.join(", ");
        selIndicesList.setAttribute(
          "aria-label",
          "Selected sketch volume indices: " + listLabel
        );
        arr.forEach(function(v, i) {
          var li = document.createElement("li");
          li.textContent = volSketchFilenameIndex(v) + " (" + montageLabelAt(i) + ")";
          selIndicesList.appendChild(li);
        });
      }
    }
  }

  function volFromPoint(pt) {
    if (!pt) return NaN;
    if (pt.id != null && String(pt.id) !== "") {
      var idv = parseInt(String(pt.id), 10);
      if (!isNaN(idv)) return idv;
    }
    var trace = gd && gd.data && gd.data[0];
    if (trace && pt.pointNumber != null && pt.pointNumber >= 0) {
      var vi = volIdAtPointIndex(trace, pt.pointNumber);
      if (!isNaN(vi)) return vi;
    }
    var row = customDataRowAsArray(pt.customdata);
    return row.length ? parseInt(String(row[0]), 10) : NaN;
  }

  function _hexToRgbVs(hex) {
    var m = /^#([0-9a-f]{3}|[0-9a-f]{6})$/i.exec(String(hex || "").trim());
    if (!m) return null;
    var s = m[1];
    if (s.length === 3) {
      return {
        r: parseInt(s[0] + s[0], 16),
        g: parseInt(s[1] + s[1], 16),
        b: parseInt(s[2] + s[2], 16),
      };
    }
    return {
      r: parseInt(s.slice(0, 2), 16),
      g: parseInt(s.slice(2, 4), 16),
      b: parseInt(s.slice(4, 6), 16),
    };
  }

  function _parseRgbLikeVs(s) {
    var m = /^rgba?\(\s*(\d+)[\s,]+(\d+)[\s,]+(\d+)/i.exec(String(s || "").trim());
    if (!m) return null;
    return { r: +m[1], g: +m[2], b: +m[3] };
  }

  function _fillRelLumVs(css) {
    var rgb = _hexToRgbVs(css) || _parseRgbLikeVs(css);
    if (!rgb) return 0.5;
    function lin(u) {
      u /= 255;
      return u <= 0.03928 ? u / 12.92 : Math.pow((u + 0.055) / 1.055, 2.4);
    }
    var R = lin(rgb.r);
    var G = lin(rgb.g);
    var B = lin(rgb.b);
    return 0.2126 * R + 0.7152 * G + 0.0722 * B;
  }

  function markerFillCssAtIndexVs(trace, gd, i) {
    var full = gd._fullData && gd._fullData[0];
    var expanded = full && full.marker && full.marker.color;
    if (Array.isArray(expanded) && expanded[i] != null && typeof expanded[i] === "string") {
      return expanded[i];
    }
    var mk = trace.marker || {};
    var mc = mk.color;
    if (typeof mc === "string") return mc;
    if (Array.isArray(mc) && mc[i] != null && typeof mc[i] === "string") return mc[i];
    return "#4a5568";
  }

  function selectionOutlineForFillVs(fillCss) {
    return _fillRelLumVs(fillCss) < 0.34 ? "#f1f5f9" : "#0f172a";
  }

  /** Match ``latent3d_landscape_vol_animations.js`` — keep selection overlay sizes tied to API trace. */
  function referenceScatter3dBaseMarkerSize(trace, graphDiv) {
    var full = graphDiv && graphDiv._fullData && graphDiv._fullData[0];
    var marker = (full && full.marker) || (trace && trace.marker);
    if (!marker || marker.size == null) return null;
    var sz = marker.size;
    if (typeof sz === "number" && isFinite(sz) && sz > 0) return sz;
    if (Array.isArray(sz) && sz.length) {
      var minS = Infinity;
      var si;
      for (si = 0; si < sz.length; si++) {
        var one = sz[si];
        if (typeof one === "number" && isFinite(one) && one > 0 && one < minS) minS = one;
      }
      if (minS < Infinity) return minS;
    }
    return null;
  }

  function refreshSelectionHighlight() {
    if (!gd || !gd.data || !gd.data[0]) return;
    var trace = gd.data[0];
    var xs = trace.x;
    if (!xs || !xs.length) return;
    var n = xs.length;
    var hasIds = trace.ids && trace.ids.length === n;
    var cd = trace.customdata;
    if (!hasIds && (!cd || cd.length !== n)) return;
    var mkFallback = 9 * (1 - 0.13);
    var refSz = referenceScatter3dBaseMarkerSize(trace, gd);
    var baseSize = refSz != null ? Math.max(2, refSz) : Math.max(2, mkFallback);
    var hiSize = refSz != null ? Math.max(3, refSz * (11 / 9)) : Math.max(3, mkFallback * (11 / 9));
    var baseOp = 0.75 * 0.8;
    var hiOp = 0.86 * 0.8;
    var dimOp = 0.66 * 0.8;
    var hasSel = selectedVols.size > 0;
    var sizes = new Array(n);
    var opacities = new Array(n);
    var lineWidths = new Array(n);
    var lineColors = new Array(n);
    var sortedSel = Array.from(selectedVols).sort(function(a, b) { return a - b; });
    var labelByVol = {};
    for (var li = 0; li < sortedSel.length; li++) {
      labelByVol[sortedSel[li]] = montageLabelAt(li);
    }
    var texts = new Array(n);
    for (var i = 0; i < n; i++) {
      var v = volIdAtPointIndex(trace, i);
      var on = !isNaN(v) && selectedVols.has(v);
      texts[i] = on ? labelByVol[v] : "";
      if (on) {
        sizes[i] = hiSize;
        opacities[i] = hiOp;
        lineWidths[i] = 1;
        lineColors[i] = selectionOutlineForFillVs(markerFillCssAtIndexVs(trace, gd, i));
      } else {
        sizes[i] = baseSize;
        opacities[i] = hasSel ? dimOp : baseOp;
        lineWidths[i] = 0;
        lineColors[i] = "rgba(0,0,0,0)";
      }
    }
    var tf = {
      size: 12,
      color: "#1a1a1a",
      family: "system-ui, Segoe UI, sans-serif",
    };
    var snap = snapshotVolsketchAxesView(gd);
    var upd = {
      "marker.size": [sizes],
      "marker.opacity": [opacities],
      "marker.line.width": [lineWidths],
      "marker.line.color": [lineColors],
      mode: "markers+text",
      text: [texts],
      textposition: "top center",
      textfont: tf,
    };
    var fallback = {
      "marker.size": [sizes],
      "marker.opacity": [opacities],
      "marker.line.width": [lineWidths],
      "marker.line.color": [lineColors],
      mode: "markers+text",
      text: [texts],
      textposition: "top center",
      textfont: tf,
    };
    var p = Plotly.restyle(gd, upd, [0]);
    function finishAxisRestore() {
      return restoreVolsketchAxesView(gd, snap);
    }
    if (p && typeof p.then === "function") {
      p.catch(function() {
        return Plotly.restyle(gd, fallback, [0]);
      }).then(function() {
        requestAnimationFrame(finishAxisRestore);
      });
    } else {
      requestAnimationFrame(finishAxisRestore);
    }
  }

  function toggleVol(v) {
    if (landscapeAnimInFlight) {
      showAnimBusySelectionMsg();
      return;
    }
    if (selectedVols.has(v)) {
      selectedVols.delete(v);
    } else {
      if (selectedVols.size >= maxSelectable()) {
        setAnimateStatus("Selection capped at the number of sketch volumes.", false, true);
        return;
      }
      selectedVols.add(v);
    }
    setPlotStatus("", false);
    setAnimateStatus("", false, false);
    updateSelLabel();
    refreshSelectionHighlight();
    scheduleAutoGif("selection");
  }

  function scheduleAutoGif(reason) {
    reason = reason || "selection";
    debouncedGifReason = reason;
    if (debounceGifTimer) clearTimeout(debounceGifTimer);
    if (!animationsEnabled()) {
      debounceGifTimer = setTimeout(function() {
        debounceGifTimer = null;
        finishSelectionWithoutAnimation(reason);
      }, 600);
      return;
    }
    debounceGifTimer = setTimeout(function() {
      debounceGifTimer = null;
      var r = debouncedGifReason;
      debouncedGifReason = "selection";
      runAutoGif(r);
    }, 600);
  }

  function runAutoGif(reason) {
    reason = reason || "selection";
    if (!animationsEnabled()) {
      finishSelectionWithoutAnimation(reason);
      return;
    }
    var vols = Array.from(selectedVols).sort(function(a,b){ return a-b; });
    if (!vols.length) {
      landscapeAnimInFlight = false;
      clearAnimPreviewState();
      setAnimateStatus("", false, false);
      syncSaveGifButton();
      return;
    }
    var myGen = ++gifReqGen;
    landscapeAnimInFlight = true;
    lastChimeraxViewMatrix = "";
    chimeraxViewMatrixUnavailable = false;
    syncSaveGifButton();
    syncViewMatrixField();
    var cpus = META && META.chimerax_cpus != null ? Number(META.chimerax_cpus) : 1;
    var cpuPhrase = cpus === 1 ? "1 CPU" : cpus + " CPUs";
    var keepPreviews = previewGrid && previewGrid.children.length > 0 && (
      reason === "color" || reason === "gif_mode" || reason === "gif_frames"
      || reason === "view_rotation" || reason === "view_matrix"
    );
    var statusMsg;
    if (keepPreviews && reason === "color") {
      statusMsg = (
        "Re-rendering volume previews with the new plot colors — current GIFs keep playing until "
        + "the update finishes (ChimeraX, "
        + cpuPhrase + ")…"
      );
    } else if (keepPreviews && (reason === "gif_mode" || reason === "gif_frames")) {
      statusMsg = (
        "Re-rendering volume previews for the new animation settings — current GIFs keep playing until "
        + "the update finishes (ChimeraX, "
        + cpuPhrase + ")…"
      );
    } else if (keepPreviews && (reason === "view_rotation" || reason === "view_matrix")) {
      statusMsg = (
        "Re-rendering volume previews for the new view angle — current GIFs keep playing until "
        + "the update finishes (ChimeraX, "
        + cpuPhrase + ")…"
      );
    } else {
      statusMsg = "Rendering GIFs with ChimeraX (using " + cpuPhrase + ")…";
    }
    setAnimateStatus(statusMsg, true, false);
    if (!keepPreviews) {
      teardownVolsketchPreviews();
      previewGrid.innerHTML = "";
    }
    var animPayload = {
      vol_indices: vols,
      mode: currentGifMode(),
      gif_frames: parseInt(gifFrames.value, 10) || 20,
      chimerax_cpus: cpus,
      color_mode: colorSel.value,
    };
    if (appliedViewMatrix) {
      animPayload.view_matrix = appliedViewMatrix;
    } else {
      animPayload.view_rotations = currentViewRotationPayload();
    }
    if (
      reason !== "color"
      && currentGifMode() === "cycle"
      && lastAnimToken
      && lastBatchMode === "rotate_each"
    ) {
      animPayload.reuse_rotate_keyframes_token = lastAnimToken;
    }
    postJson("/api/landscape_volpca/generate_animations", animPayload).then(function(j) {
      if (myGen !== gifReqGen) return;
      landscapeAnimInFlight = false;
      lastAnimToken = j.token;
      lastBatchMode = j.batch_mode != null ? j.batch_mode : null;
      lastChimeraxViewMatrix = j.view_matrix != null ? String(j.view_matrix).trim() : "";
      if (!lastChimeraxViewMatrix && j.items && j.items.length && j.items[0].view_matrix != null) {
        lastChimeraxViewMatrix = String(j.items[0].view_matrix).trim();
      }
      chimeraxViewMatrixUnavailable = !lastChimeraxViewMatrix;
      viewMatrixInputDirty = false;
      syncViewMatrixField();
      var ds = j.duration_s != null ? Number(j.duration_s) : NaN;
      var doneMsg = !isNaN(ds)
        ? ("Finished animating in " + Math.round(ds) + " s.")
        : "Finished animating.";
      setAnimateStatus(doneMsg, false, false);
      if (!j.items || !j.items.length) {
        syncSaveGifButton();
        return;
      }
      teardownVolsketchPreviews();
      previewGrid.innerHTML = "";
      j.items.forEach(function(it) {
        var fig = document.createElement("figure");
        var wrap = document.createElement("div");
        wrap.className = "volsketch-anim-wrap";
        var badge = document.createElement("span");
        badge.className = "volsketch-gif-overlay";
        var covBadge = document.createElement("span");
        covBadge.className = "volsketch-gif-overlay volsketch-gif-overlay--cov";
        var covVarEl = document.createElement("span");
        covVarEl.className = "volsketch-gif-cov-var";
        covVarEl.hidden = true;
        var covValEl = document.createElement("span");
        covValEl.className = "volsketch-gif-cov-val";
        covBadge.appendChild(covVarEl);
        covBadge.appendChild(covValEl);
        var img = document.createElement("img");
        img.alt = "";
        img.src = "data:image/gif;base64," + it.gif_b64;
        wrap.appendChild(img);
        wrap.appendChild(badge);
        wrap.appendChild(covBadge);
        fig.appendChild(wrap);
        mountVolsketchPreviewOverlay(wrap, badge, covBadge, covVarEl, covValEl, img, it.preview_overlay);
        previewGrid.appendChild(fig);
      });
      syncPreviewGridLayout();
      if (j.rendered_vol_indices && j.rendered_vol_indices.length) {
        var synced = new Set();
        j.rendered_vol_indices.forEach(function(x) {
          var v = parseInt(x, 10);
          if (!isNaN(v)) synced.add(v);
        });
        selectedVols = synced;
        updateSelLabel();
        refreshSelectionHighlight();
      }
      syncSaveGifButton();
    }).catch(function(e) {
      if (myGen !== gifReqGen) return;
      landscapeAnimInFlight = false;
      chimeraxViewMatrixUnavailable = false;
      syncViewMatrixField();
      setAnimateStatus(String(e), false, true);
      syncSaveGifButton();
    });
  }

  function allVolsketchAxisValues() {
    var out = [];
    var j;
    var nPc = META && META.n_pc != null ? META.n_pc : 0;
    var nUm = META && META.n_umap != null ? META.n_umap : 0;
    for (j = 0; j < nPc; j++) out.push("pc:" + j);
    for (j = 0; j < nUm; j++) out.push("umap:" + j);
    return out;
  }

  function bumpVolsketchYIfSameAxis() {
    var vals = allVolsketchAxisValues();
    if (vals.length < 2 || pcx.value !== pcy.value) return;
    var ix = vals.indexOf(pcx.value);
    var next = (ix + 1) % vals.length;
    pcy.value = vals[next];
  }

  function fillAxisSelects() {
    if (!META) return;
    var nPc = META.n_pc || 0;
    var nUm = META.n_umap || 0;
    var evr = META.explained_variance_ratio || null;
    var prevX = pcx.value;
    var prevY = pcy.value;
    pcx.innerHTML = "";
    pcy.innerHTML = "";
    function appendOptions(sel) {
      var j;
      for (j = 0; j < nPc; j++) {
        var ox = document.createElement("option");
        ox.value = "pc:" + j;
        var lab = "Vol PC" + (j + 1);
        if (evr && j < evr.length) {
          lab += " (" + (100 * Number(evr[j])).toFixed(1) + "%)";
        }
        ox.textContent = lab;
        sel.appendChild(ox);
      }
      for (j = 0; j < nUm; j++) {
        var ou = document.createElement("option");
        ou.value = "umap:" + j;
        ou.textContent = "UMAP" + (j + 1);
        sel.appendChild(ou);
      }
    }
    appendOptions(pcx);
    appendOptions(pcy);
    var vals = allVolsketchAxisValues();
    if (!vals.length) return;
    pcx.value = vals.indexOf(prevX) >= 0 ? prevX : vals[0];
    pcy.value = vals.indexOf(prevY) >= 0 ? prevY : vals[Math.min(1, vals.length - 1)];
    bumpVolsketchYIfSameAxis();
  }

  function fillColorSelect(meta) {
    var opts = meta.color_options || [];
    var prev = colorSel.value;
    colorSel.innerHTML = "";
    opts.forEach(function(o) {
      var el = document.createElement("option");
      el.value = o.value;
      el.textContent = o.label;
      if (o.value === "state" && !meta.has_state_color) el.disabled = true;
      colorSel.appendChild(el);
    });
    var want = prev;
    if (want === "state" && !meta.has_state_color) want = "none";
    var ok = false;
    for (var i = 0; i < colorSel.options.length; i++) {
      var opt = colorSel.options[i];
      if (opt.value === want && !opt.disabled) { ok = true; break; }
    }
    colorSel.value = ok ? want : "none";
  }

  function loadMeta(cb) {
    fetch("/api/landscape_volpca/meta")
      .then(function(r) { return r.json(); })
      .then(function(j) {
        if (!j.ok) {
          setPlotStatus(j.error || "Meta failed", true);
          if (randomSelBtn) {
            randomSelBtn.disabled = true;
            randomSelBtn.textContent = "Choose volumes at random";
          }
          return;
        }
        META = j;
        fillAxisSelects();
        fillColorSelect(j);
        syncRandomSelButton();
        updateSelLabel();
        if (cb) cb();
      })
      .catch(function(e) {
        setPlotStatus(String(e), true);
        if (randomSelBtn) {
          randomSelBtn.disabled = true;
          randomSelBtn.textContent = "Choose volumes at random";
        }
      });
  }

  function bindPlotEvents() {
    gd.on("plotly_selected", function(ev) {
      if (landscapeAnimInFlight) {
        showAnimBusySelectionMsg();
        return;
      }
      if (!ev || !ev.points || !ev.points.length) return;
      selectedVols.clear();
      ev.points.forEach(function(pt) {
        var v = volFromPoint(pt);
        if (!isNaN(v)) selectedVols.add(v);
      });
      if (selectedVols.size > maxSelectable()) {
        var arr = Array.from(selectedVols).sort(function(a,b){ return a-b; });
        selectedVols = new Set(arr.slice(0, maxSelectable()));
        setAnimateStatus("Selection capped at the number of sketch volumes.", false, true);
      } else {
        setAnimateStatus("", false, false);
      }
      updateSelLabel();
      refreshSelectionHighlight();
      scheduleAutoGif("selection");
    });
    gd.on("plotly_click", function(ev) {
      if (!ev || !ev.points || !ev.points.length) return;
      var v = volFromPoint(ev.points[0]);
      if (isNaN(v)) return;
      if (clickTimer !== null && clickLastVol === v) {
        clearTimeout(clickTimer);
        clickTimer = null;
        clickLastVol = null;
        if (landscapeAnimInFlight) {
          showAnimBusySelectionMsg();
          return;
        }
        toggleVol(v);
        return;
      }
      clickLastVol = v;
      if (clickTimer) clearTimeout(clickTimer);
      clickTimer = setTimeout(function() {
        clickTimer = null;
        clickLastVol = null;
      }, 320);
    });
  }

  function reloadScatter(scatterOpts) {
    scatterOpts = scatterOpts || {};
    if (!META) return;
    var x = pcx.value;
    var y = pcy.value;
    var c = colorSel.value;
    var url = "/api/landscape_volpca/scatter?axis_x=" + encodeURIComponent(x)
      + "&axis_y=" + encodeURIComponent(y) + "&color=" + encodeURIComponent(c);
    setPlotStatus("Loading scatter…", false);
    fetch(url)
      .then(function(r) {
        if (!r.ok) return r.json().then(function(j) { throw new Error(j.error || r.status); });
        return r.json();
      })
      .then(function(fig) {
        var plotOpts = { responsive: true, displayModeBar: true, doubleClick: false };
        var afterColor = !!scatterOpts.afterColorChange;
        var hasPlot = !!(gd && gd.data && gd.data.length);
        var prevSelections = null;
        if (afterColor && gd && gd.layout && gd.layout.selections && gd.layout.selections.length) {
          try {
            prevSelections = JSON.parse(JSON.stringify(gd.layout.selections));
          } catch (err) {
            prevSelections = null;
          }
        }
        var plotPromise;
        function redrawFull() {
          if (afterColor && hasPlot && typeof Plotly.react === "function") {
            return Plotly.react(gd, fig.data, fig.layout, plotOpts);
          }
          Plotly.purge(gd);
          return Plotly.newPlot(gd, fig.data, fig.layout, plotOpts);
        }

        var usedFullRedraw = false;
        if (afterColor && hasPlot && fig.data && fig.data[0] && gd.data && gd.data[0]) {
          var g0 = gd.data[0];
          var d0 = fig.data[0];
          var nOld = g0.x ? g0.x.length : 0;
          var nNew = d0.x ? d0.x.length : 0;
          if (nOld === nNew && typeof Plotly.restyle === "function") {
            plotPromise = Plotly.restyle(
              gd,
              {
                marker: d0.marker,
                hovertemplate: d0.hovertemplate,
                customdata: d0.customdata,
              },
              [0]
            ).catch(function() {
              usedFullRedraw = true;
              return redrawFull();
            });
          } else {
            usedFullRedraw = true;
            plotPromise = redrawFull();
          }
        } else {
          usedFullRedraw = true;
          plotPromise = redrawFull();
        }
        return plotPromise
          .then(function() {
            if (usedFullRedraw && prevSelections && prevSelections.length) {
              return Plotly.relayout(gd, { selections: prevSelections }).catch(
                function() { /* older Plotly or incompatible selection payload */ }
              );
            }
          })
          .then(function() {
            setPlotStatus("", false);
            bindPlotEvents();
            refreshSelectionHighlight();
            Plotly.Plots.resize(gd);
            if (scatterOpts.afterColorChange && selectedVols.size > 0) {
              scheduleAutoGif("color");
            }
          });
      })
      .catch(function(e) { setPlotStatus(String(e), true); });
  }

  pcx.addEventListener("change", function() {
    if (!META) return;
    bumpVolsketchYIfSameAxis();
    reloadScatter();
  });
  pcy.addEventListener("change", function() {
    if (!META) return;
    if (pcx.value === pcy.value) {
      var vals = allVolsketchAxisValues();
      if (vals.length >= 2) {
        var iy = vals.indexOf(pcy.value);
        pcx.value = vals[(iy + 1) % vals.length];
      }
    }
    reloadScatter();
  });
  colorSel.addEventListener("change", function() {
    reloadScatter({ afterColorChange: true });
  });

  document.getElementById("volsketch-clear-sel").addEventListener("click", function() {
    selectedVols.clear();
    updateSelLabel();
    refreshSelectionHighlight();
    setAnimateStatus("", false, false);
    scheduleAutoGif("selection");
  });

  function postJson(url, body) {
    return fetch(url, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body || {}),
    }).then(function(r) {
      return r.json().then(function(j) {
        if (!r.ok) throw new Error(j.error || r.status);
        return j;
      });
    });
  }

  document.querySelectorAll('input[name="volsketch-gif-mode"]').forEach(function(radio) {
    radio.addEventListener("change", function() {
      syncAnimOutputControls();
      scheduleAutoGif("gif_mode");
    });
  });
  gifFrames.addEventListener("change", function() {
    scheduleAutoGif("gif_frames");
  });
  document.querySelectorAll("[data-volsketch-view-axis]").forEach(function(btn) {
    btn.addEventListener("click", function() {
      var axis = String(btn.getAttribute("data-volsketch-view-axis") || "").toLowerCase();
      if (axis === "x" || axis === "y" || axis === "z") {
        applyViewRotation(axis);
      }
    });
  });
  if (viewMatrixInputEl) {
    viewMatrixInputEl.addEventListener("focus", function() {
      viewMatrixFieldFocused = true;
    });
    viewMatrixInputEl.addEventListener("blur", function() {
      viewMatrixFieldFocused = false;
      syncViewMatrixField();
    });
    viewMatrixInputEl.addEventListener("input", function() {
      viewMatrixInputDirty = true;
    });
  }
  if (viewMatrixApplyBtn) {
    viewMatrixApplyBtn.addEventListener("click", applyViewMatrixFromField);
  }
  syncViewMatrixField();

  if (randomSelBtn) {
    randomSelBtn.addEventListener("click", function() {
      if (landscapeAnimInFlight) {
        showAnimBusySelectionMsg();
        return;
      }
      var ids = allPlotVolIds();
      if (!ids.length) {
        setAnimateStatus("Wait for the scatter plot to finish loading.", false, true);
        return;
      }
      if (debounceGifTimer) {
        clearTimeout(debounceGifTimer);
        debounceGifTimer = null;
      }
      gifReqGen++;
      landscapeAnimInFlight = false;
      clearVolsketchPreviewGrid();
      lastAnimToken = null;
      lastBatchMode = null;
      setAnimateStatus("", false, false);
      syncSaveGifButton();
      function applyRandomPickAndHighlight() {
        var k = META && META.chimerax_cpus != null ? Number(META.chimerax_cpus) : 1;
        k = Math.max(1, Math.min(k, ids.length, maxSelectable()));
        selectedVols.clear();
        shufflePick(ids, k).forEach(function(v) { selectedVols.add(v); });
        updateSelLabel();
        refreshSelectionHighlight();
        scheduleAutoGif("selection");
      }
      if (gd && typeof Plotly !== "undefined" && Plotly.relayout) {
        var rp = Plotly.relayout(gd, { selections: [] });
        var p = rp && typeof rp.then === "function" ? rp.catch(function() {}) : Promise.resolve();
        p.then(applyRandomPickAndHighlight);
      } else {
        applyRandomPickAndHighlight();
      }
    });
  }

  document.getElementById("volsketch-save-gif").addEventListener("click", function() {
    if (!lastAnimToken) {
      setAnimateStatus("Select volumes and wait for previews to finish rendering.", false, true);
      return;
    }
    var custom = window.prompt(
      "Save folder (leave empty for default: this run’s landscape k-means directory):",
      ""
    );
    if (custom === null) return;
    setAnimateStatus("Saving…", true, false);
    var body = { token: lastAnimToken };
    if (custom && custom.trim()) body.out_dir = custom.trim();
    postJson("/api/landscape_volpca/save_animations", body)
      .then(function(j) {
        setAnimateStatus("Saved: " + (j.paths || []).join(", "), false, false);
      })
      .catch(function(e) { setAnimateStatus(String(e), false, true); });
  });

  window.addEventListener("resize", function() {
    if (gd && gd.layout) Plotly.Plots.resize(gd);
  });

  loadMeta(reloadScatter);
  updateSelLabel();
  syncViewMatrixField();
  syncSaveGifButton();
})();
