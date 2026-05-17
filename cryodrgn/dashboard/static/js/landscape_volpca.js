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
  var viewRotateSummaryEl = document.getElementById("volsketch-view-rotate-summary");
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

  /** Volume color → tinted badge background; letters stay black. */
  function volsketchApplyBadgeBackground(badge, color) {
    badge.style.color = "#1a1a1a";
    badge.classList.remove("volsketch-gif-overlay--on-tint");
    if (typeof color !== "string" || !color.trim()) {
      badge.style.backgroundColor = "";
      badge.style.borderColor = "";
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
    if (rgba) {
      badge.style.backgroundColor = rgba;
      badge.style.borderColor = "rgba(27, 31, 36, 0.28)";
      badge.classList.add("volsketch-gif-overlay--on-tint");
      return;
    }
    badge.style.backgroundColor = "";
    badge.style.borderColor = "";
  }

  function volsketchResetBadgeChrome(badge) {
    badge.style.backgroundColor = "";
    badge.style.borderColor = "";
    badge.style.color = "#1a1a1a";
    badge.classList.remove("volsketch-gif-overlay--on-tint");
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
    if (spec.style !== "cycle_segments") {
      badge.hidden = true;
      volsketchResetBadgeChrome(badge);
      hideCov();
      return;
    }
    var labels = spec.segment_labels || [];
    var segBgs = spec.segment_backgrounds || [];
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
        volsketchApplyBadgeBackground(covBadge, bg);
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
      var loopMs = totalF * fd;
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
  };

  function currentGifMode() {
    var el = document.querySelector('input[name="volsketch-gif-mode"]:checked');
    return el ? el.value : "cycle";
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
  syncGifModeDescription();
  syncGifFramesVisibility();

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
  /** True while a generate_animations request is in flight (disable save until rerender completes). */
  var landscapeAnimInFlight = false;

  function syncSaveGifButton() {
    if (!saveGifBtn) return;
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
    syncViewRotationControls();
  }

  function syncViewRotationControls() {
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

  function syncViewRotationSummary() {
    if (!viewRotateSummaryEl) return;
    if (lastChimeraxViewMatrix) {
      viewRotateSummaryEl.textContent = "ChimeraX view matrix:\n" + lastChimeraxViewMatrix;
    } else if (chimeraxViewMatrixUnavailable) {
      viewRotateSummaryEl.textContent = "ChimeraX view matrix: not reported by ChimeraX for this render.";
    } else if (landscapeAnimInFlight) {
      viewRotateSummaryEl.textContent = "ChimeraX view matrix: rendering…";
    } else {
      viewRotateSummaryEl.textContent = "ChimeraX view matrix: available after animation loads.";
    }
  }

  function applyViewRotation(axis) {
    var deg = viewRotateAngleEl ? Number(viewRotateAngleEl.value) : NaN;
    if (!isFinite(deg)) {
      setAnimateStatus("Enter a finite rotation angle in degrees.", false, true);
      return;
    }
    viewRotations[axis] = normalizeViewDegrees(Number(viewRotations[axis] || 0) + deg);
    syncViewRotationSummary();
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

  function refreshSelectionHighlight() {
    if (!gd || !gd.data || !gd.data[0]) return;
    var trace = gd.data[0];
    var xs = trace.x;
    if (!xs || !xs.length) return;
    var n = xs.length;
    var hasIds = trace.ids && trace.ids.length === n;
    var cd = trace.customdata;
    if (!hasIds && (!cd || cd.length !== n)) return;
    var mk = 9 * (1 - 0.13);
    var baseSize = Math.max(2, mk);
    var hiSize = Math.max(3, mk * (11 / 9));
    var baseOp = 0.75;
    var hiOp = 0.86;
    var dimOp = 0.66;
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
    debouncedGifReason = reason || "selection";
    if (debounceGifTimer) clearTimeout(debounceGifTimer);
    debounceGifTimer = setTimeout(function() {
      debounceGifTimer = null;
      var r = debouncedGifReason;
      debouncedGifReason = "selection";
      runAutoGif(r);
    }, 600);
  }

  function runAutoGif(reason) {
    reason = reason || "selection";
    var vols = Array.from(selectedVols).sort(function(a,b){ return a-b; });
    if (!vols.length) {
      landscapeAnimInFlight = false;
      clearVolsketchPreviewGrid();
      lastAnimToken = null;
      lastBatchMode = null;
      lastChimeraxViewMatrix = "";
      chimeraxViewMatrixUnavailable = false;
      syncViewRotationSummary();
      setAnimateStatus("", false, false);
      syncSaveGifButton();
      return;
    }
    var myGen = ++gifReqGen;
    landscapeAnimInFlight = true;
    lastChimeraxViewMatrix = "";
    chimeraxViewMatrixUnavailable = false;
    syncSaveGifButton();
    syncViewRotationSummary();
    var cpus = META && META.chimerax_cpus != null ? Number(META.chimerax_cpus) : 1;
    var cpuPhrase = cpus === 1 ? "1 CPU" : cpus + " CPUs";
    var keepPreviews = previewGrid && previewGrid.children.length > 0 && (
      reason === "color" || reason === "gif_mode" || reason === "gif_frames"
      || reason === "view_rotation"
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
    } else if (keepPreviews && reason === "view_rotation") {
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
      view_rotations: currentViewRotationPayload(),
    };
    if (
      currentGifMode() === "cycle"
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
      syncViewRotationSummary();
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
      syncViewRotationSummary();
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
      syncGifModeDescription();
      syncGifFramesVisibility();
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

  if (randomSelBtn) {
    randomSelBtn.addEventListener("click", function() {
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
  syncViewRotationSummary();
  syncSaveGifButton();
})();
