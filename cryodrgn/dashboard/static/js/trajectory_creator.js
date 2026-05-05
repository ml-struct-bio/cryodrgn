(function() {
  var b = window.CRYO_TRAJECTORY_BOOT;
  var trajAxisCols = b.trajAxisCols;
  var colorCols = b.colorCols;
  var covariateDisplayMap = b.covariateDisplayMap;
  var dx = b.defaultX;
  var dy = b.defaultY;
  var zdim = b.zdim;
  var expWorkdir = b.expWorkdir;
  var trajChimeraxCpus = b.trajChimeraxCpus;

  var trajectoryMode = "direct";
  var anchorIndicesActive = null;
  var anchorTrajXY = null;
  var editableTrajXY = null;
  var startXY = null;
  var endXY = null;
  var scatterLoadGeneration = 0;
  var volumeGeneration = 0;
  var dragEndpoint = null;
  var dragPointIndex = -1;
  var dragListenersWired = false;

  var gd = document.getElementById("scatter");
  var overlay = document.getElementById("scatter-rendering-overlay");
  var paletteFieldset = document.getElementById("traj-palette-radios");
  var plotStatus = document.getElementById("scatter-plot-status");
  var trajStatus = document.getElementById("traj-status");
  var trajProgress = document.getElementById("traj-progress");
  var trajVolColumn = document.getElementById("traj-vol-column");
  var btnModeDirectEl = document.getElementById("btn-mode-direct");
  var btnModeNearestEl = document.getElementById("btn-mode-nearest");
  var btnGenerateVolumesEl = document.getElementById("btn-generate-volumes");
  var btnResetEl = document.getElementById("btn-reset");
  var nPointsEl = document.getElementById("n-points");
  var trajPointsRowEl = document.getElementById("traj-points-row");
  var trajNeighborRowEl = document.getElementById("traj-neighbor-row");
  var maxNeighborsEl = document.getElementById("max-neighbors");
  var avgNeighborsEl = document.getElementById("avg-neighbors");
  var trajPointsLabelEl = document.getElementById("traj-points-label");
  var trajSettingsHintEl = document.getElementById("traj-settings-hint");
  var trajZPanel = document.getElementById("traj-z-panel");
  var trajZNums = document.getElementById("traj-z-nums");
  var trajZPre = document.getElementById("traj-z-pre");
  var btnSaveZpathWorkdir = document.getElementById("btn-save-zpath-workdir");
  var btnSaveZpathAs = document.getElementById("btn-save-zpath-as");
  var trajAnchorPick = document.getElementById("traj-anchor-pick");
  var trajAnchorManualRow = document.getElementById("traj-anchor-manual-row");
  var btnAnchorImportPkl = document.getElementById("btn-anchor-import-pkl");
  var btnAnchorKmeans = document.getElementById("btn-anchor-kmeans");
  var btnAnchorRandom = document.getElementById("btn-anchor-random");
  var anchorIndicesManualEl = document.getElementById("anchor-indices-manual");
  var btnAnchorManualLoad = document.getElementById("btn-anchor-manual-load");
  var btnAnchorRemove = document.getElementById("btn-anchor-remove");
  var trajFileBrowserPanel = document.getElementById("traj-file-browser-panel");
  var serverFileBrowser = document.getElementById("server-file-browser");
  var fbList = document.getElementById("fb-list");
  var fbPath = document.getElementById("fb-path");
  var fbUp = document.getElementById("fb-up");
  var fbCancel = document.getElementById("fb-cancel");
  var fbSaveRow = document.getElementById("fb-save-row");
  var fbSaveName = document.getElementById("fb-save-name");
  var fbSaveBtn = document.getElementById("fb-save-btn");
  var fbCurrentDir = null;
  var fbMode = "import";
  var lastZPathTxt = null;
  var coordsRequestId = 0;
  var coordsDragTimer = 0;
  var lastTrajectoryVolumeCacheId = "";
  var generatedTrajectoryVolumeCount = 0;

  var scatterControlsCard = document.getElementById("scatter-controls-card");
  var scatterControlsSide = scatterControlsCard ? scatterControlsCard.closest(".cryo-dash-side") : null;
  var scatterPlotStack = document.getElementById("scatter-plot-stack");
  var plotlyStackEl = gd && gd.parentElement;

  var TRAJ_PALETTE_KEYS = { Viridis: 1, Plasma: 1, Inferno: 1, Magma: 1, Cividis: 1, Turbo: 1 };
  function selectedTrajPalette() {
    var r = document.querySelector("input[name=\"traj_palette\"]:checked");
    var v = r && r.value ? String(r.value) : "Viridis";
    return TRAJ_PALETTE_KEYS[v] ? v : "Viridis";
  }

  function fillSelect(sel, values, includeNone) {
    sel.innerHTML = "";
    if (includeNone) {
      var o = document.createElement("option");
      o.value = "none"; o.textContent = "None";
      sel.appendChild(o);
    }
    values.forEach(function(c) {
      var o = document.createElement("option");
      o.value = c;
      o.textContent = covariateDisplayMap[c] || c;
      sel.appendChild(o);
    });
  }

  var sx = document.getElementById("sx");
  var sy = document.getElementById("sy");
  var sc = document.getElementById("sc");
  fillSelect(sx, trajAxisCols, false);
  fillSelect(sy, trajAxisCols, false);
  fillSelect(sc, colorCols, true);
  sx.value = trajAxisCols.indexOf(dx) >= 0 ? dx : trajAxisCols[0];
  sy.value = trajAxisCols.indexOf(dy) >= 0 ? dy : trajAxisCols[Math.min(1, trajAxisCols.length - 1)];

  /* Matches scatter_json / particle explorer: only ``labels`` is discrete (no sequential palette). */
  function trajColorByIsDiscrete() {
    return sc.value === "labels";
  }

  function syncTrajPaletteFieldset() {
    if (!paletteFieldset) return;
    var rendering = !!(overlay && overlay.classList.contains("cryo-plot-rendering-overlay--show"));
    var hasColor = sc.value && sc.value !== "none";
    var suppress = rendering || !hasColor || trajColorByIsDiscrete();
    paletteFieldset.classList.toggle("traj-palette-radios--suppressed", suppress);
  }

  function setRendering(on) {
    if (overlay) {
      overlay.classList.toggle("cryo-plot-rendering-overlay--show", on);
      overlay.setAttribute("aria-hidden", on ? "false" : "true");
    }
    syncTrajPaletteFieldset();
  }

  function setTrajStatus(msg, showProgress) {
    if (trajStatus) trajStatus.textContent = msg || "";
    if (trajProgress) trajProgress.hidden = !showProgress;
  }

  function currentNPoints() {
    var n = parseInt(nPointsEl.value, 10);
    if (!Number.isFinite(n)) n = 4;
    var minN = (hasAnchorIndices() && currentAnchorTraversalMode() === "direct") ? 0 : 2;
    return Math.max(minN, Math.min(20, n));
  }

  function syncPointsSelectorRange() {
    if (!nPointsEl) return;
    var minN = (hasAnchorIndices() && currentAnchorTraversalMode() === "direct") ? 0 : 2;
    var current = parseInt(nPointsEl.value, 10);
    if (!Number.isFinite(current)) current = 4;
    var next = Math.max(minN, Math.min(20, current));
    if (String(nPointsEl.options[0] && nPointsEl.options[0].value) === String(minN)) {
      nPointsEl.value = String(next);
      return;
    }
    nPointsEl.innerHTML = "";
    for (var i = minN; i <= 20; i++) {
      var o = document.createElement("option");
      o.value = String(i);
      o.textContent = String(i);
      nPointsEl.appendChild(o);
    }
    nPointsEl.value = String(next);
  }

  function currentMaxNeighbors() {
    var n = parseInt(maxNeighborsEl && maxNeighborsEl.value, 10);
    if (!Number.isFinite(n)) n = 10;
    return Math.max(2, Math.min(20, n));
  }

  function currentAvgNeighbors() {
    var n = parseInt(avgNeighborsEl && avgNeighborsEl.value, 10);
    if (!Number.isFinite(n)) n = 5;
    return Math.max(2, Math.min(20, n));
  }

  function requestedAnchorNPoints(mode) {
    return mode === "direct" ? 0 : currentNPoints();
  }

  function hasAnchorIndices() {
    return !!(anchorIndicesActive && anchorIndicesActive.length >= 2);
  }

  function currentAnchorTraversalMode() {
    return trajectoryMode === "graph" ? "graph" : "direct";
  }

  function syncTraversalSettingsText() {
    var hasAnchors = hasAnchorIndices();
    if (btnModeNearestEl) {
      btnModeNearestEl.textContent = hasAnchors ? "Graph traversal" : "Nearest neighbor";
    }
    if (!trajPointsLabelEl || !trajSettingsHintEl) return;
    var graphMode = hasAnchors && currentAnchorTraversalMode() === "graph";
    syncPointsSelectorRange();
    if (trajPointsRowEl) trajPointsRowEl.style.display = graphMode ? "none" : "flex";
    if (trajNeighborRowEl) {
      trajNeighborRowEl.classList.toggle("cryo-traj-neighbor-row--show", graphMode);
    }
    if (!hasAnchors) {
      trajPointsLabelEl.textContent = "Points";
      trajSettingsHintEl.textContent = "Drag the orange and navy handles on the plot to edit the trajectory.";
      return;
    }
    if (graphMode) {
      trajSettingsHintEl.textContent =
        "Anchors are fixed. Max and Avg. neighbors control graph traversal between each anchor pair.";
      return;
    }
    trajPointsLabelEl.textContent = "Interpolation points between anchors";
    trajSettingsHintEl.textContent =
      "Anchors are fixed. Points sets interpolation points between each consecutive anchor pair (0-20).";
  }

  function currentVolumeCount() {
    if (hasAnchorIndices()) {
      if (currentAnchorTraversalMode() === "direct") {
        var segs = Math.max(0, anchorIndicesActive.length - 1);
        return segs * (currentNPoints() + 1) + 1;
      }
      if (anchorTrajXY && anchorTrajXY.length) return anchorTrajXY.length;
      return anchorIndicesActive.length;
    }
    return currentNPoints();
  }

  function updateGenerateVolumesButtonLabel() {
    if (!btnGenerateVolumesEl) return;
    var n = currentVolumeCount();
    if (lastTrajectoryVolumeCacheId) {
      var saveCount = generatedTrajectoryVolumeCount > 0 ? generatedTrajectoryVolumeCount : n;
      btnGenerateVolumesEl.textContent = "Save " + saveCount + " volumes";
      return;
    }
    btnGenerateVolumesEl.textContent = "Generate " + n + " volumes";
  }

  function setGeneratedVolumesState(cacheId, nGenerated) {
    lastTrajectoryVolumeCacheId = String(cacheId || "");
    generatedTrajectoryVolumeCount = Number.isFinite(nGenerated) ? Math.max(0, nGenerated) : 0;
    if (!lastTrajectoryVolumeCacheId && trajVolColumn) {
      trajVolColumn.innerHTML = "";
    }
    updateGenerateVolumesButtonLabel();
  }

  function redrawTrajectoryOverlay() {
    if (!gd || !gd.data || !gd.data.length) return;
    if (anchorIndicesActive && anchorTrajXY && anchorTrajXY.length >= 2) {
      drawAnchorPolyline(anchorTrajXY);
      return;
    }
    if (!hasAnchorIndices() && editableTrajXY && editableTrajXY.length >= 2) {
      drawEditableTrajectory(editableTrajXY);
      return;
    }
    if (!startXY || !endXY) return;
    drawLineAndEndpoints(startXY[0], startXY[1], endXY[0], endXY[1]);
  }

  function clearTrajectoryZPanel() {
    lastZPathTxt = null;
    if (trajZNums) trajZNums.textContent = "";
    if (trajZPre) trajZPre.textContent = "";
    if (trajZPanel) trajZPanel.hidden = true;
  }

  function clearVolumeMontage() {
    volumeGeneration++;
    coordsRequestId++;
    if (trajVolColumn) trajVolColumn.innerHTML = "";
    setGeneratedVolumesState("", 0);
    editableTrajXY = null;
    clearTrajectoryZPanel();
    if (gd && gd.data && gd.data.length > 3) {
      try { Plotly.deleteTraces(gd, [3]); } catch (e) {}
    }
  }

  /** Drop numbered trajectory markers without touching thumbnails or in-flight decode. */
  function removeTrajectorySampleMarkersOnly() {
    if (gd && gd.data && gd.data.length > 3) {
      try { Plotly.deleteTraces(gd, [3]); } catch (e) {}
    }
  }

  var dragRestyleRaf = 0;
  var editableDragRaf = 0;

  function cancelDragRestyleRaf() {
    if (dragRestyleRaf) {
      cancelAnimationFrame(dragRestyleRaf);
      dragRestyleRaf = 0;
    }
  }

  /** Update line + endpoint traces (indices 1 and 2) without touching sample markers. */
  function restyleTrajectoryLineOnly() {
    if (!gd || !gd.data || gd.data.length < 3) return;
    try {
      var xup = [[startXY[0], endXY[0]]];
      var yup = [[startXY[1], endXY[1]]];
      Plotly.restyle(gd, { x: xup, y: yup }, [1]);
      Plotly.restyle(gd, { x: xup, y: yup }, [2]);
    } catch (e) {
      console.error(e);
    }
  }

  function scheduleTrajectoryRestyle() {
    if (dragRestyleRaf) return;
    dragRestyleRaf = requestAnimationFrame(function() {
      dragRestyleRaf = 0;
      restyleTrajectoryLineOnly();
    });
  }

  function cancelEditableDragRaf() {
    if (editableDragRaf) {
      cancelAnimationFrame(editableDragRaf);
      editableDragRaf = 0;
    }
  }

  function restyleEditableTrajectoryFromState() {
    if (!gd || !gd.data || !editableTrajXY || editableTrajXY.length < 2) return;
    var xs = editableTrajXY.map(function(p) { return p[0]; });
    var ys = editableTrajXY.map(function(p) { return p[1]; });
    if (gd.data.length > 1 && gd.data[1] && typeof gd.data[1].mode === "string" &&
        gd.data[1].mode.indexOf("lines+markers+text") >= 0) {
      try {
        Plotly.restyle(gd, { x: [xs], y: [ys] }, [1]);
        return;
      } catch (e) {}
    }
    redrawTrajectoryOverlay();
  }

  function scheduleEditableTrajectoryRestyle() {
    if (editableDragRaf) return;
    editableDragRaf = requestAnimationFrame(function() {
      editableDragRaf = 0;
      restyleEditableTrajectoryFromState();
    });
  }

  function fmtCoord(v) {
    if (v == null) return "?";
    var a = Math.abs(v);
    if (a === 0) return "0";
    if (a >= 100) return v.toFixed(0);
    if (a >= 1) return v.toFixed(1);
    return v.toPrecision(3);
  }
  function fmtXY(xy) {
    if (!xy) return "(?, ?)";
    return "(" + fmtCoord(xy[0]) + ", " + fmtCoord(xy[1]) + ")";
  }

  function formatZPathDisplayFromZTraj(zTraj) {
    if (!zTraj || !zTraj.length) return "";
    var nCols = Array.isArray(zTraj[0]) ? zTraj[0].length : 0;
    var zHead = [];
    for (var j = 0; j < nCols; j++) zHead.push("z" + j);
    var body = zTraj.map(function(row) {
      return row.map(function(v) {
        var x = Number(v);
        if (!Number.isFinite(x)) return String(v);
        return x.toFixed(3);
      }).join(" ");
    }).join("\n");
    return (zHead.length ? zHead.join(" ") + "\n" : "") + body;
  }

  function formatZPathGutterFromPayload(j) {
    if (!j) return "";
    var traj_xy = (j.traj_xy && j.traj_xy.length) ? j.traj_xy : [];
    var pidx = (j.traj_particle_indices && j.traj_particle_indices.length)
      ? j.traj_particle_indices
      : [];
    var n = 0;
    if (j.z_traj && j.z_traj.length) n = j.z_traj.length;
    else if (traj_xy.length) n = traj_xy.length;
    if (!n) return "";
    var nDigits = String(n).length;
    var lines = [];
    var xLabel = j.xcol || "x";
    var yLabel = j.ycol || "y";
    lines.push("pt  particle  (" + xLabel + ", " + yLabel + ")");
    for (var i = 0; i < n; i++) {
      var idxTxt = String(i + 1);
      while (idxTxt.length < nDigits) idxTxt = " " + idxTxt;
      var p = (pidx[i] != null) ? String(pidx[i]) : "—";
      var xy = traj_xy[i];
      var xyTxt = "(?, ?)";
      if (xy && xy.length >= 2) {
        xyTxt = "(" + fmtCoord(xy[0]) + ", " + fmtCoord(xy[1]) + ")";
      }
      lines.push(idxTxt + "  #" + p + "  " + xyTxt);
    }
    return lines.join("\n");
  }

  function axisDisplayName(col) {
    return covariateDisplayMap[col] || col;
  }

  function applyVolumeCellAnnotations(entries, j) {
    var mode = (j && j.mode) || trajectoryMode;
    var traj_xy = (j && j.traj_xy) || [];
    var xcol = (j && j.xcol) || sx.value;
    var ycol = (j && j.ycol) || sy.value;
    var xdn = axisDisplayName(xcol);
    var ydn = axisDisplayName(ycol);
    for (var i = 0; i < entries.length; i++) {
      var ent = entries[i];
      var xy = traj_xy[i];
      ent.inset.classList.remove("cryo-traj-vol-inset--show");
      ent.insetImg.src = "";
      if (!xy || xy.length < 2) {
        ent.meta.textContent = "";
        continue;
      }
      if ((mode === "nearest" || mode === "graph") && j.traj_particle_indices && j.traj_particle_indices[i] != null) {
        var pidx = j.traj_particle_indices[i];
        var row = j.traj_rows && j.traj_rows[i] != null ? j.traj_rows[i] : "";
        ent.meta.textContent = "#" + pidx + (row !== "" ? " · row " + row : "")
          + "\n" + xdn + ": " + fmtCoord(xy[0]) + "\n" + ydn + ": " + fmtCoord(xy[1]);
        var thumb = j.particle_thumbs && j.particle_thumbs[i];
        if (thumb) {
          ent.insetImg.src = "data:image/jpeg;base64," + thumb;
          ent.inset.classList.add("cryo-traj-vol-inset--show");
        }
      } else {
        ent.meta.textContent = xdn + ": " + fmtCoord(xy[0]) + "\n" + ydn + ": " + fmtCoord(xy[1]);
      }
    }
  }

  function syncScatterControlsAlignment() {
    if (!scatterControlsCard || !scatterControlsSide || !scatterPlotStack) return;
    var minGapPx = 8;
    var stackRect = scatterPlotStack.getBoundingClientRect();
    var sideRect = scatterControlsSide.getBoundingClientRect();
    if (!stackRect || stackRect.height <= 0 || !sideRect) {
      scatterControlsCard.style.marginTop = "0px";
      return;
    }
    var previousMargin = scatterControlsCard.style.marginTop;
    scatterControlsCard.style.marginTop = "0px";
    var baseTop = scatterControlsCard.offsetTop;
    var cardHalf = scatterControlsCard.offsetHeight * 0.5;
    scatterControlsCard.style.marginTop = previousMargin;
    var plotCenterLocal = (stackRect.top + stackRect.height * 0.5) - sideRect.top;
    var desiredTop = plotCenterLocal - cardHalf;
    var px = Math.max(minGapPx, desiredTop - baseTop);
    scatterControlsCard.style.marginTop = Math.round(px) + "px";
  }

  if (plotlyStackEl && typeof ResizeObserver !== "undefined") {
    var plotStackResizeObserver = new ResizeObserver(function() {
      window.requestAnimationFrame(function() {
        syncScatterControlsAlignment();
        if (!gd || !gd.data || !gd.data.length) return;
        try { Plotly.Plots.resize(gd); } catch (e) {}
      });
    });
    plotStackResizeObserver.observe(plotlyStackEl);
  }

  /* ── Convert a DOM click on the plot to data coordinates ── */
  function plotClickToData(evt, clamp) {
    var dragLayer = gd.querySelector(".nsewdrag");
    if (!dragLayer) return null;
    var rect = dragLayer.getBoundingClientRect();
    var fracX = (evt.clientX - rect.left) / rect.width;
    var fracY = (evt.clientY - rect.top) / rect.height;
    if (clamp) {
      fracX = Math.max(0, Math.min(1, fracX));
      fracY = Math.max(0, Math.min(1, fracY));
    } else if (fracX < 0 || fracX > 1 || fracY < 0 || fracY > 1) {
      return null;
    }
    var fl = gd._fullLayout;
    if (!fl || !fl.xaxis || !fl.yaxis) return null;
    var xr = fl.xaxis.range;
    var yr = fl.yaxis.range;
    var dataX = xr[0] + fracX * (xr[1] - xr[0]);
    var dataY = yr[1] - fracY * (yr[1] - yr[0]);
    return [dataX, dataY];
  }

  function dataToPixelXY(dataX, dataY) {
    var fl = gd._fullLayout;
    if (!fl || !fl.xaxis || !fl.yaxis) return null;
    var dragLayer = gd.querySelector(".nsewdrag");
    if (!dragLayer) return null;
    var rect = dragLayer.getBoundingClientRect();
    var xr = fl.xaxis.range;
    var yr = fl.yaxis.range;
    if (xr[1] === xr[0] || yr[1] === yr[0]) return null;
    var fracX = (dataX - xr[0]) / (xr[1] - xr[0]);
    var fracY = (yr[1] - dataY) / (yr[1] - yr[0]);
    return [rect.left + fracX * rect.width, rect.top + fracY * rect.height];
  }

  function pixelDist(ax, ay, bx, by) {
    var ddx = ax - bx, ddy = ay - by;
    return Math.sqrt(ddx * ddx + ddy * ddy);
  }

  var DRAG_HANDLE_PX = 18;

  function restyleLineAndHandles() {
    if (!gd || !gd.data || gd.data.length < 3) return;
    try {
      if (gd.data.length > 3) Plotly.deleteTraces(gd, [3]);
    } catch (e) {}
    restyleTrajectoryLineOnly();
  }

  function setTrajectoryModesLockedForAnchors(on) {
    // Keep signature for call sites; mode buttons now remain active with anchors loaded.
    syncTraversalModeButtonStates();
  }

  function directTraversalAllowedForAxes(xc, yc) {
    function zCol(c) { return /^z[0-9]+$/.test(c); }
    function pcCol(c) { return /^PC[0-9]+$/.test(c); }
    return (zCol(xc) || pcCol(xc)) && (zCol(yc) || pcCol(yc));
  }

  function syncTraversalModeButtonStates() {
    var hasAnchors = hasAnchorIndices();
    syncTraversalSettingsText();
    if (hasAnchors) {
      if (btnModeDirectEl) btnModeDirectEl.disabled = false;
      if (btnModeNearestEl) btnModeNearestEl.disabled = false;
      if (trajectoryMode !== "direct" && trajectoryMode !== "graph") {
        trajectoryMode = "direct";
      }
      if (btnModeDirectEl) {
        btnModeDirectEl.classList.toggle("active", trajectoryMode === "direct");
        btnModeDirectEl.classList.toggle("btn-secondary", trajectoryMode !== "direct");
      }
      if (btnModeNearestEl) {
        btnModeNearestEl.classList.toggle("active", trajectoryMode === "graph");
        btnModeNearestEl.classList.toggle("btn-secondary", trajectoryMode !== "graph");
      }
      return;
    }
    var allowDirect = directTraversalAllowedForAxes(sx.value, sy.value);
    if (btnModeDirectEl) {
      btnModeDirectEl.disabled = !allowDirect;
    }
    if (btnModeNearestEl) {
      btnModeNearestEl.disabled = false;
    }
    if (!allowDirect && trajectoryMode === "direct") {
      trajectoryMode = "nearest";
      if (btnModeDirectEl) {
        btnModeDirectEl.classList.remove("active");
        btnModeDirectEl.classList.add("btn-secondary");
      }
      if (btnModeNearestEl) {
        btnModeNearestEl.classList.add("active");
        btnModeNearestEl.classList.remove("btn-secondary");
      }
      if (startXY && endXY) {
        fetchTrajectoryCoords();
      }
    } else {
      if (btnModeDirectEl) {
        btnModeDirectEl.classList.toggle("active", trajectoryMode === "direct");
        btnModeDirectEl.classList.toggle("btn-secondary", trajectoryMode !== "direct");
      }
      if (btnModeNearestEl) {
        btnModeNearestEl.classList.toggle("active", trajectoryMode === "nearest");
        btnModeNearestEl.classList.toggle("btn-secondary", trajectoryMode !== "nearest");
      }
    }
  }

  function syncAnchorIndicesUI() {
    var on = anchorIndicesActive && anchorIndicesActive.length >= 2;
    if (trajAnchorPick) trajAnchorPick.hidden = !!on;
    if (trajAnchorManualRow) trajAnchorManualRow.hidden = !!on;
    if (on) {
      closeFileBrowser();
    }
    if (btnAnchorRemove) {
      btnAnchorRemove.classList.toggle("cryo-traj-remove--visible", !!on);
      btnAnchorRemove.setAttribute("aria-hidden", on ? "false" : "true");
    }
  }

  function clearAnchorMode() {
    anchorIndicesActive = null;
    anchorTrajXY = null;
    editableTrajXY = null;
    if (nPointsEl) nPointsEl.disabled = false;
    setTrajectoryModesLockedForAnchors(false);
    updateGenerateVolumesButtonLabel();
    syncAnchorIndicesUI();
    closeFileBrowser();
  }

  function onTrajectoryIndicesLoaded(res) {
    anchorIndicesActive = res.anchor_indices ? res.anchor_indices.slice() : null;
    anchorTrajXY = res.traj_xy ? res.traj_xy.slice() : null;
    editableTrajXY = null;
    if (res && (res.mode === "direct" || res.mode === "graph")) {
      trajectoryMode = res.mode;
    }
    if (nPointsEl) nPointsEl.disabled = false;
    setTrajectoryModesLockedForAnchors(true);
    syncPointsSelectorRange();
    if (nPointsEl && currentAnchorTraversalMode() === "direct") {
      nPointsEl.value = "0";
    }
    setGeneratedVolumesState("", 0);
    updateGenerateVolumesButtonLabel();
    syncAnchorIndicesUI();
    applyTrajectoryCoordsPayload(res);
    redrawTrajectoryOverlay();
    setTrajStatus("Index trajectory ready — press Generate to decode volumes.", false);
  }

  function makeLinearEditablePath(firstPt, lastPt, nPoints) {
    if (!firstPt || !lastPt) return null;
    var n = Math.max(2, parseInt(nPoints, 10) || 2);
    var out = [];
    var x0 = Number(firstPt[0]), y0 = Number(firstPt[1]);
    var x1 = Number(lastPt[0]), y1 = Number(lastPt[1]);
    if (!Number.isFinite(x0) || !Number.isFinite(y0) || !Number.isFinite(x1) || !Number.isFinite(y1)) {
      return null;
    }
    for (var i = 0; i < n; i++) {
      var t = (n <= 1) ? 0 : (i / (n - 1));
      out.push([x0 + t * (x1 - x0), y0 + t * (y1 - y0)]);
    }
    return out;
  }

  /* ── Mode toggle (also refreshes latent z via coords API) ── */
  function setMode(mode) {
    var prevMode = trajectoryMode;
    var hasAnchors = hasAnchorIndices();
    if (hasAnchors) {
      if (mode !== "direct" && mode !== "graph") return;
    } else if (mode !== "direct" && mode !== "nearest") {
      return;
    }
    if (!hasAnchors && mode === "direct" && btnModeDirectEl && btnModeDirectEl.disabled) return;
    trajectoryMode = mode;
    if (!hasAnchors && prevMode === "nearest" && mode === "direct") {
      var nLinear = (editableTrajXY && editableTrajXY.length >= 2)
        ? editableTrajXY.length
        : currentNPoints();
      var firstPt = (editableTrajXY && editableTrajXY.length >= 2)
        ? editableTrajXY[0]
        : startXY;
      var lastPt = (editableTrajXY && editableTrajXY.length >= 2)
        ? editableTrajXY[editableTrajXY.length - 1]
        : endXY;
      var linearPath = makeLinearEditablePath(firstPt, lastPt, nLinear);
      if (linearPath && linearPath.length >= 2) {
        editableTrajXY = linearPath;
      } else {
        editableTrajXY = null;
      }
    }
    syncTraversalModeButtonStates();
    updateGenerateVolumesButtonLabel();
    if (hasAnchors) {
      fetchTrajectoryCoords();
      return;
    }
    if (startXY && endXY) {
      redrawTrajectoryOverlay();
      fetchTrajectoryCoords();
    }
  }
  btnModeDirectEl.addEventListener("click", function() { setMode("direct"); });
  btnModeNearestEl.addEventListener("click", function() {
    setMode(hasAnchorIndices() ? "graph" : "nearest");
  });
  if (btnGenerateVolumesEl) {
    btnGenerateVolumesEl.addEventListener("click", function() {
      if (lastTrajectoryVolumeCacheId) {
        saveTrajectoryVolumes();
      } else {
        generateTrajectory();
      }
    });
  }
  if (btnSaveZpathAs) {
    btnSaveZpathAs.addEventListener("click", function() {
      if (!lastZPathTxt) return;
      if (fbSaveName) fbSaveName.value = "z-path.txt";
      openFileBrowser("", "save-zpath");
    });
  }
  if (btnSaveZpathWorkdir) {
    btnSaveZpathWorkdir.addEventListener("click", function() {
      if (!lastZPathTxt) return;
      fetch(b.urls.apiTrajectorySaveZpath, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ z_path_txt: lastZPathTxt })
      })
        .then(function(r) { return r.json().then(function(j) { return { ok: r.ok, j: j }; }); })
        .then(function(res) {
          if (!res.ok) {
            setTrajStatus(res.j.error || "Could not save z-path.txt.", false);
            return;
          }
          var p = res.j.path || (expWorkdir ? expWorkdir + "/z-path.txt" : "z-path.txt");
          setTrajStatus("Saved " + p, false);
        })
        .catch(function(err) {
          console.error(err);
          setTrajStatus("Save failed.", false);
        });
    });
  }
  function postTrajectoryIndices(url, statusMsg, errPrefix) {
    closeFileBrowser();
    setTrajStatus(statusMsg, false);
    var anchorMode = currentAnchorTraversalMode();
    var nPointsReq = requestedAnchorNPoints(anchorMode);
    fetch(url, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        x: sx.value,
        y: sy.value,
        mode: anchorMode,
        n_points: nPointsReq,
        max_neighbors: currentMaxNeighbors(),
        avg_neighbors: currentAvgNeighbors()
      })
    })
      .then(function(r) { return r.json().then(function(j) { return { ok: r.ok, j: j }; }); })
      .then(function(res) {
        if (!res.ok) {
          setTrajStatus(res.j.error || (errPrefix + " failed."), false);
          return;
        }
        onTrajectoryIndicesLoaded(res.j);
      })
      .catch(function(err) {
        console.error(err);
        setTrajStatus(errPrefix + " request failed.", false);
      });
  }

  function parseManualAnchorIndices(raw) {
    var txt = String(raw || "").trim();
    if (!txt) return [];
    txt = txt.replace(/[\[\]\(\)]/g, " ");
    var toks = txt.split(/[\s,;]+/).filter(Boolean);
    var out = [];
    for (var i = 0; i < toks.length; i++) {
      var t = toks[i];
      if (!/^-?[0-9]+$/.test(t)) return null;
      out.push(parseInt(t, 10));
    }
    return out;
  }

  function loadManualAnchorIndices() {
    closeFileBrowser();
    var parsed = parseManualAnchorIndices(anchorIndicesManualEl && anchorIndicesManualEl.value);
    if (!parsed || parsed.length < 2) {
      setTrajStatus("Enter at least two integer anchor indices (comma or space separated).", false);
      return;
    }
    anchorIndicesActive = parsed.slice();
    anchorTrajXY = null;
    if (nPointsEl) nPointsEl.disabled = false;
    if (nPointsEl && currentAnchorTraversalMode() === "direct") {
      nPointsEl.value = "0";
    }
    setTrajectoryModesLockedForAnchors(true);
    updateGenerateVolumesButtonLabel();
    syncAnchorIndicesUI();
    fetchTrajectoryCoords();
  }

  function syncFileBrowserModeUI() {
    var saveMode = fbMode === "save-zpath";
    if (fbSaveRow) fbSaveRow.hidden = !saveMode;
  }

  function openFileBrowser(startDir, mode) {
    if (!trajFileBrowserPanel || !serverFileBrowser) return;
    fbMode = (mode === "save-zpath") ? "save-zpath" : "import";
    trajFileBrowserPanel.hidden = false;
    trajFileBrowserPanel.setAttribute("aria-hidden", "false");
    syncFileBrowserModeUI();
    loadFileBrowserDir(startDir || "");
  }

  function closeFileBrowser() {
    if (trajFileBrowserPanel) {
      trajFileBrowserPanel.hidden = true;
      trajFileBrowserPanel.setAttribute("aria-hidden", "true");
    }
    fbCurrentDir = null;
    fbMode = "import";
    syncFileBrowserModeUI();
  }

  function normalizedSaveFilename() {
    var name = String((fbSaveName && fbSaveName.value) || "").trim();
    if (!name) name = "z-path.txt";
    if (name.indexOf("/") >= 0 || name.indexOf("\\") >= 0) return null;
    if (!/\.txt$/i.test(name)) name = name + ".txt";
    return name;
  }

  function saveZPathViaBrowser() {
    if (!lastZPathTxt) {
      setTrajStatus("No latent z path to save yet.", false);
      return;
    }
    if (!fbCurrentDir) {
      setTrajStatus("Choose a destination folder first.", false);
      return;
    }
    var fname = normalizedSaveFilename();
    if (!fname) {
      setTrajStatus("Filename cannot contain path separators.", false);
      return;
    }
    var outPath = fbCurrentDir + "/" + fname;
    setTrajStatus("Saving latent z path…", false);
    fetch(b.urls.apiTrajectorySaveZpath, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        z_path_txt: lastZPathTxt,
        out_path: outPath
      })
    })
      .then(function(r) { return r.json().then(function(j) { return { ok: r.ok, j: j }; }); })
      .then(function(res) {
        if (!res.ok) {
          setTrajStatus(res.j.error || "Could not save z-path.txt.", false);
          return;
        }
        closeFileBrowser();
        setTrajStatus("Saved " + (res.j.path || outPath), false);
      })
      .catch(function(err) {
        console.error(err);
        setTrajStatus("Save request failed.", false);
      });
  }

  function loadFileBrowserDir(dir) {
    if (!fbList || !fbPath) return;
    fbList.innerHTML = "<li class='cryo-file-browser-empty'>Loading…</li>";
    var q = dir ? ("?dir=" + encodeURIComponent(dir)) : "";
    fetch(b.urls.apiListServerFiles + q)
      .then(function(r) { return r.json(); })
      .then(function(j) {
        if (!j.ok) {
          fbList.innerHTML = "<li class='cryo-file-browser-empty'>" + (j.error || "Error") + "</li>";
          return;
        }
        fbCurrentDir = j.dir;
        fbPath.textContent = j.dir;
        fbPath.title = j.dir;
        fbList.innerHTML = "";
        if (!j.entries || !j.entries.length) {
          fbList.innerHTML = "<li class='cryo-file-browser-empty'>No entries here</li>";
          return;
        }
        j.entries.forEach(function(ent) {
          var li = document.createElement("li");
          var icon = document.createElement("span");
          icon.className = "fb-icon";
          icon.textContent = ent.type === "dir" ? "\uD83D\uDCC1" : "\uD83D\uDCC4";
          var name = document.createElement("span");
          name.className = "fb-name";
          name.textContent = ent.name;
          li.appendChild(icon);
          li.appendChild(name);
          if (ent.type === "dir") {
            li.addEventListener("click", function() {
              loadFileBrowserDir(j.dir + "/" + ent.name);
            });
          } else {
            li.addEventListener("click", function() {
              var fullPath = j.dir + "/" + ent.name;
              if (fbMode === "save-zpath") {
                if (fbSaveName) fbSaveName.value = ent.name;
              } else {
                closeFileBrowser();
                importServerPklFile(fullPath);
              }
            });
          }
          fbList.appendChild(li);
        });
      })
      .catch(function() {
        fbList.innerHTML = "<li class='cryo-file-browser-empty'>Could not list directory</li>";
      });
  }

  function importServerPklFile(serverPath) {
    setTrajStatus("Loading anchor indices…", false);
    var anchorMode = currentAnchorTraversalMode();
    var nPointsReq = requestedAnchorNPoints(anchorMode);
    fetch(b.urls.apiTrajectoryImportAnchors, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        server_path: serverPath,
        x: sx.value,
        y: sy.value,
        mode: anchorMode,
        n_points: nPointsReq,
        max_neighbors: currentMaxNeighbors(),
        avg_neighbors: currentAvgNeighbors()
      })
    })
      .then(function(r) { return r.json().then(function(j) { return { ok: r.ok, j: j }; }); })
      .then(function(res) {
        if (!res.ok) {
          setTrajStatus(res.j.error || "Import failed.", false);
          return;
        }
        onTrajectoryIndicesLoaded(res.j);
      })
      .catch(function(err) {
        console.error(err);
        setTrajStatus("Import request failed.", false);
      });
  }

  if (btnAnchorImportPkl) {
    btnAnchorImportPkl.addEventListener("click", function() {
      openFileBrowser("", "import");
    });
  }
  if (fbUp) {
    fbUp.addEventListener("click", function() {
      if (fbCurrentDir) {
        var parent = fbCurrentDir.replace(/\/[^\/]*\/?$/, "") || "/";
        loadFileBrowserDir(parent);
      }
    });
  }
  if (fbCancel) {
    fbCancel.addEventListener("click", function() {
      closeFileBrowser();
    });
  }
  if (fbSaveBtn) {
    fbSaveBtn.addEventListener("click", function() {
      saveZPathViaBrowser();
    });
  }
  if (fbSaveName) {
    fbSaveName.addEventListener("keydown", function(evt) {
      if (evt.key === "Enter") {
        evt.preventDefault();
        saveZPathViaBrowser();
      }
    });
  }
  if (btnAnchorKmeans) {
    btnAnchorKmeans.addEventListener("click", function() {
      postTrajectoryIndices(
        b.urls.apiTrajectoryKmeansCenters,
        "Loading k-means centers…",
        "K-means centers"
      );
    });
  }
  if (btnAnchorRandom) {
    btnAnchorRandom.addEventListener("click", function() {
      postTrajectoryIndices(
        b.urls.apiTrajectoryRandomIndices,
        "Sampling random indices…",
        "Random indices"
      );
    });
  }
  if (btnAnchorManualLoad) {
    btnAnchorManualLoad.addEventListener("click", loadManualAnchorIndices);
  }
  if (anchorIndicesManualEl) {
    anchorIndicesManualEl.addEventListener("keydown", function(evt) {
      if (evt.key === "Enter") {
        evt.preventDefault();
        loadManualAnchorIndices();
      }
    });
  }
  if (btnAnchorRemove) {
    btnAnchorRemove.addEventListener("click", function() {
      clearAnchorMode();
      clearTrajectoryZPanel();
      initDefaultTrajectory(scatterLoadGeneration);
    });
  }

  /* ── Before reloading scatter (axes / color): cancel in-flight requests, clear stale 2D state ── */
  function preparePlotReload() {
    volumeGeneration++;
    coordsRequestId++;
    dragEndpoint = null;
    dragPointIndex = -1;
    startXY = null;
    endXY = null;
    anchorTrajXY = null;
    editableTrajXY = null;
    setTrajStatus("", false);
    if (!anchorIndicesActive) {
      trajVolColumn.innerHTML = "";
      clearTrajectoryZPanel();
    }
  }

  function resetTrajectory() {
    clearAnchorMode();
    initDefaultTrajectory(scatterLoadGeneration);
  }
  btnResetEl.addEventListener("click", resetTrajectory);

  /* ── Overlay traces (markers + line) on the scatter ── */
  function removeOverlayTraces() {
    if (!gd || !gd.data) return;
    var toRemove = [];
    for (var i = gd.data.length - 1; i >= 1; i--) toRemove.push(i);
    if (toRemove.length) {
      try { Plotly.deleteTraces(gd, toRemove); } catch (e) {}
    }
  }

  function drawAnchorPolyline(traj_xy) {
    if (!gd || !gd.data || !traj_xy || traj_xy.length < 2) return;
    removeOverlayTraces();
    var xs = traj_xy.map(function(p) { return p[0]; });
    var ys = traj_xy.map(function(p) { return p[1]; });
    var showLabels = traj_xy.length <= 40;
    var labels = showLabels ? traj_xy.map(function(_, i) { return "<b>" + (i + 1) + "</b>"; }) : [];
    Plotly.addTraces(gd, [{
      x: xs,
      y: ys,
      mode: showLabels ? "lines+markers+text" : "lines+markers",
      line: { color: "#c4703a", width: 2.5, dash: "dot" },
      marker: { size: 11, color: "#243b53", line: { width: 2, color: "#fff" } },
      text: labels,
      textposition: "top center",
      textfont: { size: 19, color: "#243b53", family: "Barlow, sans-serif" },
      hoverinfo: "skip",
      showlegend: false
    }]);
  }

  function drawLineAndEndpoints(sx, sy, ex, ey) {
    removeOverlayTraces();
    Plotly.addTraces(gd, [
      {
        x: [sx, ex], y: [sy, ey],
        mode: "lines",
        line: { color: "#c4703a", width: 2.5, dash: "dot" },
        hoverinfo: "skip",
        showlegend: false
      },
      {
        x: [sx, ex], y: [sy, ey],
        mode: "markers",
        marker: { size: 12, color: ["#c4703a", "#243b53"],
                  line: { width: 2, color: "#fff" } },
        hoverinfo: "skip",
        showlegend: false
      }
    ]);
  }

  function drawTrajectoryPoints(points) {
    if (!gd || !gd.data) return;
    var xs = points.map(function(p) { return p[0]; });
    var ys = points.map(function(p) { return p[1]; });
    var labels = [];
    for (var i = 0; i < points.length; i++) labels.push("<b>" + (i + 1) + "</b>");
    Plotly.addTraces(gd, [{
      x: xs, y: ys,
      mode: "markers+text",
      text: labels,
      textposition: "top center",
      textfont: { size: 19, color: "#243b53", family: "Barlow, sans-serif" },
      marker: { size: 10, color: "#243b53", symbol: "diamond",
                line: { width: 1, color: "#fff" } },
      hoverinfo: "skip",
      showlegend: false
    }]);
  }

  function drawEditableTrajectory(points) {
    if (!gd || !gd.data || !points || points.length < 2) return;
    removeOverlayTraces();
    var xs = points.map(function(p) { return p[0]; });
    var ys = points.map(function(p) { return p[1]; });
    var labels = [];
    for (var i = 0; i < points.length; i++) labels.push("<b>" + (i + 1) + "</b>");
    Plotly.addTraces(gd, [{
      x: xs,
      y: ys,
      mode: "lines+markers+text",
      line: { color: "#c4703a", width: 2.2, dash: "dot" },
      text: labels,
      textposition: "top center",
      textfont: { size: 17, color: "#243b53", family: "Barlow, sans-serif" },
      marker: { size: 10, color: "#243b53", symbol: "diamond", line: { width: 1, color: "#fff" } },
      hoverinfo: "skip",
      showlegend: false
    }]);
  }

  function applyTrajectoryCoordsPayload(j) {
    if (!j || !j.ok) return;
    removeTrajectorySampleMarkersOnly();
    if (hasAnchorIndices() && j.traj_xy) {
      anchorTrajXY = j.traj_xy;
      editableTrajXY = null;
    } else if (j.traj_xy && j.traj_xy.length) {
      editableTrajXY = j.traj_xy.map(function(p) { return [Number(p[0]), Number(p[1])]; });
      redrawTrajectoryOverlay();
    }
    if (j.z_path_txt && trajZPre && trajZPanel) {
      lastZPathTxt = j.z_path_txt;
      var hasZTraj = j.z_traj && j.z_traj.length;
      trajZPre.textContent = hasZTraj
        ? formatZPathDisplayFromZTraj(j.z_traj)
        : j.z_path_txt;
      if (trajZNums) trajZNums.textContent = hasZTraj
        ? formatZPathGutterFromPayload(j)
        : "";
      trajZPanel.hidden = false;
    } else {
      clearTrajectoryZPanel();
    }
  }

  function scheduleCoordsAfterDrag() {
    if (anchorIndicesActive && anchorIndicesActive.length >= 2) return;
    if (coordsDragTimer) {
      clearTimeout(coordsDragTimer);
      coordsDragTimer = 0;
    }
    coordsDragTimer = setTimeout(function() {
      coordsDragTimer = 0;
      fetchTrajectoryCoords();
    }, 300);
  }

  function fetchTrajectoryCoords() {
    setGeneratedVolumesState("", 0);
    if (hasAnchorIndices()) {
      var myIdA = ++coordsRequestId;
      setTrajStatus("Computing latent z…", false);
      fetch(b.urls.apiTrajectoryCoords, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          anchor_indices: anchorIndicesActive,
          mode: currentAnchorTraversalMode(),
          n_points: currentNPoints(),
          max_neighbors: currentMaxNeighbors(),
          avg_neighbors: currentAvgNeighbors(),
          x: sx.value,
          y: sy.value
        })
      })
        .then(function(r) { return r.json().then(function(j) { return { ok: r.ok, j: j }; }); })
        .then(function(res) {
          if (myIdA !== coordsRequestId) return;
          if (!res.ok) {
            setTrajStatus(res.j.error || "Could not compute latent z.", false);
            return;
          }
          applyTrajectoryCoordsPayload(res.j);
          redrawTrajectoryOverlay();
          setTrajStatus("Latent z ready — press Generate to decode volumes.", false);
        })
        .catch(function(err) {
          if (myIdA !== coordsRequestId) return;
          console.error(err);
          setTrajStatus("Latent z request failed.", false);
        });
      return;
    }
    if (!startXY || !endXY) return;
    var myId = ++coordsRequestId;
    setTrajStatus("Computing latent z…", false);
    var payload = {
      mode: trajectoryMode,
      start: startXY,
      end: endXY,
      x: sx.value,
      y: sy.value,
      n_points: currentNPoints()
    };
    if (editableTrajXY && editableTrajXY.length >= 2) {
      payload.traj_xy = editableTrajXY;
    }
    fetch(b.urls.apiTrajectoryCoords, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload)
    })
      .then(function(r) { return r.json().then(function(j) { return { ok: r.ok, j: j }; }); })
      .then(function(res) {
        if (myId !== coordsRequestId) return;
        if (!res.ok) {
          setTrajStatus(res.j.error || "Could not compute latent z.", false);
          return;
        }
        applyTrajectoryCoordsPayload(res.j);
        setTrajStatus("Latent z ready — press Generate to decode volumes.", false);
      })
      .catch(function(err) {
        if (myId !== coordsRequestId) return;
        console.error(err);
        setTrajStatus("Latent z request failed.", false);
      });
  }

  /* ── Volume column ── */
  function buildVolumeColumn(n) {
    trajVolColumn.innerHTML = "";
    var entries = [];
    for (var i = 0; i < n; i++) {
      var cell = document.createElement("div");
      cell.className = "cryo-traj-vol-cell";
      var mainImg = document.createElement("img");
      mainImg.className = "cryo-traj-vol-main";
      mainImg.src = "";
      mainImg.alt = "Volume " + (i + 1);
      var lbl = document.createElement("span");
      lbl.className = "cryo-traj-vol-label";
      lbl.textContent = String(i + 1);
      var meta = document.createElement("div");
      meta.className = "cryo-traj-vol-meta";
      meta.setAttribute("aria-hidden", "true");
      var inset = document.createElement("div");
      inset.className = "cryo-traj-vol-inset";
      var insetImg = document.createElement("img");
      insetImg.alt = "";
      inset.appendChild(insetImg);
      cell.appendChild(mainImg);
      cell.appendChild(lbl);
      cell.appendChild(meta);
      cell.appendChild(inset);
      trajVolColumn.appendChild(cell);
      entries.push({
        mainImg: mainImg,
        meta: meta,
        inset: inset,
        insetImg: insetImg
      });
    }
    return entries;
  }

  /* ── Generate trajectory volumes ── */
  function generateTrajectory() {
    var useAnchors = hasAnchorIndices();
    if (!useAnchors && (!startXY || !endXY)) return;
    var myGen = ++volumeGeneration;
    var n = currentVolumeCount();
    var volumeEntries = buildVolumeColumn(n);
    setTrajStatus("Decoding " + n + " volumes and rendering (ChimeraX)…", true);
    if (plotStatus) plotStatus.textContent = "Generating trajectory…";

    redrawTrajectoryOverlay();

    var payload = useAnchors
      ? {
          anchor_indices: anchorIndicesActive,
          mode: currentAnchorTraversalMode(),
          n_points: currentNPoints(),
          max_neighbors: currentMaxNeighbors(),
          avg_neighbors: currentAvgNeighbors(),
          x: sx.value,
          y: sy.value
        }
      : {
          mode: trajectoryMode,
          start: startXY,
          end: endXY,
          x: sx.value,
          y: sy.value,
          n_points: currentNPoints()
        };
    if (!useAnchors && editableTrajXY && editableTrajXY.length >= 2) {
      payload.traj_xy = editableTrajXY;
    }
    payload.chimerax_cpus = trajChimeraxCpus;
    fetch(b.urls.apiTrajectoryVolumes, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload)
    })
      .then(function(r) { return r.json().then(function(j) { return { ok: r.ok, j: j }; }); })
      .then(function(res) {
        if (myGen !== volumeGeneration) return;
        if (!res.ok) {
          if (res.j.need_chimerax) {
            window.alert(res.j.error || "Set CHIMERAX_PATH and try again.");
          }
          setTrajStatus(res.j.error || "Volume generation failed.", false);
          if (plotStatus) plotStatus.textContent = "Volume generation failed.";
          return;
        }
        var imgs = res.j.images || [];
        for (var ci = 0; ci < volumeEntries.length && ci < imgs.length; ci++) {
          volumeEntries[ci].mainImg.src = "data:image/png;base64," + imgs[ci];
        }
        applyVolumeCellAnnotations(volumeEntries, res.j);
        applyTrajectoryCoordsPayload(res.j);
        setGeneratedVolumesState(String(res.j.volume_cache_id || ""), imgs.length);
        if (useAnchors && res.j.traj_xy) {
          redrawTrajectoryOverlay();
        }
        setTrajStatus(imgs.length + " volumes generated.", false);
        if (plotStatus) {
          if (useAnchors) {
            plotStatus.textContent = "Index trajectory — click Remove indices or change axes to edit.";
          } else {
            plotStatus.textContent = "Showing trajectory between "
              + fmtXY(startXY) + " and " + fmtXY(endXY)
              + " — drag the Start / End handles to adjust.";
          }
        }
      })
      .catch(function(err) {
        if (myGen !== volumeGeneration) return;
        console.error(err);
        setTrajStatus("Request failed.", false);
        if (plotStatus) plotStatus.textContent = "Request failed.";
      });
  }

  function saveTrajectoryVolumes() {
    if (!lastTrajectoryVolumeCacheId) {
      setTrajStatus("Generate volumes first.", false);
      return;
    }
    var n = generatedTrajectoryVolumeCount > 0 ? generatedTrajectoryVolumeCount : currentVolumeCount();
    var suggestedDir = expWorkdir ? String(expWorkdir).replace(/\/+$/, "") + "/trajectory_volumes" : "";
    var outDir = window.prompt("Choose server folder for saving " + n + " .mrc files:", suggestedDir);
    if (outDir == null) return;
    outDir = String(outDir || "").trim();
    if (!outDir) {
      setTrajStatus("Choose an output folder.", false);
      return;
    }
    setTrajStatus("Saving " + n + " volumes…", true);
    fetch(b.urls.apiTrajectorySaveVolumes, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        volume_cache_id: lastTrajectoryVolumeCacheId,
        out_dir: outDir
      })
    })
      .then(function(r) { return r.json().then(function(j) { return { ok: r.ok, j: j }; }); })
      .then(function(res) {
        if (!res.ok) {
          var msg = res.j && res.j.error ? String(res.j.error) : "Could not save volumes.";
          if (/expired|unknown/i.test(msg)) {
            setGeneratedVolumesState("", 0);
          }
          setTrajStatus(msg, false);
          return;
        }
        var savedN = Number(res.j.n_saved || 0);
        var out = res.j.out_dir || outDir;
        setTrajStatus("Saved " + savedN + " volumes to " + out, false);
      })
      .catch(function(err) {
        console.error(err);
        setTrajStatus("Save request failed.", false);
      });
  }

  function initDefaultTrajectory(gen) {
    var q = "x=" + encodeURIComponent(sx.value) + "&y=" + encodeURIComponent(sy.value);
    fetch(b.urls.apiDefaultTrajectoryEndpoints + "?" + q)
      .then(function(r) { return r.json().then(function(j) { return { ok: r.ok, j: j }; }); })
      .then(function(res) {
        if (gen !== scatterLoadGeneration) return;
        if (!res.ok || !res.j.ok || !res.j.start || !res.j.end) {
          if (plotStatus) {
            plotStatus.textContent = (res.j && res.j.error)
              ? String(res.j.error)
              : "Could not compute default trajectory for these axes.";
          }
          return;
        }
        /* New default path (load / axis change / reset): drop stale thumbnails; drag/mode/n alone do not. */
        clearVolumeMontage();
        startXY = res.j.start.slice();
        endXY = res.j.end.slice();
        editableTrajXY = null;
        redrawTrajectoryOverlay();
        fetchTrajectoryCoords();
        if (plotStatus) {
          plotStatus.textContent = "Default trajectory on long axis — choose mode for z, then Generate for volumes.";
        }
      })
      .catch(function() {
        if (gen !== scatterLoadGeneration) return;
        if (plotStatus) plotStatus.textContent = "Could not load default trajectory.";
      });
  }

  /* ── Scatter plot load ── */
  function loadPlot(opts) {
    var preserveTrajectoryState = !!(opts && opts.preserveTrajectoryState);
    var gen = ++scatterLoadGeneration;
    dragEndpoint = null;
    setRendering(true);
    if (plotStatus) plotStatus.textContent = "";
    var x = sx.value, y = sy.value, c = sc.value;
    var q = "x=" + encodeURIComponent(x)
      + "&y=" + encodeURIComponent(y)
      + "&color=" + encodeURIComponent(c)
      + "&marker_size=2"
      + "&max_points=60000"
      + "&palette=" + encodeURIComponent(selectedTrajPalette());
    fetch(b.urls.apiScatter + "?" + q)
      .then(function(r) {
        return r.text().then(function(text) {
          var j;
          try { j = JSON.parse(text); }
          catch (e) { throw new Error("Scatter API did not return JSON"); }
          if (!r.ok) throw new Error(j.error || r.statusText);
          return j;
        });
      })
      .then(function(fig) {
        if (gen !== scatterLoadGeneration) return;
        fig.layout.title = "";
        fig.layout.margin = Object.assign(fig.layout.margin || {}, {t: 16});
        fig.layout.dragmode = "pan";
        var axFont = {size: 15, family: "Barlow, sans-serif", color: "#243b53"};
        fig.layout.xaxis = Object.assign(fig.layout.xaxis || {}, {title: {text: x, font: axFont}});
        fig.layout.yaxis = Object.assign(fig.layout.yaxis || {}, {title: {text: y, font: axFont}});
        if (/^umap/i.test(x)) fig.layout.xaxis.showticklabels = false;
        if (/^umap/i.test(y)) fig.layout.yaxis.showticklabels = false;
        var opts = {responsive: true};
        var p;
        if (!gd.data || gd.data.length === 0) {
          p = Plotly.newPlot(gd, fig.data, fig.layout, opts);
        } else {
          p = Plotly.react(gd, fig.data, fig.layout, opts);
        }
        function afterPlot() {
          if (gen !== scatterLoadGeneration) return;
          setRendering(false);
          wireUpPlotDrag();
          syncScatterControlsAlignment();
          var canReuseTrajectory = false;
          if (preserveTrajectoryState) {
            if (anchorIndicesActive && anchorIndicesActive.length >= 2) {
              canReuseTrajectory = !!(anchorTrajXY && anchorTrajXY.length >= 2);
            } else {
              canReuseTrajectory = !!(
                (editableTrajXY && editableTrajXY.length >= 2) ||
                (startXY && endXY)
              );
            }
          }
          if (canReuseTrajectory) {
            redrawTrajectoryOverlay();
            requestAnimationFrame(function() {
              try { Plotly.Plots.resize(gd); } catch (e) {}
            });
            return;
          }
          if (anchorIndicesActive && anchorIndicesActive.length >= 2) {
            fetchTrajectoryCoords();
          } else {
            initDefaultTrajectory(gen);
          }
          requestAnimationFrame(function() {
            try { Plotly.Plots.resize(gd); } catch (e) {}
          });
        }
        if (p && typeof p.then === "function") {
          p.then(afterPlot).catch(function(err) {
            setRendering(false);
            if (plotStatus) plotStatus.textContent = err.message || "Could not load scatter.";
          });
        } else {
          afterPlot();
        }
      })
      .catch(function(e) {
        console.error(e);
        setRendering(false);
        if (plotStatus) plotStatus.textContent = e.message || "Could not load scatter.";
      });
  }

  /* ── Drag handles (Start / End) ── */
  function wireUpPlotDrag() {
    if (dragListenersWired) return;
    dragListenersWired = true;
    gd.addEventListener("mousedown", function(evt) {
      if (anchorIndicesActive && anchorIndicesActive.length >= 2) return;
      if (editableTrajXY && editableTrajXY.length >= 2) {
        var pClickM = [evt.clientX, evt.clientY];
        var bestI = -1;
        var bestD = Infinity;
        for (var mi = 0; mi < editableTrajXY.length; mi++) {
          var pp = dataToPixelXY(editableTrajXY[mi][0], editableTrajXY[mi][1]);
          if (!pp) continue;
          var dd = pixelDist(pClickM[0], pClickM[1], pp[0], pp[1]);
          if (dd < bestD) { bestD = dd; bestI = mi; }
        }
        if (bestI >= 0 && bestD <= DRAG_HANDLE_PX) {
          dragPointIndex = bestI;
          evt.preventDefault();
          evt.stopPropagation();
          return;
        }
      }
      if (!startXY || !endXY) return;
      var pClick = [evt.clientX, evt.clientY];
      var pS = dataToPixelXY(startXY[0], startXY[1]);
      var pE = dataToPixelXY(endXY[0], endXY[1]);
      if (!pS || !pE) return;
      var ds = pixelDist(pClick[0], pClick[1], pS[0], pS[1]);
      var de = pixelDist(pClick[0], pClick[1], pE[0], pE[1]);
      if (ds <= DRAG_HANDLE_PX && ds <= de) {
        dragEndpoint = "start";
        removeTrajectorySampleMarkersOnly();
        evt.preventDefault();
        evt.stopPropagation();
      } else if (de <= DRAG_HANDLE_PX) {
        dragEndpoint = "end";
        removeTrajectorySampleMarkersOnly();
        evt.preventDefault();
        evt.stopPropagation();
      }
    }, true);
    window.addEventListener("mousemove", function(evt) {
      if (dragPointIndex >= 0) {
        var ptd = plotClickToData(evt, true);
        if (!ptd || !editableTrajXY || !editableTrajXY.length) return;
        if (lastTrajectoryVolumeCacheId) setGeneratedVolumesState("", 0);
        var lastIdx = editableTrajXY.length - 1;
        if (trajectoryMode === "direct" && (dragPointIndex === 0 || dragPointIndex === lastIdx)) {
          var fixedFirst = editableTrajXY[0];
          var fixedLast = editableTrajXY[lastIdx];
          if (dragPointIndex === 0) {
            editableTrajXY = makeLinearEditablePath([ptd[0], ptd[1]], fixedLast, editableTrajXY.length) ||
              editableTrajXY;
          } else {
            editableTrajXY = makeLinearEditablePath(fixedFirst, [ptd[0], ptd[1]], editableTrajXY.length) ||
              editableTrajXY;
          }
        } else {
          editableTrajXY[dragPointIndex] = [ptd[0], ptd[1]];
        }
        startXY = editableTrajXY[0].slice();
        endXY = editableTrajXY[lastIdx].slice();
        scheduleEditableTrajectoryRestyle();
        return;
      }
      if (!dragEndpoint) return;
      var pt = plotClickToData(evt, true);
      if (!pt) return;
      if (lastTrajectoryVolumeCacheId) setGeneratedVolumesState("", 0);
      if (dragEndpoint === "start") {
        startXY = [pt[0], pt[1]];
      } else {
        endXY = [pt[0], pt[1]];
      }
      scheduleTrajectoryRestyle();
    });
    window.addEventListener("mouseup", function() {
      if (dragPointIndex >= 0) {
        cancelEditableDragRaf();
        restyleEditableTrajectoryFromState();
        dragPointIndex = -1;
        scheduleCoordsAfterDrag();
        return;
      }
      if (!dragEndpoint) return;
      var hadDrag = dragEndpoint;
      dragEndpoint = null;
      cancelDragRestyleRaf();
      restyleTrajectoryLineOnly();
      if (!anchorIndicesActive && (hadDrag === "start" || hadDrag === "end") && startXY && endXY) {
        scheduleCoordsAfterDrag();
      }
    });
  }

  /* ── Axis/color change ── */
  function triggerRedraw() {
    preparePlotReload();
    syncTraversalModeButtonStates();
    loadPlot();
  }
  [sx, sy].forEach(function(s) {
    s.addEventListener("change", triggerRedraw);
  });
  sc.addEventListener("change", function() {
    loadPlot({ preserveTrajectoryState: true });
  });
  document.querySelectorAll("input[name=\"traj_palette\"]").forEach(function(inp) {
    inp.addEventListener("change", function() {
      loadPlot({ preserveTrajectoryState: true });
    });
  });
  nPointsEl.addEventListener("change", function() {
    syncPointsSelectorRange();
    updateGenerateVolumesButtonLabel();
    removeTrajectorySampleMarkersOnly();
    if (hasAnchorIndices() && currentAnchorTraversalMode() === "graph") return;
    if (hasAnchorIndices()) {
      fetchTrajectoryCoords();
      return;
    }
    editableTrajXY = null;
    if (startXY && endXY) {
      fetchTrajectoryCoords();
    }
  });
  if (maxNeighborsEl) {
    maxNeighborsEl.addEventListener("change", function() {
      if (hasAnchorIndices() && currentAnchorTraversalMode() === "graph") {
        fetchTrajectoryCoords();
      }
    });
  }
  if (avgNeighborsEl) {
    avgNeighborsEl.addEventListener("change", function() {
      if (hasAnchorIndices() && currentAnchorTraversalMode() === "graph") {
        fetchTrajectoryCoords();
      }
    });
  }

  updateGenerateVolumesButtonLabel();
  syncAnchorIndicesUI();
  syncTraversalModeButtonStates();
  syncScatterControlsAlignment();
  loadPlot();
})();
