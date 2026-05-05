(function() {
  var b = window.CRYO_EXPLORER_BOOT;
  var cols = b.cols;
  var covariateDisplayMap = b.covariateDisplayMap;
  var dx = b.defaultX;
  var dy = b.defaultY;
  var initialRows = b.initialRows;
  var totalParticles = b.totalParticles;
  var workdirForSaveHint = b.workdirForSaveHint;
  var showVolumeExplorer = b.showVolumeExplorer;
  var preloadCpus = b.preloadCpus;
  /** Match ``api_scatter`` when ``explorer_scatter=1`` (``--filter-max`` or 200k default). */
  var explorerScatterMaxPoints = b.explorerScatterMaxPoints;
  var explorerScatterCapFromEnv = b.explorerScatterCapFromEnv;

  var preloadRequestedCacheTarget = null;

  var preselectRowsForNextFetch = (initialRows && initialRows.length && initialRows.length <= 5000)
    ? initialRows.slice() : [];
  var saveSelectionRows = [];
  var selectionTooltipText = "";
  var imagesViewEnabled = false;
  var particleSelFs = document.getElementById("particle-sel-fieldset");
  var clearExplorerSelectionBtn = document.getElementById("clear-explorer-selection");
  var selCountRowEl = document.getElementById("sel-count-row");
  var selCountEl = document.getElementById("sel-count");
  var selCountFootnoteEl = document.getElementById("sel-count-footnote");
  var selPieEl = document.getElementById("sel-pie");
  var selPiePctEl = document.getElementById("sel-pie-pct");
  var imageCacheProgressEl = document.getElementById("image-cache-progress");
  var imageCacheProgressBarEl = document.getElementById("image-cache-progress-bar");
  var imageCacheProgressLabelEl = document.getElementById("image-cache-progress-label");
  var montageInner = document.getElementById("montage-inner");
  var montagePreloadOverlayEl = document.getElementById("montage-preload-overlay");
  var montagePreloadOverlayLabel = document.getElementById("montage-preload-overlay-label");
  var btnExpandCache = document.getElementById("btn-expand-cache");
  var btnViewImages = document.getElementById("btn-view-images");
  var btnVolumeGenerate = document.getElementById("btn-volume-generate");
  var btnVolumeAnimate = document.getElementById("btn-volume-animate");
  var volumeStatusEl = document.getElementById("volume-status");
  var volumeStatusProgressEl = document.getElementById("volume-status-progress");
  var gridSizeSelect = document.getElementById("grid-size");
  var imageGridMenuToggle = document.getElementById("image-grid-menu-toggle");
  var imageGridSizeInline = document.getElementById("image-grid-size-inline");
  var imageGridMenuBody = document.getElementById("image-grid-menu-body");
  var imageGridInactiveNote = document.getElementById("image-grid-inactive-note");
  var scatterControlsCard = document.getElementById("scatter-controls-card");
  var scatterControlsSide = scatterControlsCard ? scatterControlsCard.closest(".cryo-dash-side") : null;
  var scatterPlotStack = document.getElementById("scatter-plot-stack");
  var boxSelectTooltipEl = document.getElementById("box-select-tooltip");
  var montageActionsEl = document.getElementById("montage-actions");

  var montageDisplayMode = "images";
  var volumeStaticDataUrls = null;
  var volumeStaticKey = null;
  var volumeGifDataUrls = null;
  var volumeGifsKey = null;
  var EXPLORER_GIF_FRAMES = 40;
  var EXPLORER_CHIMERAX_CPUS = 16;
  var lastMontageRows = [];
  var volumeCacheId = null;
  var volumeAnimateGeneration = 0;

  function imageGridMenuIsOpen() {
    return !!(imageGridMenuToggle && imageGridMenuToggle.getAttribute("aria-expanded") === "true");
  }

  function setImageGridMenuOpen(open) {
    var isOpen = !!open;
    if (imageGridMenuToggle) imageGridMenuToggle.setAttribute("aria-expanded", isOpen ? "true" : "false");
    if (imageGridSizeInline) imageGridSizeInline.hidden = !isOpen;
    if (imageGridMenuBody) imageGridMenuBody.hidden = !isOpen;
    if (!isOpen) {
      updateMontage([]);
    } else if (imagesViewEnabled && preloaded && preloaded.traceMap) {
      showRandomPreloaded();
    }
    syncHighlightTraceAnnotations();
  }

  function syncImageGridInactiveNote() {
    if (!imageGridInactiveNote) return;
    var inactiveNoCache = !preloaded && !imagesViewEnabled;
    imageGridInactiveNote.hidden = !inactiveNoCache;
  }

  function montageRowsOrderKey() {
    return lastMontageRows.length ? lastMontageRows.join(",") : "";
  }

  function montagePreloadOverlayIsShown() {
    return !!(montagePreloadOverlayEl &&
      montagePreloadOverlayEl.classList.contains("cryo-plot-rendering-overlay--show"));
  }

  function setMontagePreloadOverlay(on, label) {
    if (!montagePreloadOverlayEl) return;
    if (label && montagePreloadOverlayLabel) {
      montagePreloadOverlayLabel.textContent = label;
    }
    montagePreloadOverlayEl.classList.toggle("cryo-plot-rendering-overlay--show", !!on);
    montagePreloadOverlayEl.setAttribute("aria-hidden", on ? "false" : "true");
    if (!on) preloadRequestedCacheTarget = null;
  }

  function setVolumeStatus(msg, showProgress) {
    if (volumeStatusEl) volumeStatusEl.textContent = msg || "";
    if (volumeStatusProgressEl) volumeStatusProgressEl.hidden = !showProgress;
  }

  function invalidateVolumeArtifacts() {
    if (!showVolumeExplorer) return;
    volumeStaticDataUrls = null;
    volumeStaticKey = null;
    volumeGifDataUrls = null;
    volumeGifsKey = null;
    volumeCacheId = null;
    volumeAnimateGeneration++;
    montageDisplayMode = "images";
    setVolumeStatus("", false);
    syncVolumeExploreButtons();
  }

  function syncVolumeExploreButtons() {
    if (!showVolumeExplorer) return;
    var hasGrid = lastMontageRows && lastMontageRows.length > 0;
    var key = montageRowsOrderKey();
    var hasStatic = volumeStaticDataUrls && volumeStaticKey === key && volumeStaticDataUrls.length;
    var hasCache = !!(volumeCacheId && volumeStaticKey === key);
    if (btnVolumeGenerate) {
      btnVolumeGenerate.disabled = !plotHasScatterData();
    }
    if (btnVolumeAnimate) {
      btnVolumeAnimate.disabled = !hasGrid || !hasStatic || !hasCache;
      btnVolumeAnimate.title = btnVolumeAnimate.disabled
        ? "Generate volumes for the current image grid before animating."
        : "";
    }
  }

  function showImagesInMontageFromCache() {
    if (!preloaded || !lastMontageRows.length) return;
    for (var i = 0; i < montageCells.length; i++) {
      if (i < lastMontageRows.length) {
        var r = lastMontageRows[i];
        var src = preloaded.rowToSrc.get(r);
        if (src) montageCells[i].img.src = src;
      }
    }
    montageDisplayMode = "images";
    syncVolumeExploreButtons();
  }

  function showVolumesInMontage() {
    if (!lastMontageRows.length) return;
    var key = montageRowsOrderKey();
    if (!(volumeStaticDataUrls && volumeStaticKey === key && volumeStaticDataUrls.length)) return;
    var n = Math.min(montageCells.length, lastMontageRows.length, volumeStaticDataUrls.length);
    for (var i = 0; i < montageCells.length; i++) {
      if (i < n) {
        var g = volumeGifDataUrls && volumeGifDataUrls[i];
        var s = volumeStaticDataUrls[i];
        montageCells[i].img.src = g || s || "";
      } else {
        montageCells[i].img.src = "";
      }
    }
    montageDisplayMode = "volumes";
    syncVolumeExploreButtons();
  }

  function updateParticleSelFieldset() {
    var n = saveSelectionRows.length;
    var nFmt = Number(n).toLocaleString("en-CA");
    var totalFmt = Number(totalParticles).toLocaleString("en-CA");
    var plotCap = scatterSubsetRowCount();
    var plotCapFmt = Number(plotCap).toLocaleString("en-CA");
    var curtailed = totalParticles > plotCap;
    if (selCountEl) {
      selCountEl.textContent = "Selected: " + nFmt + "/" + totalFmt + " particles";
    }
    if (selCountFootnoteEl) {
      if (curtailed) {
        var foot = explorerScatterCapFromEnv
          ? "Scatter only shows " + plotCapFmt + " particles due to `--filter-max`."
          : "Scatter only shows " + plotCapFmt + " particles (default cap).";
        selCountFootnoteEl.textContent = foot;
        selCountFootnoteEl.hidden = false;
      } else {
        selCountFootnoteEl.textContent = "";
        selCountFootnoteEl.hidden = true;
      }
    }
    var frac = 0;
    var pctStr = "0.0%";
    if (totalParticles > 0) {
      frac = Math.min(1, Math.max(0, n / totalParticles));
      pctStr = (100 * frac).toFixed(1) + "%";
    }
    if (selPiePctEl) selPiePctEl.textContent = pctStr;
    if (selPieEl) {
      selPieEl.style.setProperty("--sel-frac", String(frac));
      selPieEl.title = pctStr + " of dataset particles selected";
    }
    if (selCountRowEl) {
      var scatterHint = "";
      if (curtailed) {
        scatterHint = "Lasso and box apply only to plotted points.";
      }
      if (selectionTooltipText && scatterHint) {
        selCountRowEl.title = selectionTooltipText + " — " + scatterHint;
      } else {
        selCountRowEl.title = selectionTooltipText || scatterHint;
      }
      var aria = "Selected " + nFmt + " of " + totalFmt + " particles, " + pctStr;
      if (curtailed) {
        aria += "; scatter " + plotCapFmt + " of " + totalFmt;
      }
      selCountRowEl.setAttribute("aria-label", aria);
    }
    if (particleSelFs) particleSelFs.disabled = (n === 0);
    if (clearExplorerSelectionBtn) clearExplorerSelectionBtn.disabled = (n === 0);
    syncCacheSelectionUncachedButton();
  }

  function countUncachedSelectionRows() {
    if (!saveSelectionRows || !saveSelectionRows.length) return 0;
    var n = 0;
    for (var i = 0; i < saveSelectionRows.length; i++) {
      var r = saveSelectionRows[i];
      if (!preloaded || !preloaded.rows || !preloaded.rows.has(r)) n++;
    }
    return n;
  }

  function syncCacheSelectionUncachedButton() {
    if (!btnCacheSelectionUncached) return;
    var busy = montagePreloadOverlayIsShown();
    var nUnc = countUncachedSelectionRows();
    var can = !!(
      saveSelectionRows
      && saveSelectionRows.length
      && nUnc > 0
      && totalParticles
      && plotHasScatterData()
      && !busy
    );
    btnCacheSelectionUncached.disabled = !can;
    btnCacheSelectionUncached.textContent = "Add " + String(nUnc) + " selection images to cache";
    btnCacheSelectionUncached.setAttribute("aria-label", btnCacheSelectionUncached.textContent);
    btnCacheSelectionUncached.title = !saveSelectionRows || !saveSelectionRows.length
      ? "Select particles on the plot or via color controls first."
      : nUnc === 0
        ? "All selected particles are already in the image cache."
        : "Fetch thumbnails for " + nUnc + " selected particle"
          + (nUnc === 1 ? "" : "s") + " not yet cached.";
  }

  function loadUncachedSelectionIntoCache() {
    if (!saveSelectionRows || !saveSelectionRows.length) return;
    var uncached = [];
    for (var i = 0; i < saveSelectionRows.length; i++) {
      var r = saveSelectionRows[i];
      if (!preloaded || !preloaded.rows || !preloaded.rows.has(r)) uncached.push(r);
    }
    if (!uncached.length) return;
    var have = preloaded ? cachedImageCount() : 0;
    var target = have + uncached.length;
    var overlayLabel = "Caching " + uncached.length + " selected image"
      + (uncached.length === 1 ? "" : "s") + "…";
    fetchPreload(sx.value, sy.value, null, {
      cacheSize: target,
      delta: !!preloaded,
      enableImages: true,
      restrictToScatterPlot: true,
      selectedRows: uncached,
      overlayLabel: overlayLabel
    }).catch(function(e) {
      console.error("cache selection failed", e);
    });
  }

  /** A–Z omitting I, O, U (easier to read vs 1/0 and avoids ambiguous glyphs). */
  var SAFE_LETTERS = (function() {
    var out = [];
    for (var c = 65; c <= 90; c++) {
      if (c !== 73 && c !== 79 && c !== 85) out.push(String.fromCharCode(c));
    }
    return out;
  })();

  /**
   * Generate labels for the current grid.
   * ≤26 images: first labels are single safe letters; any past 23 use pairs (AA, AB, …).
   * >26 images: row letter from safe alphabet × column number (e.g. A1…An, B1…Bn).
   */
  function makeLabels(gridSz) {
    var k = gridSz * gridSz;
    function labelAt(idx) {
      if (idx < SAFE_LETTERS.length) return SAFE_LETTERS[idx];
      var j = idx - SAFE_LETTERS.length;
      return SAFE_LETTERS[Math.floor(j / SAFE_LETTERS.length)]
        + SAFE_LETTERS[j % SAFE_LETTERS.length];
    }
    if (k <= 26) {
      var out = [];
      for (var i = 0; i < k; i++) out.push(labelAt(i));
      return out;
    }
    var out2 = [];
    for (var row = 0; row < gridSz; row++) {
      var letter = SAFE_LETTERS[row];
      for (var col = 0; col < gridSz; col++) {
        out2.push(letter + (col + 1));
      }
    }
    return out2;
  }

  /* RGB stops for montage cell tinting (aligned with Plotly sequential scales). */
  var CONTINUOUS_PALETTE_RGB = {
    Viridis:
      [[68,1,84],[72,39,119],[63,74,138],[49,104,142],[38,130,142],[31,157,138],[53,183,121],[109,205,89],[253,231,37]],
    Plasma:
      [[13,8,135],[76,2,161],[126,3,168],[170,35,149],[204,71,120],[230,108,92],[248,149,64],[253,197,39],[240,249,33]],
    Inferno: [[0,0,4],[41,10,108],[108,2,112],[187,55,84],[249,142,8],[252,208,74],[252,254,164]],
    Magma: [[0,0,4],[59,15,112],[122,4,111],[203,62,65],[252,141,89],[252,208,136],[252,253,191]],
    Cividis: [[0,34,78],[0,57,112],[68,90,129],[115,128,131],[159,161,135],[206,186,122],[253,231,37]],
    Turbo: [[48,18,59],[61,99,221],[26,152,223],[25,189,114],[208,230,28],[250,85,8],[144,12,2]],
    Blues: [[247,251,255],[198,219,239],[107,174,214],[33,113,181],[8,48,107]],
    Greens: [[247,252,245],[199,233,192],[116,196,118],[35,139,69],[0,68,27]],
    Greys: [[255,255,255],[217,217,217],[150,150,150],[82,82,82],[0,0,0]],
    Oranges: [[255,245,235],[253,208,162],[253,141,60],[217,72,1],[127,39,4]],
    Purples: [[252,251,253],[218,218,235],[158,154,200],[106,81,163],[63,0,125]],
    Reds: [[255,245,240],[252,187,161],[251,106,74],[203,24,29],[103,0,13]],
    YlGnBu: [[255,255,217],[199,233,180],[65,182,196],[44,127,184],[8,29,88]],
    YlOrRd: [[255,255,204],[255,237,160],[254,178,76],[240,59,32],[128,0,38]],
    RdBu: [[103,0,31],[214,96,77],[247,247,247],[67,147,195],[5,48,97]],
    Portland: [[12,51,131],[10,136,186],[242,211,56],[242,143,56],[217,30,30]],
    Jet: [[0,0,127],[0,0,255],[0,255,255],[255,255,0],[255,0,0],[127,0,0]],
    Hot: [[0,0,0],[176,0,0],[255,140,0],[255,255,102],[255,255,255]],
    Blackbody: [[0,0,0],[91,12,107],[184,50,137],[242,109,75],[249,248,113]],
    Electric: [[0,0,0],[30,0,109],[141,23,165],[217,79,213],[255,255,255]],
    Rainbow: [[110,64,170],[71,118,230],[26,199,194],[123,217,76],[244,208,63],[242,95,92]],
    Earth: [[0,0,130],[0,181,200],[125,190,66],[200,182,74],[139,69,19]]
  };
  function selectedScatterPalette() {
    var r = document.querySelector("input[name=\"scatter_palette\"]:checked");
    var v = r && r.value ? String(r.value) : "Viridis";
    return CONTINUOUS_PALETTE_RGB[v] ? v : "Viridis";
  }
  function paletteScaleCSS(paletteName, t) {
    t = Math.max(0, Math.min(1, t));
    var pal = CONTINUOUS_PALETTE_RGB[paletteName] || CONTINUOUS_PALETTE_RGB.Viridis;
    var n = pal.length - 1, i = Math.min(Math.floor(t * n), n - 1), f = t * n - i;
    var a = pal[i], b = pal[i + 1];
    return "rgb(" + Math.round(a[0]+(b[0]-a[0])*f) + ","
                  + Math.round(a[1]+(b[1]-a[1])*f) + ","
                  + Math.round(a[2]+(b[2]-a[2])*f) + ")";
  }
  function covariateScaleCSS(t) {
    return paletteScaleCSS(selectedScatterPalette(), t);
  }
  var ACCENT = "#c4703a";

  var gridSize = 3;
  var neighborK = gridSize * gridSize;
  var labels = makeLabels(gridSize);

  /* ---------- montage grid DOM ---------- */
  var gridEl = document.getElementById("montage-grid");
  var montageCells = [];  /* [{img, idx, covar, cell}] */

  function rebuildGrid() {
    neighborK = gridSize * gridSize;
    labels = makeLabels(gridSize);
    gridEl.innerHTML = "";
    gridEl.style.gridTemplateColumns = "repeat(" + gridSize + ", 1fr)";
    montageCells = [];
    var labelSize = gridSize <= 3 ? "0.82rem" : gridSize <= 5 ? "0.65rem" : "0.52rem";
    var metaSize = gridSize <= 3 ? "0.52rem" : gridSize <= 5 ? "0.46rem" : "0.4rem";
    for (var i = 0; i < neighborK; i++) {
      var cell = document.createElement("div");
      cell.className = "cryo-montage-cell";

      var hdr = document.createElement("div");
      hdr.className = "cryo-montage-header";
      var lbl = document.createElement("span");
      lbl.className = "cryo-montage-label";
      lbl.textContent = labels[i];
      lbl.style.fontSize = labelSize;
      var covar = document.createElement("span");
      covar.className = "cryo-montage-covar";
      covar.style.fontSize = metaSize;
      hdr.appendChild(lbl);
      hdr.appendChild(covar);

      var wrap = document.createElement("div");
      wrap.className = "cryo-montage-img-wrap";
      var img = document.createElement("img");
      img.src = ""; img.alt = "";
      wrap.appendChild(img);

      var ftr = document.createElement("div");
      ftr.className = "cryo-montage-footer";
      var idxSpan = document.createElement("span");
      idxSpan.className = "cryo-montage-idx";
      idxSpan.style.fontSize = metaSize;
      ftr.appendChild(idxSpan);

      cell.appendChild(hdr);
      cell.appendChild(wrap);
      cell.appendChild(ftr);
      gridEl.appendChild(cell);
      montageCells.push({img: img, idx: idxSpan, covar: covar, cell: cell});
    }
  }
  rebuildGrid();

  gridSizeSelect.addEventListener("change", function() {
    gridSize = parseInt(this.value) || 3;
    rebuildGrid();
    refillMontageForNewGridSize();
  });

  function fillSelect(sel, includeNone) {
    sel.innerHTML = "";
    if (includeNone) {
      var o = document.createElement("option");
      o.value = "none"; o.textContent = "None";
      sel.appendChild(o);
    }
    cols.forEach(function(c) {
      var o = document.createElement("option");
      o.value = c;
      o.textContent = covariateDisplayMap[c] || c;
      sel.appendChild(o);
    });
  }

  var sx = document.getElementById("sx");
  var sy = document.getElementById("sy");
  var sc = document.getElementById("sc");
  fillSelect(sx, false);
  fillSelect(sy, false);
  fillSelect(sc, true);
  sx.value = cols.indexOf(dx) >= 0 ? dx : cols[0];
  sy.value = cols.indexOf(dy) >= 0 ? dy : cols[Math.min(1, cols.length - 1)];
  var renderedScatterAxes = {x: sx.value, y: sy.value};

  var gd = document.getElementById("scatter");
  var colorHistPanel = document.getElementById("color-hist-panel");
  var colorHistDiv = document.getElementById("color-hist");
  var colorContinuousWrap = document.getElementById("color-continuous-hist-wrap");
  var colorDiscreteWrap = document.getElementById("color-discrete-wrap");
  var colorDiscreteSwitches = document.getElementById("color-discrete-switches");
  var btnCacheSelectionUncached = document.getElementById("btn-cache-selection-uncached");
  var suppressDiscreteSwitchProgrammatic = false;
  var colorThresholdUseMax = document.getElementById("color-threshold-use-max");
  var colorThresholdStatus = document.getElementById("color-threshold-status");
  var scatterPaletteToggle = document.getElementById("scatter-palette-toggle");
  var scatterPaletteOptions = document.getElementById("scatter-palette-options");
  var colorHistEventsWired = false;
  var colorThresholdLevel = null;
  var colorThresholdColorCol = null;
  var colorHistHoverLevel = null;
  var colorHistHoverRaf = null;
  var colorViolinBaseShapes = [];
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

  function equalizeMontageActionButtonWidths() {
    if (!montageActionsEl) return;
    function runGroup(selector) {
      var btns = montageActionsEl.querySelectorAll(selector);
      if (!btns || !btns.length) return;
      for (var i = 0; i < btns.length; i++) {
        btns[i].style.width = "auto";
      }
      var maxW = 0;
      for (var j = 0; j < btns.length; j++) {
        var w = btns[j].getBoundingClientRect().width;
        if (w > maxW) maxW = w;
      }
      if (!maxW) return;
      var px = Math.ceil(maxW) + "px";
      for (var k = 0; k < btns.length; k++) {
        btns[k].style.width = px;
      }
    }
    runGroup(".cryo-explorer-load-images-toolbar .cryo-explorer-action-row__btn .btn");
    runGroup(".cryo-explorer-action-row .cryo-explorer-action-row__btn .btn");
  }

  var plotlyStackEl = gd && gd.parentElement;
  if (plotlyStackEl && typeof ResizeObserver !== "undefined") {
    var plotStackResizeObserver = new ResizeObserver(function() {
      window.requestAnimationFrame(function() {
        syncScatterControlsAlignment();
        equalizeMontageActionButtonWidths();
        if (!gd || !gd.data || !gd.data.length) return;
        try { Plotly.Plots.resize(gd); } catch (e) {}
        if (colorHistPanel && !colorHistPanel.hidden && colorHistDiv) {
          try { Plotly.Plots.resize(colorHistDiv); } catch (e2) {}
        }
      });
    });
    plotStackResizeObserver.observe(plotlyStackEl);
  }
  var overlay = document.getElementById("scatter-rendering-overlay");
  var paletteFieldset = document.getElementById("scatter-palette-radios");
  var plotStatus = document.getElementById("scatter-plot-status");
  var preloadStatus = document.getElementById("preload-status");
  var scatterPlotWatchdog = null;
  var SCATTER_PLOT_TIMEOUT_MS = 90000;
  var scatterLoadGeneration = 0;
  var pendingScatterAfterPlot = null;
  var highlightTraceAdded = false;
  var suppressPlotGridHighlights = false;
  var suppressSelectionEvents = false;

  var selFileBrowserPanel = document.getElementById("sel-file-browser-panel");
  var selFbList = document.getElementById("sel-fb-list");
  var selFbPath = document.getElementById("sel-fb-path");
  var selFbUp = document.getElementById("sel-fb-up");
  var selFbCancel = document.getElementById("sel-fb-cancel");
  var selFbSaveName = document.getElementById("sel-fb-save-name");
  var selFbSaveBtn = document.getElementById("sel-fb-save-btn");
  var selFbCurrentDir = null;

  function plotHasScatterData() {
    return !!(gd && gd.data && gd.data[0] && gd.data[0].customdata && gd.data[0].customdata.length);
  }

  function shuffleInPlace(arr) {
    for (var i = arr.length - 1; i > 0; i--) {
      var j = Math.floor(Math.random() * (i + 1));
      var t = arr[i];
      arr[i] = arr[j];
      arr[j] = t;
    }
    return arr;
  }

  function allPlotRowsFromTrace() {
    var cd = gd.data[0].customdata;
    var out = [];
    for (var i = 0; i < cd.length; i++) {
      out.push(cd[i][1]);
    }
    return out;
  }

  function scatterPlottedRowBudget() {
    var n = Math.max(0, Number(totalParticles) || 0);
    return Math.min(n, explorerScatterMaxPoints);
  }

  /** Rows in the scatter subsample (matches server ``restrict_to_scatter_plot`` pool). */
  function scatterSubsetRowCount() {
    if (plotHasScatterData()) return gd.data[0].customdata.length;
    return scatterPlottedRowBudget();
  }

  /** Largest power of 10 ≤ 5% of plotted count (matches ``explorer_cache_size_power10_step``). */
  function cacheSizePower10StepFromCap(cap) {
    cap = Math.max(0, Math.floor(Number(cap) || 0));
    if (cap === 0) return 1;
    var x = 0.05 * cap;
    if (x < 1) return 1;
    var k = Math.floor(Math.log10(x));
    var step = Math.pow(10, k);
    return Math.max(1, Math.floor(step));
  }

  /**
   * ``<input type="number">`` only steps correctly when (value − min) is a multiple of step.
   * We keep min/step/value on that lattice so browser spinners use the power-of-ten step.
   */
  function montageCacheSizeStepperGrid(cap) {
    cap = Math.max(0, Math.floor(Number(cap) || 0));
    if (cap === 0) return { minG: 0, step: 1 };
    var step = cacheSizePower10StepFromCap(cap);
    if (!preloaded) {
      var stepUse = Math.min(step, cap);
      return { minG: 0, step: Math.max(1, stepUse) };
    }
    var have = cachedImageCount();
    var floor = Math.min(have + 1, cap);
    if (have >= cap) return { minG: cap, step: 1 };
    var minAligned = Math.ceil(floor / step) * step;
    if (minAligned > cap) return { minG: floor, step: 1 };
    return { minG: minAligned, step: step };
  }

  function snapMontageCacheSizeToGridNearest(val, cap, grid) {
    var minG = grid.minG;
    var step = grid.step;
    cap = Math.max(0, Math.floor(Number(cap) || 0));
    if (cap === 0) return 0;
    if (!isFinite(val)) val = minG;
    val = Math.floor(Number(val));
    val = Math.min(cap, Math.max(minG, val));
    if (step <= 1) return val;
    var k = Math.round((val - minG) / step);
    if (!isFinite(k)) k = 0;
    var out = minG + k * step;
    if (out > cap) out = minG + Math.floor((cap - minG) / step) * step;
    if (out < minG) out = minG;
    return out;
  }

  /** Default first-build target: step above, capped at plotted count (matches server ``preload_image_limit``). */
  function firstBuildCacheSizeDefault() {
    var cap = scatterSubsetRowCount();
    if (cap === 0) return 1;
    var step = cacheSizePower10StepFromCap(cap);
    return Math.min(cap, Math.max(1, step));
  }

  /** Suggested new cache total after a cache exists: min(2 × current, cap), always &gt; current when possible. */
  function defaultMontageNewCacheSizeSuggestion() {
    if (!preloaded) return firstBuildCacheSizeDefault();
    var have = cachedImageCount();
    var cap = scatterSubsetRowCount();
    if (have >= cap) return cap;
    var twiceCapped = Math.min(2 * have, cap);
    return Math.max(twiceCapped, have + 1);
  }

  function syncMontageCacheSizeLabelEl() {
    var sp = document.getElementById("montage-cache-size-label-text");
    var inp = document.getElementById("montage-cache-size-input");
    if (sp) sp.textContent = preloaded ? "New cache size" : "Cache size";
    if (inp) {
      inp.title = preloaded
        ? "Total cached images to grow toward (must be greater than the current cache; capped at plotted points)"
        : "Number of thumbnails to fetch when you click Build cache (capped at plotted scatter points)";
    }
  }

  /** Keep min/max/value in range; during preload overlay show requested target (disabled). */
  function syncMontageNewCacheSizeField() {
    var el = document.getElementById("montage-cache-size-input");
    if (!el) return;
    var cap = scatterSubsetRowCount();
    var grid = montageCacheSizeStepperGrid(cap);
    el.setAttribute("min", String(grid.minG));
    el.setAttribute("step", String(grid.step));
    el.setAttribute("max", String(cap));
    syncMontageCacheSizeLabelEl();
    if (montagePreloadOverlayIsShown()) {
      el.disabled = true;
      if (preloadRequestedCacheTarget != null) {
        var tgt = snapMontageCacheSizeToGridNearest(preloadRequestedCacheTarget, cap, grid);
        el.value = String(tgt);
      }
      return;
    }
    var have = cachedImageCount();
    if (preloaded && have >= cap) {
      el.disabled = true;
      el.setAttribute("min", String(cap));
      el.setAttribute("step", "1");
      el.value = String(cap);
      return;
    }
    el.disabled = false;
    if (!preloaded) {
      var d0 = firstBuildCacheSizeDefault();
      var raw0 = parseInt(String(el.value).replace(/\s/g, "").replace(/,/g, ""), 10);
      if (!isFinite(raw0) || raw0 < 0 || raw0 > cap) raw0 = d0;
      var snapped0 = snapMontageCacheSizeToGridNearest(raw0, cap, grid);
      if (snapped0 < 1) snapped0 = snapMontageCacheSizeToGridNearest(d0, cap, grid);
      if (snapped0 < 1) snapped0 = 1;
      el.value = String(Math.min(cap, snapped0));
      return;
    }
    var suggestion = defaultMontageNewCacheSizeSuggestion();
    var raw = parseInt(String(el.value).replace(/\s/g, "").replace(/,/g, ""), 10);
    var floor = Math.min(have + 1, cap);
    if (!isFinite(raw) || raw <= have || raw > cap) raw = suggestion;
    var snapped = snapMontageCacheSizeToGridNearest(raw, cap, grid);
    if (snapped < floor) snapped = snapMontageCacheSizeToGridNearest(suggestion, cap, grid);
    el.value = String(Math.min(cap, Math.max(floor, snapped)));
  }

  /** Target cache size from field (first build: [1, cap]; expand: (have, cap]). */
  function readMontageCacheSizeInput() {
    var cap = scatterSubsetRowCount();
    var el = document.getElementById("montage-cache-size-input");
    var raw = el ? parseInt(String(el.value).replace(/\s/g, "").replace(/,/g, ""), 10) : NaN;
    var grid = montageCacheSizeStepperGrid(cap);
    if (!preloaded) {
      if (!isFinite(raw)) raw = firstBuildCacheSizeDefault();
      var v = snapMontageCacheSizeToGridNearest(raw, cap, grid);
      if (v < 1) v = snapMontageCacheSizeToGridNearest(firstBuildCacheSizeDefault(), cap, grid);
      return Math.min(cap, Math.max(1, v));
    }
    var have = cachedImageCount();
    var floor = have >= cap ? cap : have + 1;
    if (!isFinite(raw)) raw = defaultMontageNewCacheSizeSuggestion();
    var v2 = snapMontageCacheSizeToGridNearest(raw, cap, grid);
    return Math.min(cap, Math.max(floor, v2));
  }

  function rowToTraceIndex(row) {
    var cd = gd.data[0].customdata;
    for (var i = 0; i < cd.length; i++) {
      if (cd[i][1] === row) return i;
    }
    return -1;
  }

  function montageItemForRow(r) {
    var ti = rowToTraceIndex(r);
    if (ti < 0) return null;
    var src = "";
    if (preloaded && preloaded.rowToSrc) {
      src = preloaded.rowToSrc.get(r) || "";
    }
    return { ti: ti, row: r, src: src };
  }

  /** Random grid: prefer lasso/box selection, pad from full plot if needed. */
  function pickRandomMontageItemsFilled(want) {
    if (!plotHasScatterData()) return [];
    var primary = [];
    if (saveSelectionRows && saveSelectionRows.length) {
      primary = saveSelectionRows.slice();
      shuffleInPlace(primary);
    } else {
      primary = allPlotRowsFromTrace();
      shuffleInPlace(primary);
    }
    var seen = new Set();
    var out = [];
    for (var a = 0; a < primary.length && out.length < want; a++) {
      var r = primary[a];
      if (seen.has(r)) continue;
      var it = montageItemForRow(r);
      if (!it) continue;
      seen.add(r);
      out.push(it);
    }
    if (out.length < want) {
      var fallback = allPlotRowsFromTrace();
      shuffleInPlace(fallback);
      for (var b = 0; b < fallback.length && out.length < want; b++) {
        var r2 = fallback[b];
        if (seen.has(r2)) continue;
        var it2 = montageItemForRow(r2);
        if (!it2) continue;
        seen.add(r2);
        out.push(it2);
      }
    }
    return out;
  }

  /* ---------- preloaded image cache ---------- */
  var preloaded = null;

  function cachedImageCount() {
    return preloaded && preloaded.rows ? preloaded.rows.size : 0;
  }

  /** Cached rows that also appear on the current scatter trace (subsample may omit most dataset rows). */
  function traceMapCachedCount() {
    return preloaded && preloaded.traceMap ? preloaded.traceMap.length : 0;
  }

  function formatCachedSelectionStatusCore(activeInSelection) {
    var nPlot = traceMapCachedCount();
    var nTot = cachedImageCount();
    if (!nTot) return "";
    if (nPlot === nTot) {
      return activeInSelection + " of " + nTot + " cached images in selection";
    }
    return activeInSelection + " of " + nPlot + " on-scatter cached in selection ("
      + nTot + " thumbnails loaded)";
  }

  function setPreloadStatusCachedTotalOnly() {
    if (!preloadStatus || !preloaded) return;
    var n = cachedImageCount();
    preloadStatus.textContent = n ? n + " images cached." : "";
  }

  function allImagesCached() {
    return cachedImageCount() >= scatterSubsetRowCount();
  }

  function syncExpandCacheButton(loading) {
    if (!btnExpandCache) return;
    syncMontageCacheSizeLabelEl();
    if (loading) {
      btnExpandCache.disabled = true;
      return;
    }
    if (!preloaded) {
      btnExpandCache.textContent = "Build cache";
      var capB = scatterSubsetRowCount();
      var canBuild = !!(totalParticles && plotHasScatterData() && capB > 0);
      btnExpandCache.disabled = !canBuild;
      btnExpandCache.title = !totalParticles
        ? "No particles in this experiment."
        : !plotHasScatterData()
          ? "Wait for the scatter plot to finish loading."
          : "Fetch up to " + Number(capB).toLocaleString("en-CA")
            + " thumbnails (see cache size).";
      return;
    }
    btnExpandCache.textContent = "Expand cache";
    var have = cachedImageCount();
    var want = readMontageCacheSizeInput();
    var cap = scatterSubsetRowCount();
    var canExpand = !!(want > have && have < cap);
    btnExpandCache.disabled = !canExpand;
    if (have >= cap) {
      btnExpandCache.title = "All plotted particle images are already cached.";
    } else if (want <= have) {
      btnExpandCache.title = "Set new cache size above " + Number(have).toLocaleString("en-CA") + " to expand.";
    } else {
      btnExpandCache.title = "Fetch thumbnails up to the new cache size (now " + Number(have).toLocaleString("en-CA") +
        " cached).";
    }
  }

  function syncImageCacheButton(loading) {
    btnViewImages.textContent = "Load all images";
    btnViewImages.disabled = !!loading || !preloaded || allImagesCached();
    btnViewImages.title = allImagesCached()
      ? "All particle images on the current scatter are loaded."
      : (btnViewImages.disabled ? "Build the image cache first." : "");
    syncMontageNewCacheSizeField();
    syncExpandCacheButton(loading);
    syncCacheSelectionUncachedButton();
    syncImageGridInactiveNote();
  }

  function setImageCacheProgress(on, current, total) {
    if (!imageCacheProgressEl) return;
    if (!on) {
      imageCacheProgressEl.hidden = true;
      imageCacheProgressEl.setAttribute("aria-valuenow", "0");
      if (imageCacheProgressBarEl) imageCacheProgressBarEl.style.width = "0%";
      if (imageCacheProgressLabelEl) imageCacheProgressLabelEl.textContent = "0%";
      return;
    }
    total = Math.max(1, Number(total) || 1);
    current = Math.max(0, Math.min(Number(current) || 0, total));
    var pct = Math.round((current / total) * 100);
    imageCacheProgressEl.hidden = false;
    imageCacheProgressEl.setAttribute("aria-valuenow", String(pct));
    if (imageCacheProgressBarEl) imageCacheProgressBarEl.style.width = pct + "%";
    if (imageCacheProgressLabelEl) {
      imageCacheProgressLabelEl.textContent = pct + "% loaded ("
        + Number(current).toLocaleString("en-CA") + " / "
        + Number(total).toLocaleString("en-CA") + ")";
    }
  }

  function imageCacheChunkSize(total) {
    total = Math.max(0, Number(total) || 0);
    if (!total) return 0;
    var cpus = Math.max(1, Number(preloadCpus) || 1);
    var parallelThreshold = Math.max(128, cpus * 32);
    var minChunk = cpus > 1 ? parallelThreshold : 250;
    var targetChunk = Math.ceil(total / 200);
    return Math.min(5000, Math.max(minChunk, targetChunk));
  }

  function fetchPreload(xcol, ycol, onDone, opts) {
    opts = opts || {};
    var cacheSize = opts.cacheSize || scatterSubsetRowCount();
    preloadRequestedCacheTarget = cacheSize;
    var plottedCap = scatterSubsetRowCount();
    var overlayLabel = opts.overlayLabel
      ? opts.overlayLabel
      : (opts.fullLoad
        ? "Loading all " + Number(plottedCap).toLocaleString("en-CA") + " plotted particle images..."
        : "Loading image cache...");
    setMontagePreloadOverlay(true, overlayLabel);
    syncImageCacheButton(true);
    var body = { x: xcol, y: ycol, cache_size: cacheSize };
    if (opts.delta) {
      body.response_mode = "delta";
    }
    if (opts.restrictToScatterPlot) {
      body.restrict_to_scatter_plot = true;
      body.scatter_max_points = opts.scatterMaxPoints != null
        ? opts.scatterMaxPoints
        : explorerScatterMaxPoints;
    }
    if (opts.initialRows && opts.initialRows.length) {
      body.initial_rows = opts.initialRows.slice();
    }
    if (opts.selectedRows && opts.selectedRows.length) {
      body.selected_rows = opts.selectedRows.map(Number).filter(function(x) {
        return isFinite(x);
      });
    }
    return fetch(b.urls.apiPreloadImages, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    })
      .then(function(r) {
        return r.json().then(function(j) { return { ok: r.ok, j: j }; });
      })
      .then(function(res) {
        if (!res.ok) {
          var msg = (res.j && res.j.error) ? res.j.error : "Image cache request failed.";
          if (preloadStatus) preloadStatus.textContent = msg;
          setMontagePreloadOverlay(false);
          if (onDone) onDone(new Error(msg));
          throw new Error(msg);
        }
        var data = res.j;
        if (!data.rows || !data.images || (!data.rows.length && !opts.delta)) {
          if (preloadStatus) preloadStatus.textContent = "No images available.";
          setMontagePreloadOverlay(false);
          if (onDone) onDone(new Error("no images"));
          throw new Error("no images");
        }
        var rowToSrc = (opts.delta && preloaded && preloaded.rowToSrc)
          ? new Map(preloaded.rowToSrc)
          : new Map();
        var rowsSet = (opts.delta && preloaded && preloaded.rows)
          ? new Set(preloaded.rows)
          : new Set();
        data.rows.forEach(function(r, i) {
          rowToSrc.set(r, "data:image/jpeg;base64," + data.images[i]);
          rowsSet.add(r);
        });
        preloaded = {rows: rowsSet, rowToSrc: rowToSrc, traceMap: null};
        if (opts.enableImages) {
          imagesViewEnabled = true;
          montageInner.classList.remove("cryo-dash-montage-inner--disabled");
          gridSizeSelect.disabled = false;
          setImageGridMenuOpen(true);
        }
        buildTraceMap();
        /* Lasso/box before image display only set saveSelectionRows — rebuild activePool for montage + clicks. */
        if (saveSelectionRows && saveSelectionRows.length) {
          var sset = new Set(saveSelectionRows);
          var ap = preloaded.traceMap.filter(function(item) {
            return sset.has(item.row);
          });
          if (ap.length) {
            activePool = ap;
            lassoSelectedRowSet = new Set();
            for (var li = 0; li < ap.length; li++) {
              lassoSelectedRowSet.add(ap[li].row);
            }
          }
        }
        var totalNow = preloaded.rows.size;
        var msg;
        if (opts.delta && data.rows && data.rows.length && data.batch_elapsed != null) {
          var nNew = data.rows.length;
          msg = nNew + (nNew === 1 ? " new image cached in " : " new images cached in ")
            + data.batch_elapsed + "s — " + totalNow + " total in cache";
        } else {
          var cachedTotal = data.total_cached != null ? data.total_cached : totalNow;
          msg = cachedTotal + " images cached";
          if (data.elapsed != null) msg += " in " + data.elapsed + "s";
          msg += " — " + totalNow + " total in cache";
        }
        if (preloadStatus && !opts.suppressStatus) {
          if (activePool && activePool.length) {
            preloadStatus.textContent = formatCachedSelectionStatusCore(activePool.length)
              + " — " + msg + ".";
          } else {
            preloadStatus.textContent = msg + ".";
          }
        }
        if (!opts.keepOverlayAfterSuccess) {
          setMontagePreloadOverlay(false);
        }
        if (!opts.suppressMontageUpdate) {
          showRandomPreloaded();
        }
        syncImageCacheButton(!!opts.keepOverlayAfterSuccess);
        if (onDone) onDone(null);
        return data;
      })
      .catch(function(e) {
        console.error("preload failed", e);
        if (preloadStatus) preloadStatus.textContent = "Image cache failed.";
        setMontagePreloadOverlay(false);
        syncImageCacheButton(false);
        if (onDone) onDone(e);
        if (onDone) return null;
        throw e;
      });
  }

  function expandImageCacheToTarget() {
    if (!preloaded) return;
    var target = readMontageCacheSizeInput();
    var before = cachedImageCount();
    if (target <= before) return;
    var label = "Expanding cache to " + Number(target).toLocaleString("en-CA") + " images…";
    fetchPreload(sx.value, sy.value, null, {
      cacheSize: target,
      delta: true,
      enableImages: true,
      restrictToScatterPlot: true,
      overlayLabel: label
    });
  }

  function buildInitialImageCache() {
    if (preloaded || !totalParticles || !plotHasScatterData()) return;
    if (scatterSubsetRowCount() < 1) return;
    var sz = readMontageCacheSizeInput();
    fetchPreload(sx.value, sy.value, null, {
      cacheSize: sz,
      enableImages: true,
      restrictToScatterPlot: true,
      overlayLabel: "Building cache (" + Number(sz).toLocaleString("en-CA") + " images)…"
    }).catch(function(e) {
      console.error("build cache failed", e);
    });
  }

  function initializeExplorerView() {
    loadPlot(true);
  }

  function loadAllImagesInChunks() {
    var total = scatterSubsetRowCount();
    setMontagePreloadOverlay(
      true,
      "Loading all " + Number(total).toLocaleString("en-CA") + " plotted particle images..."
    );
    if (preloadStatus) preloadStatus.textContent = "";
    setImageCacheProgress(true, cachedImageCount(), total);
    suppressPlotGridHighlights = true;
    syncHighlightTraceAnnotations();
    syncImageCacheButton(true);

    function nextChunk() {
      var before = cachedImageCount();
      if (before >= total) {
        setMontagePreloadOverlay(false);
        setImageCacheProgress(true, total, total);
        suppressPlotGridHighlights = false;
        syncHighlightTraceAnnotations();
        syncImageCacheButton(false);
        syncVolumeExploreButtons();
        return Promise.resolve(null);
      }
      var step = imageCacheChunkSize(total);
      var userTarget = readMontageCacheSizeInput();
      var roomTowardUser = userTarget - before;
      if (roomTowardUser > 0) {
        step = Math.max(1, Math.min(step, roomTowardUser));
      } else {
        step = Math.max(1, Math.min(step, total - before));
      }
      var target = Math.min(total, before + step);
      return fetchPreload(sx.value, sy.value, null, {
        cacheSize: target,
        enableImages: true,
        delta: true,
        fullLoad: true,
        keepOverlayAfterSuccess: true,
        suppressStatus: true,
        suppressMontageUpdate: true,
        restrictToScatterPlot: true
      }).then(function() {
        if (cachedImageCount() <= before) {
          throw new Error("Image cache did not advance.");
        }
        setImageCacheProgress(true, cachedImageCount(), total);
        return nextChunk();
      });
    }

    return nextChunk()
      .catch(function(e) {
        console.error("load all images failed", e);
        setImageCacheProgress(false);
        suppressPlotGridHighlights = false;
        syncHighlightTraceAnnotations();
        if (preloadStatus) preloadStatus.textContent = e.message || "Image cache failed.";
        setMontagePreloadOverlay(false);
      })
      .then(function() {
        syncImageCacheButton(false);
        syncVolumeExploreButtons();
      });
  }

  if (btnExpandCache) {
    btnExpandCache.addEventListener("click", function() {
      if (!preloaded) {
        buildInitialImageCache();
      } else {
        expandImageCacheToTarget();
      }
    });
  }
  if (btnCacheSelectionUncached) {
    btnCacheSelectionUncached.addEventListener("click", function() {
      loadUncachedSelectionIntoCache();
    });
  }

  var montageCacheSizeInput = document.getElementById("montage-cache-size-input");
  if (montageCacheSizeInput) {
    var montageCacheOverlayBusy = function() {
      return !!(montagePreloadOverlayEl &&
        montagePreloadOverlayEl.classList.contains("cryo-plot-rendering-overlay--show"));
    };
    montageCacheSizeInput.addEventListener("input", function() {
      syncExpandCacheButton(montageCacheOverlayBusy());
    });
    montageCacheSizeInput.addEventListener("change", function() {
      syncMontageNewCacheSizeField();
      syncExpandCacheButton(montageCacheOverlayBusy());
    });
    montageCacheSizeInput.addEventListener("blur", function() {
      syncMontageNewCacheSizeField();
      syncExpandCacheButton(montageCacheOverlayBusy());
    });
  }

  btnViewImages.addEventListener("click", function() {
    if (imagesViewEnabled && preloaded && montageDisplayMode === "volumes" && allImagesCached()) {
      showImagesInMontageFromCache();
      return;
    }
    if (imagesViewEnabled && preloaded && allImagesCached()) {
      montageInner.classList.remove("cryo-dash-montage-inner--disabled");
      gridSizeSelect.disabled = false;
      setImageGridMenuOpen(true);
      buildTraceMap();
      if (montageDisplayMode === "images") {
        restoreMontageAfterScatterReload();
      }
      if (preloadStatus && preloaded.rows) {
        preloadStatus.textContent = preloaded.rows.size + " images in cache (no reload).";
      }
      syncVolumeExploreButtons();
      return;
    }
    btnViewImages.disabled = true;
    loadAllImagesInChunks();
  });

  function postVolumeStatic() {
    return fetch(b.urls.apiExplorerVolumeMedia, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        rows: lastMontageRows.slice(),
        mode: "static",
        chimerax_cpus: EXPLORER_CHIMERAX_CPUS,
      }),
    })
      .then(function(r) { return r.json().then(function(j) { return { ok: r.ok, j: j }; }); })
      .then(function(res) {
        if (!res.ok) {
          if (res.j.need_chimerax) {
            window.alert(res.j.error || "Set CHIMERAX_PATH to your ChimeraX executable and try again.");
            setVolumeStatus("", false);
          } else {
            setVolumeStatus(res.j.error || "Volume request failed.", false);
          }
          return null;
        }
        return res.j;
      });
  }

  function postVolumeAnimateCell(cacheId, cellIndex) {
    return fetch(b.urls.apiExplorerVolumeMedia, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        mode: "animate",
        volume_cache_id: cacheId,
        cell_index: cellIndex,
        rows: lastMontageRows.slice(),
        gif_frames: EXPLORER_GIF_FRAMES,
        chimerax_cpus: EXPLORER_CHIMERAX_CPUS,
      }),
    })
      .then(function(r) { return r.json().then(function(j) { return { ok: r.ok, j: j }; }); })
      .then(function(res) {
        if (!res.ok) {
          if (res.j.need_chimerax) {
            window.alert(res.j.error || "Set CHIMERAX_PATH to your ChimeraX executable and try again.");
            setVolumeStatus("", false);
          } else {
            setVolumeStatus(res.j.error || "Animate request failed.", false);
          }
          return null;
        }
        return res.j;
      });
  }

  function runVolumeAnimateSequential() {
    var key = montageRowsOrderKey();
    if (!volumeCacheId || !lastMontageRows.length || volumeStaticKey !== key) return;
    var n = lastMontageRows.length;
    volumeAnimateGeneration++;
    var myGen = volumeAnimateGeneration;
    btnVolumeAnimate.disabled = true;
    if (btnVolumeGenerate) btnVolumeGenerate.disabled = true;
    volumeGifDataUrls = new Array(n);
    volumeGifsKey = null;
    montageDisplayMode = "volumes";
    showVolumesInMontage();
    var i = 0;
    function step() {
      if (myGen !== volumeAnimateGeneration) return;
      if (i >= n) {
        volumeGifsKey = key;
        setVolumeStatus("", false);
        syncVolumeExploreButtons();
        return;
      }
      setVolumeStatus("GIF " + (i + 1) + " / " + n + "…", false);
      postVolumeAnimateCell(volumeCacheId, i)
        .then(function(j) {
          if (myGen !== volumeAnimateGeneration) return;
          if (!j || !j.image) {
            syncVolumeExploreButtons();
            return;
          }
          var url = "data:image/gif;base64," + j.image;
          volumeGifDataUrls[i] = url;
          if (i < montageCells.length) montageCells[i].img.src = url;
          i++;
          step();
        })
        .catch(function() {
          if (myGen === volumeAnimateGeneration) {
            setVolumeStatus("Request failed.", false);
          }
          syncVolumeExploreButtons();
        });
    }
    step();
  }

  if (showVolumeExplorer && btnVolumeGenerate) {
    btnVolumeGenerate.addEventListener("click", function() {
      if (!plotHasScatterData()) return;
      if (!lastMontageRows.length) {
        ensureHighlightTrace();
        var picked = pickRandomMontageItemsFilled(neighborK);
        if (!picked.length) return;
        montageInner.classList.remove("cryo-dash-montage-inner--disabled");
        gridSizeSelect.disabled = false;
        updateMontage(picked);
      }
      if (!lastMontageRows.length) return;
      var key = montageRowsOrderKey();
      btnVolumeGenerate.disabled = true;
      if (btnVolumeAnimate) btnVolumeAnimate.disabled = true;
      setVolumeStatus("Decoding volumes and rendering (ChimeraX)…", true);
      postVolumeStatic()
        .then(function(j) {
          if (!j || !j.images || !j.volume_cache_id) {
            syncVolumeExploreButtons();
            return;
          }
          volumeStaticDataUrls = j.images.map(function(b) { return "data:image/png;base64," + b; });
          volumeStaticKey = key;
          volumeGifDataUrls = null;
          volumeGifsKey = null;
          volumeCacheId = j.volume_cache_id;
          showVolumesInMontage();
          setVolumeStatus("", false);
        })
        .catch(function() {
          setVolumeStatus("Request failed.", false);
        })
        .then(function() { syncVolumeExploreButtons(); });
    });
  }

  if (showVolumeExplorer && btnVolumeAnimate) {
    btnVolumeAnimate.addEventListener("click", function() {
      runVolumeAnimateSequential();
    });
  }

  function scatterTraceColorMode() {
    var meta = gd.layout && gd.layout.meta;
    var m = meta && meta.cdrgn_color_mode;
    if (m === "discrete" || m === "continuous") return m;
    var mc = gd.data && gd.data[0] && gd.data[0].marker && gd.data[0].marker.color;
    if (Array.isArray(mc) && mc.length && typeof mc[0] === "string" && /^#/.test(String(mc[0]))) {
      return "discrete";
    }
    return "continuous";
  }

  function formatColorThresholdValue(v) {
    if (typeof v !== "number" || !isFinite(v)) return "";
    var av = Math.abs(v);
    if (av !== 0 && (av < 0.001 || av >= 10000)) return v.toExponential(3);
    return String(Math.round(v * 1000) / 1000);
  }

  function escapeHtmlText(s) {
    return String(s).replace(/[&<>"']/g, function(ch) {
      return {"&": "&amp;", "<": "&lt;", ">": "&gt;", "\"": "&quot;", "'": "&#39;"}[ch];
    });
  }

  function traceColorNumericValue(ti) {
    if (!gd || !gd.data || !gd.data[0]) return null;
    var trace = gd.data[0];
    var markerColors = trace.marker ? trace.marker.color : null;
    if (Array.isArray(markerColors)) {
      var mv = Number(markerColors[ti]);
      if (isFinite(mv)) return mv;
    }
    var cd = trace.customdata;
    if (cd && cd[ti] && cd[ti].length > 2) {
      var dv = Number(cd[ti][2]);
      if (isFinite(dv)) return dv;
    }
    return null;
  }

  function currentColorValueItems() {
    if (!gd || !gd.data || !gd.data[0] || !gd.data[0].customdata) return [];
    if (!sc.value || sc.value === "none") return [];
    var cd = gd.data[0].customdata;
    var out = [];
    for (var i = 0; i < cd.length; i++) {
      var row = cd[i] && cd[i].length > 1 ? cd[i][1] : null;
      var val = traceColorNumericValue(i);
      if (row != null && val != null) out.push({ti: i, row: row, value: val});
    }
    return out;
  }

  function syncParticleExplorerRowExpandForColorHist() {
    window.requestAnimationFrame(function() {
      syncScatterControlsAlignment();
      try {
        if (gd && gd.data && gd.data.length) Plotly.Plots.resize(gd);
      } catch (e) {}
      if (colorHistPanel && !colorHistPanel.hidden && colorHistDiv) {
        try {
          Plotly.Plots.resize(colorHistDiv);
        } catch (e2) {}
      }
      equalizeMontageActionButtonWidths();
    });
  }

  function hideColorHistogram() {
    if (colorHistPanel) {
      colorHistPanel.hidden = true;
      colorHistPanel.setAttribute("aria-hidden", "true");
    }
    colorHistHoverLevel = null;
    colorViolinBaseShapes = [];
    if (colorThresholdStatus) colorThresholdStatus.textContent = "";
    if (colorContinuousWrap) colorContinuousWrap.hidden = false;
    if (colorDiscreteWrap) {
      colorDiscreteWrap.hidden = true;
      colorDiscreteWrap.classList.remove("cryo-color-discrete-wrap--show");
    }
    if (colorDiscreteSwitches) colorDiscreteSwitches.innerHTML = "";
    syncParticleExplorerRowExpandForColorHist();
  }

  function discreteCategoryKeyFromCustomRow(cdrow) {
    if (!cdrow || cdrow.length < 3) return "__na__";
    var v = cdrow[2];
    if (v == null || v === "") return "__na__";
    if (typeof v === "number" && !isFinite(v)) return "__na__";
    return String(v);
  }

  function discreteCategoryLabelForKey(key) {
    if (key === "__na__") return "(missing)";
    var rawColorCol = sc && sc.value ? String(sc.value) : "";
    var dispColorCol = covariateDisplayMap && sc && sc.value ? String(covariateDisplayMap[sc.value] || "") : "";
    var isKmeans = /k[\s_-]*means/i.test(rawColorCol)
      || /k[\s_-]*means/i.test(dispColorCol)
      || /^labels$/i.test(rawColorCol);
    if (isKmeans) return "kmeans=" + String(key);
    return key;
  }

  function measureDiscreteSwitchColumnWidthPx(labelItems) {
    if (!labelItems || !labelItems.length) return 0;
    var probe = document.createElement("label");
    probe.className = "cryo-color-discrete-switch";
    probe.style.position = "absolute";
    probe.style.visibility = "hidden";
    probe.style.left = "-9999px";
    probe.style.top = "-9999px";
    probe.style.width = "auto";
    var inp = document.createElement("input");
    inp.type = "checkbox";
    var textWrap = document.createElement("span");
    textWrap.className = "cryo-color-discrete-switch-text";
    var labelSpan = document.createElement("span");
    labelSpan.className = "cryo-color-discrete-switch-label";
    var countSpan = document.createElement("span");
    countSpan.className = "cryo-color-discrete-switch-count";
    textWrap.appendChild(labelSpan);
    textWrap.appendChild(countSpan);
    probe.appendChild(inp);
    probe.appendChild(textWrap);
    document.body.appendChild(probe);
    var maxW = 0;
    for (var i = 0; i < labelItems.length; i++) {
      var item = labelItems[i] || {};
      labelSpan.textContent = item.label || "";
      countSpan.textContent = item.count || "";
      var w = Math.ceil(probe.getBoundingClientRect().width);
      if (w > maxW) maxW = w;
    }
    probe.remove();
    return maxW;
  }

  function sortDiscreteCategoryKeys(keys) {
    var rest = keys.filter(function(k) { return k !== "__na__"; });
    var allNum = rest.length > 0 && rest.every(function(k) {
      var n = Number(k);
      return isFinite(n) && String(n) === String(k).trim();
    });
    if (allNum) {
      rest.sort(function(a, b) { return Number(a) - Number(b); });
    } else {
      rest.sort(function(a, b) {
        return String(a).localeCompare(String(b), undefined, {numeric: true});
      });
    }
    if (keys.indexOf("__na__") >= 0) rest.push("__na__");
    return rest;
  }

  function resetDiscreteCategorySwitches() {
    if (!colorDiscreteSwitches) return;
    suppressDiscreteSwitchProgrammatic = true;
    try {
      var inputs = colorDiscreteSwitches.querySelectorAll("input[type=\"checkbox\"]");
      for (var i = 0; i < inputs.length; i++) {
        inputs[i].checked = false;
        inputs[i].setAttribute("aria-checked", "false");
      }
    } finally {
      suppressDiscreteSwitchProgrammatic = false;
    }
  }

  function applyDiscreteCategorySelectionFromSwitches() {
    if (!gd || !gd.data || !gd.data[0] || scatterTraceColorMode() !== "discrete") return;
    if (!colorDiscreteSwitches) return;
    var inputs = colorDiscreteSwitches.querySelectorAll("input[type=\"checkbox\"]");
    var keySet = new Set();
    for (var i = 0; i < inputs.length; i++) {
      if (inputs[i].checked && inputs[i].dataset.catKey != null) {
        keySet.add(inputs[i].dataset.catKey);
      }
    }
    clearScatterGeometricSelection();
    colorThresholdLevel = null;
    colorThresholdColorCol = null;
    if (colorThresholdStatus) colorThresholdStatus.textContent = "";
    if (!keySet.size) {
      applyRowsSelection([], "Discrete");
      return;
    }
    var cd = gd.data[0].customdata;
    var rows = [];
    for (var j = 0; j < cd.length; j++) {
      var k = discreteCategoryKeyFromCustomRow(cd[j]);
      if (keySet.has(k)) rows.push(cd[j][1]);
    }
    applyRowsSelection(rows, "Discrete");
  }

  function renderDiscreteColorSwitches() {
    if (!colorDiscreteSwitches || !gd || !gd.data || !gd.data[0] || !gd.data[0].customdata) return;
    colorDiscreteSwitches.innerHTML = "";
    var colLabel = covariateDisplayMap[sc.value] || sc.value;
    colorDiscreteSwitches.setAttribute("aria-label", "Filter by " + colLabel);
    var cd = gd.data[0].customdata;
    var markerColors = (gd.data[0].marker && gd.data[0].marker.color) || null;
    var seen = new Set();
    var keyToColor = new Map();
    for (var si = 0; si < cd.length; si++) {
      var kSeen = discreteCategoryKeyFromCustomRow(cd[si]);
      seen.add(kSeen);
      if (!keyToColor.has(kSeen) && Array.isArray(markerColors)) {
        var cval = markerColors[si];
        if (typeof cval === "string" && /^#/.test(cval)) keyToColor.set(kSeen, cval);
      }
    }
    var keys = sortDiscreteCategoryKeys(Array.from(seen));
    var keyCounts = new Map();
    for (var ci = 0; ci < cd.length; ci++) {
      var kCount = discreteCategoryKeyFromCustomRow(cd[ci]);
      keyCounts.set(kCount, (keyCounts.get(kCount) || 0) + 1);
    }
    var labelItems = [];
    for (var ki = 0; ki < keys.length; ki++) {
      var key = keys[ki];
      var safeId = "disc-sw-" + ki + "-" + String(key).replace(/[^a-zA-Z0-9_-]/g, "_");
      var lab = document.createElement("label");
      lab.className = "cryo-color-discrete-switch";
      lab.setAttribute("for", safeId);
      var inp = document.createElement("input");
      inp.type = "checkbox";
      inp.id = safeId;
      inp.setAttribute("role", "switch");
      inp.setAttribute("aria-checked", "false");
      inp.dataset.catKey = key;
      var textWrap = document.createElement("span");
      textWrap.className = "cryo-color-discrete-switch-text";
      var span = document.createElement("span");
      span.className = "cryo-color-discrete-switch-label";
      var labelText = discreteCategoryLabelForKey(key);
      span.textContent = labelText;
      var countSpan = document.createElement("span");
      countSpan.className = "cryo-color-discrete-switch-count";
      var nForKey = keyCounts.get(key) || 0;
      var countText = String(nForKey) + (nForKey === 1 ? " particle" : " particles");
      countSpan.textContent = countText;
      labelItems.push({label: labelText, count: countText});
      var keyHex = keyToColor.get(key);
      if (keyHex) {
        var keyStyles = discreteLabelMontageStyles(String(keyHex));
        if (keyStyles.bg) lab.style.backgroundColor = keyStyles.bg;
        if (keyStyles.covarColor) lab.style.color = keyStyles.covarColor;
      }
      textWrap.appendChild(span);
      textWrap.appendChild(countSpan);
      lab.appendChild(inp);
      lab.appendChild(textWrap);
      colorDiscreteSwitches.appendChild(lab);
      inp.addEventListener("change", (function(input) {
        return function() {
          input.setAttribute("aria-checked", input.checked ? "true" : "false");
          if (suppressDiscreteSwitchProgrammatic) return;
          applyDiscreteCategorySelectionFromSwitches();
        };
      })(inp));
    }
    var colW = measureDiscreteSwitchColumnWidthPx(labelItems);
    if (colW > 0) colorDiscreteSwitches.style.setProperty("--cryo-discrete-col-w", String(colW) + "px");
    else colorDiscreteSwitches.style.removeProperty("--cryo-discrete-col-w");
    syncDiscreteSwitchesFromSelection();
  }

  function syncDiscreteSwitchesFromSelection() {
    if (!colorDiscreteSwitches || !gd || !gd.data || !gd.data[0] || scatterTraceColorMode() !== "discrete") return;
    var cd = gd.data[0].customdata;
    var sel = new Set(saveSelectionRows || []);
    var keyToRows = new Map();
    for (var i = 0; i < cd.length; i++) {
      var k = discreteCategoryKeyFromCustomRow(cd[i]);
      var row = cd[i][1];
      if (!keyToRows.has(k)) keyToRows.set(k, []);
      keyToRows.get(k).push(row);
    }
    suppressDiscreteSwitchProgrammatic = true;
    try {
      var inputs = colorDiscreteSwitches.querySelectorAll("input[type=\"checkbox\"]");
      for (var j = 0; j < inputs.length; j++) {
        var inp = inputs[j];
        var key = inp.dataset.catKey;
        var rowsForK = keyToRows.get(key);
        var full = rowsForK && rowsForK.length
          && rowsForK.every(function(r) { return sel.has(r); });
        inp.checked = !!full;
        inp.setAttribute("aria-checked", full ? "true" : "false");
      }
    } finally {
      suppressDiscreteSwitchProgrammatic = false;
    }
  }

  function syncColorCovariateSidePanel() {
    if (!colorHistPanel) return;
    if (!plotHasScatterData() || !sc.value || sc.value === "none") {
      hideColorHistogram();
      return;
    }
    colorHistPanel.hidden = false;
    colorHistPanel.setAttribute("aria-hidden", "false");
    if (scatterTraceColorMode() === "discrete") {
      if (colorContinuousWrap) colorContinuousWrap.hidden = true;
      if (colorDiscreteWrap) {
        colorDiscreteWrap.hidden = false;
        colorDiscreteWrap.classList.add("cryo-color-discrete-wrap--show");
      }
      renderDiscreteColorSwitches();
    } else {
      if (colorDiscreteWrap) {
        colorDiscreteWrap.hidden = true;
        colorDiscreteWrap.classList.remove("cryo-color-discrete-wrap--show");
      }
      if (colorDiscreteSwitches) colorDiscreteSwitches.innerHTML = "";
      if (colorContinuousWrap) colorContinuousWrap.hidden = false;
      renderColorHistogram();
    }
    syncParticleExplorerRowExpandForColorHist();
  }

  function colorHistogramSelectionBandShape(xRangeOpt) {
    if (colorThresholdLevel == null || colorThresholdColorCol !== sc.value) return null;
    var rng = xRangeOpt;
    if (!rng || rng.length < 2) {
      if (colorHistDiv && colorHistDiv._fullLayout && colorHistDiv._fullLayout.xaxis) {
        rng = colorHistDiv._fullLayout.xaxis.range;
      }
    }
    if (!rng || rng.length < 2) return null;
    var xLo = Math.min(rng[0], rng[1]);
    var xHi = Math.max(rng[0], rng[1]);
    var t = colorThresholdLevel;
    var useMax = !!(colorThresholdUseMax && colorThresholdUseMax.checked);
    var bx0;
    var bx1;
    if (useMax) {
      bx0 = xLo;
      bx1 = Math.min(Math.max(t, xLo), xHi);
    } else {
      bx0 = Math.max(Math.min(t, xHi), xLo);
      bx1 = xHi;
    }
    if (!(bx1 > bx0)) return null;
    return {
      type: "rect",
      xref: "x",
      yref: "paper",
      x0: bx0,
      x1: bx1,
      y0: 0,
      y1: 1,
      fillcolor: "rgba(196, 112, 58, 0.14)",
      line: {width: 0},
      layer: "above"
    };
  }

  function colorHistogramGuideShapes(xRangeOpt) {
    var shapes = colorViolinBaseShapes.slice();
    var band = colorHistogramSelectionBandShape(xRangeOpt);
    if (band) shapes.push(band);
    if (colorThresholdLevel != null && colorThresholdColorCol === sc.value) {
      shapes.push({
        type: "line",
        xref: "x",
        yref: "paper",
        x0: colorThresholdLevel,
        x1: colorThresholdLevel,
        y0: 0,
        y1: 1,
        line: {color: ACCENT, width: 2},
        layer: "above"
      });
    }
    if (colorHistHoverLevel != null) {
      shapes.push({
        type: "line",
        xref: "x",
        yref: "paper",
        x0: colorHistHoverLevel,
        x1: colorHistHoverLevel,
        y0: 0,
        y1: 1,
        line: {color: "#243b53", width: 1.5, dash: "dot"},
        layer: "above"
      });
    }
    return shapes;
  }

  function colorHistogramGuideAnnotations() {
    if (colorHistHoverLevel == null) return [];
    return [{
      xref: "x",
      yref: "paper",
      x: colorHistHoverLevel,
      y: 0,
      yanchor: "top",
      text: formatColorThresholdValue(colorHistHoverLevel),
      showarrow: false,
      bgcolor: "rgba(255,255,255,0.94)",
      bordercolor: "rgba(36,59,83,0.35)",
      borderwidth: 1,
      borderpad: 3,
      font: {size: 10, color: "#243b53"}
    }];
  }

  function relayoutColorHistogramGuides() {
    if (!colorHistDiv || !colorHistDiv.data) return;
    try {
      Plotly.relayout(colorHistDiv, {
        shapes: colorHistogramGuideShapes(),
        annotations: colorHistogramGuideAnnotations()
      });
    } catch (err) {}
  }

  /**
   * Lasso/box selection invalidates color-threshold mode. Drop the threshold guide
   * and status text without rebuilding the violin (avoids O(n log n) histogram work).
   */
  function clearColorHistogramThresholdForGeometricSelection() {
    var hadThresholdGuide = colorThresholdLevel != null && colorThresholdColorCol === sc.value;
    colorThresholdLevel = null;
    colorThresholdColorCol = null;
    if (colorThresholdStatus) colorThresholdStatus.textContent = "";
    resetDiscreteCategorySwitches();
    if (hadThresholdGuide && colorHistDiv && colorHistDiv.data) {
      relayoutColorHistogramGuides();
    }
  }

  function scheduleColorHistogramGuideRelayout(level) {
    colorHistHoverLevel = level;
    if (colorHistHoverRaf != null) return;
    colorHistHoverRaf = requestAnimationFrame(function() {
      colorHistHoverRaf = null;
      relayoutColorHistogramGuides();
    });
  }

  function colorHistogramLevelFromPointerEvent(evt) {
    if (!colorHistDiv || !colorHistDiv._fullLayout) return null;
    var xaxis = colorHistDiv._fullLayout.xaxis;
    var sz = colorHistDiv._fullLayout._size;
    if (!xaxis || !sz || typeof xaxis.p2c !== "function") return null;
    var rect = colorHistDiv.getBoundingClientRect();
    var px = evt.clientX - rect.left - sz.l;
    if (px < 0 || px > sz.w) return null;
    var level = xaxis.p2c(px);
    return typeof level === "number" && isFinite(level) ? level : null;
  }

  function handleColorHistogramMouseMove(evt) {
    scheduleColorHistogramGuideRelayout(colorHistogramLevelFromPointerEvent(evt));
  }

  function handleColorHistogramMouseLeave() {
    scheduleColorHistogramGuideRelayout(null);
  }

  function colorThresholdDescription(nSelected, nTotal) {
    if (colorThresholdLevel == null || colorThresholdColorCol !== sc.value) return "";
    var op = colorThresholdUseMax && colorThresholdUseMax.checked ? "≤" : "≥";
    return "Current threshold: " + op + " " + formatColorThresholdValue(colorThresholdLevel);
  }

  function quantileSorted(sorted, q) {
    if (!sorted.length) return NaN;
    var pos = (sorted.length - 1) * q;
    var lo = Math.floor(pos);
    var hi = Math.ceil(pos);
    if (lo === hi) return sorted[lo];
    var f = pos - lo;
    return sorted[lo] * (1 - f) + sorted[hi] * f;
  }

  function violinDistribution(xs) {
    var n = xs.length;
    var xmin = Infinity, xmax = -Infinity;
    for (var i = 0; i < n; i++) {
      if (xs[i] < xmin) xmin = xs[i];
      if (xs[i] > xmax) xmax = xs[i];
    }
    if (!isFinite(xmin) || !isFinite(xmax)) {
      return {outlineX: [], outlineY: [], fillTraces: [], outlierX: [], outlierY: [], outlierColors: [], range: null};
    }
    var sorted = xs.slice().sort(function(a, b) { return a - b; });
    var q1 = quantileSorted(sorted, 0.25);
    var q3 = quantileSorted(sorted, 0.75);
    var iqr = q3 - q1;
    var loFence = q1 - 3.0 * iqr;
    var hiFence = q3 + 3.0 * iqr;
    var outlierX = [], outlierY = [], vals = [];
    for (var oi = 0; oi < n; oi++) {
      var xv = xs[oi];
      if (iqr > 0 && (xv < loFence || xv > hiFence)) {
        outlierX.push(xv);
        outlierY.push(0);
      } else {
        vals.push(xv);
      }
    }
    if (!vals.length) vals = xs.slice();
    if (xmax <= xmin) {
      var oneColor = paletteScaleCSS(selectedScatterPalette(), 0.5);
      var singleOutlierColors = [];
      for (var soc = 0; soc < outlierX.length; soc++) singleOutlierColors.push(oneColor);
      return {
        outlineX: [xmin, xmin, xmin],
        outlineY: [-0.38, 0.38, -0.38],
        fillTraces: [{
          type: "scatter",
          mode: "lines",
          x: [xmin - 0.5, xmin + 0.5, xmin + 0.5, xmin - 0.5, xmin - 0.5],
          y: [-0.38, -0.38, 0.38, 0.38, -0.38],
          fill: "toself",
          fillcolor: oneColor,
          line: {width: 0, color: oneColor},
          hoverinfo: "skip",
          showlegend: false
        }],
        outlierX: outlierX,
        outlierY: outlierY,
        outlierColors: singleOutlierColors,
        range: [xmin - 0.6, xmax + 0.6]
      };
    }
    var gridN = 96;
    var span = xmax - xmin;
    var mean = 0;
    for (var mi = 0; mi < vals.length; mi++) mean += vals[mi];
    mean /= vals.length;
    var variance = 0;
    for (var vi = 0; vi < vals.length; vi++) {
      var dv = vals[vi] - mean;
      variance += dv * dv;
    }
    var std = Math.sqrt(variance / Math.max(1, vals.length - 1));
    var bandwidth = 1.06 * (std || span / 6) * Math.pow(vals.length, -0.2);
    if (!isFinite(bandwidth) || bandwidth <= 0) bandwidth = span / 24;
    bandwidth = Math.max(bandwidth, span / 200);
    var gridX = [], dens = [], maxD = 0;
    var invTwoBw2 = 1 / (2 * bandwidth * bandwidth);
    for (var g = 0; g < gridN; g++) {
      var gx = xmin + span * g / (gridN - 1);
      var dsum = 0;
      for (var k = 0; k < vals.length; k++) {
        var dx = gx - vals[k];
        dsum += Math.exp(-(dx * dx) * invTwoBw2);
      }
      var d = dsum / vals.length;
      if (d > maxD) maxD = d;
      gridX.push(gx);
      dens.push(d);
    }
    var half = [];
    for (var hi = 0; hi < dens.length; hi++) {
      half.push(maxD > 0 ? 0.42 * dens[hi] / maxD : 0);
    }
    var outlineX = [], outlineY = [];
    for (var ux = 0; ux < gridX.length; ux++) {
      outlineX.push(gridX[ux]);
      outlineY.push(half[ux]);
    }
    for (var lx = gridX.length - 1; lx >= 0; lx--) {
      outlineX.push(gridX[lx]);
      outlineY.push(-half[lx]);
    }
    outlineX.push(gridX[0]);
    outlineY.push(half[0]);
    var fillTraces = [];
    var palName = selectedScatterPalette();
    var outlierColors = [];
    for (var oc = 0; oc < outlierX.length; oc++) {
      outlierColors.push(paletteScaleCSS(palName, (outlierX[oc] - xmin) / span));
    }
    for (var s = 0; s < gridX.length - 1; s++) {
      var x0 = gridX[s], x1 = gridX[s + 1];
      var y = Math.max(half[s], half[s + 1]);
      if (y <= 0) continue;
      var t = ((x0 + x1) * 0.5 - xmin) / span;
      var fillColor = paletteScaleCSS(palName, t);
      fillTraces.push({
        type: "scatter",
        mode: "lines",
        x: [x0, x1, x1, x0, x0],
        y: [-y, -y, y, y, -y],
        fill: "toself",
        fillcolor: fillColor,
        line: {width: 0, color: fillColor},
        hoverinfo: "skip",
        showlegend: false
      });
    }
    var pad = span * 0.04;
    return {
      outlineX: outlineX,
      outlineY: outlineY,
      fillTraces: fillTraces,
      outlierX: outlierX,
      outlierY: outlierY,
      outlierColors: outlierColors,
      range: [xmin - pad, xmax + pad]
    };
  }

  function renderColorHistogram() {
    if (!colorHistPanel || !colorHistDiv) return;
    if (gd && gd.data && gd.data[0] && scatterTraceColorMode() === "discrete") return;
    var items = currentColorValueItems();
    if (!items.length) {
      hideColorHistogram();
      return;
    }
    colorHistPanel.hidden = false;
    colorHistPanel.setAttribute("aria-hidden", "false");
    if (!colorHistEventsWired) {
      colorHistEventsWired = true;
      colorHistDiv.addEventListener("click", handleColorHistogramClick);
      colorHistDiv.addEventListener("mousemove", handleColorHistogramMouseMove);
      colorHistDiv.addEventListener("mouseleave", handleColorHistogramMouseLeave);
    }
    var xs = items.map(function(item) { return item.value; });
    var violin = violinDistribution(xs);
    colorViolinBaseShapes = [];
    var colorColLabel = covariateDisplayMap[sc.value] || sc.value;
    var layout = {
      template: "plotly_white",
      paper_bgcolor: "rgba(250,248,244,0)",
      plot_bgcolor: "rgba(255,255,255,0.62)",
      margin: {l: 42, r: 14, t: 4, b: 28},
      height: 80,
      dragmode: false,
      showlegend: false,
      font: {family: "Barlow, sans-serif", color: "#243b53", size: 11},
      xaxis: {title: {text: "<b>" + escapeHtmlText(colorColLabel) + "</b>", standoff: 2, font: {size: 15.7, weight:
        700}}, zeroline: false, fixedrange: true, range: violin.range},
      yaxis: {visible: false, fixedrange: true, range: [-0.55, 0.55]},
      shapes: colorHistogramGuideShapes(violin.range),
      annotations: colorHistogramGuideAnnotations()
    };
    var data = violin.fillTraces.slice();
    data.push(
      {
        type: "scatter",
        mode: "lines",
        x: violin.outlineX,
        y: violin.outlineY,
        line: {color: "rgba(36,59,83,0.72)", width: 1},
        hoverinfo: "skip"
      },
      {
        type: "scatter",
        mode: "markers",
        x: violin.outlierX,
        y: violin.outlierY,
        marker: {symbol: "circle-open", color: violin.outlierColors, size: 5.76, line: {color: violin.outlierColors,
          width: 1.44}},
        hovertemplate: colorColLabel + " %{x}<br>outlier<extra></extra>"
      }
    );
    var opts = {
      responsive: true,
      displayModeBar: false,
      doubleClick: false,
      scrollZoom: false,
      modeBarButtonsToRemove: ["zoom2d", "pan2d", "select2d", "lasso2d", "autoScale2d", "resetScale2d"]
    };
    try {
      Plotly.react(colorHistDiv, data, layout, opts);
    } catch (err) {
      Plotly.newPlot(colorHistDiv, data, layout, opts);
    }
    if (colorThresholdStatus) {
      if (colorThresholdLevel != null && colorThresholdColorCol === sc.value) {
        colorThresholdStatus.textContent = colorThresholdDescription(saveSelectionRows.length, items.length);
      } else {
        colorThresholdStatus.textContent = "";
      }
    }
    syncParticleExplorerRowExpandForColorHist();
  }

  function traceCustomdataColorDisp(ti) {
    var cd = gd.data[0].customdata;
    if (!cd || ti < 0 || ti >= cd.length) return null;
    var row = cd[ti];
    if (row && row.length > 2 && row[2] != null && row[2] !== "") return row[2];
    return null;
  }

  /** Parse #rgb / #rrggbb → {r,g,b} or null. */
  function parseMarkerHexRgb(hex) {
    if (typeof hex !== "string") return null;
    var h = hex.trim();
    if (h.charAt(0) !== "#") return null;
    h = h.slice(1);
    var r;
    var g;
    var b;
    if (h.length === 3) {
      r = parseInt(h.charAt(0) + h.charAt(0), 16);
      g = parseInt(h.charAt(1) + h.charAt(1), 16);
      b = parseInt(h.charAt(2) + h.charAt(2), 16);
    } else if (h.length === 6) {
      r = parseInt(h.slice(0, 2), 16);
      g = parseInt(h.slice(2, 4), 16);
      b = parseInt(h.slice(4, 6), 16);
    } else {
      return null;
    }
    if (!isFinite(r) || !isFinite(g) || !isFinite(b)) return null;
    return { r: r, g: g, b: b };
  }

  /**
   * k-means / discrete ChimeraX colors: lightly pastel cell fill (small blend toward white)
   * plus a darker saturated tint for the covariate line so it stays readable on the chip.
   */
  function discreteLabelMontageStyles(hex) {
    var rgb = parseMarkerHexRgb(hex);
    if (!rgb) return { bg: "", covarColor: "" };
    var r = rgb.r;
    var g = rgb.g;
    var b = rgb.b;
    var towardWhite = 0.3;
    var rf = Math.round(r * (1 - towardWhite) + 255 * towardWhite);
    var gf = Math.round(g * (1 - towardWhite) + 255 * towardWhite);
    var bf = Math.round(b * (1 - towardWhite) + 255 * towardWhite);
    rf = Math.max(0, Math.min(255, rf));
    gf = Math.max(0, Math.min(255, gf));
    bf = Math.max(0, Math.min(255, bf));
    var textDim = 0.46;
    var rt = Math.max(0, Math.min(255, Math.round(r * textDim)));
    var gt = Math.max(0, Math.min(255, Math.round(g * textDim)));
    var bt = Math.max(0, Math.min(255, Math.round(b * textDim)));
    return {
      bg: "rgb(" + rf + "," + gf + "," + bf + ")",
      covarColor: "rgb(" + rt + "," + gt + "," + bt + ")"
    };
  }

  function discreteMontageSortKey(ti) {
    var d = traceCustomdataColorDisp(ti);
    if (d == null) return { n: NaN, s: "" };
    if (typeof d === "number" && isFinite(d)) return { n: d, s: String(d) };
    var sn = String(d);
    var pn = parseFloat(sn);
    if (sn !== "" && isFinite(pn) && String(pn) === sn) return { n: pn, s: sn };
    return { n: NaN, s: sn };
  }

  function buildTraceMap() {
    if (!preloaded || !gd.data || !gd.data[0]) return;
    var cd = gd.data[0].customdata;
    var tm = [];
    for (var i = 0; i < cd.length; i++) {
      var row = cd[i][1];
      if (preloaded.rows.has(row)) {
        tm.push({ti: i, src: preloaded.rowToSrc.get(row), row: row});
      }
    }
    preloaded.traceMap = tm;
  }

  function showRandomPreloaded() {
    if (!preloaded || !preloaded.traceMap) return;
    /* Only sample from cached rows so the image grid never picks uncached (blank) cells. */
    var pool = (activePool && activePool.length) ? activePool : preloaded.traceMap;
    if (!pool.length) return;
    var shuffled = pool.slice();
    shuffleInPlace(shuffled);
    updateMontage(shuffled.slice(0, Math.min(neighborK, shuffled.length)));
  }

  function updateMontage(nbs) {
    nbs = nbs.slice();
    if (!showGridHighlightsEnabled()) {
      nbs = [];
    }
    var newOrdered = nbs.map(function(x) { return x.row; }).join(",");
    var oldOrdered = lastMontageRows.join(",");
    var newStable = nbs.map(function(x) { return x.row; }).slice().sort(function(a, b) { return a - b; }).join(",");
    var oldStable = lastMontageRows.length
      ? lastMontageRows.slice().sort(function(a, b) { return a - b; }).join(",")
      : "";
    if (newOrdered !== oldOrdered || newStable !== oldStable) {
      invalidateVolumeArtifacts();
    }
    var colorCol = sc.value;
    var colorColLabel = covariateDisplayMap[colorCol] || colorCol;
    var markerColors = gd.data && gd.data[0] && gd.data[0].marker ? gd.data[0].marker.color : null;
    var hasColor = colorCol !== "none" && Array.isArray(markerColors);
    var colorMode = hasColor ? scatterTraceColorMode() : "continuous";
    var cmin = Infinity;
    var cmax = -Infinity;
    if (hasColor && colorMode === "continuous") {
      for (var ci = 0; ci < markerColors.length; ci++) {
        var vci = markerColors[ci];
        if (typeof vci === "number" && isFinite(vci)) {
          if (vci < cmin) cmin = vci;
          if (vci > cmax) cmax = vci;
        }
      }
    }
    if (hasColor && colorMode === "discrete") {
      nbs = nbs.slice().sort(function(a, b) {
        var ka = discreteMontageSortKey(a.ti);
        var kb = discreteMontageSortKey(b.ti);
        if (isFinite(ka.n) && isFinite(kb.n) && ka.n !== kb.n) return ka.n - kb.n;
        if (ka.s !== kb.s) return ka.s < kb.s ? -1 : ka.s > kb.s ? 1 : 0;
        return a.row - b.row;
      });
    } else if (hasColor) {
      nbs = nbs.slice().sort(function(a, b) {
        return markerColors[a.ti] - markerColors[b.ti];
      });
    } else {
      var xs = gd.data && gd.data[0] ? gd.data[0].x : null;
      nbs = nbs.slice().sort(function(a, b) {
        var ax = xs ? xs[a.ti] : 0;
        var bx = xs ? xs[b.ti] : 0;
        if (ax !== bx) return ax - bx;
        return a.row - b.row;
      });
    }
    for (var i = 0; i < montageCells.length; i++) {
      if (i < nbs.length) {
        montageCells[i].img.src = nbs[i].src || "";
        montageCells[i].idx.textContent = "idx " + nbs[i].row;
        if (hasColor) {
          var ti = nbs[i].ti;
          var dispVal = traceCustomdataColorDisp(ti);
          if (dispVal == null && colorMode === "continuous") {
            var vmc = markerColors[ti];
            dispVal = typeof vmc === "number" && isFinite(vmc) ? vmc : vmc;
          }
          var dispStr = dispVal == null ? "—" : (
            typeof dispVal === "number" && isFinite(dispVal) ? dispVal.toFixed(2) : String(dispVal)
          );
          montageCells[i].covar.textContent = colorCol === "labels"
            ? ("k-means cluster = " + dispStr)
            : (colorColLabel + "=" + dispStr);
          if (colorMode === "discrete") {
            var hexCol = markerColors[ti];
            var dStyles = discreteLabelMontageStyles(String(hexCol));
            montageCells[i].cell.style.background = dStyles.bg;
            montageCells[i].covar.style.color = dStyles.covarColor || "";
          } else {
            montageCells[i].covar.style.color = "";
            var val = markerColors[ti];
            var t = cmax > cmin && typeof val === "number" && isFinite(val)
              ? (val - cmin) / (cmax - cmin)
              : 0.5;
            montageCells[i].cell.style.background = covariateScaleCSS(t);
          }
        } else {
          montageCells[i].covar.textContent = "";
          montageCells[i].covar.style.color = "";
          montageCells[i].cell.style.background = "";
        }
      } else {
        montageCells[i].img.src = "";
        montageCells[i].idx.textContent = "";
        montageCells[i].covar.textContent = "";
        montageCells[i].covar.style.color = "";
        montageCells[i].cell.style.background = "";
      }
    }
    if (highlightTraceAdded && gd.data && gd.data[0]) {
      var nx = [], ny = [], txt = [];
      var showPlotHighlights = showGridHighlightsOnPlot();
      for (var k = 0; k < nbs.length; k++) {
        nx.push(gd.data[0].x[nbs[k].ti]);
        ny.push(gd.data[0].y[nbs[k].ti]);
        if (showPlotHighlights) txt.push(labels[k]);
      }
      var restyleData = {
        x: [nx],
        y: [ny],
        mode: [showPlotHighlights ? "markers+text" : "markers"],
        text: [txt],
        visible: [showPlotHighlights]
      };
      if (hasColor) {
        var colorVals = [];
        for (var k = 0; k < nbs.length; k++) colorVals.push(markerColors[nbs[k].ti]);
        restyleData["marker.color"] = [colorVals];
        if (colorMode === "continuous") {
          restyleData["marker.colorscale"] = [selectedScatterPalette()];
          restyleData["marker.cmin"] = cmin;
          restyleData["marker.cmax"] = cmax;
        }
        restyleData["marker.line.color"] = ["#fff"];
        restyleData["marker.line.width"] = [1.5];
      } else {
        restyleData["marker.color"] = [ACCENT];
        restyleData["marker.line.width"] = [0];
      }
      queueHighlightRestyle(restyleData);
    }
    lastMontageRows = [];
    for (var mi = 0; mi < nbs.length; mi++) {
      lastMontageRows.push(nbs[mi].row);
    }
    if (montageDisplayMode === "volumes") {
      var vkey = montageRowsOrderKey();
      if (volumeStaticDataUrls && volumeStaticKey === vkey && volumeStaticDataUrls.length) {
        showVolumesInMontage();
      } else {
        syncVolumeExploreButtons();
      }
    } else {
      syncVolumeExploreButtons();
    }
  }

  /** After color-only scatter reload: same particles, pad from lasso pool or full cache if grid needs more. */
  function restoreMontageAfterScatterReload() {
    if (!preloaded || !preloaded.traceMap) return;
    var rowToItem = new Map();
    for (var ri = 0; ri < preloaded.traceMap.length; ri++) {
      var it0 = preloaded.traceMap[ri];
      rowToItem.set(it0.row, it0);
    }
    var want = neighborK;
    var items = [];
    if (lastMontageRows && lastMontageRows.length) {
      for (var i = 0; i < lastMontageRows.length && items.length < want; i++) {
        var found = rowToItem.get(lastMontageRows[i]);
        if (found) items.push(found);
      }
    }
    var poolItems = (activePool && activePool.length) ? activePool : preloaded.traceMap;
    if (items.length < want && poolItems.length) {
      var used = new Set();
      for (var u = 0; u < items.length; u++) used.add(items[u].row);
      var extras = [];
      for (var pi = 0; pi < poolItems.length; pi++) {
        if (!used.has(poolItems[pi].row)) extras.push(poolItems[pi]);
      }
      for (var ei = extras.length - 1; ei > 0; ei--) {
        var j = Math.floor(Math.random() * (ei + 1));
        var tmp = extras[ei]; extras[ei] = extras[j]; extras[j] = tmp;
      }
      for (var k = 0; k < extras.length && items.length < want; k++) items.push(extras[k]);
    }
    if (items.length) {
      updateMontage(items);
    } else if (activePool && activePool.length >= neighborK) {
      showRandomFromPool();
    } else {
      showRandomPreloaded();
    }
  }

  /** After grid size change: keep current rows where possible; add more only from lasso pool if active, else full
    cache. */
  function refillMontageForNewGridSize() {
    var want = neighborK;
    if (!preloaded || !preloaded.traceMap) {
      if (!plotHasScatterData() || !lastMontageRows.length) return;
      ensureHighlightTrace();
      var picked = pickRandomMontageItemsFilled(want);
      if (picked.length) updateMontage(picked);
      return;
    }
    var poolItems = (activePool && activePool.length) ? activePool : preloaded.traceMap;
    if (!poolItems.length) return;
    var byRow = new Map();
    for (var pi = 0; pi < poolItems.length; pi++) {
      byRow.set(poolItems[pi].row, poolItems[pi]);
    }
    var items = [];
    if (lastMontageRows && lastMontageRows.length) {
      for (var i = 0; i < lastMontageRows.length && items.length < want; i++) {
        var it = byRow.get(lastMontageRows[i]);
        if (it) items.push(it);
      }
    }
    if (items.length < want) {
      var used2 = new Set();
      for (var u2 = 0; u2 < items.length; u2++) used2.add(items[u2].row);
      var extras2 = [];
      for (var ex = 0; ex < poolItems.length; ex++) {
        if (!used2.has(poolItems[ex].row)) extras2.push(poolItems[ex]);
      }
      for (var ei2 = extras2.length - 1; ei2 > 0; ei2--) {
        var j2 = Math.floor(Math.random() * (ei2 + 1));
        var t2 = extras2[ei2]; extras2[ei2] = extras2[j2]; extras2[j2] = t2;
      }
      for (var kk = 0; kk < extras2.length && items.length < want; kk++) {
        items.push(extras2[kk]);
      }
    }
    if (items.length) {
      updateMontage(items);
    } else if (poolItems.length >= want) {
      if (activePool && activePool.length) showRandomFromPool();
      else showRandomPreloaded();
    }
  }

  /* ---------- rendering overlay ---------- */
  function cancelPendingScatterAfterPlot() {
    if (pendingScatterAfterPlot) {
      gd.removeListener("plotly_afterplot", pendingScatterAfterPlot);
      pendingScatterAfterPlot = null;
    }
  }
  function clearScatterPlotWatchdog() {
    if (scatterPlotWatchdog) { clearTimeout(scatterPlotWatchdog); scatterPlotWatchdog = null; }
  }
  function scatterColorByIsDiscrete() {
    /* Matches ``scatter_json``: only ``labels`` uses ChimeraX discrete colors (no sequential palette). */
    return sc.value === "labels";
  }

  function syncScatterPaletteFieldset() {
    if (!paletteFieldset) return;
    var rendering = !!(overlay && overlay.classList.contains("cryo-plot-rendering-overlay--show"));
    var hasColor = sc.value && sc.value !== "none";
    var suppress = rendering || !hasColor || scatterColorByIsDiscrete();
    paletteFieldset.classList.toggle("scatter-palette-radios--suppressed", suppress);
    syncScatterPaletteOptions();
  }
  function syncScatterPaletteOptions() {
    if (!scatterPaletteOptions) return;
    var showOptions = !!(
      scatterPaletteToggle
      && scatterPaletteToggle.getAttribute("aria-expanded") === "true"
      && paletteFieldset
      && !paletteFieldset.classList.contains("scatter-palette-radios--suppressed")
    );
    scatterPaletteOptions.classList.toggle("cryo-palette-select__options--collapsed", !showOptions);
  }
  function setRendering(on) {
    if (overlay) {
      overlay.classList.toggle("cryo-plot-rendering-overlay--show", on);
      overlay.setAttribute("aria-hidden", on ? "false" : "true");
    }
    syncScatterPaletteFieldset();
  }

  /* ---------- lasso / selection state (shared by scatter load + events) ---------- */
  var activePool = null;
  var lassoSelectedRowSet = null;
  var lassoSelectionDebounceTimer = null;
  var LASSO_SELECTION_DEBOUNCE_MS = 200;
  var highlightRestyleRaf = null;
  var pendingHighlightRestyle = null;

  function clearScatterGeometricSelection() {
    if (!gd) return;
    try {
      Plotly.relayout(gd, { selections: [] });
    } catch (e) {
      try {
        Plotly.relayout(gd, { selections: [], selectionrevision: Date.now() });
      } catch (e2) {}
    }
  }

  function applyPreselectedTraceIndices() {
    var meta = gd.layout && gd.layout.meta;
    var pre = meta && meta.cdrgn_preselected;
    if (!pre || !pre.length) return;
    try {
      Plotly.restyle(gd, { selectedpoints: [pre] }, [0]);
    } catch (err) {
      console.error(err);
    }
  }

  function applySavedSelectionToCurrentTrace() {
    if (!gd || !gd.data || !gd.data[0] || !gd.data[0].customdata) return;
    var rows = saveSelectionRows && saveSelectionRows.length ? saveSelectionRows : [];
    if (!rows.length) {
      try { Plotly.restyle(gd, { selectedpoints: [null] }, [0]); } catch (err) {}
      lassoSelectedRowSet = null;
      return;
    }
    var rowSet = new Set(rows);
    var cd = gd.data[0].customdata;
    var selected = [];
    for (var i = 0; i < cd.length; i++) {
      var row = cd[i] && cd[i].length > 1 ? cd[i][1] : null;
      if (rowSet.has(row)) selected.push(i);
    }
    try {
      Plotly.restyle(gd, { selectedpoints: [selected] }, [0]);
    } catch (err) {
      console.error(err);
    }
    lassoSelectedRowSet = new Set();
    for (var si = 0; si < selected.length; si++) {
      var r = cd[selected[si]][1];
      lassoSelectedRowSet.add(r);
    }
  }

  /* ---------- highlight trace + montage grid coupling ---------- */
  function showGridHighlightsEnabled() {
    return imageGridMenuIsOpen();
  }

  function showGridHighlightsOnPlot() {
    return !suppressPlotGridHighlights && showGridHighlightsEnabled();
  }

  function syncHighlightTraceAnnotations() {
    if (!highlightTraceAdded || !gd.data || !gd.data[1]) return;
    var showHighlights = showGridHighlightsOnPlot();
    var text = [];
    if (showHighlights && lastMontageRows && lastMontageRows.length) {
      for (var i = 0; i < lastMontageRows.length && i < labels.length; i++) {
        text.push(labels[i]);
      }
    }
    try {
      Plotly.restyle(gd, {
        mode: [showHighlights ? "markers+text" : "markers"],
        text: [text],
        visible: [showHighlights]
      }, [1]);
    } catch (err) {
      console.error(err);
    }
  }

  function ensureHighlightTrace() {
    if (highlightTraceAdded) return;
    highlightTraceAdded = true;
    Plotly.addTraces(gd, [{
      x: [], y: [],
      mode: showGridHighlightsOnPlot() ? "markers+text" : "markers",
      text: [],
      textposition: "middle center",
      textfont: {size: 9, color: "#fff", family: "Georgia, 'Times New Roman', serif"},
      marker: {size: 15, color: "#c4703a", line: {width: 0}},
      visible: showGridHighlightsOnPlot(),
      hoverinfo: "skip",
      showlegend: false
    }]);
  }

  /* ---------- scatter load ---------- */
  function loadPlot(axesChanged) {
    var gen = ++scatterLoadGeneration;
    var preservedDragMode = (
      gd && gd.layout && gd.layout.dragmode
    ) ? String(gd.layout.dragmode) : "";
    var preservedSelections = null;
    if (!axesChanged && gd && gd.layout && gd.layout.selections && gd.layout.selections.length) {
      try {
        preservedSelections = JSON.parse(JSON.stringify(gd.layout.selections));
      } catch (e) {
        preservedSelections = null;
      }
    }
    var hasActiveSelection = !!(saveSelectionRows && saveSelectionRows.length);
    suppressSelectionEvents = true;
    highlightTraceAdded = false;
    cancelLassoSelectionDebounce();
    cancelPendingHighlightRestyle();
    cancelPendingScatterAfterPlot();
    clearScatterPlotWatchdog();
    if (plotStatus) plotStatus.textContent = "";
    hideColorHistogram();
    setRendering(true);
    scatterPlotWatchdog = setTimeout(function() {
      scatterPlotWatchdog = null;
      if (!overlay || !overlay.classList.contains("cryo-plot-rendering-overlay--show")) return;
      cancelPendingScatterAfterPlot();
      setRendering(false);
      suppressSelectionEvents = false;
      if (plotStatus) plotStatus.textContent = "Plot took too long — try a hard refresh.";
    }, SCATTER_PLOT_TIMEOUT_MS);
    var x = sx.value, y = sy.value, c = sc.value;
    var q = "x=" + encodeURIComponent(x)
      + "&y=" + encodeURIComponent(y)
      + "&color=" + encodeURIComponent(c)
      + "&explorer_scatter=1"
      + "&marker_size=2"
      + "&palette=" + encodeURIComponent(selectedScatterPalette());
    if (preselectRowsForNextFetch.length) {
      q += "&preselect_rows=" + preselectRowsForNextFetch.map(Number).join(",");
    }
    requestAnimationFrame(function() {
      requestAnimationFrame(function() {
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
            fig.layout.title = "";
            fig.layout.margin = Object.assign({}, fig.layout.margin || {}, {
              l: 38,
              r: 10,
              t: 4,
              b: 38,
              pad: 0
            });
            var nextDragMode = preservedDragMode;
            if (hasActiveSelection) {
              if (nextDragMode !== "lasso" && nextDragMode !== "select") {
                nextDragMode = "select";
              }
            }
            if (nextDragMode) {
              fig.layout.dragmode = nextDragMode;
            }
            var axFont = {size: 21.45, family: "Barlow, sans-serif", color: "#243b53", weight: 700};
            fig.layout.xaxis = Object.assign(fig.layout.xaxis || {}, {
              title: { text: "<b>" + escapeHtmlText(x) + "</b>", font: axFont }
            });
            fig.layout.yaxis = Object.assign(fig.layout.yaxis || {}, {
              title: { text: "<b>" + escapeHtmlText(y) + "</b>", font: axFont }
            });
            if (/^umap/i.test(x)) fig.layout.xaxis.showticklabels = false;
            if (/^umap/i.test(y)) fig.layout.yaxis.showticklabels = false;
            function finishScatterDraw() {
              if (gen !== scatterLoadGeneration) return;
              cancelPendingScatterAfterPlot();
              wireUpPlotlyEvents();
              suppressSelectionEvents = true;
              renderedScatterAxes = {x: x, y: y};
              clearScatterPlotWatchdog();
              setRendering(false);
              buildTraceMap();
              updateParticleSelFieldset();
              if (preselectRowsForNextFetch.length) {
                saveSelectionRows = preselectRowsForNextFetch.slice();
                updateParticleSelFieldset();
                applyPreselectedTraceIndices();
                preselectRowsForNextFetch = [];
              }
              applySavedSelectionToCurrentTrace();
              if (preservedSelections && preservedSelections.length) {
                try {
                  var srev = gd.layout && gd.layout.selectionrevision != null
                    ? gd.layout.selectionrevision
                    : 0;
                  Plotly.relayout(gd, {
                    selections: preservedSelections,
                    selectionrevision: srev + 1
                  });
                } catch (e) {}
              }
              syncColorCovariateSidePanel();
              if (imagesViewEnabled && preloaded && preloaded.traceMap) {
                if (lassoSelectedRowSet && lassoSelectedRowSet.size) {
                  activePool = preloaded.traceMap.filter(function(item) {
                    return lassoSelectedRowSet.has(item.row);
                  });
                  if (preloadStatus) {
                    if (activePool.length) {
                      preloadStatus.textContent = formatCachedSelectionStatusCore(activePool.length) + ".";
                    } else {
                      setPreloadStatusCachedTotalOnly();
                    }
                  }
                }
                if (lastMontageRows.length) {
                  restoreMontageAfterScatterReload();
                } else {
                  showRandomPreloaded();
                }
              } else if (!imagesViewEnabled) {
                if (lastMontageRows.length) {
                  ensureHighlightTrace();
                  var keepItems = [];
                  for (var li = 0; li < lastMontageRows.length; li++) {
                    var itKeep = montageItemForRow(lastMontageRows[li]);
                    if (itKeep) keepItems.push(itKeep);
                  }
                  if (keepItems.length) updateMontage(keepItems);
                } else {
                  updateMontage([]);
                }
              }
              if (showVolumeExplorer) syncVolumeExploreButtons();
              var montageOverlayBusy = montagePreloadOverlayEl
                && montagePreloadOverlayEl.classList.contains("cryo-plot-rendering-overlay--show");
              if (!montageOverlayBusy) syncImageCacheButton(false);
              syncScatterControlsAlignment();
              requestAnimationFrame(function() {
                try { Plotly.Plots.resize(gd); } catch (e) {}
                if (colorHistPanel && !colorHistPanel.hidden && colorHistDiv) {
                  try { Plotly.Plots.resize(colorHistDiv); } catch (e2) {}
                }
              });
              setTimeout(function() { suppressSelectionEvents = false; }, 0);
            }
            var opts = {responsive: true};
            var p;
            try {
              if (!gd.data || gd.data.length === 0) {
                p = Plotly.newPlot(gd, fig.data, fig.layout, opts);
              } else {
                p = Plotly.react(gd, fig.data, fig.layout, opts);
              }
            } catch (err) { throw err; }
            pendingScatterAfterPlot = finishScatterDraw;
            gd.on("plotly_afterplot", pendingScatterAfterPlot);
            setTimeout(function() {
              if (gen !== scatterLoadGeneration) return;
              if (pendingScatterAfterPlot && gd.data && gd.data.length) finishScatterDraw();
            }, 800);
            if (p && typeof p.then === "function") {
              return p.catch(function(err) { cancelPendingScatterAfterPlot(); throw err; });
            }
            return Promise.resolve();
          })
          .catch(function(e) {
            console.error(e);
            cancelPendingScatterAfterPlot();
            clearScatterPlotWatchdog();
            setRendering(false);
            suppressSelectionEvents = false;
            if (plotStatus) plotStatus.textContent = e.message || "Could not load scatter.";
          });
      });
    });
  }

  function cancelLassoSelectionDebounce() {
    if (lassoSelectionDebounceTimer) {
      clearTimeout(lassoSelectionDebounceTimer);
      lassoSelectionDebounceTimer = null;
    }
  }

  function cancelPendingHighlightRestyle() {
    if (highlightRestyleRaf != null) {
      cancelAnimationFrame(highlightRestyleRaf);
      highlightRestyleRaf = null;
    }
    pendingHighlightRestyle = null;
  }

  /** Copy indices immediately — Plotly may reuse the event object on the next tick. */
  function snapshotSelectedPoints(ev) {
    if (!ev || !ev.points || !ev.points.length) return [];
    var pts = ev.points;
    var out = new Array(pts.length);
    for (var i = 0; i < pts.length; i++) {
      var p = pts[i];
      out[i] = {curveNumber: p.curveNumber, pointIndex: p.pointIndex};
    }
    return out;
  }

  function snapshotRowsFromPoints(ev) {
    var set = new Set();
    if (ev && ev.points) {
      for (var i = 0; i < ev.points.length; i++) {
        var p = ev.points[i];
        if (p.curveNumber !== 0) continue;
        if (p.customdata != null && p.customdata.length > 1) {
          set.add(p.customdata[1]);
        }
      }
    }
    return Array.from(set);
  }

  function selectionRowsToTraceIndices(rows) {
    var rowSet = new Set(rows || []);
    var selected = [];
    if (!gd || !gd.data || !gd.data[0] || !gd.data[0].customdata) return selected;
    var cd = gd.data[0].customdata;
    for (var i = 0; i < cd.length; i++) {
      var row = cd[i] && cd[i].length > 1 ? cd[i][1] : null;
      if (rowSet.has(row)) selected.push(i);
    }
    return selected;
  }

  function applyRowsSelection(rows, statusPrefix) {
    rows = rows || [];
    var rowSet = new Set(rows);
    saveSelectionRows = Array.from(rowSet);
    selectionTooltipText = "";
    updateParticleSelFieldset();
    var selectedTi = selectionRowsToTraceIndices(saveSelectionRows);
    suppressSelectionEvents = true;
    try {
      Plotly.restyle(gd, { selectedpoints: [selectedTi.length ? selectedTi : null] }, [0]);
    } catch (err) {
      console.error(err);
    }
    setTimeout(function() { suppressSelectionEvents = false; }, 0);

    if (!saveSelectionRows.length) {
      activePool = null;
      lassoSelectedRowSet = null;
      if (imagesViewEnabled && preloaded && preloadStatus) {
        setPreloadStatusCachedTotalOnly();
      }
      invalidateVolumeArtifacts();
      return;
    }

    lassoSelectedRowSet = new Set(saveSelectionRows);
    if (imagesViewEnabled && preloaded && preloaded.traceMap) {
      activePool = preloaded.traceMap.filter(function(item) {
        return rowSet.has(item.row);
      });
      if (preloadStatus) {
        preloadStatus.textContent = (statusPrefix || "Selection") + ": "
          + formatCachedSelectionStatusCore(activePool.length) + ".";
      }
      if (activePool.length) {
        showRandomFromPool();
      } else {
        invalidateVolumeArtifacts();
      }
    } else if (lastMontageRows.length || montageDisplayMode === "volumes") {
      ensureHighlightTrace();
      var picked = pickRandomMontageItemsFilled(neighborK);
      if (picked.length) updateMontage(picked);
      else invalidateVolumeArtifacts();
    } else {
      invalidateVolumeArtifacts();
    }
  }

  function clearExplorerSelection() {
    if (!gd) return;
    cancelLassoSelectionDebounce();
    hideBoxSelectTooltip();
    selectionTooltipText = "";
    colorThresholdLevel = null;
    colorThresholdColorCol = null;
    if (colorThresholdStatus) colorThresholdStatus.textContent = "";
    resetDiscreteCategorySwitches();
    if (colorHistDiv && colorHistDiv.data) {
      try { relayoutColorHistogramGuides(); } catch (e) {}
    }
    clearScatterGeometricSelection();
    applyRowsSelection([]);
    if (imagesViewEnabled && preloaded) {
      showRandomPreloaded();
    }
  }

  function rowsMatchingColorThresholdFromTrace(items, level, useMax) {
    var rows = [];
    for (var i = 0; i < items.length; i++) {
      if (useMax ? items[i].value <= level : items[i].value >= level) rows.push(items[i].row);
    }
    return rows;
  }

  function applyColorThresholdSelection(level) {
    var items = currentColorValueItems();
    if (!items.length || typeof level !== "number" || !isFinite(level)) return;
    cancelLassoSelectionDebounce();
    resetDiscreteCategorySwitches();
    if (colorDiscreteWrap) {
      colorDiscreteWrap.hidden = true;
      colorDiscreteWrap.classList.remove("cryo-color-discrete-wrap--show");
    }
    if (colorContinuousWrap) colorContinuousWrap.hidden = false;
    clearScatterGeometricSelection();
    var useMax = !!(colorThresholdUseMax && colorThresholdUseMax.checked);
    colorThresholdLevel = level;
    colorThresholdColorCol = sc.value;
    fetch(b.urls.apiCovariateThresholdRows, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        column: sc.value,
        level: level,
        use_max: useMax
      })
    })
      .then(function(r) {
        return r.json().then(function(j) { return { ok: r.ok, j: j }; });
      })
      .then(function(res) {
        var rows = (res.ok && res.j && Array.isArray(res.j.rows))
          ? res.j.rows
          : rowsMatchingColorThresholdFromTrace(items, level, useMax);
        applyRowsSelection(rows, "Color threshold");
        syncColorCovariateSidePanel();
      })
      .catch(function() {
        var rows = rowsMatchingColorThresholdFromTrace(items, level, useMax);
        applyRowsSelection(rows, "Color threshold");
        syncColorCovariateSidePanel();
      });
  }

  function handleColorHistogramClick(evt) {
    evt.preventDefault();
    var level = colorHistogramLevelFromPointerEvent(evt);
    if (level == null) return;
    applyColorThresholdSelection(level);
  }

  function queueHighlightRestyle(restyleData) {
    pendingHighlightRestyle = restyleData;
    if (highlightRestyleRaf != null) return;
    highlightRestyleRaf = requestAnimationFrame(function() {
      highlightRestyleRaf = null;
      var payload = pendingHighlightRestyle;
      pendingHighlightRestyle = null;
      if (!payload || !highlightTraceAdded || !gd.data || !gd.data[0]) return;
      try {
        Plotly.restyle(gd, payload, [1]);
      } catch (err) {
        console.error(err);
      }
    });
  }

  function formatSelectionBound(v) {
    var n = Number(v);
    if (!isFinite(n)) return String(v);
    return Number(n.toPrecision(6)).toString();
  }

  function boxSelectionTooltip(ev) {
    if (!ev || !ev.range || !ev.range.x || !ev.range.y) return "";
    var xr = ev.range.x.slice().sort(function(a, b) { return a - b; });
    var yr = ev.range.y.slice().sort(function(a, b) { return a - b; });
    if (xr.length < 2 || yr.length < 2) return "";
    return "Box filter: "
      + formatSelectionBound(xr[0]) + " <= val_x <= " + formatSelectionBound(xr[1])
      + "; "
      + formatSelectionBound(yr[0]) + " <= val_y <= " + formatSelectionBound(yr[1])
      + " (val_x = " + sx.value + ", val_y = " + sy.value + ")";
  }

  function hideBoxSelectTooltip() {
    if (!boxSelectTooltipEl) return;
    boxSelectTooltipEl.hidden = true;
    boxSelectTooltipEl.textContent = "";
  }

  function showBoxSelectTooltip(ev) {
    if (!boxSelectTooltipEl || !scatterPlotStack) return;
    var text = boxSelectionTooltip(ev);
    if (!text) {
      hideBoxSelectTooltip();
      return;
    }
    var rect = scatterPlotStack.getBoundingClientRect();
    var src = ev && ev.event ? ev.event : null;
    var left = src && typeof src.clientX === "number" ? src.clientX - rect.left : 12;
    var top = src && typeof src.clientY === "number" ? src.clientY - rect.top : 12;
    left = Math.max(8, Math.min(left, rect.width - 24));
    top = Math.max(8, Math.min(top, rect.height - 8));
    boxSelectTooltipEl.textContent = text;
    boxSelectTooltipEl.style.left = left + "px";
    boxSelectTooltipEl.style.top = top + "px";
    boxSelectTooltipEl.hidden = false;
  }

  function findKNearestInPool(px, py, k) {
    var items = activePool || (preloaded && preloaded.traceMap) || [];
    if (!items.length) return [];
    var xs = gd.data[0].x, ys = gd.data[0].y;
    var best = [];
    for (var j = 0; j < items.length; j++) {
      var ti = items[j].ti;
      var ddx = xs[ti] - px, ddy = ys[ti] - py;
      var d2 = ddx * ddx + ddy * ddy;
      if (best.length < k) {
        best.push({j: j, d: d2});
        if (best.length === k) best.sort(function(a, b) { return a.d - b.d; });
      } else if (d2 < best[k - 1].d) {
        best[k - 1] = {j: j, d: d2};
        best.sort(function(a, b) { return a.d - b.d; });
      }
    }
    return best.map(function(b) { return items[b.j]; });
  }

  /* ---------- plot events ---------- */
  var plotlyEventsWired = false;
  function wireUpPlotlyEvents() {
    if (plotlyEventsWired) return;
    plotlyEventsWired = true;
    ensureHighlightTrace();
    gd.on("plotly_click", function(ev) {
      if (!imagesViewEnabled || !preloaded || !preloaded.traceMap) return;
      if (!showGridHighlightsEnabled()) return;
      if (!ev.points || !ev.points.length) return;
      var pt = ev.points[0];
      if (pt.curveNumber !== 0) return;
      var nbs = findKNearestInPool(pt.x, pt.y, neighborK);
      if (nbs.length) updateMontage(nbs);
    });
    gd.on("plotly_selecting", function(ev) {
      if (suppressSelectionEvents) return;
      showBoxSelectTooltip(ev);
    });
    function applyLassoSelectionFromSnapshot(snap) {
      if (!preloaded || !preloaded.traceMap) return;
      if (!snap.length) {
        activePool = null;
        lassoSelectedRowSet = null;
        if (preloadStatus) setPreloadStatusCachedTotalOnly();
        return;
      }
      var selectedTi = new Set();
      for (var i = 0; i < snap.length; i++) {
        if (snap[i].curveNumber === 0) selectedTi.add(snap[i].pointIndex);
      }
      activePool = preloaded.traceMap.filter(function(item) {
        return selectedTi.has(item.ti);
      });
      lassoSelectedRowSet = new Set();
      for (var lp = 0; lp < activePool.length; lp++) {
        lassoSelectedRowSet.add(activePool[lp].row);
      }
      if (preloadStatus) {
        preloadStatus.textContent = formatCachedSelectionStatusCore(activePool.length) + ".";
      }
      if (activePool.length) {
        showRandomFromPool();
      }
    }

    gd.on("plotly_selected", function(ev) {
      if (suppressSelectionEvents) return;
      clearColorHistogramThresholdForGeometricSelection();
      var snap = snapshotSelectedPoints(ev);
      var rowsSnap = snapshotRowsFromPoints(ev);
      selectionTooltipText = boxSelectionTooltip(ev);
      hideBoxSelectTooltip();
      cancelLassoSelectionDebounce();
        lassoSelectionDebounceTimer = setTimeout(function() {
        lassoSelectionDebounceTimer = null;
        saveSelectionRows = rowsSnap;
        updateParticleSelFieldset();
        if (imagesViewEnabled && preloaded && preloaded.traceMap) {
          applyLassoSelectionFromSnapshot(snap);
        }
      }, LASSO_SELECTION_DEBOUNCE_MS);
    });
    gd.on("plotly_deselect", function() {
      if (suppressSelectionEvents) return;
      clearColorHistogramThresholdForGeometricSelection();
      hideBoxSelectTooltip();
      cancelLassoSelectionDebounce();
      saveSelectionRows = [];
      selectionTooltipText = "";
      updateParticleSelFieldset();
      activePool = null;
      lassoSelectedRowSet = null;
      if (imagesViewEnabled && preloaded && preloadStatus) {
        setPreloadStatusCachedTotalOnly();
      }
    });
    gd.addEventListener("mouseleave", hideBoxSelectTooltip);
  }

  function showRandomFromPool() {
    var pool = (activePool || (preloaded && preloaded.traceMap) || []).slice();
    if (!pool.length) return;
    for (var i = pool.length - 1; i > 0; i--) {
      var j = Math.floor(Math.random() * (i + 1));
      var tmp = pool[i]; pool[i] = pool[j]; pool[j] = tmp;
    }
    updateMontage(pool.slice(0, neighborK));
  }

  function triggerRedraw(preserveSelectionAndMontage) {
    var axesChanged = sx.value !== renderedScatterAxes.x || sy.value !== renderedScatterAxes.y;
    var keepState = !!preserveSelectionAndMontage && !axesChanged;
    if (keepState) {
      /* Re-apply selection highlight after color-only redraw (keep URI size bounded). */
      preselectRowsForNextFetch = (saveSelectionRows && saveSelectionRows.length && saveSelectionRows.length <= 5000)
        ? saveSelectionRows.slice()
        : [];
    } else {
      saveSelectionRows = [];
      preselectRowsForNextFetch = [];
      activePool = null;
      lassoSelectedRowSet = null;
      selectionTooltipText = "";
      lastMontageRows = [];
      colorThresholdLevel = null;
      colorThresholdColorCol = null;
    }
    updateParticleSelFieldset();
    invalidateVolumeArtifacts();
    if (axesChanged) {
      setImageCacheProgress(false);
      if (preloaded) {
        syncImageCacheButton(false);
        if (preloadStatus) setPreloadStatusCachedTotalOnly();
      } else {
        imagesViewEnabled = false;
        montageInner.classList.add("cryo-dash-montage-inner--disabled");
        gridSizeSelect.disabled = true;
        setImageGridMenuOpen(false);
        syncImageCacheButton(false);
        if (preloadStatus) preloadStatus.textContent = "";
      }
    }
    cancelLassoSelectionDebounce();
    cancelPendingHighlightRestyle();
    plotlyEventsWired = false;
    if (axesChanged) initializeExplorerView();
    else loadPlot(axesChanged);
  }

  function saveInverseChecked() {
    var el = document.getElementById("save-inverse");
    return el && el.checked;
  }

  function openSelectionFileBrowser(startDir) {
    if (!selFileBrowserPanel || !selFbList || !selFbPath) return;
    if (selFbSaveName && !String(selFbSaveName.value || "").trim()) {
      selFbSaveName.value = "indices.pkl";
    }
    selFileBrowserPanel.hidden = false;
    selFileBrowserPanel.setAttribute("aria-hidden", "false");
    loadSelectionFileBrowserDir(startDir || "");
  }

  function closeSelectionFileBrowser() {
    if (!selFileBrowserPanel) return;
    selFileBrowserPanel.hidden = true;
    selFileBrowserPanel.setAttribute("aria-hidden", "true");
    selFbCurrentDir = null;
  }

  function normalizedSelectionBasename() {
    var raw = String((selFbSaveName && selFbSaveName.value) || "").trim();
    if (!raw) raw = "indices.pkl";
    if (raw.indexOf("/") >= 0 || raw.indexOf("\\") >= 0) return null;
    if (!/\.pkl$/i.test(raw)) raw += ".pkl";
    var base = raw.replace(/\.pkl$/i, "");
    return base || "indices";
  }

  function loadSelectionFileBrowserDir(dir) {
    if (!selFbList || !selFbPath) return;
    selFbList.innerHTML = "<li class='cryo-file-browser-empty'>Loading…</li>";
    var q = dir ? ("?dir=" + encodeURIComponent(dir)) : "";
    fetch(b.urls.apiListServerFiles + q)
      .then(function(r) { return r.json(); })
      .then(function(j) {
        if (!j.ok) {
          selFbList.innerHTML = "<li class='cryo-file-browser-empty'>" + (j.error || "Error") + "</li>";
          return;
        }
        selFbCurrentDir = j.dir;
        selFbPath.textContent = j.dir;
        selFbPath.title = j.dir;
        selFbList.innerHTML = "";
        if (!j.entries || !j.entries.length) {
          selFbList.innerHTML = "<li class='cryo-file-browser-empty'>No entries here</li>";
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
              loadSelectionFileBrowserDir(j.dir + "/" + ent.name);
            });
          } else {
            li.addEventListener("click", function() {
              if (selFbSaveName) selFbSaveName.value = ent.name;
            });
          }
          selFbList.appendChild(li);
        });
      })
      .catch(function() {
        selFbList.innerHTML = "<li class='cryo-file-browser-empty'>Could not list directory</li>";
      });
  }

  function saveSelectionViaBrowser() {
    var status = document.getElementById("save-status");
    if (!selFbCurrentDir) {
      status.textContent = "Choose a destination folder first.";
      return;
    }
    var basename = normalizedSelectionBasename();
    if (!basename) {
      status.textContent = "Filename cannot contain path separators.";
      return;
    }
    closeSelectionFileBrowser();
    postSaveSelection({ force: false, sel_dir: selFbCurrentDir, basename: basename });
  }

  function postSaveSelection(payload) {
    var status = document.getElementById("save-status");
    status.textContent = "";
    status.title = "";
    if (!saveSelectionRows.length) {
      status.textContent = "Select particles (plot, color histogram, or discrete toggles) before saving.";
      status.title = status.textContent;
      return;
    }
    payload.rows = saveSelectionRows;
    payload.save_inverse = saveInverseChecked();
    fetch(b.urls.apiSaveSelection, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload)
    })
      .then(function(r) { return r.json().then(function(j) { return { ok: r.ok, j: j }; }); })
      .then(function(res) {
        if (!res.ok) {
          status.textContent = res.j.error || "Save failed.";
          status.title = status.textContent;
          return;
        }
        var msg = "Saved " + res.j.n_selected + " indices to " + res.j.path;
        if (res.j.inverse_path) {
          msg += "; inverse to " + res.j.inverse_path;
        }
        status.textContent = msg + ".";
        status.title = status.textContent;
      })
      .catch(function() {
        status.textContent = "Request failed.";
        status.title = status.textContent;
      });
  }

  document.getElementById("save-indices-pkl").addEventListener("click", function() {
    postSaveSelection({ force: true, sel_dir: "", basename: "indices" });
  });
  if (clearExplorerSelectionBtn) {
    clearExplorerSelectionBtn.addEventListener("click", function() {
      clearExplorerSelection();
    });
  }
  document.getElementById("save-indices-custom").addEventListener("click", function() {
    if (selFbSaveName) selFbSaveName.value = "indices.pkl";
    openSelectionFileBrowser(workdirForSaveHint || "");
  });
  if (selFbUp) {
    selFbUp.addEventListener("click", function() {
      if (selFbCurrentDir) {
        var parent = selFbCurrentDir.replace(/\/[^\/]*\/?$/, "") || "/";
        loadSelectionFileBrowserDir(parent);
      }
    });
  }
  if (selFbCancel) {
    selFbCancel.addEventListener("click", function() {
      closeSelectionFileBrowser();
    });
  }
  if (selFbSaveBtn) {
    selFbSaveBtn.addEventListener("click", function() {
      saveSelectionViaBrowser();
    });
  }
  if (selFbSaveName) {
    selFbSaveName.addEventListener("keydown", function(evt) {
      if (evt.key === "Enter") {
        evt.preventDefault();
        saveSelectionViaBrowser();
      }
    });
  }

  [sx, sy].forEach(function(s) {
    s.addEventListener("change", function() { triggerRedraw(false); });
  });
  sc.addEventListener("change", function() {
    colorThresholdLevel = null;
    colorThresholdColorCol = null;
    triggerRedraw(true);
  });
  document.querySelectorAll("input[name=\"scatter_palette\"]").forEach(function(inp) {
    inp.addEventListener("change", function() { triggerRedraw(true); });
  });
  if (scatterPaletteToggle) {
    scatterPaletteToggle.addEventListener("click", function() {
      var expanded = scatterPaletteToggle.getAttribute("aria-expanded") === "true";
      scatterPaletteToggle.setAttribute("aria-expanded", expanded ? "false" : "true");
      syncScatterPaletteOptions();
    });
  }
  if (imageGridMenuToggle) {
    imageGridMenuToggle.addEventListener("click", function() {
      setImageGridMenuOpen(!imageGridMenuIsOpen());
    });
  }
  if (colorThresholdUseMax) {
    colorThresholdUseMax.addEventListener("change", function() {
      if (colorThresholdLevel != null && colorThresholdColorCol === sc.value) {
        applyColorThresholdSelection(colorThresholdLevel);
      } else {
        syncColorCovariateSidePanel();
      }
    });
  }
  updateParticleSelFieldset();
  if (showVolumeExplorer) syncVolumeExploreButtons();
  syncScatterControlsAlignment();
  equalizeMontageActionButtonWidths();
  window.addEventListener("resize", equalizeMontageActionButtonWidths);
  syncScatterPaletteFieldset();
  initializeExplorerView();
})();
