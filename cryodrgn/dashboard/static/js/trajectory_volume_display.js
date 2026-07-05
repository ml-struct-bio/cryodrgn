/**
 * Trajectory creator integrated volume panel: VTK 3D, 2D slice, or ChimeraX PNGs.
 */
(function (global) {
  "use strict";

  var VTK_BUNDLE_URL = null;
  var TRAJ_VSLICE_PAN_STEP = 3;

  function TrajectoryVolumeDisplay(options) {
    this.asideHostEl = options.asideHostEl;
    this.viewportEl = options.viewportEl;
    this.displayRowEl = options.displayRowEl;
    this.padColumnEl = options.padColumnEl;
    this.canvasEl = options.canvasEl;
    this.vtkContainerEl = options.vtkContainerEl;
    this.controlsDockEl = options.controlsDockEl;
    this.vtkSliceControlsEl = options.vtkSliceControlsEl;
    this.chimeraxViewControlsEl = options.chimeraxViewControlsEl;
    this.progressEl = options.progressEl;
    this.renderingOverlayEl = options.renderingOverlayEl;
    this.statusEl = options.statusEl;
    this.viewMode3dEl = options.viewMode3dEl;
    this.isoControlsEl = options.isoControlsEl;
    this.isoSliderEl = options.isoSliderEl;
    this.sliceControlsRowEl = options.sliceControlsRowEl;
    this.rotationLockToolbarEl = options.rotationLockToolbarEl;
    this.volumeNavEl = options.volumeNavEl;
    this.volumeSliderEl = options.volumeSliderEl;
    this.volumeSliderTicksEl = options.volumeSliderTicksEl;
    this.btnVolPrev = options.btnVolPrev;
    this.btnVolNext = options.btnVolNext;
    this.btnResetView = options.btnResetView;
    this.expandedHostEl = options.expandedHostEl;
    this.asideShellEl = options.asideShellEl;
    this.backend = "slice";
    this.volumes = [];
    this.chimeraxImages = [];
    this.sliceViewer = null;
    this.raycastView = null;
    this.vtkModulePromise = null;
    this.lastVtkViewMatrix = "";
    this.sharedViewMatrix = "";
    this._pendingApplyViewMatrixToVtk = false;
    this._pendingApplyViewTurnsToVtk = null;
    this.vtkCameraUserAdjusted = false;
    this.vtkSessionViewMatrix = "";
    this.chimeraxRenderedViewMatrix = "";
    this._vtkViewMatrixBeforeChimerax = "";
    this._vtkViewTurnsBeforeChimerax = [];
    this._vtkChimeraxSyncTurns = [];
    this._vtkChimeraxSyncMatrix = "";
    this._suppressVtkCameraCapture = false;
    this.vtkFocusIndex = 0;
    this.chimeraxFocusIndex = 0;
    this.chimeraxPreviewEl = null;
    this.raycastVolIndex = null;
    this.expandedBelow = false;
    this.chimeraxRendering = false;
    this.chimeraxRerenderInFlight = false;
    this.incompleteVolumeOverlay = false;
    this.expectedVolumeCount = null;
    this.chimeraxIsoLevel = null;
    this.isoPercentileSamples = null;
    this.isoSliderRange = null;
    this.onChimeraxIsoChange = options.onChimeraxIsoChange || null;
    this.onResetViewClick = options.onResetViewClick || null;
    this.canResetView = options.canResetView || null;
    this._chimeraxIsoRerenderTimer = null;
    this.onFocusChange = options.onFocusChange || null;
    this.getVolumeNavLabels = options.getVolumeNavLabels || null;
    this._volumeSliderTickLabelFontPx = null;
    this._viewportHomeParent = null;
    this._viewportHomeNext = null;
    this._displayRowHomeParent = null;
    this._displayRowHomeNext = null;
    this._controlsHomeParent = null;
    this._controlsHomeNext = null;
    this.stableBackendChrome = !!options.stableBackendChrome;

    if (this.canvasEl && global.CryoVolumeSliceCanvas) {
      this.sliceViewer = new global.CryoVolumeSliceCanvas({ canvas: this.canvasEl });
    }
    this._wireControls(options);
  }

  TrajectoryVolumeDisplay.prototype.setVtkBundleUrl = function (url) {
    VTK_BUNDLE_URL = url;
  };

  TrajectoryVolumeDisplay.prototype._wireControls = function (options) {
    var self = this;
    if (options.btnPanUp) options.btnPanUp.addEventListener("click", function () { self._pan(0, -1); });
    if (options.btnPanDown) options.btnPanDown.addEventListener("click", function () { self._pan(0, 1); });
    if (options.btnPanLeft) options.btnPanLeft.addEventListener("click", function () { self._pan(-1, 0); });
    if (options.btnPanRight) options.btnPanRight.addEventListener("click", function () { self._pan(1, 0); });
    if (options.btnZoomIn) options.btnZoomIn.addEventListener("click", function () { self._zoomSlice(1); });
    if (options.btnZoomOut) options.btnZoomOut.addEventListener("click", function () { self._zoomSlice(-1); });
    if (options.btnResetView) {
      options.btnResetView.addEventListener("click", function () {
        if (typeof self.onResetViewClick === "function") {
          self.onResetViewClick();
          return;
        }
        self.resetInteractiveView();
      });
    }
    if (this.isoSliderEl) {
      this.isoSliderEl.addEventListener("input", function () {
        self._onIsoInput();
      });
    }
    if (options.rotationLockEl && this.sliceViewer) {
      options.rotationLockEl.addEventListener("change", function () {
        if (options.rotationLockEl.checked) self.sliceViewer.setRotationLocked(true);
        else self.sliceViewer.setRotationLocked(false);
      });
    }
    if (options.sliceContrastEl && this.sliceViewer) {
      options.sliceContrastEl.addEventListener("input", function () {
        self.sliceViewer.setContrast(Number(options.sliceContrastEl.value));
      });
    }
    if (this.volumeSliderEl) {
      this.volumeSliderEl.addEventListener("input", function () {
        self._applySliderFocus();
      });
      this.volumeSliderEl.addEventListener("change", function () {
        self._applySliderFocus();
      });
    }
    if (this.volumeSliderTicksEl) {
      this.volumeSliderTicksEl.addEventListener("click", function (ev) {
        var tick = ev.target.closest("[data-vol-index]");
        if (!tick || tick.disabled) return;
        var idx = Number(tick.getAttribute("data-vol-index"));
        if (!Number.isFinite(idx)) return;
        self.setFocusIndex(idx);
      });
      if (typeof ResizeObserver !== "undefined") {
        var tickFitTimer = 0;
        var sliderWrap = this.volumeSliderEl && this.volumeSliderEl.parentElement;
        this._volumeSliderTicksResizeObserver = new ResizeObserver(function () {
          if (tickFitTimer) clearTimeout(tickFitTimer);
          tickFitTimer = setTimeout(function () {
            tickFitTimer = 0;
            self._fitVolumeSliderTickLabelFont();
          }, 40);
        });
        if (sliderWrap) {
          this._volumeSliderTicksResizeObserver.observe(sliderWrap);
        } else if (this.volumeSliderTicksEl) {
          this._volumeSliderTicksResizeObserver.observe(this.volumeSliderTicksEl);
        }
      }
    }
    if (this.btnVolPrev) {
      this.btnVolPrev.addEventListener("click", function () { self._cycleVtkFocus(-1); });
    }
    if (this.btnVolNext) {
      this.btnVolNext.addEventListener("click", function () { self._cycleVtkFocus(1); });
    }
    if (options.btnDockBelow) {
      options.btnDockBelow.addEventListener("click", function () {
        if (typeof options.onDockBelowClick === "function") {
          options.onDockBelowClick();
        } else {
          self.toggleExpandedBelow();
        }
      });
    }
    if (options.backendRadios) {
      options.backendRadios.forEach(function (radio) {
        radio.addEventListener("change", function () {
          if (!radio.checked) return;
          if (typeof options.onBackendChange === "function") {
            options.onBackendChange(radio.value);
          } else {
            self.setBackend(radio.value);
          }
        });
      });
    }
  };

  TrajectoryVolumeDisplay.prototype._ensureChimeraxPreview = function () {
    if (this.chimeraxPreviewEl || !this.viewportEl) return;
    var img = document.createElement("img");
    img.className = "cryo-vslice-chimerax-preview";
    img.id = "vslice-chimerax-preview";
    img.alt = "ChimeraX volume";
    img.hidden = true;
    this.viewportEl.appendChild(img);
    this.chimeraxPreviewEl = img;
  };

  TrajectoryVolumeDisplay.prototype._chimeraxImageSrc = function (b64) {
    if (!b64) return "";
    var s = String(b64).trim();
    if (!s) return "";
    if (s.indexOf("data:image") === 0) return s;
    return "data:image/png;base64," + s;
  };

  TrajectoryVolumeDisplay.prototype._chimeraxImageReadyAt = function (index) {
    var i = Math.floor(Number(index));
    return i >= 0 && i < this.chimeraxImages.length && !!this.chimeraxImages[i];
  };

  TrajectoryVolumeDisplay.prototype._volumeReadyAt = function (index) {
    if (this.backend === "chimerax") return this._chimeraxImageReadyAt(index);
    var i = Math.floor(Number(index));
    return i >= 0
      && i < this.volumes.length
      && !!(this.volumes[i] && this.volumes[i].volume_b64);
  };

  TrajectoryVolumeDisplay.prototype._countRenderedChimeraxImages = function () {
    var total = this._volumeNavCount();
    var n = 0;
    for (var i = 0; i < total; i++) {
      if (this._chimeraxImageReadyAt(i)) n++;
    }
    return n;
  };

  TrajectoryVolumeDisplay.prototype._nearestRenderedChimeraxIndex = function (target) {
    var total = this._volumeNavCount();
    if (total < 1) return 0;
    var t = Math.max(0, Math.min(total - 1, Math.floor(Number(target))));
    if (this._chimeraxImageReadyAt(t)) return t;
    for (var d = 1; d < total; d++) {
      if (t - d >= 0 && this._chimeraxImageReadyAt(t - d)) return t - d;
      if (t + d < total && this._chimeraxImageReadyAt(t + d)) return t + d;
    }
    return t;
  };

  TrajectoryVolumeDisplay.prototype._nearestReadyVolumeIndex = function (target) {
    var total = this._volumeNavCount();
    if (total < 1) return 0;
    var t = Math.max(0, Math.min(total - 1, Math.floor(Number(target))));
    if (this._volumeReadyAt(t)) return t;
    for (var d = 1; d < total; d++) {
      if (t - d >= 0 && this._volumeReadyAt(t - d)) return t - d;
      if (t + d < total && this._volumeReadyAt(t + d)) return t + d;
    }
    return t;
  };

  TrajectoryVolumeDisplay.prototype._snapFocusIndexToReady = function (target) {
    var total = this._volumeNavCount();
    if (total < 1) return -1;
    var t = Math.max(0, Math.min(total - 1, Math.floor(Number(target))));
    if (!this._volumeReadyAt(t)) {
      t = this.backend === "chimerax"
        ? this._nearestRenderedChimeraxIndex(t)
        : this._nearestReadyVolumeIndex(t);
    }
    return this._volumeReadyAt(t) ? t : -1;
  };

  TrajectoryVolumeDisplay.prototype._stepReadyFocus = function (delta) {
    var total = this._volumeNavCount();
    if (total < 1) return;
    var idx = this.backend === "chimerax" ? this.chimeraxFocusIndex : this.vtkFocusIndex;
    for (var attempt = 0; attempt < total; attempt++) {
      idx = (idx + delta + total) % total;
      if (this._volumeReadyAt(idx)) {
        this.setFocusIndex(idx);
        return;
      }
    }
  };

  TrajectoryVolumeDisplay.prototype.setChimeraxImageAt = function (index, b64) {
    var i = Math.floor(Number(index));
    if (!Number.isFinite(i) || i < 0) return;
    if (!b64) {
      this.clearVolumeAt(i);
      return;
    }
    var total = this.expectedVolumeCount != null && this.expectedVolumeCount > 0
      ? this.expectedVolumeCount
      : Math.max(this.chimeraxImages.length, i + 1);
    while (this.chimeraxImages.length < total) this.chimeraxImages.push(null);
    this.chimeraxImages[i] = b64;
    if (this.backend !== "chimerax") return;
    if (!this._chimeraxImageReadyAt(this.chimeraxFocusIndex)) {
      this.chimeraxFocusIndex = this._nearestRenderedChimeraxIndex(this.chimeraxFocusIndex);
    }
    this._renderChimerax();
    this._syncChrome();
    this._syncChimeraxRenderingOverlay();
  };

  TrajectoryVolumeDisplay.prototype.clearVolumeAt = function (index) {
    var i = Math.floor(Number(index));
    if (!Number.isFinite(i) || i < 0) return;
    var total = this._expectedVolumeCount();
    if (!total) {
      total = Math.max(this.chimeraxImages.length, this.volumes.length, i + 1);
    }
    while (this.chimeraxImages.length < total) this.chimeraxImages.push(null);
    while (this.volumes.length < total) this.volumes.push(null);
    if (i < this.chimeraxImages.length) this.chimeraxImages[i] = null;
    if (i < this.volumes.length) this.volumes[i] = null;
    if (this.backend === "chimerax" && !this._chimeraxImageReadyAt(this.chimeraxFocusIndex)) {
      this.chimeraxFocusIndex = this._nearestRenderedChimeraxIndex(this.chimeraxFocusIndex);
    }
    if (this.backend === "vtk" && !this._volumeReadyAt(this.vtkFocusIndex)) {
      var snapped = this._snapFocusIndexToReady(this.vtkFocusIndex);
      if (snapped >= 0) this.vtkFocusIndex = snapped;
    }
    this._renderCurrent();
    this._syncChrome();
    this._syncChimeraxRenderingOverlay();
  };

  TrajectoryVolumeDisplay.prototype.getFocusIndex = function () {
    if (this.backend === "chimerax") return this.chimeraxFocusIndex;
    if (this.backend === "vtk") return this.vtkFocusIndex;
    return 0;
  };

  TrajectoryVolumeDisplay.prototype._notifyFocusChange = function () {
    if (typeof this.onFocusChange === "function") {
      this.onFocusChange(this);
    }
  };

  TrajectoryVolumeDisplay.prototype._renderChimerax = function () {
    this._repatriateDisplayFromExpandedHost();
    if (this.chimeraxPreviewEl && this.chimeraxPreviewEl.parentElement !== this.viewportEl) {
      this.chimeraxPreviewEl = null;
    }
    this._ensureChimeraxPreview();
    this._destroyRaycast();
    if (this.sliceViewer) this.sliceViewer.setVolumes([], false);
    var hasImages = this._countRenderedChimeraxImages() > 0;
    if (this.chimeraxFocusIndex < 0) this.chimeraxFocusIndex = 0;
    var displayIndex = hasImages
      ? this._nearestRenderedChimeraxIndex(this.chimeraxFocusIndex)
      : this.chimeraxFocusIndex;
    if (this.canvasEl) this.canvasEl.hidden = true;
    if (this.vtkContainerEl) this.vtkContainerEl.hidden = true;
    if (this.chimeraxPreviewEl) {
      if (this._chimeraxImageReadyAt(displayIndex)) {
        this.chimeraxPreviewEl.src = this._chimeraxImageSrc(
          this.chimeraxImages[displayIndex]
        );
        this.chimeraxPreviewEl.hidden = false;
      } else {
        this.chimeraxPreviewEl.src = "";
        this.chimeraxPreviewEl.hidden = true;
      }
    }
    this._syncChrome();
    this._syncChimeraxRenderingOverlay();
  };

  TrajectoryVolumeDisplay.prototype.setStatus = function (msg, busy) {
    if (busy && !this._hasDisplayableVolumes()) {
      busy = false;
      msg = "";
    }
    var label = busy ? (msg || "Rendering…") : "";
    if (this.statusEl) this.statusEl.textContent = label;
    if (this.renderingOverlayEl) {
      this.renderingOverlayEl.hidden = !busy;
      this.renderingOverlayEl.setAttribute("aria-hidden", busy ? "false" : "true");
      this.renderingOverlayEl.classList.toggle("cryo-plot-rendering-overlay--show", !!busy);
      var labelEl = this.renderingOverlayEl.querySelector(".cryo-plot-rendering-overlay__label");
      if (labelEl) labelEl.textContent = label || "Rendering…";
    } else if (this.progressEl) {
      this.progressEl.hidden = !busy;
    }
  };

  TrajectoryVolumeDisplay.prototype._expectedVolumeCount = function () {
    var n = Number(this.expectedVolumeCount);
    return Number.isFinite(n) && n > 0 ? Math.floor(n) : 0;
  };

  TrajectoryVolumeDisplay.prototype._countReadyVolumes = function () {
    var total = this._volumeNavCount();
    var n = 0;
    for (var i = 0; i < total; i++) {
      if (this._volumeReadyAt(i)) n++;
    }
    return n;
  };

  TrajectoryVolumeDisplay.prototype._hasDisplayableVolumes = function () {
    if (this.backend === "chimerax") {
      return this._countRenderedChimeraxImages() > 0;
    }
    return this._countReadyVolumes() > 0;
  };

  TrajectoryVolumeDisplay.prototype.setIncompleteVolumeOverlay = function (on) {
    this.incompleteVolumeOverlay = !!on;
    this._syncChimeraxRenderingOverlay();
  };

  TrajectoryVolumeDisplay.prototype._syncChimeraxRenderingOverlay = function () {
    var anyReady = this._hasDisplayableVolumes();
    var rerenderBusy = this.backend === "chimerax"
      && !!this.chimeraxRendering
      && !!this.chimeraxRerenderInFlight
      && anyReady;
    var busy = rerenderBusy;
    this.setStatus(busy ? "Re-rendering…" : "", busy);
  };

  TrajectoryVolumeDisplay.prototype.setChimeraxRendering = function (on, opts) {
    opts = opts || {};
    this.chimeraxRendering = !!on;
    if (!on) {
      this.chimeraxRerenderInFlight = false;
    } else if (opts.rerender) {
      this.chimeraxRerenderInFlight = true;
    } else if (opts.rerender === false) {
      this.chimeraxRerenderInFlight = false;
    }
    if (this.backend === "chimerax") {
      this._syncChrome();
      if (this._hasDisplayableVolumes()) {
        this._renderChimerax();
      }
      this._syncChimeraxRenderingOverlay();
    } else if (!on) {
      this.setStatus("", false);
    }
  };

  TrajectoryVolumeDisplay.prototype.setBackend = function (backend) {
    backend = String(backend || "slice").toLowerCase();
    if (backend !== "vtk" && backend !== "slice" && backend !== "chimerax") backend = "slice";
    if (this.backend === "chimerax" && backend === "vtk") {
      this.raycastVolIndex = null;
      this._pendingApplyViewTurnsToVtk = null;
      this._pendingApplyViewMatrixToVtk = false;
      if (this._vtkViewMatrixBeforeChimerax) {
        this.setSharedViewMatrix(this._vtkViewMatrixBeforeChimerax);
        this._pendingApplyViewMatrixToVtk = true;
      } else if (this._vtkViewTurnsBeforeChimerax && this._vtkViewTurnsBeforeChimerax.length) {
        this._pendingApplyViewTurnsToVtk = this._vtkViewTurnsBeforeChimerax.slice();
      }
    } else if (this.backend === "vtk" && backend === "chimerax" && this.volumes.length > 0) {
      if (this.raycastView) {
        this.getChimeraxViewMatrix();
        this._vtkViewMatrixBeforeChimerax = this.getSharedViewMatrix();
        if (typeof this.raycastView.getChimeraxViewTurns === "function") {
          var vtkTurns = this.raycastView.getChimeraxViewTurns();
          if (vtkTurns && vtkTurns.length) {
            this._vtkViewTurnsBeforeChimerax = vtkTurns.slice();
          }
        }
      }
    }
    if (backend === "chimerax" && this.expandedBelow) {
      this.setExpandedBelow(false);
    }
    this.backend = backend;
    if (backend !== "chimerax") {
      this.chimeraxRendering = false;
      this.chimeraxRerenderInFlight = false;
    }
    if (this.raycastView && typeof this.raycastView.setInteractionEnabled === "function") {
      this.raycastView.setInteractionEnabled(backend === "vtk");
    }
    if (backend === "chimerax") {
      this._repatriateDisplayFromExpandedHost();
    }
    this._syncChrome();
    this._renderCurrent();
  };

  TrajectoryVolumeDisplay.prototype.setChimeraxRenderedViewMatrix = function (vm) {
    var text = vm && isValidChimeraxViewMatrixText(vm) ? String(vm).trim() : "";
    this.chimeraxRenderedViewMatrix = text;
  };

  function isValidChimeraxViewMatrixText(vm) {
    if (!vm) return false;
    var raw = String(vm).trim();
    if (!raw) return false;
    if (raw.toLowerCase().indexOf("camera") === 0) raw = raw.slice(6).trim();
    var nums = raw.match(/[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?/g);
    if (!nums || nums.length < 12) return false;
    for (var ni = 0; ni < 12; ni++) {
      if (!isFinite(Number(nums[ni]))) return false;
    }
    return true;
  }

  TrajectoryVolumeDisplay.prototype.setSharedViewMatrix = function (vm) {
    var text = vm ? String(vm).trim() : "";
    this.sharedViewMatrix = text;
    this.lastVtkViewMatrix = text;
  };

  TrajectoryVolumeDisplay.prototype.getSharedViewMatrix = function () {
    return this.sharedViewMatrix || this.lastVtkViewMatrix || "";
  };

  TrajectoryVolumeDisplay.prototype.getChimeraxViewMatrix = function () {
    if (this.raycastView) {
      if (typeof this.raycastView.getChimeraxViewMatrixCamera === "function") {
        var live = this.raycastView.getChimeraxViewMatrixCamera();
        if (live) {
          this.setSharedViewMatrix(live);
          return live;
        }
      }
      var renderer = this.raycastView.renderer;
      var cam = renderer && renderer.getActiveCamera && renderer.getActiveCamera();
      if (cam && cam.getViewMatrix) {
        var mat = new Float64Array(16);
        cam.getViewMatrix(mat);
        var out = [];
        for (var row = 0; row < 3; row++) {
          for (var col = 0; col < 4; col++) {
            out.push(mat[col * 4 + row]);
          }
        }
        if (out.every(function (v) { return isFinite(v); })) {
          this.setSharedViewMatrix(out.map(function (v) {
            return Number(v).toPrecision(8);
          }).join(","));
        }
      }
    }
    return this.getSharedViewMatrix();
  };

  TrajectoryVolumeDisplay.prototype.resetVtkNavigationCamera = function () {
    this.vtkSessionViewMatrix = "";
    this.vtkCameraUserAdjusted = false;
    this.sharedViewMatrix = "";
    this.lastVtkViewMatrix = "";
    this.raycastVolIndex = null;
  };

  TrajectoryVolumeDisplay.prototype._captureVtkSessionViewMatrix = function (view) {
    if (view && typeof view.getChimeraxViewMatrixCamera === "function") {
      var live = view.getChimeraxViewMatrixCamera();
      if (live) {
        this.setSharedViewMatrix(live);
        this.vtkSessionViewMatrix = live;
        return;
      }
    }
    this.getChimeraxViewMatrix();
    var vm = this.getSharedViewMatrix();
    if (vm) this.vtkSessionViewMatrix = vm;
  };

  TrajectoryVolumeDisplay.prototype._applyVtkNavigationViewMatrix = function (view) {
    if (!view) return;
    var vm = this.vtkCameraUserAdjusted
      ? this.getSharedViewMatrix()
      : (this.vtkSessionViewMatrix || this.getSharedViewMatrix());
    if (!vm || typeof view.setChimeraxViewMatrixCamera !== "function") return;
    view.setChimeraxViewMatrixCamera(vm);
  };

  TrajectoryVolumeDisplay.prototype._markVtkCameraUserAdjusted = function () {
    this.vtkCameraUserAdjusted = true;
    this.getChimeraxViewMatrix();
    var vm = this.getSharedViewMatrix();
    if (vm) this.vtkSessionViewMatrix = vm;
  };

  TrajectoryVolumeDisplay.prototype._applyPendingViewTurnsToVtk = function (view) {
    var turns = this._pendingApplyViewTurnsToVtk;
    this._pendingApplyViewTurnsToVtk = null;
    if (!view || !turns || !turns.length) return;
    if (typeof view.applyChimeraxViewTurns === "function") {
      view.applyChimeraxViewTurns(turns);
    }
  };

  TrajectoryVolumeDisplay.prototype._applySharedViewMatrixToVtk = function (view) {
    var vm = this.sharedViewMatrix || this.lastVtkViewMatrix;
    if (!view || !vm) return;
    if (typeof view.setChimeraxViewMatrixCamera === "function") {
      view.setChimeraxViewMatrixCamera(vm);
    }
    this._pendingApplyViewMatrixToVtk = false;
  };

  TrajectoryVolumeDisplay.prototype.loadPayload = function (payload, opts) {
    opts = opts || {};
    this.volumes = (payload && payload.volumes) ? payload.volumes.slice() : [];
    this.chimeraxImages = (payload && payload.images) ? payload.images.slice() : [];
    if (payload && payload.expectedVolumeCount != null) {
      this.expectedVolumeCount = Number(payload.expectedVolumeCount);
    } else if (!this.volumes.length && !(payload && payload.images && payload.images.length)) {
      this.expectedVolumeCount = null;
    }
    var expected = this.expectedVolumeCount;
    if (expected != null && expected > 0) {
      while (this.volumes.length < expected) this.volumes.push(null);
      while (this.chimeraxImages.length < expected) this.chimeraxImages.push(null);
    }
    if (this.vtkFocusIndex >= this.volumes.length) this.vtkFocusIndex = 0;
    if (this.vtkFocusIndex < 0) this.vtkFocusIndex = 0;
    if (this.chimeraxFocusIndex >= this.chimeraxImages.length) this.chimeraxFocusIndex = 0;
    if (this.chimeraxFocusIndex < 0) this.chimeraxFocusIndex = 0;
    this._refreshIsoSamplesFromVolumes();
    if (!opts.deferRender) this._renderCurrent();
  };

  TrajectoryVolumeDisplay.prototype._refreshIsoSamplesFromVolumes = function () {
    if (!global.CryoVolume3dUtils) return;
    var vol = null;
    for (var i = 0; i < this.volumes.length; i++) {
      if (this.volumes[i] && this.volumes[i].volume_b64) {
        vol = this.volumes[i];
        break;
      }
    }
    if (!vol) return;
    try {
      var U = global.CryoVolume3dUtils;
      var values = U.decodeFloat32Volume(vol.volume_b64, vol.D);
      this.isoPercentileSamples = U.volumePercentileSamples(values);
      this.isoSliderRange = U.isoSliderDataRange(this.isoPercentileSamples);
      if (this.chimeraxIsoLevel == null) {
        this.chimeraxIsoLevel = U.suggestIsoDataValue(this.isoPercentileSamples);
      }
      this._syncIsoSliderFromState();
    } catch (err) {
      // Keep prior iso state if volume decode fails.
    }
  };

  TrajectoryVolumeDisplay.prototype.setChimeraxIsoRange = function (min, max, level) {
    min = Number(min);
    max = Number(max);
    if (isFinite(min) && isFinite(max) && max > min) {
      this.isoSliderRange = { min: min, max: max };
    }
    if (level != null && isFinite(Number(level))) {
      this.chimeraxIsoLevel = Number(level);
    }
    this._syncIsoSliderFromState();
  };

  TrajectoryVolumeDisplay.prototype.getChimeraxIsoLevel = function () {
    if (this.chimeraxIsoLevel != null && isFinite(this.chimeraxIsoLevel)) {
      return this.chimeraxIsoLevel;
    }
    if (this.raycastView && typeof this.raycastView.getIsoLevel === "function") {
      return this.raycastView.getIsoLevel();
    }
    return null;
  };

  TrajectoryVolumeDisplay.prototype.hasInteractiveVolumes = function () {
    return this.volumes.some(function (v) { return v && v.volume_b64; });
  };

  TrajectoryVolumeDisplay.prototype.hasChimeraxImages = function () {
    return this._countRenderedChimeraxImages() > 0;
  };

  TrajectoryVolumeDisplay.prototype.setFocusIndex = function (index) {
    var i = Number(index);
    if (!Number.isFinite(i)) return;
    i = this._snapFocusIndexToReady(i);
    if (i < 0) return;
    if (this.backend === "chimerax") {
      var cxTotal = this._volumeNavCount();
      if (cxTotal < 1) return;
      if (this.chimeraxFocusIndex === i && this._chimeraxImageReadyAt(i)) {
        this._syncVolumeNavChrome();
        return;
      }
      this.chimeraxFocusIndex = i;
      this._renderChimerax();
      this._notifyFocusChange();
      this._syncVolumeNavChrome();
      return;
    }
    if (this.backend !== "vtk") return;
    var n = this._volumeNavCount();
    if (n < 1) return;
    if (this.vtkFocusIndex === i && this.raycastVolIndex === i) {
      this._syncVolumeNavChrome();
      return;
    }
    this.vtkFocusIndex = i;
    this.raycastVolIndex = null;
    this._syncVolumeNavChrome();
    this._renderVtk();
    this._notifyFocusChange();
  };

  TrajectoryVolumeDisplay.prototype._applySliderFocus = function () {
    if (!this.volumeSliderEl) return;
    var snapped = this._snapFocusIndexToReady(Number(this.volumeSliderEl.value));
    if (snapped < 0) return;
    if (String(this.volumeSliderEl.value) !== String(snapped)) {
      this.volumeSliderEl.value = String(snapped);
    }
    this.setFocusIndex(snapped);
  };

  TrajectoryVolumeDisplay.prototype._volumeNavCount = function () {
    var isVtk = this.backend === "vtk";
    var isChimeraX = this.backend === "chimerax";
    if (isChimeraX) {
      if (this.expectedVolumeCount != null && this.expectedVolumeCount > 0) {
        return this.expectedVolumeCount;
      }
      return this.chimeraxImages.length;
    }
    if (isVtk) {
      if (this.expectedVolumeCount != null && this.expectedVolumeCount > 0) {
        return this.expectedVolumeCount;
      }
      return this.volumes.length;
    }
    return 0;
  };

  TrajectoryVolumeDisplay.prototype._volumeNavLabels = function () {
    var total = this._volumeNavCount();
    if (typeof this.getVolumeNavLabels === "function") {
      var custom = this.getVolumeNavLabels(total);
      if (Array.isArray(custom) && custom.length) return custom;
    }
    var labels = [];
    for (var i = 0; i < total; i++) labels.push(String(i + 1));
    return labels;
  };

  TrajectoryVolumeDisplay.prototype._volumeSliderTickLabelMeasurer = function () {
    if (!this._volumeSliderTickLabelMeasurerEl) {
      var el = document.createElement("span");
      el.className = "cryo-vslice-volume-slider-tick-label";
      el.setAttribute("aria-hidden", "true");
      el.style.position = "absolute";
      el.style.left = "-10000px";
      el.style.top = "0";
      el.style.visibility = "hidden";
      el.style.maxWidth = "none";
      el.style.whiteSpace = "pre-line";
      el.style.pointerEvents = "none";
      document.body.appendChild(el);
      this._volumeSliderTickLabelMeasurerEl = el;
    }
    return this._volumeSliderTickLabelMeasurerEl;
  };

  TrajectoryVolumeDisplay.prototype._measureVolumeSliderTickLabelWidth = function (
    measurer,
    text,
    fontSizePx,
    fontWeight,
    fontFamily
  ) {
    measurer.style.fontSize = String(fontSizePx) + "px";
    measurer.style.fontWeight = fontWeight || "600";
    measurer.style.fontFamily = fontFamily || '"Barlow", ui-sans-serif, system-ui, sans-serif';
    measurer.style.whiteSpace = "pre-line";
    var raw = String(text || "");
    var lines = raw.split("\n");
    var maxW = 0;
    for (var li = 0; li < lines.length; li++) {
      measurer.textContent = lines[li];
      maxW = Math.max(maxW, measurer.getBoundingClientRect().width);
    }
    return maxW;
  };

  TrajectoryVolumeDisplay.prototype._volumeSliderTickLabelWidthsAtSize = function (
    labels,
    fontSizePx
  ) {
    var measurer = this._volumeSliderTickLabelMeasurer();
    var widths = [];
    for (var i = 0; i < labels.length; i++) {
      var cs = window.getComputedStyle(labels[i]);
      widths.push(this._measureVolumeSliderTickLabelWidth(
        measurer,
        labels[i].textContent,
        fontSizePx,
        cs.fontWeight,
        cs.fontFamily
      ));
    }
    return widths;
  };

  TrajectoryVolumeDisplay.prototype._volumeSliderThumbPx = function () {
    var wrap = this.volumeSliderEl && this.volumeSliderEl.parentElement;
    if (!wrap) return 16;
    var raw = window.getComputedStyle(wrap).getPropertyValue("--cryo-vol-slider-thumb").trim();
    if (raw) {
      var parsed = parseFloat(raw);
      if (isFinite(parsed) && parsed > 0) return parsed;
    }
    return 16;
  };

  TrajectoryVolumeDisplay.prototype._volumeSliderTrackWidthPx = function () {
    var wrap = this.volumeSliderEl && this.volumeSliderEl.parentElement;
    if (!wrap) {
      return this.volumeSliderTicksEl ? this.volumeSliderTicksEl.clientWidth : 0;
    }
    return wrap.clientWidth;
  };

  TrajectoryVolumeDisplay.prototype._volumeSliderTickCentersPx = function (total) {
    total = Number(total);
    if (!isFinite(total) || total < 1) return [];
    var trackW = this._volumeSliderTrackWidthPx();
    if (trackW <= 0) return [];
    var thumb = this._volumeSliderThumbPx();
    var travel = Math.max(0, trackW - thumb);
    var inset = thumb * 0.5;
    if (total === 1) return [inset];
    var centers = [];
    for (var i = 0; i < total; i++) {
      centers.push(inset + (travel * i) / (total - 1));
    }
    return centers;
  };

  TrajectoryVolumeDisplay.prototype._positionVolumeSliderTicks = function (total) {
    if (!this.volumeSliderTicksEl || total < 1) return;
    var ticks = this.volumeSliderTicksEl.querySelectorAll(".cryo-vslice-volume-slider-tick");
    for (var i = 0; i < ticks.length; i++) {
      var frac = total > 1 ? i / (total - 1) : 0;
      ticks[i].style.setProperty("--tick-pos", String(frac));
    }
  };

  TrajectoryVolumeDisplay.prototype._volumeSliderTickLabelMaxWidthAt = function (
    index,
    centersPx
  ) {
    if (!centersPx.length) return 0;
    var leftHalf = index > 0
      ? (centersPx[index] - centersPx[index - 1]) * 0.5
      : (centersPx.length > 1 ? (centersPx[1] - centersPx[0]) * 0.5 : centersPx[0]);
    var rightHalf = index < centersPx.length - 1
      ? (centersPx[index + 1] - centersPx[index]) * 0.5
      : leftHalf;
    return Math.max(0, Math.min(leftHalf, rightHalf) * 2);
  };

  TrajectoryVolumeDisplay.prototype._volumeSliderTickLabelsFitAtSize = function (
    widths,
    centersPx
  ) {
    if (!widths.length) return true;
    var lMax = 0;
    for (var i = 0; i < widths.length; i++) {
      if (widths[i] > lMax) lMax = widths[i];
    }
    var minGap = 0.2 * lMax;
    for (var j = 0; j < widths.length; j++) {
      var maxW = this._volumeSliderTickLabelMaxWidthAt(j, centersPx);
      if (widths[j] > maxW + 0.5) return false;
    }
    if (widths.length < 2) return true;
    for (var k = 0; k < widths.length - 1; k++) {
      var gap = (centersPx[k + 1] - centersPx[k]) - ((widths[k] + widths[k + 1]) * 0.5);
      if (gap + 0.5 < minGap) return false;
    }
    return true;
  };

  TrajectoryVolumeDisplay.prototype._markVolumeSliderTickLabelsPending = function () {
    if (!this.volumeSliderTicksEl) return;
    this.volumeSliderTicksEl.classList.remove("cryo-vslice-volume-slider-ticks--ready");
    this.volumeSliderTicksEl.classList.add("cryo-vslice-volume-slider-ticks--pending-fit");
  };

  TrajectoryVolumeDisplay.prototype._applyVolumeSliderTickLabelFont = function (
    fontSizePx
  ) {
    if (!this.volumeSliderTicksEl || fontSizePx == null) return;
    var labels = this.volumeSliderTicksEl.querySelectorAll(
      ".cryo-vslice-volume-slider-tick-label"
    );
    var size = String(fontSizePx) + "px";
    var centersPx = this._volumeSliderTickCentersPx(labels.length);
    for (var li = 0; li < labels.length; li++) {
      labels[li].style.fontSize = size;
      if (centersPx.length) {
        labels[li].style.maxWidth =
          String(this._volumeSliderTickLabelMaxWidthAt(li, centersPx)) + "px";
      }
    }
    this._volumeSliderTickLabelFontPx = fontSizePx;
    var ticksEl = this.volumeSliderTicksEl;
    requestAnimationFrame(function () {
      requestAnimationFrame(function () {
        if (!ticksEl) return;
        ticksEl.classList.remove("cryo-vslice-volume-slider-ticks--pending-fit");
        ticksEl.classList.add("cryo-vslice-volume-slider-ticks--ready");
      });
    });
  };

  TrajectoryVolumeDisplay.prototype._fitVolumeSliderTickLabelFont = function () {
    if (!this.volumeSliderTicksEl) return;
    var labels = this.volumeSliderTicksEl.querySelectorAll(
      ".cryo-vslice-volume-slider-tick-label"
    );
    if (!labels.length) return;
    this._markVolumeSliderTickLabelsPending();

    var total = labels.length;
    this._positionVolumeSliderTicks(total);
    var centersPx = this._volumeSliderTickCentersPx(total);
    if (!centersPx.length) {
      this._scheduleVolumeSliderTickLabelFit();
      return;
    }

    var lo = 6;
    var multiline = false;
    for (var mi = 0; mi < labels.length; mi++) {
      if (String(labels[mi].textContent || "").indexOf("\n") >= 0) {
        multiline = true;
        break;
      }
    }
    var hi = multiline ? 11 : 40;
    var best = lo;
    while (lo <= hi) {
      var mid = Math.ceil((lo + hi) * 0.5);
      var widths = this._volumeSliderTickLabelWidthsAtSize(labels, mid);
      if (this._volumeSliderTickLabelsFitAtSize(widths, centersPx)) {
        best = mid;
        lo = mid + 1;
      } else {
        hi = mid - 1;
      }
    }

    this._applyVolumeSliderTickLabelFont(best);
  };

  TrajectoryVolumeDisplay.prototype._scheduleVolumeSliderTickLabelFit = function () {
    var self = this;
    this._markVolumeSliderTickLabelsPending();
    if (this._volumeSliderTickFitRaf) return;
    this._volumeSliderTickFitRaf = requestAnimationFrame(function () {
      requestAnimationFrame(function () {
        self._volumeSliderTickFitRaf = 0;
        self._fitVolumeSliderTickLabelFont();
      });
    });
  };

  TrajectoryVolumeDisplay.prototype._syncVolumeSliderTicks = function (focus, total) {
    if (!this.volumeSliderTicksEl) return;
    var labels = this._volumeNavLabels();
    if (total < 1) {
      this.volumeSliderTicksEl.innerHTML = "";
      this.volumeSliderTicksEl.classList.remove(
        "cryo-vslice-volume-slider-ticks--ready",
        "cryo-vslice-volume-slider-ticks--pending-fit"
      );
      return;
    }
    var needsLabelReflow = this.volumeSliderTicksEl.childElementCount !== total;
    var labelTextChanged = false;
    if (needsLabelReflow) {
      this._markVolumeSliderTickLabelsPending();
      this.volumeSliderTicksEl.innerHTML = "";
      for (var ti = 0; ti < total; ti++) {
        var tick = document.createElement("button");
        tick.type = "button";
        tick.className = "cryo-vslice-volume-slider-tick";
        tick.setAttribute("data-vol-index", String(ti));
        tick.setAttribute("aria-label", "Volume " + String(labels[ti] || ti + 1));
        var mark = document.createElement("span");
        mark.className = "cryo-vslice-volume-slider-tick-mark";
        mark.setAttribute("aria-hidden", "true");
        tick.appendChild(mark);
        var lbl = document.createElement("span");
        lbl.className = "cryo-vslice-volume-slider-tick-label";
        lbl.textContent = String(labels[ti] != null ? labels[ti] : ti + 1);
        if (this._volumeSliderTickLabelFontPx != null) {
          lbl.style.fontSize = String(this._volumeSliderTickLabelFontPx) + "px";
        }
        tick.appendChild(lbl);
        this.volumeSliderTicksEl.appendChild(tick);
      }
      labelTextChanged = true;
    } else {
      var tickBtns = this.volumeSliderTicksEl.querySelectorAll(".cryo-vslice-volume-slider-tick");
      for (var uj = 0; uj < tickBtns.length; uj++) {
        var lblEl = tickBtns[uj].querySelector(".cryo-vslice-volume-slider-tick-label");
        if (!lblEl) continue;
        var nextText = String(labels[uj] != null ? labels[uj] : uj + 1);
        if (lblEl.textContent !== nextText) {
          labelTextChanged = true;
          this._markVolumeSliderTickLabelsPending();
          lblEl.textContent = nextText;
        }
      }
    }
    var ticks = this.volumeSliderTicksEl.querySelectorAll(".cryo-vslice-volume-slider-tick");
    var isChimeraX = this.backend === "chimerax";
    var pipelineBusy = isChimeraX && this.chimeraxRendering;
    var rerenderBusy = isChimeraX && this.chimeraxRerenderInFlight;
    for (var tj = 0; tj < ticks.length; tj++) {
      var tickReady = isChimeraX
        ? this._chimeraxImageReadyAt(tj)
        : this._volumeReadyAt(tj);
      ticks[tj].classList.toggle("cryo-vslice-volume-slider-tick--active", tj === focus);
      var isPending = !!((pipelineBusy && !tickReady) || rerenderBusy);
      ticks[tj].classList.toggle("cryo-vslice-volume-slider-tick--pending", isPending);
      ticks[tj].classList.toggle(
        "cryo-vslice-volume-slider-tick--inactive",
        !tickReady && !isPending
      );
      ticks[tj].disabled = focus < 0 || !tickReady;
    }
    this._positionVolumeSliderTicks(total);
    if (labelTextChanged || this._volumeSliderTickLabelFontPx == null) {
      this._scheduleVolumeSliderTickLabelFit();
    }
  };

  TrajectoryVolumeDisplay.prototype._syncVolumeNavChrome = function () {
    var isVtk = this.backend === "vtk";
    var isChimeraX = this.backend === "chimerax";
    var total = this._volumeNavCount();
    var readyCount = this._countReadyVolumes();
    var multi = total > 1;
    var active = multi;
    if (isChimeraX && this.chimeraxRendering && readyCount < 1) {
      active = false;
    }
    if (this.volumeNavEl) {
      this.volumeNavEl.classList.toggle("cryo-vslice-volume-nav--inactive", !active);
      this.volumeNavEl.setAttribute("aria-disabled", active ? "false" : "true");
      this.volumeNavEl.hidden = !(isVtk || isChimeraX) || total < 1;
    }
    var focus = isChimeraX ? this.chimeraxFocusIndex : this.vtkFocusIndex;
    if (active) {
      var snapped = this._snapFocusIndexToReady(focus);
      if (snapped >= 0 && snapped !== focus) {
        focus = snapped;
        if (isChimeraX) this.chimeraxFocusIndex = snapped;
        else this.vtkFocusIndex = snapped;
      }
    }
    var navStepEnabled = readyCount > 1;
    if (this.btnVolPrev) this.btnVolPrev.disabled = !navStepEnabled;
    if (this.btnVolNext) this.btnVolNext.disabled = !navStepEnabled;
    if (this.volumeSliderEl) {
      this.volumeSliderEl.disabled = !active;
      this.volumeSliderEl.min = "0";
      this.volumeSliderEl.max = String(Math.max(0, total - 1));
      this.volumeSliderEl.step = "1";
      if (active && focus >= 0) {
        var sliderIdx = this._snapFocusIndexToReady(focus);
        if (sliderIdx < 0) sliderIdx = Math.max(0, Math.min(total - 1, focus));
        if (String(this.volumeSliderEl.value) !== String(sliderIdx)) {
          this.volumeSliderEl.value = String(sliderIdx);
        }
        var navLabels = this._volumeNavLabels();
        var navText = navLabels[sliderIdx] != null
          ? String(navLabels[sliderIdx])
          : String(sliderIdx + 1);
        this.volumeSliderEl.setAttribute("aria-valuetext", navText);
        focus = sliderIdx;
      } else {
        this.volumeSliderEl.setAttribute("aria-valuetext", "—");
      }
    }
    this._syncVolumeSliderTicks(active ? focus : -1, total);
  };

  TrajectoryVolumeDisplay.prototype.resetInteractiveView = function () {
    if (this.backend === "vtk" && this.raycastView) {
      this._suppressVtkCameraCapture = true;
      this.raycastView.resetCamera();
      this._suppressVtkCameraCapture = false;
      this.vtkCameraUserAdjusted = false;
      this.setSharedViewMatrix("");
      return;
    }
    if (this.backend === "slice" && this.sliceViewer) {
      this.sliceViewer.resetView();
    }
  };

  TrajectoryVolumeDisplay.prototype.canResetInteractiveView = function () {
    var hasVol = this.volumes.some(function (v) { return v && v.volume_b64; });
    if (this.backend === "vtk") {
      return !!(this.raycastView && hasVol);
    }
    if (this.backend === "slice") {
      return !!(this.sliceViewer && hasVol);
    }
    return false;
  };

  TrajectoryVolumeDisplay.prototype._resetViewEnabled = function () {
    if (typeof this.canResetView === "function") {
      return !!this.canResetView();
    }
    return this.canResetInteractiveView();
  };

  TrajectoryVolumeDisplay.prototype._setLayoutStablePanelVisible = function (el, visible) {
    if (!el) return;
    if (!this.stableBackendChrome) {
      el.classList.remove("cryo-traj-vol-backend-panel--inactive");
      el.classList.remove("cryo-traj-vol-backend-panel--reserved");
      el.hidden = !visible;
      return;
    }
    el.removeAttribute("hidden");
    el.classList.remove("cryo-traj-vol-backend-panel--reserved");
    el.classList.toggle("cryo-traj-vol-backend-panel--inactive", !visible);
    el.setAttribute("aria-hidden", visible ? "false" : "true");
  };

  TrajectoryVolumeDisplay.prototype._setLayoutStablePanelReserved = function (el, visible) {
    if (!el) return;
    if (!this.stableBackendChrome) {
      el.classList.remove("cryo-traj-vol-backend-panel--inactive");
      el.classList.remove("cryo-traj-vol-backend-panel--reserved");
      el.hidden = !visible;
      return;
    }
    el.removeAttribute("hidden");
    el.classList.remove("cryo-traj-vol-backend-panel--inactive");
    el.classList.toggle("cryo-traj-vol-backend-panel--reserved", !visible);
    el.setAttribute("aria-hidden", visible ? "false" : "true");
  };

  TrajectoryVolumeDisplay.prototype.syncResetViewButton = function () {
    if (!this.btnResetView) return;
    var isChimeraX = this.backend === "chimerax";
    var enabled = this._resetViewEnabled();
    this.btnResetView.hidden = false;
    this.btnResetView.textContent = "Reset view";
    this.btnResetView.setAttribute("aria-label", "Reset volume view");
    this.btnResetView.disabled = !enabled;
    this.btnResetView.title = enabled
      ? (isChimeraX
        ? "Clear rotations and re-render with the default view."
        : "Reset pan, zoom, and camera.")
      : (isChimeraX
        ? "Render ChimeraX images before resetting the view."
        : "Load volumes before resetting the view.");
  };

  TrajectoryVolumeDisplay.prototype._syncChrome = function () {
    var isChimeraX = this.backend === "chimerax";
    var isVtk = this.backend === "vtk";
    var isSlice = this.backend === "slice";
    var hasVol = this.volumes.some(function (v) { return v && v.volume_b64; });
    var hasCxImages = this._countRenderedChimeraxImages() > 0;
    var showChimeraxPanel = isChimeraX;
    var showInteractive = !isChimeraX && hasVol;
    if (this.asideShellEl) {
      this.asideShellEl.hidden = isChimeraX ? !showChimeraxPanel : false;
    }
    if (this.displayRowEl) {
      this.displayRowEl.hidden = !(showInteractive || showChimeraxPanel);
    }
    if (this.viewportEl && isChimeraX) {
      this.viewportEl.hidden = !hasCxImages;
    } else if (this.viewportEl) {
      this.viewportEl.hidden = !hasVol;
    }
    var showControlsDock = (showInteractive || showChimeraxPanel)
      && (hasVol || hasCxImages);
    if (this.controlsDockEl) {
      this.controlsDockEl.hidden = !showControlsDock;
    }
    if (this.vtkSliceControlsEl) {
      if (this.stableBackendChrome && showControlsDock) {
        this._setLayoutStablePanelVisible(this.vtkSliceControlsEl, true);
      } else {
        this.vtkSliceControlsEl.hidden = isChimeraX;
        this.vtkSliceControlsEl.classList.remove("cryo-traj-vol-backend-panel--inactive");
      }
    }
    if (this.chimeraxViewControlsEl) {
      if (this.stableBackendChrome && showControlsDock) {
        this._setLayoutStablePanelVisible(this.chimeraxViewControlsEl, isChimeraX);
      } else {
        this.chimeraxViewControlsEl.hidden = !isChimeraX;
        this.chimeraxViewControlsEl.classList.remove("cryo-traj-vol-backend-panel--inactive");
      }
    }
    if (this.padColumnEl) {
      var showPad = !isChimeraX && hasVol;
      if (this.stableBackendChrome && hasVol) {
        this._setLayoutStablePanelReserved(this.padColumnEl, showPad);
      } else {
        this.padColumnEl.hidden = !showPad;
        this.padColumnEl.classList.remove("cryo-traj-vol-backend-panel--inactive");
        this.padColumnEl.classList.remove("cryo-traj-vol-backend-panel--reserved");
      }
    }
    if (this.isoControlsEl) {
      this.isoControlsEl.hidden = !(isVtk && hasVol) && !(isChimeraX && hasCxImages);
    }
    if (this.sliceControlsRowEl) {
      this.sliceControlsRowEl.hidden = !isSlice || !hasVol;
    }
    if (this.canvasEl) this.canvasEl.hidden = isVtk || isChimeraX;
    if (this.vtkContainerEl) this.vtkContainerEl.hidden = !isVtk;
    if (this.chimeraxPreviewEl) {
      var showCxPreview = isChimeraX && hasCxImages
        && this._chimeraxImageReadyAt(this._nearestRenderedChimeraxIndex(this.chimeraxFocusIndex));
      this.chimeraxPreviewEl.hidden = !showCxPreview;
    }
    if (this.rotationLockToolbarEl) {
      this.rotationLockToolbarEl.hidden = isVtk || isChimeraX
        || !(this.sliceViewer && this.sliceViewer.layers && this.sliceViewer.layers.length > 1);
    }
    this.syncResetViewButton();
    this._syncVolumeNavChrome();
    if (isChimeraX || isVtk) this._syncChimeraxRenderingOverlay();
    if (isVtk && this.raycastView) this._scheduleVtkResize();
  };

  TrajectoryVolumeDisplay.prototype._renderCurrent = function () {
    if (!this._hasDisplayableVolumes()) {
      if (this.chimeraxPreviewEl) {
        this.chimeraxPreviewEl.hidden = true;
        this.chimeraxPreviewEl.src = "";
      }
      this.setStatus("", false);
      this._syncChrome();
      return;
    }
    if (this.backend === "chimerax") {
      this._renderChimerax();
      return;
    }
    if (this.chimeraxPreviewEl) {
      this.chimeraxPreviewEl.hidden = true;
      this.chimeraxPreviewEl.src = "";
    }
    this._renderInteractive();
  };

  TrajectoryVolumeDisplay.prototype._stepChimeraxFocus = function (delta) {
    this._stepReadyFocus(delta);
  };

  TrajectoryVolumeDisplay.prototype._cycleVtkFocus = function (delta) {
    if (this.backend === "chimerax" || this.backend === "vtk") {
      this._stepReadyFocus(delta);
    }
  };

  TrajectoryVolumeDisplay.prototype._renderInteractive = function () {
    var self = this;
    if (!this._hasDisplayableVolumes()) {
      if (this.viewportEl) this.viewportEl.hidden = true;
      if (this.controlsDockEl) this.controlsDockEl.hidden = true;
      this._destroyRaycast();
      if (this.sliceViewer) this.sliceViewer.setVolumes([], false);
      this.setStatus("", false);
      this._syncChrome();
      return;
    }
    if (!this.volumes.length) {
      if (this.viewportEl) this.viewportEl.hidden = true;
      if (this.controlsDockEl) this.controlsDockEl.hidden = true;
      this._destroyRaycast();
      if (this.sliceViewer) this.sliceViewer.setVolumes([], false);
      this._syncChrome();
      return;
    }
    if (this.backend === "vtk") {
      if (this.sliceViewer) this.sliceViewer.setVolumes([], false);
      this._renderVtk();
      return;
    }
    this._destroyRaycast();
    var layers = [];
    for (var i = 0; i < this.volumes.length; i++) {
      var v = this.volumes[i];
      if (!v || !v.volume_b64) continue;
      layers.push({
        id: "traj:" + i,
        b64: v.volume_b64,
        d: v.D,
        label: String(i + 1),
        manual: false
      });
    }
    if (this.viewportEl) this.viewportEl.hidden = false;
    this._syncChrome();
    if (this.sliceViewer) {
      this.sliceViewer.setVolumes(layers, true);
    }
  };

  TrajectoryVolumeDisplay.prototype._scheduleVtkResize = function () {
    var self = this;
    if (!this.raycastView || typeof this.raycastView.resize !== "function") return;
    requestAnimationFrame(function () {
      requestAnimationFrame(function () {
        if (self.raycastView) self.raycastView.resize();
      });
    });
  };

  TrajectoryVolumeDisplay.prototype._vtkResizeTarget = function () {
    if (this.expandedBelow && this.expandedHostEl) return this.expandedHostEl;
    if (this.displayRowEl) return this.displayRowEl;
    return this.viewportEl;
  };

  TrajectoryVolumeDisplay.prototype._ensureVtkModule = function () {
    if (global.CryoVolumeRaycastView) return Promise.resolve();
    if (!VTK_BUNDLE_URL) return Promise.reject(new Error("VTK bundle URL not set."));
    if (!this.vtkModulePromise) {
      var self = this;
      this.vtkModulePromise = new Promise(function (resolve, reject) {
        var script = document.createElement("script");
        script.src = VTK_BUNDLE_URL;
        script.async = true;
        script.onload = function () { resolve(); };
        script.onerror = function () { reject(new Error("Failed to load VTK bundle.")); };
        document.head.appendChild(script);
      });
    }
    return this.vtkModulePromise;
  };

  TrajectoryVolumeDisplay.prototype._ensureRaycastView = function () {
    var self = this;
    return this._ensureVtkModule().then(function () {
      if (!self.vtkContainerEl || !global.CryoVolumeRaycastView) return null;
      if (!self.raycastView) {
        self.raycastView = new global.CryoVolumeRaycastView(self.vtkContainerEl);
        self.raycastView.onCameraChanged = function () {
          if (self._suppressVtkCameraCapture || self.backend !== "vtk") return;
          if (typeof self.raycastView.syncChimeraxTurnsFromCamera === "function") {
            self.raycastView.syncChimeraxTurnsFromCamera();
          }
          self._markVtkCameraUserAdjusted();
        };
        var resizeTarget = self._vtkResizeTarget();
        if (resizeTarget) self.raycastView.observeResize(resizeTarget);
      }
      return self.raycastView;
    });
  };

  TrajectoryVolumeDisplay.prototype._renderVtk = function () {
    var self = this;
    if (!this._hasDisplayableVolumes()) {
      this.setStatus("", false);
      if (this.viewportEl) this.viewportEl.hidden = true;
      this._syncChrome();
      return;
    }
    var navTotal = this._volumeNavCount();
    if (navTotal < 1) return;
    var idx = Math.max(0, Math.min(this.vtkFocusIndex, navTotal - 1));
    if (!this._volumeReadyAt(idx)) {
      if (this._countReadyVolumes() > 0) {
        idx = this._nearestReadyVolumeIndex(idx);
        this.vtkFocusIndex = idx;
      }
    }
    while (this.volumes.length < navTotal) this.volumes.push(null);
    var vol = this.volumes[idx];
    if (!vol || !vol.volume_b64) {
      this.setStatus("", false);
      this._syncChrome();
      return;
    }
    this.vtkFocusIndex = idx;
    this.setStatus("Loading 3D viewer…", true);
    if (this.viewportEl) this.viewportEl.hidden = false;
    this._syncChrome();
    this._ensureRaycastView().then(function (view) {
      if (!view) throw new Error("3D viewer unavailable.");
      if (self.raycastVolIndex === idx && self.raycastView === view) {
        self._syncIsoSlider(view);
        if (self._pendingApplyViewTurnsToVtk) self._applyPendingViewTurnsToVtk(view);
        else if (self._pendingApplyViewMatrixToVtk) self._applySharedViewMatrixToVtk(view);
        if (self.vtkContainerEl) {
          self.vtkContainerEl.hidden = false;
        }
        if (self.viewportEl) self.viewportEl.hidden = false;
        self._syncChrome();
        self.setStatus("", false);
        self._scheduleVtkResize();
        return;
      }
      var preserveView = !!(self.raycastView && self.raycastView.volume
        && (self.vtkSessionViewMatrix || self.getSharedViewMatrix()));
      var pendingTurns = !!(self._pendingApplyViewTurnsToVtk && self._pendingApplyViewTurnsToVtk.length);
      var pendingMatrix = !!self._pendingApplyViewMatrixToVtk;
      self._suppressVtkCameraCapture = true;
      view.setVolumeFromB64(vol.volume_b64, vol.D, { skipDefaultCamera: preserveView });
      if (typeof view.setInteractionEnabled === "function") view.setInteractionEnabled(true);
      if (self.vtkContainerEl) self.vtkContainerEl.hidden = false;
      self.raycastVolIndex = idx;
      if (pendingMatrix && self.sharedViewMatrix
        && typeof view.setChimeraxViewMatrixCamera === "function") {
        var matrixOk = view.setChimeraxViewMatrixCamera(self.sharedViewMatrix);
        self._pendingApplyViewMatrixToVtk = false;
        if (!matrixOk && self._vtkViewTurnsBeforeChimerax
          && self._vtkViewTurnsBeforeChimerax.length
          && typeof view.applyChimeraxViewTurns === "function") {
          view.applyChimeraxViewTurns(self._vtkViewTurnsBeforeChimerax);
        }
        self._captureVtkSessionViewMatrix(view);
      } else if (pendingTurns) {
        self._applyPendingViewTurnsToVtk(view);
        self._captureVtkSessionViewMatrix(view);
      } else if (preserveView) {
        self._applyVtkNavigationViewMatrix(view);
      } else {
        self._captureVtkSessionViewMatrix(view);
      }
      self._suppressVtkCameraCapture = false;
      if (self.vtkContainerEl) self.vtkContainerEl.hidden = false;
      self.getChimeraxViewMatrix();
      self._syncIsoSlider(view);
      self.setStatus("", false);
      self._syncVolumeNavChrome();
      self._scheduleVtkResize();
    }).catch(function (err) {
      self.setStatus(err.message || "Failed to load 3D viewer.", false);
    });
  };

  TrajectoryVolumeDisplay.prototype._destroyRaycast = function () {
    if (this.raycastView) {
      this.raycastView.destroy();
      this.raycastView = null;
    }
    this.raycastVolIndex = null;
  };

  TrajectoryVolumeDisplay.prototype._getIsoDataRange = function () {
    if (this.raycastView && typeof this.raycastView.getIsoDataRange === "function") {
      return this.raycastView.getIsoDataRange();
    }
    if (this.isoSliderRange) return this.isoSliderRange;
    return { min: 0, max: 1 };
  };

  TrajectoryVolumeDisplay.prototype._syncIsoSliderFromState = function () {
    if (!this.isoSliderEl || !global.CryoVolume3dUtils) return;
    var range = this._getIsoDataRange();
    var level = this.chimeraxIsoLevel;
    if (this.backend === "vtk" && this.raycastView && typeof this.raycastView.getIsoLevel === "function") {
      level = this.raycastView.getIsoLevel();
    }
    if (level == null || !isFinite(Number(level))) return;
    var slider = global.CryoVolume3dUtils.isoDataValueToSlider(
      Number(level), range.min, range.max
    );
    this.isoSliderEl.value = String(Math.round(slider));
  };

  TrajectoryVolumeDisplay.prototype._syncIsoSlider = function (view) {
    if (!view || !this.isoSliderEl || !global.CryoVolume3dUtils) return;
    var range = view.getIsoDataRange();
    var slider = global.CryoVolume3dUtils.isoDataValueToSlider(
      view.getIsoLevel(), range.min, range.max
    );
    this.isoSliderEl.value = String(Math.round(slider));
    if (this.backend === "chimerax") {
      this.chimeraxIsoLevel = view.getIsoLevel();
    }
  };

  TrajectoryVolumeDisplay.prototype._onIsoInput = function () {
    if (!this.isoSliderEl || !global.CryoVolume3dUtils) return;
    var range = this._getIsoDataRange();
    var level = global.CryoVolume3dUtils.isoSliderToDataValue(
      Number(this.isoSliderEl.value), range.min, range.max
    );
    if (this.backend === "chimerax") {
      this.chimeraxIsoLevel = level;
      if (typeof this.onChimeraxIsoChange === "function") {
        var self = this;
        if (self._chimeraxIsoRerenderTimer) clearTimeout(self._chimeraxIsoRerenderTimer);
        self._chimeraxIsoRerenderTimer = setTimeout(function () {
          self._chimeraxIsoRerenderTimer = null;
          self.onChimeraxIsoChange(level);
        }, 400);
      }
      return;
    }
    if (!this.raycastView) return;
    this.raycastView.setIsoLevel(level);
    this.chimeraxIsoLevel = level;
  };

  TrajectoryVolumeDisplay.prototype._pan = function (dx, dy) {
    if (this.backend === "vtk" && this.raycastView) {
      this.raycastView.panBy(-dx * 5, dy * 5);
      this._markVtkCameraUserAdjusted();
      return;
    }
    if (this.sliceViewer && this.sliceViewer.layers && this.sliceViewer.layers.length) {
      this.sliceViewer.panBy(-dx * TRAJ_VSLICE_PAN_STEP, dy * TRAJ_VSLICE_PAN_STEP);
    }
  };

  TrajectoryVolumeDisplay.prototype._zoomSlice = function (dz) {
    if (this.backend === "vtk" && this.raycastView) {
      if (dz > 0) this.raycastView.dollyBy(1.12);
      else this.raycastView.dollyBy(1 / 1.12);
      this._markVtkCameraUserAdjusted();
      return;
    }
    if (this.sliceViewer && this.sliceViewer.layers && this.sliceViewer.layers.length) {
      if (dz > 0) this.sliceViewer.zoomIn();
      else this.sliceViewer.zoomOut();
    }
  };

  TrajectoryVolumeDisplay.prototype._rememberDisplayRowHome = function () {
    if (!this.displayRowEl || this._displayRowHomeParent) return;
    var parent = this.displayRowEl.parentElement;
    if (!parent || parent === this.expandedHostEl) return;
    this._displayRowHomeParent = parent;
    this._displayRowHomeNext = this.displayRowEl.nextSibling;
  };

  TrajectoryVolumeDisplay.prototype._repatriateDisplayFromExpandedHost = function () {
    if (!this.expandedHostEl || !this.displayRowEl) return;
    if (this.displayRowEl.parentElement !== this.expandedHostEl) return;
    this._rememberDisplayRowHome();
    this._rememberControlsHome();
    if (this._displayRowHomeParent) {
      this._displayRowHomeParent.insertBefore(this.displayRowEl, this._displayRowHomeNext);
    } else if (this.asideShellEl) {
      var anchor = this.controlsDockEl;
      if (anchor && anchor.parentElement === this.asideShellEl) {
        this.asideShellEl.insertBefore(this.displayRowEl, anchor);
      } else {
        this.asideShellEl.appendChild(this.displayRowEl);
      }
    }
    if (this.controlsDockEl && this._controlsHomeParent) {
      this._controlsHomeParent.insertBefore(this.controlsDockEl, this._controlsHomeNext);
    }
    this.expandedBelow = false;
    this.expandedHostEl.hidden = true;
  };
  TrajectoryVolumeDisplay.prototype._rememberViewportHome = function () {
    if (!this.viewportEl || this._viewportHomeParent) return;
    this._viewportHomeParent = this.viewportEl.parentElement;
    this._viewportHomeNext = this.viewportEl.nextSibling;
  };

  TrajectoryVolumeDisplay.prototype._rememberControlsHome = function () {
    if (!this.controlsDockEl || this._controlsHomeParent) return;
    this._controlsHomeParent = this.controlsDockEl.parentElement;
    this._controlsHomeNext = this.controlsDockEl.nextSibling;
  };

  TrajectoryVolumeDisplay.prototype.toggleExpandedBelow = function () {
    this.setExpandedBelow(!this.expandedBelow);
  };

  TrajectoryVolumeDisplay.prototype.setExpandedBelow = function (expanded) {
    if (!this.expandedHostEl) return;
    if (this.backend === "chimerax" && expanded) return;
    this.expandedBelow = !!expanded;
    this._rememberDisplayRowHome();
    this._rememberControlsHome();
    if (this.expandedBelow) {
      if (this.displayRowEl) this.expandedHostEl.appendChild(this.displayRowEl);
      if (this.controlsDockEl) this.expandedHostEl.appendChild(this.controlsDockEl);
      this.expandedHostEl.hidden = false;
    } else {
      if (this.displayRowEl && this._displayRowHomeParent) {
        this._displayRowHomeParent.insertBefore(this.displayRowEl, this._displayRowHomeNext);
      }
      if (this.controlsDockEl && this._controlsHomeParent) {
        this._controlsHomeParent.insertBefore(this.controlsDockEl, this._controlsHomeNext);
      }
      this.expandedHostEl.hidden = true;
    }
    if (this.raycastView) {
      var resizeTarget = this.expandedBelow ? this.expandedHostEl : this._vtkResizeTarget();
      if (resizeTarget && this.raycastView.observeResize) {
        this.raycastView.observeResize(resizeTarget);
      }
      this.raycastView.resize();
    }
    if (this.sliceViewer && this.sliceViewer.resize) this.sliceViewer.resize();
    if (typeof this.onExpandedBelowChange === "function") {
      this.onExpandedBelowChange(this.expandedBelow);
    }
  };

  global.CryoTrajectoryVolumeDisplay = TrajectoryVolumeDisplay;
})(typeof window !== "undefined" ? window : this);
