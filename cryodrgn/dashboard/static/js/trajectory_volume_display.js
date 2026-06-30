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
    this.volumeNavLabelEl = options.volumeNavLabelEl;
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
    this.expectedVolumeCount = null;
    this.chimeraxIsoLevel = null;
    this.isoPercentileSamples = null;
    this.isoSliderRange = null;
    this.onChimeraxIsoChange = options.onChimeraxIsoChange || null;
    this._chimeraxIsoRerenderTimer = null;
    this.onFocusChange = options.onFocusChange || null;
    this.getVolumeNavLabels = options.getVolumeNavLabels || null;
    this._viewportHomeParent = null;
    this._viewportHomeNext = null;
    this._displayRowHomeParent = null;
    this._displayRowHomeNext = null;
    this._controlsHomeParent = null;
    this._controlsHomeNext = null;

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
        if (self.backend === "vtk" && self.raycastView) {
          self._suppressVtkCameraCapture = true;
          self.raycastView.resetCamera();
          self._suppressVtkCameraCapture = false;
          self.vtkCameraUserAdjusted = false;
          self.setSharedViewMatrix("");
        } else if (self.sliceViewer) self.sliceViewer.resetView();
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
    }
    if (this.volumeSliderTicksEl) {
      this.volumeSliderTicksEl.addEventListener("click", function (ev) {
        var tick = ev.target.closest("[data-vol-index]");
        if (!tick) return;
        var idx = Number(tick.getAttribute("data-vol-index"));
        if (!Number.isFinite(idx)) return;
        self.setFocusIndex(idx);
      });
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

  TrajectoryVolumeDisplay.prototype.setChimeraxFocusIndex = function (index) {
    var n = this.chimeraxImages.length;
    if (!n) return;
    var i = Number(index);
    if (!Number.isFinite(i)) i = 0;
    this.chimeraxFocusIndex = Math.max(0, Math.min(n - 1, Math.floor(i)));
    if (this.backend === "chimerax") this._renderChimerax();
    this._notifyFocusChange();
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
    var hasImages = this.chimeraxImages.length > 0;
    if (this.chimeraxFocusIndex >= this.chimeraxImages.length) {
      this.chimeraxFocusIndex = 0;
    }
    if (this.chimeraxFocusIndex < 0) this.chimeraxFocusIndex = 0;
    if (this.canvasEl) this.canvasEl.hidden = true;
    if (this.vtkContainerEl) this.vtkContainerEl.hidden = true;
    if (this.chimeraxPreviewEl) {
      if (hasImages) {
        this.chimeraxPreviewEl.src = this._chimeraxImageSrc(
          this.chimeraxImages[this.chimeraxFocusIndex]
        );
        this.chimeraxPreviewEl.hidden = false;
      } else {
        this.chimeraxPreviewEl.src = "";
        this.chimeraxPreviewEl.hidden = true;
      }
    }
    if (this.viewportEl) this.viewportEl.hidden = !hasImages && !this.chimeraxRendering;
    this._syncChrome();
  };

  TrajectoryVolumeDisplay.prototype.setStatus = function (msg, busy) {
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

  TrajectoryVolumeDisplay.prototype.setChimeraxRendering = function (on) {
    this.chimeraxRendering = !!on;
    if (on) this.setStatus("Rendering…", true);
    else if (this.backend === "chimerax") this.setStatus("", false);
    if (this.backend === "chimerax") {
      this._renderChimerax();
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
    if (backend !== "chimerax") this.chimeraxRendering = false;
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
    } else if (!this.volumes.length) {
      this.expectedVolumeCount = null;
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
    return this.chimeraxImages.length > 0;
  };

  TrajectoryVolumeDisplay.prototype.setFocusIndex = function (index) {
    var i = Number(index);
    if (!Number.isFinite(i)) return;
    i = Math.floor(i);
    if (this.backend === "chimerax") {
      var cxTotal = this._volumeNavCount();
      if (cxTotal < 1) return;
      i = Math.max(0, Math.min(cxTotal - 1, i));
      if (this.chimeraxFocusIndex === i && this.chimeraxImages.length > i) {
        this._syncVolumeNavChrome();
        return;
      }
      this.chimeraxFocusIndex = i;
      if (this.chimeraxImages.length > i) {
        this._ensureChimeraxPreview();
        if (this.chimeraxPreviewEl) {
          this.chimeraxPreviewEl.src = this._chimeraxImageSrc(
            this.chimeraxImages[this.chimeraxFocusIndex]
          );
          this.chimeraxPreviewEl.hidden = false;
        }
        if (this.viewportEl) this.viewportEl.hidden = false;
      } else {
        this._renderChimerax();
      }
      this._notifyFocusChange();
      this._syncVolumeNavChrome();
      return;
    }
    if (this.backend !== "vtk" || !this.volumes.length) return;
    var n = this.expectedVolumeCount != null && this.expectedVolumeCount > 0
      ? this.expectedVolumeCount
      : this.volumes.length;
    if (n < 1) return;
    i = Math.max(0, Math.min(n - 1, i));
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
    this.setFocusIndex(Number(this.volumeSliderEl.value));
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

  TrajectoryVolumeDisplay.prototype._syncVolumeSliderTicks = function (focus, total) {
    if (!this.volumeSliderTicksEl) return;
    var labels = this._volumeNavLabels();
    if (total < 1) {
      this.volumeSliderTicksEl.innerHTML = "";
      return;
    }
    if (this.volumeSliderTicksEl.childElementCount !== total) {
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
        tick.appendChild(lbl);
        this.volumeSliderTicksEl.appendChild(tick);
      }
    } else {
      var tickBtns = this.volumeSliderTicksEl.querySelectorAll(".cryo-vslice-volume-slider-tick");
      for (var uj = 0; uj < tickBtns.length; uj++) {
        var lblEl = tickBtns[uj].querySelector(".cryo-vslice-volume-slider-tick-label");
        if (lblEl) lblEl.textContent = String(labels[uj] != null ? labels[uj] : uj + 1);
      }
    }
    var ticks = this.volumeSliderTicksEl.querySelectorAll(".cryo-vslice-volume-slider-tick");
    for (var tj = 0; tj < ticks.length; tj++) {
      ticks[tj].classList.toggle("cryo-vslice-volume-slider-tick--active", tj === focus);
      ticks[tj].disabled = focus < 0;
    }
  };

  TrajectoryVolumeDisplay.prototype._syncVolumeNavChrome = function () {
    var isVtk = this.backend === "vtk";
    var isChimeraX = this.backend === "chimerax";
    var total = this._volumeNavCount();
    var multi = total > 1;
    var active = multi;
    if (isChimeraX && this.chimeraxRendering) active = false;
    if (this.volumeNavEl) {
      this.volumeNavEl.classList.toggle("cryo-vslice-volume-nav--inactive", !active);
      this.volumeNavEl.setAttribute("aria-disabled", active ? "false" : "true");
      this.volumeNavEl.hidden = !(isVtk || isChimeraX) || total < 1;
    }
    var focus = isChimeraX ? this.chimeraxFocusIndex : this.vtkFocusIndex;
    var total = this._volumeNavCount();
    if (this.btnVolPrev) this.btnVolPrev.disabled = !active;
    if (this.btnVolNext) this.btnVolNext.disabled = !active;
    if (this.volumeSliderEl) {
      this.volumeSliderEl.disabled = !active;
      this.volumeSliderEl.min = "0";
      this.volumeSliderEl.max = String(Math.max(0, total - 1));
      this.volumeSliderEl.step = "1";
      if (active) {
        var clamped = Math.max(0, Math.min(total - 1, focus));
        if (String(this.volumeSliderEl.value) !== String(clamped)) {
          this.volumeSliderEl.value = String(clamped);
        }
      }
    }
    this._syncVolumeSliderTicks(active ? focus : -1, total);
    if (this.volumeNavLabelEl) {
      if (!isVtk && !(isChimeraX && total > 0)) {
        this.volumeNavLabelEl.textContent = "—";
        this.volumeNavLabelEl.classList.add("cryo-vslice-volume-nav-label--inactive");
      } else if (!active) {
        this.volumeNavLabelEl.textContent = "—";
        this.volumeNavLabelEl.classList.add("cryo-vslice-volume-nav-label--inactive");
      } else {
        this.volumeNavLabelEl.classList.remove("cryo-vslice-volume-nav-label--inactive");
        var labels = this._volumeNavLabels();
        var labelText = labels[focus] != null ? String(labels[focus]) : String(focus + 1);
        this.volumeNavLabelEl.textContent =
          labelText + " (" + String(focus + 1) + " / " + String(total) + ")";
      }
    }
  };

  TrajectoryVolumeDisplay.prototype._syncChrome = function () {
    var isChimeraX = this.backend === "chimerax";
    var isVtk = this.backend === "vtk";
    var isSlice = this.backend === "slice";
    var hasVol = this.volumes.some(function (v) { return v && v.volume_b64; });
    var hasCxImages = this.chimeraxImages.length > 0;
    var showChimeraxLoading = isChimeraX && this.chimeraxRendering;
    var showInteractive = !isChimeraX && hasVol;
    var showChimeraxPanel = isChimeraX && (hasCxImages || showChimeraxLoading);
    if (this.asideShellEl) {
      this.asideShellEl.hidden = isChimeraX ? !showChimeraxPanel : false;
    }
    if (this.displayRowEl) {
      this.displayRowEl.hidden = !(showInteractive || showChimeraxPanel);
    }
    if (this.viewportEl && !isChimeraX) {
      this.viewportEl.hidden = !hasVol;
    }
    if (this.controlsDockEl) {
      var showControls = (showInteractive || showChimeraxPanel)
        && (hasVol || hasCxImages || showChimeraxLoading);
      this.controlsDockEl.hidden = !showControls;
    }
    if (this.vtkSliceControlsEl) {
      this.vtkSliceControlsEl.hidden = isChimeraX;
    }
    if (this.chimeraxViewControlsEl) {
      this.chimeraxViewControlsEl.hidden = !isChimeraX;
    }
    if (this.padColumnEl) {
      this.padColumnEl.hidden = isChimeraX || !hasVol || showChimeraxLoading;
    }
    if (this.isoControlsEl) {
      this.isoControlsEl.hidden = !(isVtk || (isChimeraX && (hasCxImages || showChimeraxLoading)));
    }
    if (this.sliceControlsRowEl) {
      this.sliceControlsRowEl.hidden = !isSlice || !hasVol;
    }
    if (this.canvasEl) this.canvasEl.hidden = isVtk || isChimeraX;
    if (this.vtkContainerEl) this.vtkContainerEl.hidden = !isVtk;
    if (this.chimeraxPreviewEl) {
      this.chimeraxPreviewEl.hidden = !isChimeraX || !hasCxImages;
    }
    if (this.rotationLockToolbarEl) {
      this.rotationLockToolbarEl.hidden = isVtk || isChimeraX
        || !(this.sliceViewer && this.sliceViewer.layers && this.sliceViewer.layers.length > 1);
    }
    if (this.btnResetView) {
      this.btnResetView.hidden = isChimeraX;
      this.btnResetView.textContent = isVtk ? "Reset camera" : "Reset view";
      this.btnResetView.setAttribute(
        "aria-label",
        isVtk ? "Reset 3D camera" : "Reset slice view"
      );
    }
    this._syncVolumeNavChrome();
    if (isVtk && this.raycastView) this._scheduleVtkResize();
  };

  TrajectoryVolumeDisplay.prototype._renderCurrent = function () {
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

  TrajectoryVolumeDisplay.prototype._cycleVtkFocus = function (delta) {
    if (this.backend === "chimerax") {
      if (this.chimeraxImages.length < 2) return;
      var cn = this.chimeraxImages.length;
      this.chimeraxFocusIndex = (this.chimeraxFocusIndex + delta + cn) % cn;
      this._renderChimerax();
      this._notifyFocusChange();
      return;
    }
    if (this.backend !== "vtk" || this.volumes.length < 2) return;
    var n = this.volumes.length;
    this.vtkFocusIndex = (this.vtkFocusIndex + delta + n) % n;
    this.raycastVolIndex = null;
    this._syncVolumeNavChrome();
    this._renderVtk();
    this._notifyFocusChange();
  };

  TrajectoryVolumeDisplay.prototype._renderInteractive = function () {
    var self = this;
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
    var idx = Math.max(0, Math.min(this.vtkFocusIndex, this.volumes.length - 1));
    var vol = this.volumes[idx];
    if (!vol || !vol.volume_b64) {
      this.setStatus("Loading volume…", true);
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
