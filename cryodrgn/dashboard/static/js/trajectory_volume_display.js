/**
 * Trajectory creator integrated volume panel: VTK 3D (default for manual analyze volumes), 2D slice, or ChimeraX PNGs.
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
    this.onFocusChange = options.onFocusChange || null;
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

  TrajectoryVolumeDisplay.prototype._captureVtkViewMatrix = function () {
    if (!this.vtkCameraUserAdjusted) return;
    this.getChimeraxViewMatrix();
  };

  TrajectoryVolumeDisplay.prototype._markVtkCameraUserAdjusted = function () {
    this.vtkCameraUserAdjusted = true;
    this.getChimeraxViewMatrix();
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
    if (this.vtkFocusIndex >= this.volumes.length) this.vtkFocusIndex = 0;
    if (this.vtkFocusIndex < 0) this.vtkFocusIndex = 0;
    if (this.chimeraxFocusIndex >= this.chimeraxImages.length) this.chimeraxFocusIndex = 0;
    if (this.chimeraxFocusIndex < 0) this.chimeraxFocusIndex = 0;
    if (!opts.deferRender) this._renderCurrent();
  };

  TrajectoryVolumeDisplay.prototype.hasInteractiveVolumes = function () {
    return this.volumes.length > 0;
  };

  TrajectoryVolumeDisplay.prototype.hasChimeraxImages = function () {
    return this.chimeraxImages.length > 0;
  };

  TrajectoryVolumeDisplay.prototype._syncVolumeNavChrome = function () {
    var isVtk = this.backend === "vtk";
    var isChimeraX = this.backend === "chimerax";
    var multi = (isVtk && this.volumes.length > 1)
      || (isChimeraX && this.chimeraxImages.length > 1);
    var active = multi;
    if (this.volumeNavEl) {
      this.volumeNavEl.classList.toggle("cryo-vslice-volume-nav--inactive", !active);
      this.volumeNavEl.setAttribute("aria-disabled", active ? "false" : "true");
      this.volumeNavEl.hidden = !(isVtk || (isChimeraX && this.chimeraxImages.length > 0));
    }
    if (this.btnVolPrev) this.btnVolPrev.disabled = !active;
    if (this.btnVolNext) this.btnVolNext.disabled = !active;
    if (this.volumeNavLabelEl) {
      if (!isVtk && !(isChimeraX && this.chimeraxImages.length > 0)) {
        this.volumeNavLabelEl.textContent = "—";
        this.volumeNavLabelEl.classList.add("cryo-vslice-volume-nav-label--inactive");
      } else if (!active) {
        this.volumeNavLabelEl.textContent = "—";
        this.volumeNavLabelEl.classList.add("cryo-vslice-volume-nav-label--inactive");
      } else {
        this.volumeNavLabelEl.classList.remove("cryo-vslice-volume-nav-label--inactive");
        var focus = isChimeraX ? this.chimeraxFocusIndex : this.vtkFocusIndex;
        var total = isChimeraX ? this.chimeraxImages.length : this.volumes.length;
        this.volumeNavLabelEl.textContent =
          String(focus + 1) + " / " + String(total);
      }
    }
  };

  TrajectoryVolumeDisplay.prototype._syncChrome = function () {
    var isChimeraX = this.backend === "chimerax";
    var isVtk = this.backend === "vtk";
    var isSlice = this.backend === "slice";
    var hasVol = this.volumes.length > 0;
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
    if (this.isoControlsEl) this.isoControlsEl.hidden = !isVtk;
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
    var layers = this.volumes.map(function (v, i) {
      return {
        id: "traj:" + i,
        b64: v.volume_b64,
        d: v.D,
        label: String(i + 1),
        manual: false
      };
    });
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
    if (!vol || !vol.volume_b64) return;
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
      var pendingTurns = !!(self._pendingApplyViewTurnsToVtk && self._pendingApplyViewTurnsToVtk.length);
      var pendingMatrix = !!self._pendingApplyViewMatrixToVtk;
      self._suppressVtkCameraCapture = true;
      view.setVolumeFromB64(vol.volume_b64, vol.D, { skipDefaultCamera: false });
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
      } else if (pendingTurns) {
        self._applyPendingViewTurnsToVtk(view);
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

  TrajectoryVolumeDisplay.prototype._syncIsoSlider = function (view) {
    if (!view || !this.isoSliderEl || !global.CryoVolume3dUtils) return;
    var range = view.getIsoDataRange();
    var slider = global.CryoVolume3dUtils.isoDataValueToSlider(
      view.getIsoLevel(), range.min, range.max
    );
    this.isoSliderEl.value = String(Math.round(slider));
  };

  TrajectoryVolumeDisplay.prototype._onIsoInput = function () {
    if (!this.raycastView || !this.isoSliderEl || !global.CryoVolume3dUtils) return;
    var range = this.raycastView.getIsoDataRange();
    var level = global.CryoVolume3dUtils.isoSliderToDataValue(
      Number(this.isoSliderEl.value), range.min, range.max
    );
    this.raycastView.setIsoLevel(level);
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
