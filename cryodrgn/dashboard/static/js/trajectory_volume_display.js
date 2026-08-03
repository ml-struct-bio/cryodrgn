/**
 * Trajectory creator integrated volume panel: VTK 3D, 2D slice, or ChimeraX PNGs.
 *
 * Invariants with the Decode/Render button (via Session.decodeRenderDebts):
 *   - ``renderDebtCount()`` === ``inactiveTickCount()`` (both use ``_slotTickReadyAt``).
 *   - Decode debt is the inactive-tick set intersected with undecoded slots, so
 *     decode count ≤ inactive tick count by set inclusion. Prefer Session for
 *     the joint ``{decode, render}`` read — Display alone owns the inactive set.
 *   - While a ChimeraX render batch is active, ``chimeraxViewBatchSnapshot()`` is
 *     the frozen view/iso for every frame in that batch (rotation consistency).
 *
 * Two distinct notions of "ready" are used on purpose and must not be mixed:
 *   - ``_slotTickReadyAt`` (nav-ready): the single predicate behind tick color/
 *     click-enablement, Decode/Render debt, nav enablement, and every focus
 *     entry point (slider drag/change, Prev/Next, keyboard, tick clicks) via
 *     ``_nearestNavReadyIndex`` / ``_snapFocusIndexToReady`` / ``_stepReadyFocus``.
 *     A tick is a valid stopping point if and only if it renders active.
 *   - ``_volumeReadyAt`` / ``_chimeraxImageReadyAt`` (paintable): raw "do we
 *     hold bytes locally" checks used only by the paint-fallback family
 *     (``_volumeDisplayIndexForFocus`` and the ``_render*`` methods) to pick
 *     what to actually draw once a focus index has been chosen.
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
    this.underToolsEl = options.underToolsEl || null;
    this.sliceSliderColumnEl = options.sliceSliderColumnEl || null;
    this.isoResetRowEl = options.isoResetRowEl || null;
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
    this._vtkStickyViewTurns = null;
    this.vtkSessionViewMatrix = "";
    this.chimeraxRenderedViewMatrix = "";
    this.chimeraxRenderedViewTurns = [];
    this._vtkViewMatrixBeforeChimerax = "";
    this._vtkViewTurnsBeforeChimerax = [];
    this._vtkChimeraxSyncTurns = [];
    this._vtkChimeraxSyncMatrix = "";
    this._suppressVtkCameraCapture = false;
    this.vtkFocusIndex = 0;
    this.chimeraxFocusIndex = 0;
    this._lastChimeraxDisplayIndex = -1;
    this.chimeraxPreviewEl = null;
    this.raycastVolIndex = null;
    this.expandedBelow = false;
    this.chimeraxRendering = false;
    this.chimeraxRerenderInFlight = false;
    this._pendingChimeraxViewSyncClear = false;
    this.volumeGenerationBusy = false;
    this._explicitJobStatusBusy = false;
    this._chimeraxRerenderReadyMask = null;
    this.incompleteVolumeOverlay = false;
    this.expectedVolumeCount = null;
    this.chimeraxIsoLevel = null;
    /**
     * 2D-slice density floor in map units. Null keeps the full sliceMin–sliceMax
     * window (correct default greyscale). Do not reuse ChimeraX/VTK contour
     * suggestions here — those crush 2D contrast on load.
     */
    this.sliceIsoLevel = null;
    this.isoPercentileSamples = null;
    this.isoSliderRange = null;
    this.onChimeraxIsoChange = options.onChimeraxIsoChange || null;
    this.onResetViewClick = options.onResetViewClick || null;
    this.canResetView = options.canResetView || null;
    this._chimeraxIsoRerenderTimer = null;
    this.onFocusChange = options.onFocusChange || null;
    this.getVolumeNavLabels = options.getVolumeNavLabels || null;
    this.allowUnreadyVolumeNav = options.allowUnreadyVolumeNav || null;
    // Optional: page-level readiness (Session / heap / debt) so slider ticks
    // match Decode/Render accounting when display arrays briefly lag.
    this.slotReadyAt = options.slotReadyAt || null;
    this._volumeSliderTickLabelFontPx = null;
    this._viewportHomeParent = null;
    this._viewportHomeNext = null;
    this._displayRowHomeParent = null;
    this._displayRowHomeNext = null;
    this._controlsHomeParent = null;
    this._controlsHomeNext = null;
    this.stableBackendChrome = !!options.stableBackendChrome;
    /** Frozen ChimeraX view/iso for the in-flight render batch (or null). */
    this._chimeraxViewBatch = null;
    /** Optional page hook: () => { view_matrix, view_turns, viewKey, isoKey, iso_level } */
    this.captureChimeraxViewSnapshot = options.captureChimeraxViewSnapshot || null;
    /**
     * Optional page hook when leaving VTK with a rendered volume: push the live
     * camera into ChimeraX applied-view state so the next PNG batch matches.
     * Receives ``{ view_matrix, view_turns }``.
     */
    this.onSyncVtkViewToChimerax = options.onSyncVtkViewToChimerax || null;

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
        if (typeof self.sliceViewer.setContrastLevel === "function") {
          self.sliceViewer.setContrastLevel(Number(options.sliceContrastEl.value));
        } else if (typeof self.sliceViewer.setContrast === "function") {
          self.sliceViewer.setContrast(Number(options.sliceContrastEl.value));
        }
      });
    }
    if (this.volumeSliderEl) {
      this.volumeSliderEl.addEventListener("input", function (ev) {
        self._applySliderFocus(ev);
      });
      this.volumeSliderEl.addEventListener("change", function (ev) {
        self._applySliderFocus(ev);
      });
    }
    if (this.volumeSliderTicksEl) {
      this.volumeSliderTicksEl.addEventListener("click", function (ev) {
        if (self._volumeNavSuspended()) return;
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
    var overlay = this.renderingOverlayEl;
    if (overlay && overlay.parentElement === this.viewportEl) {
      this.viewportEl.insertBefore(img, overlay);
    } else {
      this.viewportEl.appendChild(img);
    }
    this.chimeraxPreviewEl = img;
  };

  TrajectoryVolumeDisplay.prototype._chimeraxImageSrc = function (b64) {
    if (!b64) return "";
    var s = String(b64).trim();
    if (!s) return "";
    if (s.indexOf("data:image") === 0) return s;
    return "data:image/png;base64," + s;
  };

  /**
   * Raw local-data checks: does this display instance already hold the actual
   * bytes to paint at ``index``? These answer "what can we render right now"
   * and must stay local-only (a hook can never supply pixels). Reserved for
   * the paint-fallback family (``_volumeDisplayIndexForFocus`` and the
   * ``_render*`` methods it feeds) — never use these for navigation/tick
   * readiness decisions; use ``_slotTickReadyAt`` for those instead so nav
   * never disagrees with the tick chrome it drives.
   */
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

  /**
   * Canonical slot-readiness predicate. Every other "is this slot ready"
   * decision in the class — tick chrome, Decode/Render debt, nav enablement,
   * and slider/keyboard/tick stopping points — must call this (directly or
   * via ``inactiveTick*``/``_nearestNavReadyIndex``) so they can never
   * disagree about which slots are active. Prefers the page-level
   * Session/debt-aware ``slotReadyAt`` hook, without recursing through
   * ``_chimeraxImageReadyAt`` (used by those helpers).
   */
  TrajectoryVolumeDisplay.prototype._slotTickReadyAt = function (index) {
    if (this.backend === "chimerax" && this.chimeraxRerenderInFlight) {
      return this._chimeraxTickReadyAt(index);
    }
    if (typeof this.slotReadyAt === "function") return !!this.slotReadyAt(index);
    return this.backend === "chimerax"
      ? this._chimeraxImageReadyAt(index)
      : this._volumeReadyAt(index);
  };

  /**
   * Number of volume-slider ticks (live path length on the panel).
   */
  TrajectoryVolumeDisplay.prototype.tickCount = function () {
    return this._volumeNavCount();
  };

  /**
   * Indices of inactive volume-slider ticks (same chrome as grey ticks).
   */
  TrajectoryVolumeDisplay.prototype.inactiveTickIndices = function () {
    var total = this.tickCount();
    var out = [];
    for (var i = 0; i < total; i++) {
      if (!this._slotTickReadyAt(i)) out.push(i);
    }
    return out;
  };

  TrajectoryVolumeDisplay.prototype.inactiveTickCount = function () {
    return this.inactiveTickIndices().length;
  };

  /**
   * Render debt for the Decode/Render button. Always equals inactiveTickCount
   * (class invariant with the volume-slider chrome).
   */
  TrajectoryVolumeDisplay.prototype.renderDebtCount = function () {
    return this.inactiveTickCount();
  };

  /**
   * Decode debt among inactive ticks. ``isDecodedFn(i)`` should return true when
   * slot ``i`` already has decode material. Result ⊆ inactiveTickIndices, so
   * length ≤ inactiveTickCount by construction.
   */
  TrajectoryVolumeDisplay.prototype.decodeDebtIndices = function (isDecodedFn) {
    var inactive = this.inactiveTickIndices();
    if (typeof isDecodedFn !== "function") return inactive.slice();
    var out = [];
    for (var i = 0; i < inactive.length; i++) {
      var idx = inactive[i];
      if (!isDecodedFn(idx)) out.push(idx);
    }
    return out;
  };

  TrajectoryVolumeDisplay.prototype.decodeDebtCount = function (isDecodedFn) {
    return this.decodeDebtIndices(isDecodedFn).length;
  };

  /**
   * Clone a ChimeraX view/iso snapshot used for one render batch.
   */
  TrajectoryVolumeDisplay.prototype._cloneChimeraxViewSnapshot = function (snap) {
    snap = snap || {};
    return {
      view_matrix: snap.view_matrix ? String(snap.view_matrix) : "",
      view_turns: Array.isArray(snap.view_turns)
        ? snap.view_turns.map(function (t) {
            return {
              axis: String((t && t.axis) || "").toLowerCase(),
              degrees: Number(t && t.degrees)
            };
          }).filter(function (t) {
            return t.axis && Number.isFinite(t.degrees) && Math.abs(t.degrees) > 1e-9;
          })
        : [],
      viewKey: snap.viewKey != null ? String(snap.viewKey) : "",
      isoKey: snap.isoKey != null ? String(snap.isoKey) : "",
      iso_level: snap.iso_level != null && Number.isFinite(Number(snap.iso_level))
        ? Number(snap.iso_level)
        : null
    };
  };

  /**
   * Live ChimeraX view/iso for the next render (or the frozen batch snapshot).
   */
  TrajectoryVolumeDisplay.prototype.chimeraxViewSnapshot = function () {
    if (this._chimeraxViewBatch) {
      return this._cloneChimeraxViewSnapshot(this._chimeraxViewBatch);
    }
    if (typeof this.captureChimeraxViewSnapshot === "function") {
      try {
        var live = this.captureChimeraxViewSnapshot();
        if (live) return this._cloneChimeraxViewSnapshot(live);
      } catch (errSnap) { /* fall through */ }
    }
    return this._cloneChimeraxViewSnapshot({
      view_matrix: "",
      view_turns: [],
      viewKey: "default",
      isoKey: this.chimeraxIsoLevel != null && Number.isFinite(Number(this.chimeraxIsoLevel))
        ? "iso:" + Number(this.chimeraxIsoLevel).toFixed(4)
        : "iso:auto",
      iso_level: this.chimeraxIsoLevel
    });
  };

  TrajectoryVolumeDisplay.prototype.chimeraxViewBatchSnapshot = function () {
    return this._chimeraxViewBatch
      ? this._cloneChimeraxViewSnapshot(this._chimeraxViewBatch)
      : null;
  };

  /**
   * Freeze view/iso for every frame in the current ChimeraX render batch.
   */
  TrajectoryVolumeDisplay.prototype.beginChimeraxViewBatch = function (snap) {
    this._chimeraxViewBatch = this._cloneChimeraxViewSnapshot(
      snap || this.chimeraxViewSnapshot()
    );
    return this.chimeraxViewBatchSnapshot();
  };

  TrajectoryVolumeDisplay.prototype.endChimeraxViewBatch = function () {
    this._chimeraxViewBatch = null;
    return this;
  };

  /**
   * Apply the (batch) view snapshot onto a ChimeraX API payload.
   */
  TrajectoryVolumeDisplay.prototype.applyChimeraxViewToPayload = function (payload, snap) {
    payload = payload || {};
    snap = this._cloneChimeraxViewSnapshot(snap || this.chimeraxViewSnapshot());
    if (snap.view_matrix) {
      payload.view_matrix = snap.view_matrix;
      delete payload.view_turns;
    } else if (snap.view_turns && snap.view_turns.length) {
      payload.view_turns = snap.view_turns.slice();
      delete payload.view_matrix;
    }
    if (snap.iso_level != null && Number.isFinite(snap.iso_level)) {
      payload.iso_level = snap.iso_level;
    }
    return payload;
  };

  TrajectoryVolumeDisplay.prototype._countRenderedChimeraxImages = function () {
    var total = this._volumeNavCount();
    var n = 0;
    for (var i = 0; i < total; i++) {
      if (this._chimeraxImageReadyAt(i)) n++;
    }
    return n;
  };

  /**
   * Nearest index with actual local bytes to paint — the chimerax/vtk-specific
   * halves of the raw paint-fallback used only by ``_volumeDisplayIndexForFocus``.
   */
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

  TrajectoryVolumeDisplay.prototype._allowsSparseInteriorFocus = function () {
    return typeof this.allowUnreadyVolumeNav === "function" && !!this.allowUnreadyVolumeNav();
  };

  /**
   * Which local slot to actually paint for a given focus index. Backend-specific
   * and raw-data-based by necessity — a hook can say a slot is "ready" (decoded
   * server-side) before its bytes have been mirrored into this display instance,
   * but we can only ever paint bytes we already hold. Nav/tick-stop decisions
   * must not use this; they use ``_nearestNavReadyIndex`` instead.
   */
  TrajectoryVolumeDisplay.prototype._volumeDisplayIndexForFocus = function (focus) {
    var t = Math.max(0, Math.min(this._volumeNavCount() - 1, Math.floor(Number(focus))));
    if (this.backend === "chimerax") {
      return this._countRenderedChimeraxImages() > 0
        ? this._nearestRenderedChimeraxIndex(t)
        : t;
    }
    if (this._volumeReadyAt(t)) return t;
    return this._nearestReadyVolumeIndex(t);
  };

  /**
   * Nearest slot that satisfies ``_slotTickReadyAt`` — the same predicate the
   * ticks render active/clickable with (see ``_syncVolumeSliderTicks``). One
   * backend-agnostic search suffices because ``_slotTickReadyAt`` already
   * dispatches on ``this.backend`` internally. Returns -1 when nothing on the
   * path is nav-ready.
   */
  TrajectoryVolumeDisplay.prototype._nearestNavReadyIndex = function (target) {
    var total = this._volumeNavCount();
    if (total < 1) return -1;
    var t = Math.max(0, Math.min(total - 1, Math.floor(Number(target))));
    if (this._slotTickReadyAt(t)) return t;
    for (var d = 1; d < total; d++) {
      if (t - d >= 0 && this._slotTickReadyAt(t - d)) return t - d;
      if (t + d < total && this._slotTickReadyAt(t + d)) return t + d;
    }
    return -1;
  };

  /**
   * Snap a candidate focus index onto the nearest nav-ready slot (or -1 when
   * none are ready). Every navigation entry point — slider drag/change,
   * keyboard stepping, programmatic focus — funnels through this so none of
   * them can stop somewhere the tick chrome disagrees with.
   */
  TrajectoryVolumeDisplay.prototype._snapFocusIndexToReady = function (target) {
    return this._nearestNavReadyIndex(target);
  };

  TrajectoryVolumeDisplay.prototype._stepReadyFocus = function (delta) {
    var total = this._volumeNavCount();
    if (total < 1) return;
    var idx = this.backend === "chimerax" ? this.chimeraxFocusIndex : this.vtkFocusIndex;
    for (var attempt = 0; attempt < total; attempt++) {
      idx = (idx + delta + total) % total;
      if (this._slotTickReadyAt(idx)) {
        this.setFocusIndex(idx);
        return;
      }
    }
  };

  TrajectoryVolumeDisplay.prototype._holdChimeraxImagesForRerender = function () {
    var n = this._volumeNavCount();
    if (n < 1) return;
    while (this.chimeraxImages.length < n) this.chimeraxImages.push(null);
    this._chimeraxRerenderReadyMask = new Array(n);
    for (var hi = 0; hi < n; hi++) {
      this._chimeraxRerenderReadyMask[hi] = this._chimeraxImageReadyAt(hi);
      this.chimeraxImages[hi] = null;
    }
    if (this.chimeraxImages.length > n) {
      this.chimeraxImages.length = n;
    } else {
      while (this.chimeraxImages.length < n) this.chimeraxImages.push(null);
    }
    if (this.chimeraxPreviewEl) {
      this.chimeraxPreviewEl.src = "";
      this.chimeraxPreviewEl.hidden = true;
    }
  };

  TrajectoryVolumeDisplay.prototype._chimeraxTickReadyAt = function (index) {
    if (this.chimeraxRerenderInFlight) {
      if (this._chimeraxRerenderReadyMask
          && index >= 0
          && index < this._chimeraxRerenderReadyMask.length
          && !this._chimeraxRerenderReadyMask[index]) {
        return false;
      }
      return false;
    }
    return this._chimeraxImageReadyAt(index);
  };

  TrajectoryVolumeDisplay.prototype.setChimeraxImageAt = function (index, b64) {
    if (this.chimeraxRerenderInFlight) return;
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
      ? this._volumeDisplayIndexForFocus(this.chimeraxFocusIndex)
      : this.chimeraxFocusIndex;
    this._lastChimeraxDisplayIndex = displayIndex;
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

  TrajectoryVolumeDisplay.prototype.setJobStatus = function (msg, busy) {
    this._explicitJobStatusBusy = !!busy;
    this.setStatus(msg, busy);
    this._syncChrome();
  };

  TrajectoryVolumeDisplay.prototype._volumeDisplayPlaceholderActive = function () {
    if (this.backend === "slice"
        && this._sliceFixedGridN() > 0
        && !this._hasDisplayableVolumes()) {
      return true;
    }
    if (this._hasDisplayableVolumes()
        && !this._explicitJobStatusBusy
        && !this.volumeGenerationBusy
        && !(this.backend === "chimerax" && this.chimeraxRendering)) {
      return false;
    }
    return !!this._explicitJobStatusBusy
      || !!this.volumeGenerationBusy
      || (this.backend === "chimerax" && !!this.chimeraxRendering);
  };

  TrajectoryVolumeDisplay.prototype.setStatus = function (msg, busy) {
    if (busy && !this._hasDisplayableVolumes()
        && !(this.backend === "chimerax" && this.chimeraxRendering)
        && !this.volumeGenerationBusy
        && !this._explicitJobStatusBusy) {
      busy = false;
      msg = "";
    }
    if (!busy && !msg) {
      this._explicitJobStatusBusy = false;
    }
    var label = busy ? (msg || "Rendering…") : "";
    if (this.statusEl) this.statusEl.textContent = label;
    if (this.renderingOverlayEl) {
      this.renderingOverlayEl.hidden = !busy;
      this.renderingOverlayEl.setAttribute("aria-hidden", busy ? "false" : "true");
      this.renderingOverlayEl.classList.toggle("cryo-plot-rendering-overlay--show", !!busy);
      // 2D slice: keep the montage visible under a lower-right corner badge.
      // Other backends: corner badge once something is already on screen —
      // except VTK focus-pending (new slider tick still fetching).
      var vtkFocusPending = this.backend === "vtk"
        && typeof this._volumeReadyAt === "function"
        && !this._volumeReadyAt(this.vtkFocusIndex);
      var nonblocking = !!busy && !vtkFocusPending && (
        this.backend === "slice" || this._hasDisplayableVolumes()
      );
      this.renderingOverlayEl.classList.toggle(
        "cryo-plot-rendering-overlay--nonblocking",
        nonblocking
      );
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

  /**
   * n×n side length for the 2D slice montage, pinned to the final expected
   * volume count (or the number of already-decoded layers, whichever is larger)
   * so cell size does not jump as volumes stream in.
   */
  TrajectoryVolumeDisplay.prototype._sliceFixedGridN = function () {
    var expected = this._expectedVolumeCount();
    var ready = 0;
    for (var i = 0; i < this.volumes.length; i++) {
      if (this.volumes[i] && this.volumes[i].volume_b64) ready++;
    }
    var count = Math.max(expected, ready);
    if (count < 1) return 0;
    return Math.max(1, Math.ceil(Math.sqrt(count)));
  };

  TrajectoryVolumeDisplay.prototype._syncSliceFixedGrid = function () {
    if (!this.sliceViewer || typeof this.sliceViewer.setFixedGridN !== "function") return;
    this.sliceViewer.setFixedGridN(this.backend === "slice" ? this._sliceFixedGridN() : 0);
  };

  /**
   * Raw count of locally-paintable VTK slots — feeds ``_hasDisplayableVolumes``
   * (busy-overlay/placeholder gating) only. Nav enablement uses the nav-ready
   * count in ``_syncVolumeNavChrome`` instead.
   */
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
    var pipelineBusy = this.backend === "chimerax" && !!this.chimeraxRendering && !anyReady;
    var rerenderBusy = this.backend === "chimerax"
      && !!this.chimeraxRendering
      && !!this.chimeraxRerenderInFlight
      && anyReady;
    var showOverlay = pipelineBusy || rerenderBusy
      || this._explicitJobStatusBusy || this.volumeGenerationBusy;
    if (!this.renderingOverlayEl) return;
    if (showOverlay) {
      this.renderingOverlayEl.hidden = false;
      this.renderingOverlayEl.setAttribute("aria-hidden", "false");
      this.renderingOverlayEl.classList.add("cryo-plot-rendering-overlay--show");
      return;
    }
    if (!this._explicitJobStatusBusy) {
      this.renderingOverlayEl.hidden = true;
      this.renderingOverlayEl.setAttribute("aria-hidden", "true");
      this.renderingOverlayEl.classList.remove("cryo-plot-rendering-overlay--show");
    }
  };

  TrajectoryVolumeDisplay.prototype.setVolumeGenerationBusy = function (on) {
    on = !!on;
    if (this.volumeGenerationBusy === on) return;
    this.volumeGenerationBusy = on;
    // Nav chrome only — a full `_syncChrome()` here can toggle overlay
    // attributes and re-enter `syncTrajectoryUiBusy()` via MutationObserver.
    this._syncVolumeNavChrome();
  };

  TrajectoryVolumeDisplay.prototype._volumeNavSuspended = function () {
    if (this.volumeGenerationBusy) return true;
    // Do not lock the volume slider during VTK paint: disabling the range
    // input mid-drag aborts the gesture after a single tick. Overlapping
    // paints are discarded via _vtkPaintGen + vtkFocusIndex checks in
    // _renderVtk instead.
    if (this.backend === "chimerax") {
      if (this.chimeraxRerenderInFlight) return false;
      if (this.chimeraxRendering) return true;
    }
    return false;
  };

  TrajectoryVolumeDisplay.prototype.setChimeraxRendering = function (on, opts) {
    opts = opts || {};
    var startingRerender = !!(on && opts.rerender);
    this.chimeraxRendering = !!on;
    if (!on) {
      this.chimeraxRerenderInFlight = false;
      this._chimeraxRerenderReadyMask = null;
      this.endChimeraxViewBatch();
    } else if (opts.rerender) {
      this.chimeraxRerenderInFlight = true;
      this._holdChimeraxImagesForRerender();
    } else if (opts.rerender === false) {
      this.chimeraxRerenderInFlight = false;
    }
    if (this.backend === "chimerax") {
      this._syncChrome();
      if (!startingRerender || this._hasDisplayableVolumes()) {
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
      this._prepareVtkViewFromChimerax();
    } else if (this.backend === "vtk" && backend === "chimerax") {
      this._prepareChimeraxViewFromVtk();
    }
    if (backend === "chimerax" && this.expandedBelow) {
      this.setExpandedBelow(false);
    }
    this.backend = backend;
    if (backend !== "chimerax") {
      this.chimeraxRendering = false;
      this.chimeraxRerenderInFlight = false;
      this._pendingChimeraxViewSyncClear = false;
    }
    if (this.raycastView && typeof this.raycastView.setInteractionEnabled === "function") {
      this.raycastView.setInteractionEnabled(backend === "vtk");
    }
    if (backend === "chimerax") {
      this._repatriateDisplayFromExpandedHost();
    }
    if (backend === "chimerax" && this._pendingChimeraxViewSyncClear) {
      this._pendingChimeraxViewSyncClear = false;
      this.setChimeraxRendering(true, { rerender: true });
      return;
    }
    this._syncChrome();
    this._renderCurrent();
  };

  /**
   * Leaving ChimeraX for VTK: adopt the viewing angle of the rendered ChimeraX
   * frame when one is on screen.
   *
   * Canonical sync format is ChimeraX ``turn`` commands (resetCamera + orbit).
   * Never apply ChimeraX's absolute ``view matrix`` (includes camera translation
   * that blanks VTK). When a ChimeraX frame is showing, do not restore a stale
   * pre-ChimeraX VTK matrix either — that fought the ChimeraX view on repeated
   * switches and could re-apply a degenerate camera.
   */
  TrajectoryVolumeDisplay.prototype._prepareVtkViewFromChimerax = function () {
    this.raycastVolIndex = null;
    this._pendingApplyViewTurnsToVtk = null;
    this._pendingApplyViewMatrixToVtk = false;
    var hasCxFrame = this._countRenderedChimeraxImages() > 0;
    var cxTurns = [];
    if (hasCxFrame) {
      if (this.chimeraxRenderedViewTurns && this.chimeraxRenderedViewTurns.length) {
        cxTurns = this.chimeraxRenderedViewTurns.slice();
      } else {
        try {
          var snap = this.chimeraxViewSnapshot();
          if (snap && snap.view_turns && snap.view_turns.length) {
            cxTurns = snap.view_turns.slice();
          } else if (snap && snap.view_matrix && isValidChimeraxViewMatrixText(snap.view_matrix)) {
            cxTurns = chimeraxMatrixToViewTurns(snap.view_matrix);
          }
        } catch (errCxSnap) { /* keep empty */ }
      }
      if (!cxTurns.length && this.chimeraxRenderedViewMatrix
          && isValidChimeraxViewMatrixText(this.chimeraxRenderedViewMatrix)) {
        cxTurns = chimeraxMatrixToViewTurns(this.chimeraxRenderedViewMatrix);
      }
      if (cxTurns.length) {
        this._pendingApplyViewTurnsToVtk = cxTurns.slice();
        this.vtkCameraUserAdjusted = true;
      }
      // Identity / default ChimeraX view → VTK resetCamera (no pending matrix).
      return;
    }
    // No ChimeraX frame: restore last VTK camera if we have a sane stash.
    if (this._vtkViewMatrixBeforeChimerax
        && !isDegenerateChimeraxViewMatrix(this._vtkViewMatrixBeforeChimerax)) {
      this.setSharedViewMatrix(this._vtkViewMatrixBeforeChimerax);
      this._pendingApplyViewMatrixToVtk = true;
    } else if (this._vtkViewTurnsBeforeChimerax && this._vtkViewTurnsBeforeChimerax.length) {
      this._pendingApplyViewTurnsToVtk = this._vtkViewTurnsBeforeChimerax.slice();
    }
  };

  /**
   * Map a ChimeraX ``view matrix camera`` rotation (translation ignored) to
   * approximate ``turn y`` / ``turn x`` degrees relative to identity/orient.
   */
  function chimeraxMatrixToViewTurns(vm) {
    if (!vm) return [];
    var raw = String(vm).trim();
    if (raw.toLowerCase().indexOf("camera") === 0) raw = raw.slice(6).trim();
    var parts = raw.match(/[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?/g);
    if (!parts || parts.length < 12) return [];
    var nums = parts.slice(0, 12).map(Number);
    if (!nums.every(function (v) { return isFinite(v); })) return [];
    // Row-major 3x3 rotation from ChimeraX camera matrix layout.
    var r0 = nums[0], r1 = nums[1], r2 = nums[2];
    var r4 = nums[5], r8 = nums[10];
    var eps = 1e-3;
    if (Math.abs(r0 - 1) < eps && Math.abs(r4 - 1) < eps && Math.abs(r8 - 1) < eps
        && Math.abs(nums[1]) < eps && Math.abs(nums[2]) < eps && Math.abs(nums[4]) < eps
        && Math.abs(nums[6]) < eps && Math.abs(nums[8]) < eps && Math.abs(nums[9]) < eps) {
      return [];
    }
    // VTK-style Euler from the ChimeraX camera matrix, then negate to scene-turn
    // degrees (ChimeraX ``turn`` / applyChimeraxViewTurns convention).
    var y = Math.atan2(r2, r0) * (180 / Math.PI);
    var x = -Math.asin(Math.max(-1, Math.min(1, r1))) * (180 / Math.PI);
    function normDeg(d) {
      if (!isFinite(d)) return 0;
      d = d % 360;
      if (d > 180) d -= 360;
      if (d <= -180) d += 360;
      return d;
    }
    y = normDeg(-y);
    x = normDeg(-x);
    var out = [];
    if (Math.abs(y) > 1e-4) out.push({ axis: "y", degrees: y });
    if (Math.abs(x) > 1e-4) out.push({ axis: "x", degrees: x });
    return out;
  }

  function isDegenerateChimeraxViewMatrix(vm) {
    if (!isValidChimeraxViewMatrixText(vm)) return true;
    var raw = String(vm).trim();
    if (raw.toLowerCase().indexOf("camera") === 0) raw = raw.slice(6).trim();
    var parts = raw.match(/[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?/g);
    if (!parts || parts.length < 12) return true;
    var mag = 0;
    for (var i = 0; i < 12; i++) mag += Math.abs(Number(parts[i]));
    return !(mag > 1e-3);
  }

  /**
   * Leaving VTK for ChimeraX: push the live VTK orientation as ChimeraX turns
   * (preferred) so round-trips share one representation.
   */
  TrajectoryVolumeDisplay.prototype._prepareChimeraxViewFromVtk = function () {
    this._vtkChimeraxSyncMatrix = "";
    this._vtkChimeraxSyncTurns = [];
    var hasVtkVolume = !!(this.raycastView && this.raycastView.volume
      && this._countReadyVolumes() > 0);
    if (!hasVtkVolume) return;
    this.getChimeraxViewMatrix();
    var vm = this.getSharedViewMatrix();
    if (vm && !isDegenerateChimeraxViewMatrix(vm)) {
      this._vtkViewMatrixBeforeChimerax = vm;
    } else {
      this._vtkViewMatrixBeforeChimerax = "";
    }
    var vtkTurns = [];
    if (typeof this.raycastView.getChimeraxViewTurns === "function") {
      vtkTurns = this.raycastView.getChimeraxViewTurns() || [];
      if (vtkTurns.length) {
        this._vtkViewTurnsBeforeChimerax = vtkTurns.slice();
      } else {
        this._vtkViewTurnsBeforeChimerax = [];
      }
    }
    this._vtkChimeraxSyncTurns = vtkTurns.slice();
    this._vtkChimeraxSyncMatrix = this._vtkViewMatrixBeforeChimerax;
    var syncPayload = { view_matrix: "", view_turns: [] };
    // Prefer orient-relative ``turn y`` / ``turn x`` for VTK→ChimeraX. ChimeraX
    // runs ``volume center`` then ``view orient`` + turns; a VTK view-matrix with
    // zeroed translation places the camera at the origin and yields a different
    // wrong angle on every orbit. Turns share the same resetCamera/orbit path
    // used for ChimeraX→VTK.
    if (vtkTurns.length) {
      syncPayload.view_turns = vtkTurns.slice();
    }
    this._pendingChimeraxViewSyncClear = !!(
      syncPayload.view_turns.length || syncPayload.view_matrix
    );
    if (typeof this.onSyncVtkViewToChimerax !== "function") return;
    try {
      this.onSyncVtkViewToChimerax(syncPayload);
    } catch (errSync) { /* page hook is best-effort */ }
  };

  TrajectoryVolumeDisplay.prototype.setChimeraxRenderedViewMatrix = function (vm) {
    var text = vm && isValidChimeraxViewMatrixText(vm) ? String(vm).trim() : "";
    this.chimeraxRenderedViewMatrix = text;
  };

  TrajectoryVolumeDisplay.prototype.setChimeraxRenderedViewTurns = function (turns) {
    this.chimeraxRenderedViewTurns = Array.isArray(turns)
      ? turns.map(function (t) {
          return {
            axis: String((t && t.axis) || "").toLowerCase(),
            degrees: Number(t && t.degrees)
          };
        }).filter(function (t) {
          return t.axis && Number.isFinite(t.degrees) && Math.abs(t.degrees) > 1e-9;
        })
      : [];
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
        if (live && !isDegenerateChimeraxViewMatrix(live)) {
          this.setSharedViewMatrix(live);
          return live;
        }
      }
      var renderer = this.raycastView.renderer;
      var cam = renderer && renderer.getActiveCamera && renderer.getActiveCamera();
      if (cam && cam.getViewMatrix) {
        // vtk.js Camera.getViewMatrix() returns a mat4 (no out-arg).
        var mat = cam.getViewMatrix();
        if (mat && mat.length >= 16) {
          var out = [];
          for (var row = 0; row < 3; row++) {
            for (var col = 0; col < 4; col++) {
              out.push(mat[col * 4 + row]);
            }
          }
          out[3] = 0;
          out[7] = 0;
          out[11] = 0;
          if (out.every(function (v) { return isFinite(v); })) {
            var text = out.map(function (v) {
              return Number(v).toPrecision(8);
            }).join(",");
            if (!isDegenerateChimeraxViewMatrix(text)) {
              this.setSharedViewMatrix(text);
            }
          }
        }
      }
    }
    return this.getSharedViewMatrix();
  };

  TrajectoryVolumeDisplay.prototype.resetVtkNavigationCamera = function () {
    this.vtkSessionViewMatrix = "";
    this.vtkCameraUserAdjusted = false;
    this._vtkStickyViewTurns = null;
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
    // Freeze the orbit once per user gesture. Re-deriving axis-angle turns from
    // the live camera on every slider hop drifts, so volumes look out of sync.
    if (this.raycastView
        && typeof this.raycastView.getChimeraxViewTurns === "function") {
      var liveTurns = this.raycastView.getChimeraxViewTurns() || [];
      if (liveTurns.length) {
        this._vtkStickyViewTurns = liveTurns.map(function (t) {
          return { axis: t.axis, degrees: t.degrees };
        });
      }
    }
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
    this._syncVolumeNavChrome();
  };

  TrajectoryVolumeDisplay.prototype.reorderVolumeSlots = function (permutation, opts) {
    opts = opts || {};
    if (!permutation || permutation.length < 2) return;
    var n = permutation.length;
    var prevFocus = this.backend === "chimerax" ? this.chimeraxFocusIndex : this.vtkFocusIndex;
    var permuteArray = global.CryoTrajectoryVolumeStateUtils
      && global.CryoTrajectoryVolumeStateUtils.permuteArray;
    if (typeof permuteArray !== "function") return;

    this.volumes = permuteArray(this.volumes, permutation);
    this.chimeraxImages = permuteArray(this.chimeraxImages, permutation);
    if (Array.isArray(this._chimeraxRerenderReadyMask)) {
      this._chimeraxRerenderReadyMask = permuteArray(this._chimeraxRerenderReadyMask, permutation);
    }
    this.raycastVolIndex = null;

    var nextFocus = opts.focusIndex;
    if (nextFocus == null && prevFocus >= 0 && prevFocus < n) {
      for (var fi = 0; fi < n; fi++) {
        if (permutation[fi] === prevFocus) {
          nextFocus = fi;
          break;
        }
      }
    }
    if (nextFocus != null && Number.isFinite(Number(nextFocus))) {
      nextFocus = Math.max(0, Math.min(n - 1, Math.floor(Number(nextFocus))));
      if (this.backend === "chimerax") this.chimeraxFocusIndex = nextFocus;
      else this.vtkFocusIndex = nextFocus;
    }

    if (this.backend === "chimerax") {
      this._renderChimerax();
    } else if (this.backend === "vtk") {
      this._renderVtk();
    } else if (this.backend === "slice") {
      this._renderInteractive();
    }
    this._notifyFocusChange();
    this._syncVolumeNavChrome();
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
    // Skip full float32 decode when the iso source blob is unchanged (e.g. waypoint
    // append remounts the slider with the same decoded volumes plus empty slots).
    if (this._isoSampleSourceB64 === vol.volume_b64
        && this.isoPercentileSamples
        && this.isoPercentileSamples.length) {
      this._syncIsoSliderFromState();
      this._applySliceIsoLevel();
      return;
    }
    try {
      var U = global.CryoVolume3dUtils;
      var values = U.decodeFloat32Volume(vol.volume_b64, vol.D);
      this.isoPercentileSamples = U.volumePercentileSamples(values);
      this.isoSliderRange = U.isoSliderDataRange(this.isoPercentileSamples);
      if (this.chimeraxIsoLevel == null) {
        this.chimeraxIsoLevel = U.suggestIsoDataValue(this.isoPercentileSamples);
      }
      this._isoSampleSourceB64 = vol.volume_b64;
      this._syncIsoSliderFromState();
      this._applySliceIsoLevel();
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

  TrajectoryVolumeDisplay.prototype.snapFocusToReadyTick = function (preferredIndex) {
    var snapped = this._snapFocusIndexToReady(
      preferredIndex != null ? preferredIndex : this.getFocusIndex()
    );
    if (snapped < 0) return -1;
    this.setFocusIndex(snapped);
    return snapped;
  };

  TrajectoryVolumeDisplay.prototype.setFocusIndex = function (index, opts) {
    opts = opts || {};
    var i = Number(index);
    if (!Number.isFinite(i)) return;
    var totalOpen = this._volumeNavCount();
    if (totalOpen < 1) return;
    if (opts.allowInterior && this._allowsSparseInteriorFocus()) {
      i = Math.max(0, Math.min(totalOpen - 1, Math.floor(i)));
    } else {
      i = this._snapFocusIndexToReady(i);
      if (i < 0) return;
    }
    if (this.backend === "chimerax") {
      var cxTotal = this._volumeNavCount();
      if (cxTotal < 1) return;
      var cxDisplay = this._volumeDisplayIndexForFocus(i);
      if (this.chimeraxFocusIndex === i
          && this._lastChimeraxDisplayIndex === cxDisplay) {
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
    var vtkDisplay = this._volumeDisplayIndexForFocus(i);
    if (this.vtkFocusIndex === i && this.raycastVolIndex === vtkDisplay) {
      this._syncVolumeNavChrome();
      return;
    }
    this.vtkFocusIndex = i;
    this.raycastVolIndex = null;
    this._syncVolumeNavChrome();
    this._renderVtk();
    this._notifyFocusChange();
  };

  TrajectoryVolumeDisplay.prototype._applySliderFocus = function (ev) {
    if (!this.volumeSliderEl) return;
    if (this._volumeNavSuspended()) return;
    var raw = Number(this.volumeSliderEl.value);
    var dragging = !!(ev && ev.type === "input");
    if (dragging && this._allowsSparseInteriorFocus()) {
      this.setFocusIndex(raw, { allowInterior: true });
      return;
    }
    var snapped = this._snapFocusIndexToReady(raw);
    if (snapped < 0) return;
    if (String(this.volumeSliderEl.value) !== String(snapped)) {
      this.volumeSliderEl.value = String(snapped);
    }
    this.setFocusIndex(snapped);
  };

  TrajectoryVolumeDisplay.prototype._volumeNavCount = function () {
    var isChimeraX = this.backend === "chimerax";
    var isInteractive = this.backend === "vtk" || this.backend === "slice";
    if (isChimeraX) {
      if (this.expectedVolumeCount != null && this.expectedVolumeCount > 0) {
        return this.expectedVolumeCount;
      }
      return this.chimeraxImages.length;
    }
    if (isInteractive) {
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
    var navSuspended = this._volumeNavSuspended();
    var pipelineBusy = isChimeraX && this.chimeraxRendering;
    var rerenderBusy = isChimeraX && this.chimeraxRerenderInFlight;
    var allowUnreadyNav = this._allowsSparseInteriorFocus();
    for (var tj = 0; tj < ticks.length; tj++) {
      var tickReady = this._slotTickReadyAt(tj);
      ticks[tj].classList.toggle("cryo-vslice-volume-slider-tick--active", !navSuspended && tj === focus);
      if (navSuspended) {
        ticks[tj].classList.remove("cryo-vslice-volume-slider-tick--pending");
        ticks[tj].classList.add("cryo-vslice-volume-slider-tick--inactive");
        ticks[tj].disabled = true;
        continue;
      }
      var isPending = false;
      if (rerenderBusy) {
        isPending = !!(this._chimeraxRerenderReadyMask
          && tj >= 0
          && tj < this._chimeraxRerenderReadyMask.length
          && this._chimeraxRerenderReadyMask[tj]);
      } else if (pipelineBusy) {
        isPending = !tickReady;
      }
      ticks[tj].classList.toggle("cryo-vslice-volume-slider-tick--pending", isPending);
      ticks[tj].classList.toggle(
        "cryo-vslice-volume-slider-tick--inactive",
        !tickReady && !isPending && tj !== focus
      );
      ticks[tj].disabled = focus < 0 || (!allowUnreadyNav && !tickReady);
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
    // Nav-ready count via the same predicate as the tick chrome below
    // (_slotTickReadyAt, through inactiveTickCount) — nav enablement,
    // Prev/Next stepping, and tick color/click-enablement must always agree.
    var readyCount = total - this.inactiveTickCount();
    var multi = total > 1;
    var navSuspended = this._volumeNavSuspended();
    var active = multi && readyCount > 0 && !navSuspended;
    if (this.volumeNavEl) {
      this.volumeNavEl.classList.toggle("cryo-vslice-volume-nav--inactive", !active);
      this.volumeNavEl.setAttribute("aria-disabled", active ? "false" : "true");
      this.volumeNavEl.hidden = !(isVtk || isChimeraX) || total < 1;
    }
    var focus = isChimeraX ? this.chimeraxFocusIndex : this.vtkFocusIndex;
    var navStepEnabled = active && (this._allowsSparseInteriorFocus() ? total > 1 : readyCount > 1);
    if (this.btnVolPrev) this.btnVolPrev.disabled = !navStepEnabled;
    if (this.btnVolNext) this.btnVolNext.disabled = !navStepEnabled;
    if (this.volumeSliderEl) {
      this.volumeSliderEl.disabled = !active;
      this.volumeSliderEl.min = "0";
      this.volumeSliderEl.max = String(Math.max(0, total - 1));
      this.volumeSliderEl.step = "1";
      if (active && focus >= 0) {
        var sliderIdx = Math.max(0, Math.min(total - 1, focus));
        if (String(this.volumeSliderEl.value) !== String(sliderIdx)) {
          this.volumeSliderEl.value = String(sliderIdx);
        }
        var navLabels = this._volumeNavLabels();
        var navText = navLabels[sliderIdx] != null
          ? String(navLabels[sliderIdx])
          : String(sliderIdx + 1);
        this.volumeSliderEl.setAttribute("aria-valuetext", navText);
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

  /**
   * Host pan/zoom + iso/reset (and slice contrast) beside each other under the
   * viewport for VTK and 2D slice. ChimeraX keeps iso+reset in the controls dock.
   */
  TrajectoryVolumeDisplay.prototype._syncSliceControlHosts = function () {
    var under = this.underToolsEl;
    var sliderCol = this.sliceSliderColumnEl;
    var isoResetRow = this.isoResetRowEl;
    var controls = this.vtkSliceControlsEl;
    if (!under || !sliderCol || !isoResetRow) return;
    var isSlice = this.backend === "slice";
    var isVtk = this.backend === "vtk";
    var besideUnder = isSlice || isVtk;

    if (besideUnder) {
      if (this.padColumnEl && this.padColumnEl.parentElement !== under) {
        under.insertBefore(this.padColumnEl, sliderCol);
      }
      if (isSlice && this.rotationLockToolbarEl) {
        sliderCol.appendChild(this.rotationLockToolbarEl);
      } else if (this.rotationLockToolbarEl && controls) {
        controls.insertBefore(this.rotationLockToolbarEl, controls.firstChild);
      }
      if (this.isoControlsEl) sliderCol.appendChild(this.isoControlsEl);
      if (isSlice && this.sliceControlsRowEl) {
        sliderCol.appendChild(this.sliceControlsRowEl);
      } else if (this.sliceControlsRowEl && controls) {
        controls.appendChild(this.sliceControlsRowEl);
      }
      if (this.btnResetView) sliderCol.appendChild(this.btnResetView);
      sliderCol.hidden = false;
      under.classList.add("cryo-vslice-under-tools--beside");
    } else {
      if (this.rotationLockToolbarEl && controls) {
        controls.insertBefore(this.rotationLockToolbarEl, controls.firstChild);
      }
      if (this.isoControlsEl) {
        if (this.btnResetView && isoResetRow.contains(this.btnResetView)) {
          isoResetRow.insertBefore(this.isoControlsEl, this.btnResetView);
        } else {
          isoResetRow.insertBefore(this.isoControlsEl, isoResetRow.firstChild);
        }
      }
      if (this.btnResetView) isoResetRow.appendChild(this.btnResetView);
      if (this.sliceControlsRowEl && controls) {
        controls.appendChild(this.sliceControlsRowEl);
      }
      sliderCol.hidden = true;
      under.classList.remove("cryo-vslice-under-tools--beside");
    }
  };

  TrajectoryVolumeDisplay.prototype._applySliceIsoLevel = function () {
    if (this.backend !== "slice" || !this.sliceViewer) return;
    if (typeof this.sliceViewer.setIsoLevel !== "function") return;
    // null → full dynamic range (legacy/default 2D appearance).
    this.sliceViewer.setIsoLevel(this.sliceIsoLevel);
  };

  TrajectoryVolumeDisplay.prototype._syncChrome = function () {
    var isChimeraX = this.backend === "chimerax";
    var isVtk = this.backend === "vtk";
    var isSlice = this.backend === "slice";
    if (this.asideShellEl) {
      this.asideShellEl.setAttribute("data-vol-backend", this.backend);
    }
    this._syncSliceControlHosts();
    var hasVol = this.volumes.some(function (v) { return v && v.volume_b64; });
    var hasCxImages = this._countRenderedChimeraxImages() > 0;
    var placeholder = this._volumeDisplayPlaceholderActive();
    var showChimeraxPanel = isChimeraX;
    var showInteractive = !isChimeraX && hasVol;
    if (this.asideShellEl) {
      this.asideShellEl.hidden = isChimeraX ? !showChimeraxPanel : false;
    }
    if (this.displayRowEl) {
      this.displayRowEl.hidden = !(showInteractive || showChimeraxPanel || placeholder);
    }
    if (this.viewportEl && isChimeraX) {
      this.viewportEl.hidden = !(hasCxImages || placeholder);
    } else if (this.viewportEl) {
      this.viewportEl.hidden = !(hasVol || placeholder);
    }
    // VTK / slice host iso (+ contrast) beside the pan pad under the viewer;
    // the dock is only needed for ChimeraX control blocks.
    var showControlsDock = isChimeraX
      && showChimeraxPanel
      && (hasVol || hasCxImages);
    if (this.controlsDockEl) {
      this.controlsDockEl.hidden = !showControlsDock;
    }
    if (this.underToolsEl) {
      // Pan / iso / contrast stay hidden while volumes are still loading.
      var showUnder = !isChimeraX && hasVol;
      this.underToolsEl.hidden = !showUnder;
      this.underToolsEl.classList.toggle(
        "cryo-vslice-under-tools--beside",
        (isSlice || isVtk) && showUnder
      );
    }
    if (this.vtkSliceControlsEl) {
      if (isSlice || isVtk) {
        this.vtkSliceControlsEl.hidden = true;
        this.vtkSliceControlsEl.classList.remove("cryo-traj-vol-backend-panel--inactive");
      } else if (this.stableBackendChrome && showControlsDock) {
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
      this.isoControlsEl.hidden = !(
        (isVtk && hasVol) || (isChimeraX && hasCxImages) || (isSlice && hasVol)
      );
    }
    if (this.sliceControlsRowEl) {
      this.sliceControlsRowEl.hidden = !isSlice || !hasVol;
    }
    if (this.sliceSliderColumnEl) {
      this.sliceSliderColumnEl.hidden = !((isSlice || isVtk) && hasVol);
    }
    if (this.canvasEl) this.canvasEl.hidden = isVtk || isChimeraX;
    if (this.vtkContainerEl) this.vtkContainerEl.hidden = !isVtk;
    if (isSlice) {
      this._syncSliceFixedGrid();
      this._applySliceIsoLevel();
      this._syncIsoSliderFromState();
    }
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
    if (this.renderingOverlayEl
        && this.renderingOverlayEl.classList.contains("cryo-plot-rendering-overlay--show")) {
      // Slice always uses the lower-right badge; other backends only once content exists.
      // VTK focus-pending keeps a full blocking overlay while the new tick loads.
      var vtkFocusPendingChrome = isVtk
        && typeof this._volumeReadyAt === "function"
        && !this._volumeReadyAt(this.vtkFocusIndex);
      this.renderingOverlayEl.classList.toggle(
        "cryo-plot-rendering-overlay--nonblocking",
        !vtkFocusPendingChrome && (isSlice || this._hasDisplayableVolumes())
      );
    }
    if (isVtk && this.raycastView) this._scheduleVtkResize();
  };

  TrajectoryVolumeDisplay.prototype._renderCurrent = function () {
    if (!this._hasDisplayableVolumes()) {
      if (this.chimeraxPreviewEl) {
        this.chimeraxPreviewEl.hidden = true;
        this.chimeraxPreviewEl.src = "";
      }
      // Keep an in-flight job/generation message; only clear when truly idle.
      if (!this._explicitJobStatusBusy && !this.volumeGenerationBusy
          && !(this.backend === "chimerax" && this.chimeraxRendering)) {
        this.setStatus("", false);
      }
      // Slice: still paint the final-sized empty montage while volumes load.
      if (this.backend === "slice") {
        this._renderInteractive();
        return;
      }
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

  TrajectoryVolumeDisplay.prototype._cycleVtkFocus = function (delta) {
    if (this._volumeNavSuspended()) return;
    if (this.backend === "chimerax" || this.backend === "vtk") {
      this._stepReadyFocus(delta);
    }
  };

  TrajectoryVolumeDisplay.prototype._renderInteractive = function () {
    if (this.backend === "vtk") {
      if (this.sliceViewer) {
        this._syncSliceFixedGrid();
        this.sliceViewer.setVolumes([], false);
      }
      if (!this._hasDisplayableVolumes()) {
        if (this.viewportEl) this.viewportEl.hidden = true;
        if (this.controlsDockEl) this.controlsDockEl.hidden = true;
        this._destroyRaycast();
        if (!this._explicitJobStatusBusy && !this.volumeGenerationBusy) {
          this.setStatus("", false);
        }
        this._syncChrome();
        return;
      }
      this._renderVtk();
      return;
    }
    // 2D slice montage
    this._destroyRaycast();
    this._syncSliceFixedGrid();
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
    var showEmptyGrid = !layers.length && this._sliceFixedGridN() > 0;
    if (!layers.length && !showEmptyGrid) {
      if (this.viewportEl) this.viewportEl.hidden = true;
      if (this.controlsDockEl) this.controlsDockEl.hidden = true;
      if (this.sliceViewer) this.sliceViewer.setVolumes([], false);
      if (!this._explicitJobStatusBusy && !this.volumeGenerationBusy) {
        this.setStatus("", false);
      }
      this._syncChrome();
      return;
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
        if (!self.raycastView) return;
        self.raycastView.resize();
        if (self.raycastView.renderWindow
            && typeof self.raycastView.renderWindow.render === "function") {
          self.raycastView.renderWindow.render();
        }
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
    var navTotal = this._volumeNavCount();
    if (navTotal < 1) return Promise.resolve(false);
    var focus = Math.max(0, Math.min(navTotal - 1, Math.floor(Number(this.vtkFocusIndex)) || 0));
    this.vtkFocusIndex = focus;
    while (this.volumes.length < navTotal) this.volumes.push(null);
    // One volume at a time: never paint a neighbour while the focused slot
    // is still fetching — keep the loading overlay instead.
    if (!this._volumeReadyAt(focus)) {
      if (!this._explicitJobStatusBusy && !this.volumeGenerationBusy) {
        this.setJobStatus("Loading volume\u2026", true);
      }
      if (this.viewportEl) this.viewportEl.hidden = false;
      this._syncChrome();
      return Promise.resolve(false);
    }
    var idx = focus;
    var vol = this.volumes[idx];
    if (!vol || !vol.volume_b64) {
      if (!this._explicitJobStatusBusy && !this.volumeGenerationBusy) {
        this.setJobStatus("", false);
      }
      this._syncChrome();
      return Promise.resolve(false);
    }
    // Keep an explicit job busy flag through setVolumeFromB64's clear→add gap
    // so the overlay covers the blank frame between volumes. Slider stays
    // enabled; stale paints abort when vtkFocusIndex / _vtkPaintGen diverge.
    this._vtkPaintGen = (this._vtkPaintGen || 0) + 1;
    var paintGen = this._vtkPaintGen;
    this.setJobStatus("Loading volume\u2026", true);
    if (this.viewportEl) this.viewportEl.hidden = false;
    this._syncChrome();
    return this._ensureRaycastView().then(function (view) {
      if (!view) throw new Error("3D viewer unavailable.");
      if (self.backend !== "vtk" || self.vtkFocusIndex !== idx || paintGen !== self._vtkPaintGen) {
        return false;
      }
      if (self.raycastVolIndex === idx && self.raycastView === view) {
        self._syncIsoSlider(view);
        if (self._pendingApplyViewTurnsToVtk) self._applyPendingViewTurnsToVtk(view);
        else if (self._pendingApplyViewMatrixToVtk) self._applySharedViewMatrixToVtk(view);
        if (self.vtkContainerEl) self.vtkContainerEl.hidden = false;
        if (self.viewportEl) self.viewportEl.hidden = false;
        self.setJobStatus("", false);
        self._syncChrome();
        self._scheduleVtkResize();
        return true;
      }
      var prevRaycastIdx = self.raycastVolIndex;
      var switchingVolume = prevRaycastIdx !== idx;
      var hadVolume = !!(self.raycastView && self.raycastView.volume);
      var savedTurns = null;
      // Only reuse a frozen user orbit. Never re-read live turns on switch —
      // that round-trips axis-angle and desyncs the initial default path and
      // post-rotate slider hops alike. No sticky → pure default camera.
      if (switchingVolume && hadVolume
          && self.vtkCameraUserAdjusted
          && self._vtkStickyViewTurns
          && self._vtkStickyViewTurns.length) {
        savedTurns = self._vtkStickyViewTurns;
      }
      var preserveView = !!(hadVolume
        && !switchingVolume
        && (self.vtkSessionViewMatrix || self.getSharedViewMatrix()));
      var pendingTurns = !!(self._pendingApplyViewTurnsToVtk && self._pendingApplyViewTurnsToVtk.length);
      var pendingMatrix = !!self._pendingApplyViewMatrixToVtk;
      self._suppressVtkCameraCapture = true;
      // Always reset the default camera when the focused slot changes —
      // preserving the prior camera across setVolumeFromB64 left a blank
      // viewport after the first volume (clipping / framing race).
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
        // Freeze ChimeraX→VTK pending orbit so later slider hops stay aligned.
        if (self.vtkCameraUserAdjusted
            && typeof view.getChimeraxViewTurns === "function") {
          var pendingFrozen = view.getChimeraxViewTurns() || [];
          if (pendingFrozen.length) {
            self._vtkStickyViewTurns = pendingFrozen.map(function (t) {
              return { axis: t.axis, degrees: t.degrees };
            });
          }
        }
      } else if (switchingVolume && savedTurns && savedTurns.length
          && typeof view.applyChimeraxViewTurns === "function") {
        // Re-apply frozen orbit after this slot's default frame. Full reset
        // inside apply keeps the turn basis identical for every volume.
        view.applyChimeraxViewTurns(savedTurns);
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
      self.setJobStatus("", false);
      self._syncChrome();
      self._syncVolumeNavChrome();
      if (typeof view.resize === "function") view.resize();
      if (view.renderWindow && typeof view.renderWindow.render === "function") {
        view.renderWindow.render();
      }
      self._scheduleVtkResize();
      return true;
    }).catch(function (err) {
      self.setJobStatus(err.message || "Failed to load 3D viewer.", false);
      return false;
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
    if (this.backend === "slice") {
      // Default (null) shows the bottom of the iso window so 2D greyscale is
      // not clipped by a ChimeraX/VTK contour suggestion.
      level = this.sliceIsoLevel != null ? this.sliceIsoLevel : range.min;
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
    if (this.backend === "slice") {
      this.sliceIsoLevel = level;
      this._applySliceIsoLevel();
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
    // Pop-out expansion for 2D slice / VTK 3D is disabled; ChimeraX uses the
    // gallery column in the modal instead of relocating this display row.
    if (expanded || !this.expandedBelow) return;
    this.expandedBelow = false;
    this._rememberDisplayRowHome();
    this._rememberControlsHome();
    if (this.displayRowEl && this._displayRowHomeParent) {
      this._displayRowHomeParent.insertBefore(this.displayRowEl, this._displayRowHomeNext);
    }
    if (this.controlsDockEl && this._controlsHomeParent) {
      this._controlsHomeParent.insertBefore(this.controlsDockEl, this._controlsHomeNext);
    }
    this.expandedHostEl.hidden = true;
    if (this.raycastView) {
      var resizeTarget = this._vtkResizeTarget();
      if (resizeTarget && this.raycastView.observeResize) {
        this.raycastView.observeResize(resizeTarget);
      }
      this.raycastView.resize();
    }
    if (this.sliceViewer && this.sliceViewer.resize) this.sliceViewer.resize();
    if (typeof this.onExpandedBelowChange === "function") {
      this.onExpandedBelowChange(false);
    }
  };

  global.CryoTrajectoryVolumeDisplay = TrajectoryVolumeDisplay;
})(typeof window !== "undefined" ? window : this);
