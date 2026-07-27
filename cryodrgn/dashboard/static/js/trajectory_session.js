/**
 * TrajectorySession — facade owning path, volume state, and volume pipeline.
 *
 * Page controls should dispatch PathMutation commands through this object so
 * geometry and per-slot decode/render status stay aligned.
 *
 * Decode/Render button debts (by construction, not post-hoc clamps):
 *   - ``renderDebtCount()`` === inactive volume-slider ticks (display hook or
 *     VolumeState).
 *   - ``decodeDebtIndices()`` ⊆ those inactive ticks (intersected with slots
 *     that still need decode), so ``decodeDebtCount() ≤ renderDebtCount()``.
 * Prefer ``decodeRenderDebts()`` for the button label — one joint read.
 *
 * ChimeraX view/iso for an in-flight render batch is frozen via
 * ``chimeraxViewSnapshot()`` / Display ``beginChimeraxViewBatch`` so every
 * frame in the batch shares one rotation.
 */
(function (global) {
  "use strict";

  function TrajectorySession(opts) {
    opts = opts || {};
    this.hooks = opts.hooks || {};
    this._store = opts.store || null;
    this._pathHooks = opts.pathHooks || {};
    this._volumes = opts.volumes
      || (global.CryoTrajectoryVolumeState ? new global.CryoTrajectoryVolumeState() : null);
    this._pipeline = opts.pipeline
      || (global.CryoTrajectoryVolumePipeline
        ? new global.CryoTrajectoryVolumePipeline({
            volumes: this._volumes,
            hooks: opts.pipelineHooks || {}
          })
        : null);
    if (this._pipeline && this._volumes) this._pipeline.setVolumeState(this._volumes);
    this._visitOrder = "preserve";
    this._listeners = [];
    if (this._volumes && typeof this._volumes.onChange === "function") {
      var self = this;
      this._volumes.onChange(function (state, reason) {
        self._emit("volumes", { reason: reason, state: state });
        if (typeof self.hooks.onVolumesChanged === "function") {
          self.hooks.onVolumesChanged(state, reason);
        }
      });
    }
  }

  TrajectorySession.prototype.on = function (fn) {
    if (typeof fn === "function") this._listeners.push(fn);
    return this;
  };

  TrajectorySession.prototype._emit = function (type, detail) {
    for (var i = 0; i < this._listeners.length; i++) {
      try { this._listeners[i](type, detail || {}); } catch (err) { /* ignore */ }
    }
  };

  TrajectorySession.prototype.path = function () {
    if (!this._store) return null;
    if (global.CryoTrajectoryPath && typeof CryoTrajectoryPath.ensureCurrentPath === "function") {
      return CryoTrajectoryPath.ensureCurrentPath(this._store, this._pathHooks);
    }
    return this._store.currentPath || null;
  };

  TrajectorySession.prototype.volumes = function () {
    return this._volumes;
  };

  TrajectorySession.prototype.pipeline = function () {
    return this._pipeline;
  };

  /**
   * Inactive volume-slider tick indices (render-debt set).
   * Prefer the live display hook so chrome and button stay aligned.
   */
  TrajectorySession.prototype.inactiveVolumeTickIndices = function () {
    if (typeof this.hooks.inactiveTickIndices === "function") {
      try {
        var hooked = this.hooks.inactiveTickIndices();
        if (Array.isArray(hooked)) {
          return hooked.map(function (i) { return Math.floor(Number(i)); })
            .filter(function (i) { return Number.isFinite(i) && i >= 0; });
        }
      } catch (errIdx) { /* fall through */ }
    }
    var vols = this._volumes;
    if (vols && typeof vols.inactiveTickIndices === "function") {
      return vols.inactiveTickIndices("chimerax");
    }
    var n = this.renderDebtCount();
    var fallback = [];
    for (var i = 0; i < n; i++) fallback.push(i);
    return fallback;
  };

  /**
   * Whether slot ``i`` already has decode material for the button / pipeline.
   * Prefer the page VTK/catalog hook when present.
   */
  TrajectorySession.prototype.isSlotDecoded = function (i) {
    i = Math.floor(Number(i));
    if (!Number.isFinite(i) || i < 0) return false;
    if (typeof this.hooks.isSlotDecoded === "function") {
      try {
        return !!this.hooks.isSlotDecoded(i);
      } catch (errDec) { /* fall through */ }
    }
    var vols = this._volumes;
    return !!(vols && typeof vols.isDecoded === "function" && vols.isDecoded(i));
  };

  /**
   * Render debt for the Decode/Render button.
   * Equals inactive volume-slider ticks when a display is hooked; otherwise
   * Session VolumeState inactive ChimeraX ticks.
   */
  TrajectorySession.prototype.renderDebtCount = function () {
    if (typeof this.hooks.renderDebtCount === "function") {
      try {
        var hooked = this.hooks.renderDebtCount();
        if (Number.isFinite(hooked)) return Math.max(0, Math.floor(hooked));
      } catch (errHook) { /* fall through */ }
    }
    if (typeof this.hooks.inactiveTickCount === "function") {
      try {
        var inactive = this.hooks.inactiveTickCount();
        if (Number.isFinite(inactive)) return Math.max(0, Math.floor(inactive));
      } catch (errInactive) { /* fall through */ }
    }
    var indices = null;
    if (typeof this.hooks.inactiveTickIndices === "function") {
      try {
        indices = this.hooks.inactiveTickIndices();
      } catch (errIdx) { indices = null; }
    }
    if (Array.isArray(indices)) return indices.length;
    var vols = this._volumes;
    if (vols && typeof vols.renderDebtCount === "function") {
      return vols.renderDebtCount();
    }
    if (vols && typeof vols.inactiveTickCount === "function") {
      return vols.inactiveTickCount("chimerax");
    }
    return 0;
  };

  /**
   * Alias of renderDebtCount — documents the slider ↔ button invariant.
   */
  TrajectorySession.prototype.inactiveVolumeTickCount = function () {
    return this.renderDebtCount();
  };

  /**
   * Decode debt indices: inactive ticks that still need decode.
   * Always derived as a subset of ``inactiveVolumeTickIndices()`` so
   * ``decodeDebtCount() ≤ renderDebtCount()`` by set inclusion. Hooks may
   * supply the inactive set and per-slot decode readiness — not a separate
   * decode count that could exceed inactive ticks.
   */
  TrajectorySession.prototype.decodeDebtIndices = function () {
    var inactive = this.inactiveVolumeTickIndices();
    var out = [];
    for (var i = 0; i < inactive.length; i++) {
      var idx = inactive[i];
      if (!this.isSlotDecoded(idx)) out.push(idx);
    }
    return out;
  };

  TrajectorySession.prototype.decodeDebtCount = function () {
    return this.decodeDebtIndices().length;
  };

  /**
   * Joint Decode/Render button debts. ``decode ≤ render`` by set inclusion
   * (``decodeIndices`` is filtered from ``renderIndices``).
   */
  TrajectorySession.prototype.decodeRenderDebts = function () {
    var renderIndices = this.inactiveVolumeTickIndices();
    var decodeIndices = [];
    for (var i = 0; i < renderIndices.length; i++) {
      var idx = renderIndices[i];
      if (!this.isSlotDecoded(idx)) decodeIndices.push(idx);
    }
    return {
      decode: decodeIndices.length,
      render: renderIndices.length,
      decodeIndices: decodeIndices,
      renderIndices: renderIndices
    };
  };

  /**
   * ChimeraX view/iso for renders. Prefers the frozen in-flight batch snapshot
   * from Display / page hooks so every frame shares one rotation.
   */
  TrajectorySession.prototype.chimeraxViewSnapshot = function () {
    if (typeof this.hooks.chimeraxViewSnapshot === "function") {
      try {
        var hooked = this.hooks.chimeraxViewSnapshot();
        if (hooked) return hooked;
      } catch (errSnap) { /* fall through */ }
    }
    return null;
  };

  /**
   * Freeze the current ChimeraX view/iso for the render batch about to start.
   */
  TrajectorySession.prototype.beginChimeraxViewBatch = function (snap) {
    if (typeof this.hooks.beginChimeraxViewBatch === "function") {
      try {
        return this.hooks.beginChimeraxViewBatch(snap) || snap || null;
      } catch (errBegin) { /* fall through */ }
    }
    return snap || this.chimeraxViewSnapshot();
  };

  TrajectorySession.prototype.endChimeraxViewBatch = function () {
    if (typeof this.hooks.endChimeraxViewBatch === "function") {
      try { this.hooks.endChimeraxViewBatch(); } catch (errEnd) { /* ignore */ }
    }
    return this;
  };

  TrajectorySession.prototype.waypointRows = function () {
    var path = this.path();
    if (path && typeof path.getAnchorRows === "function") {
      var rows = path.getAnchorRows();
      if (rows && rows.length >= 2) return rows.slice();
    }
    if (path && typeof path.getSamplePlotRows === "function") {
      var samples = path.getSamplePlotRows();
      if (samples && samples.length >= 2) return samples.slice();
    }
    if (path && typeof path.anchorRowsFromSelection === "function") {
      var selection = path.anchorRowsFromSelection({ includeCustom: true });
      if (selection && selection.rows && selection.rows.length) {
        return selection.rows.slice();
      }
    }
    return [];
  };

  TrajectorySession.prototype.waypointCount = function () {
    var rows = this.waypointRows();
    if (rows && rows.length >= 2) return rows.length;
    var path = this.path();
    if (path && typeof path.getSelectedVolumeIds === "function") {
      var ids = path.getSelectedVolumeIds();
      if (ids && ids.length >= 2) return ids.length;
    }
    return rows ? rows.length : 0;
  };

  TrajectorySession.prototype.sampleCount = function () {
    var path = this.path();
    var hasSamples = !!(path && (
      path.hasInteriorSamples
      || (typeof path.interpolationArmed === "function" && path.interpolationArmed())
    ));
    if (hasSamples && path && typeof path.sampleCount === "function") {
      var n = Math.floor(Number(path.sampleCount()));
      if (Number.isFinite(n) && n >= 2) return n;
    }
    return this.waypointCount();
  };

  TrajectorySession.prototype.visitOrder = function () {
    return this._visitOrder;
  };

  TrajectorySession.prototype.setVisitOrder = function (order) {
    this._visitOrder = String(order || "preserve");
    return this;
  };

  TrajectorySession.prototype.setHooks = function (hooks) {
    this.hooks = Object.assign({}, this.hooks, hooks || {});
    return this;
  };

  TrajectorySession.prototype.setPathHooks = function (hooks) {
    this._pathHooks = Object.assign({}, this._pathHooks, hooks || {});
    return this;
  };

  TrajectorySession.prototype.setPipelineHooks = function (hooks) {
    if (this._pipeline) {
      this._pipeline.hooks = Object.assign({}, this._pipeline.hooks, hooks || {});
    }
    return this;
  };

  /**
   * Apply a PathMutation (or any object with ``apply(session)``).
   */
  TrajectorySession.prototype.prepareForMutation = function () {
    return this;
  };

  TrajectorySession.prototype.dispatch = function (mutation, opts) {
    opts = opts || {};
    if (!mutation || typeof mutation.apply !== "function") {
      return { ok: false, reason: "invalid-mutation" };
    }
    var result;
    try {
      result = mutation.apply(this) || { ok: false };
    } catch (err) {
      result = { ok: false, error: err, reason: "exception" };
    }
    if (result.ok && !opts.silent) {
      this._emit("mutation", { mutation: mutation.name || "mutation", result: result });
      if (typeof this.hooks.afterMutation === "function") {
        try { this.hooks.afterMutation(mutation, result); } catch (err2) { /* ignore */ }
      }
    }
    return result;
  };

  TrajectorySession.prototype.reverse = function () {
    var M = global.CryoTrajectoryPathMutations;
    if (!M) return { ok: false, reason: "no-mutations" };
    return this.dispatch(M.reverse());
  };

  TrajectorySession.prototype.setVisitOrderPath = function (order, orderedRows) {
    var M = global.CryoTrajectoryPathMutations;
    if (!M) return { ok: false, reason: "no-mutations" };
    this.setVisitOrder(order);
    return this.dispatch(M.visitOrder(order, orderedRows));
  };

  TrajectorySession.prototype.appendPlotRows = function (rows, opts) {
    var M = global.CryoTrajectoryPathMutations;
    if (!M) return { ok: false, reason: "no-mutations" };
    return this.dispatch(M.appendPlotRows(rows), opts || {});
  };

  TrajectorySession.prototype.rebuildFromSelection = function (opts) {
    var M = global.CryoTrajectoryPathMutations;
    if (!M) return { ok: false, reason: "no-mutations" };
    return this.dispatch(M.rebuildFromSelection(opts));
  };

  TrajectorySession.prototype.alignVolumesToPath = function () {
    var M = global.CryoTrajectoryPathMutations;
    if (!M) return { ok: false, reason: "no-mutations" };
    return this.dispatch(M.alignVolumes());
  };

  /**
   * Align volume slots to the current path using the SlotModel layout contract.
   */
  TrajectorySession.prototype.syncVolumesToPath = function (opts) {
    opts = opts || {};
    var volumes = this.volumes();
    if (!volumes) return { ok: false, reason: "no-volumes" };
    var path = this.path();
    var layout = {};
    if (typeof this.hooks.volumePathLayout === "function") {
      try { layout = this.hooks.volumePathLayout(path, opts) || {}; } catch (err) { layout = {}; }
    }
    if (opts.pathN != null) layout.pathN = opts.pathN;
    if (opts.nAnchors != null) layout.nAnchors = opts.nAnchors;
    if (opts.nInterp != null) layout.nInterp = opts.nInterp;
    if (opts.compactIds) layout.compactIds = opts.compactIds;
    if (opts.ids) layout.ids = opts.ids;
    if (opts.pathXY) layout.pathXY = opts.pathXY;
    if (opts.pathSamples) layout.pathSamples = opts.pathSamples;
    if (opts.remapByPathXy != null) layout.remapByPathXy = opts.remapByPathXy;
    if (opts.rematchPath != null) layout.rematchPath = opts.rematchPath;
    if (!layout.pathXY && path && typeof path.getPathXY === "function") {
      try { layout.pathXY = path.getPathXY(); } catch (errXy) { /* ignore */ }
    }
    if (!layout.ids && !layout.compactIds && typeof this.hooks.slotIdsForPath === "function") {
      var ids = this.hooks.slotIdsForPath(path);
      if (ids && ids.length) {
        if (layout.pathN && ids.length === layout.pathN) layout.ids = ids;
        else layout.compactIds = ids;
      }
    }
    if (layout.ids && layout.ids.length >= 2) {
      var allowGeomRematch = layout.pathXY && layout.pathXY.length >= 2
        && (layout.remapByPathXy || layout.rematchPath)
        && layout.ids.length !== volumes.slotCount();
      if (!allowGeomRematch) {
        volumes.alignToIds(layout.ids);
        if (layout.pathXY && layout.pathXY.length === layout.ids.length
            && typeof volumes.stampPathXy === "function") {
          volumes.stampPathXy(layout.pathXY);
        }
        if (layout.pathSamples && layout.pathSamples.length === layout.ids.length
            && typeof volumes.rematchToPath === "function" && layout.rematchPath) {
          volumes.rematchToPath(layout.pathSamples, { ids: layout.ids });
        }
        return { ok: true, count: volumes.slotCount() };
      }
    }
    // No full per-slot ``ids`` layout available. ``alignToPath`` below treats
    // ``compactIds`` (or, absent that, ``path.getSelectedVolumeIds()``) as
    // catalog anchors spread across the path — a same-length reorder/resync
    // where compactIds already covers every slot degenerates safely to plain
    // alignToIds there, but a compactIds list *shorter* than the target
    // length is only safe when the path is actually resizing (densify/shrink
    // or bootstrap from empty). If nothing is resizing, a short compactIds
    // list would silently collapse an already-aligned mixed catalog+Other
    // path down to catalog-only anchors, orphaning every Other slot's
    // decode/render identity at its old position — never do that.
    var targetPathN = layout.pathN;
    var compactIds = layout.compactIds && layout.compactIds.length
      ? layout.compactIds
      : (path && typeof path.getSelectedVolumeIds === "function" ? path.getSelectedVolumeIds() : []);
    var compactLen = (compactIds || []).length;
    if (targetPathN == null) {
      var sampleN = (path && typeof path.sampleCount === "function")
        ? Math.floor(Number(path.sampleCount())) || 0
        : 0;
      targetPathN = sampleN >= 2 ? sampleN : Math.max(compactLen, volumes.slotCount());
    }
    var spreadNeeded = compactLen !== targetPathN;
    var resizeOrBootstrap = volumes.slotCount() < 2 || targetPathN !== volumes.slotCount();
    if (spreadNeeded && !resizeOrBootstrap) {
      return { ok: true, count: volumes.slotCount(), preserved: true };
    }
    if (typeof volumes.alignToPath === "function") {
      volumes.alignToPath(path, layout);
    } else if (layout.compactIds) {
      volumes.alignToIds(layout.compactIds);
    }
    return { ok: true, count: volumes.slotCount() };
  };

  TrajectorySession.prototype.clear = function () {
    var M = global.CryoTrajectoryPathMutations;
    if (!M) return { ok: false, reason: "no-mutations" };
    return this.dispatch(M.clear());
  };

  TrajectorySession.prototype.decodeVolumes = function (opts) {
    if (!this._pipeline) return Promise.resolve({ ok: false, reason: "no-pipeline" });
    return this._pipeline.decode(opts || {});
  };

  /**
   * ChimeraX render for outstanding slots. Without explicit ``indices`` /
   * ``forceAll``, the batch is exactly ``inactiveVolumeTickIndices()`` — the
   * same set as Decode/Render button render debt — so progress totals and work
   * stay one number.
   */
  TrajectorySession.prototype.renderVolumes = function (opts) {
    if (!this._pipeline) return Promise.resolve({ ok: false, reason: "no-pipeline" });
    opts = opts || {};
    if (!opts.forceAll) {
      var batch = Array.isArray(opts.indices) ? opts.indices.slice() : null;
      if (!batch || !batch.length) {
        batch = this.inactiveVolumeTickIndices();
      }
      opts = Object.assign({}, opts, { indices: batch });
    }
    return this._pipeline.render(opts);
  };

  /**
   * Mirror VolumeState into the page cache (``lastVolumePayload``, stale map, etc.).
   */
  TrajectorySession.prototype.syncVolumeCache = function () {
    if (!this._volumes || typeof this.hooks.syncVolumeCache !== "function") return this;
    try { this.hooks.syncVolumeCache(this._volumes.snapshot()); } catch (err) { /* ignore */ }
    return this;
  };

  /**
   * Load a volume snapshot into VolumeState (bootstrap / external writes).
   */
  TrajectorySession.prototype.loadVolumeSnapshot = function (snap) {
    if (!this._volumes) return this;
    this._volumes.loadSnapshot(snap || {});
    return this;
  };

  global.CryoTrajectorySession = TrajectorySession;
})(typeof window !== "undefined" ? window : this);
