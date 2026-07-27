/**
 * Direct-trace / nearest UI state for the trajectory creator.
 *
 * Owns tick invalidation maps (stale / detached / drag) and the pure helpers
 * that keep Trajectory controls, path length, and volume chrome aligned after
 * path mutations. Page templates supply hooks for Session / slider sync.
 *
 * Hierarchy (with existing classes):
 *   DirectTracePath     — latent polyline geometry + endpoint catalog ids
 *   TrajectoryVolumeState / Session — decode/render + traj_xy / path_t identity
 *   DirectTraceUiState  — control ↔ path ↔ volume-chrome consistency
 */
(function (global) {
  "use strict";

  function finiteInt(v, fallback) {
    var n = Math.floor(Number(v));
    return Number.isFinite(n) ? n : (fallback || 0);
  }

  function reversePermutation(n) {
    n = Math.max(0, finiteInt(n, 0));
    var perm = new Array(n);
    for (var i = 0; i < n; i++) perm[i] = n - 1 - i;
    return perm;
  }

  function remapIndexFlags(flags, perm) {
    flags = flags || {};
    perm = perm || [];
    var n = perm.length;
    if (n < 2) return {};
    var next = {};
    Object.keys(flags).forEach(function (key) {
      if (!flags[key]) return;
      var src = finiteInt(key, -1);
      if (src < 0 || src >= n) return;
      for (var di = 0; di < n; di++) {
        if (perm[di] === src) {
          next[di] = true;
          break;
        }
      }
    });
    return next;
  }

  /**
   * @param {object} [opts]
   * @param {object} [opts.hooks] Page callbacks (syncVolumeChrome, …).
   */
  function DirectTraceUiState(opts) {
    opts = opts || {};
    this.hooks = opts.hooks || {};
    this._staleAt = {};
    this._detachedAt = {};
    this._dragInvalidated = {};
  }

  DirectTraceUiState.prototype.setHooks = function (hooks) {
    this.hooks = hooks || {};
    return this;
  };

  /* ── Invalidation maps ─────────────────────────────────────────────── */

  DirectTraceUiState.prototype.staleAt = function () {
    return this._staleAt;
  };

  DirectTraceUiState.prototype.detachedAt = function () {
    return this._detachedAt;
  };

  DirectTraceUiState.prototype.dragInvalidated = function () {
    return this._dragInvalidated;
  };

  DirectTraceUiState.prototype.isStale = function (index) {
    index = finiteInt(index, -1);
    return index >= 0 && !!this._staleAt[index];
  };

  DirectTraceUiState.prototype.isDetached = function (index) {
    index = finiteInt(index, -1);
    return index >= 0 && !!this._detachedAt[index];
  };

  DirectTraceUiState.prototype.isDragInvalidated = function (index) {
    index = finiteInt(index, -1);
    return index >= 0 && !!this._dragInvalidated[index];
  };

  DirectTraceUiState.prototype.markStale = function (index) {
    index = finiteInt(index, -1);
    if (index < 0) return this;
    this._staleAt[index] = true;
    return this;
  };

  DirectTraceUiState.prototype.markDetached = function (index) {
    index = finiteInt(index, -1);
    if (index < 0) return this;
    // Label-only: slider shows traj-vol text instead of catalog particle
    // labels. Does not block tick readiness after Decode+Render.
    this._detachedAt[index] = true;
    return this;
  };

  DirectTraceUiState.prototype.markDragInvalidated = function (index) {
    index = finiteInt(index, -1);
    if (index < 0) return this;
    this._dragInvalidated[index] = true;
    return this;
  };

  DirectTraceUiState.prototype.clearDragInvalidated = function () {
    this._dragInvalidated = {};
    return this;
  };

  DirectTraceUiState.prototype.clearDetached = function () {
    this._detachedAt = {};
    return this;
  };

  DirectTraceUiState.prototype.clearStale = function () {
    this._staleAt = {};
    return this;
  };

  DirectTraceUiState.prototype.clearAllInvalidation = function () {
    this._staleAt = {};
    this._detachedAt = {};
    this._dragInvalidated = {};
    return this;
  };

  DirectTraceUiState.prototype.clearInvalidationAtIndices = function (indices) {
    if (!indices || !indices.length) return this;
    for (var i = 0; i < indices.length; i++) {
      var idx = finiteInt(indices[i], -1);
      if (idx < 0) continue;
      delete this._staleAt[idx];
      delete this._dragInvalidated[idx];
      // Keep detachedAt — freely moved endpoints keep traj-vol tick labels
      // through Decode+Render until path reset / direct-trace handoff.
    }
    return this;
  };

  DirectTraceUiState.prototype.replaceStaleMap = function (map) {
    this._staleAt = Object.assign({}, map || {});
    return this;
  };

  /**
   * Remap stale / detached / drag flags through a reverse permutation of length n.
   */
  DirectTraceUiState.prototype.reverseInvalidationMaps = function (n) {
    n = Math.max(0, finiteInt(n, 0));
    if (n < 2) return this;
    var perm = reversePermutation(n);
    this._staleAt = remapIndexFlags(this._staleAt, perm);
    this._detachedAt = remapIndexFlags(this._detachedAt, perm);
    this._dragInvalidated = remapIndexFlags(this._dragInvalidated, perm);
    return this;
  };

  /**
   * Remap stale / detached / drag flags through a densify/shrink src→dest map
   * (from VolumeState.lastRematchMap). Unmapped sources are dropped.
   */
  DirectTraceUiState.prototype.remapFlagsBySlotMap = function (srcToDest) {
    srcToDest = srcToDest || {};
    function remap(flags) {
      var next = {};
      Object.keys(flags || {}).forEach(function (key) {
        if (!flags[key]) return;
        var dest = srcToDest[key];
        if (dest == null) dest = srcToDest[String(key)];
        if (dest == null) return;
        dest = finiteInt(dest, -1);
        if (dest < 0) return;
        next[dest] = true;
      });
      return next;
    }
    this._staleAt = remap(this._staleAt);
    this._detachedAt = remap(this._detachedAt);
    this._dragInvalidated = remap(this._dragInvalidated);
    return this;
  };

  DirectTraceUiState.prototype.hasAnyStale = function () {
    var keys = Object.keys(this._staleAt);
    for (var i = 0; i < keys.length; i++) {
      if (this._staleAt[keys[i]]) return true;
    }
    return false;
  };

  /* ── Path ↔ controls alignment ─────────────────────────────────────── */

  /**
   * Resolve the authoritative live path length for direct-trace UI.
   * @param {object} ctx
   * @param {Array|null} ctx.editableTrajXY
   * @param {number} [ctx.controlNPoints]
   * @param {number} [ctx.latentCount]
   * @param {number} [ctx.displayCount]
   */
  DirectTraceUiState.prototype.resolvePathLength = function (ctx) {
    ctx = ctx || {};
    var editableN = (ctx.editableTrajXY && ctx.editableTrajXY.length >= 2)
      ? ctx.editableTrajXY.length
      : 0;
    return Math.max(
      2,
      editableN,
      finiteInt(ctx.controlNPoints, 0),
      finiteInt(ctx.latentCount, 0),
      finiteInt(ctx.displayCount, 0)
    );
  };

  /**
   * Whether volume display length disagrees with the live path.
   */
  DirectTraceUiState.prototype.pathDisplayLengthMismatch = function (pathN, displayN) {
    pathN = finiteInt(pathN, 0);
    displayN = finiteInt(displayN, 0);
    return pathN >= 2 && displayN >= 1 && pathN !== displayN;
  };

  /**
   * Align Trajectory-controls nPoints to the live editable path when the path
   * owns length (e.g. after double-click insert).
   *
   * Only bumps the control upward. Shrinking must come from Trajectory-controls
   * (nPoints); pulling the control down here raced applyResampled and left the
   * volume slider on the old tick count.
   *
   * @returns {{ pathN: number, controlUpdated: boolean }}
   */
  DirectTraceUiState.prototype.alignControlsToLivePath = function (ctx) {
    ctx = ctx || {};
    var pathN = this.resolvePathLength(ctx);
    var controlN = finiteInt(ctx.controlNPoints, 0);
    var controlUpdated = false;
    var editableN = (ctx.editableTrajXY && ctx.editableTrajXY.length >= 2)
      ? ctx.editableTrajXY.length
      : 0;
    // Path-owned growth (insert): push the control up to match.
    if (editableN >= 2 && editableN > controlN && typeof ctx.setControlNPoints === "function") {
      ctx.setControlNPoints(editableN);
      controlUpdated = true;
      pathN = editableN;
    } else if (controlN >= 2 && controlN >= editableN) {
      // Control-owned length (nPoints change): path follows the control.
      pathN = Math.max(pathN, controlN);
    }
    return { pathN: pathN, controlUpdated: controlUpdated };
  };

  /**
   * Decide whether →direct should rebuild at live editable length vs control.
   */
  DirectTraceUiState.prototype.rebuildPointCountForDirectMode = function (ctx) {
    ctx = ctx || {};
    var editableN = (ctx.editableTrajXY && ctx.editableTrajXY.length >= 2)
      ? ctx.editableTrajXY.length
      : 0;
    var controlN = Math.max(2, finiteInt(ctx.controlNPoints, 0));
    return Math.max(controlN, editableN);
  };

  /**
   * Orchestrate UI sync after a direct-trace path mutation.
   * Hooks (all optional):
   *   alignControls(ctx) → { pathN, controlUpdated }
   *   ensureEndpointIds()
   *   syncVolumeChrome(opts)
   *   syncTickLabels()
   *   updateDebtLabels()
   *   syncOverlay()
   *   syncPendingOverlay()
   */
  DirectTraceUiState.prototype.afterPathMutation = function (opts) {
    opts = opts || {};
    var hooks = this.hooks || {};
    var align = { pathN: Math.max(0, finiteInt(opts.pathN, 0)), controlUpdated: false };
    // nPoints / control-owned mutations pass alignControls: false so we never
    // rewrite Trajectory-controls from a stale editable polyline.
    if (opts.alignControls !== false && typeof hooks.alignControls === "function") {
      try { align = hooks.alignControls(opts) || align; } catch (err) { /* ignore */ }
    }
    if (opts.ensureEndpointIds !== false && typeof hooks.ensureEndpointIds === "function") {
      try { hooks.ensureEndpointIds(); } catch (err2) { /* ignore */ }
    }
    if (typeof hooks.syncVolumeChrome === "function") {
      try {
        hooks.syncVolumeChrome({
          forceExpand: opts.forceExpand !== false,
          pathN: align.pathN || opts.pathN
        });
      } catch (err3) { /* ignore */ }
    }
    if (opts.syncTicks !== false && typeof hooks.syncTickLabels === "function") {
      try { hooks.syncTickLabels(); } catch (err4) { /* ignore */ }
    }
    if (typeof hooks.updateDebtLabels === "function") {
      try { hooks.updateDebtLabels(); } catch (err5) { /* ignore */ }
    }
    if (opts.redraw !== false && typeof hooks.syncOverlay === "function") {
      try { hooks.syncOverlay(); } catch (err6) { /* ignore */ }
    }
    if (typeof hooks.syncPendingOverlay === "function") {
      try { hooks.syncPendingOverlay(); } catch (err7) { /* ignore */ }
    }
    return align;
  };

  DirectTraceUiState.reversePermutation = reversePermutation;
  DirectTraceUiState.remapIndexFlags = remapIndexFlags;

  global.CryoDirectTraceUiState = DirectTraceUiState;
})(typeof window !== "undefined" ? window : this);
