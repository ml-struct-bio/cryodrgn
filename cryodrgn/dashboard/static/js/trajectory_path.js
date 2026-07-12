/**
 * Trajectory path model for the cryoDRGN trajectory creator.
 *
 * Path geometry and mode-specific mutation live here. The dashboard page
 * supplies a mutable ``store`` adapter that mirrors the legacy IIFE globals
 * (``anchorIndicesActive``, ``anchorTrajXY``, ``editableTrajXY``, …) so existing
 * call sites keep working while new logic goes through ``currentPath``.
 *
 * Hierarchy:
 *   TrajectoryPath
 *     ├── DirectTracePath  — latent-space polyline (start/end + editable samples)
 *     └── WaypointPath     — particle / catalog volume waypoints
 */
(function (global) {
  "use strict";

  function asFiniteNumber(v) {
    var n = Number(v);
    return Number.isFinite(n) ? n : NaN;
  }

  function cloneXY(pt) {
    if (!pt || pt.length < 2) return null;
    var x = asFiniteNumber(pt[0]);
    var y = asFiniteNumber(pt[1]);
    if (!Number.isFinite(x) || !Number.isFinite(y)) return null;
    return [x, y];
  }

  function cloneXYList(list) {
    if (!list || !list.length) return null;
    var out = [];
    for (var i = 0; i < list.length; i++) {
      var pt = cloneXY(list[i]);
      if (!pt) return null;
      out.push(pt);
    }
    return out;
  }

  function cloneRowList(list) {
    if (!list || !list.length) return null;
    return list.map(function (r) { return Number(r); }).filter(function (r) {
      return Number.isFinite(r);
    });
  }

  /**
   * @param {object} store  Mutable adapter with getters/setters for path fields.
   * @param {object} [hooks] Optional UI/API callbacks supplied by the page.
   */
  function TrajectoryPath(store, hooks) {
    this.store = store || {};
    this.hooks = hooks || {};
  }

  TrajectoryPath.prototype.kind = "base";

  TrajectoryPath.prototype.getMode = function () {
    return String(this.store.trajectoryMode || "manual");
  };

  TrajectoryPath.prototype.setMode = function (mode) {
    this.store.trajectoryMode = String(mode || "manual");
  };

  TrajectoryPath.prototype.getAnchorRows = function () {
    return this.store.anchorIndicesActive || null;
  };

  TrajectoryPath.prototype.setAnchorRows = function (rows) {
    this.store.anchorIndicesActive = rows && rows.length ? cloneRowList(rows) : null;
  };

  TrajectoryPath.prototype.getPathXY = function () {
    var anchors = this.store.anchorTrajXY;
    if (anchors && anchors.length >= 2) return anchors;
    var editable = this.store.editableTrajXY;
    if (editable && editable.length >= 2) return editable;
    return null;
  };

  TrajectoryPath.prototype.setAnchorPathXY = function (xy) {
    this.store.anchorTrajXY = xy && xy.length ? cloneXYList(xy) : null;
  };

  TrajectoryPath.prototype.getSamplePlotRows = function () {
    return this.store.trajPlotRows || null;
  };

  TrajectoryPath.prototype.setSamplePlotRows = function (rows) {
    this.store.trajPlotRows = rows && rows.length ? cloneRowList(rows) : null;
  };

  TrajectoryPath.prototype.getEditableXY = function () {
    return this.store.editableTrajXY || null;
  };

  TrajectoryPath.prototype.setEditableXY = function (xy) {
    this.store.editableTrajXY = xy && xy.length ? cloneXYList(xy) : null;
  };

  TrajectoryPath.prototype.getStart = function () {
    return cloneXY(this.store.startXY);
  };

  TrajectoryPath.prototype.getEnd = function () {
    return cloneXY(this.store.endXY);
  };

  TrajectoryPath.prototype.setEndpoints = function (start, end) {
    this.store.startXY = cloneXY(start);
    this.store.endXY = cloneXY(end);
  };

  TrajectoryPath.prototype.getInterpolatedCount = function () {
    return Math.max(0, parseInt(this.store.manualInterpTrajectoryCount, 10) || 0);
  };

  TrajectoryPath.prototype.setInterpolatedCount = function (n) {
    this.store.manualInterpTrajectoryCount = Math.max(0, parseInt(n, 10) || 0);
  };

  TrajectoryPath.prototype.hasAnchors = function () {
    var rows = this.getAnchorRows();
    return !!(rows && rows.length >= 2);
  };

  TrajectoryPath.prototype.anchorCount = function () {
    var rows = this.getAnchorRows();
    return rows && rows.length ? rows.length : 0;
  };

  TrajectoryPath.prototype.sampleCount = function () {
    var samples = this.getSamplePlotRows();
    if (samples && samples.length) return samples.length;
    var xy = this.getPathXY();
    if (xy && xy.length) return xy.length;
    var n = this.getInterpolatedCount();
    return n >= 2 ? n : 0;
  };

  /**
   * True when the path has samples denser than its defining anchors
   * (snap / graph / interpolated interiors).
   */
  Object.defineProperty(TrajectoryPath.prototype, "hasInteriorSamples", {
    configurable: true,
    enumerable: true,
    get: function () {
      var nAnchors = this.anchorCount();
      if (nAnchors < 2) {
        // Endpoint-only editable polylines are not "interior relative to anchors".
        return false;
      }
      var samples = this.getSamplePlotRows();
      if (samples && samples.length > nAnchors) return true;
      var xy = this.store.anchorTrajXY;
      if (xy && xy.length > nAnchors) return true;
      if (this.getInterpolatedCount() > nAnchors) return true;
      return false;
    }
  });

  TrajectoryPath.prototype.isWaypointMode = function () {
    return this.getMode() === "manual";
  };

  TrajectoryPath.prototype.isDirectTraceMode = function () {
    var mode = this.getMode();
    return mode === "direct" || mode === "nearest";
  };

  TrajectoryPath.prototype.endpoints = function () {
    var xy = this.getPathXY();
    if (xy && xy.length >= 2) {
      return { start: xy[0].slice(), end: xy[xy.length - 1].slice() };
    }
    var start = this.getStart();
    var end = this.getEnd();
    if (start && end) return { start: start, end: end };
    return null;
  };

  TrajectoryPath.prototype.clearGeometry = function () {
    this.setAnchorRows(null);
    this.setAnchorPathXY(null);
    this.setSamplePlotRows(null);
    this.setEditableXY(null);
    this.setInterpolatedCount(0);
    this.setEndpoints(null, null);
  };

  TrajectoryPath.prototype.applyCoordsPayload = function (payload) {
    payload = payload || {};
    if (payload.traj_rows && payload.traj_rows.length) {
      this.setSamplePlotRows(payload.traj_rows);
    } else if (this.getMode() !== "nearest") {
      this.setSamplePlotRows(null);
    }
    if (payload.traj_marker_colors && payload.traj_marker_colors.length) {
      this.store.trajMarkerColors = payload.traj_marker_colors.map(function (c) {
        return (c == null || c === "") ? null : String(c);
      });
    }
    if (this.hasAnchors() && payload.traj_xy) {
      this.setAnchorPathXY(payload.traj_xy);
      if (payload.traj_rows && payload.traj_rows.length >= 2) {
        this.setSamplePlotRows(payload.traj_rows);
        if (payload.traj_xy.length !== payload.traj_rows.length
            && typeof this.hooks.trajXYFromPlotRows === "function") {
          var rebuilt = this.hooks.trajXYFromPlotRows(payload.traj_rows);
          if (rebuilt && rebuilt.length === payload.traj_rows.length) {
            this.setAnchorPathXY(rebuilt);
          }
        }
      }
      this.setInterpolatedCount(this.store.anchorTrajXY ? this.store.anchorTrajXY.length : 0);
      this.setEditableXY(null);
    } else if (payload.traj_xy && payload.traj_xy.length) {
      this.setEditableXY(payload.traj_xy);
      if (this.store.editableTrajXY && this.store.editableTrajXY.length >= 2) {
        this.setEndpoints(
          this.store.editableTrajXY[0],
          this.store.editableTrajXY[this.store.editableTrajXY.length - 1]
        );
      }
      this.setInterpolatedCount(this.store.editableTrajXY ? this.store.editableTrajXY.length : 0);
    }
    if (payload.anchor_indices && payload.anchor_indices.length >= 2) {
      this.setAnchorRows(payload.anchor_indices);
    }
    if (payload.kmeans_volume_ids && payload.kmeans_volume_ids.length >= 2) {
      this.store.anchorKmeansVolumeIds = payload.kmeans_volume_ids.slice();
    }
  };

  TrajectoryPath.prototype.toCoordsRequestFields = function () {
    if (this.hasAnchors()) {
      return {
        anchor_indices: this.getAnchorRows().slice(),
        has_anchors: true
      };
    }
    var start = this.getStart();
    var end = this.getEnd();
    var fields = { has_anchors: false };
    if (start && end) {
      fields.start = start;
      fields.end = end;
    }
    var editable = this.getEditableXY();
    if (editable && editable.length >= 2) fields.traj_xy = editable.slice();
    return fields;
  };

  TrajectoryPath.prototype.reverse = function () {
    var rows = this.getAnchorRows();
    if (rows && rows.length >= 2) this.setAnchorRows(rows.slice().reverse());
    var axy = this.store.anchorTrajXY;
    if (axy && axy.length >= 2) this.setAnchorPathXY(axy.slice().reverse());
    var samples = this.getSamplePlotRows();
    if (samples && samples.length >= 2) this.setSamplePlotRows(samples.slice().reverse());
    var editable = this.getEditableXY();
    if (editable && editable.length >= 2) {
      this.setEditableXY(editable.slice().reverse());
      this.setEndpoints(
        this.store.editableTrajXY[0],
        this.store.editableTrajXY[this.store.editableTrajXY.length - 1]
      );
    } else {
      var start = this.getStart();
      var end = this.getEnd();
      if (start && end) this.setEndpoints(end, start);
    }
    var ids = this.store.anchorKmeansVolumeIds;
    if (ids && ids.length >= 2) {
      this.store.anchorKmeansVolumeIds = ids.slice().reverse();
    }
  };

  /**
   * Double-click / pick entry point. Subclasses override.
   * @returns {{ok:boolean, reason?:string}}
   */
  TrajectoryPath.prototype.addPointFromScatter = function (/* detail */) {
    return { ok: false, reason: "unsupported" };
  };

  /**
   * Rebuild anchors from currently selected catalog volumes + active custom rows.
   * Used in waypoint mode and when direct-trace mode is driven by particle sets.
   */
  TrajectoryPath.prototype.rebuildFromSelection = function (opts) {
    opts = opts || {};
    if (this.kind === "waypoint" && this.interpolationArmed && this.interpolationArmed() && !opts.force) {
      return { ok: true, preserved: true };
    }
    var rowXY = this.hooks.rowXYForPlotRow;
    var markersByVolId = this.store.manualMarkersByVolId || {};
    var selected = this.store.manualSelectedVolIds || [];
    var activeCustom = this.store.manualActiveCustomPlotRows || [];
    var rows = [];
    var xy = [];
    var i;
    for (i = 0; i < selected.length; i++) {
      var marker = markersByVolId[selected[i]];
      if (!marker || marker.plot_row == null) continue;
      var pr = Number(marker.plot_row);
      if (!Number.isFinite(pr) || rows.indexOf(pr) >= 0) continue;
      rows.push(pr);
      if (typeof rowXY === "function") {
        var pt = rowXY(pr);
        if (pt) xy.push(pt.slice());
      }
    }
    for (i = 0; i < activeCustom.length; i++) {
      var crow = Number(activeCustom[i]);
      if (!Number.isFinite(crow) || rows.indexOf(crow) >= 0) continue;
      rows.push(crow);
      if (typeof rowXY === "function") {
        var cpt = rowXY(crow);
        if (cpt) xy.push(cpt.slice());
      }
    }
    if (rows.length < 2) {
      this.setAnchorRows(null);
      this.setAnchorPathXY(null);
      this.setSamplePlotRows(null);
      this.setEditableXY(null);
      this.setInterpolatedCount(0);
      this.store.anchorKmeansVolumeIds = [];
      this.store.manualTrajectoryVolumeIds = [];
      return { ok: true, cleared: true, rows: rows };
    }
    this.setAnchorRows(rows);
    if (xy.length === rows.length) {
      this.setAnchorPathXY(xy);
    } else if (typeof this.hooks.trajXYFromPlotRows === "function") {
      var rebuilt = this.hooks.trajXYFromPlotRows(rows);
      this.setAnchorPathXY(rebuilt && rebuilt.length === rows.length ? rebuilt : null);
    } else {
      this.setAnchorPathXY(null);
    }
    this.setEditableXY(null);
    this.setSamplePlotRows(null);
    this.setInterpolatedCount(0);
    this.store.manualTrajectoryVolumeIds = [];
    return { ok: true, rows: rows.slice(), xy: this.store.anchorTrajXY };
  };

  /* ── Direct latent-space trace ─────────────────────────────────────── */

  function DirectTracePath(store, hooks) {
    TrajectoryPath.call(this, store, hooks);
  }

  DirectTracePath.prototype = Object.create(TrajectoryPath.prototype);
  DirectTracePath.prototype.constructor = DirectTracePath;
  DirectTracePath.prototype.kind = "direct";

  DirectTracePath.prototype.isNearest = function () {
    return this.getMode() === "nearest";
  };

  DirectTracePath.prototype.ensureEditablePolyline = function (nPoints) {
    nPoints = Math.max(2, parseInt(nPoints, 10) || 2);
    var editable = this.getEditableXY();
    if (editable && editable.length >= 2) return editable;
    var start = this.getStart();
    var end = this.getEnd();
    if (!start || !end) return null;
    var out = [];
    for (var i = 0; i < nPoints; i++) {
      var t = (nPoints <= 1) ? 0 : (i / (nPoints - 1));
      out.push([
        start[0] + t * (end[0] - start[0]),
        start[1] + t * (end[1] - start[1])
      ]);
    }
    this.setEditableXY(out);
    return out;
  };

  DirectTracePath.prototype.buildLinearPath = function (nPoints) {
    var start = this.getStart();
    var end = this.getEnd();
    if (!start || !end) return false;
    nPoints = Math.max(2, parseInt(nPoints, 10) || 2);
    var out = [];
    for (var i = 0; i < nPoints; i++) {
      var t = (nPoints <= 1) ? 0 : (i / (nPoints - 1));
      out.push([
        start[0] + t * (end[0] - start[0]),
        start[1] + t * (end[1] - start[1])
      ]);
    }
    this.setEditableXY(out);
    this.setAnchorRows(null);
    this.setAnchorPathXY(null);
    this.setSamplePlotRows(null);
    this.setInterpolatedCount(out.length);
    return true;
  };

  DirectTracePath.prototype.setEndpointVolumeIds = function (ids) {
    this.store.directEndpointVolumeIds = (ids && ids.length >= 2) ? ids.slice() : [];
  };

  DirectTracePath.prototype.getEndpointVolumeIds = function () {
    return this.store.directEndpointVolumeIds || [];
  };

  /**
   * Insert a free latent-space coordinate (double-click in tracing mode).
   * - No endpoints yet → becomes start
   * - Start only → becomes end and builds a 2-point path
   * - Existing polyline → inserts before the end endpoint (or at nearest segment)
   */
  DirectTracePath.prototype.addCoordinate = function (xy, opts) {
    opts = opts || {};
    var pt = cloneXY(xy);
    if (!pt) return { ok: false, reason: "invalid_xy" };
    var start = this.getStart();
    var end = this.getEnd();
    if (!start) {
      this.setEndpoints(pt, null);
      return { ok: true, action: "set_start" };
    }
    if (!end) {
      this.setEndpoints(start, pt);
      var n = Math.max(2, parseInt(opts.nPoints, 10) || 2);
      this.buildLinearPath(n);
      return { ok: true, action: "set_end" };
    }
    // Anchored catalog path in direct mode: convert to editable polyline first.
    if (this.hasAnchors() && !this.getEditableXY()) {
      var axy = this.store.anchorTrajXY;
      if (axy && axy.length >= 2) this.setEditableXY(axy);
      this.setAnchorRows(null);
      this.setSamplePlotRows(null);
    }
    var poly = this.ensureEditablePolyline(opts.nPoints || (this.getEditableXY() || []).length || 4);
    if (!poly || poly.length < 2) return { ok: false, reason: "no_polyline" };
    var insertAt = poly.length - 1;
    if (opts.insertIndex != null && Number.isFinite(Number(opts.insertIndex))) {
      insertAt = Math.max(1, Math.min(poly.length - 1, Math.floor(Number(opts.insertIndex))));
    } else if (opts.nearestSegment) {
      insertAt = nearestSegmentInsertIndex(poly, pt);
    }
    poly = poly.slice();
    poly.splice(insertAt, 0, pt);
    this.setEditableXY(poly);
    this.setEndpoints(poly[0], poly[poly.length - 1]);
    this.setInterpolatedCount(poly.length);
    this.setEndpointVolumeIds([]);
    return { ok: true, action: "insert", index: insertAt };
  };

  DirectTracePath.prototype.addPointFromScatter = function (detail) {
    detail = detail || {};
    return this.addCoordinate(detail.xy, {
      nPoints: detail.nPoints,
      nearestSegment: detail.nearestSegment !== false
    });
  };

  DirectTracePath.prototype.onEndpointDrag = function (which, xy) {
    var pt = cloneXY(xy);
    if (!pt) return false;
    var editable = this.getEditableXY();
    if (which === "start") {
      this.store.startXY = pt;
      if (editable && editable.length) editable[0] = pt.slice();
    } else if (which === "end") {
      this.store.endXY = pt;
      if (editable && editable.length) editable[editable.length - 1] = pt.slice();
    } else {
      return false;
    }
    if (editable) this.setEditableXY(editable);
    return true;
  };

  /* ── Waypoint / particle-set path ──────────────────────────────────── */

  function WaypointPath(store, hooks) {
    TrajectoryPath.call(this, store, hooks);
  }

  WaypointPath.prototype = Object.create(TrajectoryPath.prototype);
  WaypointPath.prototype.constructor = WaypointPath;
  WaypointPath.prototype.kind = "waypoint";

  WaypointPath.prototype.getSelectedVolumeIds = function () {
    return this.store.manualSelectedVolIds || [];
  };

  WaypointPath.prototype.setSelectedVolumeIds = function (ids) {
    this.store.manualSelectedVolIds = (ids || []).map(String);
  };

  WaypointPath.prototype.interpolationArmed = function () {
    return !!(this.store.manualVolumeSnapActive
      || this.store.manualParticleSnapPathActive
      || this.store.graphTraversalActive);
  };

  WaypointPath.prototype.armInterpolation = function (kind) {
    var volume = kind === "volume";
    this.store.manualVolumeSnapActive = volume;
    this.store.manualParticleSnapPathActive = volume;
    this.store.graphTraversalActive = !volume;
    this.setEditableXY(null);
  };

  WaypointPath.prototype.disarmInterpolation = function () {
    this.store.manualVolumeSnapActive = false;
    this.store.manualParticleSnapPathActive = false;
    this.store.graphTraversalActive = false;
    this.setInterpolatedCount(0);
    this.setSamplePlotRows(null);
    // Clear dense geometry so scatter/preview rebuild from the current selection
    // rather than treating a stale longer polyline as still-valid interiors.
    this.setAnchorPathXY(null);
    this.store.lastLatentTrajectoryPointCount = 0;
  };

  WaypointPath.prototype.anchorTraversalMode = function () {
    return this.store.graphTraversalActive ? "graph" : "direct";
  };

  /**
   * Add the hover / clicked particle as a custom waypoint (double-click in
   * choosing-waypoints mode).
   */
  WaypointPath.prototype.addParticleRow = function (plotRow, opts) {
    opts = opts || {};
    plotRow = Number(plotRow);
    if (!Number.isFinite(plotRow) || plotRow < 0) {
      return { ok: false, reason: "invalid_row" };
    }
    if (typeof this.hooks.registerCustomPlotRows === "function") {
      var changed = this.hooks.registerCustomPlotRows([plotRow], opts.activate !== false);
      return { ok: !!changed, action: "register_custom", plotRow: plotRow };
    }
    var custom = this.store.manualCustomPlotRows || [];
    var active = this.store.manualActiveCustomPlotRows || [];
    if (custom.indexOf(plotRow) < 0) {
      custom = custom.concat([plotRow]).sort(function (a, b) { return a - b; });
      this.store.manualCustomPlotRows = custom;
    }
    if (active.indexOf(plotRow) < 0) {
      active = active.concat([plotRow]).sort(function (a, b) { return a - b; });
      this.store.manualActiveCustomPlotRows = active;
    }
    this.disarmInterpolation();
    return this.rebuildFromSelection();
  };

  WaypointPath.prototype.addPointFromScatter = function (detail) {
    detail = detail || {};
    var row = detail.plotRow;
    if (row == null && detail.hoverPlotRow != null) row = detail.hoverPlotRow;
    return this.addParticleRow(row, detail);
  };

  WaypointPath.prototype.setAnchorsFromIndices = function (indices) {
    var rows = cloneRowList(indices);
    if (!rows || rows.length < 2) {
      this.setAnchorRows(null);
      return { ok: false, reason: "need_two_anchors" };
    }
    this.setAnchorRows(rows);
    this.setEditableXY(null);
    this.setSamplePlotRows(null);
    return { ok: true, rows: rows };
  };

  /* ── helpers ───────────────────────────────────────────────────────── */

  function nearestSegmentInsertIndex(poly, pt) {
    var bestI = poly.length - 1;
    var bestD = Infinity;
    for (var i = 0; i < poly.length - 1; i++) {
      var a = poly[i];
      var b = poly[i + 1];
      var d = pointSegmentDistanceSq(pt, a, b);
      if (d < bestD) {
        bestD = d;
        bestI = i + 1;
      }
    }
    return Math.max(1, Math.min(poly.length - 1, bestI));
  }

  function pointSegmentDistanceSq(p, a, b) {
    var abx = b[0] - a[0];
    var aby = b[1] - a[1];
    var apx = p[0] - a[0];
    var apy = p[1] - a[1];
    var ab2 = abx * abx + aby * aby;
    var t = ab2 > 0 ? (apx * abx + apy * aby) / ab2 : 0;
    t = Math.max(0, Math.min(1, t));
    var cx = a[0] + t * abx;
    var cy = a[1] + t * aby;
    var dx = p[0] - cx;
    var dy = p[1] - cy;
    return dx * dx + dy * dy;
  }

  /**
   * Build the appropriate path instance for the current store mode.
   */
  function createForStore(store, hooks) {
    store = store || {};
    var mode = String(store.trajectoryMode || "manual");
    if (mode === "manual") return new WaypointPath(store, hooks);
    return new DirectTracePath(store, hooks);
  }

  /**
   * Ensure ``host.currentPath`` matches ``host`` store mode, recreating when
   * the mode family changes (waypoint ↔ direct).
   */
  function ensureCurrentPath(host, hooks) {
    if (!host) return null;
    var mode = String(host.trajectoryMode || "manual");
    var wantWaypoint = mode === "manual";
    var cur = host.currentPath;
    if (cur) {
      var isWaypoint = cur.kind === "waypoint";
      if (wantWaypoint === isWaypoint) {
        cur.hooks = hooks || cur.hooks || {};
        return cur;
      }
    }
    host.currentPath = createForStore(host, hooks || (cur && cur.hooks) || {});
    return host.currentPath;
  }

  global.CryoTrajectoryPath = {
    TrajectoryPath: TrajectoryPath,
    DirectTracePath: DirectTracePath,
    WaypointPath: WaypointPath,
    createForStore: createForStore,
    ensureCurrentPath: ensureCurrentPath
  };
})(typeof window !== "undefined" ? window : this);
