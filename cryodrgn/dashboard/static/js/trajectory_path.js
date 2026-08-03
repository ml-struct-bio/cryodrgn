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

  function linearXYFromEndpoints(start, end, nPoints) {
    start = cloneXY(start);
    end = cloneXY(end);
    if (!start || !end) return null;
    nPoints = Math.max(2, parseInt(nPoints, 10) || 2);
    var out = [];
    for (var i = 0; i < nPoints; i++) {
      var t = (nPoints <= 1) ? 0 : (i / (nPoints - 1));
      out.push([
        start[0] + t * (end[0] - start[0]),
        start[1] + t * (end[1] - start[1])
      ]);
    }
    return out;
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

  TrajectoryPath.prototype.applyCoordsPayload = function (payload, opts) {
    payload = payload || {};
    opts = opts || {};
    var preserveWaypointGeometry = !!opts.preserveWaypointGeometry;
    if (payload.traj_marker_colors && payload.traj_marker_colors.length) {
      this.store.trajMarkerColors = payload.traj_marker_colors.map(function (c) {
        return (c == null || c === "") ? null : String(c);
      });
    }
    if (preserveWaypointGeometry) {
      if (payload.kmeans_volume_ids && payload.kmeans_volume_ids.length >= 2) {
        this.store.anchorKmeansVolumeIds = payload.kmeans_volume_ids.slice();
      }
      return;
    }
    if (payload.traj_rows && payload.traj_rows.length) {
      this.setSamplePlotRows(payload.traj_rows);
    } else if (this.getMode() !== "nearest") {
      this.setSamplePlotRows(null);
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

  TrajectoryPath.prototype.reverse = function () {
    var rows = this.getAnchorRows();
    if (rows && rows.length >= 2) this.setAnchorRows(rows.slice().reverse());
    var axy = this.store.anchorTrajXY;
    if (axy && axy.length >= 2) this.setAnchorPathXY(axy.slice().reverse());
    var samples = this.getSamplePlotRows();
    if (samples && samples.length >= 2) this.setSamplePlotRows(samples.slice().reverse());
    var colors = this.store.trajMarkerColors;
    if (colors && colors.length >= 2) {
      this.store.trajMarkerColors = colors.slice().reverse();
    }
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
    var directIds = this.store.directEndpointVolumeIds;
    if (directIds && directIds.length >= 2) {
      this.store.directEndpointVolumeIds = directIds.slice().reverse();
    }
  };

  function rowMultisetKey(rows) {
    return (rows || []).map(Number).filter(function (r) {
      return Number.isFinite(r);
    }).slice().sort(function (a, b) { return a - b; }).join(",");
  }

  function sameRowMultiset(a, b) {
    if (!a || !b || a.length !== b.length) return false;
    return rowMultisetKey(a) === rowMultisetKey(b);
  }

  TrajectoryPath.prototype.anchorRowsFromSelection = function (opts) {
    opts = opts || {};
    var includeCustom = opts.includeCustom !== false;
    var rowXY = this.hooks.rowXYForPlotRow;
    var markersByVolId = this.store.manualMarkersByVolId || {};
    var selected = this.store.manualSelectedVolIds || [];
    var activeCustom = this.store.manualActiveCustomPlotRows || [];
    var rows = [];
    var xy = [];
    var volumeRows = [];
    var i;

    function markerXY(marker, plotRow) {
      if (!marker) return null;
      if (marker.xy && marker.xy.length === 2) {
        var mx = Number(marker.xy[0]);
        var my = Number(marker.xy[1]);
        if (Number.isFinite(mx) && Number.isFinite(my)) return [mx, my];
      }
      if (!Number.isFinite(plotRow) || typeof rowXY !== "function") return null;
      var pt = rowXY(plotRow);
      return pt ? pt.slice() : null;
    }

    function xyForRowList(rowList) {
      // One slot per selected catalog volume (queue by plot_row), matching the
      // catalog loop below. Reusing the first marker for a shared plot_row would
      // collapse distinct PC sample positions onto one particle (e.g. 3 under 4).
      var volsByRow = {};
      for (var si = 0; si < selected.length; si++) {
        var mk = markersByVolId[selected[si]];
        if (!mk || mk.plot_row == null) continue;
        var pr = Number(mk.plot_row);
        if (!Number.isFinite(pr)) continue;
        if (!volsByRow[pr]) volsByRow[pr] = [];
        volsByRow[pr].push(mk);
      }
      var out = [];
      for (var ri = 0; ri < rowList.length; ri++) {
        var want = Number(rowList[ri]);
        if (!Number.isFinite(want)) return null;
        var list = volsByRow[want];
        var found = null;
        if (list && list.length) {
          found = markerXY(list.shift(), want);
        } else if (typeof rowXY === "function") {
          found = rowXY(want);
          if (found) found = found.slice();
        }
        if (!found) return null;
        out.push(found);
      }
      return out;
    }

    for (i = 0; i < selected.length; i++) {
      var marker = markersByVolId[selected[i]];
      // Trajectory APIs require real dataset plot rows. Do not invent negative
      // indices for latent-only markers — those cannot densify via anchors.
      if (!marker || marker.plot_row == null) continue;
      var pr = Number(marker.plot_row);
      if (!Number.isFinite(pr)) continue;
      // Keep one path slot per selected volume, even when two volumes share a
      // nearest-particle plot_row (common for dense PC samples).
      rows.push(pr);
      volumeRows.push(pr);
      var pt = markerXY(marker, pr);
      if (pt) xy.push(pt);
    }
    if (includeCustom) {
      for (i = 0; i < activeCustom.length; i++) {
        var crow = Number(activeCustom[i]);
        if (!Number.isFinite(crow)
            || volumeRows.indexOf(crow) >= 0
            || rows.indexOf(crow) >= 0) {
          continue;
        }
        rows.push(crow);
        if (typeof rowXY === "function") {
          var cpt = rowXY(crow);
          if (cpt) xy.push(cpt.slice());
        }
      }
      // Exact / Inexact visit order interleaves catalog + Other waypoints.
      // Prefer the live anchor permutation when it covers the same multiset;
      // otherwise rebuildFromSelection would collapse back to catalog-then-custom.
      var active = this.store.anchorIndicesActive;
      if (opts.preferActiveOrder !== false
          && active && active.length >= 2
          && sameRowMultiset(active, rows)) {
        var ordered = active.map(Number).filter(function (r) {
          return Number.isFinite(r);
        });
        if (ordered.length === rows.length) {
          var orderedXY = xyForRowList(ordered);
          return {
            rows: ordered,
            xy: orderedXY && orderedXY.length === ordered.length ? orderedXY : null
          };
        }
      }
    }
    return {
      rows: rows,
      xy: xy.length === rows.length ? xy : null
    };
  };

  /**
   * Install an explicit waypoint order (e.g. Exact / Inexact tour) without
   * regrouping catalog volumes ahead of Other / random indices.
   */
  TrajectoryPath.prototype.setAnchorsFromOrderedRows = function (rows) {
    rows = cloneRowList(rows).map(Number).filter(function (r) {
      return Number.isFinite(r);
    });
    if (rows.length < 2) {
      this.setAnchorRows(null);
      this.setAnchorPathXY(null);
      this.setSamplePlotRows(null);
      this.setEditableXY(null);
      this.setInterpolatedCount(0);
      this.store.manualTrajectoryVolumeIds = [];
      return { ok: true, cleared: true, rows: [] };
    }
    this.setAnchorRows(rows);
    var xy = null;
    // Prefer catalog marker.xy (true PC / kmeans latent coords) when the page
    // supplies pathXYForOrderedRows — trajXYFromPlotRows alone snaps to nearest
    // particles and misplaces PC samples after Original-order rewrites.
    if (typeof this.hooks.pathXYForOrderedRows === "function") {
      xy = this.hooks.pathXYForOrderedRows(rows);
    }
    if ((!xy || xy.length !== rows.length)
        && typeof this.hooks.trajXYFromPlotRows === "function") {
      xy = this.hooks.trajXYFromPlotRows(rows);
    } else if ((!xy || xy.length !== rows.length)
        && typeof this.hooks.rowXYForPlotRow === "function") {
      xy = [];
      for (var i = 0; i < rows.length; i++) {
        var pt = this.hooks.rowXYForPlotRow(rows[i]);
        if (!pt) {
          xy = null;
          break;
        }
        xy.push(pt.slice());
      }
    }
    this.setAnchorPathXY(xy && xy.length === rows.length ? xy : null);
    this.setEditableXY(null);
    this.setSamplePlotRows(null);
    this.setInterpolatedCount(0);
    this.store.manualTrajectoryVolumeIds = [];
    return { ok: true, rows: rows.slice(), xy: this.store.anchorTrajXY };
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
    var selection = this.anchorRowsFromSelection({
      includeCustom: true,
      preferActiveOrder: opts.preferActiveOrder
    });
    var rows = selection.rows;
    var xy = selection.xy || [];
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
    var nextXY = null;
    if (xy.length === rows.length) {
      nextXY = xy;
    }
    // Prefer catalog marker.xy (true PC / kmeans latent coords) when available —
    // trajXYFromPlotRows alone snaps to nearest particles and also fails for
    // plot rows absent from the scatter subsample (common after Add-random).
    if ((!nextXY || nextXY.length !== rows.length)
        && typeof this.hooks.pathXYForOrderedRows === "function") {
      nextXY = this.hooks.pathXYForOrderedRows(rows);
    }
    if ((!nextXY || nextXY.length !== rows.length)
        && typeof this.hooks.trajXYFromPlotRows === "function") {
      nextXY = this.hooks.trajXYFromPlotRows(rows);
    }
    this.setAnchorPathXY(nextXY && nextXY.length === rows.length ? nextXY : null);
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

  DirectTracePath.prototype.ensureEditablePolyline = function (nPoints) {
    nPoints = Math.max(2, parseInt(nPoints, 10) || 2);
    var editable = this.getEditableXY();
    if (editable && editable.length >= 2) return editable;
    var out = linearXYFromEndpoints(this.getStart(), this.getEnd(), nPoints);
    if (!out) return null;
    this.setEditableXY(out);
    return out;
  };

  DirectTracePath.prototype.buildLinearPath = function (nPoints) {
    var out = linearXYFromEndpoints(this.getStart(), this.getEnd(), nPoints);
    if (!out) return false;
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
   * Densified volume-slot ids for Session / VolumeState: catalog endpoints at
   * 0 and pathN-1, null interiors (same role as WaypointPath densify anchors).
   */
  DirectTracePath.prototype.getVolumeSlotIds = function (pathN) {
    var endpoints = this.getEndpointVolumeIds();
    if (!endpoints || endpoints.length < 2) return [];
    pathN = Math.max(
      0,
      Math.floor(Number(pathN)) || 0,
      this.sampleCount() || 0
    );
    if (pathN < 2) pathN = 2;
    var out = new Array(pathN);
    for (var i = 0; i < pathN; i++) out[i] = null;
    out[0] = String(endpoints[0]);
    out[pathN - 1] = String(endpoints[1]);
    return out;
  };

  /**
   * Compact catalog ids for seedCompactCatalog (two endpoints).
   */
  DirectTracePath.prototype.getSelectedVolumeIds = function () {
    return this.getEndpointVolumeIds().map(String);
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
    // Keep endpoint catalog ids — insert changes interiors only; clearing ids
    // dropped PC1 tick identity and broke Decode/Render debt accounting.
    return { ok: true, action: "insert", index: insertAt, pathN: poly.length };
  };

  DirectTracePath.prototype._xyForPlotRow = function (row) {
    var hooks = this.hooks || {};
    if (typeof hooks.rowXYForPlotRow === "function") {
      var direct = hooks.rowXYForPlotRow(row);
      if (direct) {
        var pt = cloneXY(direct);
        if (pt) return pt;
      }
    }
    if (typeof hooks.trajXYFromPlotRows === "function") {
      var rebuilt = hooks.trajXYFromPlotRows([row]);
      if (rebuilt && rebuilt[0]) return cloneXY(rebuilt[0]);
    }
    return null;
  };

  /**
   * Map a direct-trace latent polyline to dataset plot rows when only XY is
   * stored (e.g. endpoint line before densify row materialization).
   */
  DirectTracePath.prototype._derivePlotRowsForXY = function (xy) {
    xy = cloneXYList(xy);
    if (!xy || xy.length < 2) return null;
    var hooks = this.hooks || {};
    if (typeof hooks.plotRowsForPathXY === "function") {
      var batch = hooks.plotRowsForPathXY(xy);
      if (batch && batch.length === xy.length) {
        return cloneRowList(batch).map(Number).filter(function (r) {
          return Number.isFinite(r) && r >= 0;
        }).length === xy.length ? cloneRowList(batch) : null;
      }
    }
    if (typeof hooks.nearestPlotRowForXY === "function") {
      var rows = [];
      for (var i = 0; i < xy.length; i++) {
        var pt = xy[i];
        if (!pt || pt.length < 2) return null;
        var row = Number(hooks.nearestPlotRowForXY(Number(pt[0]), Number(pt[1])));
        if (!Number.isFinite(row) || row < 0) return null;
        rows.push(row);
      }
      return rows;
    }
    return null;
  };

  /**
   * Snapshot the live direct-trace polyline (densified samples or anchors).
   */
  DirectTracePath.prototype._livePathRowsAndXY = function () {
    var rows = null;
    var xy = null;
    var usedEditable = false;

    var samples = this.getSamplePlotRows();
    if (samples && samples.length >= 2) {
      rows = samples.slice();
      xy = this.getEditableXY();
      if (!xy || xy.length !== rows.length) xy = this.store.anchorTrajXY;
      if ((!xy || xy.length !== rows.length)
          && typeof this.hooks.trajXYFromPlotRows === "function") {
        xy = this.hooks.trajXYFromPlotRows(rows);
      }
      if (xy && xy.length === rows.length) {
        return { rows: rows, xy: cloneXYList(xy), usedEditable: !!this.getEditableXY() };
      }
    }

    var anchors = this.getAnchorRows();
    if (anchors && anchors.length >= 2) {
      rows = anchors.slice();
      xy = this.store.anchorTrajXY;
      if ((!xy || xy.length !== rows.length)
          && typeof this.hooks.trajXYFromPlotRows === "function") {
        xy = this.hooks.trajXYFromPlotRows(rows);
      }
      if (xy && xy.length === rows.length) {
        return { rows: rows, xy: cloneXYList(xy), usedEditable: false };
      }
    }

    var editable = this.getEditableXY();
    if (editable && editable.length >= 2) {
      rows = (samples && samples.length === editable.length)
        ? samples.slice()
        : ((anchors && anchors.length === editable.length) ? anchors.slice() : null);
      if (!rows || rows.length !== editable.length) {
        rows = this._derivePlotRowsForXY(editable);
      }
      if (rows && rows.length === editable.length) {
        return {
          rows: rows,
          xy: cloneXYList(editable),
          usedEditable: true
        };
      }
      return {
        rows: null,
        xy: cloneXYList(editable),
        usedEditable: true
      };
    }

    return null;
  };

  /**
   * Append custom particle rows to the end of the current direct-trace path
   * without collapsing densified geometry back to catalog anchors.
   */
  DirectTracePath.prototype.appendParticleRows = function (plotRows) {
    plotRows = cloneRowList(plotRows).map(Number).filter(function (r) {
      return Number.isFinite(r) && r >= 0;
    });
    if (!plotRows.length) return { ok: false, reason: "empty" };

    var live = this._livePathRowsAndXY();
    if (!live || !live.xy || live.xy.length < 2) {
      return { ok: false, reason: "no-path", useRebuild: true };
    }

    var rows = live.rows ? live.rows.slice() : null;
    var xy = live.xy.slice();
    if (rows && rows.length !== xy.length) {
      if (typeof this.hooks.trajXYFromPlotRows === "function") {
        var aligned = this.hooks.trajXYFromPlotRows(rows);
        if (aligned && aligned.length === rows.length) {
          xy = cloneXYList(aligned);
        } else {
          return { ok: false, reason: "row-xy-mismatch", useRebuild: true };
        }
      } else {
        return { ok: false, reason: "row-xy-mismatch", useRebuild: true };
      }
    } else if (!rows) {
      rows = this._derivePlotRowsForXY(xy);
      if (!rows || rows.length !== xy.length) {
        return { ok: false, reason: "no-row-map", useRebuild: true };
      }
    }

    var seen = {};
    for (var i = 0; i < rows.length; i++) seen[rows[i]] = true;

    var added = 0;
    for (var j = 0; j < plotRows.length; j++) {
      var row = plotRows[j];
      if (seen[row]) continue;
      var pt = this._xyForPlotRow(row);
      if (!pt) continue;
      rows.push(row);
      xy.push(pt);
      seen[row] = true;
      added++;
    }

    if (!added) return { ok: false, reason: "no-new-rows" };
    if (rows.length < 2) return { ok: false, reason: "too-few" };

    this.setSamplePlotRows(rows);
    this.setAnchorRows(rows);
    this.setInterpolatedCount(rows.length);
    this.store.lastLatentTrajectoryPointCount = rows.length;
    this.setEndpoints(xy[0], xy[xy.length - 1]);
    if (live.usedEditable) {
      this.setEditableXY(xy);
      this.setAnchorPathXY(null);
    } else {
      this.setAnchorPathXY(xy);
      this.setEditableXY(null);
    }

    return { ok: true, added: added, pathN: rows.length };
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

  WaypointPath.prototype.interpolationArmed = function () {
    return !!(this.store.manualVolumeSnapActive
      || this.store.manualParticleSnapPathActive
      || this.store.graphTraversalActive);
  };

  WaypointPath.prototype.armInterpolation = function (kind) {
    // "volume" / direct-line: densify along the current waypoint polyline.
    // Do not arm particle-midpoint snap — that relocates catalog/PC paths onto
    // the particle cloud and destroys the displayed trajectory.
    var volume = kind === "volume";
    this.store.manualVolumeSnapActive = volume;
    this.store.manualParticleSnapPathActive = false;
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

  /**
   * Subset key for catalog volumes (PC1 / PC2 / kmeans). Returns null for
   * Other / unknown ids so Original-order sorting can ignore them.
   */
  function catalogVolumeSubsetKey(volId, catalogById) {
    volId = String(volId || "");
    var entry = catalogById && catalogById[volId];
    if (entry) {
      if (entry.kind === "kmeans") return "kmeans";
      if (entry.kind === "pc") return "pc" + String(entry.pc);
      return null;
    }
    if (/^kmeans:/.test(volId)) return "kmeans";
    var pc = /^pc(\d+):/.exec(volId);
    return pc ? ("pc" + pc[1]) : null;
  }

  /**
   * Native catalog index used by Original order: sample_index for PC volumes,
   * cluster_label for kmeans.
   */
  function catalogVolumeOriginalIndex(volId, catalogById) {
    volId = String(volId || "");
    var entry = catalogById && catalogById[volId];
    if (entry) {
      if (entry.kind === "kmeans") return Number(entry.cluster_label);
      if (entry.kind === "pc") return Number(entry.sample_index);
      return Infinity;
    }
    var km = /^kmeans:(\d+)$/.exec(volId);
    if (km) return Number(km[1]);
    var pc = /^pc(\d+):(\d+)$/.exec(volId);
    if (pc) return Number(pc[2]);
    return Infinity;
  }

  function catalogVolumeSubsetRank(volId, catalogById) {
    var key = catalogVolumeSubsetKey(volId, catalogById);
    if (key === "kmeans") return 0;
    if (key === "pc1") return 1;
    if (key === "pc2") return 2;
    return 9;
  }

  /**
   * Sort selected catalog volume ids by subset then native index. Non-catalog
   * ids keep relative order at the end (Original order ignores them as keys).
   */
  function sortCatalogVolumeIdsByOriginalIndex(ids, catalogById) {
    var catalog = [];
    var other = [];
    (ids || []).forEach(function (id) {
      id = String(id);
      if (catalogVolumeSubsetKey(id, catalogById)) catalog.push(id);
      else other.push(id);
    });
    catalog.sort(function (a, b) {
      var ra = catalogVolumeSubsetRank(a, catalogById);
      var rb = catalogVolumeSubsetRank(b, catalogById);
      if (ra !== rb) return ra - rb;
      return catalogVolumeOriginalIndex(a, catalogById)
        - catalogVolumeOriginalIndex(b, catalogById);
    });
    return catalog.concat(other);
  }

  /**
   * Original-order path rewrite: within each PC1 / PC2 / kmeans subset, replace
   * catalog slots so those volumes appear in native index order. Other / custom
   * slots keep their path positions and are ignored when ordering subsets.
   *
   * ``slots`` entries are ``{ kind: "catalog", volId }`` or ``{ kind: "custom", row }``.
   * Returns ``{ slots, volIds, customRows }`` with catalog slots rewritten.
   */
  function reorderSlotsByOriginalCatalogIndex(slots, catalogById) {
    slots = (slots || []).map(function (slot) {
      return slot && typeof slot === "object" ? Object.assign({}, slot) : slot;
    });
    var bySubset = {};
    var i;
    for (i = 0; i < slots.length; i++) {
      var slot = slots[i];
      if (!slot || slot.kind !== "catalog") continue;
      var key = catalogVolumeSubsetKey(slot.volId, catalogById);
      if (!key) continue;
      if (!bySubset[key]) bySubset[key] = [];
      bySubset[key].push(String(slot.volId));
    }
    Object.keys(bySubset).forEach(function (key) {
      bySubset[key].sort(function (a, b) {
        return catalogVolumeOriginalIndex(a, catalogById)
          - catalogVolumeOriginalIndex(b, catalogById);
      });
    });
    var cursors = {};
    var volIds = [];
    var customRows = [];
    for (i = 0; i < slots.length; i++) {
      slot = slots[i];
      if (!slot) continue;
      if (slot.kind === "custom") {
        customRows.push(Number(slot.row));
        continue;
      }
      if (slot.kind !== "catalog") continue;
      var subset = catalogVolumeSubsetKey(slot.volId, catalogById);
      if (!subset || !bySubset[subset]) {
        volIds.push(String(slot.volId));
        continue;
      }
      if (cursors[subset] == null) cursors[subset] = 0;
      var nextId = bySubset[subset][cursors[subset]++];
      slot.volId = nextId;
      volIds.push(nextId);
    }
    return { slots: slots, volIds: volIds, customRows: customRows };
  }

  global.CryoTrajectoryPath = {
    TrajectoryPath: TrajectoryPath,
    DirectTracePath: DirectTracePath,
    WaypointPath: WaypointPath,
    linearXYFromEndpoints: linearXYFromEndpoints,
    createForStore: createForStore,
    ensureCurrentPath: ensureCurrentPath,
    catalogVolumeSubsetKey: catalogVolumeSubsetKey,
    catalogVolumeOriginalIndex: catalogVolumeOriginalIndex,
    sortCatalogVolumeIdsByOriginalIndex: sortCatalogVolumeIdsByOriginalIndex,
    reorderSlotsByOriginalCatalogIndex: reorderSlotsByOriginalCatalogIndex
  };
})(typeof window !== "undefined" ? window : this);
