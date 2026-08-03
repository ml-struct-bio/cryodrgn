/**
 * Path mutation command objects for the trajectory creator.
 *
 * Each command mutates TrajectoryPath geometry and optionally aligns
 * TrajectoryVolumeState. TrajectorySession.dispatch runs them uniformly so
 * Reverse / visit-order / add-points share one post-mutation chrome path.
 */
(function (global) {
  "use strict";

  function cloneRows(rows) {
    return (rows || []).map(Number).filter(function (r) { return Number.isFinite(r); });
  }

  function PathMutation(name) {
    this.name = String(name || "mutation");
  }

  PathMutation.prototype.apply = function (/* session */) {
    return { ok: false, reason: "unimplemented" };
  };

  function ReverseMutation() {
    PathMutation.call(this, "reverse");
  }
  ReverseMutation.prototype = Object.create(PathMutation.prototype);
  ReverseMutation.prototype.constructor = ReverseMutation;
  ReverseMutation.prototype.apply = function (session) {
    var path = session.path();
    var volumes = session.volumes();
    if (!path) return { ok: false, reason: "no-path" };
    var n = path.sampleCount();
    if (n < 2 && !(path.getSelectedVolumeIds && path.getSelectedVolumeIds().length >= 2)) {
      var rows = path.getAnchorRows();
      if (!(rows && rows.length >= 2)) return { ok: false, reason: "too-short" };
    }
    // Pure permutation: reverse path geometry and reverse the current volume
    // view. Do NOT re-derive slot ids from selection after reverse — that can
    // assign different catalog ids to the same plot rows and disturb Decode /
    // Render button counts even though media is unchanged in the id map.
    path.reverse();
    if (volumes && volumes.slotCount() >= 2) {
      var usedCompact = false;
      if (typeof session.hooks.shouldReverseCompactCatalog === "function"
          && session.hooks.shouldReverseCompactCatalog(session)) {
        if (typeof session.hooks.reverseCompactCatalogVolumes === "function") {
          usedCompact = !!session.hooks.reverseCompactCatalogVolumes(session);
        }
      }
      if (!usedCompact) {
        volumes.reverse();
      }
    }
    return { ok: true, mutation: this.name };
  };

  function VisitOrderMutation(order, orderedRows) {
    PathMutation.call(this, "visit-order");
    this.order = String(order || "preserve");
    this.orderedRows = cloneRows(orderedRows);
  }
  VisitOrderMutation.prototype = Object.create(PathMutation.prototype);
  VisitOrderMutation.prototype.constructor = VisitOrderMutation;
  VisitOrderMutation.prototype.apply = function (session) {
    var path = session.path();
    var volumes = session.volumes();
    if (!path) return { ok: false, reason: "no-path" };
    var rows = this.orderedRows;
    if (rows.length < 2) {
      return { ok: true, mutation: this.name, order: this.order, noop: true };
    }
    if (typeof session.hooks.reorderSelectionFromRows === "function") {
      session.hooks.reorderSelectionFromRows(rows);
    }
    // Do not rebuildFromSelection here: that helper historically emitted
    // catalog volumes then Other indices, undoing an interleaved Exact /
    // Inexact tour over PC1 + random waypoints.
    if (typeof path.setAnchorsFromOrderedRows === "function") {
      path.setAnchorsFromOrderedRows(rows);
    } else {
      path.setAnchorRows(rows);
      if (path.kind === "waypoint") {
        path.setSamplePlotRows(null);
        path.setEditableXY(null);
        path.setInterpolatedCount(0);
      }
    }
    if (volumes && typeof session.syncVolumesToPath === "function") {
      session.syncVolumesToPath();
    } else {
      var nextIds = typeof session.hooks.slotIdsForPath === "function"
        ? session.hooks.slotIdsForPath(path)
        : null;
      if (volumes && nextIds && nextIds.length >= 2) {
        volumes.alignToIds(nextIds);
      }
    }
    if (typeof session.hooks.setVisitOrderMode === "function") {
      session.hooks.setVisitOrderMode(this.order);
    }
    return { ok: true, mutation: this.name, order: this.order };
  };

  function AppendPlotRowsMutation(plotRows) {
    PathMutation.call(this, "append-plot-rows");
    this.plotRows = cloneRows(plotRows);
  }
  AppendPlotRowsMutation.prototype = Object.create(PathMutation.prototype);
  AppendPlotRowsMutation.prototype.constructor = AppendPlotRowsMutation;
  AppendPlotRowsMutation.prototype.apply = function (session) {
    var path = session.path();
    var volumes = session.volumes();
    if (!path) return { ok: false, reason: "no-path" };
    if (!this.plotRows.length) return { ok: false, reason: "empty" };
    var hasAppendHook = typeof session.hooks.appendPlotRows === "function";
    if (hasAppendHook) {
      var ok = session.hooks.appendPlotRows(this.plotRows);
      if (!ok) return { ok: false, reason: "append-failed" };
    } else if (path.kind === "waypoint" && typeof path.addParticleRow === "function") {
      for (var i = 0; i < this.plotRows.length; i++) {
        path.addParticleRow(this.plotRows[i], { activate: true });
      }
    }
    if (path.kind === "direct" && hasAppendHook
        && typeof path.appendParticleRows === "function") {
      var directAppend = path.appendParticleRows(this.plotRows);
      if (directAppend && directAppend.ok) {
        /* geometry extended in place */
      } else if (directAppend && directAppend.useRebuild) {
        if (typeof path.rebuildFromSelection === "function") {
          path.rebuildFromSelection({ force: true });
        }
      } else if (directAppend && !directAppend.useRebuild) {
        return { ok: false, reason: directAppend.reason || "append-failed" };
      }
    } else if (typeof path.rebuildFromSelection === "function"
        && path.kind === "waypoint") {
      path.rebuildFromSelection({ force: true });
    }
    if (volumes && typeof session.syncVolumesToPath === "function") {
      // Waypoint / direct append grows by durable slot ids — never rematch by
      // latent XY here (that misplaces catalog ChimeraX frames before expand).
      session.syncVolumesToPath({ remapByPathXy: false, rematchPath: false });
    } else {
      var nextIds = typeof session.hooks.slotIdsForPath === "function"
        ? session.hooks.slotIdsForPath(path)
        : null;
      if (volumes && nextIds && nextIds.length >= 2) {
        volumes.alignToIds(nextIds);
      }
    }
    if (typeof session.hooks.resetVisitOrderToOriginal === "function") {
      session.hooks.resetVisitOrderToOriginal();
    }
    return { ok: true, mutation: this.name, added: this.plotRows.length };
  };

  function RebuildFromSelectionMutation(opts) {
    PathMutation.call(this, "rebuild-from-selection");
    this.opts = opts || {};
  }
  RebuildFromSelectionMutation.prototype = Object.create(PathMutation.prototype);
  RebuildFromSelectionMutation.prototype.constructor = RebuildFromSelectionMutation;
  RebuildFromSelectionMutation.prototype.apply = function (session) {
    var path = session.path();
    var volumes = session.volumes();
    if (!path || typeof path.rebuildFromSelection !== "function") {
      return { ok: false, reason: "no-path" };
    }
    var result = path.rebuildFromSelection(this.opts);
    if (volumes && typeof session.syncVolumesToPath === "function") {
      session.syncVolumesToPath();
    } else {
      var nextIds = typeof session.hooks.slotIdsForPath === "function"
        ? session.hooks.slotIdsForPath(path)
        : null;
      if (volumes && nextIds && nextIds.length >= 2) {
        volumes.alignToIds(nextIds);
      } else if (volumes && (!nextIds || nextIds.length < 2)) {
        volumes.clear();
      }
    }
    return Object.assign({ ok: true, mutation: this.name }, result || {});
  };

  function ClearPathMutation() {
    PathMutation.call(this, "clear");
  }
  ClearPathMutation.prototype = Object.create(PathMutation.prototype);
  ClearPathMutation.prototype.constructor = ClearPathMutation;
  ClearPathMutation.prototype.apply = function (session) {
    var path = session.path();
    var volumes = session.volumes();
    if (path && typeof path.clearGeometry === "function") path.clearGeometry();
    if (volumes) volumes.clear();
    return { ok: true, mutation: this.name };
  };

  function AlignVolumesMutation(ids) {
    PathMutation.call(this, "align-volumes");
    this.ids = (ids || []).map(function (id) {
      return id == null || id === "" ? null : String(id);
    });
  }
  AlignVolumesMutation.prototype = Object.create(PathMutation.prototype);
  AlignVolumesMutation.prototype.constructor = AlignVolumesMutation;
  AlignVolumesMutation.prototype.apply = function (session) {
    var volumes = session.volumes();
    if (!volumes) return { ok: false, reason: "no-volumes" };
    var ids = this.ids.length
      ? this.ids
      : (typeof session.hooks.slotIdsForPath === "function"
        ? session.hooks.slotIdsForPath(session.path())
        : []);
    if (!ids || ids.length < 2) {
      volumes.clear();
      return { ok: true, mutation: this.name, cleared: true };
    }
    volumes.alignToIds(ids);
    return { ok: true, mutation: this.name, count: ids.length };
  };

  global.CryoTrajectoryPathMutations = {
    PathMutation: PathMutation,
    ReverseMutation: ReverseMutation,
    VisitOrderMutation: VisitOrderMutation,
    AppendPlotRowsMutation: AppendPlotRowsMutation,
    RebuildFromSelectionMutation: RebuildFromSelectionMutation,
    ClearPathMutation: ClearPathMutation,
    AlignVolumesMutation: AlignVolumesMutation,
    reverse: function () { return new ReverseMutation(); },
    visitOrder: function (order, rows) { return new VisitOrderMutation(order, rows); },
    appendPlotRows: function (rows) { return new AppendPlotRowsMutation(rows); },
    rebuildFromSelection: function (opts) { return new RebuildFromSelectionMutation(opts); },
    clear: function () { return new ClearPathMutation(); },
    alignVolumes: function (ids) { return new AlignVolumesMutation(ids); }
  };
})(typeof window !== "undefined" ? window : this);
