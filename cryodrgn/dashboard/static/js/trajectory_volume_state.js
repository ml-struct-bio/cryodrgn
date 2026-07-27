/**
 * Per-slot decode / render status for a trajectory path.
 *
 * Owns catalog / waypoint ids, VTK blobs, ChimeraX images, and path-match
 * metadata used by the volume slider. Slot identity is primarily a durable
 * catalog / custom id. For null-id interiors (direct-trace densify), latent
 * ``traj_xy`` / ``path_t`` on each decoded volume is the stable link to the
 * scatter polyline and slider tick under Trajectory-controls length changes —
 * never bare slot index alone.
 *
 * Invariants (ChimeraX / Decode/Render button):
 *   - ``renderDebtCount()`` === ``inactiveTickCount("chimerax")`` (alias).
 *   - ``decodeDebtIndices()`` ⊆ ``inactiveTickIndices("chimerax")``, so
 *     ``decodeDebtCount()`` ≤ ``inactiveTickCount("chimerax")`` by set inclusion
 *     (decode debt is never larger than inactive slider ticks).
 */
(function (global) {
  "use strict";

  function asId(v) {
    if (v == null || v === "") return null;
    return String(v);
  }

  function cloneVol(vol) {
    if (!vol || typeof vol !== "object") return vol || null;
    return Object.assign({}, vol);
  }

  function cloneImage(img) {
    return img || null;
  }

  function finitePlotRow(v) {
    if (v == null || v === "") return null;
    var n = Math.floor(Number(v));
    return Number.isFinite(n) && n >= 0 ? n : null;
  }

  function cloneMeta(meta, id, vol) {
    meta = meta || {};
    vol = vol || null;
    var membership = asId(
      meta.membership != null ? meta.membership
        : (meta.volId != null ? meta.volId
          : (meta.id != null ? meta.id
            : (id != null ? id
              : (vol && (vol.catalog_id || vol.vol_id || vol.id)))))
    );
    var plotRow = finitePlotRow(
      meta.plot_row != null ? meta.plot_row
        : (meta.plotRow != null ? meta.plotRow
          : (meta.particle_index != null ? meta.particle_index
            : (meta.particleIndex != null ? meta.particleIndex
              : (vol && (vol.plot_row != null ? vol.plot_row : vol.particle_index)))))
    );
    return {
      plot_row: plotRow,
      membership: membership
    };
  }

  function cloneMatch(match) {
    match = match || {};
    return {
      matched: match.matched !== false,
      labelId: asId(match.labelId != null ? match.labelId : match.membership),
      plot_row: finitePlotRow(match.plot_row != null ? match.plot_row : match.plotRow)
    };
  }

  function normalizePathSamples(samples, ids) {
    if (!samples || !samples.length) return [];
    var out = [];
    for (var i = 0; i < samples.length; i++) {
      var src = samples[i];
      var sample;
      if (Array.isArray(src)) {
        sample = { xy: cloneXy(src) };
      } else if (src && typeof src === "object") {
        sample = {
          xy: cloneXy(src.xy || src.traj_xy || src.path_xy),
          plot_row: finitePlotRow(src.plot_row != null ? src.plot_row : src.plotRow),
          membership: asId(src.membership != null ? src.membership
            : (src.volId != null ? src.volId : (src.id != null ? src.id : null)))
        };
      } else {
        sample = { xy: null };
      }
      if (!sample.membership && ids && i < ids.length) sample.membership = asId(ids[i]);
      out.push(sample);
    }
    return out;
  }

  function finiteInt(v, fallback) {
    var n = Math.floor(Number(v));
    return Number.isFinite(n) ? n : (fallback || 0);
  }

  function anchorSlotIndices(nAnchors, nInterp, pathN) {
    nAnchors = Math.max(0, finiteInt(nAnchors, 0));
    nInterp = Math.max(0, finiteInt(nInterp, 0));
    pathN = Math.max(0, finiteInt(pathN, 0));
    if (nAnchors < 2) return [];
    if (nInterp <= 0) {
      var only = [];
      for (var ai = 0; ai < nAnchors; ai++) {
        if (!pathN || ai < pathN) only.push(ai);
      }
      return only;
    }
    var perSeg = nInterp + 1;
    var slots = [];
    for (var i = 0; i < nAnchors; i++) {
      var slot = i * perSeg;
      if (!pathN || slot < pathN) slots.push(slot);
    }
    return slots;
  }

  function compactCatalogPackedPrefix(images, pathN, nAnchors, nInterp) {
    images = images || [];
    pathN = Math.max(0, finiteInt(pathN, 0));
    nAnchors = Math.max(0, finiteInt(nAnchors, 0));
    if (nAnchors < 2 || pathN <= nAnchors || images.length <= nAnchors) return false;
    var slots = anchorSlotIndices(nAnchors, nInterp, pathN);
    if (slots.length < nAnchors) return false;
    var prefixReady = 0;
    for (var i = 0; i < nAnchors && i < images.length; i++) {
      if (images[i]) prefixReady++;
    }
    var upperTotal = 0;
    var upperReady = 0;
    for (var si = 0; si < slots.length; si++) {
      if (slots[si] < nAnchors) continue;
      upperTotal++;
      if (slots[si] < images.length && images[slots[si]]) upperReady++;
    }
    return prefixReady >= 2 && upperTotal > 0 && upperReady === 0;
  }

  function pathSampleCount(path) {
    if (path && typeof path.sampleCount === "function") {
      var n = finiteInt(path.sampleCount(), 0);
      if (n >= 2) return n;
    }
    if (path && typeof path.getPathXY === "function") {
      var xy = path.getPathXY();
      if (xy && xy.length >= 2) return xy.length;
    }
    if (path && typeof path.getSamplePlotRows === "function") {
      var rows = path.getSamplePlotRows();
      if (rows && rows.length >= 2) return rows.length;
    }
    if (path && typeof path.getAnchorRows === "function") {
      var anchors = path.getAnchorRows();
      if (anchors && anchors.length >= 2) return anchors.length;
    }
    return 0;
  }

  function permuteArray(arr, permutation) {
    permutation = permutation || [];
    var n = permutation.length;
    var padded = (arr || []).slice();
    while (padded.length < n) padded.push(null);
    if (padded.length > n) padded.length = n;
    var out = new Array(n);
    for (var i = 0; i < n; i++) {
      var src = permutation[i];
      out[i] = (src >= 0 && src < padded.length) ? padded[src] : null;
    }
    return out;
  }

  function cloneXy(xy) {
    if (!xy || xy.length < 2) return null;
    var x = Number(xy[0]);
    var y = Number(xy[1]);
    if (!Number.isFinite(x) || !Number.isFinite(y)) return null;
    return [x, y];
  }

  function normalizePathXY(pathXY) {
    if (!pathXY || !pathXY.length) return [];
    var out = [];
    for (var i = 0; i < pathXY.length; i++) {
      var pt = cloneXy(pathXY[i]);
      if (pt) out.push(pt);
    }
    return out;
  }

  function pathParamT(index, n) {
    n = Math.max(0, finiteInt(n, 0));
    index = Math.max(0, finiteInt(index, 0));
    if (n <= 1) return 0;
    return index / (n - 1);
  }

  /**
   * Squared-distance helper (tests / diagnostics). Matching uses per-axis
   * ``latentXyMatchAtol`` to mirror server ``np.allclose`` rounding tolerance.
   */
  function xyDist2(a, b) {
    a = cloneXy(a);
    b = cloneXy(b);
    if (!a || !b) return Infinity;
    var dx = a[0] - b[0];
    var dy = a[1] - b[1];
    return dx * dx + dy * dy;
  }

  /**
   * Half-ulp absolute tolerance matching server ``_round_direct_mode_traj_xy``
   * (3 decimals below |v|=100, else 2).
   */
  function latentXyMatchAtol(a, b) {
    a = cloneXy(a);
    b = cloneXy(b);
    var maxAbs = 0;
    if (a) maxAbs = Math.max(maxAbs, Math.abs(a[0]), Math.abs(a[1]));
    if (b) maxAbs = Math.max(maxAbs, Math.abs(b[0]), Math.abs(b[1]));
    return maxAbs >= 100 ? 5e-3 : 5e-4;
  }

  function latentXyMatchEps2(a, b) {
    var atol = latentXyMatchAtol(a, b);
    return atol * atol;
  }

  function xyMatchesLatent(a, b) {
    a = cloneXy(a);
    b = cloneXy(b);
    if (!a || !b) return false;
    var atol = latentXyMatchAtol(a, b);
    return Math.abs(a[0] - b[0]) <= atol && Math.abs(a[1] - b[1]) <= atol;
  }

  /**
   * Attach latent path geometry onto a decoded volume object so identity survives
   * Trajectory-controls densify / shrink (index alone is not stable).
   */
  function stampVolGeometry(vol, xy, t) {
    if (!vol || typeof vol !== "object") return vol || null;
    var out = cloneVol(vol);
    var pt = cloneXy(xy);
    if (pt) out.traj_xy = pt;
    if (Number.isFinite(Number(t))) out.path_t = Number(t);
    return out;
  }

  function TrajectoryVolumeState(opts) {
    opts = opts || {};
    this._ids = [];
    this._volumes = [];
    this._images = [];
    this._xy = [];
    this._meta = [];
    this._match = [];
    this._stale = {};
    this._lastRematchMap = {};
    // Durable id → media map. Path order is always a view via alignToIds;
    // reverse / visit-order must never drop readiness that still belongs to an id.
    this._byId = Object.create(null);
    this._generated = false;
    this._cacheId = "";
    this._subsetSlotIndices = null;
    this._layout = {
      kind: "empty",
      pathN: 0,
      nAnchors: 0,
      nInterp: 0,
      anchorSlotIndices: []
    };
    this._listeners = [];
    if (opts.ids) this.replaceSlots(opts.ids, opts.volumes, opts.images, opts.xy);
  }

  TrajectoryVolumeState.prototype.onChange = function (fn) {
    if (typeof fn === "function") this._listeners.push(fn);
    return this;
  };

  TrajectoryVolumeState.prototype._emit = function (reason) {
    for (var i = 0; i < this._listeners.length; i++) {
      try { this._listeners[i](this, reason || "change"); } catch (err) { /* ignore */ }
    }
  };

  TrajectoryVolumeState.prototype._rememberId = function (id, vol, img, stale) {
    id = asId(id);
    if (!id) return;
    var prev = this._byId[id] || { volume: null, image: null, stale: false };
    var nextVol = prev.volume;
    var nextImg = prev.image;
    if (vol !== undefined) {
      nextVol = vol ? cloneVol(vol) : null;
    }
    if (img !== undefined) {
      nextImg = img || null;
    }
    // Prefer keeping a stronger prior payload when the incoming slot is empty.
    if (!nextVol && prev.volume && (prev.volume.volume_b64 || prev.volume.decoded === true)
        && vol !== null) {
      nextVol = prev.volume;
    }
    if (!nextImg && prev.image && img !== null) {
      nextImg = prev.image;
    }
    this._byId[id] = {
      volume: nextVol || null,
      image: nextImg || null,
      stale: stale == null ? !!prev.stale : !!stale
    };
  };

  TrajectoryVolumeState.prototype._ensureSlotArrays = function (n) {
    n = Math.max(0, finiteInt(n, 0));
    while (this._ids.length < n) this._ids.push(null);
    while (this._volumes.length < n) this._volumes.push(null);
    while (this._images.length < n) this._images.push(null);
    while (this._xy.length < n) this._xy.push(null);
    while (this._meta.length < n) this._meta.push(cloneMeta(null, this._ids[this._meta.length], null));
    while (this._match.length < n) this._match.push({ matched: true, labelId: null, plot_row: null });
    if (this._ids.length > n) this._ids.length = n;
    if (this._volumes.length > n) this._volumes.length = n;
    if (this._images.length > n) this._images.length = n;
    if (this._xy.length > n) this._xy.length = n;
    if (this._meta.length > n) this._meta.length = n;
    if (this._match.length > n) this._match.length = n;
  };

  TrajectoryVolumeState.prototype._pathIdentityForSlot = function (i) {
    var meta = cloneMeta(this._meta[i], this._ids[i], this._volumes[i]);
    var labelId = meta.membership || this._ids[i] || null;
    return {
      matched: !(this._match[i] && this._match[i].matched === false),
      labelId: labelId,
      plot_row: meta.plot_row
    };
  };

  TrajectoryVolumeState.prototype._rememberSlots = function () {
    for (var i = 0; i < this._ids.length; i++) {
      this._rememberId(
        this._ids[i],
        this._volumes[i] === undefined ? undefined : (this._volumes[i] || null),
        this._images[i] === undefined ? undefined : (this._images[i] || null),
        !!this._stale[i]
      );
    }
  };

  TrajectoryVolumeState.prototype._materializeFromById = function (targetIds) {
    targetIds = (targetIds || []).map(asId);
    var vols = [];
    var imgs = [];
    var stale = {};
    for (var i = 0; i < targetIds.length; i++) {
      var id = targetIds[i];
      var hit = id && this._byId[id] ? this._byId[id] : null;
      vols.push(hit && hit.volume ? cloneVol(hit.volume) : null);
      imgs.push(hit && hit.image ? hit.image : null);
      if (hit && hit.stale) stale[i] = true;
    }
    return { volumes: vols, images: imgs, stale: stale };
  };

  TrajectoryVolumeState.prototype.slotCount = function () {
    return this._ids.length;
  };

  TrajectoryVolumeState.prototype.slotIdAt = function (i) {
    i = Math.floor(Number(i));
    if (!Number.isFinite(i) || i < 0 || i >= this._ids.length) return null;
    return this._ids[i];
  };

  TrajectoryVolumeState.prototype.slotXyAt = function (i) {
    i = Math.floor(Number(i));
    if (!Number.isFinite(i) || i < 0 || i >= this._xy.length) return null;
    return cloneXy(this._xy[i]);
  };

  TrajectoryVolumeState.prototype.pathXy = function () {
    return this._xy.map(function (pt) { return cloneXy(pt); });
  };

  TrajectoryVolumeState.prototype.setSlotXy = function (i, xy) {
    i = Math.floor(Number(i));
    if (!Number.isFinite(i) || i < 0) return this;
    while (this._xy.length <= i) this._xy.push(null);
    while (this._ids.length <= i) this._ids.push(null);
    while (this._volumes.length <= i) this._volumes.push(null);
    while (this._images.length <= i) this._images.push(null);
    this._xy[i] = cloneXy(xy);
    if (this._volumes[i]) {
      this._volumes[i] = stampVolGeometry(
        this._volumes[i],
        this._xy[i],
        pathParamT(i, this._ids.length)
      );
    }
    return this;
  };

  /**
   * Latent XY where slot media was decoded (vol.traj_xy), if any.
   * Distinct from live path ``_xy`` after nearest-snap / drag.
   */
  TrajectoryVolumeState.prototype.decodedXyAt = function (i) {
    i = Math.floor(Number(i));
    if (!Number.isFinite(i) || i < 0 || i >= this._volumes.length) return null;
    var vol = this._volumes[i];
    return vol ? cloneXy(vol.traj_xy) : null;
  };

  /**
   * Indices whose decoded / rendered media is not on the live polyline sample.
   * Tick activity should follow this: active iff scatter point sits on a
   * decoded location.
   */
  TrajectoryVolumeState.prototype.indicesOffDecodedPath = function (pathXY, eps2) {
    pathXY = normalizePathXY(pathXY);
    var useFixedEps = eps2 != null && Number.isFinite(Number(eps2));
    if (useFixedEps) eps2 = Number(eps2);
    var n = Math.max(this._ids.length, pathXY.length);
    var out = [];
    for (var i = 0; i < n; i++) {
      var hasMedia = !!(this._volumes[i] || this._images[i]);
      if (!hasMedia) continue;
      var decoded = this.decodedXyAt(i) || cloneXy(this._xy[i]);
      var live = i < pathXY.length ? pathXY[i] : null;
      if (!decoded || !live) {
        out.push(i);
        continue;
      }
      if (useFixedEps) {
        if (xyDist2(decoded, live) > eps2) out.push(i);
      } else if (!xyMatchesLatent(decoded, live)) {
        out.push(i);
      }
    }
    return out;
  };

  /**
   * Indices with media whose sample moved between two polylines (nearest snap).
   */
  TrajectoryVolumeState.prototype.indicesMovedBetweenPaths = function (preXY, postXY) {
    preXY = normalizePathXY(preXY);
    postXY = normalizePathXY(postXY);
    var n = Math.max(this._ids.length, preXY.length, postXY.length);
    var out = [];
    for (var i = 0; i < n; i++) {
      if (!(this._volumes[i] || this._images[i])) continue;
      var pre = i < preXY.length ? preXY[i] : null;
      var post = i < postXY.length ? postXY[i] : null;
      if (!pre || !post || !xyMatchesLatent(pre, post)) out.push(i);
    }
    return out;
  };

  /**
   * Stamp the live polyline onto slots without moving media.
   * When a slot already has decoded media at a different latent XY, keep
   * ``vol.traj_xy`` as decode identity (nearest snap must not pretend the
   * volume was decoded at the snapped particle). Matching or unset geometry
   * is rebound so densify pre-stamp / post-realign still work.
   */
  TrajectoryVolumeState.prototype.stampPathXy = function (pathXY) {
    pathXY = normalizePathXY(pathXY);
    var n = Math.max(this._ids.length, pathXY.length);
    this._ensureSlotArrays(n);
    for (var i = 0; i < pathXY.length; i++) {
      this._xy[i] = pathXY[i];
      this._meta[i] = cloneMeta(this._meta[i], this._ids[i], this._volumes[i]);
      if (!this._match[i]) this._match[i] = this._pathIdentityForSlot(i);
      if (this._volumes[i]) {
        var existing = cloneXy(this._volumes[i].traj_xy);
        if (!existing || xyMatchesLatent(existing, pathXY[i])) {
          this._volumes[i] = stampVolGeometry(
            this._volumes[i],
            pathXY[i],
            pathParamT(i, Math.max(pathXY.length, 2))
          );
          this._rememberId(this._ids[i], this._volumes[i], this._images[i], !!this._stale[i]);
        }
      }
    }
    this._emit("stamp-xy");
    return this;
  };

  TrajectoryVolumeState.prototype._slotPathT = function (i) {
    i = Math.floor(Number(i));
    var vol = (i >= 0 && i < this._volumes.length) ? this._volumes[i] : null;
    if (vol && Number.isFinite(Number(vol.path_t))) return Number(vol.path_t);
    return pathParamT(i, this._ids.length);
  };

  TrajectoryVolumeState.prototype._harvestGeomRecords = function () {
    var n = this._ids.length;
    var records = [];
    for (var i = 0; i < n; i++) {
      var vol = this._volumes[i] || null;
      var xy = cloneXy(this._xy[i])
        || (vol && cloneXy(vol.traj_xy))
        || null;
      var meta = cloneMeta(this._meta[i], this._ids[i], vol);
      records.push({
        src: i,
        id: this._ids[i],
        vol: vol ? cloneVol(vol) : null,
        img: this._images[i] || null,
        stale: !!this._stale[i],
        xy: xy,
        meta: meta,
        match: cloneMatch(this._match[i] || {
          matched: true,
          labelId: meta.membership || this._ids[i] || null,
          plot_row: meta.plot_row
        }),
        t: this._slotPathT(i),
        used: false
      });
    }
    return records;
  };

  TrajectoryVolumeState.prototype.rematchToPath = function (pathSamples, opts) {
    opts = opts || {};
    var samples = normalizePathSamples(
      pathSamples,
      opts.ids || null
    );
    var pathXY = samples.map(function (sample) { return cloneXy(sample.xy); });
    var n = samples.length;
    if (n < 2) return this;
    if (!this._ids.length) {
      this._ids = new Array(n);
      this._volumes = new Array(n);
      this._images = new Array(n);
      this._xy = pathXY.map(cloneXy);
      this._meta = new Array(n);
      this._match = new Array(n);
      if (opts.ids && opts.ids.length === n) {
        this._ids = opts.ids.map(asId);
      }
      for (var zi = 0; zi < n; zi++) {
        if (!this._ids[zi] && samples[zi].membership) this._ids[zi] = asId(samples[zi].membership);
        this._volumes[zi] = null;
        this._images[zi] = null;
        this._meta[zi] = cloneMeta(samples[zi], this._ids[zi], null);
        this._match[zi] = {
          matched: true,
          labelId: this._meta[zi].membership || this._ids[zi] || null,
          plot_row: this._meta[zi].plot_row
        };
      }
      this._stale = {};
      this._emit("rematch-path");
      return this;
    }

    this._rememberSlots();
    var records = this._harvestGeomRecords();
    var oldN = records.length;
    var targetIds = (opts.ids && opts.ids.length === n)
      ? opts.ids.map(asId)
      : null;
    if (!targetIds) {
      targetIds = new Array(n);
      for (var ti = 0; ti < n; ti++) targetIds[ti] = null;
      // Preserve catalog endpoint ids at the new ends when present.
      if (oldN >= 2 && records[0].id) targetIds[0] = records[0].id;
      if (oldN >= 2 && records[oldN - 1].id) targetIds[n - 1] = records[oldN - 1].id;
      if (opts.compactIds && opts.compactIds.length >= 2) {
        targetIds[0] = asId(opts.compactIds[0]);
        targetIds[n - 1] = asId(opts.compactIds[opts.compactIds.length - 1]);
      }
    }
    for (var si = 0; si < n; si++) {
      if (!targetIds[si] && samples[si].membership) targetIds[si] = asId(samples[si].membership);
    }

    var nextVols = new Array(n);
    var nextImgs = new Array(n);
    var nextMeta = new Array(n);
    var nextMatch = new Array(n);
    var nextStale = {};
    var rematchMap = {};
    for (var j = 0; j < n; j++) {
      nextVols[j] = null;
      nextImgs[j] = null;
      nextMeta[j] = cloneMeta(samples[j], targetIds[j], null);
      nextMatch[j] = {
        matched: true,
        labelId: nextMeta[j].membership || targetIds[j] || null,
        plot_row: nextMeta[j].plot_row
      };
    }

    function claim(rec, dest) {
      if (!rec || rec.used) return;
      rec.used = true;
      if (rec.src != null && Number.isFinite(Number(rec.src))) {
        rematchMap[String(rec.src)] = dest;
      }
      nextVols[dest] = rec.vol
        ? stampVolGeometry(rec.vol, pathXY[dest], pathParamT(dest, n))
        : null;
      nextImgs[dest] = rec.img || null;
      if (rec.stale) nextStale[dest] = true;
      // Never promote a harvested catalog id onto a free sample — that labelled
      // densified interiors as "PC1 vol10" when an endpoint record snapped near
      // the end. Path identity comes only from targetIds / Pass 1.
      nextMeta[dest] = cloneMeta(samples[dest], targetIds[dest], rec.vol);
      if (rec.meta) {
        if (nextMeta[dest].plot_row == null && rec.meta.plot_row != null) {
          nextMeta[dest].plot_row = rec.meta.plot_row;
        }
      }
      if (targetIds[dest] != null) {
        nextMeta[dest].membership = targetIds[dest];
      } else {
        nextMeta[dest].membership = null;
      }
      nextMatch[dest] = {
        matched: true,
        labelId: targetIds[dest] || null,
        plot_row: nextMeta[dest].plot_row
      };
    }

    // Pass 1: durable catalog / waypoint ids.
    for (var di = 0; di < n; di++) {
      var want = targetIds[di];
      if (!want) continue;
      for (var ri = 0; ri < records.length; ri++) {
        if (records[ri].used) continue;
        if (records[ri].id != null && String(records[ri].id) === String(want)) {
          claim(records[ri], di);
          break;
        }
      }
    }

    // Pass 2: particle identity, when present, is stronger than projected XY.
    for (var pr = 0; pr < n; pr++) {
      if (nextVols[pr] || nextImgs[pr]) continue;
      if (samples[pr].plot_row == null) continue;
      for (var rr = 0; rr < records.length; rr++) {
        if (records[rr].used) continue;
        if (!records[rr].vol && !records[rr].img) continue;
        if (!records[rr].meta || records[rr].meta.plot_row == null) continue;
        if (Number(records[rr].meta.plot_row) === Number(samples[pr].plot_row)) {
          claim(records[rr], pr);
          break;
        }
      }
    }

    // Pass 3: exact latent XY (within server rounding atol).
    for (var nj = 0; nj < n; nj++) {
      if (nextVols[nj] || nextImgs[nj]) continue;
      if (!pathXY[nj]) continue;
      var best = -1;
      var bestDist = Infinity;
      for (var rk = 0; rk < records.length; rk++) {
        if (records[rk].used) continue;
        if (!records[rk].vol && !records[rk].img) continue;
        if (!records[rk].xy) continue;
        if (!xyMatchesLatent(records[rk].xy, pathXY[nj])) continue;
        var d = xyDist2(records[rk].xy, pathXY[nj]);
        if (d < bestDist) {
          bestDist = d;
          best = rk;
        }
      }
      if (best < 0) continue;
      claim(records[best], nj);
    }

    // Pass 4: densify / shrink leftovers. Arc-length resampling (4→12) rarely
    // keeps the old sample XY on the new grid, so exact match alone wrongly
    // inactivated already-rendered interiors. Snap each leftover media record
    // to its nearest empty sample when that sample is within half the local
    // neighbour spacing (Voronoi cell) — never by bare path_t / index.
    function localHalfSpacing2(j) {
      var dPrev = j > 0 ? Math.sqrt(xyDist2(pathXY[j], pathXY[j - 1])) : Infinity;
      var dNext = j < n - 1 ? Math.sqrt(xyDist2(pathXY[j], pathXY[j + 1])) : Infinity;
      var half = 0.5 * Math.min(dPrev, dNext);
      if (!Number.isFinite(half) || half <= 0) {
        return latentXyMatchEps2(pathXY[j], pathXY[j]);
      }
      return half * half;
    }
    var snapPairs = [];
    for (var rr = 0; rr < records.length; rr++) {
      if (records[rr].used) continue;
      if (!records[rr].vol && !records[rr].img) continue;
      if (!records[rr].xy) continue;
      var snapBest = -1;
      var snapDist = Infinity;
      for (var sj = 0; sj < n; sj++) {
        if (nextVols[sj] || nextImgs[sj]) continue;
        if (!pathXY[sj]) continue;
        var sd = xyDist2(records[rr].xy, pathXY[sj]);
        if (sd < snapDist) {
          snapDist = sd;
          snapBest = sj;
        }
      }
      if (snapBest < 0) continue;
      if (snapDist > localHalfSpacing2(snapBest)) continue;
      snapPairs.push({ rec: rr, dest: snapBest, dist: snapDist });
    }
    snapPairs.sort(function (a, b) { return a.dist - b.dist; });
    for (var sp = 0; sp < snapPairs.length; sp++) {
      var pair = snapPairs[sp];
      if (records[pair.rec].used) continue;
      if (nextVols[pair.dest] || nextImgs[pair.dest]) continue;
      claim(records[pair.rec], pair.dest);
    }

    this._ids = targetIds;
    this._volumes = nextVols;
    this._images = nextImgs;
    this._meta = nextMeta;
    this._match = nextMatch;
    this._stale = nextStale;
    this._xy = pathXY.map(cloneXy);
    while (this._xy.length < n) this._xy.push(null);
    this._lastRematchMap = rematchMap;
    this._rememberSlots();
    this._layout = Object.assign({}, this._layout, {
      kind: this._layout.kind === "empty" ? "path" : this._layout.kind,
      pathN: n
    });
    this._emit("rematch-path");
    return this;
  };

  TrajectoryVolumeState.prototype.lastRematchMap = function () {
    return Object.assign({}, this._lastRematchMap || {});
  };

  /**
   * Remap slot media onto a new polyline by durable id first, then by latent
   * XY (exact, then nearest sample within half local spacing). Index and bare
   * path_t alone are not stable under Trajectory-controls densify / shrink.
   */
  TrajectoryVolumeState.prototype.realignToPathXy = function (pathXY, opts) {
    opts = opts || {};
    pathXY = normalizePathXY(pathXY);
    var samples = pathXY.map(function (xy, i) {
      return {
        xy: xy,
        membership: opts.ids && i < opts.ids.length ? opts.ids[i] : null
      };
    });
    return this.rematchToPath(samples, opts);
  };

  TrajectoryVolumeState.prototype.ids = function () {
    return this._ids.slice();
  };

  TrajectoryVolumeState.prototype.volumes = function () {
    return this._volumes.map(cloneVol);
  };

  TrajectoryVolumeState.prototype.images = function () {
    return this._images.slice();
  };

  TrajectoryVolumeState.prototype.meta = function () {
    return this._meta.map(function (m, i) {
      return cloneMeta(m, this._ids[i], this._volumes[i]);
    }, this);
  };

  TrajectoryVolumeState.prototype.match = function () {
    return this._match.map(function (m, i) {
      var meta = cloneMeta(this._meta[i], this._ids[i], this._volumes[i]);
      var out = cloneMatch(m || {
        matched: true,
        labelId: meta.membership || this._ids[i] || null,
        plot_row: meta.plot_row
      });
      if (!out.labelId) out.labelId = meta.membership || this._ids[i] || null;
      if (out.plot_row == null) out.plot_row = meta.plot_row;
      return out;
    }, this);
  };

  TrajectoryVolumeState.prototype.matchAt = function (i) {
    i = Math.floor(Number(i));
    if (!Number.isFinite(i) || i < 0 || i >= this._ids.length) {
      return { matched: false, labelId: null, plot_row: null };
    }
    return this.match()[i];
  };

  TrajectoryVolumeState.prototype.tickLabelAt = function (i) {
    var m = this.matchAt(i);
    return m.labelId || null;
  };

  TrajectoryVolumeState.prototype.isDecoded = function (i) {
    i = Math.floor(Number(i));
    if (i < 0 || i >= this._volumes.length) return false;
    if (this._stale[i]) return false;
    var vol = this._volumes[i];
    if (!vol) return false;
    // Full VTK blob, or cache-only Generate marker ({ decoded: true }).
    return !!(vol.volume_b64 || vol.decoded === true);
  };

  TrajectoryVolumeState.prototype.isVtkHydrated = function (i) {
    i = Math.floor(Number(i));
    return i >= 0 && i < this._volumes.length
      && !!(this._volumes[i] && this._volumes[i].volume_b64);
  };

  TrajectoryVolumeState.prototype.isRendered = function (i) {
    i = Math.floor(Number(i));
    if (i < 0 || i >= this._images.length) return false;
    if (this._stale[i]) return false;
    return !!this._images[i];
  };

  TrajectoryVolumeState.prototype.isReady = function (i, backend) {
    backend = String(backend || "chimerax").toLowerCase();
    if (backend === "chimerax") return this.isRendered(i);
    return this.isDecoded(i);
  };

  TrajectoryVolumeState.prototype.tickReadyAt = function (i, backend) {
    i = Math.floor(Number(i));
    if (!Number.isFinite(i) || i < 0 || i >= this._ids.length) return false;
    if (this._stale[i]) return false;
    var match = this.matchAt(i);
    if (match && match.matched === false) return false;
    return this.isReady(i, backend);
  };

  /**
   * Indices whose slider ticks are inactive for ``backend``.
   * For ChimeraX this is the render-debt set (see class invariant).
   */
  TrajectoryVolumeState.prototype.inactiveTickIndices = function (backend) {
    backend = String(backend || "chimerax").toLowerCase();
    var out = [];
    for (var i = 0; i < this._ids.length; i++) {
      if (!this.tickReadyAt(i, backend)) out.push(i);
    }
    return out;
  };

  TrajectoryVolumeState.prototype.inactiveTickCount = function (backend) {
    return this.inactiveTickIndices(backend).length;
  };

  /**
   * Render debt for Decode/Render: inactive ChimeraX ticks on this path.
   */
  TrajectoryVolumeState.prototype.renderDebtCount = function () {
    return this.inactiveTickCount("chimerax");
  };

  TrajectoryVolumeState.prototype.isStale = function (i) {
    return !!this._stale[Math.floor(Number(i))];
  };

  TrajectoryVolumeState.prototype.decodedCount = function () {
    var n = 0;
    for (var i = 0; i < this._volumes.length; i++) if (this.isDecoded(i)) n++;
    return n;
  };

  TrajectoryVolumeState.prototype.renderedCount = function () {
    var n = 0;
    for (var i = 0; i < this._images.length; i++) if (this.isRendered(i)) n++;
    return n;
  };

  /**
   * Every slot lacking a decode blob, regardless of slider activity.
   * Prefer ``decodeDebtIndices`` for Decode/Render and the pipeline.
   */
  TrajectoryVolumeState.prototype.undecodedIndices = function () {
    var out = [];
    for (var i = 0; i < this._ids.length; i++) {
      if (!this.isDecoded(i)) out.push(i);
    }
    return out;
  };

  /**
   * Decode debt: inactive ChimeraX ticks that still need decode.
   * Defined as a subset of ``inactiveTickIndices("chimerax")``, so
   * ``decodeDebtCount() ≤ inactiveTickCount("chimerax")`` by construction.
   */
  TrajectoryVolumeState.prototype.decodeDebtIndices = function () {
    var inactive = this.inactiveTickIndices("chimerax");
    var out = [];
    for (var i = 0; i < inactive.length; i++) {
      var idx = inactive[i];
      if (!this.isDecoded(idx)) out.push(idx);
    }
    return out;
  };

  TrajectoryVolumeState.prototype.decodeDebtCount = function () {
    return this.decodeDebtIndices().length;
  };

  /**
   * Decode/Render debts from this path alone.
   * ``decode ≤ render`` holds because decode indices ⊆ inactive indices.
   */
  TrajectoryVolumeState.prototype.decodeRenderDebts = function () {
    var renderIndices = this.inactiveTickIndices("chimerax");
    var decodeIndices = [];
    for (var i = 0; i < renderIndices.length; i++) {
      var idx = renderIndices[i];
      if (!this.isDecoded(idx)) decodeIndices.push(idx);
    }
    return {
      decode: decodeIndices.length,
      render: renderIndices.length,
      decodeIndices: decodeIndices,
      renderIndices: renderIndices
    };
  };

  TrajectoryVolumeState.prototype.missingDecodeIndices = function () {
    // Pipeline / button: only inactive undecoded slots (class invariant).
    return this.decodeDebtIndices();
  };

  TrajectoryVolumeState.prototype.missingRenderIndices = function () {
    // Align with slider inactivity (tickReadyAt), not raw image presence alone.
    return this.inactiveTickIndices("chimerax");
  };

  TrajectoryVolumeState.prototype.hasGeneratedTrajectory = function () {
    return !!this._generated && (this.decodedCount() >= 2 || !!this._cacheId);
  };

  TrajectoryVolumeState.prototype.setGenerated = function (on, cacheId) {
    this._generated = !!on;
    this._cacheId = on ? String(cacheId || this._cacheId || "") : "";
    this._emit("generated");
  };

  TrajectoryVolumeState.prototype.cacheId = function () {
    return this._cacheId || "";
  };

  TrajectoryVolumeState.prototype.subsetSlotIndices = function () {
    return this._subsetSlotIndices ? this._subsetSlotIndices.slice() : null;
  };

  TrajectoryVolumeState.prototype.setSubsetSlotIndices = function (indices) {
    if (!indices || !indices.length) {
      this._subsetSlotIndices = null;
      this._emit("subset");
      return this;
    }
    this._subsetSlotIndices = indices.map(function (i) {
      return Math.floor(Number(i));
    }).filter(function (i) {
      return Number.isFinite(i) && i >= 0;
    });
    if (!this._subsetSlotIndices.length) this._subsetSlotIndices = null;
    this._emit("subset");
    return this;
  };

  TrajectoryVolumeState.prototype.layout = function () {
    return {
      kind: this._layout.kind,
      pathN: this._layout.pathN,
      nAnchors: this._layout.nAnchors,
      nInterp: this._layout.nInterp,
      anchorSlotIndices: this._layout.anchorSlotIndices.slice()
    };
  };

  TrajectoryVolumeState.prototype.toPayload = function () {
    var n = this._ids.length;
    var vols = this.volumes();
    var imgs = this.images();
    while (vols.length < n) vols.push(null);
    while (imgs.length < n) imgs.push(null);
    return {
      volumes: vols,
      slots: vols.slice(),
      images: imgs,
      ids: this.ids(),
      meta: this.meta(),
      match: this.match(),
      expectedVolumeCount: n,
      expected_volume_count: n,
      volume_cache_id: this._cacheId || undefined,
      slot_indices: this._subsetSlotIndices ? this._subsetSlotIndices.slice() : undefined
    };
  };

  TrajectoryVolumeState.prototype.replaceSlots = function (ids, volumes, images, xy) {
    ids = (ids || []).map(asId);
    volumes = (volumes || []).slice();
    images = (images || []).slice();
    var pathXY = normalizePathXY(xy);
    var n = ids.length;
    while (volumes.length < n) volumes.push(null);
    while (images.length < n) images.push(null);
    if (volumes.length > n) volumes.length = n;
    if (images.length > n) images.length = n;
    this._ids = ids;
    this._volumes = volumes.map(function (vol, idx) {
      return stampVolGeometry(vol, pathXY[idx] || (vol && vol.traj_xy), pathParamT(idx, n));
    });
    this._images = images.map(cloneImage);
    this._xy = [];
    this._meta = [];
    this._match = [];
    for (var xi = 0; xi < n; xi++) {
      this._xy[xi] = pathXY[xi]
        || (this._volumes[xi] && cloneXy(this._volumes[xi].traj_xy))
        || null;
      this._meta[xi] = cloneMeta(null, this._ids[xi], this._volumes[xi]);
      this._match[xi] = {
        matched: true,
        labelId: this._meta[xi].membership || this._ids[xi] || null,
        plot_row: this._meta[xi].plot_row
      };
    }
    this._stale = {};
    this._rememberSlots();
    this._layout = {
      kind: ids.length ? "path" : "empty",
      pathN: ids.length,
      nAnchors: 0,
      nInterp: 0,
      anchorSlotIndices: []
    };
    this._emit("replace");
    return this;
  };

  TrajectoryVolumeState.prototype.setDecoded = function (i, vol) {
    i = Math.floor(Number(i));
    if (!Number.isFinite(i) || i < 0) return this;
    while (this._volumes.length <= i) this._volumes.push(null);
    while (this._ids.length <= i) this._ids.push(null);
    while (this._images.length <= i) this._images.push(null);
    while (this._xy.length <= i) this._xy.push(null);
    while (this._meta.length <= i) this._meta.push(cloneMeta(null, this._ids[this._meta.length], null));
    while (this._match.length <= i) this._match.push({ matched: true, labelId: null, plot_row: null });
    var xy = this._xy[i] || (vol && vol.traj_xy) || null;
    this._volumes[i] = stampVolGeometry(vol, xy, pathParamT(i, Math.max(this._ids.length, i + 1)));
    if (xy) this._xy[i] = cloneXy(xy);
    this._meta[i] = cloneMeta(this._meta[i], this._ids[i], this._volumes[i]);
    this._match[i] = {
      matched: true,
      labelId: this._meta[i].membership || this._ids[i] || null,
      plot_row: this._meta[i].plot_row
    };
    // Never clear stale on decode alone. Freely moved direct-trace slots must
    // stay in ChimeraX render debt until setRendered writes a frame at the new
    // XY. Clearing stale here let seedCompactCatalog re-apply old catalog
    // endpoint PNGs and zero Decode/Render debt while the page display still
    // nulls invalidated ticks (slider falls back to two interior frames).
    this._rememberId(this._ids[i], this._volumes[i], this._images[i], !!this._stale[i]);
    this._emit("decode");
    return this;
  };

  TrajectoryVolumeState.prototype.setRendered = function (i, imageB64) {
    i = Math.floor(Number(i));
    if (!Number.isFinite(i) || i < 0) return this;
    while (this._images.length <= i) this._images.push(null);
    while (this._ids.length <= i) this._ids.push(null);
    while (this._volumes.length <= i) this._volumes.push(null);
    while (this._xy.length <= i) this._xy.push(null);
    while (this._meta.length <= i) this._meta.push(cloneMeta(null, this._ids[this._meta.length], null));
    while (this._match.length <= i) this._match.push({ matched: true, labelId: null, plot_row: null });
    this._images[i] = imageB64 || null;
    this._meta[i] = cloneMeta(this._meta[i], this._ids[i], this._volumes[i]);
    this._match[i] = {
      matched: true,
      labelId: this._meta[i].membership || this._ids[i] || null,
      plot_row: this._meta[i].plot_row
    };
    if (imageB64) delete this._stale[i];
    this._rememberId(this._ids[i], this._volumes[i], this._images[i], !!this._stale[i]);
    this._emit("render");
    return this;
  };

  TrajectoryVolumeState.prototype.markStale = function (i) {
    i = Math.floor(Number(i));
    if (Number.isFinite(i) && i >= 0) {
      this._stale[i] = true;
      this._rememberId(this._ids[i], this._volumes[i], this._images[i], true);
    }
    this._emit("stale");
    return this;
  };

  TrajectoryVolumeState.prototype.clearStaleAt = function (i) {
    i = Math.floor(Number(i));
    if (Number.isFinite(i) && i >= 0 && this._stale[i]) {
      delete this._stale[i];
      this._rememberId(this._ids[i], this._volumes[i], this._images[i], false);
      this._emit("stale");
    }
    return this;
  };

  TrajectoryVolumeState.prototype.clearStale = function () {
    this._stale = {};
    this._emit("stale");
    return this;
  };

  /**
   * Drop durable id→media entries that are not in ``keepIds`` so a mode handoff
   * cannot rematerialize stale catalog frames onto new path interiors.
   */
  TrajectoryVolumeState.prototype.forgetIdsExcept = function (keepIds) {
    var keep = Object.create(null);
    (keepIds || []).forEach(function (id) {
      id = asId(id);
      if (id) keep[id] = true;
    });
    Object.keys(this._byId).forEach(function (id) {
      if (!keep[id]) delete this._byId[id];
    }, this);
    return this;
  };

  /**
   * Reorder / expand / shrink slots to match targetIds. Media is looked up by
   * durable id (and the live arrays), so reverse / visit-order only change the
   * view — never the decode/render status of a volume identity.
   */
  TrajectoryVolumeState.prototype.alignToIds = function (targetIds) {
    targetIds = (targetIds || []).map(asId);
    // Harvest current slot media into the id map before rebuilding the view.
    this._rememberSlots();
    var oldIds = this._ids.slice();
    var oldVols = this._volumes.slice();
    var oldImgs = this._images.slice();
    var oldMeta = this._meta.slice();
    var oldMatch = this._match.slice();
    var oldStale = this._stale;
    for (var oi = 0; oi < oldIds.length; oi++) {
      this._rememberId(
        oldIds[oi],
        oldVols[oi] || null,
        oldImgs[oi] || null,
        !!oldStale[oi]
      );
    }
    var built = this._materializeFromById(targetIds);
    var builtMeta = new Array(targetIds.length);
    var builtMatch = new Array(targetIds.length);
    // Fall back to positional old arrays when an id was never remembered
    // (null ids / brand-new slots).
    var used = {};
    var lengthStable = targetIds.length === oldIds.length;
    var oldXy = this._xy.slice();
    for (var i = 0; i < targetIds.length; i++) {
      var want = targetIds[i];
      builtMeta[i] = cloneMeta(null, want, built.volumes[i]);
      builtMatch[i] = {
        matched: true,
        labelId: builtMeta[i].membership || want || null,
        plot_row: builtMeta[i].plot_row
      };
      if (want && (built.volumes[i] || built.images[i])) continue;
      if (!want) {
        // Null-id interiors: preserve by index ONLY when path length is stable.
        // Densify (4→10) must use realignToPathXy instead — index meaning moves.
        if (lengthStable && i < oldIds.length && oldIds[i] == null) {
          if (!built.volumes[i] && oldVols[i]) built.volumes[i] = oldVols[i];
          if (!built.images[i] && oldImgs[i]) built.images[i] = oldImgs[i];
          if (oldStale[i]) built.stale[i] = true;
          builtMeta[i] = cloneMeta(oldMeta[i], targetIds[i], built.volumes[i]);
          builtMatch[i] = cloneMatch(oldMatch[i] || {
            matched: true,
            labelId: builtMeta[i].membership || targetIds[i] || null,
            plot_row: builtMeta[i].plot_row
          });
        }
        continue;
      }
      for (var j = 0; j < oldIds.length; j++) {
        if (used[j]) continue;
        if (oldIds[j] != null && String(oldIds[j]) === want) {
          used[j] = true;
          if (!built.volumes[i] && oldVols[j]) built.volumes[i] = oldVols[j];
          if (!built.images[i] && oldImgs[j]) built.images[i] = oldImgs[j];
          if (oldStale[j]) built.stale[i] = true;
          builtMeta[i] = cloneMeta(oldMeta[j], targetIds[i], built.volumes[i]);
          builtMatch[i] = cloneMatch(oldMatch[j] || {
            matched: true,
            labelId: builtMeta[i].membership || targetIds[i] || null,
            plot_row: builtMeta[i].plot_row
          });
          break;
        }
      }
    }
    this._ids = targetIds.slice();
    this._volumes = built.volumes;
    this._images = built.images;
    this._meta = builtMeta;
    this._match = builtMatch;
    this._stale = built.stale;
    if (lengthStable) {
      this._xy = oldXy.slice(0, targetIds.length);
      while (this._xy.length < targetIds.length) this._xy.push(null);
    } else {
      // Length changed without geometry remap — clear stale XY; caller should
      // stampPathXy / realignToPathXy next.
      this._xy = new Array(targetIds.length);
      for (var xi = 0; xi < targetIds.length; xi++) this._xy[xi] = null;
    }
    this._rememberSlots();
    this._layout = {
      kind: targetIds.length ? "path" : "empty",
      pathN: targetIds.length,
      nAnchors: this._layout.nAnchors || 0,
      nInterp: this._layout.nInterp || 0,
      anchorSlotIndices: (this._layout.anchorSlotIndices || []).slice()
    };
    this._emit("align");
    return this;
  };

  /**
   * Align volume slots to the current trajectory path. Once interpolation is
   * armed, slots are always path-indexed; compact catalog media is a seed only.
   */
  TrajectoryVolumeState.prototype.alignToPath = function (path, opts) {
    opts = opts || {};
    var explicitPathN = finiteInt(opts.pathN, 0);
    // Explicit pathN is authoritative (Trajectory-controls shrink 10→4). Do not
    // floor on the previous slot count — that ratcheted the slider at the old N.
    var pathN = explicitPathN >= 2
      ? explicitPathN
      : Math.max(pathSampleCount(path), this._ids.length);
    var compactIds = opts.compactIds || opts.ids || null;
    if (!compactIds && path && typeof path.getSelectedVolumeIds === "function") {
      compactIds = path.getSelectedVolumeIds();
    }
    compactIds = (compactIds || []).map(asId);
    var nAnchors = Math.max(
      finiteInt(opts.nAnchors, 0),
      compactIds.length
    );
    var nInterp = Math.max(0, finiteInt(opts.nInterp, 0));
    if (pathN < 2 && nAnchors >= 2) {
      pathN = nInterp > 0 ? ((nAnchors - 1) * (nInterp + 1) + 1) : nAnchors;
    }
    if (pathN < 1) return this.clear();
    var slots = anchorSlotIndices(nAnchors, nInterp, pathN);
    var targetIds = new Array(pathN);
    for (var i = 0; i < pathN; i++) targetIds[i] = null;
    if (compactIds.length === pathN && (!nAnchors || pathN === nAnchors)) {
      targetIds = compactIds.slice();
    } else if (compactIds.length) {
      for (var ai = 0; ai < slots.length && ai < compactIds.length; ai++) {
        targetIds[slots[ai]] = compactIds[ai];
      }
    } else {
      for (var ti = 0; ti < pathN && ti < this._ids.length; ti++) {
        targetIds[ti] = this._ids[ti];
      }
    }

    var pathXY = normalizePathXY(opts.pathXY);
    if (pathXY.length < 2 && path && typeof path.getPathXY === "function") {
      pathXY = normalizePathXY(path.getPathXY());
    }
    var pathSamples = normalizePathSamples(opts.pathSamples || null, null);
    if ((!pathSamples || pathSamples.length < 2) && pathXY.length >= 2) {
      pathSamples = pathXY.map(function (xy, idx) {
        return { xy: xy, membership: targetIds[idx] || null };
      });
    } else if (pathSamples && pathSamples.length >= 2) {
      for (var ps = 0; ps < pathSamples.length && ps < targetIds.length; ps++) {
        if (!pathSamples[ps].membership && targetIds[ps]) {
          pathSamples[ps].membership = targetIds[ps];
        }
      }
    }
    var lengthChanged = this._ids.length >= 2 && this._ids.length !== pathN;
    if (pathSamples.length === pathN && (lengthChanged || opts.remapByPathXy || opts.rematchPath)) {
      // Geometry-owned remap: keep media only at samples whose latent XY still
      // matches a previously decoded coordinate (densify / shrink).
      this.rematchToPath(pathSamples, {
        ids: targetIds,
        compactIds: compactIds
      });
    } else {
      var oldIds = this._ids.slice();
      var oldVols = this._volumes.slice();
      var oldImgs = this._images.slice();
      var oldStale = Object.assign({}, this._stale || {});
      this.alignToIds(targetIds);
      // Length-stable: always preserve prior stale flags by index. Catalog-id
      // endpoints (pc1:*) were previously skipped here, so seedCompactCatalog
      // could drop free-drag stale and reseed old PNGs into VolumeState only —
      // zeroing ChimeraX debt while the page display stayed blank.
      if (oldIds.length === targetIds.length) {
        for (var pi = 0; pi < targetIds.length; pi++) {
          if (oldStale[pi]) this._stale[pi] = true;
          if (targetIds[pi] != null) continue;
          if (pi >= oldIds.length || oldIds[pi] != null) continue;
          if (!this._volumes[pi] && oldVols[pi]) this._volumes[pi] = cloneVol(oldVols[pi]);
          if (!this._images[pi] && oldImgs[pi]) this._images[pi] = cloneImage(oldImgs[pi]);
        }
      }
      if (pathXY.length === pathN) this.stampPathXy(pathXY);
      if (pathSamples.length === pathN) {
        this._ensureSlotArrays(pathN);
        for (var mi = 0; mi < pathN; mi++) {
          this._meta[mi] = cloneMeta(pathSamples[mi], this._ids[mi], this._volumes[mi]);
          this._match[mi] = {
            matched: !(this._match[mi] && this._match[mi].matched === false),
            labelId: this._meta[mi].membership || this._ids[mi] || null,
            plot_row: this._meta[mi].plot_row
          };
        }
      }
    }
    this._layout = {
      kind: nAnchors >= 2 && pathN > nAnchors ? "sparse-path" : "path",
      pathN: pathN,
      nAnchors: nAnchors,
      nInterp: nInterp,
      anchorSlotIndices: slots.slice()
    };
    this._emit("layout");
    return this;
  };

  TrajectoryVolumeState.prototype.seedCompactCatalog = function (compactIds, compactVolumes, compactImages, opts) {
    opts = opts || {};
    compactIds = (compactIds || []).map(asId);
    compactVolumes = (compactVolumes || []).map(cloneVol);
    compactImages = (compactImages || []).map(cloneImage);
    var requestedAnchors = finiteInt(opts.nAnchors, 0);
    var nInterp = Math.max(0, finiteInt(opts.nInterp, this._layout.nInterp || 0));
    var pathN = Math.max(
      finiteInt(opts.pathN, 0),
      this._layout.pathN || 0,
      this._ids.length || 0
    );
    // Dense path-length arrays are not compact catalogs. For two endpoints,
    // take geometric ends before nAnchors is derived from compactIds.length
    // (otherwise nAnchors becomes pathN and compact[1] lands on the path end).
    if (requestedAnchors === 2 && pathN > 2) {
      if (compactImages.length === pathN) {
        compactImages = [compactImages[0], compactImages[pathN - 1]];
      }
      if (compactVolumes.length === pathN) {
        compactVolumes = [compactVolumes[0], compactVolumes[pathN - 1]];
      }
      if (compactIds.length === pathN) {
        compactIds = [compactIds[0], compactIds[pathN - 1]];
      }
    }
    var nAnchors = Math.max(requestedAnchors, compactIds.length);
    if (pathN < 2 && nAnchors >= 2) {
      pathN = nInterp > 0 ? ((nAnchors - 1) * (nInterp + 1) + 1) : nAnchors;
    }
    if (compactIds.length === 0 && nAnchors > 0) {
      compactIds = new Array(nAnchors);
      for (var ci = 0; ci < nAnchors; ci++) compactIds[ci] = null;
    }
    var staleBeforeAlign = Object.assign({}, this._stale || {});
    this.alignToPath(null, {
      pathN: pathN,
      nAnchors: nAnchors,
      nInterp: nInterp,
      compactIds: compactIds,
      pathSamples: opts.pathSamples || null
    });
    // alignToPath / alignToIds can drop index-stable stale on catalog-id
    // endpoints (pc1:*). Restore any pre-align stale flags so freely moved
    // endpoints stay in ChimeraX debt and are not reseeded below.
    Object.keys(staleBeforeAlign).forEach(function (key) {
      if (!staleBeforeAlign[key]) return;
      var si = Math.floor(Number(key));
      if (Number.isFinite(si) && si >= 0) this._stale[si] = true;
    }, this);
    var slots = this._layout.anchorSlotIndices;
    for (var i = 0; i < slots.length && i < nAnchors; i++) {
      var slot = slots[i];
      if (slot < 0 || slot >= this._ids.length) continue;
      // Freely dragged endpoints stay stale until Decode+Render at the new XY.
      // Reseeding catalog media here previously cleared stale and zeroed debt.
      if (this._stale[slot]) continue;
      // After a real VTK decode at a new XY (volume_b64) with no ChimeraX yet,
      // do not paint old catalog PNGs — that zeroed render debt while the page
      // display still nulled invalidated ticks.
      if (this._volumes[slot] && this._volumes[slot].volume_b64 && !this._images[slot]) {
        continue;
      }
      if (i < compactIds.length) this._ids[slot] = compactIds[i];
      if (i < compactVolumes.length && compactVolumes[i]) this._volumes[slot] = cloneVol(compactVolumes[i]);
      if (i < compactImages.length && compactImages[i]) this._images[slot] = cloneImage(compactImages[i]);
      this._meta[slot] = cloneMeta(this._meta[slot], this._ids[slot], this._volumes[slot]);
      this._match[slot] = {
        matched: true,
        labelId: this._meta[slot].membership || this._ids[slot] || null,
        plot_row: this._meta[slot].plot_row
      };
      delete this._stale[slot];
    }
    this._emit("seed-catalog");
    return this;
  };

  TrajectoryVolumeState.prototype.applyDecodeResult = function (res) {
    res = res || {};
    if (res.cacheId || res.volume_cache_id) {
      this._cacheId = String(res.cacheId || res.volume_cache_id || "");
      this._generated = true;
    }
    var slotIndices = res.slot_indices || res.slotIndices || null;
    if (slotIndices && slotIndices.length) this.setSubsetSlotIndices(slotIndices);
    if (res.volumesByIndex) {
      Object.keys(res.volumesByIndex).forEach(function (key) {
        var idx = finiteInt(key, -1);
        if (idx >= 0) this.setDecoded(idx, res.volumesByIndex[key]);
      }, this);
    }
    if (res.volumesById) {
      var ids = this._ids.slice();
      Object.keys(res.volumesById).forEach(function (id) {
        var idx = ids.indexOf(String(id));
        if (idx >= 0) this.setDecoded(idx, res.volumesById[id]);
      }, this);
    }
    var vols = res.volumes || null;
    if (Array.isArray(vols)) {
      // slot_indices maps a *compact* volumes list onto path slots. A full-path
      // (or otherwise length-mismatched) array must stay positional — otherwise
      // toPayload()'s leftover subset slot_indices remaps volumes[0] onto slot 1.
      var useVolSlotMap = Array.isArray(slotIndices)
        && slotIndices.length > 0
        && slotIndices.length === vols.length;
      for (var i = 0; i < vols.length; i++) {
        var vol = vols[i];
        if (!vol) continue;
        var idx = vol.index != null ? finiteInt(vol.index, i) : (
          useVolSlotMap ? finiteInt(slotIndices[i], i) : i
        );
        if (idx >= 0) {
          this.setDecoded(idx, vol);
        }
      }
    }
    if (this.decodedCount() >= 2 || this._cacheId) this._generated = true;
    this._emit("decode-result");
    return this;
  };

  TrajectoryVolumeState.prototype.applyRenderResult = function (res) {
    res = res || {};
    var slotIndices = res.slot_indices || res.slotIndices || null;
    var pathSlots = this._ids.length;
    if (res.imagesByIndex) {
      Object.keys(res.imagesByIndex).forEach(function (key) {
        var idx = finiteInt(key, -1);
        if (idx < 0) return;
        if (pathSlots >= 2 && idx >= pathSlots) return;
        this.setRendered(idx, res.imagesByIndex[key]);
      }, this);
    }
    if (res.imagesById) {
      var ids = this._ids.slice();
      Object.keys(res.imagesById).forEach(function (id) {
        var idx = ids.indexOf(String(id));
        if (idx >= 0) this.setRendered(idx, res.imagesById[id]);
      }, this);
    }
    var imgs = res.images || null;
    if (Array.isArray(imgs)) {
      // Compact PNG lists align with slot_indices 1:1. Full-path image arrays
      // (length === path, or any length mismatch) are positional — reusing a
      // stale subset slot_indices list would copy images[0] onto slot 1 and
      // make adjacent direct-trace ticks show the same ChimeraX frame.
      var useImgSlotMap = Array.isArray(slotIndices)
        && slotIndices.length > 0
        && slotIndices.length === imgs.length;
      for (var i = 0; i < imgs.length; i++) {
        var img = cloneImage(imgs[i]);
        if (!img) continue;
        var idx = useImgSlotMap ? finiteInt(slotIndices[i], i) : i;
        if (idx < 0) continue;
        // Positional dense arrays must not grow past the current path (stale
        // PC1×10 finishes after nPoints shrink). Explicit slot_indices may
        // still target indices within the current path only here; sync the
        // path length first when densifying.
        if (pathSlots >= 2 && idx >= pathSlots) continue;
        this.setRendered(idx, img);
      }
    }
    this._emit("render-result");
    return this;
  };

  TrajectoryVolumeState.prototype.permute = function (permutation) {
    if (!permutation || permutation.length < 2) return this;
    this._rememberSlots();
    var n = permutation.length;
    this._ids = permuteArray(this._ids, permutation);
    this._volumes = permuteArray(this._volumes, permutation).map(cloneVol);
    this._images = permuteArray(this._images, permutation).map(cloneImage);
    this._xy = permuteArray(this._xy, permutation).map(cloneXy);
    this._meta = permuteArray(this._meta, permutation).map(function (m, i) {
      return cloneMeta(m, this._ids[i], this._volumes[i]);
    }, this);
    this._match = permuteArray(this._match, permutation).map(function (m, i) {
      var meta = cloneMeta(this._meta[i], this._ids[i], this._volumes[i]);
      var out = cloneMatch(m || {
        matched: true,
        labelId: meta.membership || this._ids[i] || null,
        plot_row: meta.plot_row
      });
      if (!out.labelId) out.labelId = meta.membership || this._ids[i] || null;
      if (out.plot_row == null) out.plot_row = meta.plot_row;
      return out;
    }, this);
    var nextStale = {};
    for (var key in this._stale) {
      if (!this._stale[key]) continue;
      var src = Math.floor(Number(key));
      if (!Number.isFinite(src) || src < 0 || src >= n) continue;
      for (var di = 0; di < n; di++) {
        if (permutation[di] === src) {
          nextStale[di] = true;
          break;
        }
      }
    }
    this._stale = nextStale;
    // path_t / traj_xy must follow the new slot index after permute / reverse.
    for (var i = 0; i < n; i++) {
      if (this._volumes[i]) {
        this._volumes[i] = stampVolGeometry(
          this._volumes[i],
          this._xy[i],
          pathParamT(i, n)
        );
      }
    }
    this._rememberSlots();
    this._layout = Object.assign({}, this._layout, { pathN: n });
    this._emit("permute");
    return this;
  };

  /**
   * Reverse a densified compact catalog: extract anchor media from a pre-reverse
   * snapshot (or current anchor slots), reverse the compact arrays, and
   * re-seed onto the path layout. Used when undecoded interpolated catalogs
   * cannot use a simple slot permutation.
   */
  TrajectoryVolumeState.prototype.reverseCompactCatalog = function (opts) {
    opts = opts || {};
    var nAnchors = finiteInt(opts.nAnchors, this._layout.nAnchors || 0);
    var nInterp = Math.max(0, finiteInt(opts.nInterp, this._layout.nInterp || 0));
    var pathN = Math.max(this._ids.length, finiteInt(opts.pathN, 0));
    if (nAnchors < 2 || pathN < 2) return this;

    var snap = opts.snapshot || null;
    var vols = snap ? (snap.volumes || snap.slots || []) : this._volumes;
    var imgs = snap ? (snap.images || []) : this._images;
    var snapIds = snap && snap.ids ? snap.ids : null;
    var formulaN = nAnchors < 2 ? 0 : ((nAnchors - 1) * (nInterp + 1) + 1);
    var compactVols = [];
    var compactImgs = [];
    var compactIds = [];

    if (pathN > nAnchors && (vols.length === pathN || imgs.length === pathN
        || vols.length === formulaN || imgs.length === formulaN)) {
      var slots = anchorSlotIndices(nAnchors, nInterp, pathN);
      for (var ai = 0; ai < slots.length && ai < nAnchors; ai++) {
        var slot = slots[ai];
        compactVols.push(slot >= 0 && slot < vols.length ? cloneVol(vols[slot]) : null);
        compactImgs.push(slot >= 0 && slot < imgs.length ? cloneImage(imgs[slot]) : null);
        if (snapIds && slot >= 0 && slot < snapIds.length) {
          compactIds.push(asId(snapIds[slot]));
        } else if (slot >= 0 && slot < this._ids.length) {
          compactIds.push(this._ids[slot]);
        } else {
          compactIds.push(null);
        }
      }
    } else if (vols.length === nAnchors) {
      compactVols = vols.map(cloneVol);
      compactImgs = (imgs || []).slice(0, nAnchors).map(cloneImage);
      compactIds = snapIds
        ? snapIds.slice(0, nAnchors).map(asId)
        : this._ids.slice(0, nAnchors);
    } else {
      compactVols = vols.slice(0, nAnchors).map(cloneVol);
      compactImgs = (imgs || []).slice(0, nAnchors).map(cloneImage);
      compactIds = snapIds
        ? snapIds.slice(0, nAnchors).map(asId)
        : this._ids.slice(0, nAnchors);
    }

    compactVols.reverse();
    compactImgs.reverse();
    compactIds.reverse();

    return this.seedCompactCatalog(compactIds, compactVols, compactImgs, {
      nAnchors: nAnchors,
      nInterp: nInterp,
      pathN: pathN,
      pathSamples: opts.pathSamples || null
    });
  };

  /**
   * Reverse path order. Slot arrays reverse together so null-id interiors keep
   * their media↔geometry pairing; path_t / traj_xy are restamped to the new
   * indices (alignToIds alone would leave interiors fixed while XY flipped).
   */
  TrajectoryVolumeState.prototype.reverse = function () {
    var n = this._ids.length;
    if (n < 2) return this;
    var perm = new Array(n);
    for (var i = 0; i < n; i++) perm[i] = n - 1 - i;
    return this.permute(perm);
  };

  TrajectoryVolumeState.prototype.insertSlot = function (index, id) {
    index = Math.max(0, Math.min(this._ids.length, Math.floor(Number(index)) || 0));
    this._ids.splice(index, 0, asId(id));
    this._volumes.splice(index, 0, null);
    this._images.splice(index, 0, null);
    this._xy.splice(index, 0, null);
    this._meta.splice(index, 0, cloneMeta(null, id, null));
    this._match.splice(index, 0, {
      matched: true,
      labelId: asId(id),
      plot_row: null
    });
    var nextStale = {};
    for (var key in this._stale) {
      if (!this._stale[key]) continue;
      var src = Math.floor(Number(key));
      nextStale[src >= index ? src + 1 : src] = true;
    }
    this._stale = nextStale;
    this._emit("insert");
    return this;
  };

  TrajectoryVolumeState.prototype.removeSlots = function (indices) {
    var remove = {};
    (indices || []).forEach(function (i) {
      i = Math.floor(Number(i));
      if (Number.isFinite(i) && i >= 0) remove[i] = true;
    });
    var ids = [];
    var vols = [];
    var imgs = [];
    var xys = [];
    var metas = [];
    var matches = [];
    var stale = {};
    var newIdx = 0;
    for (var i = 0; i < this._ids.length; i++) {
      if (remove[i]) continue;
      ids.push(this._ids[i]);
      vols.push(this._volumes[i] || null);
      imgs.push(this._images[i] || null);
      xys.push(this._xy[i] || null);
      metas.push(cloneMeta(this._meta[i], this._ids[i], this._volumes[i]));
      matches.push(cloneMatch(this._match[i] || {
        matched: true,
        labelId: this._ids[i] || null,
        plot_row: null
      }));
      if (this._stale[i]) stale[newIdx] = true;
      newIdx++;
    }
    this._ids = ids;
    this._volumes = vols;
    this._images = imgs;
    this._xy = xys;
    this._meta = metas;
    this._match = matches;
    this._stale = stale;
    this._emit("remove");
    return this;
  };

  TrajectoryVolumeState.prototype.clear = function () {
    this._ids = [];
    this._volumes = [];
    this._images = [];
    this._xy = [];
    this._meta = [];
    this._match = [];
    this._stale = {};
    this._byId = Object.create(null);
    this._generated = false;
    this._cacheId = "";
    this._subsetSlotIndices = null;
    this._layout = {
      kind: "empty",
      pathN: 0,
      nAnchors: 0,
      nInterp: 0,
      anchorSlotIndices: []
    };
    this._emit("clear");
    return this;
  };

  /**
   * Snapshot for bridging legacy page globals.
   */
  TrajectoryVolumeState.prototype.snapshot = function () {
    return {
      ids: this.ids(),
      volumes: this.volumes(),
      images: this.images(),
      xy: this.pathXy(),
      meta: this.meta(),
      match: this.match(),
      stale: Object.assign({}, this._stale),
      generated: this._generated,
      cacheId: this._cacheId,
      subsetSlotIndices: this.subsetSlotIndices(),
      layout: this.layout(),
      ready: this.decodedCount() > 0 || this.renderedCount() > 0
    };
  };

  TrajectoryVolumeState.prototype.loadSnapshot = function (snap) {
    if (!snap) return this.clear();
    this.replaceSlots(snap.ids, snap.volumes, snap.images, snap.xy);
    if (snap.meta && snap.meta.length) {
      this._meta = snap.meta.map(function (m, i) {
        return cloneMeta(m, this._ids[i], this._volumes[i]);
      }, this);
    }
    if (snap.match && snap.match.length) {
      this._match = snap.match.map(function (m, i) {
        var meta = cloneMeta(this._meta[i], this._ids[i], this._volumes[i]);
        var out = cloneMatch(m);
        if (!out.labelId) out.labelId = meta.membership || this._ids[i] || null;
        if (out.plot_row == null) out.plot_row = meta.plot_row;
        return out;
      }, this);
    }
    this._stale = Object.assign({}, snap.stale || {});
    this._generated = !!snap.generated;
    this._cacheId = String(snap.cacheId || "");
    this._subsetSlotIndices = snap.subsetSlotIndices ? snap.subsetSlotIndices.slice() : null;
    if (snap.layout) {
      this._layout = {
        kind: snap.layout.kind || "path",
        pathN: finiteInt(snap.layout.pathN, this._ids.length),
        nAnchors: finiteInt(snap.layout.nAnchors, 0),
        nInterp: finiteInt(snap.layout.nInterp, 0),
        anchorSlotIndices: (snap.layout.anchorSlotIndices || []).slice()
      };
    }
    this._emit("load");
    return this;
  };

  global.CryoTrajectoryVolumeState = TrajectoryVolumeState;
  global.CryoTrajectoryVolumeStateUtils = {
    permuteArray: permuteArray,
    anchorSlotIndices: anchorSlotIndices,
    compactCatalogPackedPrefix: compactCatalogPackedPrefix,
    pathParamT: pathParamT,
    normalizePathXY: normalizePathXY,
    xyDist2: xyDist2,
    latentXyMatchEps2: latentXyMatchEps2,
    latentXyMatchAtol: latentXyMatchAtol,
    xyMatchesLatent: xyMatchesLatent
  };
})(typeof window !== "undefined" ? window : this);
