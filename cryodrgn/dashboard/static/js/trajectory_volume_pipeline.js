/**
 * Decode / render job orchestration for trajectory volumes.
 *
 * Owns generation counters, in-flight flags, and the request lifecycle.
 * Slot payloads are written into TrajectoryVolumeState; UI status is reported
 * via hooks so the page templates stay thin.
 */
(function (global) {
  "use strict";

  function VolumePipeline(opts) {
    opts = opts || {};
    this.hooks = opts.hooks || {};
    this.volumes = opts.volumes || null;
    this._decodeGen = 0;
    this._renderGen = 0;
    this._decodeInFlight = false;
    this._renderInFlight = false;
    this._activeJobId = null;
  }

  VolumePipeline.prototype.setVolumeState = function (volumes) {
    this.volumes = volumes || null;
    return this;
  };

  VolumePipeline.prototype.isBusy = function () {
    return this._decodeInFlight || this._renderInFlight;
  };

  VolumePipeline.prototype.isDecoding = function () {
    return this._decodeInFlight;
  };

  VolumePipeline.prototype.isRendering = function () {
    return this._renderInFlight;
  };

  VolumePipeline.prototype._notify = function (eventName, detail) {
    var fn = this.hooks[eventName];
    if (typeof fn === "function") {
      try { fn(detail || {}); } catch (err) { /* ignore */ }
    }
  };

  VolumePipeline.prototype.cancel = function () {
    this._decodeGen++;
    this._renderGen++;
    this._decodeInFlight = false;
    this._renderInFlight = false;
    this._activeJobId = null;
    this._notify("onCancel", {});
    return this;
  };

  /**
   * Decode missing (or all) slots. ``fetchDecode`` must return a Promise that
   * resolves to ``{ ok, volumesById|volumesByIndex, cacheId?, jobId? }``.
   */
  VolumePipeline.prototype.decode = function (opts) {
    opts = opts || {};
    var self = this;
    if (!this.volumes) {
      return Promise.resolve({ ok: false, reason: "no-volumes" });
    }
    if (typeof this.hooks.fetchDecode !== "function") {
      return Promise.resolve({ ok: false, reason: "no-fetch" });
    }
    var gen = ++this._decodeGen;
    this._decodeInFlight = true;
    var indices = opts.indices;
    if (!indices) {
      indices = opts.forceAll
        ? this.volumes.ids().map(function (_, i) { return i; })
        : (typeof this.volumes.decodeDebtIndices === "function"
          ? this.volumes.decodeDebtIndices()
          : this.volumes.missingDecodeIndices());
    }
    var ids = indices.map(function (i) { return self.volumes.slotIdAt(i); }).filter(Boolean);
    this._notify("onDecodeStart", { indices: indices, ids: ids });
    return Promise.resolve(this.hooks.fetchDecode({
      indices: indices,
      ids: ids,
      forceAll: !!opts.forceAll
    })).then(function (res) {
      if (gen !== self._decodeGen) return { ok: false, reason: "stale" };
      self._decodeInFlight = false;
      if (!res || !res.ok) {
        self._notify("onDecodeError", res || {});
        return res || { ok: false };
      }
      if (res.jobId) self._activeJobId = res.jobId;
      if (typeof self.volumes.applyDecodeResult === "function") {
        self.volumes.applyDecodeResult(res);
      } else {
        if (res.volumesByIndex) {
          Object.keys(res.volumesByIndex).forEach(function (key) {
            var idx = Math.floor(Number(key));
            if (Number.isFinite(idx)) self.volumes.setDecoded(idx, res.volumesByIndex[key]);
          });
        }
        if (res.volumesById) {
          var slotIds = self.volumes.ids();
          Object.keys(res.volumesById).forEach(function (id) {
            var idx = slotIds.indexOf(String(id));
            if (idx >= 0) self.volumes.setDecoded(idx, res.volumesById[id]);
          });
        }
        if (res.cacheId) self.volumes.setGenerated(true, res.cacheId);
        else if (self.volumes.decodedCount() >= 2) self.volumes.setGenerated(true, self.volumes.cacheId());
      }
      self._notify("onDecodeComplete", res);
      return res;
    }).catch(function (err) {
      if (gen !== self._decodeGen) return { ok: false, reason: "stale" };
      self._decodeInFlight = false;
      self._notify("onDecodeError", { ok: false, error: err });
      return { ok: false, error: err };
    });
  };

  /**
   * Render ChimeraX images for missing (or all) slots.
   */
  VolumePipeline.prototype.render = function (opts) {
    opts = opts || {};
    var self = this;
    if (!this.volumes) {
      return Promise.resolve({ ok: false, reason: "no-volumes" });
    }
    if (typeof this.hooks.fetchRender !== "function") {
      return Promise.resolve({ ok: false, reason: "no-fetch" });
    }
    var gen = ++this._renderGen;
    this._renderInFlight = true;
    var indices = opts.indices;
    if (!indices || !indices.length) {
      indices = opts.forceAll
        ? this.volumes.ids().map(function (_, i) { return i; })
        : this.volumes.missingRenderIndices();
    }
    var ids = indices.map(function (i) { return self.volumes.slotIdAt(i); }).filter(Boolean);
    this._notify("onRenderStart", { indices: indices, ids: ids });
    return Promise.resolve(this.hooks.fetchRender({
      indices: indices,
      ids: ids,
      forceAll: !!opts.forceAll
    })).then(function (res) {
      if (gen !== self._renderGen) return { ok: false, reason: "stale" };
      self._renderInFlight = false;
      if (!res || !res.ok) {
        self._notify("onRenderError", res || {});
        return res || { ok: false };
      }
      if (typeof self.volumes.applyRenderResult === "function") {
        self.volumes.applyRenderResult(res);
      } else {
        if (res.imagesByIndex) {
          Object.keys(res.imagesByIndex).forEach(function (key) {
            var idx = Math.floor(Number(key));
            if (Number.isFinite(idx)) self.volumes.setRendered(idx, res.imagesByIndex[key]);
          });
        }
        if (res.imagesById) {
          var slotIds = self.volumes.ids();
          Object.keys(res.imagesById).forEach(function (id) {
            var idx = slotIds.indexOf(String(id));
            if (idx >= 0) self.volumes.setRendered(idx, res.imagesById[id]);
          });
        }
      }
      self._notify("onRenderComplete", res);
      return res;
    }).catch(function (err) {
      if (gen !== self._renderGen) return { ok: false, reason: "stale" };
      self._renderInFlight = false;
      self._notify("onRenderError", { ok: false, error: err });
      return { ok: false, error: err };
    });
  };

  global.CryoTrajectoryVolumePipeline = VolumePipeline;
})(typeof window !== "undefined" ? window : this);
