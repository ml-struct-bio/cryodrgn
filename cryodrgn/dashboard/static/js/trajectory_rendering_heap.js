/**
 * Byte-budget LRU heap for inactive ChimeraX trajectory renderings (PNG b64).
 *
 * Keeps frames that fall out of the live slider so a later return to the same
 * volume identity + view/iso can avoid a ChimeraX round-trip.
 */
(function (global) {
  "use strict";

  var DEFAULT_MAX_BYTES = 100 * 1024 * 1024;

  function stripDataUrl(b64) {
    var s = String(b64 || "");
    var i = s.indexOf(",");
    if (i >= 0 && /^data:/i.test(s.slice(0, i + 1))) return s.slice(i + 1);
    return s;
  }

  function estimateImageBytes(b64) {
    var s = stripDataUrl(b64);
    if (!s) return 0;
    // Approximate decoded PNG size from base64; add small per-entry overhead.
    return Math.ceil(s.length * 0.75) + 64;
  }

  function RenderingHeap(opts) {
    opts = opts || {};
    this._maxBytes = Math.max(1024 * 1024, Math.floor(Number(opts.maxBytes) || DEFAULT_MAX_BYTES));
    this._entries = new Map(); // key -> { b64, bytes }
    this._usedBytes = 0;
  }

  RenderingHeap.prototype.usedBytes = function () {
    return this._usedBytes;
  };

  RenderingHeap.prototype.size = function () {
    return this._entries.size;
  };

  RenderingHeap.prototype.has = function (key) {
    return this._entries.has(String(key || ""));
  };

  RenderingHeap.prototype.get = function (key) {
    key = String(key || "");
    if (!key || !this._entries.has(key)) return null;
    var entry = this._entries.get(key);
    // LRU: refresh recency.
    this._entries.delete(key);
    this._entries.set(key, entry);
    return entry.b64;
  };

  RenderingHeap.prototype.put = function (key, b64) {
    key = String(key || "");
    b64 = stripDataUrl(b64);
    if (!key || !b64) return false;
    var bytes = estimateImageBytes(b64);
    if (bytes > this._maxBytes) return false;
    if (this._entries.has(key)) {
      var prev = this._entries.get(key);
      this._usedBytes -= prev.bytes;
      this._entries.delete(key);
    }
    this._evictToFit(bytes);
    this._entries.set(key, { b64: b64, bytes: bytes });
    this._usedBytes += bytes;
    return true;
  };

  RenderingHeap.prototype.remove = function (key) {
    key = String(key || "");
    if (!this._entries.has(key)) return false;
    var entry = this._entries.get(key);
    this._usedBytes -= entry.bytes;
    this._entries.delete(key);
    return true;
  };

  /**
   * Drop every frame keyed to a volume identity (any view/iso). Used when a
   * path sample is freely dragged so old PNGs cannot satisfy Render debt.
   */
  RenderingHeap.prototype.removeForVolume = function (volumeKey) {
    volumeKey = String(volumeKey || "");
    if (!volumeKey) return 0;
    var prefix = "cx|" + volumeKey + "|";
    var toRemove = [];
    this._entries.forEach(function (_entry, key) {
      if (key.indexOf(prefix) === 0) toRemove.push(key);
    });
    var n = 0;
    for (var i = 0; i < toRemove.length; i++) {
      if (this.remove(toRemove[i])) n++;
    }
    return n;
  };

  RenderingHeap.prototype.clear = function () {
    this._entries.clear();
    this._usedBytes = 0;
    return this;
  };

  /**
   * Look up by exact key, else the most recently used frame for the same
   * volume identity (any view/iso). Keys are "cx|<volId>|<view>|<iso>".
   */
  RenderingHeap.prototype.getForVolume = function (volumeKey, preferredKey) {
    volumeKey = String(volumeKey || "");
    preferredKey = preferredKey ? String(preferredKey) : "";
    if (preferredKey && this._entries.has(preferredKey)) {
      return this.get(preferredKey);
    }
    if (!volumeKey) return null;
    var prefix = "cx|" + volumeKey + "|";
    var match = null;
    this._entries.forEach(function (_entry, key) {
      if (key.indexOf(prefix) === 0) match = key;
    });
    return match ? this.get(match) : null;
  };

  RenderingHeap.prototype._evictToFit = function (incomingBytes) {
    incomingBytes = Math.max(0, Math.floor(Number(incomingBytes)) || 0);
    while (this._entries.size && this._usedBytes + incomingBytes > this._maxBytes) {
      var oldest = this._entries.keys().next().value;
      var entry = this._entries.get(oldest);
      this._entries.delete(oldest);
      this._usedBytes -= entry.bytes;
    }
    if (this._usedBytes < 0) this._usedBytes = 0;
  };

  global.CryoTrajectoryRenderingHeap = RenderingHeap;
})(typeof window !== "undefined" ? window : this);
