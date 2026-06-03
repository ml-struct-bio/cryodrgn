/**
 * Small shared utilities for ResizeObserver usage in the dashboard.
 *
 * Exposes `window.CryoResizeObserverUtils`:
 * - `create(callback)` -> ResizeObserver | null
 * - `observeResize(targets, callback)` -> ResizeObserver | null
 * - `observeMany(ro, targets)` -> observes elements into an existing observer
 */
(function (global) {
  "use strict";

  function asArray(targets) {
    if (!targets) return [];
    if (Array.isArray(targets)) return targets;
    if (typeof targets.length === "number") {
      return Array.prototype.slice.call(targets);
    }
    return [targets];
  }

  function create(callback) {
    if (!global || typeof global.ResizeObserver === "undefined") return null;
    if (typeof callback !== "function") return null;
    return new global.ResizeObserver(function (entries, observer) {
      callback(entries, observer);
    });
  }

  function observeMany(ro, targets) {
    if (!ro || typeof ro.observe !== "function") return;
    var arr = asArray(targets);
    for (var i = 0; i < arr.length; i++) {
      var t = arr[i];
      if (!t || t.nodeType !== 1) continue; // only Elements
      try {
        ro.observe(t);
      } catch (e) {
        // Ignore non-observable elements.
      }
    }
  }

  function observeResize(targets, callback) {
    var ro = create(callback);
    if (!ro) return null;
    observeMany(ro, targets);
    return ro;
  }

  global.CryoResizeObserverUtils = {
    create: create,
    observeMany: observeMany,
    observeResize: observeResize,
  };
})(typeof window !== "undefined" ? window : this);
