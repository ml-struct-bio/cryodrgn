/**
 * Client side of the compact "columns" scatter format
 * (server: ``cryodrgn/dashboard/figure_columns.py``).
 *
 * ``CryoFigureColumns.fetchFigure(url, init, slot)`` asks a scatter endpoint for
 * ``wire=columns`` and resolves to the same Plotly figure object the JSON response
 * would have produced: ``fig.data[0].x`` etc. are typed arrays, ``customdata`` is an
 * array of per-point rows, per-point colours are hex strings. Page code that reads
 * the figure does not need to know which format was used.
 *
 * ``slot`` names the plot (one per page is enough). The coordinate and index
 * columns of the last response are kept per slot; the next request tells the
 * server which ones it holds (``have``) so a colour-only change downloads only
 * the colour columns.
 */
(function () {
  "use strict";

  var MAGIC = "CDFC1\n";
  var MIMETYPE = "application/vnd.cryodrgn.figure-columns";
  var TYPED = {
    f4: Float32Array, f8: Float64Array,
    i1: Int8Array, i2: Int16Array, i4: Int32Array,
    u1: Uint8Array, u2: Uint16Array, u4: Uint32Array
  };
  /** slot -> { key: base_key, cols: { name: decoded column } } */
  var baseCache = {};

  /**
   * C ``%.5g`` of a finite double, including round-half-to-even on exact ties
   * (JavaScript's own toPrecision rounds ties away from zero).
   */
  function formatG5(v) {
    if (v === 0) return (1 / v < 0) ? "-0" : "0";
    var neg = v < 0;
    var a = neg ? -v : v;
    var e = a.toExponential(25);
    var ePos = e.indexOf("e");
    var exp = parseInt(e.slice(ePos + 1), 10);
    var digits = e.slice(0, ePos).replace(".", "");
    var head = parseInt(digits.slice(0, 5), 10);
    var rest = digits.slice(5);
    var first = rest.charCodeAt(0) - 48;
    var up;
    if (first > 5) up = true;
    else if (first < 5) up = false;
    else if (/[1-9]/.test(rest.slice(1))) up = true;
    else up = (head % 2) === 1;
    if (up) {
      head += 1;
      if (head === 100000) { head = 10000; exp += 1; }
    }
    var d = String(head);
    var out;
    if (exp < -4 || exp >= 5) {
      var frac = d.slice(1).replace(/0+$/, "");
      var ae = exp < 0 ? -exp : exp;
      out = d.charAt(0) + (frac ? "." + frac : "") + "e" + (exp < 0 ? "-" : "+") + (ae < 10 ? "0" + ae : String(ae));
    } else if (exp >= 0) {
      var ip = d.slice(0, exp + 1);
      var fp = d.slice(exp + 1).replace(/0+$/, "");
      out = ip + (fp ? "." + fp : "");
    } else {
      var body = (new Array(-exp).join("0") + d).replace(/0+$/, "");
      out = "0." + body;
    }
    return neg ? "-" + out : out;
  }

  function isColumnRef(obj) {
    if (!obj || typeof obj !== "object" || Array.isArray(obj)) return false;
    var keys = Object.keys(obj);
    return keys.length === 1 && keys[0] === "$col";
  }

  function decode(buffer, slot) {
    var bytes = new Uint8Array(buffer);
    for (var m = 0; m < MAGIC.length; m++) {
      if (bytes[m] !== MAGIC.charCodeAt(m)) throw new Error("Not a figure-columns response");
    }
    var view = new DataView(buffer);
    var hlen = view.getUint32(MAGIC.length, true);
    var h0 = MAGIC.length + 4;
    var header = JSON.parse(new TextDecoder("utf-8").decode(bytes.subarray(h0, h0 + hlen)));
    var data0 = h0 + hlen;
    data0 += (8 - (data0 % 8)) % 8;
    var specs = header.cols;
    var decoded = new Array(specs.length);
    var held = (slot && baseCache[slot] && baseCache[slot].key === header.base_key)
      ? baseCache[slot].cols : null;
    var baseOut = {};

    function typed(spec) {
      var Ctor = TYPED[spec.dtype];
      if (!Ctor) throw new Error("Unsupported column dtype " + spec.dtype);
      return new Ctor(buffer, data0 + spec.offset, spec.n);
    }

    function col(i) {
      if (decoded[i] !== undefined) return decoded[i];
      var s = specs[i];
      var out, k, n, src, table, codes;
      switch (s.kind) {
        case "cached":
          if (!held || held[s.name] === undefined) {
            var miss = new Error("figure-columns base cache miss: " + s.name);
            miss.cryoBaseCacheMiss = true;
            throw miss;
          }
          out = held[s.name];
          baseOut[s.name] = out;
          break;
        case "arange":
          out = new Int32Array(s.n);
          for (k = 0; k < s.n; k++) out[k] = k;
          break;
        case "num":
          /* Copy out of the response buffer so it can be released. */
          out = typed(s).slice();
          break;
        case "dict":
          table = s.table;
          codes = typed(s);
          n = codes.length;
          out = new Array(n);
          for (k = 0; k < n; k++) out[k] = table[codes[k]];
          break;
        case "const":
          out = new Array(s.n);
          for (k = 0; k < s.n; k++) out[k] = s.value;
          break;
        case "json":
          out = s.values;
          break;
        case "g5":
          src = col(s.src);
          n = src.length;
          out = new Array(n);
          for (k = 0; k < n; k++) {
            var v = src[k];
            out[k] = (v === null || !isFinite(v)) ? s.missing : s.prefix + formatG5(v);
          }
          break;
        case "rows":
          var members = [];
          for (k = 0; k < s.cols.length; k++) members.push(col(s.cols[k]));
          var nm = members.length;
          n = nm ? members[0].length : 0;
          out = new Array(n);
          for (var r = 0; r < n; r++) {
            var row = new Array(nm);
            for (var c = 0; c < nm; c++) {
              var x = members[c][r];
              if (s.nan_null[c] && typeof x === "number" && !isFinite(x)) x = null;
              row[c] = x;
            }
            out[r] = row;
          }
          break;
        default:
          throw new Error("Unknown figure column kind " + s.kind);
      }
      if (s.base && s.name) baseOut[s.name] = out;
      decoded[i] = out;
      return out;
    }

    function hydrate(obj, topLevel) {
      if (isColumnRef(obj)) {
        var out = col(obj.$col);
        /* Hand Plotly its own copy of cached coordinate arrays: the figure object is
           mutated by page code and must not alias what the cache holds. */
        var s = specs[obj.$col];
        if ((s.base || s.kind === "cached") && out && typeof out.slice === "function") {
          return out.slice();
        }
        return out;
      }
      if (Array.isArray(obj)) {
        for (var i = 0; i < obj.length; i++) obj[i] = hydrate(obj[i]);
        return obj;
      }
      if (obj && typeof obj === "object") {
        for (var key in obj) {
          if (Object.prototype.hasOwnProperty.call(obj, key)) obj[key] = hydrate(obj[key]);
        }
      }
      return obj;
    }

    var fig = hydrate(header.fig);
    if (slot) baseCache[slot] = { key: header.base_key, cols: baseOut };
    return fig;
  }

  function withWireParams(url, init, have) {
    init = init || {};
    var method = String(init.method || "GET").toUpperCase();
    if (method === "GET") {
      var sep = url.indexOf("?") >= 0 ? "&" : "?";
      url = url + sep + "wire=columns" + (have ? "&have=" + encodeURIComponent(have) : "");
      return { url: url, init: init };
    }
    var body = {};
    try { body = init.body ? JSON.parse(init.body) : {}; } catch (e) { body = {}; }
    body.wire = "columns";
    if (have) body.have = have; else delete body.have;
    var next = {};
    for (var k in init) { if (Object.prototype.hasOwnProperty.call(init, k)) next[k] = init[k]; }
    next.body = JSON.stringify(body);
    return { url: url, init: next };
  }

  /**
   * Fetch a scatter figure in the columns format.
   *
   * Resolves to the figure object. Rejects with ``Error(message)`` where message is
   * the server's JSON ``error`` (or the HTTP status text), like the JSON fetches
   * this replaces. A server without the columns format (plain JSON reply) works too.
   *
   * @param {string} url      endpoint (query string allowed)
   * @param {object} [init]   fetch options; a POST ``body`` must be a JSON string
   * @param {string} [slot]   cache name for this plot's coordinate columns
   * @param {string} [notJsonMessage]  error text when an error reply is not JSON
   */
  function fetchFigure(url, init, slot, notJsonMessage) {
    function attempt(useHave) {
      var have = (useHave && slot && baseCache[slot]) ? baseCache[slot].key : null;
      var req = withWireParams(url, init, have);
      return fetch(req.url, req.init).then(function (r) {
        var ctype = String(r.headers.get("Content-Type") || "");
        if (r.ok && ctype.indexOf(MIMETYPE) >= 0) {
          return r.arrayBuffer().then(function (buf) {
            try {
              return decode(buf, slot);
            } catch (err) {
              if (err && err.cryoBaseCacheMiss && useHave) {
                if (slot) delete baseCache[slot];
                return attempt(false);
              }
              throw err;
            }
          });
        }
        return r.text().then(function (text) {
          var j;
          try { j = JSON.parse(text); }
          catch (e) { throw new Error(notJsonMessage || "API did not return JSON"); }
          if (!r.ok) throw new Error(j.error || r.statusText);
          return j;
        });
      });
    }
    return attempt(true);
  }

  window.CryoFigureColumns = {
    fetchFigure: fetchFigure,
    decode: decode,
    formatG5: formatG5,
    clearCache: function (slot) {
      if (slot) delete baseCache[slot]; else baseCache = {};
    }
  };
})();
