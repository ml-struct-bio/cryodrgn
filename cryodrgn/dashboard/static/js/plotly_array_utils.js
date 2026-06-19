/**
 * Helpers for Plotly trace arrays after server serialization.
 *
 * Plotly 5 emits plain JSON lists; Plotly 6 may binary-encode numpy-backed arrays as
 * ``{dtype, bdata[, shape], _inputArray}`` objects or TypedArray views without a
 * reliable ``.length``. ``shape`` may be a list or a comma-separated string.
 */
(function(global) {
  "use strict";

  var _binaryDecodeCache = typeof WeakMap !== "undefined" ? new WeakMap() : null;

  function isPlotlyBinary(arr) {
    return !!(arr && (arr.bdata || arr.dtype));
  }

  function normalizeDtype(dtype) {
    var d = String(dtype || "f8").trim();
    if (d.charAt(0) === "<" || d.charAt(0) === ">" || d.charAt(0) === "|" || d.charAt(0) === "=") {
      d = d.slice(1);
    }
    return d;
  }

  function parseShape(shape) {
    if (shape == null) return null;
    if (Array.isArray(shape)) return shape;
    if (typeof shape === "string") {
      var parts = shape.split(",").map(function(s) {
        return parseInt(String(s).trim(), 10);
      });
      if (parts.length && parts.every(function(n) { return isFinite(n); })) return parts;
      return null;
    }
    if (typeof shape === "number" && isFinite(shape)) return [shape];
    if (typeof shape.length === "number" && shape.length > 0) {
      try {
        return Array.from(shape);
      } catch (e) { /* ignore */ }
    }
    return null;
  }

  function shapeDims(arr) {
    if (!arr) return null;
    return parseShape(arr.shape);
  }

  function shapeRowCount(arr) {
    var shape = shapeDims(arr);
    if (!shape || !shape.length) return 0;
    return shape[0];
  }

  function shapeColCount(arr) {
    var shape = shapeDims(arr);
    if (!shape || shape.length < 2) return 0;
    return shape[1];
  }

  function inputArrayRowCount(input) {
    if (!input || typeof input.length !== "number") return 0;
    if (input.length !== 1) return input.length;
    var inner = input[0];
    if (inner == null) return 1;
    if (Array.isArray(inner)) return inner.length;
    if (typeof inner.length === "number") return inner.length;
    return 1;
  }

  function length(arr) {
    if (!arr) return 0;
    if (typeof arr === "string") return 0;
    var fromShape = shapeRowCount(arr);
    if (fromShape > 0) return fromShape;
    if (isPlotlyBinary(arr)) {
      var inputLen = inputArrayRowCount(arr._inputArray);
      if (inputLen > 0) return inputLen;
      var flat = decodedBinaryFlat(arr);
      if (flat) return flat.length;
      return 0;
    }
    if (typeof arr.length === "number") return arr.length;
    return inputArrayRowCount(arr._inputArray);
  }

  function isArray(arr) {
    if (!arr) return false;
    if (typeof arr === "string") return false;
    if (Array.isArray(arr)) return true;
    if (shapeRowCount(arr) > 0) return true;
    if (isPlotlyBinary(arr)) return true;
    if (arr._inputArray && typeof arr._inputArray.length === "number") return true;
    if (typeof arr.length === "number" && typeof arr[0] !== "undefined") return true;
    return false;
  }

  function base64ToUint8Array(b64) {
    var binary = atob(String(b64));
    var len = binary.length;
    var bytes = new Uint8Array(len);
    for (var i = 0; i < len; i++) bytes[i] = binary.charCodeAt(i);
    return bytes;
  }

  function typedArrayFromBytes(bytes, dtype) {
    dtype = normalizeDtype(dtype);
    try {
      if (dtype === "f8") return new Float64Array(bytes.buffer, bytes.byteOffset, Math.floor(bytes.byteLength / 8));
      if (dtype === "f4") return new Float32Array(bytes.buffer, bytes.byteOffset, Math.floor(bytes.byteLength / 4));
      if (dtype === "f2" && typeof Float16Array !== "undefined") {
        return new Float16Array(bytes.buffer, bytes.byteOffset, Math.floor(bytes.byteLength / 2));
      }
      if (dtype === "i8" && typeof BigInt64Array !== "undefined") {
        return new BigInt64Array(bytes.buffer, bytes.byteOffset, Math.floor(bytes.byteLength / 8));
      }
      if (dtype === "u8" && typeof BigUint64Array !== "undefined") {
        return new BigUint64Array(bytes.buffer, bytes.byteOffset, Math.floor(bytes.byteLength / 8));
      }
      if (dtype === "i4") return new Int32Array(bytes.buffer, bytes.byteOffset, Math.floor(bytes.byteLength / 4));
      if (dtype === "u4") return new Uint32Array(bytes.buffer, bytes.byteOffset, Math.floor(bytes.byteLength / 4));
      if (dtype === "i2") return new Int16Array(bytes.buffer, bytes.byteOffset, Math.floor(bytes.byteLength / 2));
      if (dtype === "u2") return new Uint16Array(bytes.buffer, bytes.byteOffset, Math.floor(bytes.byteLength / 2));
      if (dtype === "i1") return new Int8Array(bytes.buffer, bytes.byteOffset, bytes.byteLength);
      if (dtype === "u1") return new Uint8Array(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    } catch (e) { /* ignore */ }
    return null;
  }

  function decodedBinaryFlat(arr) {
    if (!arr || !arr.bdata) return null;
    if (_binaryDecodeCache && _binaryDecodeCache.has(arr)) {
      return _binaryDecodeCache.get(arr);
    }
    var bytes = base64ToUint8Array(arr.bdata);
    var typed = typedArrayFromBytes(bytes, arr.dtype);
    if (typed && _binaryDecodeCache) _binaryDecodeCache.set(arr, typed);
    return typed;
  }

  function numericFromTyped(val) {
    if (typeof val === "number" && isFinite(val)) return val;
    if (typeof val === "bigint") return Number(val);
    return val;
  }

  function slice(arr) {
    if (!arr) return [];
    if (Array.isArray(arr)) return arr.slice();
    if (isArray(arr)) {
      try {
        return Array.prototype.slice.call(arr);
      } catch (e) {
        try {
          return Array.from(arr);
        } catch (e2) {
          var flat = decodedBinaryFlat(arr);
          if (flat) return Array.prototype.slice.call(flat);
          return [];
        }
      }
    }
    return [];
  }

  function rowAsArray(row) {
    if (row == null) return [];
    if (Array.isArray(row)) return row;
    if (typeof row === "object" && row.length !== undefined) {
      try {
        return Array.from(row);
      } catch (e) {
        return [];
      }
    }
    return [row];
  }

  function rowAt(matrix, rowIdx) {
    if (matrix == null || rowIdx == null || rowIdx < 0) return [];
    try {
      if (matrix[rowIdx] != null) return rowAsArray(matrix[rowIdx]);
    } catch (e) { /* ignore */ }
    var input = matrix._inputArray;
    if (input) {
      try {
        if (input[rowIdx] != null) return rowAsArray(input[rowIdx]);
        if (input.length === 1 && input[0] != null && input[0][rowIdx] != null) {
          return rowAsArray(input[0][rowIdx]);
        }
      } catch (e2) { /* ignore */ }
    }
    var flat = decodedBinaryFlat(matrix);
    if (flat) {
      var cols = shapeColCount(matrix);
      if (cols > 0) {
        var start = rowIdx * cols;
        if (start + cols <= flat.length) {
          var row = [];
          for (var c = 0; c < cols; c++) row.push(numericFromTyped(flat[start + c]));
          return row;
        }
      }
      if (rowIdx < flat.length) return [numericFromTyped(flat[rowIdx])];
    }
    return [];
  }

  function rowLength(row) {
    return rowAsArray(row).length;
  }

  function rowValue(source, idx, colIdx) {
    if (colIdx !== undefined) {
      return rowValue(rowAt(source, idx), colIdx);
    }
    var arr = rowAsArray(source);
    return idx >= 0 && idx < arr.length ? arr[idx] : undefined;
  }

  function valueAt(arr, idx) {
    if (!arr || idx == null || idx < 0) return undefined;
    try {
      var direct = arr[idx];
      if (typeof direct === "number" && isFinite(direct)) return direct;
      if (typeof direct === "bigint") return Number(direct);
      if (direct !== undefined && typeof direct !== "object") return direct;
    } catch (e) { /* ignore */ }
    var input = arr._inputArray;
    if (input) {
      try {
        var fromInput = input[idx];
        if (typeof fromInput === "number" && isFinite(fromInput)) return fromInput;
        if (typeof fromInput === "bigint") return Number(fromInput);
        if (fromInput !== undefined && typeof fromInput !== "object") return fromInput;
      } catch (e2) { /* ignore */ }
    }
    var flat = decodedBinaryFlat(arr);
    if (flat && idx < flat.length) return numericFromTyped(flat[idx]);
    return undefined;
  }

  /** Prefer decoded ``_fullData`` arrays, then ``data`` trace arrays. */
  function traceValueAt(gd, traceIndex, key, pointIndex) {
    var sources = [];
    if (gd && gd._fullData && gd._fullData[traceIndex]) {
      sources.push(gd._fullData[traceIndex][key]);
    }
    if (gd && gd.data && gd.data[traceIndex]) {
      sources.push(gd.data[traceIndex][key]);
    }
    for (var si = 0; si < sources.length; si++) {
      var v = valueAt(sources[si], pointIndex);
      if (typeof v === "number" && isFinite(v)) return v;
    }
    return undefined;
  }

  function rowsEqualLength(a, b) {
    return length(a) === length(b);
  }

  global.CryoPlotlyArrays = {
    length: length,
    isArray: isArray,
    slice: slice,
    rowAsArray: rowAsArray,
    rowAt: rowAt,
    rowLength: rowLength,
    rowValue: rowValue,
    valueAt: valueAt,
    traceValueAt: traceValueAt,
    rowsEqualLength: rowsEqualLength,
    customdataLength: length,
  };
})(typeof window !== "undefined" ? window : this);
