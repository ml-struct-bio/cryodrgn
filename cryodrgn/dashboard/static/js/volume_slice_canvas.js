/**
 * Interactive single-slice viewer for decoded cryoDRGN volumes.
 *
 * Left-drag rotates the cutting plane. Pan translates within the plane.
 * Zoom (+/−) steps the plane along its normal through the volume.
 */
(function (global) {
  "use strict";

  function clamp(v, lo, hi) {
    return Math.max(lo, Math.min(hi, v));
  }

  function trilinearSample(vol, d, x, y, z) {
    if (x < 0 || y < 0 || z < 0 || x > d - 1 || y > d - 1 || z > d - 1) {
      return 0;
    }
    var x0 = Math.floor(x);
    var y0 = Math.floor(y);
    var z0 = Math.floor(z);
    var x1 = Math.min(x0 + 1, d - 1);
    var y1 = Math.min(y0 + 1, d - 1);
    var z1 = Math.min(z0 + 1, d - 1);
    var xd = x - x0;
    var yd = y - y0;
    var zd = z - z0;
    var i000 = vol[x0 * d * d + y0 * d + z0];
    var i100 = vol[x1 * d * d + y0 * d + z0];
    var i010 = vol[x0 * d * d + y1 * d + z0];
    var i110 = vol[x1 * d * d + y1 * d + z0];
    var i001 = vol[x0 * d * d + y0 * d + z1];
    var i101 = vol[x1 * d * d + y0 * d + z1];
    var i011 = vol[x0 * d * d + y1 * d + z1];
    var i111 = vol[x1 * d * d + y1 * d + z1];
    var c00 = i000 * (1 - xd) + i100 * xd;
    var c01 = i001 * (1 - xd) + i101 * xd;
    var c10 = i010 * (1 - xd) + i110 * xd;
    var c11 = i011 * (1 - xd) + i111 * xd;
    var c0 = c00 * (1 - yd) + c10 * yd;
    var c1 = c01 * (1 - yd) + c11 * yd;
    return c0 * (1 - zd) + c1 * zd;
  }

  function rotationMatrix(yaw, pitch) {
    var cy = Math.cos(yaw);
    var sy = Math.sin(yaw);
    var cp = Math.cos(pitch);
    var sp = Math.sin(pitch);
    return [
      [cy, 0, sy],
      [sp * sy, cp, -sp * cy],
      [-cp * sy, sp, cp * cy],
    ];
  }

  function decodeFloat32Volume(b64, d) {
    var binary = atob(b64);
    var len = binary.length;
    var bytes = new Uint8Array(len);
    var chunk = 0x8000;
    for (var offset = 0; offset < len; offset += chunk) {
      var end = Math.min(offset + chunk, len);
      for (var i = offset; i < end; i++) {
        bytes[i] = binary.charCodeAt(i);
      }
    }
    return new Float32Array(bytes.buffer);
  }

  /** Match particle explorer image grid meta band (see _particle_explorer_scatterjs.html). */
  var VSLICE_META_TOP_FRAC = 0.19;
  var VSLICE_META_IMG_GAP = 1;
  var VSLICE_CELL_PAD_Y = 2;
  var VSLICE_LABEL_META_HEIGHT_FRAC = 0.7885;
  var VSLICE_GEORGIA_CAP_HEIGHT_EM = 0.715;
  var VSLICE_MONTAGE_PAPER = "#faf8f4";
  var VSLICE_MONTAGE_LABEL_COLOR = "#243b53";
  var VSLICE_DEPTH_STEP = 3;

  function VolumeSliceCanvas(opts) {
    this.canvas = opts.canvas;
    this.onViewChange = opts.onViewChange || null;
    this.ctx = this.canvas.getContext("2d");
    /** @type {Array<{id:string,d:number,vol:Float32Array,sliceData:Float32Array|null,sliceMin:number,sliceMax:number}>} */
    this.layers = [];
    this.vol = null;
    this.d = 0;
    this.sliceData = null;
    this.sliceMin = 0;
    this.sliceMax = 1;
    this.yaw = 0;
    this.pitch = 0;
    /** In-plane translation of the cutting plane (voxel units, plane-local axes). */
    this.slicePanU = 0;
    this.slicePanV = 0;
    /** Offset of the slice centre along the plane normal (voxel units). */
    this.sliceDepth = 0;
    this.dragMode = null;
    this.lastX = 0;
    this.lastY = 0;
    /** Display contrast slider level (0–100, 50 = neutral). */
    this.contrastLevel = 50;
    this.recomputePending = false;
    this._wirePointer();
  }

  VolumeSliceCanvas.prototype._resetSliceView = function () {
    this.yaw = 0;
    this.pitch = 0;
    this.slicePanU = 0;
    this.slicePanV = 0;
    this.sliceDepth = 0;
  };

  VolumeSliceCanvas.prototype._syncPrimaryLayer = function () {
    if (this.layers.length) {
      var primary = this.layers[0];
      this.vol = primary.vol;
      this.d = primary.d;
      this.sliceData = primary.sliceData;
      this.sliceMin = primary.sliceMin;
      this.sliceMax = primary.sliceMax;
    } else {
      this.vol = null;
      this.d = 0;
      this.sliceData = null;
      this.sliceMin = 0;
      this.sliceMax = 1;
    }
  };

  VolumeSliceCanvas.prototype._gridLayout = function (layerCount) {
    var n = Math.max(1, Math.ceil(Math.sqrt(Math.max(1, layerCount))));
    return { cols: n, rows: n, n: n, slots: n * n };
  };

  /** Height / width of one montage-style slice card (meta band + square image). */
  VolumeSliceCanvas.prototype._cardBlockHeight = function (cellW) {
    var metaH = cellW * VSLICE_META_TOP_FRAC;
    return metaH + VSLICE_META_IMG_GAP + cellW + VSLICE_CELL_PAD_Y;
  };

  VolumeSliceCanvas.prototype.cardAspectRatio = function () {
    return this._cardBlockHeight(1) / 1;
  };

  VolumeSliceCanvas.prototype._gridGap = function () {
    return 6;
  };

  VolumeSliceCanvas.prototype._metaRowExtra = function () {
    return VSLICE_CELL_PAD_Y + VSLICE_META_IMG_GAP;
  };

  /** n×n grid of slice cards; slots may exceed layerCount for empty cells. */
  VolumeSliceCanvas.prototype._computeGridCells = function (cw, ch, layerCount) {
    var grid = this._gridLayout(layerCount);
    var n = grid.n;
    var gap = this._gridGap();
    var cellW = (cw - gap * (n + 1)) / n;
    var rowExtra = this._metaRowExtra();
    var blockAspect = 1 + VSLICE_META_TOP_FRAC;
    var cellWFromH = (ch - gap * (n + 1) - n * rowExtra) / (n * blockAspect);
    cellW = Math.min(cellW, cellWFromH);
    if (!isFinite(cellW) || cellW <= 0) cellW = Math.max(1, Math.min(cw, ch) / n);

    var metaH = cellW * VSLICE_META_TOP_FRAC;
    var imgSize = cellW;
    var blockH = this._cardBlockHeight(cellW);
    var cells = [];
    for (var slot = 0; slot < grid.slots; slot++) {
      var col = slot % n;
      var row = Math.floor(slot / n);
      var x0 = gap + col * (cellW + gap);
      var y0 = gap + row * (blockH + gap);
      cells.push({
        slot: slot,
        layerIndex: slot < layerCount ? slot : -1,
        x0: x0,
        y0: y0,
        cellW: cellW,
        metaH: metaH,
        imgX: x0,
        imgY: y0 + metaH + VSLICE_META_IMG_GAP,
        imgSize: imgSize,
        blockH: blockH,
      });
    }
    return { grid: grid, gap: gap, cellW: cellW, cells: cells };
  };

  VolumeSliceCanvas.prototype._drawCellChrome = function (ctx, cell, layer) {
    ctx.fillStyle = layer && layer.cellBg ? layer.cellBg : VSLICE_MONTAGE_PAPER;
    ctx.fillRect(cell.x0, cell.y0, cell.cellW, cell.blockH);
    if (layer && layer.cellBorder) {
      ctx.strokeStyle = layer.cellBorder;
      ctx.lineWidth = 2;
      ctx.strokeRect(cell.x0 + 1, cell.y0 + 1, cell.cellW - 2, cell.blockH - 2);
    }
    ctx.fillStyle = "#e8ecf1";
    ctx.fillRect(cell.imgX, cell.imgY, cell.imgSize, cell.imgSize);
    if (layer) {
      this._drawMontageLabel(
        ctx,
        cell.x0,
        cell.y0,
        cell.cellW,
        cell.metaH,
        layer.label,
        layer.manual
      );
      if (layer.sliceData) {
        var off = this._drawLayerImage(layer);
        ctx.drawImage(off, cell.imgX, cell.imgY, cell.imgSize, cell.imgSize);
      }
    }
  };

  VolumeSliceCanvas.prototype._drawMontageLabel = function (ctx, x, y, cellW, metaH, label, manual) {
    var labStr = String(label || "");
    if (!labStr) return;
    var letterBandPx = metaH * VSLICE_LABEL_META_HEIGHT_FRAC;
    var letterFontPx = letterBandPx / VSLICE_GEORGIA_CAP_HEIGHT_EM;
    if (labStr.length > 1) {
      letterFontPx *= Math.min(0.92, 1.22 / labStr.length);
    }
    var fontFamily = manual
      ? 'Palatino, "Palatino Linotype", "Book Antiqua", serif'
      : 'Georgia, "Times New Roman", serif';
    var fontStyle = manual ? "italic " : "";
    ctx.font = fontStyle + "bold " + letterFontPx.toFixed(2) + "px " + fontFamily;
    ctx.fillStyle = VSLICE_MONTAGE_LABEL_COLOR;
    ctx.textBaseline = "middle";
    ctx.textAlign = "left";
    ctx.fillText(labStr, x + letterFontPx * 0.11, y + metaH * 0.5);
  };

  VolumeSliceCanvas.prototype._planeGeometry = function (layerD) {
    var d = Number(layerD) || (this.layers.length ? this.layers[0].d : 0);
    if (!d || !this.layers.length) return null;
    var center = (d - 1) / 2;
    var R = rotationMatrix(this.yaw, this.pitch);
    var t0x = R[0][0];
    var t0y = R[1][0];
    var t0z = R[2][0];
    var t1x = R[0][1];
    var t1y = R[1][1];
    var t1z = R[2][1];
    var offX = this.slicePanU * t0x + this.slicePanV * t1x;
    var offY = this.slicePanU * t0y + this.slicePanV * t1y;
    var offZ = this.slicePanU * t0z + this.slicePanV * t1z;
    var nx = R[0][2];
    var ny = R[1][2];
    var nz = R[2][2];
    var nlen = Math.sqrt(nx * nx + ny * ny + nz * nz) || 1;
    nx /= nlen;
    ny /= nlen;
    nz /= nlen;
    var cx = center + offX + this.sliceDepth * nx;
    var cy = center + offY + this.sliceDepth * ny;
    var cz = center + offZ + this.sliceDepth * nz;
    var ix = Math.round(cx);
    var iy = Math.round(cy);
    var iz = Math.round(cz);
    var ax = Math.abs(nx);
    var ay = Math.abs(ny);
    var az = Math.abs(nz);
    var orth = null;
    var orthThreshold = 0.92;
    if (ax >= orthThreshold && ax >= ay && ax >= az) {
      orth = { label: "Y–Z (fixed X)", fixedAxis: "X", fixedIndex: ix };
    } else if (ay >= orthThreshold && ay >= ax && ay >= az) {
      orth = { label: "X–Z (fixed Y)", fixedAxis: "Y", fixedIndex: iy };
    } else if (az >= orthThreshold && az >= ax && az >= ay) {
      orth = { label: "X–Y (fixed Z)", fixedAxis: "Z", fixedIndex: iz };
    }
    return {
      d: d,
      center: center,
      cx: cx,
      cy: cy,
      cz: cz,
      ix: ix,
      iy: iy,
      iz: iz,
      nx: nx,
      ny: ny,
      nz: nz,
      span: Math.max(1, Math.round(2 * center)),
      orthogonal: orth,
      R: R,
      offX: offX,
      offY: offY,
      offZ: offZ,
    };
  };

  VolumeSliceCanvas.prototype._clampSliceDepth = function () {
    if (!this.layers.length) return;
    var d = this.layers[0].d;
    var center = (d - 1) / 2;
    var R = rotationMatrix(this.yaw, this.pitch);
    var offX = this.slicePanU * R[0][0] + this.slicePanV * R[0][1];
    var offY = this.slicePanU * R[1][0] + this.slicePanV * R[1][1];
    var offZ = this.slicePanU * R[2][0] + this.slicePanV * R[2][1];
    var nx = R[0][2];
    var ny = R[1][2];
    var nz = R[2][2];
    var lo = -Infinity;
    var hi = Infinity;
    function axisInterval(off, n) {
      if (Math.abs(n) < 1e-8) return;
      var dLo = (0 - center - off) / n;
      var dHi = ((d - 1) - center - off) / n;
      lo = Math.max(lo, Math.min(dLo, dHi));
      hi = Math.min(hi, Math.max(dLo, dHi));
    }
    axisInterval(offX, nx);
    axisInterval(offY, ny);
    axisInterval(offZ, nz);
    if (!isFinite(lo) || !isFinite(hi) || lo > hi) {
      this.sliceDepth = 0;
      return;
    }
    this.sliceDepth = clamp(this.sliceDepth, lo, hi);
  };

  VolumeSliceCanvas.prototype.getViewInfo = function () {
    var geom = this._planeGeometry();
    if (!geom) {
      return { d: 0 };
    }
    return {
      d: geom.d,
      ix: geom.ix,
      iy: geom.iy,
      iz: geom.iz,
      centre: [geom.ix, geom.iy, geom.iz],
      normal: [
        Math.round(geom.nx * 1000) / 1000,
        Math.round(geom.ny * 1000) / 1000,
        Math.round(geom.nz * 1000) / 1000,
      ],
      span: geom.span,
      orthogonal: geom.orthogonal,
    };
  };

  VolumeSliceCanvas.prototype._notifyViewChange = function () {
    if (typeof this.onViewChange === "function") {
      this.onViewChange(this.getViewInfo());
    }
  };

  /** Pan the cutting plane within the volume (plane-local u/v in voxels). */
  VolumeSliceCanvas.prototype.panBy = function (du, dv) {
    if (!this.layers.length) return;
    this.slicePanU += Number(du) || 0;
    this.slicePanV += Number(dv) || 0;
    this.scheduleRecompute();
  };

  /** Step the cutting plane along its normal (+factor moves forward). */
  VolumeSliceCanvas.prototype.zoomBy = function (factor) {
    if (!this.layers.length) return;
    var f = Number(factor);
    if (!isFinite(f) || f === 0) return;
    var sign = f >= 1 ? 1 : -1;
    this.sliceDepth += sign * VSLICE_DEPTH_STEP;
    this._clampSliceDepth();
    this.scheduleRecompute();
  };

  VolumeSliceCanvas.prototype.zoomIn = function () {
    this.zoomBy(1.12);
  };

  VolumeSliceCanvas.prototype.zoomOut = function () {
    this.zoomBy(1 / 1.12);
  };

  VolumeSliceCanvas.prototype.loadVolume = function (b64, d) {
    this.setVolumes([{ b64: b64, d: d, id: "" }], true);
  };

  VolumeSliceCanvas.prototype.setVolumes = function (volumeSpecs, resetView) {
    var specs = volumeSpecs || [];
    var newLayers = [];
    for (var i = 0; i < specs.length; i++) {
      var spec = specs[i];
      var d = Number(spec.d);
      if (!spec.b64 || !isFinite(d) || d < 1) continue;
      newLayers.push({
        id: String(spec.id || ""),
        label: String(spec.label || ""),
        manual: !!spec.manual,
        cellBg: spec.cellBg || null,
        cellBorder: spec.cellBorder || null,
        d: d,
        vol: decodeFloat32Volume(spec.b64, d),
        sliceData: null,
        sliceMin: 0,
        sliceMax: 1,
      });
    }
    this.layers = newLayers;
    if (resetView !== false) this._resetSliceView();
    var self = this;
    window.requestAnimationFrame(function () {
      self.recomputeSlice();
      self.draw();
      self._notifyViewChange();
    });
  };

  VolumeSliceCanvas.prototype._contrastFactor = function () {
    return Math.pow(2, (this.contrastLevel - 50) / 20);
  };

  VolumeSliceCanvas.prototype.setContrastLevel = function (level) {
    var n = Number(level);
    if (!isFinite(n)) n = 50;
    this.contrastLevel = clamp(n, 0, 100);
    this.draw();
  };

  VolumeSliceCanvas.prototype.resetView = function () {
    this._resetSliceView();
    this.recomputeSlice();
    this.draw();
    this._notifyViewChange();
  };

  VolumeSliceCanvas.prototype.recomputeSlice = function () {
    if (!this.layers.length) return;
    for (var li = 0; li < this.layers.length; li++) {
      var layer = this.layers[li];
      var geom = this._planeGeometry(layer.d);
      if (!geom) continue;
      var d = geom.d;
      var center = geom.center;
      var scale = center;
      var R = geom.R;
      var cx = geom.cx;
      var cy = geom.cy;
      var cz = geom.cz;
      if (!layer.sliceData || layer.sliceData.length !== d * d) {
        layer.sliceData = new Float32Array(d * d);
      }
      var minV = Infinity;
      var maxV = -Infinity;
      for (var j = 0; j < d; j++) {
        for (var i = 0; i < d; i++) {
          var u = (i - center) / scale;
          var v = (j - center) / scale;
          var px = cx + scale * (R[0][0] * u + R[0][1] * v);
          var py = cy + scale * (R[1][0] * u + R[1][1] * v);
          var pz = cz + scale * (R[2][0] * u + R[2][1] * v);
          var val = trilinearSample(layer.vol, d, px, py, pz);
          var idx = j * d + i;
          layer.sliceData[idx] = val;
          if (val < minV) minV = val;
          if (val > maxV) maxV = val;
        }
      }
      if (maxV <= minV) {
        minV = 0;
        maxV = 1;
      }
      layer.sliceMin = minV;
      layer.sliceMax = maxV;
    }
    this._syncPrimaryLayer();
  };

  VolumeSliceCanvas.prototype._drawLayerImage = function (layer) {
    var d = layer.d;
    var img = this.ctx.createImageData(d, d);
    var span = layer.sliceMax - layer.sliceMin;
    var data = img.data;
    var factor = this._contrastFactor();
    var mid = 0.5;
    for (var k = 0; k < d * d; k++) {
      var norm = span > 0 ? (layer.sliceData[k] - layer.sliceMin) / span : 0;
      norm = clamp(mid + (norm - mid) * factor, 0, 1);
      var g = Math.round(norm * 255);
      var p = k * 4;
      data[p] = g;
      data[p + 1] = g;
      data[p + 2] = g;
      data[p + 3] = 255;
    }
    var off = document.createElement("canvas");
    off.width = d;
    off.height = d;
    off.getContext("2d").putImageData(img, 0, 0);
    return off;
  };

  VolumeSliceCanvas.prototype.draw = function () {
    if (!this.layers.length) return;
    var canvas = this.canvas;
    var ctx = this.ctx;
    var cw = canvas.clientWidth || 256;
    var ch = canvas.clientHeight || 256;
    if (canvas.width !== cw || canvas.height !== ch) {
      canvas.width = cw;
      canvas.height = ch;
    }
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, cw, ch);

    var layerCount = this.layers.length;
    var layout = this._computeGridCells(cw, ch, layerCount);

    for (var si = 0; si < layout.cells.length; si++) {
      var cell = layout.cells[si];
      var layer = cell.layerIndex >= 0 ? this.layers[cell.layerIndex] : null;
      this._drawCellChrome(ctx, cell, layer);
    }
    this._syncPrimaryLayer();
  };

  VolumeSliceCanvas.prototype.scheduleRecompute = function () {
    var self = this;
    if (this.recomputePending) return;
    this.recomputePending = true;
    window.requestAnimationFrame(function () {
      self.recomputePending = false;
      if (!self.layers.length) return;
      self.recomputeSlice();
      self.draw();
      self._notifyViewChange();
    });
  };

  VolumeSliceCanvas.prototype._wirePointer = function () {
    var self = this;
    var canvas = this.canvas;

    canvas.addEventListener("contextmenu", function (e) {
      e.preventDefault();
    });

    canvas.addEventListener("pointerdown", function (e) {
      if (!self.layers.length) return;
      if (e.button !== 0) return;
      canvas.setPointerCapture(e.pointerId);
      self.lastX = e.clientX;
      self.lastY = e.clientY;
      self.dragMode = "rotate";
      e.preventDefault();
    });

    canvas.addEventListener("pointermove", function (e) {
      if (self.dragMode !== "rotate" || !self.layers.length) return;
      var dx = e.clientX - self.lastX;
      var dy = e.clientY - self.lastY;
      self.lastX = e.clientX;
      self.lastY = e.clientY;
      self.yaw += dx * 0.012;
      self.pitch = clamp(self.pitch + dy * 0.012, -1.45, 1.45);
      self.scheduleRecompute();
    });

    function endDrag(e) {
      if (self.dragMode) {
        try {
          canvas.releasePointerCapture(e.pointerId);
        } catch (err) { /* already released */ }
      }
      self.dragMode = null;
    }
    canvas.addEventListener("pointerup", endDrag);
    canvas.addEventListener("pointercancel", endDrag);

    canvas.addEventListener("dblclick", function () {
      if (!self.layers.length) return;
      self.resetView();
    });
  };

  global.CryoVolumeSliceCanvas = VolumeSliceCanvas;
})(typeof window !== "undefined" ? window : this);
