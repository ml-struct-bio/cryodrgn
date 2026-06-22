/**
 * Interactive single-slice viewer for decoded cryoDRGN volumes.
 *
 * Left-drag rotates the cutting plane. Pan/zoom buttons translate and magnify
 * the slice within the volume cube (not the 2-D canvas viewport).
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

  /** 3×3 rotation: Ry(yaw) then Rx(pitch). */
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

  function VolumeSliceCanvas(opts) {
    this.canvas = opts.canvas;
    this.onViewChange = opts.onViewChange || null;
    this.ctx = this.canvas.getContext("2d");
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
    /** Magnification within the cutting plane (>1 zooms into the volume). */
    this.sliceZoom = 1;
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
    this.sliceZoom = 1;
  };

  VolumeSliceCanvas.prototype._planeGeometry = function () {
    if (!this.vol || !this.d) return null;
    var d = this.d;
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
    var cx = center + offX;
    var cy = center + offY;
    var cz = center + offZ;
    var nx = R[0][2];
    var ny = R[1][2];
    var nz = R[2][2];
    var nlen = Math.sqrt(nx * nx + ny * ny + nz * nz) || 1;
    nx /= nlen;
    ny /= nlen;
    nz /= nlen;
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
      span: Math.max(1, Math.round((2 * center) / this.sliceZoom)),
      orthogonal: orth,
      R: R,
      invZoom: 1 / this.sliceZoom,
      offX: offX,
      offY: offY,
      offZ: offZ,
    };
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
    if (!this.vol) return;
    this.slicePanU += Number(du) || 0;
    this.slicePanV += Number(dv) || 0;
    this.scheduleRecompute();
  };

  /** Zoom the cutting plane within the volume. */
  VolumeSliceCanvas.prototype.zoomBy = function (factor) {
    if (!this.vol) return;
    this.sliceZoom = clamp(this.sliceZoom * factor, 0.35, 12);
    this.scheduleRecompute();
  };

  VolumeSliceCanvas.prototype.zoomIn = function () {
    this.zoomBy(1.12);
  };

  VolumeSliceCanvas.prototype.zoomOut = function () {
    this.zoomBy(1 / 1.12);
  };

  VolumeSliceCanvas.prototype.loadVolume = function (b64, d) {
    this.d = d;
    this.vol = decodeFloat32Volume(b64, d);
    this._resetSliceView();
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
    if (!this.vol || !this.d) return;
    var geom = this._planeGeometry();
    if (!geom) return;
    var d = geom.d;
    var center = geom.center;
    var scale = center;
    var R = geom.R;
    var invZoom = geom.invZoom;
    var offX = geom.offX;
    var offY = geom.offY;
    var offZ = geom.offZ;
    if (!this.sliceData || this.sliceData.length !== d * d) {
      this.sliceData = new Float32Array(d * d);
    }
    var minV = Infinity;
    var maxV = -Infinity;
    for (var j = 0; j < d; j++) {
      for (var i = 0; i < d; i++) {
        var u = ((i - center) / scale) * invZoom;
        var v = ((j - center) / scale) * invZoom;
        var px = center + offX + scale * (R[0][0] * u + R[0][1] * v);
        var py = center + offY + scale * (R[1][0] * u + R[1][1] * v);
        var pz = center + offZ + scale * (R[2][0] * u + R[2][1] * v);
        var val = trilinearSample(this.vol, d, px, py, pz);
        var idx = j * d + i;
        this.sliceData[idx] = val;
        if (val < minV) minV = val;
        if (val > maxV) maxV = val;
      }
    }
    if (maxV <= minV) {
      minV = 0;
      maxV = 1;
    }
    this.sliceMin = minV;
    this.sliceMax = maxV;
  };

  VolumeSliceCanvas.prototype.draw = function () {
    if (!this.sliceData || !this.d) return;
    var d = this.d;
    var canvas = this.canvas;
    var ctx = this.ctx;
    var cw = canvas.clientWidth || d;
    var ch = canvas.clientHeight || d;
    if (canvas.width !== cw || canvas.height !== ch) {
      canvas.width = cw;
      canvas.height = ch;
    }
    var img = ctx.createImageData(d, d);
    var span = this.sliceMax - this.sliceMin;
    var data = img.data;
    var factor = this._contrastFactor();
    var mid = 0.5;
    for (var k = 0; k < d * d; k++) {
      var norm = span > 0 ? (this.sliceData[k] - this.sliceMin) / span : 0;
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

    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.fillStyle = "#243b53";
    ctx.fillRect(0, 0, cw, ch);
    var drawW = Math.min(cw, ch) * 0.92;
    ctx.drawImage(off, (cw - drawW) / 2, (ch - drawW) / 2, drawW, drawW);
  };

  VolumeSliceCanvas.prototype.scheduleRecompute = function () {
    var self = this;
    if (this.recomputePending) return;
    this.recomputePending = true;
    window.requestAnimationFrame(function () {
      self.recomputePending = false;
      if (!self.vol) return;
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
      if (!self.vol) return;
      if (e.button !== 0) return;
      canvas.setPointerCapture(e.pointerId);
      self.lastX = e.clientX;
      self.lastY = e.clientY;
      self.dragMode = "rotate";
      e.preventDefault();
    });

    canvas.addEventListener("pointermove", function (e) {
      if (self.dragMode !== "rotate" || !self.vol) return;
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
      if (!self.vol) return;
      self.resetView();
    });
  };

  global.CryoVolumeSliceCanvas = VolumeSliceCanvas;
})(typeof window !== "undefined" ? window : this);
