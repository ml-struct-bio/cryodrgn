/**
 * Interactive single-slice viewer for decoded cryoDRGN volumes.
 *
 * Default: central XY slice (fixed Z). Left-drag rotates the cutting plane,
 * right-drag pans, wheel zooms.
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
    this.hintEl = opts.hintEl || null;
    this.ctx = this.canvas.getContext("2d");
    this.vol = null;
    this.d = 0;
    this.sliceData = null;
    this.sliceMin = 0;
    this.sliceMax = 1;
    this.yaw = 0;
    this.pitch = 0;
    this.viewScale = 1;
    this.viewPanX = 0;
    this.viewPanY = 0;
    this.dragMode = null;
    this.lastX = 0;
    this.lastY = 0;
    this.rafPending = false;
    this.recomputePending = false;
    this._wirePointer();
  }

  VolumeSliceCanvas.prototype.setHint = function (text) {
    if (this.hintEl) this.hintEl.textContent = text || "";
  };

  VolumeSliceCanvas.prototype.loadVolume = function (b64, d) {
    this.d = d;
    this.vol = decodeFloat32Volume(b64, d);
    this.yaw = 0;
    this.pitch = 0;
    this.viewScale = 1;
    this.viewPanX = 0;
    this.viewPanY = 0;
    this.setHint("Left-drag: rotate · Right-drag: pan · Wheel: zoom");
    var self = this;
    window.requestAnimationFrame(function () {
      self.recomputeSlice();
      self.draw();
    });
  };

  VolumeSliceCanvas.prototype.resetView = function () {
    this.yaw = 0;
    this.pitch = 0;
    this.viewScale = 1;
    this.viewPanX = 0;
    this.viewPanY = 0;
    this.recomputeSlice();
    this.draw();
  };

  VolumeSliceCanvas.prototype.recomputeSlice = function () {
    if (!this.vol || !this.d) return;
    var d = this.d;
    var center = (d - 1) / 2;
    var scale = center;
    var R = rotationMatrix(this.yaw, this.pitch);
    if (!this.sliceData || this.sliceData.length !== d * d) {
      this.sliceData = new Float32Array(d * d);
    }
    var minV = Infinity;
    var maxV = -Infinity;
    for (var j = 0; j < d; j++) {
      for (var i = 0; i < d; i++) {
        var u = (i - center) / scale;
        var v = (j - center) / scale;
        var px = center + scale * (R[0][0] * u + R[0][1] * v);
        var py = center + scale * (R[1][0] * u + R[1][1] * v);
        var pz = center + scale * (R[2][0] * u + R[2][1] * v);
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
    for (var k = 0; k < d * d; k++) {
      var g = Math.round(((this.sliceData[k] - this.sliceMin) / span) * 255);
      g = clamp(g, 0, 255);
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
    ctx.fillStyle = "#e8ecf1";
    ctx.fillRect(0, 0, cw, ch);
    ctx.save();
    ctx.translate(cw / 2 + this.viewPanX, ch / 2 + this.viewPanY);
    ctx.scale(this.viewScale, this.viewScale);
    var drawW = Math.min(cw, ch) * 0.92;
    ctx.drawImage(off, -drawW / 2, -drawW / 2, drawW, drawW);
    ctx.restore();
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
    });
  };

  VolumeSliceCanvas.prototype.scheduleDraw = function () {
    var self = this;
    if (this.rafPending) return;
    this.rafPending = true;
    requestAnimationFrame(function () {
      self.rafPending = false;
      self.draw();
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
      canvas.setPointerCapture(e.pointerId);
      self.lastX = e.clientX;
      self.lastY = e.clientY;
      if (e.button === 0) {
        self.dragMode = "rotate";
      } else if (e.button === 1 || e.button === 2) {
        self.dragMode = "pan";
      }
      e.preventDefault();
    });

    canvas.addEventListener("pointermove", function (e) {
      if (!self.dragMode || !self.vol) return;
      var dx = e.clientX - self.lastX;
      var dy = e.clientY - self.lastY;
      self.lastX = e.clientX;
      self.lastY = e.clientY;
      if (self.dragMode === "rotate") {
        self.yaw += dx * 0.012;
        self.pitch = clamp(self.pitch + dy * 0.012, -1.45, 1.45);
        self.scheduleRecompute();
      } else if (self.dragMode === "pan") {
        self.viewPanX += dx;
        self.viewPanY += dy;
        self.scheduleDraw();
      }
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

    canvas.addEventListener(
      "wheel",
      function (e) {
        if (!self.vol) return;
        e.preventDefault();
        var factor = e.deltaY < 0 ? 1.08 : 1 / 1.08;
        self.viewScale = clamp(self.viewScale * factor, 0.35, 12);
        self.scheduleDraw();
      },
      { passive: false }
    );

    canvas.addEventListener("dblclick", function () {
      if (!self.vol) return;
      self.resetView();
    });
  };

  global.CryoVolumeSliceCanvas = VolumeSliceCanvas;
})(typeof window !== "undefined" ? window : this);
