/**
 * vtk.js volume raycasting viewer for cryoDRGN dashboard (single volume).
 */
import "@kitware/vtk.js/Rendering/Profiles/Volume";
import vtkColorTransferFunction from "@kitware/vtk.js/Rendering/Core/ColorTransferFunction";
import vtkDataArray from "@kitware/vtk.js/Common/Core/DataArray";
import vtkImageData from "@kitware/vtk.js/Common/DataModel/ImageData";
import vtkPiecewiseFunction from "@kitware/vtk.js/Common/DataModel/PiecewiseFunction";
import vtkGenericRenderWindow from "@kitware/vtk.js/Rendering/Misc/GenericRenderWindow";
import vtkVolume from "@kitware/vtk.js/Rendering/Core/Volume";
import vtkVolumeMapper from "@kitware/vtk.js/Rendering/Core/VolumeMapper";
import { cross } from "@kitware/vtk.js/Common/Core/Math";

var VOLUME_BG = [0.98, 0.97, 0.96];

function utils() {
  if (typeof window !== "undefined" && window.CryoVolume3dUtils) {
    return window.CryoVolume3dUtils;
  }
  throw new Error("CryoVolume3dUtils is not loaded.");
}

function clamp(v, lo, hi) {
  return Math.max(lo, Math.min(hi, v));
}

function volumeRange(values) {
  var min = Infinity;
  var max = -Infinity;
  for (var i = 0; i < values.length; i++) {
    var v = values[i];
    if (!isFinite(v)) continue;
    if (v < min) min = v;
    if (v > max) max = v;
  }
  if (!isFinite(min) || !isFinite(max) || min === max) {
    return { min: 0, max: 1 };
  }
  return { min: min, max: max };
}

function applyTransferFunctions(ctfun, ofun, range, isoLevel) {
  var threshold = clamp(Number(isoLevel), range.min, range.max);
  var span = range.max - range.min;
  if (!(span > 0)) span = 1;
  var minSep = span * 0.0015;
  var rampLo = threshold - span * 0.018;
  var rampHi = threshold + span * 0.045;

  rampLo = Math.max(range.min, Math.min(rampLo, range.max - 3 * minSep));
  threshold = Math.max(rampLo + minSep, Math.min(threshold, range.max - 2 * minSep));
  rampHi = Math.max(threshold + minSep, Math.min(rampHi, range.max - minSep));

  ctfun.removeAllPoints();
  ctfun.addRGBPoint(range.min, 0.04, 0.04, 0.05);
  ctfun.addRGBPoint(rampLo, 0.14, 0.14, 0.16);
  ctfun.addRGBPoint(threshold, 0.58, 0.58, 0.61);
  ctfun.addRGBPoint(rampHi, 0.94, 0.94, 0.96);
  ctfun.addRGBPoint(range.max, 1.0, 1.0, 1.0);

  ofun.removeAllPoints();
  ofun.addPoint(range.min, 0.0);
  ofun.addPoint(rampLo, 0.0);
  ofun.addPoint(threshold, 0.32);
  ofun.addPoint(rampHi, 0.50);
  ofun.addPoint(range.max, 0.55);
}

export class VolumeRaycastView {
  constructor(container) {
    this.container = container;
    this.grw = vtkGenericRenderWindow.newInstance({
      background: VOLUME_BG,
      listenWindowResize: false,
    });
    this.grw.setContainer(container);
    this.renderer = this.grw.getRenderer();
    this.renderWindow = this.grw.getRenderWindow();
    this.volume = null;
    this.mapper = null;
    this.imageData = null;
    this.ctfun = vtkColorTransferFunction.newInstance();
    this.ofun = vtkPiecewiseFunction.newInstance();
    this._resizeObserver = null;
    this.volumeRange = null;
    this.isoSliderRange = null;
    this.percentileSamples = null;
    this.isoLevel = 0;
    this.onCameraChanged = null;
    this._baseAzimuth = 25;
    this._baseElevation = 12;
    /** Camera rotation captured right after ``resetCamera()`` (ChimeraX ``view orient``). */
    this._referenceRot3 = null;
    /** Live camera rotation relative to the reference; null means "at default view". */
    this._chimeraxDelta3 = null;
    this._cameraChangeFromUser = false;
    this._wireCameraChangeHook();
  }

  _wireCameraChangeHook() {
    if (!this.renderWindow) return;
    var interactor = this.renderWindow.getInteractor();
    if (!interactor) return;
    var self = this;
    if (typeof interactor.onEndAnimation === "function") {
      interactor.onEndAnimation(function () {
        if (!self._cameraChangeFromUser) return;
        self._cameraChangeFromUser = false;
        if (typeof self.onCameraChanged === "function") self.onCameraChanged();
      });
    }
    if (this.container) {
      this.container.addEventListener("pointerdown", function (e) {
        if (e.button === 0) self._cameraChangeFromUser = true;
      });
    }
  }

  _reapplyTransferFunctions() {
    if (!this.volume || !this.volumeRange || !this.ctfun || !this.ofun) {
      return;
    }
    applyTransferFunctions(
      this.ctfun,
      this.ofun,
      this.volumeRange,
      this.isoLevel
    );
    if (this.renderWindow) this.renderWindow.render();
  }

  getIsoLevel() {
    return this.isoLevel;
  }

  getIsoDataRange() {
    if (this.isoSliderRange) return this.isoSliderRange;
    if (this.volumeRange) return { min: this.volumeRange.min, max: this.volumeRange.max };
    return { min: 0, max: 1 };
  }

  setIsoLevel(level) {
    var range = this.getIsoDataRange();
    if (!this.volumeRange) {
      this.isoLevel = Number(level);
      return;
    }
    this.isoLevel = clamp(Number(level), range.min, range.max);
    this._reapplyTransferFunctions();
  }

  _clearVolume() {
    if (this.volume) {
      this.renderer.removeVolume(this.volume);
      this.volume.delete();
      this.volume = null;
    }
    if (this.mapper) {
      this.mapper.delete();
      this.mapper = null;
    }
    if (this.imageData) {
      this.imageData.delete();
      this.imageData = null;
    }
  }

  setVolumeFromB64(b64, d, options) {
    options = options || {};
    var U = utils();
    var srcDim = Number(d);
    if (!b64 || !isFinite(srcDim) || srcDim < 1) {
      this._clearVolume();
      this.renderWindow.render();
      return;
    }
    var vol = U.decodeFloat32Volume(b64, srcDim);
    if (vol.length !== srcDim * srcDim * srcDim) {
      throw new Error("Volume byte length does not match D³.");
    }
    var values = U.downsampleVolumeBoxAverage(vol, srcDim, U.PLOT3D_TARGET_D);
    var dim = U.PLOT3D_TARGET_D;
    this._clearVolume();

    this.imageData = vtkImageData.newInstance();
    this.imageData.setDimensions(dim, dim, dim);
    this.imageData.setSpacing(1, 1, 1);
    this.imageData.setOrigin(0, 0, 0);

    var scalars = vtkDataArray.newInstance({
      name: "scalars",
      numberOfComponents: 1,
      values: values,
    });
    this.imageData.getPointData().setScalars(scalars);

    this.mapper = vtkVolumeMapper.newInstance();
    this.mapper.setInputData(this.imageData);
    this.mapper.setSampleDistance(1.0);

    this.volumeRange = volumeRange(values);
    this.percentileSamples = U.volumePercentileSamples(values);
    this.isoSliderRange = U.isoSliderDataRange(this.percentileSamples);
    this.isoLevel = 0.5 * (this.isoSliderRange.min + this.isoSliderRange.max);
    applyTransferFunctions(
      this.ctfun,
      this.ofun,
      this.volumeRange,
      this.isoLevel
    );

    this.volume = vtkVolume.newInstance();
    this.volume.setMapper(this.mapper);
    var prop = this.volume.getProperty();
    prop.setRGBTransferFunction(0, this.ctfun);
    prop.setScalarOpacity(0, this.ofun);
    prop.setScalarOpacityUnitDistance(0, 2.0);
    prop.setInterpolationTypeToLinear();
    prop.setShade(true);
    prop.setAmbient(0.25);
    prop.setDiffuse(0.75);
    prop.setSpecular(0.15);

    this.renderer.addVolume(this.volume);
    if (!options.skipDefaultCamera) {
      this._applyDefaultCamera();
    }
    this.renderer.resetCameraClippingRange();
    this.resize();
    this.renderWindow.render();
  }

  resize() {
    if (this.grw) this.grw.resize();
  }

  _normDeg(d) {
    d = ((d % 360) + 360) % 360;
    if (d > 180) d -= 360;
    return d;
  }

  _rot3FromCamera(cam) {
    if (!cam || !cam.getViewMatrix) return null;
    var mat = new Float64Array(16);
    cam.getViewMatrix(mat);
    var rot = new Array(9);
    for (var row = 0; row < 3; row++) {
      for (var col = 0; col < 3; col++) {
        rot[row * 3 + col] = mat[col * 4 + row];
      }
    }
    return rot;
  }

  _mat3Transpose(a) {
    return [
      a[0], a[3], a[6],
      a[1], a[4], a[7],
      a[2], a[5], a[8],
    ];
  }

  _mat3Mul(a, b) {
    var out = new Array(9);
    for (var r = 0; r < 3; r++) {
      for (var c = 0; c < 3; c++) {
        out[r * 3 + c] = a[r * 3] * b[c]
          + a[r * 3 + 1] * b[3 + c]
          + a[r * 3 + 2] * b[6 + c];
      }
    }
    return out;
  }

  /** Exact axis-angle (deg + unit axis) of a 3x3 rotation; null if ~identity. */
  _axisAngleFromMat3(r) {
    var trace = r[0] + r[4] + r[8];
    var cosA = clamp((trace - 1) / 2, -1, 1);
    var angle = Math.acos(cosA) * (180 / Math.PI);
    if (!isFinite(angle) || Math.abs(angle) < 1e-3) return null;
    var ax = r[7] - r[5];
    var ay = r[2] - r[6];
    var az = r[3] - r[1];
    var n = Math.sqrt(ax * ax + ay * ay + az * az);
    if (n < 1e-9) {
      // Near 180°: recover axis from the largest diagonal term.
      var k = 0;
      if (r[4] > r[0]) k = 1;
      if (r[8] > r[k * 3 + k]) k = 2;
      var v = [0, 0, 0];
      v[k] = Math.sqrt(Math.max(0, (r[k * 3 + k] - cosA) / (1 - cosA)));
      ax = v[0]; ay = v[1]; az = v[2];
      n = Math.sqrt(ax * ax + ay * ay + az * az) || 1;
    }
    return { axis: [ax / n, ay / n, az / n], degrees: angle };
  }

  /** Read the live VTK camera rotation relative to the default ("orient") view. */
  syncChimeraxTurnsFromCamera() {
    if (!this.renderer || !this._referenceRot3) return;
    var cam = this.renderer.getActiveCamera();
    if (!cam) return;
    var cur = this._rot3FromCamera(cam);
    if (!cur) return;
    this._chimeraxDelta3 = this._mat3Mul(this._mat3Transpose(this._referenceRot3), cur);
  }

  _applyDefaultCamera() {
    if (!this.renderer) return;
    this.renderer.resetCamera();
    this._chimeraxDelta3 = null;
    var cam = this.renderer.getActiveCamera();
    if (cam) {
      this._referenceRot3 = this._rot3FromCamera(cam);
    }
    if (this._baseAzimuth) cam.azimuth(this._baseAzimuth);
    if (this._baseElevation) cam.elevation(this._baseElevation);
    this.syncChimeraxTurnsFromCamera();
    this.renderer.resetCameraClippingRange();
  }

  resetCamera() {
    if (!this.volume || !this.renderer) return;
    this._applyDefaultCamera();
    this.renderWindow.render();
  }

  /** Map orient-relative rotation to ChimeraX ``turn y`` then ``turn x`` degrees. */
  _chimeraxYxTurnsFromDelta(deltaRot3) {
    var y = Math.atan2(deltaRot3[2], deltaRot3[0]) * (180 / Math.PI);
    var x = -Math.asin(clamp(deltaRot3[1], -1, 1)) * (180 / Math.PI);
    return { y: this._normDeg(y), x: this._normDeg(x) };
  }

  /** ChimeraX ``view #1 orient`` + ``turn y`` / ``turn x`` matching the live VTK camera. */
  getChimeraxViewTurns() {
    this.syncChimeraxTurnsFromCamera();
    if (!this._chimeraxDelta3) return [];
    var turns = this._chimeraxYxTurnsFromDelta(this._chimeraxDelta3);
    var out = [];
    if (Math.abs(turns.y) > 1e-4) out.push({ axis: "y", degrees: turns.y });
    if (Math.abs(turns.x) > 1e-4) out.push({ axis: "x", degrees: turns.x });
    return out;
  }

  applyChimeraxViewTurns(turns) {
    if (!this.renderer || !turns || !turns.length) return;
    this.renderer.resetCamera();
    var cam = this.renderer.getActiveCamera();
    if (cam) {
      this._referenceRot3 = this._rot3FromCamera(cam);
    }
    this._chimeraxDelta3 = null;
    for (var i = 0; i < turns.length; i++) {
      var t = turns[i] || {};
      var axis = String(t.axis || "").toLowerCase();
      var deg = Number(t.degrees);
      if (!isFinite(deg)) continue;
      if (axis === "y") {
        this.orbitBy(deg, 0);
      } else if (axis === "x") {
        this.orbitBy(0, deg);
      } else if (axis === "z" && typeof this.renderer.getActiveCamera().roll === "function") {
        this.renderer.getActiveCamera().roll(deg);
        this.renderer.resetCameraClippingRange();
        this.renderWindow.render();
      }
    }
    this.syncChimeraxTurnsFromCamera();
  }

  orbitBy(azimuthDeg, elevationDeg) {
    if (!this.renderer || !this.volume) return;
    const cam = this.renderer.getActiveCamera();
    if (azimuthDeg) {
      cam.azimuth(azimuthDeg);
    }
    if (elevationDeg) {
      cam.elevation(elevationDeg);
    }
    this.renderer.resetCameraClippingRange();
    this.renderWindow.render();
    this.syncChimeraxTurnsFromCamera();
  }

  /** Translate camera and focal point in the view plane (screen-space pan). */
  panBy(dx, dy) {
    if (!this.renderer || !this.volume) return;
    const cam = this.renderer.getActiveCamera();
    if (typeof cam.orthogonalizeViewUp === "function") cam.orthogonalizeViewUp();
    const camPos = cam.getPosition();
    const fp = cam.getFocalPoint();
    const up = cam.getViewUp();
    const vpn = cam.getViewPlaneNormal();
    const right = [0, 0, 0];
    cross(vpn, up, right);
    const dist = typeof cam.getDistance === "function" ? cam.getDistance() : 1;
    const scale = Math.max(dist * 0.025, 0.05);
    const stepX = Number(dx) * scale;
    const stepY = Number(dy) * scale;
    const nx = camPos[0] + right[0] * stepX + up[0] * stepY;
    const ny = camPos[1] + right[1] * stepX + up[1] * stepY;
    const nz = camPos[2] + right[2] * stepX + up[2] * stepY;
    const fx = fp[0] + right[0] * stepX + up[0] * stepY;
    const fy = fp[1] + right[1] * stepX + up[1] * stepY;
    const fz = fp[2] + right[2] * stepX + up[2] * stepY;
    cam.setPosition(nx, ny, nz);
    cam.setFocalPoint(fx, fy, fz);
    this.renderer.resetCameraClippingRange();
    this.renderWindow.render();
  }

  dollyBy(factor) {
    if (!this.renderer || !this.volume) return;
    const cam = this.renderer.getActiveCamera();
    if (typeof cam.dolly === "function") {
      cam.dolly(factor);
    }
    this.renderer.resetCameraClippingRange();
    this.renderWindow.render();
  }

  _parseChimeraxViewMatrixCamera(cameraArg) {
    if (!cameraArg) return null;
    var raw = String(cameraArg).trim();
    if (!raw) return null;
    if (raw.toLowerCase().startsWith("camera")) raw = raw.slice(6).trim();
    var parts = raw.match(/[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?/g);
    if (!parts || parts.length < 12) return null;
    var nums = parts.slice(0, 12).map(Number);
    if (!nums.every(function (v) { return isFinite(v); })) return null;
    return nums;
  }

  /** ChimeraX ``view matrix camera`` comma-separated floats (12 numbers).

  Translation columns are zeroed: ChimeraX runs ``volume center`` first, and a VTK
  view-matrix translation often frames the map off-screen (blank PNG + default retry).
  */
  getChimeraxViewMatrixCamera() {
    if (!this.renderer) return null;
    var cam = this.renderer.getActiveCamera();
    if (!cam || !cam.getViewMatrix) return null;
    var mat = new Float64Array(16);
    cam.getViewMatrix(mat);
    var out = [];
    for (var row = 0; row < 3; row++) {
      for (var col = 0; col < 4; col++) {
        out.push(mat[col * 4 + row]);
      }
    }
    out[3] = 0;
    out[7] = 0;
    out[11] = 0;
    if (!out.every(function(v) { return isFinite(v); })) return null;
    return out.map(function(v) { return Number(v).toPrecision(8); }).join(",");
  }

  /** Apply ChimeraX ``view matrix camera`` (12 numbers) to the VTK camera. */
  setChimeraxViewMatrixCamera(cameraArg) {
    if (!this.renderer) return false;
    var nums = this._parseChimeraxViewMatrixCamera(cameraArg);
    if (!nums) return false;
    var cam = this.renderer.getActiveCamera();
    if (!cam || !cam.setViewMatrix) return false;
    var mat = new Float64Array(16);
    mat[12] = 0;
    mat[13] = 0;
    mat[14] = 0;
    mat[15] = 1;
    for (var row = 0; row < 3; row++) {
      for (var col = 0; col < 4; col++) {
        mat[col * 4 + row] = nums[row * 4 + col];
      }
    }
    cam.setViewMatrix(mat);
    if (typeof cam.orthogonalizeViewUp === "function") cam.orthogonalizeViewUp();
    this.renderer.resetCameraClippingRange();
    this.renderWindow.render();
    this.syncChimeraxTurnsFromCamera();
    return true;
  }

  setInteractionEnabled(enabled) {
    if (!this.renderWindow) return;
    var interactor = this.renderWindow.getInteractor();
    if (!interactor || !interactor.setInteractorStyle) return;
    var style = interactor.getInteractorStyle();
    if (!style) return;
    if (typeof style.setEnabled === "function") {
      style.setEnabled(!!enabled);
    }
  }

  destroy() {
    this._clearVolume();
    if (this.ctfun) {
      this.ctfun.delete();
      this.ctfun = null;
    }
    if (this.ofun) {
      this.ofun.delete();
      this.ofun = null;
    }
    if (this._resizeObserver) {
      this._resizeObserver.disconnect();
      this._resizeObserver = null;
    }
    if (this.grw) {
      this.grw.setContainer(null);
      this.grw.delete();
      this.grw = null;
    }
    this.renderer = null;
    this.renderWindow = null;
  }

  observeResize(target) {
    if (typeof ResizeObserver === "undefined" || !target) return;
    var self = this;
    if (this._resizeObserver) this._resizeObserver.disconnect();
    this._resizeObserver = new ResizeObserver(function () {
      self.resize();
    });
    this._resizeObserver.observe(target);
  }
}

if (typeof window !== "undefined") {
  window.CryoVolumeRaycastView = VolumeRaycastView;
}
