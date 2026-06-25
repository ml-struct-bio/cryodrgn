/**
 * vtk.js volume raycasting viewer for cryoDRGN dashboard (single volume, rotate only).
 */
import "@kitware/vtk.js/Rendering/Profiles/Volume";
import vtkColorTransferFunction from "@kitware/vtk.js/Rendering/Core/ColorTransferFunction";
import vtkDataArray from "@kitware/vtk.js/Common/Core/DataArray";
import vtkImageData from "@kitware/vtk.js/Common/DataModel/ImageData";
import vtkPiecewiseFunction from "@kitware/vtk.js/Common/DataModel/PiecewiseFunction";
import vtkGenericRenderWindow from "@kitware/vtk.js/Rendering/Misc/GenericRenderWindow";
import vtkVolume from "@kitware/vtk.js/Rendering/Core/Volume";
import vtkVolumeMapper from "@kitware/vtk.js/Rendering/Core/VolumeMapper";

var VOLUME_BG = [0.98, 0.97, 0.96];
var DEFAULT_ISO_LEVEL = 42;

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

function applyTransferFunctions(ctfun, ofun, range, percentileSamples, isoPercentile) {
  var U = utils();
  var pct = clamp(Number(isoPercentile), 0, 100);
  var span = range.max - range.min;
  var minSep = span > 0 ? span * 0.0015 : 1e-6;

  var threshold = U.percentileValue(percentileSamples, pct);
  var rampPctLo = Math.max(0, pct - 4 - (100 - pct) * 0.04);
  var rampPctHi = Math.min(100, pct + 8 + (100 - pct) * 0.06);
  var rampLo = U.percentileValue(percentileSamples, rampPctLo);
  var rampHi = U.percentileValue(percentileSamples, rampPctHi);

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
  ofun.addPoint(threshold, 0.22);
  ofun.addPoint(rampHi, 0.38);
  ofun.addPoint(range.max, 0.42);
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
    this.percentileSamples = null;
    this.isoLevel = DEFAULT_ISO_LEVEL;
  }

  _reapplyTransferFunctions() {
    if (!this.volume || !this.volumeRange || !this.percentileSamples || !this.ctfun || !this.ofun) {
      return;
    }
    applyTransferFunctions(
      this.ctfun,
      this.ofun,
      this.volumeRange,
      this.percentileSamples,
      this.isoLevel
    );
    if (this.renderWindow) this.renderWindow.render();
  }

  getIsoLevel() {
    return this.isoLevel;
  }

  setIsoLevel(level) {
    this.isoLevel = clamp(Number(level), 0, 100);
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

  setVolumeFromB64(b64, d) {
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
    this.isoLevel = U.suggestIsoPercentile(this.percentileSamples);
    applyTransferFunctions(
      this.ctfun,
      this.ofun,
      this.volumeRange,
      this.percentileSamples,
      this.isoLevel
    );

    this.volume = vtkVolume.newInstance();
    this.volume.setMapper(this.mapper);
    var prop = this.volume.getProperty();
    prop.setRGBTransferFunction(0, this.ctfun);
    prop.setScalarOpacity(0, this.ofun);
    prop.setScalarOpacityUnitDistance(0, 2.5);
    prop.setInterpolationTypeToLinear();
    prop.setShade(true);
    prop.setAmbient(0.25);
    prop.setDiffuse(0.75);
    prop.setSpecular(0.15);

    this.renderer.addVolume(this.volume);
    this.renderer.resetCamera();
    this.renderer.getActiveCamera().azimuth(25);
    this.renderer.getActiveCamera().elevation(12);
    this.renderer.resetCameraClippingRange();
    this.resize();
    this.renderWindow.render();
  }

  resize() {
    if (this.grw) this.grw.resize();
  }

  resetCamera() {
    if (!this.volume || !this.renderer) return;
    this.renderer.resetCamera();
    this.renderer.getActiveCamera().azimuth(25);
    this.renderer.getActiveCamera().elevation(12);
    this.renderer.resetCameraClippingRange();
    this.renderWindow.render();
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
