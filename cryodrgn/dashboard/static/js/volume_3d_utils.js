/**
 * Shared client-side volume utilities for 3D preview and slice canvas.
 */
(function (global) {
  "use strict";

  var PLOT3D_TARGET_D = 128;
  var DEFAULT_ISO_PERCENTILE = 42;

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

  function trilinearSample(vol, srcD, x, y, z) {
    if (x < 0 || y < 0 || z < 0 || x > srcD - 1 || y > srcD - 1 || z > srcD - 1) {
      return 0;
    }
    var x0 = Math.floor(x);
    var y0 = Math.floor(y);
    var z0 = Math.floor(z);
    var x1 = Math.min(x0 + 1, srcD - 1);
    var y1 = Math.min(y0 + 1, srcD - 1);
    var z1 = Math.min(z0 + 1, srcD - 1);
    var xd = x - x0;
    var yd = y - y0;
    var zd = z - z0;
    var i000 = vol[x0 * srcD * srcD + y0 * srcD + z0];
    var i100 = vol[x1 * srcD * srcD + y0 * srcD + z0];
    var i010 = vol[x0 * srcD * srcD + y1 * srcD + z0];
    var i110 = vol[x1 * srcD * srcD + y1 * srcD + z0];
    var i001 = vol[x0 * srcD * srcD + y0 * srcD + z1];
    var i101 = vol[x1 * srcD * srcD + y0 * srcD + z1];
    var i011 = vol[x0 * srcD * srcD + y1 * srcD + z1];
    var i111 = vol[x1 * srcD * srcD + y1 * srcD + z1];
    var c00 = i000 * (1 - xd) + i100 * xd;
    var c01 = i001 * (1 - xd) + i101 * xd;
    var c10 = i010 * (1 - xd) + i110 * xd;
    var c11 = i011 * (1 - xd) + i111 * xd;
    var c0 = c00 * (1 - yd) + c10 * yd;
    var c1 = c01 * (1 - yd) + c11 * yd;
    return c0 * (1 - zd) + c1 * zd;
  }

  /** Box-average downsample to ``targetD``³ (prefilter before 3D raycast). */
  function downsampleVolumeBoxAverage(vol, srcD, targetD) {
    srcD = srcD | 0;
    targetD = targetD | 0;
    if (srcD === targetD) return vol;
    if (targetD < 1) throw new Error("targetD must be positive.");
    var out = new Float32Array(targetD * targetD * targetD);
    var scale = srcD / targetD;
    for (var iz = 0; iz < targetD; iz++) {
      var z0 = Math.floor(iz * scale);
      var z1 = Math.min(srcD, Math.ceil((iz + 1) * scale));
      for (var iy = 0; iy < targetD; iy++) {
        var y0 = Math.floor(iy * scale);
        var y1 = Math.min(srcD, Math.ceil((iy + 1) * scale));
        for (var ix = 0; ix < targetD; ix++) {
          var x0 = Math.floor(ix * scale);
          var x1 = Math.min(srcD, Math.ceil((ix + 1) * scale));
          var sum = 0;
          var count = 0;
          for (var z = z0; z < z1; z++) {
            for (var y = y0; y < y1; y++) {
              for (var x = x0; x < x1; x++) {
                sum += vol[x * srcD * srcD + y * srcD + z];
                count++;
              }
            }
          }
          out[ix * targetD * targetD + iy * targetD + iz] = count ? sum / count : 0;
        }
      }
    }
    return out;
  }

  /** Trilinear resample to ``targetD``³ (smoother alternative). */
  function downsampleVolumeTrilinear(vol, srcD, targetD) {
    srcD = srcD | 0;
    targetD = targetD | 0;
    if (srcD === targetD) return vol;
    var out = new Float32Array(targetD * targetD * targetD);
    for (var iz = 0; iz < targetD; iz++) {
      var z = iz * (srcD - 1) / Math.max(1, targetD - 1);
      for (var iy = 0; iy < targetD; iy++) {
        var y = iy * (srcD - 1) / Math.max(1, targetD - 1);
        for (var ix = 0; ix < targetD; ix++) {
          var x = ix * (srcD - 1) / Math.max(1, targetD - 1);
          out[ix * targetD * targetD + iy * targetD + iz] = trilinearSample(vol, srcD, x, y, z);
        }
      }
    }
    return out;
  }

  function volumePercentileSamples(values, maxSamples) {
    maxSamples = maxSamples || 65536;
    var flat = values;
    var samples = [];
    var step = Math.max(1, Math.floor(flat.length / maxSamples));
    for (var i = 0; i < flat.length; i += step) {
      var v = flat[i];
      if (isFinite(v)) samples.push(v);
    }
    samples.sort(function(a, b) { return a - b; });
    if (!samples.length) return [0, 1];
    return samples;
  }

  function percentileValue(samples, pct) {
    if (!samples || !samples.length) return 0;
    pct = Math.max(0, Math.min(100, Number(pct)));
    var idx = (pct / 100) * (samples.length - 1);
    var lo = Math.floor(idx);
    var hi = Math.ceil(idx);
    if (lo === hi) return samples[lo];
    var frac = idx - lo;
    return samples[lo] * (1 - frac) + samples[hi] * frac;
  }

  function suggestIsoPercentile(samples) {
    if (!samples || samples.length < 8) return DEFAULT_ISO_PERCENTILE;
    var vmin = samples[0];
    var vmax = samples[samples.length - 1];
    var span = vmax - vmin;
    if (!(span > 0)) return DEFAULT_ISO_PERCENTILE;
    var p10 = percentileValue(samples, 10);
    var p25 = percentileValue(samples, 25);
    var p50 = percentileValue(samples, 50);
    var p90 = percentileValue(samples, 90);
    var p99 = percentileValue(samples, 99);
    var baseline = p10 + (p25 - p10) * 0.35;
    var signal = p99 - baseline;
    if (!(signal > 0)) return DEFAULT_ISO_PERCENTILE;
    var medianFrac = (p50 - vmin) / span;
    var tailFrac = (vmax - p90) / span;
    var alpha = 0.24;
    if (medianFrac < 0.12) alpha = 0.18;
    else if (medianFrac < 0.22) alpha = 0.21;
    else if (medianFrac > 0.4) alpha = 0.32;
    if (tailFrac < 0.08) alpha += 0.06;
    var target = baseline + signal * alpha;
    for (var p = 5; p < 99; p++) {
      if (percentileValue(samples, p) >= target) {
        return Math.max(12, Math.min(94, Math.round(p)));
      }
    }
    return 88;
  }

  function suggestIsoDataValue(samples) {
    return percentileValue(samples, suggestIsoPercentile(samples));
  }

  /** Rank window for the iso slider (linear in data values, ChimeraX map units). */
  var ISO_SLIDER_RANK_LO = 2;
  var ISO_SLIDER_RANK_HI = 99.5;

  function isoSliderDataRange(samples) {
    if (!samples || samples.length < 2) return { min: 0, max: 1 };
    var lo = percentileValue(samples, ISO_SLIDER_RANK_LO);
    var hi = percentileValue(samples, ISO_SLIDER_RANK_HI);
    if (!(hi > lo)) {
      lo = samples[0];
      hi = samples[samples.length - 1];
    }
    return { min: lo, max: hi };
  }

  /** Slider 0–100 ↔ data value, linear between slider endpoints (ChimeraX contour units). */
  function isoSliderToDataValue(slider, vmin, vmax) {
    vmin = Number(vmin);
    vmax = Number(vmax);
    var t = Math.max(0, Math.min(100, Number(slider))) / 100;
    if (!isFinite(vmin) || !isFinite(vmax)) return 0;
    if (vmax <= vmin) return vmin;
    return vmin + t * (vmax - vmin);
  }

  function isoDataValueToSlider(value, vmin, vmax) {
    vmin = Number(vmin);
    vmax = Number(vmax);
    value = Number(value);
    var span = vmax - vmin;
    if (!isFinite(span) || !(span > 0)) return 0;
    var t = Math.max(0, Math.min(1, (value - vmin) / span));
    return 100 * t;
  }

  function isoRangeFromPercentile(samples, isoPercentile) {
    if (!samples || !samples.length) return { isomin: 0, isomax: 1 };
    var pct = Math.max(0, Math.min(100, Number(isoPercentile)));
    var isomin = percentileValue(samples, Math.max(0, pct - 1));
    var isomax = percentileValue(samples, Math.min(100, Math.max(pct + 12, 95)));
    if (isomax <= isomin) isomax = samples[samples.length - 1];
    return { isomin: isomin, isomax: isomax };
  }

  global.CryoVolume3dUtils = {
    PLOT3D_TARGET_D: PLOT3D_TARGET_D,
    decodeFloat32Volume: decodeFloat32Volume,
    downsampleVolumeBoxAverage: downsampleVolumeBoxAverage,
    downsampleVolumeTrilinear: downsampleVolumeTrilinear,
    volumePercentileSamples: volumePercentileSamples,
    suggestIsoPercentile: suggestIsoPercentile,
    suggestIsoDataValue: suggestIsoDataValue,
    isoSliderDataRange: isoSliderDataRange,
    isoSliderToDataValue: isoSliderToDataValue,
    isoDataValueToSlider: isoDataValueToSlider,
    isoRangeFromPercentile: isoRangeFromPercentile,
    percentileValue: percentileValue
  };
})(typeof window !== "undefined" ? window : this);
