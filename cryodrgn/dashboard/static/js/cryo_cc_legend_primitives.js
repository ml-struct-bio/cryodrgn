/**
 * Pure helpers + palette data for ``color_covariate_legend.js``.
 * Load this file before ``color_covariate_legend.js``.
 */
(function (global) {
  "use strict";

  var PALETTE_RGB = {
    Viridis: [
      [68, 1, 84], [72, 39, 119], [63, 74, 138], [49, 104, 142], [38, 130, 142],
      [31, 157, 138], [53, 183, 121], [109, 205, 89], [253, 231, 37]
    ],
    Plasma: [
      [13, 8, 135], [76, 2, 161], [126, 3, 168], [170, 35, 149], [204, 71, 120],
      [230, 108, 92], [248, 149, 64], [253, 197, 39], [240, 249, 33]
    ],
    Inferno: [[0, 0, 4], [41, 10, 108], [108, 2, 112], [187, 55, 84], [249, 142, 8], [252, 208, 74], [252, 254, 164]],
    Magma: [[0, 0, 4], [59, 15, 112], [122, 4, 111], [203, 62, 65], [252, 141, 89], [252, 208, 136], [252, 253, 191]],
    Cividis: [
      [0, 34, 78], [0, 57, 112], [68, 90, 129], [115, 128, 131], [159, 161, 135],
      [206, 186, 122], [253, 231, 37]
    ],
    Turbo: [[48, 18, 59], [61, 99, 221], [26, 152, 223], [25, 189, 114], [208, 230, 28], [250, 85, 8], [144, 12, 2]],
    Blues: [[247, 251, 255], [198, 219, 239], [107, 174, 214], [33, 113, 181], [8, 48, 107]],
    Greens: [[247, 252, 245], [199, 233, 192], [116, 196, 118], [35, 139, 69], [0, 68, 27]],
    Greys: [[255, 255, 255], [217, 217, 217], [150, 150, 150], [82, 82, 82], [0, 0, 0]],
    Oranges: [[255, 245, 235], [253, 208, 162], [253, 141, 60], [217, 72, 1], [127, 39, 4]],
    Purples: [[252, 251, 253], [218, 218, 235], [158, 154, 200], [106, 81, 163], [63, 0, 125]],
    Reds: [[255, 245, 240], [252, 187, 161], [251, 106, 74], [203, 24, 29], [103, 0, 13]],
    YlGnBu: [[255, 255, 217], [199, 233, 180], [65, 182, 196], [44, 127, 184], [8, 29, 88]],
    YlOrRd: [[255, 255, 204], [255, 237, 160], [254, 178, 76], [240, 59, 32], [128, 0, 38]],
    RdBu: [[103, 0, 31], [214, 96, 77], [247, 247, 247], [67, 147, 195], [5, 48, 97]],
    Portland: [[12, 51, 131], [10, 136, 186], [242, 211, 56], [242, 143, 56], [217, 30, 30]],
    Jet: [[0, 0, 127], [0, 0, 255], [0, 255, 255], [255, 255, 0], [255, 0, 0], [127, 0, 0]],
    Hot: [[0, 0, 0], [176, 0, 0], [255, 140, 0], [255, 255, 102], [255, 255, 255]],
    Blackbody: [[0, 0, 0], [91, 12, 107], [184, 50, 137], [242, 109, 75], [249, 248, 113]],
    Electric: [[0, 0, 0], [30, 0, 109], [141, 23, 165], [217, 79, 213], [255, 255, 255]],
    Rainbow: [[110, 64, 170], [71, 118, 230], [26, 199, 194], [123, 217, 76], [244, 208, 63], [242, 95, 92]],
    Earth: [[0, 0, 130], [0, 181, 200], [125, 190, 66], [200, 182, 74], [139, 69, 19]]
  };

  var DISCRETE_CB_SEL = 'input[type="checkbox"][data-cat-key]';

  function discreteCategoryInputsFrom(root, checkedOnly) {
    if (!root) {
      return [];
    }
    return root.querySelectorAll(DISCRETE_CB_SEL + (checkedOnly ? ":checked" : ""));
  }

  function normalizeDiscreteLegendHex(hex) {
    if (hex == null) return "";
    var s = String(hex).trim();
    if (/^#[0-9a-fA-F]{6}$/.test(s)) {
      return s.toLowerCase();
    }
    if (/^[0-9a-fA-F]{6}$/.test(s)) {
      return ("#" + s).toLowerCase();
    }
    return "";
  }

  function discreteLegendRgbTuple(hexNorm) {
    var h = normalizeDiscreteLegendHex(hexNorm);
    if (!h || h.length !== 7) return { r: 0, g: 0, b: 0 };
    return {
      r: parseInt(h.slice(1, 3), 16),
      g: parseInt(h.slice(3, 5), 16),
      b: parseInt(h.slice(5, 7), 16)
    };
  }

  function cryoDiscreteColorWheelSvg() {
    return "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 24 24\" aria-hidden=\"true\""
      + " focusable=\"false\"><defs>"
      + "<linearGradient id=\"cryo-cc-rainbow\" x1=\"0%\" y1=\"0%\" x2=\"100%\" y2=\"100%\">"
      + "<stop offset=\"0%\" stop-color=\"#ff0000\"/>"
      + "<stop offset=\"16%\" stop-color=\"#ff8000\"/>"
      + "<stop offset=\"33%\" stop-color=\"#ffff00\"/>"
      + "<stop offset=\"50%\" stop-color=\"#00ff00\"/>"
      + "<stop offset=\"66%\" stop-color=\"#0080ff\"/>"
      + "<stop offset=\"83%\" stop-color=\"#4000ff\"/>"
      + "<stop offset=\"100%\" stop-color=\"#ff00ff\"/>"
      + "</linearGradient></defs>"
      + "<circle cx=\"12\" cy=\"12\" r=\"9\" stroke=\"rgba(36,59,83,0.35)\""
      + " stroke-width=\"1.2\" fill=\"url(#cryo-cc-rainbow)\"/>"
      + "<circle cx=\"12\" cy=\"12\" r=\"3.5\" fill=\"white\" opacity=\"0.9\"/>"
      + "</svg>";
  }

  function paletteScaleCSS(paletteName, t) {
    t = Math.max(0, Math.min(1, t));
    var pal = PALETTE_RGB[paletteName] || PALETTE_RGB.Viridis;
    var n = pal.length - 1;
    var i = Math.min(Math.floor(t * n), n - 1);
    var f = t * n - i;
    var a = pal[i];
    var b = pal[i + 1];
    return "rgb(" + Math.round(a[0] + (b[0] - a[0]) * f) + ","
      + Math.round(a[1] + (b[1] - a[1]) * f) + ","
      + Math.round(a[2] + (b[2] - a[2]) * f) + ")";
  }

  function formatThresholdValue(v) {
    if (typeof v !== "number" || !isFinite(v)) return "";
    var av = Math.abs(v);
    if (av !== 0 && (av < 0.001 || av >= 10000)) return v.toExponential(3);
    return String(Math.round(v * 1000) / 1000);
  }

  function rangeInequalityPhrase(lo, hi, label) {
    return formatThresholdValue(Math.min(lo, hi)) + " <= " + label + " <= "
      + formatThresholdValue(Math.max(lo, hi));
  }

  function quantileSorted(sorted, q) {
    if (!sorted.length) return NaN;
    var pos = (sorted.length - 1) * q;
    var lo = Math.floor(pos);
    var hi = Math.ceil(pos);
    if (lo === hi) return sorted[lo];
    var f = pos - lo;
    return sorted[lo] * (1 - f) + sorted[hi] * f;
  }

  function violinDistribution(xs, paletteName) {
    var n = xs.length;
    var xmin = Infinity;
    var xmax = -Infinity;
    for (var i = 0; i < n; i++) {
      if (xs[i] < xmin) xmin = xs[i];
      if (xs[i] > xmax) xmax = xs[i];
    }
    if (!isFinite(xmin) || !isFinite(xmax)) {
      return { outlineX: [], outlineY: [], fillTraces: [], range: [0, 1] };
    }
    var sorted = xs.slice().sort(function (a, b) { return a - b; });
    var q1 = quantileSorted(sorted, 0.25);
    var q3 = quantileSorted(sorted, 0.75);
    var iqr = q3 - q1;
    var loFence = q1 - 3.0 * iqr;
    var hiFence = q3 + 3.0 * iqr;
    var vals = [];
    for (var oi = 0; oi < n; oi++) {
      var xv = xs[oi];
      if (iqr > 0 && (xv < loFence || xv > hiFence)) continue;
      vals.push(xv);
    }
    if (!vals.length) vals = xs.slice();
    if (xmax <= xmin) {
      var oneColor = paletteScaleCSS(paletteName, 0.5);
      return {
        outlineX: [xmin, xmin, xmin],
        outlineY: [-0.38, 0.38, -0.38],
        fillTraces: [{
          type: "scatter",
          mode: "lines",
          x: [xmin - 0.5, xmin + 0.5, xmin + 0.5, xmin - 0.5, xmin - 0.5],
          y: [-0.38, -0.38, 0.38, 0.38, -0.38],
          fill: "toself",
          fillcolor: oneColor,
          line: { width: 0, color: oneColor },
          hoverinfo: "skip",
          showlegend: false
        }],
        range: [xmin - 0.6, xmax + 0.6]
      };
    }
    var gridN = 96;
    var span = xmax - xmin;
    var mean = 0;
    for (var mi = 0; mi < vals.length; mi++) mean += vals[mi];
    mean /= vals.length;
    var variance = 0;
    for (var vi = 0; vi < vals.length; vi++) {
      var dv = vals[vi] - mean;
      variance += dv * dv;
    }
    var std = Math.sqrt(variance / Math.max(1, vals.length - 1));
    var bandwidth = 1.06 * (std || span / 6) * Math.pow(vals.length, -0.2);
    if (!isFinite(bandwidth) || bandwidth <= 0) bandwidth = span / 24;
    bandwidth = Math.max(bandwidth, span / 200);
    var gridX = [];
    var dens = [];
    var maxD = 0;
    var invTwoBw2 = 1 / (2 * bandwidth * bandwidth);
    for (var g = 0; g < gridN; g++) {
      var gx = xmin + span * g / (gridN - 1);
      var dsum = 0;
      for (var k = 0; k < vals.length; k++) {
        var dx = gx - vals[k];
        dsum += Math.exp(-(dx * dx) * invTwoBw2);
      }
      var d = dsum / vals.length;
      if (d > maxD) maxD = d;
      gridX.push(gx);
      dens.push(d);
    }
    var half = [];
    for (var hi = 0; hi < dens.length; hi++) {
      half.push(maxD > 0 ? 0.42 * dens[hi] / maxD : 0);
    }
    var outlineX = [];
    var outlineY = [];
    for (var ux = 0; ux < gridX.length; ux++) {
      outlineX.push(gridX[ux]);
      outlineY.push(half[ux]);
    }
    for (var lx = gridX.length - 1; lx >= 0; lx--) {
      outlineX.push(gridX[lx]);
      outlineY.push(-half[lx]);
    }
    outlineX.push(gridX[0]);
    outlineY.push(half[0]);
    var fillTraces = [];
    var palName = paletteName;
    for (var s = 0; s < gridX.length - 1; s++) {
      var x0 = gridX[s];
      var x1 = gridX[s + 1];
      var y = Math.max(half[s], half[s + 1]);
      if (y <= 0) continue;
      var t = ((x0 + x1) * 0.5 - xmin) / span;
      var fillColor = paletteScaleCSS(palName, t);
      fillTraces.push({
        type: "scatter",
        mode: "lines",
        x: [x0, x1, x1, x0, x0],
        y: [-y, -y, y, y, -y],
        fill: "toself",
        fillcolor: fillColor,
        line: { width: 0, color: fillColor },
        hoverinfo: "skip",
        showlegend: false
      });
    }
    var pad = span * 0.04;
    return {
      outlineX: outlineX,
      outlineY: outlineY,
      fillTraces: fillTraces,
      range: [xmin - pad, xmax + pad]
    };
  }

  function violinDistributionVerticalFromHorizontal(v) {
    var fillTracesV = [];
    for (var s = 0; s < v.fillTraces.length; s++) {
      var t = v.fillTraces[s];
      fillTracesV.push({
        type: "scatter",
        mode: "lines",
        x: t.y.slice(),
        y: t.x.slice(),
        fill: t.fill,
        fillcolor: t.fillcolor,
        line: t.line,
        hoverinfo: "skip",
        showlegend: false
      });
    }
    return {
      outlineX: v.outlineY.slice(),
      outlineY: v.outlineX.slice(),
      fillTraces: fillTracesV,
      range: v.range
    };
  }

  function parseMarkerHexRgb(hex) {
    if (typeof hex !== "string") return null;
    var h = hex.trim();
    if (h.charAt(0) !== "#") return null;
    h = h.slice(1);
    var r;
    var g;
    var b;
    if (h.length === 3) {
      r = parseInt(h.charAt(0) + h.charAt(0), 16);
      g = parseInt(h.charAt(1) + h.charAt(1), 16);
      b = parseInt(h.charAt(2) + h.charAt(2), 16);
    } else if (h.length === 6) {
      r = parseInt(h.slice(0, 2), 16);
      g = parseInt(h.slice(2, 4), 16);
      b = parseInt(h.slice(4, 6), 16);
    } else {
      return null;
    }
    if (!isFinite(r) || !isFinite(g) || !isFinite(b)) return null;
    return { r: r, g: g, b: b };
  }

  function rgbChannelByte(n) {
    var x = Math.round(Number(n));
    if (!isFinite(x)) return 0;
    return Math.max(0, Math.min(255, x));
  }

  function rgbToHex6(r, g, b) {
    function hb(v) {
      var s = rgbChannelByte(v).toString(16);
      return s.length === 1 ? "0" + s : s;
    }
    return ("#" + hb(r) + hb(g) + hb(b)).toLowerCase();
  }

  function discreteCardStyles(hex) {
    var rgb = parseMarkerHexRgb(hex);
    if (!rgb) return { bg: "", borderColor: "", labelColor: "", countColor: "" };
    var r = rgb.r;
    var g = rgb.g;
    var b = rgb.b;
    var towardWhite = 0.3;
    var rf = Math.round(r * (1 - towardWhite) + 255 * towardWhite);
    var gf = Math.round(g * (1 - towardWhite) + 255 * towardWhite);
    var bf = Math.round(b * (1 - towardWhite) + 255 * towardWhite);
    rf = Math.max(0, Math.min(255, rf));
    gf = Math.max(0, Math.min(255, gf));
    bf = Math.max(0, Math.min(255, bf));
    var textDim = 0.46;
    var rt = Math.max(0, Math.min(255, Math.round(r * textDim)));
    var gt = Math.max(0, Math.min(255, Math.round(g * textDim)));
    var bt = Math.max(0, Math.min(255, Math.round(b * textDim)));
    return {
      bg: "rgb(" + rf + "," + gf + "," + bf + ")",
      borderColor: "rgba(" + r + "," + g + "," + b + ",0.38)",
      labelColor: "rgb(" + rt + "," + gt + "," + bt + ")",
      countColor: "rgba(" + rt + "," + gt + "," + bt + ",0.88)"
    };
  }

  function discretePlotMatchedChipStyles(hex, blend, plastic) {
    var rgb0 = parseMarkerHexRgb(hex);
    if (!rgb0) return { bg: "", borderColor: "", labelColor: "", countColor: "" };
    var r = rgb0.r;
    var g = rgb0.g;
    var b = rgb0.b;
    if (blend && blend.alpha != null && blend.bgHex) {
      var bgHx = normalizeDiscreteLegendHex(blend.bgHex) || blend.bgHex;
      var bg = parseMarkerHexRgb(bgHx.charAt(0) === "#" ? bgHx : "#" + bgHx);
      if (bg) {
        var a = Number(blend.alpha);
        if (isFinite(a)) {
          a = Math.max(0, Math.min(1, a));
          r = Math.round(r * a + bg.r * (1 - a));
          g = Math.round(g * a + bg.g * (1 - a));
          b = Math.round(b * a + bg.b * (1 - a));
        }
      }
    }
    var lum = (0.2126 * r + 0.7152 * g + 0.0722 * b) / 255;
    var ba = plastic ? 0.58 : 0.48;
    return {
      bg: "rgb(" + r + "," + g + "," + b + ")",
      borderColor: "rgba(" + r + "," + g + "," + b + "," + ba + ")",
      labelColor: lum > 0.52 ? "#131820" : "#f8fafc",
      countColor: lum > 0.52 ? "rgba(36,59,83,0.88)" : "rgba(248,250,252,0.82)"
    };
  }

  function formatDiscreteDatasetCount(n) {
    var v = Number(n);
    if (!isFinite(v)) v = 0;
    var nlab = v.toLocaleString("en-CA");
    return nlab + " pts";
  }

  function discreteDisplayLabel(key, rawCol, displayMap) {
    if (key === "__na__") return "(missing)";
    var disp = (displayMap && rawCol) ? String(displayMap[rawCol] || rawCol) : "";
    var rc = String(rawCol || "");
    if (rc === "landscape_vol_cluster") return "Vol cluster=" + String(key);
    var isKmeans = /k[\s_-]*means/i.test(rc) || /k[\s_-]*means/i.test(disp) || /^labels$/i.test(rc);
    if (isKmeans) return "kmeans=" + String(key);
    return String(key);
  }

  function sortDiscreteKeys(keys) {
    var rest = keys.filter(function (k) { return k !== "__na__"; });
    var allNum = rest.length > 0 && rest.every(function (k) {
      var num = Number(k);
      return isFinite(num) && String(num) === String(k).trim();
    });
    if (allNum) {
      rest.sort(function (a, b) { return Number(a) - Number(b); });
    } else {
      rest.sort(function (a, b) {
        return String(a).localeCompare(String(b), undefined, { numeric: true });
      });
    }
    if (keys.indexOf("__na__") >= 0) rest.push("__na__");
    return rest;
  }

  /**
   * Fit invert-chip label typography inside its cell (pair plot / explorer / 3-D footer).
   * @param {HTMLElement} invertCell
   */
  function fitInvertCellTypography(invertCell) {
    var fitBox = invertCell.querySelector(".cryo-cc-discrete-invert-fitbox");
    var label = invertCell.querySelector(".cryo-cc-discrete-switch-label");
    if (!fitBox || !label) return;
    if (
      invertCell.classList.contains("cryo-cc-discrete-cell--invert-heading")
      && invertCell.hidden
    ) {
      return;
    }
    if (invertCell.classList.contains("cryo-cc-discrete-cell--invert-heading")) {
      label.style.fontSize = "";
      label.style.whiteSpace = "nowrap";
      label.style.overflow = "";
      label.style.textOverflow = "";
      return;
    }
    var innerW = fitBox.clientWidth;
    var innerH = fitBox.clientHeight;
    if (typeof window !== "undefined" && window.getComputedStyle) {
      var bxs = window.getComputedStyle(fitBox);
      var pl = parseFloat(bxs.paddingLeft) || 0;
      var pr = parseFloat(bxs.paddingRight) || 0;
      var pt = parseFloat(bxs.paddingTop) || 0;
      var pb = parseFloat(bxs.paddingBottom) || 0;
      innerW -= pl + pr;
      innerH -= pt + pb;
    }
    if (!(innerW > 2 && innerH > 2)) {
      return;
    }
    label.style.fontSize = "";
    label.style.whiteSpace = "";
    var raw = (label.textContent || "").replace(/\r/g, "").trim();
    if (!raw) return;

    var cs = typeof window !== "undefined" && window.getComputedStyle
      ? window.getComputedStyle(label)
      : null;
    var probe = document.createElement("span");
    probe.setAttribute("aria-hidden", "true");
    probe.style.cssText =
      "visibility:hidden;position:absolute;left:-9999px;top:0;white-space:nowrap;padding:0;margin:0;";
    if (cs) {
      probe.style.fontFamily = cs.fontFamily;
      probe.style.fontWeight = cs.fontWeight;
      probe.style.fontStyle = cs.fontStyle;
      probe.style.letterSpacing = cs.letterSpacing;
    }
    document.body.appendChild(probe);

    var padPx = 2;
    var vertSlack = 6;

    if (invertCell.classList.contains("cryo-cc-discrete-cell--invert-footer")) {
      label.style.whiteSpace = "nowrap";
      var phrase = raw.replace(/\s+/g, " ").trim();
      var loF = 5;
      var hiF = 144;
      var bestF = loF;
      try {
        while (loF <= hiF) {
          var midF = (loF + hiF) >> 1;
          label.style.fontSize = midF + "px";
          probe.style.fontSize = midF + "px";
          probe.textContent = phrase;
          var wOkF = probe.offsetWidth <= innerW + padPx;
          var hOkF = label.scrollHeight <= innerH - slackF;
          if (wOkF && hOkF) {
            bestF = midF;
            loF = midF + 1;
          } else {
            hiF = midF - 1;
          }
        }
        label.style.fontSize = bestF + "px";
        while (bestF > 5 && label.scrollHeight > innerH - slackF) {
          bestF -= 1;
          label.style.fontSize = bestF + "px";
        }
        probe.style.fontSize = bestF + "px";
        probe.textContent = phrase;
        while (bestF > 5 && probe.offsetWidth > innerW + padPx) {
          bestF -= 1;
          label.style.fontSize = bestF + "px";
          probe.style.fontSize = bestF + "px";
        }
      } finally {
        if (probe.parentNode) probe.parentNode.removeChild(probe);
      }
      label.style.whiteSpace = "";
      return;
    }

    var lines = raw.split("\n").map(function (ln) {
      return ln.replace(/\s+/g, " ").trim();
    }).filter(Boolean);
    var words = [];
    for (var li = 0; li < lines.length; li++) {
      var parts = lines[li].split(/\s+/).filter(Boolean);
      for (var pi = 0; pi < parts.length; pi++) {
        words.push(parts[pi]);
      }
    }
    if (!words.length) return;

    var minPx = 5;
    var lo = minPx;
    var hi = 144;
    var best = lo;
    try {
      while (lo <= hi) {
        var mid = (lo + hi) >> 1;
        label.style.fontSize = mid + "px";
        probe.style.fontSize = mid + "px";
        var wordsFit = true;
        for (var wi = 0; wi < words.length; wi++) {
          probe.textContent = words[wi];
          if (probe.offsetWidth > innerW + padPx) {
            wordsFit = false;
            break;
          }
        }
        var fits =
          wordsFit &&
          label.scrollWidth <= innerW + padPx &&
          label.scrollHeight <= innerH - vertSlack;
        if (fits) {
          best = mid;
          lo = mid + 1;
        } else {
          hi = mid - 1;
        }
      }
      label.style.fontSize = best + "px";
      while (best > minPx && label.scrollHeight > innerH - vertSlack) {
        best -= 1;
        label.style.fontSize = best + "px";
      }
    } finally {
      if (probe.parentNode) probe.parentNode.removeChild(probe);
    }
  }

  global.CryoCcLegendPrimitives = {
    PALETTE_RGB: PALETTE_RGB,
    paletteScaleCSS: paletteScaleCSS,
    quantileSorted: quantileSorted,
    violinDistribution: violinDistribution,
    violinDistributionVerticalFromHorizontal: violinDistributionVerticalFromHorizontal,
    formatThresholdValue: formatThresholdValue,
    rangeInequalityPhrase: rangeInequalityPhrase,
    parseMarkerHexRgb: parseMarkerHexRgb,
    rgbChannelByte: rgbChannelByte,
    rgbToHex6: rgbToHex6,
    discreteCardStyles: discreteCardStyles,
    discretePlotMatchedChipStyles: discretePlotMatchedChipStyles,
    normalizeDiscreteLegendHex: normalizeDiscreteLegendHex,
    discreteLegendRgbTuple: discreteLegendRgbTuple,
    cryoDiscreteColorWheelSvg: cryoDiscreteColorWheelSvg,
    formatDiscreteDatasetCount: formatDiscreteDatasetCount,
    discreteDisplayLabel: discreteDisplayLabel,
    sortDiscreteKeys: sortDiscreteKeys,
    discreteCategoryInputsFrom: discreteCategoryInputsFrom,
    fitInvertCellTypography: fitInvertCellTypography
  };
})(typeof window !== "undefined" ? window : globalThis);
