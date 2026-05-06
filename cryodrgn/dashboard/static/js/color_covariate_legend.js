/**
 * Colour covariate histogram + discrete toggles (particle-explorer style) for
 * pair-grid and 3-D latent views. Optional `histPlotVertical`: covariate on y,
 * density on x (pair-plot aside). Optional `onDiscreteLayout` / `onContinuousHistogramLayout`
 * run after DOM / Plotly layout. Panel `vertical` controls CSS layout class.
 */
(function (global) {
  "use strict";

  var PALETTE_RGB = {
    Viridis: [[68, 1, 84], [72, 39, 119], [63, 74, 138], [49, 104, 142], [38, 130, 142], [31, 157, 138], [53, 183, 121], [109, 205, 89], [253, 231, 37]],
    Plasma: [[13, 8, 135], [76, 2, 161], [126, 3, 168], [170, 35, 149], [204, 71, 120], [230, 108, 92], [248, 149, 64], [253, 197, 39], [240, 249, 33]],
    Inferno: [[0, 0, 4], [41, 10, 108], [108, 2, 112], [187, 55, 84], [249, 142, 8], [252, 208, 74], [252, 254, 164]],
    Magma: [[0, 0, 4], [59, 15, 112], [122, 4, 111], [203, 62, 65], [252, 141, 89], [252, 208, 136], [252, 253, 191]],
    Cividis: [[0, 34, 78], [0, 57, 112], [68, 90, 129], [115, 128, 131], [159, 161, 135], [206, 186, 122], [253, 231, 37]],
    Turbo: [[48, 18, 59], [61, 99, 221], [26, 152, 223], [25, 189, 114], [208, 230, 28], [250, 85, 8], [144, 12, 2]]
  };

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

  /** Swap axes so covariate runs vertically (y) and density horizontally (x). */
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

  /** Parse #rgb / #rrggbb → {r,g,b} or null. */
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

  /**
   * Pastel chip fill from legend hex (matches particle explorer montage chips) plus readable text.
   */
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

  function formatParticleCount(n) {
    var v = Number(n);
    if (!isFinite(v)) v = 0;
    return v.toLocaleString("en-CA") + " pts";
  }

  function discreteDisplayLabel(key, rawCol, displayMap) {
    if (key === "__na__") return "(missing)";
    var disp = (displayMap && rawCol) ? String(displayMap[rawCol] || rawCol) : "";
    var rc = String(rawCol || "");
    var isKmeans = /k[\s_-]*means/i.test(rc) || /k[\s_-]*means/i.test(disp) || /^labels$/i.test(rc);
    if (isKmeans) return "kmeans=" + String(key);
    return String(key);
  }

  function sortDiscreteKeys(keys) {
    var rest = keys.filter(function (k) { return k !== "__na__"; });
    var allNum = rest.length > 0 && rest.every(function (k) {
      var n = Number(k);
      return isFinite(n) && String(n) === String(k).trim();
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

  function CryoColorCovariateLegend(opts) {
    this.legendContextUrl = opts.legendContextUrl;
    this.covariateDisplayMap = opts.covariateDisplayMap || {};
    this.getColorColumn = opts.getColorColumn;
    this.getPaletteName = opts.getPaletteName;
    this.panel = opts.panel;
    this.continuousWrap = opts.continuousWrap;
    this.discreteWrap = opts.discreteWrap;
    this.histDiv = opts.histDiv;
    this.thresholdUseMax = opts.thresholdUseMax;
    this.thresholdStatus = opts.thresholdStatus;
    this.discreteSwitches = opts.discreteSwitches;
    this.invertBtn = opts.invertBtn;
    this.onFilterChange = opts.onFilterChange || function () {};
    this.onModeChange = opts.onModeChange || function () {};
    this.onContinuousHistogramLayout = opts.onContinuousHistogramLayout || null;
    this.onDiscreteLayout = opts.onDiscreteLayout || null;
    this.vertical = opts.vertical !== false;
    this.histPlotVertical = opts.histPlotVertical === true;
    this._mode = null;
    this._allCatKeys = [];
    this._thresholdLevel = null;
    this._thresholdCol = null;
    this._histWired = false;
    this._suppressDiscrete = false;
    this._hoverRaf = null;
    this._hoverLevel = null;
    this._continuousValues = [];
    this._histValueRange = null;
    if (this.panel && this.vertical) {
      this.panel.classList.add("cryo-cc-legend--vertical");
    }
    var self = this;
    if (this.thresholdUseMax) {
      this.thresholdUseMax.addEventListener("change", function () {
        if (self._thresholdLevel != null && self.thresholdStatus) {
          var useMax = !!(self.thresholdUseMax && self.thresholdUseMax.checked);
          var op = useMax ? "\u2264" : "\u2265";
          self.thresholdStatus.textContent = "Threshold: " + op + " " + formatThresholdValue(self._thresholdLevel);
        }
        self._notify();
        self._relayoutHistGuides();
      });
    }
    if (this.invertBtn) {
      this.invertBtn.addEventListener("click", function () { self.invertDiscrete(); });
    }
  }

  CryoColorCovariateLegend.prototype._palette = function () {
    var n = this.getPaletteName ? this.getPaletteName() : "Viridis";
    return PALETTE_RGB[n] ? n : "Viridis";
  };

  CryoColorCovariateLegend.prototype.hide = function () {
    if (this.panel) {
      this.panel.hidden = true;
      this.panel.setAttribute("aria-hidden", "true");
      this.panel.classList.remove("cryo-cc-legend--heading-only");
    }
    this._mode = null;
    this._allCatKeys = [];
    this._thresholdLevel = null;
    this._thresholdCol = null;
    this._histValueRange = null;
  };

  /** Panel visible with only the heading (no histogram / discrete UI). Used when no colour covariate is selected. */
  CryoColorCovariateLegend.prototype.showHeadingOnly = function () {
    this._mode = null;
    this._allCatKeys = [];
    this._thresholdLevel = null;
    this._thresholdCol = null;
    this._histValueRange = null;
    this._continuousValues = [];
    this.clearThreshold();
    if (this.histDiv) {
      try {
        if (typeof Plotly !== "undefined" && Plotly && Plotly.purge) {
          Plotly.purge(this.histDiv);
        }
      } catch (e0) { /* ignore */ }
    }
    if (this.continuousWrap) {
      this.continuousWrap.hidden = true;
    }
    if (this.discreteWrap) {
      this.discreteWrap.hidden = true;
      this.discreteWrap.classList.remove("cryo-cc-discrete-wrap--show");
    }
    if (this.panel) {
      this.panel.hidden = false;
      this.panel.setAttribute("aria-hidden", "false");
      this.panel.classList.add("cryo-cc-legend--heading-only");
      this.panel.style.removeProperty("--pair-legend-middle-w");
    }
  };

  CryoColorCovariateLegend.prototype.show = function () {
    if (this.panel) {
      this.panel.hidden = false;
      this.panel.setAttribute("aria-hidden", "false");
      this.panel.classList.remove("cryo-cc-legend--heading-only");
    }
  };

  CryoColorCovariateLegend.prototype.refreshHistogramColors = function () {
    if (this._mode === "continuous" && this._continuousValues && this._continuousValues.length) {
      this._renderContinuousHistogram(this._continuousValues);
    }
  };

  CryoColorCovariateLegend.prototype._notify = function () {
    this.onFilterChange(this.getFilterForApi());
  };

  CryoColorCovariateLegend.prototype.getFilterForApi = function () {
    var col = this.getColorColumn();
    if (!col || col === "none") return null;
    if (this._mode === "continuous") {
      if (this._thresholdLevel == null || this._thresholdCol !== col) return null;
      return {
        kind: "threshold",
        level: this._thresholdLevel,
        use_max: !!(this.thresholdUseMax && this.thresholdUseMax.checked)
      };
    }
    if (this._mode === "discrete") {
      var keys = [];
      var inputs = this.discreteSwitches
        ? this.discreteSwitches.querySelectorAll("input[type=\"checkbox\"]:checked")
        : [];
      for (var i = 0; i < inputs.length; i++) {
        var k = inputs[i].getAttribute("data-cat-key");
        if (k != null) keys.push(k);
      }
      if (this._allCatKeys.length && keys.length === this._allCatKeys.length) return null;
      return { kind: "discrete", keys: keys };
    }
    return null;
  };

  CryoColorCovariateLegend.prototype.clearThreshold = function () {
    this._thresholdLevel = null;
    this._thresholdCol = null;
    if (this.thresholdStatus) this.thresholdStatus.textContent = "";
  };

  CryoColorCovariateLegend.prototype.refresh = function () {
    var self = this;
    var col = this.getColorColumn();
    if (!col || col === "none") {
      this.showHeadingOnly();
      this.onModeChange(null);
      this._notify();
      return;
    }
    fetch(this.legendContextUrl, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ column: col })
    })
      .then(function (r) {
        return r.json().then(function (j) { return { ok: r.ok, j: j }; });
      })
      .then(function (res) {
        if (!res.ok) {
          self.hide();
          self.onModeChange(null);
          return;
        }
        var data = res.j;
        self.clearThreshold();
        if (data.mode === "continuous") {
          self._mode = "continuous";
          self._continuousValues = data.values || [];
          self.show();
          if (self.continuousWrap) self.continuousWrap.hidden = false;
          if (self.discreteWrap) {
            self.discreteWrap.hidden = true;
            self.discreteWrap.classList.remove("cryo-cc-discrete-wrap--show");
          }
          self._renderContinuousHistogram(self._continuousValues);
        } else if (data.mode === "discrete") {
          self._continuousValues = [];
          self._mode = "discrete";
          self.show();
          if (self.continuousWrap) self.continuousWrap.hidden = true;
          if (self.discreteWrap) {
            self.discreteWrap.hidden = false;
            self.discreteWrap.classList.add("cryo-cc-discrete-wrap--show");
          }
          self._renderDiscrete(data.categories || []);
        } else {
          self.hide();
        }
        self.onModeChange(self._mode);
        self._notify();
      })
      .catch(function () {
        self.hide();
        self.onModeChange(null);
      });
  };

  CryoColorCovariateLegend.prototype._renderContinuousHistogram = function (values) {
    var self = this;
    if (!this.histDiv) return;
    if (!values.length) {
      Plotly.purge(this.histDiv);
      this._histValueRange = null;
      return;
    }
    var pal = this._palette();
    var violin = violinDistribution(values, pal);
    if (this.histPlotVertical) {
      violin = violinDistributionVerticalFromHorizontal(violin);
    }
    this._histValueRange = violin.range.slice();
    var layout;
    if (this.histPlotVertical) {
      layout = {
        template: "plotly_white",
        paper_bgcolor: "rgba(250,248,244,0)",
        plot_bgcolor: "rgba(255,255,255,0.62)",
        margin: { l: 4, r: 4, t: 2, b: 2 },
        xaxis: {
          visible: false,
          range: [-0.55, 0.55],
          zeroline: false,
          fixedrange: true
        },
        yaxis: {
          range: violin.range,
          showgrid: true,
          zeroline: false,
          tickfont: { size: 9, color: "#243b53" }
        },
        showlegend: false,
        hovermode: false
      };
    } else {
      layout = {
        template: "plotly_white",
        paper_bgcolor: "rgba(250,248,244,0)",
        plot_bgcolor: "rgba(255,255,255,0.62)",
        margin: { l: 4, r: 4, t: 2, b: 2 },
        xaxis: {
          range: violin.range,
          showgrid: true,
          zeroline: false,
          tickfont: { size: 9, color: "#243b53" }
        },
        yaxis: {
          visible: false,
          range: [-0.55, 0.55],
          fixedrange: true
        },
        showlegend: false,
        hovermode: false
      };
    }
    var traces = violin.fillTraces.concat([{
      type: "scatter",
      mode: "lines",
      x: violin.outlineX,
      y: violin.outlineY,
      line: { color: "#243b53", width: 1.1 },
      hoverinfo: "skip",
      showlegend: false
    }]);
    Plotly.react(this.histDiv, traces, layout, { displayModeBar: false, responsive: true });
    if (typeof this.onContinuousHistogramLayout === "function") {
      var layCb = this.onContinuousHistogramLayout;
      requestAnimationFrame(function() {
        requestAnimationFrame(function() { layCb(); });
      });
    }
    if (!this._histWired) {
      this._histWired = true;
      this.histDiv.addEventListener("click", function (evt) { self._onHistClick(evt); });
      this.histDiv.addEventListener("mousemove", function (evt) { self._onHistMove(evt); });
      this.histDiv.addEventListener("mouseleave", function () { self._onHistLeave(); });
    }
  };

  CryoColorCovariateLegend.prototype._histDataValueFromEvent = function (evt) {
    if (!this.histDiv || !this.histDiv._fullLayout) return null;
    var sz = this.histDiv._fullLayout._size;
    if (!sz) return null;
    var rect = this.histDiv.getBoundingClientRect();
    if (this.histPlotVertical) {
      var py = evt.clientY - rect.top - sz.t;
      if (py < 0 || py > sz.h) return null;
      var r = this._histValueRange;
      if (r && r.length >= 2 && r[1] !== r[0]) {
        var frac = 1 - py / sz.h;
        var level = r[0] + frac * (r[1] - r[0]);
        return isFinite(level) ? level : null;
      }
      var yaxis = this.histDiv._fullLayout.yaxis;
      if (yaxis && typeof yaxis.p2c === "function") {
        var lv = yaxis.p2c(py);
        return typeof lv === "number" && isFinite(lv) ? lv : null;
      }
      return null;
    }
    var xaxis = this.histDiv._fullLayout.xaxis;
    if (!xaxis || typeof xaxis.p2c !== "function") return null;
    var px = evt.clientX - rect.left - sz.l;
    if (px < 0 || px > sz.w) return null;
    var level = xaxis.p2c(px);
    return typeof level === "number" && isFinite(level) ? level : null;
  };

  CryoColorCovariateLegend.prototype._onHistMove = function (evt) {
    var self = this;
    var lv = this._histDataValueFromEvent(evt);
    this._hoverLevel = lv;
    if (this._hoverRaf != null) return;
    this._hoverRaf = requestAnimationFrame(function () {
      self._hoverRaf = null;
      self._relayoutHistGuides();
    });
  };

  CryoColorCovariateLegend.prototype._onHistLeave = function () {
    var self = this;
    this._hoverLevel = null;
    if (this._hoverRaf != null) return;
    this._hoverRaf = requestAnimationFrame(function () {
      self._hoverRaf = null;
      self._relayoutHistGuides();
    });
  };

  CryoColorCovariateLegend.prototype._relayoutHistGuides = function () {
    if (!this.histDiv || !this.histDiv.data) return;
    var shapes = [];
    var annotations = [];
    if (this.histPlotVertical) {
      if (this._thresholdLevel != null && this._thresholdCol === this.getColorColumn()) {
        shapes.push({
          type: "line",
          xref: "paper",
          yref: "y",
          x0: 0,
          x1: 1,
          y0: this._thresholdLevel,
          y1: this._thresholdLevel,
          line: { color: "#c4703a", width: 2 }
        });
      } else if (this._hoverLevel != null) {
        shapes.push({
          type: "line",
          xref: "paper",
          yref: "y",
          x0: 0,
          x1: 1,
          y0: this._hoverLevel,
          y1: this._hoverLevel,
          line: { color: "#243b53", width: 1.5, dash: "dot" }
        });
      }
      if (this._hoverLevel != null) {
        annotations.push({
          xref: "paper",
          yref: "y",
          x: 1,
          xanchor: "left",
          y: this._hoverLevel,
          yanchor: "middle",
          text: formatThresholdValue(this._hoverLevel),
          showarrow: false,
          bgcolor: "rgba(255,255,255,0.94)",
          bordercolor: "rgba(36,59,83,0.35)",
          borderwidth: 1,
          borderpad: 3,
          font: { size: 10, color: "#243b53" }
        });
      }
    } else {
      if (this._thresholdLevel != null && this._thresholdCol === this.getColorColumn()) {
        shapes.push({
          type: "line",
          xref: "x",
          yref: "paper",
          x0: this._thresholdLevel,
          x1: this._thresholdLevel,
          y0: 0,
          y1: 1,
          line: { color: "#c4703a", width: 2 }
        });
      } else if (this._hoverLevel != null) {
        shapes.push({
          type: "line",
          xref: "x",
          yref: "paper",
          x0: this._hoverLevel,
          x1: this._hoverLevel,
          y0: 0,
          y1: 1,
          line: { color: "#243b53", width: 1.5, dash: "dot" }
        });
      }
      if (this._hoverLevel != null) {
        annotations.push({
          xref: "x",
          yref: "paper",
          x: this._hoverLevel,
          y: 0,
          yanchor: "top",
          text: formatThresholdValue(this._hoverLevel),
          showarrow: false,
          bgcolor: "rgba(255,255,255,0.94)",
          bordercolor: "rgba(36,59,83,0.35)",
          borderwidth: 1,
          borderpad: 3,
          font: { size: 10, color: "#243b53" }
        });
      }
    }
    try {
      Plotly.relayout(this.histDiv, { shapes: shapes, annotations: annotations });
    } catch (e) { /* ignore */ }
  };

  CryoColorCovariateLegend.prototype._onHistClick = function (evt) {
    evt.preventDefault();
    var level = this._histDataValueFromEvent(evt);
    if (level == null) return;
    this._thresholdLevel = level;
    this._thresholdCol = this.getColorColumn();
    var useMax = !!(this.thresholdUseMax && this.thresholdUseMax.checked);
    var op = useMax ? "\u2264" : "\u2265";
    if (this.thresholdStatus) {
      this.thresholdStatus.textContent = "Threshold: " + op + " " + formatThresholdValue(level);
    }
    this._relayoutHistGuides();
    this._notify();
  };

  CryoColorCovariateLegend.prototype._renderDiscrete = function (categories) {
    var self = this;
    if (!this.discreteSwitches) return;
    this.discreteSwitches.innerHTML = "";
    var col = this.getColorColumn();
    var byKey = {};
    for (var j = 0; j < (categories || []).length; j++) {
      var cj = categories[j];
      byKey[String(cj.key)] = cj;
    }
    this._allCatKeys = sortDiscreteKeys(Object.keys(byKey));
    for (var i = 0; i < this._allCatKeys.length; i++) {
      var cat = byKey[this._allCatKeys[i]];
      if (!cat) continue;
      var key = String(cat.key);
      var safeId = "cc-disc-" + i + "-" + key.replace(/[^a-zA-Z0-9_-]/g, "_");
      var lab = document.createElement("label");
      lab.className = "cryo-cc-discrete-switch";
      var dSt = cat.color ? discreteCardStyles(String(cat.color)) : null;
      if (dSt && dSt.bg) {
        lab.style.backgroundColor = dSt.bg;
        lab.style.border = "1px solid " + dSt.borderColor;
      } else if (cat.color) {
        lab.style.borderLeft = "3px solid " + String(cat.color);
      }
      var inp = document.createElement("input");
      inp.type = "checkbox";
      inp.id = safeId;
      inp.setAttribute("role", "switch");
      inp.setAttribute("aria-checked", "true");
      inp.setAttribute("data-cat-key", key);
      inp.checked = true;
      var tw = document.createElement("span");
      tw.className = "cryo-cc-discrete-switch-text";
      var sp = document.createElement("span");
      sp.className = "cryo-cc-discrete-switch-label";
      sp.textContent = discreteDisplayLabel(key, col, this.covariateDisplayMap);
      if (dSt && dSt.labelColor) sp.style.color = dSt.labelColor;
      var cn = document.createElement("span");
      cn.className = "cryo-cc-discrete-switch-count";
      cn.textContent = formatParticleCount(cat.count != null ? cat.count : 0);
      if (dSt && dSt.countColor) cn.style.color = dSt.countColor;
      tw.appendChild(sp);
      tw.appendChild(cn);
      lab.appendChild(inp);
      lab.appendChild(tw);
      this.discreteSwitches.appendChild(lab);
      inp.addEventListener("change", function (input) {
        return function () {
          input.setAttribute("aria-checked", input.checked ? "true" : "false");
          if (self._suppressDiscrete) return;
          self._syncInvertBtn();
          self._notify();
        };
      }(inp));
    }
    this._syncInvertBtn();
    if (typeof this.onDiscreteLayout === "function") {
      var dLay = this.onDiscreteLayout;
      requestAnimationFrame(function () {
        requestAnimationFrame(function () { dLay(); });
      });
    }
  };

  CryoColorCovariateLegend.prototype._syncInvertBtn = function () {
    if (!this.invertBtn || !this.discreteSwitches) return;
    var inputs = this.discreteSwitches.querySelectorAll("input[type=\"checkbox\"]");
    var any = false;
    for (var i = 0; i < inputs.length; i++) {
      if (inputs[i].checked) { any = true; break; }
    }
    this.invertBtn.disabled = !any;
    this.invertBtn.textContent = any ? "Invert selection" : "no selection";
  };

  CryoColorCovariateLegend.prototype.invertDiscrete = function () {
    if (!this.discreteSwitches) return;
    this._suppressDiscrete = true;
    try {
      var inputs = this.discreteSwitches.querySelectorAll("input[type=\"checkbox\"]");
      for (var i = 0; i < inputs.length; i++) {
        inputs[i].checked = !inputs[i].checked;
        inputs[i].setAttribute("aria-checked", inputs[i].checked ? "true" : "false");
      }
    } finally {
      this._suppressDiscrete = false;
    }
    this._syncInvertBtn();
    this._notify();
  };

  global.CryoColorCovariateLegend = CryoColorCovariateLegend;
})(typeof window !== "undefined" ? window : globalThis);
