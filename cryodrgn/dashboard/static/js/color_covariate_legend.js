/**
 * Colour covariate histogram + discrete toggles (particle explorer, pair-grid, latent 3-D).
 * Continuous legends support click thresholds and drag ranges; optional ``histPlotVertical`` swaps
 * covariate onto y and density onto x. Discrete chip + popover styling lives in
 * ``static/css/cryo_cc_legend_palette_menu.css`` alongside the vertical palette menu.
 * Discrete colour picking uses one anchored panel (RGB sliders + hex + Apply / Cancel); Apply
 * commits but leaves the panel open. Plastic toggle chips default on (``discreteChipPlasticLabel``
 * opt-out). Optional ``histVerticalMargins`` shallow-merges into vertical histogram margins.
 * Panel ``vertical`` toggles the ``cryo-cc-legend--vertical`` layout class.
 * Pure helpers live in ``cryo_cc_legend_primitives.js`` (load first).
 */
(function (global) {
  "use strict";

  var CCL = global.CryoCcLegendPrimitives;
  if (!CCL) {
    throw new Error("cryo_cc_legend_primitives.js must load before color_covariate_legend.js");
  }
  var PALETTE_RGB = CCL.PALETTE_RGB;
  var paletteScaleCSS = CCL.paletteScaleCSS;
  var quantileSorted = CCL.quantileSorted;
  var violinDistribution = CCL.violinDistribution;
  var violinDistributionVerticalFromHorizontal = CCL.violinDistributionVerticalFromHorizontal;
  var formatThresholdValue = CCL.formatThresholdValue;
  var rangeInequalityPhrase = CCL.rangeInequalityPhrase;
  var parseMarkerHexRgb = CCL.parseMarkerHexRgb;
  var rgbChannelByte = CCL.rgbChannelByte;
  var rgbToHex6 = CCL.rgbToHex6;
  var discreteCardStyles = CCL.discreteCardStyles;
  var discretePlotMatchedChipStyles = CCL.discretePlotMatchedChipStyles;
  var normalizeDiscreteLegendHex = CCL.normalizeDiscreteLegendHex;
  var discreteLegendRgbTuple = CCL.discreteLegendRgbTuple;
  var cryoDiscreteColorWheelSvg = CCL.cryoDiscreteColorWheelSvg;
  var formatDiscreteDatasetCount = CCL.formatDiscreteDatasetCount;
  var discreteDisplayLabel = CCL.discreteDisplayLabel;
  var sortDiscreteKeys = CCL.sortDiscreteKeys;
  var discreteCategoryInputsFrom = CCL.discreteCategoryInputsFrom;
  var fitInvertCellTypography = CCL.fitInvertCellTypography;

  /** Shared Plotly annotation chrome for histogram threshold / range callouts. */
  function histGuideAnnotation(vertical, spec) {
    var o = {
      showarrow: false,
      bgcolor: "rgba(255,255,255,0.94)",
      bordercolor: "rgba(36,59,83,0.35)",
      borderwidth: 1,
      borderpad: 3,
      font: { size: 10, color: "#243b53" },
      text: spec.text
    };
    if (vertical) {
      o.xref = "paper";
      o.yref = "y";
      o.x = 1;
      o.xanchor = "right";
      o.y = spec.yv;
      o.yanchor = "middle";
    } else {
      o.xref = "x";
      o.yref = "paper";
      o.x = spec.xh;
      o.y = 0;
      o.yanchor = "top";
    }
    return o;
  }

  function CryoColorCovariateLegend(opts) {
    this.legendContextUrl = opts.legendContextUrl;
    this.getLegendContextData = opts.getLegendContextData || null;
    this.legendContextBodyExtra = opts.legendContextBodyExtra || null;
    this.covariateDisplayMap = opts.covariateDisplayMap || {};
    this.getColorColumn = opts.getColorColumn;
    this.getPaletteName = opts.getPaletteName;
    this.panel = opts.panel;
    this.continuousWrap = opts.continuousWrap;
    this.discreteWrap = opts.discreteWrap;
    this.histDiv = opts.histDiv;
    this.thresholdUseMax = opts.thresholdUseMax;
    this.thresholdModeCaption = opts.thresholdModeCaption || null;
    this.thresholdStatus = opts.thresholdStatus;
    this.discreteSwitches = opts.discreteSwitches;
    this.invertHeadingSlotEl = opts.invertHeadingSlotEl || null;
    this.invertHeadingRowEl = opts.invertHeadingRowEl || null;
    this._discretePanelCollapsed = false;
    this._discreteCollapseToggleBtn = null;
    /**
     * When true, category cells wrap in a grid and invert sits in a dedicated right column
     * (legacy). Ignored when invertHeadingSlotEl is set.
     */
    this.discreteInvertAsideColumn =
      opts.discreteInvertAsideColumn === true && !this.invertHeadingSlotEl;
    /**
     * When true, invert sits on its own final row spanning the toggle grid (legacy). Ignored when
     * invertHeadingSlotEl is set or discreteInvertAsideColumn is true.
     */
    this.discreteInvertFooterRow =
      opts.discreteInvertFooterRow === true
      && !this.discreteInvertAsideColumn
      && !this.invertHeadingSlotEl;
    this.invertBtn = opts.invertBtn;
    this.discreteCheckedDefault = opts.discreteCheckedDefault !== false;
    this.allDiscreteSelectedIsNull = opts.allDiscreteSelectedIsNull !== false;
    this.noDiscreteSelectionText = opts.noDiscreteSelectionText || "no selection";
    this.thresholdEmptyStatusText = opts.thresholdEmptyStatusText || "";
    this.regionHeadingEl = opts.regionHeadingEl || null;
    this.regionHeadingContinuousText = opts.regionHeadingContinuousText || "Range selection";
    this.regionHeadingDiscreteText = opts.regionHeadingDiscreteText || "Toggle selection";
    this.formatDiscreteDatasetCount =
      opts.formatDiscreteDatasetCount || formatDiscreteDatasetCount;
    this.onFilterChange = opts.onFilterChange || function () {};
    this.onModeChange = opts.onModeChange || function () {};
    this.notifyOnRefresh = opts.notifyOnRefresh !== false;
    this.onContinuousHistogramLayout = opts.onContinuousHistogramLayout || null;
    this.onDiscreteLayout = opts.onDiscreteLayout || null;
    this.vertical = opts.vertical !== false;
    this.histPlotVertical = opts.histPlotVertical === true;
    this.histVerticalMargins = opts.histVerticalMargins || null;
    this.histTickFontSize = opts.histTickFontSize != null ? Number(opts.histTickFontSize) : 9;
    this.getDiscreteColorOverrides = opts.getDiscreteColorOverrides || null;
    this.enableDiscreteColorPicker = opts.enableDiscreteColorPicker !== false;
    this.onDiscreteColorChange = opts.onDiscreteColorChange || null;
    /** When set (3-D viewer), hide the legend panel entirely if no colour covariate — skip heading-only placeholder. */
    this.hideLegendPanelWhenNoCovariate = opts.hideLegendPanelWhenNoCovariate === true;
    // Some 3-D-ish interfaces don't pass this flag explicitly; detect by panel id to keep the UI empty until
    // a colour covariate is chosen.
    if (
      !this.hideLegendPanelWhenNoCovariate
      && typeof opts.hideLegendPanelWhenNoCovariate === "undefined"
      && this.panel
      && (this.panel.id === "latent3d-color-legend-panel" || this.panel.id === "color-hist-panel")
    ) {
      this.hideLegendPanelWhenNoCovariate = true;
    }
    this.discreteChipMatplotlibBlend = opts.discreteChipMatplotlibBlend === true;
    this.discreteChipPlasticLabel = opts.discreteChipPlasticLabel !== false;
    this._discreteChipBlendAlpha = null;
    this._discreteChipBlendBg = null;
    this._discreteColorPopover = null;
    this._discreteColorPopoverState = null;
    this._docDiscretePopoverCloser = null;
    this._discretePopoverEscape = null;
    this._mode = null;
    this._allCatKeys = [];
    this._thresholdLevel = null;
    this._thresholdRange = null;
    this._thresholdCol = null;
    this._histWired = false;
    this._suppressDiscrete = false;
    this._hoverRaf = null;
    this._hoverLevel = null;
    this._histPointerDown = false;
    this._histDragMoved = false;
    this._histDragStartLevel = null;
    this._histDragStartClientCoord = null;
    this._histDragPreview = null;
    this._histSuppressClick = false;
    this._thresholdModeWasRange = false;
    this._continuousValues = [];
    this._histValueRange = null;
    this._invertTypographyFitScheduled = false;
    if (this._discreteColorPopover && !this._discreteColorPopover.hidden) {
      this._closeDiscreteColorPopover(true);
    }
    if (this.panel && this.vertical) {
      this.panel.classList.add("cryo-cc-legend--vertical");
    }
    var self = this;
    if (this.thresholdUseMax) {
      this.thresholdUseMax.addEventListener("change", function () {
        self._syncThresholdModeCaption();
        self._notify();
        self._relayoutHistGuides();
      });
    }
    if (this.invertBtn) {
      this.invertBtn.addEventListener("click", function () { self.invertDiscrete(); });
    }
    if (this.discreteSwitches) {
      this.discreteSwitches.classList.toggle(
        "cryo-color-discrete-switches--invert-aside",
        !!this.discreteInvertAsideColumn
      );
      this.discreteSwitches.classList.toggle(
        "cryo-color-discrete-switches--invert-footer-row",
        !!this.discreteInvertFooterRow
      );
    }
  }

  CryoColorCovariateLegend.prototype._getInvertCell = function () {
    if (!this.discreteSwitches) return null;
    if (this.invertHeadingSlotEl) {
      var inSlot =
        this.invertHeadingSlotEl.querySelector('[data-invert-cell="true"]')
        || this.invertHeadingSlotEl.querySelector(".cryo-cc-discrete-cell--invert");
      if (inSlot) return inSlot;
    }
    return this.discreteSwitches.querySelector('[data-invert-cell="true"]')
      || this.discreteSwitches.querySelector(".cryo-cc-discrete-cell--invert");
  };

  CryoColorCovariateLegend.prototype._clearInvertHeadingChrome = function () {
    if (!this.invertHeadingSlotEl) return;
    this.invertHeadingSlotEl.innerHTML = "";
    this.invertHeadingSlotEl.hidden = true;
    this.invertHeadingSlotEl.setAttribute("aria-hidden", "true");
    if (this.invertHeadingRowEl) {
      this.invertHeadingRowEl.hidden = true;
    }
  };

  CryoColorCovariateLegend.prototype._palette = function () {
    var n = this.getPaletteName ? this.getPaletteName() : "Viridis";
    return PALETTE_RGB[n] ? n : "Viridis";
  };

  CryoColorCovariateLegend.prototype.hide = function () {
    if (this._discreteColorPopover && !this._discreteColorPopover.hidden) {
      this._closeDiscreteColorPopover(true);
    }
    if (this.panel) {
      this.panel.hidden = true;
      this.panel.setAttribute("aria-hidden", "true");
      this.panel.classList.remove("cryo-cc-legend--heading-only");
    }
    this._mode = null;
    this._allCatKeys = [];
    this._thresholdLevel = null;
    this._thresholdRange = null;
    this._thresholdCol = null;
    this._histValueRange = null;
    this._histDragPreview = null;
    this._histPointerDown = false;
    this._discreteChipBlendAlpha = null;
    this._discreteChipBlendBg = null;
    this._syncThresholdModeCaption();
    if (this.discreteSwitches) {
      this.discreteSwitches.innerHTML = "";
    }
    this._clearInvertHeadingChrome();
    this._updateRegionHeading();
  };

  /** Panel visible with only the heading (no histogram / discrete UI). Used when no colour covariate is selected. */
  CryoColorCovariateLegend.prototype.showHeadingOnly = function () {
    if (this._discreteColorPopover && !this._discreteColorPopover.hidden) {
      this._closeDiscreteColorPopover(true);
    }
    this._mode = null;
    this._allCatKeys = [];
    this._thresholdLevel = null;
    this._thresholdRange = null;
    this._thresholdCol = null;
    this._histValueRange = null;
    this._histDragPreview = null;
    this._histPointerDown = false;
    this._continuousValues = [];
    this._discreteChipBlendAlpha = null;
    this._discreteChipBlendBg = null;
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
      this.discreteWrap.classList.remove("cryo-cc-discrete-panel--collapsed");
      this.discreteWrap.classList.remove("cryo-discrete-legend--collapsed");
    }
    if (this.discreteSwitches) {
      this.discreteSwitches.innerHTML = "";
    }
    this._discretePanelCollapsed = false;
    this._clearInvertHeadingChrome();
    if (this.panel) {
      this.panel.hidden = false;
      this.panel.setAttribute("aria-hidden", "false");
      this.panel.classList.add("cryo-cc-legend--heading-only");
      this.panel.style.removeProperty("--pair-legend-middle-w");
    }
    this._updateRegionHeading();
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
      if (this._thresholdCol !== col) return null;
      if (this._thresholdRange != null) {
        return {
          kind: "range",
          range_min: Math.min(this._thresholdRange.lo, this._thresholdRange.hi),
          range_max: Math.max(this._thresholdRange.lo, this._thresholdRange.hi),
          invert_range: !!(this.thresholdUseMax && this.thresholdUseMax.checked)
        };
      }
      if (this._thresholdLevel == null) return null;
      return {
        kind: "threshold",
        level: this._thresholdLevel,
        use_max: !!(this.thresholdUseMax && this.thresholdUseMax.checked)
      };
    }
    if (this._mode === "discrete") {
      var keys = [];
      var inputs = discreteCategoryInputsFrom(this.discreteSwitches, true);
      for (var i = 0; i < inputs.length; i++) {
        var k = inputs[i].getAttribute("data-cat-key");
        if (k != null) keys.push(k);
      }
      if (
        this.allDiscreteSelectedIsNull
        && this._allCatKeys.length
        && keys.length === this._allCatKeys.length
      ) return null;
      return { kind: "discrete", keys: keys };
    }
    return null;
  };

  CryoColorCovariateLegend.prototype.clearThreshold = function () {
    if (this._hoverRaf != null) {
      cancelAnimationFrame(this._hoverRaf);
      this._hoverRaf = null;
    }
    this._hoverLevel = null;
    this._thresholdLevel = null;
    this._thresholdRange = null;
    this._thresholdCol = null;
    this._histDragPreview = null;
    this._histPointerDown = false;
    if (this.thresholdUseMax) {
      this.thresholdUseMax.checked = false;
    }
    if (this._boundHistWindowMouseMove) {
      window.removeEventListener("mousemove", this._boundHistWindowMouseMove);
    }
    if (this._boundEndHistDrag) {
      window.removeEventListener("mouseup", this._boundEndHistDrag, true);
    }
    if (this.thresholdStatus) this.thresholdStatus.textContent = "";
    this._syncThresholdModeCaption();
  };

  /**
   * Reapply threshold/range selection onto the legend histogram after host-driven refresh
   * (e.g. particle explorer keeps parallel colour-filter state).
   */
  CryoColorCovariateLegend.prototype.syncThresholdMirrorFromHost = function (host) {
    if (!host || this._mode !== "continuous" || !this.histDiv) return;
    var col = this.getColorColumn ? this.getColorColumn() : null;
    if (!col || host.thresholdCol == null || String(host.thresholdCol) !== String(col)) return;
    var hasRange =
      host.range != null &&
      typeof host.range.lo === "number" &&
      typeof host.range.hi === "number" &&
      isFinite(host.range.lo) &&
      isFinite(host.range.hi);
    var hasLevel = host.level != null && isFinite(Number(host.level));
    if (!hasRange && !hasLevel) return;
    if (hasRange) {
      this._thresholdLevel = null;
      this._thresholdRange = {
        lo: Math.min(host.range.lo, host.range.hi),
        hi: Math.max(host.range.lo, host.range.hi)
      };
    } else {
      this._thresholdRange = null;
      this._thresholdLevel = Number(host.level);
    }
    this._thresholdCol = col;
    // Align "previous mode" with the new threshold/range before touching the checkbox so
    // `_syncThresholdModeCaption` does not treat host-driven updates as a level↔range flip and
    // force-clear `thresholdUseMax` (breaking invert / ≤ highlighting on first mirror sync).
    this._thresholdModeWasRange = this._rangeModeActive();
    if (this.thresholdUseMax && typeof host.useMaxChecked === "boolean") {
      this.thresholdUseMax.checked = !!host.useMaxChecked;
    }
    this._syncThresholdModeCaption();
    this._relayoutHistGuides();
  };

  CryoColorCovariateLegend.prototype._discreteToggleHeadingRow = function () {
    if (this.invertHeadingRowEl) return this.invertHeadingRowEl;
    if (!this.discreteWrap) return null;
    return this.discreteWrap.querySelector(".cryo-cc-discrete-toggle-heading-row");
  };

  CryoColorCovariateLegend.prototype._upgradeDiscreteCollapseTitleButton = function (row) {
    if (!row) return null;
    var btn = row.querySelector("[data-cryo-discrete-collapse-toggle]");
    var titleSpan = row.querySelector(".cryo-cc-discrete-toggle-heading-text");
    var label = this.regionHeadingDiscreteText || "Toggle selection";
    if (btn && btn.classList.contains("cryo-cc-discrete-toggle-title")) {
      if (titleSpan && titleSpan.textContent) {
        label = titleSpan.textContent;
      }
      return btn;
    }
    var chev = btn && btn.querySelector(".cryo-cc-discrete-collapse-chevron");
    var expanded = btn ? btn.getAttribute("aria-expanded") : "true";
    var controls = btn ? btn.getAttribute("aria-controls") : null;
    if (!chev) {
      chev = document.createElement("span");
      chev.className = "cryo-cc-discrete-collapse-chevron";
      chev.setAttribute("aria-hidden", "true");
      chev.textContent = "\u25bc";
    }
    var titleBtn = document.createElement("button");
    titleBtn.type = "button";
    titleBtn.className = "cryo-cc-discrete-toggle-title";
    titleBtn.setAttribute("data-cryo-discrete-collapse-toggle", "true");
    titleBtn.setAttribute("aria-expanded", expanded || "true");
    if (controls) titleBtn.setAttribute("aria-controls", controls);
    else if (this.discreteSwitches && this.discreteSwitches.id) {
      titleBtn.setAttribute("aria-controls", this.discreteSwitches.id);
    }
    var titleInner = document.createElement("span");
    titleInner.className = "cryo-cc-discrete-toggle-heading-text";
    titleInner.textContent = titleSpan
      ? titleSpan.textContent
      : label;
    titleBtn.appendChild(chev);
    titleBtn.appendChild(titleInner);
    if (btn) row.replaceChild(titleBtn, btn);
    else row.insertBefore(titleBtn, row.firstChild);
    if (titleSpan && titleSpan.parentNode === row) titleSpan.remove();
    return titleBtn;
  };

  CryoColorCovariateLegend.prototype._ensureDiscreteCollapseToggle = function () {
    var self = this;
    if (!this.discreteWrap || !this.discreteSwitches) return;
    var row = this._discreteToggleHeadingRow();
    if (!row) {
      row = document.createElement("div");
      row.className = "cryo-cc-discrete-toggle-heading-row";
      this.discreteWrap.insertBefore(row, this.discreteSwitches);
      this.invertHeadingRowEl = row;
    }
    row.hidden = false;
    var btn = this._upgradeDiscreteCollapseTitleButton(row);
    if (!btn) return;
    this._discreteCollapseToggleBtn = btn;
    btn.setAttribute("aria-expanded", this._discretePanelCollapsed ? "false" : "true");
    btn.title = this._discretePanelCollapsed
      ? "Show toggle selection"
      : "Hide toggle selection";
    if (this.discreteWrap && !this.discreteWrap.dataset.cryoDiscreteCollapseBound) {
      this.discreteWrap.dataset.cryoDiscreteCollapseBound = "1";
      this.discreteWrap.addEventListener("click", function (ev) {
        var toggle = ev.target.closest("[data-cryo-discrete-collapse-toggle]");
        if (!toggle || !self.discreteWrap.contains(toggle)) return;
        ev.preventDefault();
        self._discreteCollapseToggleBtn = toggle;
        self.setDiscretePanelCollapsed(!self._discretePanelCollapsed);
      });
    }
  };

  CryoColorCovariateLegend.prototype.setDiscretePanelCollapsed = function (collapsed) {
    this._discretePanelCollapsed = !!collapsed;
    if (!this.discreteWrap) return;
    this.discreteWrap.classList.toggle("cryo-cc-discrete-panel--collapsed", this._discretePanelCollapsed);
    this.discreteWrap.classList.toggle("cryo-discrete-legend--collapsed", this._discretePanelCollapsed);
    if (this._discreteCollapseToggleBtn) {
      this._discreteCollapseToggleBtn.setAttribute(
        "aria-expanded",
        this._discretePanelCollapsed ? "false" : "true"
      );
      this._discreteCollapseToggleBtn.title = this._discretePanelCollapsed
        ? "Show toggle selection"
        : "Hide toggle selection";
    }
    this._syncInvertBtn();
    if (typeof this.onDiscretePanelCollapsedChange === "function") {
      this.onDiscretePanelCollapsedChange(this._discretePanelCollapsed);
    }
    if (typeof this.onDiscreteLayout === "function") {
      this.onDiscreteLayout();
    }
  };

  CryoColorCovariateLegend.prototype.isDiscretePanelCollapsed = function () {
    return !!this._discretePanelCollapsed;
  };

  CryoColorCovariateLegend.prototype._updateRegionHeading = function () {
    var innerDiscreteHeading = this._mode === "discrete" && !!this._discreteToggleHeadingRow();
    if (this.regionHeadingEl) {
      if (this._mode === "continuous") {
        this.regionHeadingEl.hidden = false;
        this.regionHeadingEl.textContent = this.regionHeadingContinuousText;
      } else if (this._mode === "discrete" && !innerDiscreteHeading) {
        this.regionHeadingEl.hidden = false;
        this.regionHeadingEl.textContent = this.regionHeadingDiscreteText;
      } else {
        this.regionHeadingEl.hidden = true;
        this.regionHeadingEl.textContent = "";
      }
    }
  };

  CryoColorCovariateLegend.prototype._activeColorLabel = function () {
    var col = this.getColorColumn ? this.getColorColumn() : "";
    return (this.covariateDisplayMap && col) ? String(this.covariateDisplayMap[col] || col) : String(col || "value");
  };

  CryoColorCovariateLegend.prototype._rangeModeActive = function () {
    return this._thresholdRange != null && this._thresholdCol === this.getColorColumn();
  };

  CryoColorCovariateLegend.prototype._syncThresholdModeCaption = function () {
    var rangeMode = this._rangeModeActive();
    if (this.thresholdModeCaption) {
      this.thresholdModeCaption.textContent = rangeMode ? "Invert range selection" : "Use as maximum (\u2264)";
    }
    if (this._thresholdModeWasRange !== rangeMode && this.thresholdUseMax) {
      this.thresholdUseMax.checked = false;
    }
    this._thresholdModeWasRange = rangeMode;
    if (!this.thresholdStatus) return;
    if (this._thresholdCol !== this.getColorColumn()) {
      this.thresholdStatus.textContent = "";
      return;
    }
    if (this._thresholdRange != null) {
      var phrase = rangeInequalityPhrase(this._thresholdRange.lo, this._thresholdRange.hi, this._activeColorLabel());
      this.thresholdStatus.textContent = (this.thresholdUseMax && this.thresholdUseMax.checked)
        ? "Outside " + phrase
        : "Range: " + phrase;
      return;
    }
    if (this._thresholdLevel != null) {
      var useMax = !!(this.thresholdUseMax && this.thresholdUseMax.checked);
      var op = useMax ? "\u2264" : "\u2265";
      this.thresholdStatus.textContent = "Threshold: " + op + " " + formatThresholdValue(this._thresholdLevel);
      return;
    }
    this.thresholdStatus.textContent = this.thresholdEmptyStatusText;
  };

  CryoColorCovariateLegend.prototype.refresh = function (extra) {
    var ex = extra && typeof extra === "object" ? extra : {};
    var suppressNotify = !!ex.suppressNotify;
    var afterLayout = ex.afterLayout;
    var self = this;
    var col = this.getColorColumn();
    if (!col || col === "none") {
      if (this.hideLegendPanelWhenNoCovariate) {
        this.hide();
      } else {
        this.showHeadingOnly();
      }
      this.onModeChange(null);
      if (this.notifyOnRefresh && !suppressNotify) this._notify();
      if (typeof afterLayout === "function") {
        try {
          afterLayout(null);
        } catch (e0) { /* ignore */ }
      }
      return Promise.resolve({ ok: true, j: {} });
    }
    var contextPromise = this.getLegendContextData
      ? Promise.resolve().then(function () {
          return { ok: true, j: self.getLegendContextData(col) || {} };
        })
      : fetch(this.legendContextUrl, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify(
            (function () {
              var body = { column: col };
              var bodyExtra = self.legendContextBodyExtra;
              if (bodyExtra && typeof bodyExtra === "object") {
                Object.keys(bodyExtra).forEach(function (k) {
                  body[k] = bodyExtra[k];
                });
              }
              return body;
            })()
          )
        }).then(function (r) {
          return r.json().then(function (j) { return { ok: r.ok, j: j }; });
        });
    return contextPromise
      .then(function (res) {
        if (!res.ok) {
          self.hide();
          self.onModeChange(null);
          return res;
        }
        var data = res.j;
        var preservedContinuousThreshold = null;
        if (
          self._mode === "continuous"
          && self._thresholdCol === col
          && (self._thresholdLevel != null || self._thresholdRange != null)
        ) {
          preservedContinuousThreshold = {
            level: self._thresholdLevel,
            range: self._thresholdRange,
            thresholdCol: self._thresholdCol,
            useMaxChecked: !!(self.thresholdUseMax && self.thresholdUseMax.checked)
          };
        }
        self.clearThreshold();
        if (data.mode === "continuous") {
          if (self._discreteColorPopover && !self._discreteColorPopover.hidden) {
            self._closeDiscreteColorPopover(true);
          }
          self._mode = "continuous";
          self._continuousValues = data.values || [];
          self._discreteChipBlendAlpha = null;
          self._discreteChipBlendBg = null;
          self.show();
          if (self.continuousWrap) self.continuousWrap.hidden = false;
          if (self.discreteWrap) {
            self.discreteWrap.hidden = true;
            self.discreteWrap.classList.remove("cryo-cc-discrete-wrap--show");
            self.discreteWrap.classList.remove("cryo-cc-discrete-panel--collapsed");
            self.discreteWrap.classList.remove("cryo-discrete-legend--collapsed");
            self._discretePanelCollapsed = false;
          }
          if (self.discreteSwitches) {
            self.discreteSwitches.innerHTML = "";
          }
          self._clearInvertHeadingChrome();
          self._allCatKeys = [];
          if (
            preservedContinuousThreshold
            && preservedContinuousThreshold.thresholdCol === col
          ) {
            self._thresholdLevel = preservedContinuousThreshold.level;
            self._thresholdRange = preservedContinuousThreshold.range;
            self._thresholdCol = preservedContinuousThreshold.thresholdCol;
            if (self.thresholdUseMax) {
              self.thresholdUseMax.checked = preservedContinuousThreshold.useMaxChecked;
            }
            self._thresholdModeWasRange = preservedContinuousThreshold.range != null;
          }
          return self._renderContinuousHistogram(self._continuousValues).then(function () {
            finishLegendLayout();
            return res;
          });
        }
        if (data.mode === "discrete") {
          self._continuousValues = [];
          self._discreteChipBlendAlpha =
            typeof data.chip_blend_alpha === "number" && isFinite(data.chip_blend_alpha)
              ? data.chip_blend_alpha
              : null;
          self._discreteChipBlendBg =
            typeof data.chip_blend_bg === "string" ? data.chip_blend_bg : null;
          self._mode = "discrete";
          self.show();
          if (self.continuousWrap) self.continuousWrap.hidden = true;
          if (self.discreteWrap) {
            self.discreteWrap.hidden = false;
            self.discreteWrap.classList.add("cryo-cc-discrete-wrap");
            self.discreteWrap.classList.add("cryo-cc-discrete-wrap--show");
          }
          self._renderDiscrete(data.categories || []);
          try {
            if (self.histDiv && self.histDiv._fullLayout) {
              Plotly.relayout(self.histDiv, { shapes: [], annotations: [] });
            }
          } catch (eHist) { /* ignore */ }
        } else {
          self.hide();
        }
        finishLegendLayout();
        return res;

        function finishLegendLayout() {
          self._updateRegionHeading();
          self.onModeChange(self._mode);
          if (typeof afterLayout === "function") {
            try {
              afterLayout(self._mode);
            } catch (eLay) { /* ignore */ }
          }
          if (self.notifyOnRefresh && !suppressNotify) self._notify();
        }
      })
      .catch(function () {
        self.hide();
        self.onModeChange(null);
      });
  };

  CryoColorCovariateLegend.prototype._renderContinuousHistogram = function (values) {
    var self = this;
    if (!this.histDiv) return Promise.resolve();
    if (!values.length) {
      Plotly.purge(this.histDiv);
      this._histValueRange = null;
      return Promise.resolve();
    }
    var pal = this._palette();
    var violin = violinDistribution(values, pal);
    if (this.histPlotVertical) {
      violin = violinDistributionVerticalFromHorizontal(violin);
    }
    this._histValueRange = violin.range.slice();
    var layout;
    if (this.histPlotVertical) {
      var vMargin = { l: 4, r: 4, t: 2, b: 2 };
      if (this.histVerticalMargins && typeof this.histVerticalMargins === "object") {
        var hv = this.histVerticalMargins;
        for (var mk in hv) {
          if (Object.prototype.hasOwnProperty.call(hv, mk)) {
            vMargin[mk] = hv[mk];
          }
        }
      }
      layout = {
        template: "plotly_white",
        paper_bgcolor: "rgba(250,248,244,0)",
        plot_bgcolor: "rgba(255,255,255,0.62)",
        margin: vMargin,
        dragmode: false,
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
          fixedrange: true,
          tickfont: { size: this.histTickFontSize, color: "#243b53" }
        },
        showlegend: false,
        hovermode: false
      };
    } else {
      layout = {
        template: "plotly_white",
        paper_bgcolor: "rgba(250,248,244,0)",
        plot_bgcolor: "rgba(255,255,255,0.62)",
        margin: { l: 4, r: 4, t: 2, b: 17 },
        dragmode: false,
        xaxis: {
          range: violin.range,
          showgrid: true,
          zeroline: false,
          fixedrange: true,
          tickfont: { size: this.histTickFontSize, color: "#243b53" }
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
    var reactMaybe = Plotly.react(this.histDiv, traces, layout, {
      displayModeBar: false,
      doubleClick: false,
      scrollZoom: false,
      responsive: true
    });
    this._syncThresholdModeCaption();
    this._relayoutHistGuides();
    requestAnimationFrame(function () {
      requestAnimationFrame(function () {
        self._relayoutHistGuides();
      });
    });
    if (typeof this.onContinuousHistogramLayout === "function") {
      var layCb = this.onContinuousHistogramLayout;
      requestAnimationFrame(function() {
        requestAnimationFrame(function() { layCb(); });
      });
    }
    if (!this._histWired) {
      this._histWired = true;
      this.histDiv.addEventListener("click", function (evt) { self._onHistClick(evt); });
      this.histDiv.addEventListener("mousedown", function (evt) { self._onHistMouseDown(evt); }, true);
      this.histDiv.addEventListener("mousemove", function (evt) { self._onHistMove(evt); });
      this.histDiv.addEventListener("mouseleave", function () { self._onHistLeave(); });
    }
    var pr = Promise.resolve(reactMaybe && typeof reactMaybe.then === "function" ? reactMaybe : null);
    return pr.catch(function () {}).then(function () {
      self._relayoutHistGuides();
      return new Promise(function (resolve) {
        requestAnimationFrame(function () {
          requestAnimationFrame(function () {
            self._relayoutHistGuides();
            resolve();
          });
        });
      });
    });
  };

  CryoColorCovariateLegend.prototype._histDataValueFromEvent = function (evt, opts) {
    opts = opts || {};
    if (!evt) return null;
    var clamp = !!opts.clamp;
    if (!this.histDiv || !this.histDiv._fullLayout) return null;
    var sz = this.histDiv._fullLayout._size;
    if (!sz) return null;
    var rect = this.histDiv.getBoundingClientRect();
    if (this.histPlotVertical) {
      var py = evt.clientY - rect.top - sz.t;
      if (clamp) {
        if (sz.h <= 0) return null;
        py = Math.max(0, Math.min(sz.h, py));
      } else if (py < 0 || py > sz.h) {
        return null;
      }
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
    if (clamp) {
      if (sz.w <= 0) return null;
      px = Math.max(0, Math.min(sz.w, px));
    } else if (px < 0 || px > sz.w) {
      return null;
    }
    var level = xaxis.p2c(px);
    return typeof level === "number" && isFinite(level) ? level : null;
  };

  CryoColorCovariateLegend.prototype._onHistMove = function (evt) {
    if (this._histPointerDown) return;
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
    if (this._histPointerDown) return;
    var self = this;
    this._hoverLevel = null;
    if (this._hoverRaf != null) return;
    this._hoverRaf = requestAnimationFrame(function () {
      self._hoverRaf = null;
      self._relayoutHistGuides();
    });
  };

  CryoColorCovariateLegend.prototype._plotValueRange = function () {
    var r = this._histValueRange;
    if (r && r.length >= 2) {
      return [Math.min(r[0], r[1]), Math.max(r[0], r[1])];
    }
    if (!this.histDiv || !this.histDiv._fullLayout) return null;
    var ax = this.histPlotVertical ? this.histDiv._fullLayout.yaxis : this.histDiv._fullLayout.xaxis;
    if (!ax || !ax.range || ax.range.length < 2) return null;
    var rng = ax.range;
    return [Math.min(rng[0], rng[1]), Math.max(rng[0], rng[1])];
  };

  CryoColorCovariateLegend.prototype._oneHighlightRect = function (range, lo, hi) {
    if (!range) return null;
    var a = Math.max(range[0], Math.min(lo, hi));
    var b = Math.min(range[1], Math.max(lo, hi));
    if (!(b > a)) return null;
    if (this.histPlotVertical) {
      return {
        type: "rect",
        xref: "paper",
        yref: "y",
        x0: 0,
        x1: 1,
        y0: a,
        y1: b,
        fillcolor: "rgba(196, 112, 58, 0.14)",
        line: { width: 0 },
        layer: "above"
      };
    }
    return {
      type: "rect",
      xref: "x",
      yref: "paper",
      x0: a,
      x1: b,
      y0: 0,
      y1: 1,
      fillcolor: "rgba(196, 112, 58, 0.14)",
      line: { width: 0 },
      layer: "above"
    };
  };

  CryoColorCovariateLegend.prototype._selectionHighlightRects = function () {
    var range = this._plotValueRange();
    if (!range) return [];
    if (this._histDragPreview != null) {
      var dr = this._oneHighlightRect(range, this._histDragPreview.lo, this._histDragPreview.hi);
      return dr ? [dr] : [];
    }
    if (this._thresholdCol !== this.getColorColumn()) return [];
    if (this._thresholdRange != null) {
      var lo = Math.min(this._thresholdRange.lo, this._thresholdRange.hi);
      var hi = Math.max(this._thresholdRange.lo, this._thresholdRange.hi);
      if (!(this.thresholdUseMax && this.thresholdUseMax.checked)) {
        var rr = this._oneHighlightRect(range, lo, hi);
        return rr ? [rr] : [];
      }
      var out = [];
      var r1 = this._oneHighlightRect(range, range[0], lo);
      var r2 = this._oneHighlightRect(range, hi, range[1]);
      if (r1) out.push(r1);
      if (r2) out.push(r2);
      return out;
    }
    if (this._thresholdLevel == null) return [];
    var useMax = !!(this.thresholdUseMax && this.thresholdUseMax.checked);
    var bx0 = useMax ? range[0] : this._thresholdLevel;
    var bx1 = useMax ? this._thresholdLevel : range[1];
    var br = this._oneHighlightRect(range, bx0, bx1);
    return br ? [br] : [];
  };

  CryoColorCovariateLegend.prototype._guideLine = function (value, opts) {
    opts = opts || {};
    var color = opts.color || "#c4703a";
    var line = { color: color, width: opts.width || 2 };
    if (opts.dash) line.dash = opts.dash;
    if (this.histPlotVertical) {
      return {
        type: "line",
        xref: "paper",
        yref: "y",
        x0: 0,
        x1: 1,
        y0: value,
        y1: value,
        line: line,
        layer: "above"
      };
    }
    return {
      type: "line",
      xref: "x",
      yref: "paper",
      x0: value,
      x1: value,
      y0: 0,
      y1: 1,
      line: line,
      layer: "above"
    };
  };

  CryoColorCovariateLegend.prototype._valueAnnotation = function (value) {
    return histGuideAnnotation(this.histPlotVertical, {
      text: formatThresholdValue(value),
      yv: value,
      xh: value
    });
  };

  CryoColorCovariateLegend.prototype._rangeAnnotation = function (lo, hi) {
    var mid = (lo + hi) * 0.5;
    return histGuideAnnotation(this.histPlotVertical, {
      text: rangeInequalityPhrase(lo, hi, this._activeColorLabel()),
      yv: mid,
      xh: mid
    });
  };

  CryoColorCovariateLegend.prototype._relayoutHistGuides = function () {
    if (!this.histDiv || !this.histDiv.data) return;
    var shapes = this._selectionHighlightRects();
    var annotations = [];
    if (this._histDragPreview != null) {
      var pLo = Math.min(this._histDragPreview.lo, this._histDragPreview.hi);
      var pHi = Math.max(this._histDragPreview.lo, this._histDragPreview.hi);
      shapes.push(this._guideLine(pLo));
      shapes.push(this._guideLine(pHi));
      annotations.push(this._rangeAnnotation(pLo, pHi));
    } else if (this._thresholdCol === this.getColorColumn()) {
      if (this._thresholdRange != null) {
        var rLo = Math.min(this._thresholdRange.lo, this._thresholdRange.hi);
        var rHi = Math.max(this._thresholdRange.lo, this._thresholdRange.hi);
        shapes.push(this._guideLine(rLo));
        shapes.push(this._guideLine(rHi));
        annotations.push(this._rangeAnnotation(rLo, rHi));
      } else if (this._thresholdLevel != null) {
        shapes.push(this._guideLine(this._thresholdLevel));
        annotations.push(this._valueAnnotation(this._thresholdLevel));
      }
    }
    if (this._hoverLevel != null && !this._histPointerDown) {
      shapes.push(this._guideLine(this._hoverLevel, { color: "#243b53", width: 1.5, dash: "dot" }));
      var colNow = this.getColorColumn ? this.getColorColumn() : null;
      var showHoverValueAnnot =
        this._histDragPreview == null
        && !(this._thresholdCol === colNow && this._thresholdRange != null)
        && !(this._thresholdCol === colNow && this._thresholdLevel != null);
      if (showHoverValueAnnot) {
        annotations.push(this._valueAnnotation(this._hoverLevel));
      }
    }
    try {
      Plotly.relayout(this.histDiv, { shapes: shapes, annotations: annotations });
    } catch (e) { /* ignore */ }
  };

  CryoColorCovariateLegend.prototype._onHistClick = function (evt) {
    evt.preventDefault();
    if (this._histSuppressClick) {
      this._histSuppressClick = false;
      return;
    }
    var level = this._histDataValueFromEvent(evt);
    if (level == null) return;
    this._thresholdRange = null;
    this._thresholdLevel = level;
    this._thresholdCol = this.getColorColumn();
    this._syncThresholdModeCaption();
    this._relayoutHistGuides();
    this._notify();
  };

  CryoColorCovariateLegend.prototype._onHistWindowMouseMove = function (evt) {
    if (!this._histPointerDown || this._histDragStartLevel == null) return;
    var level = this._histDataValueFromEvent(evt, { clamp: true });
    if (level == null) return;
    this._histDragMoved = true;
    this._histDragPreview = { lo: this._histDragStartLevel, hi: level };
    this._relayoutHistGuides();
  };

  CryoColorCovariateLegend.prototype._endHistDrag = function (evt) {
    window.removeEventListener("mousemove", this._boundHistWindowMouseMove);
    window.removeEventListener("mouseup", this._boundEndHistDrag, true);
    if (!this._histPointerDown) return;
    this._histPointerDown = false;
    var appliedRange = false;
    var pxSpan = 0;
    var coord = this.histPlotVertical ? "clientY" : "clientX";
    if (this._histDragStartClientCoord != null && evt && typeof evt[coord] === "number") {
      pxSpan = Math.abs(evt[coord] - this._histDragStartClientCoord);
    }
    var levelEnd = this._histDataValueFromEvent(evt, { clamp: true });
    if (this._histDragMoved && levelEnd != null && pxSpan >= 5) {
      this._thresholdLevel = null;
      this._thresholdRange = {
        lo: Math.min(this._histDragStartLevel, levelEnd),
        hi: Math.max(this._histDragStartLevel, levelEnd)
      };
      this._thresholdCol = this.getColorColumn();
      appliedRange = true;
      this._histSuppressClick = true;
      this._syncThresholdModeCaption();
      this._notify();
    }
    this._histDragPreview = null;
    this._histDragStartLevel = null;
    this._histDragStartClientCoord = null;
    this._histDragMoved = false;
    this._relayoutHistGuides();
    if (!appliedRange) {
      this._onHistMove(evt);
    }
  };

  CryoColorCovariateLegend.prototype._onHistMouseDown = function (evt) {
    if (evt.button !== 0) return;
    var level = this._histDataValueFromEvent(evt);
    if (level == null) return;
    this._histPointerDown = true;
    this._histDragMoved = false;
    this._histDragStartLevel = level;
    this._histDragStartClientCoord = this.histPlotVertical ? evt.clientY : evt.clientX;
    this._histDragPreview = null;
    var self = this;
    this._boundHistWindowMouseMove = this._boundHistWindowMouseMove || function (e) {
      self._onHistWindowMouseMove(e);
    };
    this._boundEndHistDrag = this._boundEndHistDrag || function (e) {
      self._endHistDrag(e);
    };
    window.addEventListener("mousemove", this._boundHistWindowMouseMove);
    window.addEventListener("mouseup", this._boundEndHistDrag, true);
    evt.preventDefault();
    evt.stopPropagation();
  };

  CryoColorCovariateLegend.prototype._applyDiscreteRowPaint = function (
    labEl,
    labelSpanEl,
    countSpanEl,
    hex
  ) {
    var hx = normalizeDiscreteLegendHex(hex);
    if (!hx) return;
    var dSt = discreteCardStyles(hx);
    labEl.style.borderLeft = "";
    labEl.style.backgroundColor = "";
    labEl.style.border = "";
    if (dSt && dSt.bg) {
      labEl.style.backgroundColor = dSt.bg;
      labEl.style.border = "1px solid " + dSt.borderColor;
    }
    if (labelSpanEl) {
      if (dSt && dSt.labelColor) labelSpanEl.style.color = dSt.labelColor;
      else labelSpanEl.style.color = "";
    }
    if (countSpanEl) {
      if (dSt && dSt.countColor) countSpanEl.style.color = dSt.countColor;
      else countSpanEl.style.color = "";
    }
  };

  CryoColorCovariateLegend.prototype._closeDiscreteColorPopover = function (asCancel) {
    var pop = this._discreteColorPopover;
    var st = this._discreteColorPopoverState;
    if (!pop || pop.hidden) return;
    if (asCancel && st) {
      var comm = normalizeDiscreteLegendHex(st.commitHex) || "#aab4bf";
      this._applyDiscreteCellPaint(st.cellEl, comm);
      this._discretePopoverSetUiFromHex(comm);
    }
    pop.hidden = true;
    this._discreteColorPopoverState = null;
  };

  CryoColorCovariateLegend.prototype._discretePopoverSetUiFromHex = function (hex) {
    var pop = this._discreteColorPopover;
    if (!pop) return;
    var hx = normalizeDiscreteLegendHex(hex) || "#aab4bf";
    var rgb = discreteLegendRgbTuple(hx);
    var preview = pop.querySelector(".cryo-cc-discrete-color-popover-preview");
    var hexIn = pop.querySelector(".cryo-cc-discrete-color-popover-hex");
    var rr = pop.querySelector('input[type=\"range\"][data-channel=\"r\"]');
    var gg = pop.querySelector('input[type=\"range\"][data-channel=\"g\"]');
    var bb = pop.querySelector('input[type=\"range\"][data-channel=\"b\"]');
    if (preview) preview.style.backgroundColor = hx;
    if (hexIn) hexIn.value = hx;
    if (rr) rr.value = String(rgb.r);
    if (gg) gg.value = String(rgb.g);
    if (bb) bb.value = String(rgb.b);
  };

  CryoColorCovariateLegend.prototype._discretePopoverReadHexFromUi = function () {
    var pop = this._discreteColorPopover;
    if (!pop) return "";
    var rr = pop.querySelector('input[type=\"range\"][data-channel=\"r\"]');
    var gg = pop.querySelector('input[type=\"range\"][data-channel=\"g\"]');
    var bb = pop.querySelector('input[type=\"range\"][data-channel=\"b\"]');
    if (!rr || !gg || !bb) return "";
    var r = rgbChannelByte(rr.value);
    var g = rgbChannelByte(gg.value);
    var b = rgbChannelByte(bb.value);
    return rgbToHex6(r, g, b);
  };

  CryoColorCovariateLegend.prototype._discretePopoverPushLivePreview = function () {
    var st = this._discreteColorPopoverState;
    var pop = this._discreteColorPopover;
    if (!st || !pop) return;
    var hx = this._discretePopoverReadHexFromUi();
    if (!hx) return;
    var preview = pop.querySelector(".cryo-cc-discrete-color-popover-preview");
    var hexIn = pop.querySelector(".cryo-cc-discrete-color-popover-hex");
    if (preview) preview.style.backgroundColor = hx;
    if (hexIn) hexIn.value = hx;
    this._applyDiscreteCellPaint(st.cellEl, hx);
  };

  CryoColorCovariateLegend.prototype._ensureDiscreteColorPopover = function () {
    var self = this;
    if (this._discreteColorPopover) return this._discreteColorPopover;
    var uid =
      "ccdp-" +
      (typeof crypto !== "undefined" && crypto.randomUUID
        ? crypto.randomUUID().slice(0, 8)
        : String(Math.random()).slice(2, 10));

    var root = document.createElement("div");
    root.className = "cryo-cc-discrete-color-popover";
    root.setAttribute("role", "dialog");
    root.setAttribute("aria-label", "Choose category colour");
    root.hidden = true;

    var title = document.createElement("div");
    title.className = "cryo-cc-discrete-color-popover-title";
    title.textContent = "Category colour";

    var preview = document.createElement("div");
    preview.className = "cryo-cc-discrete-color-popover-preview";

    var slidersWrap = document.createElement("div");
    slidersWrap.className = "cryo-cc-discrete-color-popover-sliders";

    function addSlider(channel, letter) {
      var row = document.createElement("div");
      row.className = "cryo-cc-discrete-color-popover-slider-row";
      var lab = document.createElement("label");
      lab.setAttribute("for", uid + "-" + channel);
      lab.textContent = letter;
      var rg = document.createElement("input");
      rg.type = "range";
      rg.id = uid + "-" + channel;
      rg.setAttribute("data-channel", channel);
      rg.min = "0";
      rg.max = "255";
      rg.value = "128";
      rg.addEventListener("input", function () {
        self._discretePopoverPushLivePreview();
      });
      row.appendChild(lab);
      row.appendChild(rg);
      slidersWrap.appendChild(row);
    }
    addSlider("r", "R");
    addSlider("g", "G");
    addSlider("b", "B");

    var hexRow = document.createElement("div");
    hexRow.className = "cryo-cc-discrete-color-popover-hex-row";
    var hexLab = document.createElement("label");
    hexLab.setAttribute("for", uid + "-hex");
    hexLab.textContent = "Hex";
    var hexIn = document.createElement("input");
    hexIn.type = "text";
    hexIn.className = "cryo-cc-discrete-color-popover-hex";
    hexIn.id = uid + "-hex";
    hexIn.setAttribute("spellcheck", "false");
    hexIn.setAttribute("autocomplete", "off");
    hexIn.setAttribute("aria-label", "Colour as #rrggbb");
    hexIn.addEventListener("input", function () {
      var raw = String(hexIn.value || "").trim();
      var hx = normalizeDiscreteLegendHex(raw);
      if (!hx && /^[0-9a-fA-F]{6}$/.test(raw.replace(/^#/, ""))) {
        hx = normalizeDiscreteLegendHex("#" + raw.replace(/^#/, ""));
      }
      if (!hx || hx.length !== 7) return;
      self._discretePopoverSetUiFromHex(hx);
      var st = self._discreteColorPopoverState;
      if (st) self._applyDiscreteCellPaint(st.cellEl, hx);
    });

    hexRow.appendChild(hexLab);
    hexRow.appendChild(hexIn);

    var actions = document.createElement("div");
    actions.className = "cryo-cc-discrete-color-popover-actions";
    var applyBtn = document.createElement("button");
    applyBtn.type = "button";
    applyBtn.className = "cryo-cc-discrete-color-popover-apply";
    applyBtn.textContent = "Apply";
    var cancelBtn = document.createElement("button");
    cancelBtn.type = "button";
    cancelBtn.className = "cryo-cc-discrete-color-popover-cancel";
    cancelBtn.textContent = "Cancel";

    applyBtn.addEventListener("click", function () {
      var st = self._discreteColorPopoverState;
      if (!st) return;
      var hexInEl = root.querySelector(".cryo-cc-discrete-color-popover-hex");
      var hx = hexInEl ? normalizeDiscreteLegendHex(hexInEl.value) : "";
      if (!hx) hx = normalizeDiscreteLegendHex(self._discretePopoverReadHexFromUi());
      if (!hx) return;
      self._discretePopoverSetUiFromHex(hx);
      self._applyDiscreteCellPaint(st.cellEl, hx);
      if (typeof self.onDiscreteColorChange === "function") {
        self.onDiscreteColorChange({ catKey: st.catKey, hex: hx.toLowerCase() });
      }
      st.commitHex = hx;
    });

    cancelBtn.addEventListener("click", function () {
      self._closeDiscreteColorPopover(true);
    });

    actions.appendChild(cancelBtn);
    actions.appendChild(applyBtn);

    root.appendChild(title);
    root.appendChild(preview);
    root.appendChild(slidersWrap);
    root.appendChild(hexRow);
    root.appendChild(actions);

    this._docDiscretePopoverCloser = function (ev) {
      var pop = self._discreteColorPopover;
      if (!pop || pop.hidden) return;
      var st = self._discreteColorPopoverState;
      var t = ev.target;
      if (pop.contains(t)) return;
      if (st && st.anchorEl && st.anchorEl.contains(t)) return;
      self._closeDiscreteColorPopover(true);
    };
    document.addEventListener("mousedown", this._docDiscretePopoverCloser, true);

    this._discretePopoverEscape = function (ev) {
      if (ev.key !== "Escape") return;
      var pop = self._discreteColorPopover;
      if (!pop || pop.hidden) return;
      ev.preventDefault();
      ev.stopPropagation();
      self._closeDiscreteColorPopover(true);
    };
    document.addEventListener("keydown", this._discretePopoverEscape, true);

    this._discreteColorPopover = root;
    return root;
  };

  CryoColorCovariateLegend.prototype._relayoutDiscreteColorPopover = function (anchorEl) {
    var pop = this._discreteColorPopover;
    if (!pop || pop.hidden || !anchorEl) return;
    var ar = anchorEl.getBoundingClientRect();
    pop.style.left = Math.round(ar.left) + "px";
    pop.style.top = Math.round(ar.bottom + 6) + "px";
    var self = this;
    requestAnimationFrame(function () {
      if (!self._discreteColorPopover || self._discreteColorPopover.hidden) return;
      var p = self._discreteColorPopover;
      var pr = p.getBoundingClientRect();
      var vw = window.innerWidth;
      var vh = window.innerHeight;
      var pad = 10;
      var left = pr.left;
      var top = pr.top;
      if (pr.right > vw - pad) left -= pr.right - vw + pad;
      if (pr.bottom > vh - pad) top -= pr.bottom - vh + pad;
      if (left < pad) left = pad;
      if (top < pad) top = pad;
      p.style.left = Math.round(left) + "px";
      p.style.top = Math.round(top) + "px";
    });
  };

  /** Reattach popover cell/anchor after `_renderDiscrete` rebuilds the switch list. */
  CryoColorCovariateLegend.prototype._syncDiscretePopoverAnchorAfterDiscreteRender = function () {
    var st = this._discreteColorPopoverState;
    var pop = this._discreteColorPopover;
    if (!st || !pop || pop.hidden) return;
    var sw = this.discreteSwitches;
    if (!sw) {
      this._closeDiscreteColorPopover(true);
      return;
    }
    var want = String(st.catKey);
    var cells = sw.querySelectorAll(".cryo-cc-discrete-cell[data-cat-key]");
    var cell = null;
    var ci;
    for (ci = 0; ci < cells.length; ci++) {
      if (String(cells[ci].getAttribute("data-cat-key")) === want) {
        cell = cells[ci];
        break;
      }
    }
    if (!cell) {
      this._closeDiscreteColorPopover(true);
      return;
    }
    st.cellEl = cell;
    var anchor =
      cell.querySelector(".cryo-cc-discrete-colorwheel-btn") ||
      cell.querySelector(".cryo-cc-discrete-solo-btn");
    if (!anchor) {
      this._closeDiscreteColorPopover(true);
      return;
    }
    st.anchorEl = anchor;
    var hx = normalizeDiscreteLegendHex(st.commitHex) || "#aab4bf";
    this._applyDiscreteCellPaint(cell, hx);
    this._relayoutDiscreteColorPopover(anchor);
  };

  CryoColorCovariateLegend.prototype._openDiscreteColorPopover = function (
    anchorEl,
    cellEl,
    catKey,
    currentHex
  ) {
    if (this._discreteColorPopover && !this._discreteColorPopover.hidden) {
      this._closeDiscreteColorPopover(true);
    }
    var pop = this._ensureDiscreteColorPopover();
    var hx = normalizeDiscreteLegendHex(currentHex) || "#aab4bf";
    this._discreteColorPopoverState = {
      cellEl: cellEl,
      catKey: String(catKey),
      commitHex: hx,
      anchorEl: anchorEl
    };
    this._discretePopoverSetUiFromHex(hx);
    pop.hidden = false;
    if (!pop.parentNode) {
      document.body.appendChild(pop);
    }
    this._relayoutDiscreteColorPopover(anchorEl);
    var rg0 = pop.querySelector('input[type=\"range\"][data-channel=\"r\"]');
    try {
      if (rg0 && rg0.focus) rg0.focus({ preventScroll: true });
    } catch (e0) {
      try {
        if (rg0) rg0.focus();
      } catch (e1) { /* ignore */ }
    }
  };

  CryoColorCovariateLegend.prototype._renderDiscrete = function (categories) {
    var self = this;
    if (!this.discreteSwitches) return;
    var prevChecked = new Set();
    var prevKeySet = new Set();
    var prevInputs = discreteCategoryInputsFrom(this.discreteSwitches, false);
    for (var pi = 0; pi < prevInputs.length; pi++) {
      var pk = String(prevInputs[pi].getAttribute("data-cat-key"));
      prevKeySet.add(pk);
      if (prevInputs[pi].checked) {
        prevChecked.add(pk);
      }
    }
    this.discreteSwitches.innerHTML = "";
    if (this.invertHeadingSlotEl) {
      this.invertHeadingSlotEl.innerHTML = "";
    }
    if (this.invertHeadingRowEl) {
      this.invertHeadingRowEl.hidden = false;
    }
    this._ensureDiscreteCollapseToggle();
    var cellsMount = this.discreteSwitches;
    var invertMount = this.discreteSwitches;
    if (this.discreteInvertAsideColumn) {
      var cellsWrap = document.createElement("div");
      cellsWrap.className = "cryo-cc-discrete-cells-wrap";
      var invertAside = document.createElement("div");
      invertAside.className = "cryo-cc-discrete-invert-aside";
      this.discreteSwitches.appendChild(cellsWrap);
      this.discreteSwitches.appendChild(invertAside);
      cellsMount = cellsWrap;
      invertMount = invertAside;
    }
    var col = this.getColorColumn();
    var ov =
      typeof this.getDiscreteColorOverrides === "function"
        ? this.getDiscreteColorOverrides()
        : null;
    var byKey = {};
    for (var j = 0; j < (categories || []).length; j++) {
      var cj = categories[j];
      byKey[String(cj.key)] = cj;
    }
    this._allCatKeys = sortDiscreteKeys(Object.keys(byKey));
    var nextKeySet = new Set(this._allCatKeys);
    var preserveDiscreteChecks =
      prevKeySet.size > 0 &&
      prevKeySet.size === nextKeySet.size &&
      (function (a, b) {
        var ok = true;
        a.forEach(function (x) {
          if (!b.has(x)) ok = false;
        });
        return ok;
      }(prevKeySet, nextKeySet));
    for (var i = 0; i < this._allCatKeys.length; i++) {
      var cat = byKey[this._allCatKeys[i]];
      if (!cat) continue;
      var key = String(cat.key);
      var ovHex = ov && typeof ov[key] !== "undefined" && ov[key] != null ? String(ov[key]) : "";
      var hx = normalizeDiscreteLegendHex(ovHex) || normalizeDiscreteLegendHex(cat.color);
      var safeId = "cc-disc-" + i + "-" + key.replace(/[^a-zA-Z0-9_-]/g, "_");
      var initialHex =
        normalizeDiscreteLegendHex(ov && ov[key] != null ? String(ov[key]) : "") ||
        normalizeDiscreteLegendHex(hx) ||
        normalizeDiscreteLegendHex(cat.color) ||
        "#aab4bf";

      // Create cell wrapper; checkbox row uses .cryo-cc-discrete-switch flex container.
      // Keep solo/colour buttons outside <label for> so clicks don't toggle the checkbox twice.
      var cell = document.createElement("div");
      cell.className = "cryo-cc-discrete-cell";
      if (this.discreteChipPlasticLabel) {
        cell.classList.add("cryo-cc-discrete-cell--plastic");
      }
      cell.setAttribute("data-cat-key", key);

      var dispLbl = discreteDisplayLabel(key, col, this.covariateDisplayMap);
      var cntDispRaw = cat.count != null ? cat.count : 0;
      var cntDisp = Number(cntDispRaw);
      if (!isFinite(cntDisp)) cntDisp = 0;
      var countPhrase = this.formatDiscreteDatasetCount(cntDisp);
      cell.title =
        dispLbl +
        " — " +
        countPhrase +
        " in full analysis table. Checkbox toggles visibility on the plot.";
      cell.setAttribute(
        "aria-label",
        dispLbl + ". " + countPhrase + " in full analysis table."
      );

      var controlsWrap = document.createElement("div");
      controlsWrap.className = "cryo-cc-discrete-cell-controls";

      // Colour wheel opens anchored popover (preview + Apply / Cancel).
      if (this.enableDiscreteColorPicker) {
        var pickTitle = dispLbl;

        var wheelSlot = document.createElement("div");
        wheelSlot.className = "cryo-cc-discrete-wheel-slot";

        var pick = document.createElement("button");
        pick.type = "button";
        pick.className = "cryo-cc-discrete-colorwheel-btn";
        pick.title = "Choose colour…";
        pick.setAttribute("aria-label", "Choose colour for " + pickTitle);
        pick.innerHTML = cryoDiscreteColorWheelSvg();

        pick.addEventListener("mousedown", function (ev) {
          ev.preventDefault();
          ev.stopPropagation();
        });
        pick.addEventListener("click", function (ce, kk, ih) {
          return function (ev) {
            ev.preventDefault();
            ev.stopPropagation();
            self._openDiscreteColorPopover(pick, ce, kk, ih);
          };
        }(cell, key, initialHex));

        wheelSlot.appendChild(pick);
        controlsWrap.appendChild(wheelSlot);
      } else {
        controlsWrap.classList.add("cryo-cc-discrete-cell-controls--solo-only");
      }

      var soloBtn = document.createElement("button");
      soloBtn.type = "button";
      soloBtn.className = "cryo-cc-discrete-solo-btn";
      soloBtn.title = "Select only this category";
      soloBtn.setAttribute("aria-label", "Select only " + dispLbl);
      soloBtn.textContent = "1";
      soloBtn.addEventListener("mousedown", function (ev) {
        ev.preventDefault();
        ev.stopPropagation();
      });
      soloBtn.addEventListener("click", function (k) {
        return function (e) {
          e.preventDefault();
          e.stopPropagation();
          self.soloDiscrete(k);
        };
      }(key));
      controlsWrap.appendChild(soloBtn);

      var rowWrap = document.createElement("div");
      rowWrap.className = "cryo-cc-discrete-switch";

      var inp = document.createElement("input");
      inp.type = "checkbox";
      inp.id = safeId;
      inp.setAttribute("role", "switch");
      inp.setAttribute(
        "aria-checked",
        preserveDiscreteChecks
          ? (prevChecked.has(key) ? "true" : "false")
          : (this.discreteCheckedDefault ? "true" : "false")
      );
      inp.setAttribute("data-cat-key", key);
      inp.checked = preserveDiscreteChecks ? prevChecked.has(key) : this.discreteCheckedDefault;

      var lab = document.createElement("label");
      lab.setAttribute("for", safeId);

      var tw = document.createElement("span");
      tw.className = "cryo-cc-discrete-switch-text";
      var sp = document.createElement("span");
      sp.className = "cryo-cc-discrete-switch-label";
      sp.textContent = dispLbl;
      var cn = document.createElement("span");
      cn.className = "cryo-cc-discrete-switch-count";
      cn.textContent = countPhrase;
      cn.title = countPhrase + " in full analysis table";

      tw.appendChild(sp);
      tw.appendChild(cn);
      lab.appendChild(tw);

      rowWrap.appendChild(inp);
      rowWrap.appendChild(lab);
      rowWrap.appendChild(controlsWrap);
      cell.appendChild(rowWrap);
      cellsMount.appendChild(cell);

      this._applyDiscreteCellPaint(cell, initialHex);

      inp.addEventListener("change", function (input) {
        return function () {
          input.setAttribute("aria-checked", input.checked ? "true" : "false");
          if (self._suppressDiscrete) return;
          self._syncInvertBtn();
          self._notify();
        };
      }(inp));
    }

    // Add invert control as its own cell (styling: idle label vs action button — see base CSS).
    if (this.invertBtn) {
      var invertCell = document.createElement("div");
      invertCell.className = "cryo-cc-discrete-cell cryo-cc-discrete-cell--invert";
      invertCell.setAttribute("data-invert-cell", "true");
      invertCell.setAttribute("role", "status");
      invertCell.setAttribute("aria-live", "polite");
      invertCell.setAttribute("aria-disabled", "true");
      invertCell.tabIndex = -1;
      invertCell.classList.add("cryo-cc-discrete-cell--invert-idle");

      var invertInner = document.createElement("div");
      invertInner.className = "cryo-cc-discrete-switch cryo-cc-discrete-switch--invert";

      var fitBox = document.createElement("div");
      fitBox.className = "cryo-cc-discrete-invert-fitbox";

      var invertSpan = document.createElement("span");
      invertSpan.className = "cryo-cc-discrete-switch-label";
      invertSpan.textContent = "";

      fitBox.appendChild(invertSpan);
      invertInner.appendChild(fitBox);
      invertCell.appendChild(invertInner);

      invertCell.addEventListener("click", function (e) {
        if (!invertCell.classList.contains("cryo-cc-discrete-cell--invert-action")) {
          return;
        }
        e.preventDefault();
        e.stopPropagation();
        self.invertDiscrete();
      });
      invertCell.addEventListener("keydown", function (e) {
        if (!invertCell.classList.contains("cryo-cc-discrete-cell--invert-action")) {
          return;
        }
        if (e.key === "Enter" || e.key === " ") {
          e.preventDefault();
          e.stopPropagation();
          self.invertDiscrete();
        }
      });

      if (this.invertHeadingSlotEl) {
        invertCell.classList.add("cryo-cc-discrete-cell--invert-heading");
        invertCell.hidden = true;
        this.invertHeadingSlotEl.hidden = true;
        this.invertHeadingSlotEl.setAttribute("aria-hidden", "true");
        invertMount = this.invertHeadingSlotEl;
      } else if (this.discreteInvertAsideColumn) {
        invertCell.classList.add("cryo-cc-discrete-cell--invert-aside");
      } else if (this.discreteInvertFooterRow) {
        invertCell.classList.add("cryo-cc-discrete-cell--invert-footer");
      }

      invertMount.appendChild(invertCell);
    }

    this._syncInvertBtn();
    if (this._discretePanelCollapsed) {
      this.setDiscretePanelCollapsed(true);
    }
    this._syncDiscretePopoverAnchorAfterDiscreteRender();
    var afterDiscreteLay = function () {
      self.fitDiscreteSwitchColumnWidths();
      self.measureDiscreteInvertAsideColumn();
      self.scheduleInvertTypographyFit();
      if (typeof self.onDiscreteLayout === "function") {
        self.onDiscreteLayout();
      }
    };
    requestAnimationFrame(function () {
      requestAnimationFrame(afterDiscreteLay);
    });
  };

  /**
   * Measure widest toggle cell, cap at 79% of legend width, and set ``--cryo-discrete-col-w`` so
   * CSS ``repeat(auto-fill, …)`` packs as many columns as fit in the toggle-selection region.
   */
  CryoColorCovariateLegend.prototype.fitDiscreteSwitchColumnWidths = function () {
    if (!this.discreteSwitches) return;
    if (this._mode !== "discrete") {
      this.discreteSwitches.style.removeProperty("--cryo-discrete-col-w");
      this.discreteSwitches.style.removeProperty("--cryo-discrete-cell-max-w");
      this.discreteSwitches.style.removeProperty("--cryo-discrete-cell-min-h");
      return;
    }
    var gridHost =
      this.discreteSwitches.querySelector(".cryo-cc-discrete-cells-wrap") ||
      this.discreteSwitches;
    var regionRoot =
      this.discreteWrap ||
      this.discreteSwitches.closest("#color-discrete-wrap") ||
      this.discreteSwitches.closest("#latent3d-color-discrete-wrap") ||
      this.discreteSwitches.closest("#pair-color-discrete-wrap") ||
      gridHost.parentElement;
    var legendRoot =
      this.discreteSwitches.closest("#latent3d-color-legend-panel") ||
      this.discreteSwitches.closest(".pairplot-color-legend-aside") ||
      this.discreteSwitches.closest(".cryo-cc-legend") ||
      this.discreteSwitches.closest("#color-hist-panel") ||
      regionRoot ||
      this.discreteSwitches.parentElement;
    var legendW =
      legendRoot && legendRoot.getBoundingClientRect
        ? legendRoot.getBoundingClientRect().width
        : 0;
    var regionW =
      regionRoot && regionRoot.getBoundingClientRect
        ? regionRoot.getBoundingClientRect().width
        : gridHost.clientWidth || 0;
    var maxColW = legendW > 1 ? Math.floor(legendW * 0.79) : 0;
    if (maxColW > 0) {
      this.discreteSwitches.style.setProperty("--cryo-discrete-cell-max-w", maxColW + "px");
    } else {
      this.discreteSwitches.style.removeProperty("--cryo-discrete-cell-max-w");
    }
    this.discreteSwitches.style.removeProperty("--cryo-discrete-cell-min-h");

    var cells = gridHost.querySelectorAll(
      ".cryo-cc-discrete-cell:not(.cryo-cc-discrete-cell--invert):not(.cryo-cc-discrete-cell--invert-heading)"
    );
    if (!cells.length) {
      this.discreteSwitches.style.removeProperty("--cryo-discrete-col-w");
      return;
    }
    this.discreteSwitches.style.removeProperty("--cryo-discrete-col-w");
    var maxCell = 0;
    for (var ci = 0; ci < cells.length; ci++) {
      var cell = cells[ci];
      var w = cell.scrollWidth;
      if (w > maxCell) maxCell = w;
    }
    if (maxCell < 1) {
      this.discreteSwitches.style.removeProperty("--cryo-discrete-col-w");
      return;
    }
    var colW = Math.ceil(maxCell + 1);
    if (maxColW > 0) colW = Math.min(colW, maxColW);
    if (regionW > 1 && colW > regionW) colW = Math.max(1, Math.floor(regionW));
    this.discreteSwitches.style.setProperty("--cryo-discrete-col-w", colW + "px");
  };

  CryoColorCovariateLegend.prototype.measureDiscreteInvertAsideColumn = function () {
    if (!this.discreteInvertAsideColumn || !this.discreteSwitches) return;
    var wrap = this.discreteSwitches.querySelector(".cryo-cc-discrete-cells-wrap");
    if (!wrap) return;
    var probe = wrap.querySelector(".cryo-cc-discrete-cell:not(.cryo-cc-discrete-cell--invert)");
    if (probe && probe.offsetWidth > 0 && probe.offsetHeight > 0) {
      this.discreteSwitches.style.setProperty(
        "--cryo-explorer-discrete-cell-w",
        probe.offsetWidth + "px"
      );
      this.discreteSwitches.style.setProperty(
        "--cryo-explorer-discrete-cell-h",
        probe.offsetHeight + "px"
      );
    }
    this.scheduleInvertTypographyFit();
  };

  CryoColorCovariateLegend.prototype.scheduleInvertTypographyFit = function () {
    var self = this;
    if (!this.discreteSwitches || this._invertTypographyFitScheduled) return;
    this._invertTypographyFitScheduled = true;
    requestAnimationFrame(function () {
      self._invertTypographyFitScheduled = false;
      var cell = self._getInvertCell();
      if (cell) self._fitInvertCellTypography(cell);
    });
  };

  CryoColorCovariateLegend.prototype._fitInvertCellTypography = function (invertCell) {
    fitInvertCellTypography(invertCell);
  };

  CryoColorCovariateLegend.prototype._discreteChipBlendForPaint = function () {
    if (!this.discreteChipMatplotlibBlend || this._discreteChipBlendAlpha == null || !this._discreteChipBlendBg) {
      return null;
    }
    return { alpha: this._discreteChipBlendAlpha, bgHex: this._discreteChipBlendBg };
  };

  CryoColorCovariateLegend.prototype._applyDiscreteCellPaint = function (
    cellEl,
    hex
  ) {
    var hx = normalizeDiscreteLegendHex(hex);
    if (!hx) return;
    var plastic =
      this.discreteChipPlasticLabel &&
      !cellEl.classList.contains("cryo-cc-discrete-cell--invert");
    var dSt = discretePlotMatchedChipStyles(
      hx,
      this._discreteChipBlendForPaint(),
      plastic
    );
    var lab = cellEl.querySelector(".cryo-cc-discrete-switch");
    var sp = cellEl.querySelector(".cryo-cc-discrete-switch-label");
    var cn = cellEl.querySelector(".cryo-cc-discrete-switch-count");
    if (!lab) return;
    lab.style.backgroundColor = "transparent";
    lab.style.border = "none";
    if (dSt && dSt.bg) {
      cellEl.style.backgroundColor = dSt.bg;
      cellEl.style.border = "1px solid " + dSt.borderColor;
    }
    if (sp) {
      if (dSt && dSt.labelColor) sp.style.color = dSt.labelColor;
      else sp.style.color = "";
    }
    if (cn) {
      if (dSt && dSt.countColor) cn.style.color = dSt.countColor;
      else cn.style.color = "";
    }
  };

  CryoColorCovariateLegend.prototype._syncInvertBtn = function () {
    if (!this.discreteSwitches) return;
    var inputs = discreteCategoryInputsFrom(this.discreteSwitches, false);
    var any = false;
    for (var i = 0; i < inputs.length; i++) {
      if (inputs[i].checked) { any = true; break; }
    }
    /** Verbal aria label (no ⇄ character); chip shows ⇄ beside phrase for visual invert cue */
    var invertVerbPhrase = "Invert selection";
    var invertDisplayed = invertVerbPhrase + " \u21c4";

    // Update original invert button (hidden but kept for accessibility/legacy)
    if (this.invertBtn) {
      this.invertBtn.disabled = !any;
      this.invertBtn.textContent = any ? invertDisplayed : "";
    }

    var invertCell = this._getInvertCell();
    var headingInvert = !!(
      invertCell
      && invertCell.classList.contains("cryo-cc-discrete-cell--invert-heading")
    );

    if (invertCell && headingInvert) {
      var hLabel = invertCell.querySelector(".cryo-cc-discrete-switch-label");
      if (this._discretePanelCollapsed) {
        invertCell.hidden = true;
        if (this.invertHeadingSlotEl) {
          this.invertHeadingSlotEl.hidden = true;
          this.invertHeadingSlotEl.setAttribute("aria-hidden", "true");
        }
        return;
      }
      invertCell.hidden = !any;
      if (this.invertHeadingSlotEl) {
        this.invertHeadingSlotEl.hidden = !any;
        this.invertHeadingSlotEl.setAttribute("aria-hidden", any ? "false" : "true");
      }
      if (!any) {
        if (hLabel) {
          hLabel.textContent = "";
        }
        invertCell.classList.remove(
          "cryo-cc-discrete-cell--invert-idle",
          "cryo-cc-discrete-cell--invert-action"
        );
        invertCell.setAttribute("role", "presentation");
        invertCell.setAttribute("aria-hidden", "true");
        invertCell.removeAttribute("aria-label");
        invertCell.removeAttribute("tabIndex");
        invertCell.removeAttribute("aria-live");
        return;
      }
      invertCell.removeAttribute("aria-hidden");
      if (hLabel) {
        hLabel.textContent = invertDisplayed;
      }
      invertCell.classList.remove("cryo-cc-discrete-cell--invert-idle");
      invertCell.classList.add("cryo-cc-discrete-cell--invert-action");
      invertCell.setAttribute("role", "button");
      invertCell.setAttribute("aria-disabled", "false");
      invertCell.setAttribute("aria-label", invertVerbPhrase);
      invertCell.removeAttribute("aria-live");
      invertCell.tabIndex = 0;
      this.scheduleInvertTypographyFit();
      return;
    }

    if (!invertCell) return;

    var invertLabel = invertCell.querySelector(".cryo-cc-discrete-switch-label");
    if (invertLabel) {
      invertLabel.textContent = any ? invertDisplayed : (this.noDiscreteSelectionText || "");
    }
    invertCell.classList.toggle("cryo-cc-discrete-cell--invert-idle", !any);
    invertCell.classList.toggle("cryo-cc-discrete-cell--invert-action", any);
    invertCell.hidden = false;
    invertCell.removeAttribute("aria-hidden");
    invertCell.setAttribute("role", any ? "button" : "status");
    invertCell.setAttribute("aria-disabled", any ? "false" : "true");
    invertCell.setAttribute(
      "aria-label",
      any ? invertVerbPhrase : (this.noDiscreteSelectionText || "")
    );
    invertCell.tabIndex = any ? 0 : -1;
    if (any) {
      invertCell.removeAttribute("aria-live");
    } else {
      invertCell.setAttribute("aria-live", "polite");
    }
    this.scheduleInvertTypographyFit();
  };

  CryoColorCovariateLegend.prototype.setDiscreteCheckedKeys = function (keys) {
    if (!this.discreteSwitches) return;
    var keySet = new Set((keys || []).map(function (k) { return String(k); }));
    this._suppressDiscrete = true;
    try {
      var inputs = discreteCategoryInputsFrom(this.discreteSwitches, false);
      for (var i = 0; i < inputs.length; i++) {
        var k = inputs[i].getAttribute("data-cat-key");
        var checked = keySet.has(String(k));
        inputs[i].checked = checked;
        inputs[i].setAttribute("aria-checked", checked ? "true" : "false");
      }
    } finally {
      this._suppressDiscrete = false;
    }
    this._syncInvertBtn();
  };

  CryoColorCovariateLegend.prototype.invertDiscrete = function () {
    if (!this.discreteSwitches) return;
    this._suppressDiscrete = true;
    try {
      var inputs = discreteCategoryInputsFrom(this.discreteSwitches, false);
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

  CryoColorCovariateLegend.prototype.soloDiscrete = function (soloKey) {
    if (!this.discreteSwitches) return;
    this._suppressDiscrete = true;
    try {
      var inputs = discreteCategoryInputsFrom(this.discreteSwitches, false);
      for (var i = 0; i < inputs.length; i++) {
        var k = inputs[i].getAttribute("data-cat-key");
        var checked = String(k) === String(soloKey);
        inputs[i].checked = checked;
        inputs[i].setAttribute("aria-checked", checked ? "true" : "false");
      }
    } finally {
      this._suppressDiscrete = false;
    }
    this._syncInvertBtn();
    this._notify();
  };

  /** SVG markup for the discrete-legend colour-wheel icon (particle explorer region chips, etc.). */
  CryoColorCovariateLegend.discreteColorWheelSvgInnerHTML = function () {
    return cryoDiscreteColorWheelSvg();
  };

  CryoColorCovariateLegend.sortDiscreteKeys = sortDiscreteKeys;
  CryoColorCovariateLegend.paletteScaleCSS = paletteScaleCSS;
  CryoColorCovariateLegend.formatThresholdValue = formatThresholdValue;
  CryoColorCovariateLegend.normalizeDiscreteLegendHex = normalizeDiscreteLegendHex;
  CryoColorCovariateLegend.hasPalette = function (paletteName) {
    return Object.prototype.hasOwnProperty.call(PALETTE_RGB, paletteName);
  };
  CryoColorCovariateLegend.wirePaletteSelect = function (opts) {
    opts = opts || {};
    var toggle = opts.toggle;
    var options = opts.options;
    var select = opts.select;
    var hiddenClass = opts.hiddenClass;
    function sync() {
      if (!options || !toggle || !select) return;
      var show = toggle.getAttribute("aria-expanded") === "true"
        && !(hiddenClass && select.classList.contains(hiddenClass));
      options.classList.toggle("pairplot-palette-select__options--collapsed", !show);
    }
    if (toggle) {
      toggle.addEventListener("click", function () {
        var expanded = toggle.getAttribute("aria-expanded") === "true";
        toggle.setAttribute("aria-expanded", expanded ? "false" : "true");
        sync();
      });
    }
    return {
      sync: sync,
      close: function () {
        if (toggle) toggle.setAttribute("aria-expanded", "false");
        sync();
      }
    };
  };

  global.CryoColorCovariateLegend = CryoColorCovariateLegend;
})(typeof window !== "undefined" ? window : globalThis);
