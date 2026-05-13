/**
 * Colour covariate histogram + discrete toggles shared by particle explorer,
 * pair-grid, and 3-D latent views. Continuous legends support click thresholds
 * and drag ranges; optional `histPlotVertical` swaps covariate onto y and
 * density onto x (pair-plot/3-D asides). Shared vertical palette menu styling:
 * ``static/css/cryo_cc_legend_palette_menu.css`` (classes ``cryo-cc-legend-palette-shell``
 * on the select wrapper and ``cryo-cc-legend-palette-menu`` on the options pane). Discrete colour picking uses one anchored
 * panel (RGB sliders + hex field + Apply / Cancel) with no separate OS colour dialog;
 * Apply commits but leaves the panel open (Cancel, Escape, or outside click closes). Plastic styling for toggle chips is on by default
 * (`discreteChipPlasticLabel` opt-out). Optional `histVerticalMargins` shallow-
 * merges into the vertical histogram layout margin. Panel `vertical` controls CSS layout class.
 */
(function (global) {
  "use strict";

  var PALETTE_RGB = {
    Viridis: [[68, 1, 84], [72, 39, 119], [63, 74, 138], [49, 104, 142], [38, 130, 142], [31, 157, 138], [53, 183, 121], [109, 205, 89], [253, 231, 37]],
    Plasma: [[13, 8, 135], [76, 2, 161], [126, 3, 168], [170, 35, 149], [204, 71, 120], [230, 108, 92], [248, 149, 64], [253, 197, 39], [240, 249, 33]],
    Inferno: [[0, 0, 4], [41, 10, 108], [108, 2, 112], [187, 55, 84], [249, 142, 8], [252, 208, 74], [252, 254, 164]],
    Magma: [[0, 0, 4], [59, 15, 112], [122, 4, 111], [203, 62, 65], [252, 141, 89], [252, 208, 136], [252, 253, 191]],
    Cividis: [[0, 34, 78], [0, 57, 112], [68, 90, 129], [115, 128, 131], [159, 161, 135], [206, 186, 122], [253, 231, 37]],
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

  function injectCryoCcLegendStylesOnce() {
    if (document.getElementById("cryo-cc-legend-dynamic-style")) {
      return;
    }
    var st = document.createElement("style");
    st.id = "cryo-cc-legend-dynamic-style";
    st.textContent = ""
      // Single flex row: [checkbox][gap][label stack][icons flush]; chip paint on .cryo-cc-discrete-cell
      + ".cryo-cc-discrete-cell{box-sizing:border-box;min-width:0;padding:0;border-radius:2px;"
      + "border:1px solid transparent;}"
      + ".cryo-cc-discrete-switch{display:flex;flex-direction:row;flex-wrap:nowrap;align-items:center;"
      + "gap:0;margin:0;padding:0;width:100%;min-width:0;"
      + "cursor:default;user-select:none;border:none;background:transparent;border-radius:0;}"
      + ".cryo-cc-discrete-switch > input[type=\"checkbox\"]{flex-shrink:0;margin:0 0.35rem 0 0;"
      + "accent-color:#c4703a;cursor:pointer;width:0.62rem;height:0.62rem;}"
      + ".cryo-cc-discrete-switch > label[for]{flex:1 1 auto;min-width:0;margin:0;padding:0;"
      + "cursor:pointer;display:flex;flex-direction:column;align-items:flex-start;justify-content:center;}"
      + ".cryo-cc-discrete-switch-text{display:flex;flex-direction:column;align-items:flex-start;"
      + "flex:1 1 auto;min-width:0;line-height:1;}"
      + ".cryo-cc-discrete-cell-controls{display:flex;flex-direction:column;align-items:center;"
      + "justify-content:center;gap:0;flex-shrink:0;margin:0;padding:0;line-height:0;}"
      + ".cryo-cc-discrete-cell-controls--solo-only{justify-content:center;}"
      // Colour wheel slot (native picker sits over the icon button)
      + ".cryo-cc-discrete-wheel-slot{position:relative;flex-shrink:0;width:1rem;height:1rem;}"
      + ".cryo-cc-discrete-native-color{position:absolute;inset:0;width:100%;height:100%;"
      + "opacity:0;pointer-events:none;z-index:0;padding:0;border:none;margin:0;cursor:pointer;}"
      + ".cryo-cc-discrete-colorwheel-btn{position:absolute;inset:0;"
      + "display:inline-flex;align-items:center;justify-content:center;"
      + "width:100%;height:100%;margin:0;padding:0;"
      + "border:none;border-radius:999px;background:transparent;"
      + "cursor:pointer;color:#243b53;line-height:0;z-index:2;}"
      + ".cryo-cc-discrete-colorwheel-btn:focus-visible{outline:2px solid var(--accent,#c4703a);outline-offset:1px;}"
      + ".cryo-cc-discrete-colorwheel-btn:hover{transform:scale(1.15);}"
      + ".cryo-cc-discrete-colorwheel-btn svg{display:block;width:0.83rem;height:0.83rem;}"
      // Solo \"1\" button — same column as wheel, vertically separated by flex layout
      + ".cryo-cc-discrete-solo-btn{display:inline-flex;align-items:center;justify-content:center;"
      + "flex-shrink:0;width:1rem;height:1rem;margin:0;padding:0;"
      + "border:none;border-radius:2px;background:transparent;"
      + "cursor:pointer;color:#1a2332;font-size:0.58rem;font-weight:900;"
      + "font-family:Georgia,serif;font-style:italic;line-height:1;"
      + "text-shadow:-1px -1px 0 #fff,1px -1px 0 #fff,-1px 1px 0 #fff,1px 1px 0 #fff,"
      + "-1px 0 0 #fff,1px 0 0 #fff,0 -1px 0 #fff,0 1px 0 #fff,"
      + "0 0 3px rgba(255,255,255,0.95);}"
      + ".cryo-cc-discrete-solo-btn:focus-visible{outline:2px solid var(--accent,#c4703a);outline-offset:1px;}"
      + ".cryo-cc-discrete-solo-btn:hover{color:#000;transform:scale(1.2);}"
      // Invert chip matches data chips (paint applied via _applyDiscreteCellPaint)
      + ".cryo-cc-discrete-cell--invert{cursor:pointer;transition:opacity 0.15s ease;padding:0;}"
      + ".cryo-cc-discrete-cell--invert:hover{opacity:0.85;}"
      + ".cryo-cc-discrete-cell--invert .cryo-cc-discrete-switch{cursor:inherit;}"
      + ".cryo-cc-discrete-cell--plastic:not(.cryo-cc-discrete-cell--invert){"
      + "position:relative;border-radius:7px;overflow:visible;"
      + "box-shadow:"
      + "inset 0 2px 3px rgba(255,255,255,0.68),"
      + "inset 0 -3px 6px rgba(0,0,0,0.16),"
      + "inset 0 0 0 1px rgba(255,255,255,0.28),"
      + "0 2px 5px rgba(26,35,50,0.2);}"
      + ".cryo-cc-discrete-color-popover{position:fixed;z-index:2147483646;"
      + "background:var(--paper,#fafaf8);color:var(--nav-bg,#243b53);"
      + "border:1px solid rgba(36,59,83,0.26);border-radius:10px;"
      + "padding:0.65rem 0.72rem;box-shadow:0 14px 44px rgba(36,59,83,0.22);"
      + "display:flex;flex-direction:column;gap:0.55rem;min-width:14.25rem;"
      + "max-width:calc(100vw - 20px);}"
      + ".cryo-cc-discrete-color-popover[hidden]{display:none!important;}"
      + ".cryo-cc-discrete-color-popover-title{font-size:0.82rem;font-weight:700;"
      + "letter-spacing:0.02em;color:#1a2332;margin:0;line-height:1.2;}"
      + ".cryo-cc-discrete-color-popover-preview{width:100%;height:2.85rem;"
      + "border-radius:8px;border:1px solid rgba(36,59,83,0.22);"
      + "box-sizing:border-box;box-shadow:inset 0 1px 2px rgba(255,255,255,0.35);}"
      + ".cryo-cc-discrete-color-popover-sliders{display:flex;flex-direction:column;"
      + "gap:0.38rem;width:100%;}"
      + ".cryo-cc-discrete-color-popover-slider-row{display:flex;align-items:center;"
      + "gap:0.42rem;width:100%;min-width:0;}"
      + ".cryo-cc-discrete-color-popover-slider-row label{flex:0 0 1.05rem;"
      + "font-size:0.68rem;font-weight:800;color:#334e68;text-align:center;line-height:1;}"
      + ".cryo-cc-discrete-color-popover-slider-row input[type=range]{flex:1 1 auto;"
      + "min-width:0;height:0.42rem;accent-color:#c4703a;}"
      + ".cryo-cc-discrete-color-popover-hex-row{display:flex;align-items:center;"
      + "gap:0.45rem;width:100%;min-width:0;}"
      + ".cryo-cc-discrete-color-popover-hex-row label{flex:0 0 auto;font-size:0.72rem;"
      + "font-weight:700;color:#334e68;}"
      + ".cryo-cc-discrete-color-popover-hex{flex:1 1 auto;min-width:0;font-family:ui-monospace,Menlo,Consolas,monospace;"
      + "font-size:0.76rem;padding:0.28rem 0.38rem;border:1px solid rgba(36,59,83,0.22);"
      + "border-radius:6px;background:#fff;color:#1a2332;box-sizing:border-box;}"
      + ".cryo-cc-discrete-color-popover-actions{display:flex;justify-content:flex-end;"
      + "gap:0.42rem;flex-wrap:wrap;padding-top:0.15rem;"
      + "border-top:1px solid rgba(36,59,83,0.1);margin-top:0.05rem;}"
      + ".cryo-cc-discrete-color-popover-apply,.cryo-cc-discrete-color-popover-cancel{"
      + "font:inherit;font-size:0.78rem;font-weight:600;padding:0.32rem 0.72rem;"
      + "border-radius:6px;cursor:pointer;line-height:1.2;"
      + "touch-action:manipulation;}"
      + ".cryo-cc-discrete-color-popover-apply{border:1px solid rgba(196,112,58,0.55);"
      + "background:#c4703a;color:#fff;}"
      + ".cryo-cc-discrete-color-popover-apply:hover{filter:brightness(1.05);}"
      + ".cryo-cc-discrete-color-popover-cancel{border:1px solid rgba(36,59,83,0.28);"
      + "background:rgba(255,255,255,0.95);color:#243b53;}"
      + ".cryo-cc-discrete-color-popover-cancel:hover{background:#fff;"
      + "border-color:rgba(36,59,83,0.4);}"
      + ".cryo-cc-discrete-switch-checkbox-spacer{flex-shrink:0;width:0.62rem;height:0.62rem;"
      + "margin:0 0.35rem 0 0;box-sizing:border-box;pointer-events:none;}"
      // Text label and count - readable on colored backgrounds
      + ".cryo-cc-discrete-switch-label{font-weight:600;font-size:0.7rem;line-height:1.1;}"
      + ".cryo-cc-discrete-cell--invert .cryo-cc-discrete-invert-fitbox .cryo-cc-discrete-switch-label{font-weight:700;}"
      + ".cryo-cc-discrete-switch-count{font-size:0.6rem;line-height:1.1;}"
      // Legacy pick button (kept for backwards compatibility)
      + ".cryo-cc-discrete-pick-btn{display:inline-flex;align-items:center;justify-content:center;"
      + "flex-shrink:0;width:1.32rem;height:1.32rem;margin:0;padding:0;"
      + "border:1px solid rgba(36,59,83,0.28);border-radius:999px;background:rgba(255,255,255,0.92);"
      + "cursor:pointer;color:#243b53;line-height:0;"
      + "touch-action:manipulation;transition:background 0.12s ease,border-color 0.12s ease;}"
      + ".cryo-cc-discrete-pick-btn:focus-visible{outline:2px solid var(--accent,#c4703a);outline-offset:2px;}"
      + ".cryo-cc-discrete-pick-btn:hover{border-color:rgba(36,59,83,0.45);background:#fff;"
      + "box-shadow:0 1px 3px rgba(36,59,83,0.12);}"
      + ".cryo-cc-discrete-pick-btn svg{display:block;width:0.92rem;height:0.92rem;}"
      // RGB dialog styles (kept for backwards compatibility)
      + ".cryo-cc-rgb-dialog-backdrop{position:fixed;inset:0;background:rgba(26,35,50,0.38);"
      + "z-index:2147483000;display:flex;align-items:center;justify-content:center;padding:1rem;}"
      + ".cryo-cc-rgb-dialog{background:var(--paper,#fff);color:var(--nav-bg,#243b53);border-radius:8px;"
      + "box-shadow:0 12px 40px rgba(36,59,83,0.22);padding:1rem 1.1rem;width:min(22rem,calc(100vw - 2rem));"
      + "border:1px solid rgba(36,59,83,0.14);}"
      + ".cryo-cc-rgb-dialog h3{margin:0 0 0.62rem;font-size:1.02rem;line-height:1.25;color:var(--nav-bg,#243b53);}"
      + ".cryo-cc-rgb-row{display:flex;align-items:center;gap:0.5rem;margin:0 0 0.55rem;flex-wrap:wrap;}"
      + ".cryo-cc-rgb-row label{font-size:0.82rem;font-weight:600;min-width:4.25rem;color:var(--nav-bg,#243b53);}"
      + ".cryo-cc-rgb-row input[type=number]{width:4.75rem;font-size:0.88rem;padding:0.22rem 0.35rem;"
      + "border:1px solid rgba(36,59,83,0.22);border-radius:4px;}"
      + ".cryo-cc-rgb-colorwheel{width:3.05rem;height:2.05rem;padding:2px;"
      + "border:1px solid rgba(36,59,83,0.26);border-radius:6px;cursor:pointer;background:#fff;"
      + "touch-action:manipulation;}"
      + ".cryo-cc-rgb-dialog-actions{display:flex;justify-content:flex-end;gap:0.48rem;"
      + "flex-wrap:wrap;margin-top:0.82rem;padding-top:0.58rem;"
      + "border-top:1px solid rgba(36,59,83,0.08);}";
    document.head.appendChild(st);
  }

  injectCryoCcLegendStylesOnce();

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
    // Rainbow-themed color wheel icon
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

  function createDiscreteRgbPopover() {
    var backdrop = document.createElement("div");
    backdrop.className = "cryo-cc-rgb-dialog-backdrop";
    backdrop.hidden = true;
    var dlg = document.createElement("div");
    dlg.className = "cryo-cc-rgb-dialog";
    dlg.setAttribute("role", "dialog");
    dlg.setAttribute("aria-modal", "true");
    var titleEl = document.createElement("h3");
    dlg.appendChild(titleEl);

    function numRow(txt, inp) {
      var row = document.createElement("div");
      row.className = "cryo-cc-rgb-row";
      var lab = document.createElement("label");
      lab.textContent = txt;
      row.appendChild(lab);
      row.appendChild(inp);
      dlg.appendChild(row);
      return inp;
    }

    var clr = document.createElement("input");
    clr.type = "color";
    clr.className = "cryo-cc-rgb-colorwheel";
    clr.setAttribute("aria-label", "Colour wheel picker");
    numRow("Preview", clr);

    var rIn = document.createElement("input");
    var gIn = document.createElement("input");
    var bIn = document.createElement("input");
    rIn.type = "number"; gIn.type = "number"; bIn.type = "number";
    rIn.min = "0"; rIn.max = "255"; gIn.min = "0"; gIn.max = "255"; bIn.min = "0"; bIn.max = "255";
    numRow("Red", rIn);
    numRow("Green", gIn);
    numRow("Blue", bIn);

    var actions = document.createElement("div");
    actions.className = "cryo-cc-rgb-dialog-actions";
    var cancelBt = document.createElement("button");
    cancelBt.type = "button";
    cancelBt.className = "btn btn-secondary";
    cancelBt.textContent = "Cancel";
    var okBt = document.createElement("button");
    okBt.type = "button";
    okBt.className = "btn";
    okBt.textContent = "Apply";
    actions.appendChild(cancelBt);
    actions.appendChild(okBt);
    dlg.appendChild(actions);
    backdrop.appendChild(dlg);
    document.body.appendChild(backdrop);

    function clampByte(v) {
      var n = Math.round(Number(v));
      if (!isFinite(n)) return 0;
      return Math.max(0, Math.min(255, n));
    }

    function syncRgbFromClr() {
      var hx = normalizeDiscreteLegendHex(clr.value || "#808080") || "#808080";
      var t = discreteLegendRgbTuple(hx);
      rIn.value = String(t.r);
      gIn.value = String(t.g);
      bIn.value = String(t.b);
      clr.value = hx;
    }

    function syncClrFromRgb() {
      var rs = clampByte(rIn.value).toString(16);
      var gs = clampByte(gIn.value).toString(16);
      var bs = clampByte(bIn.value).toString(16);
      clr.value = "#"
        + (rs.length === 1 ? "0" : "") + rs
        + (gs.length === 1 ? "0" : "") + gs
        + (bs.length === 1 ? "0" : "") + bs;
    }

    clr.addEventListener("input", function () { syncRgbFromClr(); });
    rIn.addEventListener("input", function () { syncClrFromRgb(); });
    gIn.addEventListener("input", function () { syncClrFromRgb(); });
    bIn.addEventListener("input", function () { syncClrFromRgb(); });

    var onCommit = null;
    var trapKey = null;
    function close() {
      backdrop.hidden = true;
      backdrop.setAttribute("aria-hidden", "true");
      titleEl.textContent = "";
      onCommit = null;
      if (trapKey) {
        document.removeEventListener("keydown", trapKey, true);
        trapKey = null;
      }
    }

    backdrop.addEventListener("click", function (ev) {
      if (ev.target === backdrop) {
        close();
      }
    });
    cancelBt.addEventListener("click", function () { close(); });
    okBt.addEventListener("click", function () {
      syncClrFromRgb();
      var hxFin = normalizeDiscreteLegendHex(clr.value);
      if (onCommit && hxFin) {
        onCommit(hxFin);
      }
      close();
    });

    function open(meta) {
      meta = meta || {};
      injectCryoCcLegendStylesOnce();
      var start = normalizeDiscreteLegendHex(meta.initialHex || "#8899aa");
      if (!start) start = "#8899aa";
      titleEl.textContent = meta.title || "Choose RGB colour";
      clr.value = start;
      syncRgbFromClr();
      dlg.setAttribute("aria-label", meta.title || "RGB colour picker");
      onCommit = typeof meta.onCommit === "function" ? meta.onCommit : null;
      backdrop.hidden = false;
      backdrop.setAttribute("aria-hidden", "false");
      try {
        clr.focus();
      } catch (fe) { /* ignore */ }
      if (!trapKey) {
        trapKey = function (kev) {
          if (kev.key === "Escape") {
            kev.preventDefault();
            close();
          }
        };
        document.addEventListener("keydown", trapKey, true);
      }
    }

    return { open: open, close: close, element: backdrop };
  }

  var _cryoRgbPopoverInst = null;
  function discreteRgbPopover() {
    if (!_cryoRgbPopoverInst) _cryoRgbPopoverInst = createDiscreteRgbPopover();
    return _cryoRgbPopoverInst;
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

  /**
   * Discrete chip fill matches plotted markers: optional blend over dashboard cream
   * using the same alpha as matplotlib pair-grid lower-triangle scatter.
   */
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

  /** Full-table counts shown on discrete covariate chips (not plotted subsample). */
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
    this.getLegendContextData = opts.getLegendContextData || null;
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
    /** When true, category cells stay in a wrapping grid and invert sits in a dedicated right column (legacy). Ignored when invertHeadingSlotEl is set. */
    this.discreteInvertAsideColumn =
      opts.discreteInvertAsideColumn === true && !this.invertHeadingSlotEl;
    /** When true, invert sits on its own final row spanning the toggle grid (legacy). Ignored when invertHeadingSlotEl is set. */
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
    this.getDiscreteColorOverrides = opts.getDiscreteColorOverrides || null;
    this.enableDiscreteColorPicker = opts.enableDiscreteColorPicker !== false;
    this.onDiscreteColorChange = opts.onDiscreteColorChange || null;
    /** When set (3-D viewer), hide the legend panel entirely if no colour covariate — skip heading-only placeholder. */
    this.hideLegendPanelWhenNoCovariate = opts.hideLegendPanelWhenNoCovariate === true;
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
    }
    if (this.discreteSwitches) {
      this.discreteSwitches.innerHTML = "";
    }
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
      var inputs = this.discreteSwitches
        ? this.discreteSwitches.querySelectorAll("input[type=\"checkbox\"][data-cat-key]:checked")
        : [];
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

  CryoColorCovariateLegend.prototype._updateRegionHeading = function () {
    if (!this.regionHeadingEl) return;
    if (this._mode === "continuous") {
      this.regionHeadingEl.hidden = false;
      this.regionHeadingEl.textContent = this.regionHeadingContinuousText;
    } else if (this._mode === "discrete") {
      this.regionHeadingEl.hidden = false;
      this.regionHeadingEl.textContent = this.regionHeadingDiscreteText;
    } else {
      this.regionHeadingEl.hidden = true;
      this.regionHeadingEl.textContent = "";
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
          body: JSON.stringify({ column: col })
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
        margin: { l: 4, r: 4, t: 2, b: 17 },
        dragmode: false,
        xaxis: {
          range: violin.range,
          showgrid: true,
          zeroline: false,
          fixedrange: true,
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
    if (this.histPlotVertical) {
      return {
        xref: "paper",
        yref: "y",
        x: 1,
        xanchor: "right",
        y: value,
        yanchor: "middle",
        text: formatThresholdValue(value),
        showarrow: false,
        bgcolor: "rgba(255,255,255,0.94)",
        bordercolor: "rgba(36,59,83,0.35)",
        borderwidth: 1,
        borderpad: 3,
        font: { size: 10, color: "#243b53" }
      };
    }
    return {
      xref: "x",
      yref: "paper",
      x: value,
      y: 0,
      yanchor: "top",
      text: formatThresholdValue(value),
      showarrow: false,
      bgcolor: "rgba(255,255,255,0.94)",
      bordercolor: "rgba(36,59,83,0.35)",
      borderwidth: 1,
      borderpad: 3,
      font: { size: 10, color: "#243b53" }
    };
  };

  CryoColorCovariateLegend.prototype._rangeAnnotation = function (lo, hi) {
    var text = rangeInequalityPhrase(lo, hi, this._activeColorLabel());
    if (this.histPlotVertical) {
      return {
        xref: "paper",
        yref: "y",
        x: 1,
        xanchor: "right",
        y: (lo + hi) * 0.5,
        yanchor: "middle",
        text: text,
        showarrow: false,
        bgcolor: "rgba(255,255,255,0.94)",
        bordercolor: "rgba(36,59,83,0.35)",
        borderwidth: 1,
        borderpad: 3,
        font: { size: 10, color: "#243b53" }
      };
    }
    return {
      xref: "x",
      yref: "paper",
      x: (lo + hi) * 0.5,
      y: 0,
      yanchor: "top",
      text: text,
      showarrow: false,
      bgcolor: "rgba(255,255,255,0.94)",
      bordercolor: "rgba(36,59,83,0.35)",
      borderwidth: 1,
      borderpad: 3,
      font: { size: 10, color: "#243b53" }
    };
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
    var prevInputs = this.discreteSwitches.querySelectorAll(
      "input[type=\"checkbox\"][data-cat-key]"
    );
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
    this._syncDiscretePopoverAnchorAfterDiscreteRender();
    var afterDiscreteLay = function () {
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
    var fitBox = invertCell.querySelector(".cryo-cc-discrete-invert-fitbox");
    var label = invertCell.querySelector(".cryo-cc-discrete-switch-label");
    if (!fitBox || !label) return;
    if (
      invertCell.classList.contains("cryo-cc-discrete-cell--invert-heading")
      && invertCell.hidden
    ) {
      return;
    }
    /* Heading chip: sizing comes from CSS (auto width × title line-height); keep one readable line */
    if (invertCell.classList.contains("cryo-cc-discrete-cell--invert-heading")) {
      label.style.fontSize = "";
      label.style.whiteSpace = "nowrap";
      label.style.overflow = "";
      label.style.textOverflow = "";
      return;
    }
    /* Content box inside padding — clientWidth/Height include padding, but text lays out in the inner box */
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

    /* Pair plot / 3-D footer invert: one line, scale font to full phrase width × row height */
    if (invertCell.classList.contains("cryo-cc-discrete-cell--invert-footer")) {
      var slackF = 3;
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

    /* Multi-line: largest font so (a) every whole word fits inner width and (b) wrapped block fits box */
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
    var inputs = this.discreteSwitches.querySelectorAll("input[type=\"checkbox\"][data-cat-key]");
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
      var inputs = this.discreteSwitches.querySelectorAll("input[type=\"checkbox\"][data-cat-key]");
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
      var inputs = this.discreteSwitches.querySelectorAll("input[type=\"checkbox\"][data-cat-key]");
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
      var inputs = this.discreteSwitches.querySelectorAll("input[type=\"checkbox\"][data-cat-key]");
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
