/**
 * Volume GIF animations for the 3D volume landscapes page (nearest k-means vol per particle).
 * Reuses the same /api/landscape_volpca/* endpoints as the sketched landscape explorer.
 */
(function(global) {
  "use strict";

  function postJson(url, body) {
    return fetch(url, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body || {}),
    }).then(function(r) {
      return r.json().then(function(j) {
        if (!r.ok) throw new Error(j.error || r.status);
        return j;
      });
    });
  }

  function customDataRowAsArray(row) {
    if (row == null) return [];
    if (Array.isArray(row)) return row;
    if (typeof row === "object" && row.length !== undefined) {
      try {
        return Array.from(row);
      } catch (e) {
        return [];
      }
    }
    return [row];
  }

  function layoutHasVolAnim(gd) {
    var m = gd && gd.layout && gd.layout.meta;
    return !!(m && m.cdrgn_landscape_vol_animation);
  }

  function layoutMetaSketchCentroidCd(gd) {
    var m = gd && gd.layout && gd.layout.meta;
    return !!(m && m.cdrgn_landscape_sketch_centroid_cd);
  }

  // If Plotly drops/loses the `layout.meta` centroid flag for some reason, we still
  // want to be able to detect centroid particles via the 0/1 values carried in the
  // trace's `customdata` (second-last column for the volume-landscape plot).
  var sketchCentroidCdInferred = null; // null=unknown, bool=inferred

  function customdataForHighlight(gd, trace) {
    trace = trace || (gd && gd.data && gd.data[0]);
    if (!trace || !trace.x) return null;
    var n = trace.x.length;
    if (trace.customdata && trace.customdata.length === n) return trace.customdata;
    var full = gd && gd._fullData && gd._fullData[0];
    if (full && full.customdata && full.customdata.length === n) return full.customdata;
    return null;
  }

  function isSketchCentroidAtPointIndex(gd, trace, i) {
    var cd = customdataForHighlight(gd, trace);
    if (!cd || i < 0 || i >= cd.length || cd[i] == null) return false;
    var row = customDataRowAsArray(cd[i]);
    if (row.length < 2) return false;
    var centFlag = parseInt(String(row[row.length - 2]), 10);
    if (isNaN(centFlag)) return false;
    if (layoutMetaSketchCentroidCd(gd)) return centFlag === 1;

    // Meta flag isn't present; infer from a small customdata sample that the
    // second-last column is behaving like a 0/1 centroid flag.
    if (sketchCentroidCdInferred == null) {
      var maxJ = Math.min(cd.length, 100);
      var ok = true;
      for (var j = 0; j < maxJ; j++) {
        var rj = customDataRowAsArray(cd[j]);
        if (rj.length < 2) continue;
        var vj = parseInt(String(rj[rj.length - 2]), 10);
        // Some Plotly/serialisation paths can yield null-ish customdata entries;
        // those shouldn't poison the inference.
        if (!isNaN(vj) && (vj !== 0 && vj !== 1)) {
          ok = false;
          break;
        }
      }
      sketchCentroidCdInferred = ok;
    }
    if (!sketchCentroidCdInferred) return false;
    return centFlag === 1;
  }

  function volIdAtPointIndex(gd, trace, i) {
    if (trace.ids && trace.ids.length === trace.x.length && trace.ids[i] != null && trace.ids[i] !== "") {
      var vid = parseInt(String(trace.ids[i]), 10);
      if (!isNaN(vid)) return vid;
    }
    var cd = customdataForHighlight(gd, trace);
    if (cd && cd.length === trace.x.length && cd[i] != null) {
      var row = customDataRowAsArray(cd[i]);
      if (layoutHasVolAnim(gd) && row.length >= 2) {
        var last = row[row.length - 1];
        var vLast = parseInt(String(last), 10);
        if (!isNaN(vLast)) return vLast;
      }
      if (row.length) return parseInt(String(row[0]), 10);
    }
    return NaN;
  }

  function volFromPoint(gd, pt) {
    if (!pt) return NaN;
    if (pt.id != null && String(pt.id) !== "") {
      var idv = parseInt(String(pt.id), 10);
      if (!isNaN(idv)) return idv;
    }
    var trace = gd && gd.data && gd.data[0];
    if (trace && pt.pointNumber != null && pt.pointNumber >= 0) {
      var vi = volIdAtPointIndex(gd, trace, pt.pointNumber);
      if (!isNaN(vi)) return vi;
    }
    var row = customDataRowAsArray(pt.customdata);
    if (layoutHasVolAnim(gd) && row.length >= 2) {
      var v2 = parseInt(String(row[row.length - 1]), 10);
      if (!isNaN(v2)) return v2;
    }
    return row.length ? parseInt(String(row[0]), 10) : NaN;
  }

  function hex6ToRgba(hex, a) {
    var m = /^#([0-9a-fA-F]{6})$/.exec(hex);
    if (!m) return null;
    var x = parseInt(m[1], 16);
    var r = (x >> 16) & 255;
    var g = (x >> 8) & 255;
    var b = x & 255;
    return "rgba(" + r + "," + g + "," + b + "," + a + ")";
  }

  function isCircledMontageGlyphText(text) {
    var s = text == null ? "" : String(text);
    if (!s) return false;
    var code = s.charCodeAt(0);
    return code >= 0x24b6 && code <= 0x24cf;
  }

  function applyBadgeBackground(badge, color) {
    badge.style.color = "#1a1a1a";
    badge.style.backgroundColor = "";
    badge.classList.remove("volsketch-gif-overlay--on-tint");
    if (typeof color !== "string" || !color.trim()) {
      return;
    }
    var c = color.trim();
    var rgba = null;
    if (/^#([0-9a-fA-F]{6})$/.test(c)) {
      rgba = hex6ToRgba(c, 0.88);
    } else {
      var rm = /^rgba?\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)/i.exec(c);
      if (rm) rgba = "rgba(" + rm[1] + "," + rm[2] + "," + rm[3] + ",0.88)";
    }
    if (rgba) {
      badge.style.backgroundColor = rgba;
      badge.classList.add("volsketch-gif-overlay--on-tint");
    }
  }

  /** Montage letter on GIF previews — glyph fill matches plot covariate colour (not covariate side badge). */
  function applyBadgeLetterColor(badge, color) {
    badge.style.backgroundColor = "";
    badge.classList.remove("volsketch-gif-overlay--on-tint");
    badge.classList.remove("volsketch-gif-overlay--circle-chip");
    badge.classList.remove("volsketch-gif-overlay--letter-cov-glyph");
    if (typeof color !== "string" || !color.trim()) {
      badge.style.color = "#1a1a1a";
      return;
    }
    var c = color.trim();
    var label = badge.textContent || "";
    if (isCircledMontageGlyphText(label)) {
      var rgba = null;
      if (/^#([0-9a-fA-F]{6})$/.test(c)) {
        rgba = hex6ToRgba(c, 0.95);
      } else {
        var rm = /^rgba?\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)/i.exec(c);
        if (rm) rgba = "rgba(" + rm[1] + "," + rm[2] + "," + rm[3] + ",0.95)";
        else rgba = c;
      }
      badge.style.color = rgba;
      return;
    }
    badge.classList.add("volsketch-gif-overlay--letter-cov-glyph");
    badge.style.color = c;
  }

  function resetBadgeChrome(badge) {
    badge.style.backgroundColor = "";
    badge.style.color = "#1a1a1a";
    badge.classList.remove("volsketch-gif-overlay--on-tint");
    badge.classList.remove("volsketch-gif-overlay--circle-chip");
    badge.classList.remove("volsketch-gif-overlay--letter-cov-glyph");
  }

  function mountPreviewOverlay(wrap, badge, covBadge, covVarEl, covValEl, img, spec) {
    wrap._l3dvaCancelOverlay = null;
    function hideCov() {
      covVarEl.textContent = "";
      covValEl.textContent = "";
      covVarEl.hidden = true;
      covBadge.hidden = true;
      resetBadgeChrome(covBadge);
    }
    if (!spec || !spec.style) {
      badge.textContent = "";
      badge.hidden = true;
      resetBadgeChrome(badge);
      hideCov();
      return;
    }
    if (spec.style === "static") {
      badge.textContent = spec.text || "";
      badge.hidden = !badge.textContent;
      applyBadgeLetterColor(badge, spec.badge_background);
      wrap._l3dvaReapplyLetterColor = function() {
        applyBadgeLetterColor(badge, spec.badge_background || "");
      };
      var ctext = spec.covariate_text != null ? String(spec.covariate_text).trim() : "";
      if (ctext) {
        var vLabStatic = (
          spec.covariate_variable_label != null
            ? String(spec.covariate_variable_label).trim()
            : ""
        );
        covVarEl.textContent = vLabStatic;
        covVarEl.hidden = !vLabStatic;
        covValEl.textContent = ctext;
        covBadge.hidden = false;
        applyBadgeBackground(
          covBadge,
          spec.covariate_badge_background != null ? spec.covariate_badge_background : ""
        );
      } else {
        hideCov();
      }
      return;
    }
    if (spec.style === "rotate_frames") {
      var frameTexts = spec.frame_covariate_texts || [];
      var frameBgs = spec.frame_badge_backgrounds || [];
      var covVarRotate = (
        spec.covariate_variable_label != null
          ? String(spec.covariate_variable_label).trim()
          : ""
      );
      var fdRot = Math.max(1, parseInt(spec.frame_duration_ms, 10) || 100);
      var totalRot = Math.max(
        1,
        parseInt(spec.total_frames, 10) || frameTexts.length || 1
      );
      badge.textContent = spec.text || "";
      badge.hidden = !badge.textContent;
      applyBadgeLetterColor(badge, spec.badge_background || "");
      wrap._l3dvaReapplyLetterColor = function() {
        applyBadgeLetterColor(badge, spec.badge_background || "");
      };
      if (!frameTexts.length) {
        hideCov();
        return;
      }
      var rafRot = 0;
      var bootRot = 0;
      var startRot = null;
      function applyRotateFrameIdx(fi) {
        fi = Math.max(0, Math.min(frameTexts.length - 1, fi | 0));
        var covt = frameTexts[fi] != null ? String(frameTexts[fi]).trim() : "";
        var bg = frameBgs[fi] != null ? frameBgs[fi] : "";
        covVarEl.textContent = covVarRotate;
        covVarEl.hidden = !covVarRotate;
        if (covt) {
          covValEl.textContent = covt;
          covBadge.hidden = false;
          applyBadgeBackground(covBadge, bg);
        } else {
          hideCov();
        }
      }
      function tickRotate(now) {
        if (!wrap.isConnected) {
          if (rafRot) cancelAnimationFrame(rafRot);
          rafRot = 0;
          return;
        }
        if (startRot === null) startRot = now;
        var elapsed = now - startRot;
        var frameIdx = Math.floor(elapsed / fdRot) % totalRot;
        var covIdx = Math.min(
          frameTexts.length - 1,
          Math.floor((frameIdx * frameTexts.length) / totalRot)
        );
        applyRotateFrameIdx(covIdx);
        rafRot = requestAnimationFrame(tickRotate);
      }
      function beginRotateClock() {
        applyRotateFrameIdx(0);
        startRot = null;
        bootRot = requestAnimationFrame(function() {
          bootRot = 0;
          rafRot = requestAnimationFrame(tickRotate);
        });
      }
      wrap._l3dvaCancelOverlay = function() {
        if (bootRot) cancelAnimationFrame(bootRot);
        if (rafRot) cancelAnimationFrame(rafRot);
        bootRot = 0;
        rafRot = 0;
      };
      var decRot = img.decode && img.decode();
      if (decRot && typeof decRot.then === "function") {
        decRot.then(beginRotateClock).catch(beginRotateClock);
      } else if (img.complete) {
        beginRotateClock();
      } else {
        img.addEventListener("load", beginRotateClock, { once: true });
      }
      return;
    }
    if (spec.style !== "cycle_segments") {
      badge.hidden = true;
      resetBadgeChrome(badge);
      hideCov();
      return;
    }
    var labels = spec.segment_labels || [];
    var segBgs = spec.segment_backgrounds || [];
    var segCovBgs = spec.segment_covariate_backgrounds || [];
    var segCov = spec.segment_covariate_texts || [];
    var covVarLabelCycle = (
      spec.covariate_variable_label != null
        ? String(spec.covariate_variable_label).trim()
        : ""
    );
    var fpv = Math.max(1, parseInt(spec.frames_per_segment, 10) || 8);
    var fd = Math.max(1, parseInt(spec.frame_duration_ms, 10) || 100);
    if (!labels.length) {
      badge.hidden = true;
      resetBadgeChrome(badge);
      hideCov();
      return;
    }
    badge.hidden = false;
    var totalF = Math.max(1, parseInt(spec.total_frames, 10) || labels.length * fpv);
    var rafId = 0;
    var bootRaf = 0;
    var startT = null;
    function applySegmentIdx(idx) {
      idx = Math.max(0, Math.min(labels.length - 1, idx | 0));
      wrap._l3dvaCurrentSegIdx = idx;
      badge.textContent = labels[idx] || "";
      var bg = segBgs[idx] != null ? segBgs[idx] : "";
      applyBadgeLetterColor(badge, bg);
      var covt = "";
      if (segCov.length > idx && segCov[idx] != null) {
        covt = String(segCov[idx]).trim();
      }
      covVarEl.textContent = covVarLabelCycle;
      covVarEl.hidden = !covVarLabelCycle;
      if (covt) {
        covValEl.textContent = covt;
        covBadge.hidden = false;
        var covBg = segCovBgs.length > idx && segCovBgs[idx] != null ? segCovBgs[idx] : bg;
        applyBadgeBackground(covBadge, covBg);
      } else {
        covVarEl.hidden = true;
        covVarEl.textContent = "";
        covValEl.textContent = "";
        covBadge.hidden = true;
        resetBadgeChrome(covBadge);
      }
    }
    function tick(now) {
      if (!wrap.isConnected) {
        if (rafId) cancelAnimationFrame(rafId);
        rafId = 0;
        return;
      }
      if (startT === null) startT = now;
      var elapsed = now - startT;
      var frameIdx = Math.floor(elapsed / fd) % totalF;
      var idx = Math.min(labels.length - 1, Math.floor(frameIdx / fpv));
      applySegmentIdx(idx);
      rafId = requestAnimationFrame(tick);
    }
    wrap._l3dvaReapplyLetterColor = function() {
      applySegmentIdx(wrap._l3dvaCurrentSegIdx != null ? wrap._l3dvaCurrentSegIdx : 0);
    };
    function beginClock() {
      applySegmentIdx(0);
      startT = null;
      bootRaf = requestAnimationFrame(function() {
        bootRaf = 0;
        rafId = requestAnimationFrame(tick);
      });
    }
    wrap._l3dvaCancelOverlay = function() {
      if (bootRaf) cancelAnimationFrame(bootRaf);
      if (rafId) cancelAnimationFrame(rafId);
      bootRaf = 0;
      rafId = 0;
    };
    var dec = img.decode && img.decode();
    if (dec && typeof dec.then === "function") {
      dec.then(beginClock).catch(beginClock);
    } else if (img.complete) {
      beginClock();
    } else {
      img.addEventListener("load", beginClock, { once: true });
    }
  }

  var SAFE_LETTERS = (function() {
    var out = [];
    for (var c = 65; c <= 90; c++) {
      if (c !== 73 && c !== 79 && c !== 85) out.push(String.fromCharCode(c));
    }
    return out;
  })();

  function montageLabelAt(idx) {
    if (idx < SAFE_LETTERS.length) return SAFE_LETTERS[idx];
    var j = idx - SAFE_LETTERS.length;
    return SAFE_LETTERS[Math.floor(j / SAFE_LETTERS.length)]
      + SAFE_LETTERS[j % SAFE_LETTERS.length];
  }

  // Legacy helper: scatter3d WebGL text doesn't reliably support HTML tags like <b>/<i>,
  // but some dashboard wiring tests assert the presence of this formatter.
  function formatVolMontageLetterText(letter, isBold, isItalic) {
    var s = letter == null ? "" : String(letter);
    if (!s) return "";
    if (isBold && isItalic) return "<b><i>" + s + "</i></b>";
    if (isBold) return "<b>" + s + "</b>";
    if (isItalic) return "<i>" + s + "</i>";
    return s;
  }

  function volFilenameIndex(v) {
    var n = typeof v === "number" ? v : parseInt(String(v), 10);
    if (isNaN(n)) return String(v);
    return String(n).padStart(3, "0");
  }

  var GIF_MODE_DESC = {
    cycle: (
      "Single GIF, fixed camera: steps through your selected k-means volumes (nearest PCA match per particle). "
      + "If you select more than fit, up to <strong>50</strong> are chosen at random for the cycle."
    ),
    rotate_each: (
      "One rotating GIF per selected volume. "
      + "If you select more than fit, up to <strong>5</strong> volumes are chosen at random."
    ),
  };

  var P3S = global.CryoPlotlyScatter3dScene;
  if (!P3S || typeof P3S.snapshot !== "function") {
    throw new Error(
      "CryoPlotlyScatter3dScene is required — include plotly_scatter3d_scene.js before latent3d_landscape_vol_animations.js"
    );
  }

  function sceneSnapHasCamera(snap) {
    if (!snap || !snap.camera || typeof snap.camera !== "object") return false;
    var eye = snap.camera.eye;
    return !!(eye && typeof eye === "object");
  }

  /** Merge pointerdown snapshot (pre Plotly double-click camera reset) with a live read. */
  function mergeSnapsPreferPointerdown(pd, live) {
    if (!pd || typeof pd !== "object") return live;
    if (!live || typeof live !== "object") return pd;
    return {
      camera: pd.camera || live.camera,
      xrange: pd.xrange || live.xrange,
      yrange: pd.yrange || live.yrange,
      zrange: pd.zrange || live.zrange,
    };
  }

  function _hexToRgb(hex) {
    var m = /^#([0-9a-f]{3}|[0-9a-f]{6})$/i.exec(String(hex || "").trim());
    if (!m) return null;
    var s = m[1];
    if (s.length === 3) {
      return {
        r: parseInt(s[0] + s[0], 16),
        g: parseInt(s[1] + s[1], 16),
        b: parseInt(s[2] + s[2], 16),
      };
    }
    return {
      r: parseInt(s.slice(0, 2), 16),
      g: parseInt(s.slice(2, 4), 16),
      b: parseInt(s.slice(4, 6), 16),
    };
  }

  function _parseRgbLike(s) {
    var m = /^rgba?\(\s*(\d+)[\s,]+(\d+)[\s,]+(\d+)/i.exec(String(s || "").trim());
    if (!m) return null;
    return { r: +m[1], g: +m[2], b: +m[3] };
  }

  /** WCAG relative luminance in [0, 1]; unknown colours → 0.5 */
  function _fillRelativeLuminance(css) {
    var rgb = _hexToRgb(css) || _parseRgbLike(css);
    if (!rgb) return 0.5;
    function lin(u) {
      u /= 255;
      return u <= 0.03928 ? u / 12.92 : Math.pow((u + 0.055) / 1.055, 2.4);
    }
    var R = lin(rgb.r);
    var G = lin(rgb.g);
    var B = lin(rgb.b);
    return 0.2126 * R + 0.7152 * G + 0.0722 * B;
  }

  function plotlyColorToCss(c) {
    if (c == null) return null;
    if (typeof c === "string") return c;
    if (typeof c === "object") {
      if (c.r != null && c.g != null && c.b != null) {
        var a = c.a != null && isFinite(c.a) ? c.a : 1;
        if (a >= 0.999) {
          return "rgb(" + Math.round(c.r) + "," + Math.round(c.g) + "," + Math.round(c.b) + ")";
        }
        return "rgba(" + Math.round(c.r) + "," + Math.round(c.g) + "," + Math.round(c.b) + "," + a + ")";
      }
    }
    return null;
  }

  function plotColorModeFromGd(gd) {
    try {
      var m = gd && gd.layout && gd.layout.meta;
      if (!m) return null;
      if (m.cdrgn_color_mode) return String(m.cdrgn_color_mode);
      var leg = m.cdrgn_color_legend;
      if (leg && leg.type) return String(leg.type);
    } catch (eMode) { /* ignore */ }
    return null;
  }

  function discreteColorLegendMap(gd) {
    var out = Object.create(null);
    try {
      var leg = gd && gd.layout && gd.layout.meta && gd.layout.meta.cdrgn_color_legend;
      if (!leg || leg.type !== "discrete" || !leg.items) return out;
      for (var li = 0; li < leg.items.length; li++) {
        var it = leg.items[li];
        if (!it) continue;
        var key = it.label == null ? "" : String(it.label);
        if (it.color) out[key] = String(it.color);
      }
    } catch (eMap) { /* ignore */ }
    return out;
  }

  function colorCovariateKeyFromCustomdataRow(row) {
    var r = customDataRowAsArray(row);
    return r.length >= 3 && r[2] != null ? String(r[2]) : "";
  }

  function continuousCssFromValue(val, trace, gd, paletteName) {
    if (typeof val !== "number" || !isFinite(val)) return null;
    var mk = (trace && trace.marker) || {};
    var cmin = mk.cmin;
    var cmax = mk.cmax;
    try {
      var leg = gd && gd.layout && gd.layout.meta && gd.layout.meta.cdrgn_color_legend;
      if (leg && leg.type === "continuous") {
        if (cmin == null || !isFinite(cmin)) cmin = leg.min;
        if (cmax == null || !isFinite(cmax)) cmax = leg.max;
      }
    } catch (eLeg) { /* ignore */ }
    cmin = typeof cmin === "number" && isFinite(cmin) ? cmin : 0;
    cmax = typeof cmax === "number" && isFinite(cmax) ? cmax : 1;
    var span = cmax - cmin;
    var t = span > 0 ? (val - cmin) / span : 0.5;
    t = Math.max(0, Math.min(1, t));
    var pal = paletteName || "Viridis";
    try {
      var metaPal = gd && gd.layout && gd.layout.meta && gd.layout.meta.cdrgn_continuous_palette;
      if (metaPal && global.CryoColorCovariateLegend
          && CryoColorCovariateLegend.hasPalette(metaPal)) {
        pal = String(metaPal);
      }
    } catch (ePal) { /* ignore */ }
    var CCL = global.CryoColorCovariateLegend;
    if (CCL && typeof CCL.paletteScaleCSS === "function" && CCL.hasPalette(pal)) {
      return CCL.paletteScaleCSS(pal, t);
    }
    return null;
  }

  /**
   * Resolve plotted marker fill for trace index ``i`` — discrete hex, continuous colorscale,
   * or Plotly-expanded colours on ``_fullData``.
   */
  function markerFillCssAtIndex(trace, gd, i, paletteName) {
    if (!trace || !gd || i < 0) return "#4a5568";
    var full = gd._fullData && gd._fullData[0];
    var expanded = full && full.marker && full.marker.color;
    if (Array.isArray(expanded) && expanded[i] != null) {
      var exCss = plotlyColorToCss(expanded[i]);
      if (exCss) return exCss;
    }
    var mk = (full && full.marker) || trace.marker || {};
    var mc = mk.color;
    if (typeof mc === "string") return mc;
    if (Array.isArray(mc) && mc[i] != null) {
      var mcCss = plotlyColorToCss(mc[i]);
      if (mcCss) return mcCss;
      if (typeof mc[i] === "string") return mc[i];
    }
    var mode = plotColorModeFromGd(gd);
    var cd = customdataForHighlight(gd, trace);
    var row = cd && cd[i] != null ? cd[i] : null;
    if (mode === "discrete") {
      var fk = colorCovariateKeyFromCustomdataRow(row);
      if (fk) {
        var dmap = discreteColorLegendMap(gd);
        if (dmap[fk]) return dmap[fk];
      }
      if (Array.isArray(mc) && mc[i] != null) {
        var dCss = plotlyColorToCss(mc[i]);
        if (dCss) return dCss;
        if (typeof mc[i] === "string") return mc[i];
      }
    }
    if (mode === "continuous") {
      var cKey = colorCovariateKeyFromCustomdataRow(row);
      var cVal = cKey !== "" ? parseFloat(cKey) : NaN;
      if (!isNaN(cVal)) {
        var fromCd = continuousCssFromValue(cVal, trace, gd, paletteName);
        if (fromCd) return fromCd;
      }
      if (Array.isArray(mc) && mc[i] != null && typeof mc[i] === "number" && isFinite(mc[i])) {
        var fromMc = continuousCssFromValue(mc[i], trace, gd, paletteName);
        if (fromMc) return fromMc;
      }
    }
    if (Array.isArray(mc) && mc[i] != null) {
      var anyCss = plotlyColorToCss(mc[i]);
      if (anyCss) return anyCss;
      if (typeof mc[i] === "string") return mc[i];
    }
    return "#4a5568";
  }

  /** Dark fills → light outline; light fills → dark outline (selection ring). */
  function selectionOutlineColorForFill(fillCss) {
    return _fillRelativeLuminance(fillCss) < 0.34 ? "#f1f5f9" : "#0f172a";
  }

  function montagePlotLetterStrokeShadow(strokeCss) {
    var s = strokeCss || "#000";
    return (
      "0 0 1px " + s + ",-1px 0 1px " + s + ",1px 0 1px " + s + ","
      + "0 -1px 1px " + s + ",0 1px 1px " + s
    );
  }

  /** Plot letter border: white on dark covariate glyphs, black on light glyphs. */
  function montagePlotLetterStrokeForFill(fillCss) {
    var strokeCss = _fillRelativeLuminance(fillCss) < 0.34 ? "#ffffff" : "#000000";
    return montagePlotLetterStrokeShadow(strokeCss);
  }

  // 3D plot scene letters only (GIF preview size is separate CSS).
  var VOL_MONTAGE_PLOT_LETTER_PX = 36 * 0.8 * 0.75 * 1.5 * 0.8 * 0.8 * 0.8 * 0.8;
  /** Non-centroid, non–last-dblclick-anchor montage letters on duplicate particles per volume. */
  var VOL_MONTAGE_PLOT_LETTER_PX_SECONDARY = VOL_MONTAGE_PLOT_LETTER_PX * 0.8;
  var VOL_MONTAGE_LETTER_FONT = {
    size: VOL_MONTAGE_PLOT_LETTER_PX,
    color: "#1a1a1a",
    family: "system-ui, Segoe UI, sans-serif",
  };

  // For sketch-centroid particles: circle the montage letter glyph.
  // (We generate circled letters at runtime to keep the source ASCII-only.)
  function circledMontageLabelText(label) {
    var s = label == null ? "" : String(label);
    if (!s) return "";
    var out = "";
    for (var j = 0; j < s.length; j++) {
      var ch = s[j];
      var code = ch.charCodeAt(0);
      // Unicode: CIRCLED LATIN CAPITAL LETTER A starts at U+24B6.
      if (code >= 65 && code <= 90) {
        out += String.fromCharCode(0x24B6 + (code - 65));
      } else {
        out += ch;
      }
    }
    return out;
  }

  /**
   * Bold / italic montage letters on scatter3d: WebGL trace ``text`` ignores HTML and per-point
   * ``textfont.weight`` / ``style``; use ``layout.scene.annotations`` at anchor points instead.
   */
  function scatterMarkerOpacityFromGd(gd) {
    try {
      var tr = gd && gd.data && gd.data[0];
      var mo = tr && tr.marker && tr.marker.opacity;
      if (typeof mo === "number" && isFinite(mo)) return mo;
      if (Array.isArray(mo) && mo.length && typeof mo[0] === "number" && isFinite(mo[0])) {
        return mo[0];
      }
    } catch (eMo) { /* ignore */ }
    return null;
  }

  function volMontageSceneAnnotation(x, y, z, letter, isBold, isItalic, opts) {
    if (!letter) return null;
    opts = opts || {};
    var fontSizePx = typeof opts.fontSizePx === "number" && isFinite(opts.fontSizePx)
      ? opts.fontSizePx
      : VOL_MONTAGE_LETTER_FONT.size;
    var fontColorCss = typeof opts.fontColorCss === "string" && opts.fontColorCss.trim()
      ? opts.fontColorCss
      : VOL_MONTAGE_LETTER_FONT.color;
    var scatterOp = typeof opts.scatterOpacity === "number" && isFinite(opts.scatterOpacity)
      ? opts.scatterOpacity
      : scatterMarkerOpacityFromGd(opts.gd);
    var annOpacity = scatterOp != null
      ? Math.max(0.0, Math.min(1.0, (1 + scatterOp) / 2))
      : 0.5;
    var font = {
      size: fontSizePx,
      color: fontColorCss,
      family: VOL_MONTAGE_LETTER_FONT.family,
    };
    if (isBold) font.weight = "bold";
    if (isItalic) font.style = "italic";
    font.shadow = montagePlotLetterStrokeForFill(fontColorCss);
    return {
      x: x,
      y: y,
      z: z,
      text: String(letter),
      opacity: annOpacity,
      showarrow: false,
      xanchor: "center",
      // Move montage letters closer to the anchor point.
      yanchor: "middle",
      // Scene annotations can end up sitting on top of scatter markers.
      // Disable event capture so Plotly click/dblclick keeps working.
      captureevents: false,
      font: font,
    };
  }

  function boot(cfg) {
    var gd = cfg.gd;
    var getColorMode = cfg.getColorMode;
    var getPalette = cfg.getPalette;
    var onAfterPlot = cfg.onAfterPlot;
    var setSelectionRendering = typeof cfg.setSelectionRendering === "function"
      ? cfg.setSelectionRendering
      : null;
    var setPreviewAnimLoading = typeof cfg.setPreviewAnimLoading === "function"
      ? cfg.setPreviewAnimLoading
      : null;

    var PREFIX = "l3dva-";
    function $(id) {
      return document.getElementById(PREFIX + id);
    }

    var animateStatusEl = $("animate-status");
    var animateProgressEl = $("animate-progress");
    var selSummaryEl = $("sel-summary");
    var selIndicesDetails = $("sel-indices-details");
    var selIndicesSummaryLabel = $("sel-indices-summary-label");
    var selIndicesList = $("sel-indices-list");
    var gifModeDescEl = $("mode-desc");
    var gifModeSwitch = $("gif-mode-switch");
    var gifFrames = $("gif-frames");
    var gifFramesWrap = $("gif-frames-wrap");
    var previewGrid = $("preview-grid");
    var randomSelBtn = $("random-sel");
    var saveGifBtn = $("save-gif");
    var viewRotateRow = document.querySelector(".l3dva-view-rotate-row");
    var viewRotateAngleEl = $("view-rotate-angle");
    var viewRotateSummaryEl = $("view-rotate-summary");
    var viewRotateBtns = Array.from(document.querySelectorAll("[data-l3dva-view-axis]"));

    var META = null;
    var selectedVols = new Set();
    /** Volume ids among the current selection that contain a sketch-centroid particle. */
    var centroidVolsInSelection = new Set();
    /** Stable montage letter per k-means volume index while it stays selected (not re-sorted on add). */
    var volMontageLabel = {};
    var nextMontageLabelIdx = 0;
    var lastAnimToken = null;
    var lastBatchMode = null;
    var clickTimer = null;
    var clickLastVol = null;
    /** Per sketch volume: scatter point index last double-clicked to toggle that volume (bold letter). */
    var volDblClickPointNum = {};
    /** Last volume toggled via double-click (bold falls back to its centroid after ``Plotly.react``). */
    var lastDblClickVol = null;
    var debounceGifTimer = null;
    var gifReqGen = 0;
    var debouncedGifReason = "selection";
    var viewRotations = { x: 0, y: 0, z: 0 };
    var lastChimeraxViewMatrix = "";
    var chimeraxViewMatrixUnavailable = false;
    var landscapeAnimInFlight = false;
    /** Plotly only defines ``gd.on`` after the graph div has been drawn once. */
    var plotEventsBound = false;
    /** Scene snapshot from last ``pointerdown`` (capture); survives Plotly’s double-click camera reset. */
    var volAnimPointerdownSnap = null;
    var volAnimPointerdownAt = 0;
    var volAnimSelOverlayDepth = 0;
    /** Count of selection ``Plotly.update``/``restyle`` chains still running (only when overlay depth is on). */
    var volAnimSelHighlightRestylesInFlight = 0;
    /** Ignore extra Plotly click events right after our double-click toggle. */
    var ignorePlotlyClickUntilMs = 0;
    /** Prevent accidental double toggles (toggle on then immediately off). */
    var lastToggleVol = null;
    var lastToggleAtMs = 0;
    /** Show frosted overlay once per selection restyle (idempotent if already active). */
    function beginVolAnimSelectionOverlay() {
      if (!setSelectionRendering) return;
      if (volAnimSelOverlayDepth >= 1) return;
      volAnimSelOverlayDepth = 1;
      try {
        setSelectionRendering(true);
      } catch (eO) { /* ignore */ }
    }
    function endVolAnimSelectionOverlay() {
      if (!setSelectionRendering) return;
      if (volAnimSelOverlayDepth < 1) return;
      volAnimSelOverlayDepth = 0;
      try {
        setSelectionRendering(false);
      } catch (eO) { /* ignore */ }
    }

    function safeSetPreviewAnimLoading(on) {
      if (!setPreviewAnimLoading) return;
      try {
        setPreviewAnimLoading(!!on);
      } catch (eL) { /* ignore */ }
    }

    function teardownPreviews() {
      if (!previewGrid) return;
      previewGrid.querySelectorAll(".volsketch-anim-wrap").forEach(function(w) {
        if (w._l3dvaCancelOverlay) {
          w._l3dvaCancelOverlay();
          w._l3dvaCancelOverlay = null;
        }
      });
    }

    function clearPreviewGrid() {
      teardownPreviews();
      if (previewGrid) previewGrid.innerHTML = "";
      syncViewRotationControls();
    }

    function currentGifMode() {
      if (!gifModeSwitch) return "cycle";
      return gifModeSwitch.checked ? "rotate_each" : "cycle";
    }

    function syncGifModeDescription() {
      if (!gifModeDescEl) return;
      var m = currentGifMode();
      gifModeDescEl.innerHTML = GIF_MODE_DESC[m] || GIF_MODE_DESC.cycle;
    }

    function syncPreviewGridLayout() {
      if (!previewGrid) return;
      var cycle = currentGifMode() === "cycle";
      previewGrid.classList.toggle("volsketch-preview-grid--cycle", cycle);
    }

    function syncGifFramesVisibility() {
      if (!gifFramesWrap) return;
      var m = currentGifMode();
      gifFramesWrap.style.display = m === "rotate_each" ? "block" : "none";
      syncPreviewGridLayout();
    }
    syncGifModeDescription();
    syncGifFramesVisibility();

    function syncSaveGifButton() {
      if (!saveGifBtn) return;
      var hasPreviews = previewGrid && previewGrid.children.length > 0;
      var canSave = !!lastAnimToken && hasPreviews && !landscapeAnimInFlight;
      saveGifBtn.disabled = !canSave;
      if (canSave) {
        saveGifBtn.removeAttribute("title");
      } else if (landscapeAnimInFlight) {
        saveGifBtn.title = "Available after ChimeraX finishes rendering the current previews.";
      } else {
        saveGifBtn.title = "Select volumes and wait for GIF previews to finish rendering.";
      }
      syncViewRotationControls();
    }

    function syncViewRotationControls() {
      var hasLoadedPreview = !!(
        lastAnimToken
        && previewGrid
        && previewGrid.children.length > 0
        && !landscapeAnimInFlight
      );
      var hint = hasLoadedPreview
        ? "Rotate the loaded animation preview by the entered angle."
        : "Load an animation preview before rotating the view.";
      if (viewRotateRow) {
        viewRotateRow.classList.toggle("is-disabled", !hasLoadedPreview);
        viewRotateRow.title = hint;
      }
      if (viewRotateAngleEl) {
        viewRotateAngleEl.disabled = !hasLoadedPreview;
        viewRotateAngleEl.title = hint;
      }
      viewRotateBtns.forEach(function(btn) {
        btn.disabled = !hasLoadedPreview;
        btn.title = hint;
      });
    }

    function normalizeViewDegrees(degrees) {
      var d = Number(degrees);
      if (!isFinite(d)) return 0;
      d = d % 360;
      if (Math.abs(d) < 1e-9) return 0;
      return d;
    }

    function currentViewRotationPayload() {
      return {
        x: normalizeViewDegrees(viewRotations.x),
        y: normalizeViewDegrees(viewRotations.y),
        z: normalizeViewDegrees(viewRotations.z),
      };
    }

    function viewRotationSummary() {
      var r = currentViewRotationPayload();
      return "requested view turns X " + r.x.toFixed(1)
        + "°, Y " + r.y.toFixed(1)
        + "°, Z " + r.z.toFixed(1) + "°";
    }

    function syncViewRotationSummary() {
      if (!viewRotateSummaryEl) return;
      if (lastChimeraxViewMatrix) {
        viewRotateSummaryEl.textContent = "ChimeraX view matrix:\n" + lastChimeraxViewMatrix;
      } else if (chimeraxViewMatrixUnavailable) {
        viewRotateSummaryEl.textContent = "ChimeraX view matrix: not reported by ChimeraX for this render.";
      } else if (landscapeAnimInFlight) {
        viewRotateSummaryEl.textContent = "ChimeraX view matrix: rendering…";
      } else {
        viewRotateSummaryEl.textContent = "ChimeraX view matrix: available after animation loads.";
      }
    }

    function setAnimateStatus(msg, showProgress, isErr) {
      if (!animateStatusEl) return;
      animateStatusEl.textContent = msg || "";
      animateStatusEl.style.color = isErr ? "var(--error, #b42318)" : "";
      if (animateProgressEl) {
        animateProgressEl.hidden = !showProgress || !!isErr;
        if (showProgress && !isErr) {
          animateProgressEl.setAttribute("aria-busy", "true");
        } else {
          animateProgressEl.removeAttribute("aria-busy");
        }
      }
    }

    function maxSelectable() {
      return META && META.n_volumes ? META.n_volumes : 10000;
    }

    function syncRandomSelButton() {
      if (!randomSelBtn || !META || META.chimerax_cpus == null) return;
      var n = Math.max(1, Number(META.chimerax_cpus));
      randomSelBtn.textContent = "Choose " + n + " volume" + (n === 1 ? "" : "s") + " at random";
      randomSelBtn.disabled = false;
    }

    function applyViewRotation(axis) {
      var deg = viewRotateAngleEl ? Number(viewRotateAngleEl.value) : NaN;
      if (!isFinite(deg)) {
        setAnimateStatus("Enter a finite rotation angle in degrees.", false, true);
        return;
      }
      viewRotations[axis] = normalizeViewDegrees(Number(viewRotations[axis] || 0) + deg);
      syncViewRotationSummary();
      setAnimateStatus("Updated " + viewRotationSummary() + ".", false, false);
      if (selectedVols.size > 0) {
        scheduleAutoGif("view_rotation");
      } else {
        syncSaveGifButton();
      }
    }

    function allPlotVolIds() {
      if (!gd || !gd.data || !gd.data[0]) return [];
      var trace = gd.data[0];
      var n = trace.x ? trace.x.length : 0;
      if (!n) return [];
      var s = new Set();
      for (var i = 0; i < n; i++) {
        var v = volIdAtPointIndex(gd, trace, i);
        if (!isNaN(v)) s.add(v);
      }
      return Array.from(s).sort(function(a, b) { return a - b; });
    }

    function shufflePick(arr, k) {
      var a = arr.slice();
      for (var i = a.length - 1; i > 0; i--) {
        var j = Math.floor(Math.random() * (i + 1));
        var t = a[i];
        a[i] = a[j];
        a[j] = t;
      }
      return a.slice(0, k);
    }

    function syncStableSelectionLabels() {
      if (!selectedVols.size) {
        volMontageLabel = {};
        nextMontageLabelIdx = 0;
        return;
      }
      var k;
      var toDrop = [];
      for (k in volMontageLabel) {
        if (!Object.prototype.hasOwnProperty.call(volMontageLabel, k)) continue;
        var vk = parseInt(k, 10);
        if (!selectedVols.has(vk)) toDrop.push(k);
      }
      for (var di = 0; di < toDrop.length; di++) {
        delete volMontageLabel[toDrop[di]];
      }
      selectedVols.forEach(function(v) {
        if (volMontageLabel[v] === undefined) {
          volMontageLabel[v] = montageLabelAt(nextMontageLabelIdx++);
        }
      });
    }

    function updateSelLabel() {
      syncStableSelectionLabels();
      var arr = Array.from(selectedVols).sort(function(a, b) { return a - b; });
      var nTotal = META && META.n_volumes != null ? Number(META.n_volumes) : NaN;
      if (selSummaryEl) {
        if (!isNaN(nTotal) && nTotal >= 0) {
          selSummaryEl.textContent = arr.length + "/" + nTotal + " k-means volumes selected (nearest PCA match)";
        } else {
          selSummaryEl.textContent = arr.length
            ? (arr.length + " k-means volumes selected")
            : "No k-means volumes selected";
        }
      }
      if (selIndicesDetails && selIndicesList && selIndicesSummaryLabel) {
        selIndicesList.innerHTML = "";
        if (!arr.length) {
          selIndicesDetails.hidden = true;
          selIndicesDetails.open = false;
        } else {
          selIndicesDetails.hidden = false;
          selIndicesSummaryLabel.textContent = (
            "Selected volume indices (" + arr.length + ")"
          );
          var listParts = [];
          arr.forEach(function(v) {
            listParts.push(volFilenameIndex(v) + " (" + volMontageLabel[v] + ")");
          });
          var listLabel = listParts.join(", ");
          selIndicesList.setAttribute(
            "aria-label",
            "Selected sketch volume indices: " + listLabel
          );
          arr.forEach(function(v) {
            var li = document.createElement("li");
            li.textContent = volFilenameIndex(v) + " (" + volMontageLabel[v] + ")";
            selIndicesList.appendChild(li);
          });
        }
      }
    }

    /**
     * Unselected-point marker diameter from the live trace (Plotly may expand a scalar ``size``
     * on ``_fullData``). After a prior highlight pass, sizes are per-point — use the minimum so
     * we do not clobber the server / ``scatter3d_z_json`` glyph scale when restyling overlays.
     */
    function referenceScatter3dBaseMarkerSize(trace, graphDiv) {
      var full = graphDiv && graphDiv._fullData && graphDiv._fullData[0];
      var marker = (full && full.marker) || (trace && trace.marker);
      if (!marker || marker.size == null) return null;
      var sz = marker.size;
      if (typeof sz === "number" && isFinite(sz) && sz > 0) return sz;
      if (Array.isArray(sz) && sz.length) {
        var minS = Infinity;
        var si;
        for (si = 0; si < sz.length; si++) {
          var one = sz[si];
          if (typeof one === "number" && isFinite(one) && one > 0 && one < minS) minS = one;
        }
        if (minS < Infinity) return minS;
      }
      return null;
    }

    function sketchCentroidPointIndexForVol(gd, trace, vol) {
      var n = trace.x ? trace.x.length : 0;
      for (var ci = 0; ci < n; ci++) {
        if (volIdAtPointIndex(gd, trace, ci) !== vol) continue;
        if (isSketchCentroidAtPointIndex(gd, trace, ci)) return ci;
      }
      return null;
    }

    function pruneVolDblClickAnchors(gd, trace) {
      if (!trace || !trace.x) return;
      var toDrop = [];
      var k;
      for (k in volDblClickPointNum) {
        if (!Object.prototype.hasOwnProperty.call(volDblClickPointNum, k)) continue;
        var v = parseInt(k, 10);
        var pt = volDblClickPointNum[k];
        if (
          isNaN(v)
          || pt == null
          || pt < 0
          || pt >= trace.x.length
          || volIdAtPointIndex(gd, trace, pt) !== v
        ) {
          toDrop.push(k);
        }
      }
      for (var di = 0; di < toDrop.length; di++) {
        delete volDblClickPointNum[toDrop[di]];
      }
    }

    /**
     * Bold anchor:
     * - valid double-click point index (stored in ``volDblClickPointNum``)
     * - else centroid for the *last* double-clicked volume (after Plotly camera/scene churn)
     * - else ``null`` (centroid becomes italic for other selected volumes)
     */
    function resolvePlotPalette() {
      if (getPalette) {
        var p = getPalette();
        if (p && global.CryoColorCovariateLegend && CryoColorCovariateLegend.hasPalette(p)) {
          return String(p);
        }
      }
      return "Viridis";
    }

    function colorCovariateActive() {
      var cm = getColorMode ? String(getColorMode()) : "none";
      if (cm && cm !== "none") return true;
      var mode = plotColorModeFromGd(gd);
      return !!(mode && mode !== "none");
    }

    function markerFillCssForVol(gd, trace, vol) {
      var v = typeof vol === "number" ? vol : parseInt(String(vol), 10);
      if (isNaN(v) || !trace || !trace.x) return null;
      var pt = sketchCentroidPointIndexForVol(gd, trace, v);
      if (pt == null) pt = boldAnchorPointForVol(gd, trace, v);
      if (pt == null) {
        for (var fi = 0; fi < trace.x.length; fi++) {
          if (volIdAtPointIndex(gd, trace, fi) === v) {
            pt = fi;
            break;
          }
        }
      }
      if (pt == null) return null;
      return markerFillCssAtIndex(trace, gd, pt, resolvePlotPalette());
    }

    function montageLetterFontColor(trace, gd, pointIndex) {
      if (!colorCovariateActive()) return VOL_MONTAGE_LETTER_FONT.color;
      return markerFillCssAtIndex(trace, gd, pointIndex, resolvePlotPalette())
        || VOL_MONTAGE_LETTER_FONT.color;
    }

    function refreshPreviewOverlayLetterColors() {
      if (!previewGrid || !gd || !gd.data || !gd.data[0]) return;
      var trace = gd.data[0];
      var active = colorCovariateActive();
      previewGrid.querySelectorAll(".volsketch-anim-wrap").forEach(function(wrap) {
        var badge = wrap.querySelector(".volsketch-gif-overlay:not(.volsketch-gif-overlay--cov)");
        var spec = wrap._l3dvaOverlaySpec;
        if (!badge || !spec) return;
        if (!active) {
          if (spec.style === "static" || spec.style === "rotate_frames") {
            spec.badge_background = "";
          } else if (spec.style === "cycle_segments" && spec.segment_backgrounds) {
            spec.segment_backgrounds = spec.segment_backgrounds.map(function() { return ""; });
          }
          resetBadgeChrome(badge);
          if (typeof wrap._l3dvaReapplyLetterColor === "function") {
            wrap._l3dvaReapplyLetterColor();
          }
          return;
        }
        if (spec.style === "static" || spec.style === "rotate_frames") {
          var vol = wrap._l3dvaVol;
          spec.badge_background = markerFillCssForVol(gd, trace, vol) || "";
          applyBadgeLetterColor(badge, spec.badge_background);
          return;
        }
        if (spec.style === "cycle_segments" && spec.segment_labels) {
          var ltv = wrap._l3dvaLabelToVol || {};
          var labels = spec.segment_labels;
          if (!spec.segment_backgrounds || spec.segment_backgrounds.length !== labels.length) {
            spec.segment_backgrounds = labels.map(function() { return ""; });
          }
          for (var si = 0; si < labels.length; si++) {
            var lbl = labels[si] == null ? "" : String(labels[si]);
            var vv = ltv[lbl];
            spec.segment_backgrounds[si] = vv != null
              ? (markerFillCssForVol(gd, trace, vv) || "")
              : "";
          }
          if (typeof wrap._l3dvaReapplyLetterColor === "function") {
            wrap._l3dvaReapplyLetterColor();
          }
        }
      });
    }

    function boldAnchorPointForVol(gd, trace, vol) {
      if (Object.prototype.hasOwnProperty.call(volDblClickPointNum, vol)) {
        var anchored = volDblClickPointNum[vol];
        if (
          anchored != null
          && trace.x
          && anchored >= 0
          && anchored < trace.x.length
          && volIdAtPointIndex(gd, trace, anchored) === vol
        ) {
          return anchored;
        }
        delete volDblClickPointNum[vol];
      }
      if (lastDblClickVol != null && vol === lastDblClickVol) {
        var centLast = sketchCentroidPointIndexForVol(gd, trace, vol);
        if (centLast != null) return centLast;
      }
      return null;
    }

    /**
     * Full-size plot letters: sketch-centroid particle, or bold anchor on the volume last
     * double-clicked in the current k-means selection.
     */
    function montagePlotLetterFontSizePx(isCent, isBold, vol) {
      var primary = !!isCent || (!!isBold && lastDblClickVol != null && vol === lastDblClickVol);
      return primary ? VOL_MONTAGE_PLOT_LETTER_PX : VOL_MONTAGE_PLOT_LETTER_PX_SECONDARY;
    }

    /** Live orbit for montage relayout / GIF tail — never the double-click pointerdown snapshot. */
    function captureVolAnimScenePinLiveOnly() {
      var pin = null;
      try {
        pin = P3S.snapshot(gd);
      } catch (ePin) { /* ignore */ }
      if (!sceneSnapHasCamera(pin) && typeof P3S.getPinnedSnap === "function") {
        pin = P3S.getPinnedSnap(gd);
      }
      if (!sceneSnapHasCamera(pin) && gd._cryoLatent3dLoadScenePin) {
        pin = gd._cryoLatent3dLoadScenePin;
      }
      return sceneSnapHasCamera(pin) ? pin : null;
    }

    /**
     * Pointerdown merge only during the brief Plotly double-click window (pre-reset orbit).
     * Using it later would snap back to the view at volume selection after the user orbited.
     */
    function captureVolAnimScenePin() {
      var pin = captureVolAnimScenePinLiveOnly();
      if (
        volAnimPointerdownSnap
        && volAnimPointerdownAt
        && Date.now() - volAnimPointerdownAt < 1200
      ) {
        pin = mergeSnapsPreferPointerdown(volAnimPointerdownSnap, pin);
      }
      return sceneSnapHasCamera(pin) ? pin : null;
    }

    function scheduleRestoreVolAnimScenePin(pin, onDone) {
      if (!pin || !sceneSnapHasCamera(pin)) {
        if (typeof onDone === "function") {
          try {
            onDone();
          } catch (eDone) { /* ignore */ }
        }
        return Promise.resolve();
      }
      if (typeof P3S.scheduleEnforceCamera === "function") {
        return new Promise(function(res) {
          P3S.scheduleEnforceCamera(gd, pin, function() {
            if (typeof onDone === "function") {
              try {
                onDone();
              } catch (eDone) { /* ignore */ }
            }
            res();
          });
        });
      }
      if (typeof P3S.scheduleRestoreCameraOnly === "function") {
        return new Promise(function(res) {
          P3S.scheduleRestoreCameraOnly(gd, pin, function() {
            if (typeof onDone === "function") {
              try {
                onDone();
              } catch (eDone) { /* ignore */ }
            }
            res();
          });
        });
      }
      if (typeof P3S.scheduleRestore === "function") {
        return new Promise(function(res) {
          P3S.scheduleRestore(gd, pin, function() {
            if (typeof onDone === "function") {
              try {
                onDone();
              } catch (eDone) { /* ignore */ }
            }
            res();
          });
        });
      }
      var p = typeof P3S.restoreCameraOnly === "function"
        ? P3S.restoreCameraOnly(gd, pin)
        : P3S.restore(gd, pin);
      var chain = p && typeof p.then === "function" ? p.catch(function() {}) : Promise.resolve();
      return chain.then(function() {
        if (typeof onDone === "function") {
          try {
            onDone();
          } catch (eDone) { /* ignore */ }
        }
      });
    }

    function relayoutMontageAnnotations(sceneAnnotations, scenePin) {
      if (typeof Plotly === "undefined" || !Plotly.relayout) {
        return Promise.resolve();
      }
      var anns = sceneAnnotations || [];
      var patch = { "scene.annotations": anns };
      // Patching scene.camera when only clearing/updating letter overlays resets the main
      // scatter3d orbit after a colour-covariate redraw (especially with no volumes selected).
      var camPin = captureVolAnimScenePinLiveOnly() || scenePin;
      if (anns.length && camPin && camPin.camera && typeof camPin.camera === "object") {
        try {
          patch["scene.camera"] = JSON.parse(JSON.stringify(camPin.camera));
        } catch (eCam) {
          patch["scene.camera"] = camPin.camera;
        }
      }
      return Plotly.relayout(gd, patch).catch(function() {});
    }

    function refreshSelectionHighlight() {
      centroidVolsInSelection.clear();
      if (!gd || !gd.data || !gd.data[0]) return Promise.resolve();
      if (gd._cryoLatent3dLoadScenePin && volAnimSelOverlayDepth >= 1) {
        gd._cryoLatent3dSceneHoldUntil = Date.now() + 6000;
      }
      var scenePin = captureVolAnimScenePinLiveOnly();
      if (scenePin && typeof P3S.setPinnedSnap === "function") {
        P3S.setPinnedSnap(gd, scenePin);
      }
      syncStableSelectionLabels();
      var trace = gd.data[0];
      pruneVolDblClickAnchors(gd, trace);
      var xs = trace.x;
      if (!xs || !xs.length) return Promise.resolve();
      var n = xs.length;
      var hasIds = trace.ids && trace.ids.length === n;
      var cd = customdataForHighlight(gd, trace);
      if (!hasIds && (!cd || cd.length !== n)) return Promise.resolve();

      var sceneAnnotations = [];
      for (var i = 0; i < n; i++) {
        var v = volIdAtPointIndex(gd, trace, i);
        var on = !isNaN(v) && selectedVols.has(v);
        if (!on) continue;
        var letter = volMontageLabel[v] || "";
        var isCent = isSketchCentroidAtPointIndex(gd, trace, i);
        if (isCent) centroidVolsInSelection.add(v);
        var boldPt = boldAnchorPointForVol(gd, trace, v);
        var isBold = boldPt != null && i === boldPt;
        var labelText = isCent ? circledMontageLabelText(letter) : letter;
        var ann = volMontageSceneAnnotation(
          trace.x[i],
          trace.y[i],
          trace.z[i],
          labelText,
          isBold,
          isCent && !isBold,
          {
            fontColorCss: montageLetterFontColor(trace, gd, i),
            fontSizePx: montagePlotLetterFontSizePx(isCent, isBold, v),
            gd: gd,
          }
        );
        if (ann) sceneAnnotations.push(ann);
      }

      return new Promise(function(resolve) {
        var countThisRestyle = volAnimSelOverlayDepth >= 1 ? 1 : 0;
        if (countThisRestyle) {
          volAnimSelHighlightRestylesInFlight++;
        }
        function finishSelectionRendering() {
          if (countThisRestyle) {
            volAnimSelHighlightRestylesInFlight--;
            if (volAnimSelHighlightRestylesInFlight <= 0) {
              volAnimSelHighlightRestylesInFlight = 0;
              endVolAnimSelectionOverlay();
            }
          }
          resolve();
        }
        function finishSelectionRenderingWithCameraPin() {
          if (!scenePin || !sceneSnapHasCamera(scenePin)) {
            finishSelectionRendering();
            return;
          }
          var liveAfterRelayout = null;
          try {
            if (typeof P3S.snapshot === "function") liveAfterRelayout = P3S.snapshot(gd);
          } catch (eLive) { /* ignore */ }
          if (
            liveAfterRelayout
            && sceneSnapHasCamera(liveAfterRelayout)
            && typeof P3S.cameraEyeDiffers === "function"
            && !P3S.cameraEyeDiffers(scenePin, liveAfterRelayout)
          ) {
            finishSelectionRendering();
            return;
          }
          scheduleRestoreVolAnimScenePin(scenePin, finishSelectionRendering);
        }
        relayoutMontageAnnotations(sceneAnnotations, scenePin)
          .then(finishSelectionRenderingWithCameraPin)
          .catch(finishSelectionRenderingWithCameraPin);
      });
    }

    function toggleVol(v, fromPointNumber) {
      var now = Date.now();
      var dup = lastToggleVol === v && (now - lastToggleAtMs) < 260;
      // If a duplicate toggle arrives while the volume is already selected, we allow
      // re-anchoring (bold move). If it's a duplicate while the volume is not yet selected,
      // ignore to avoid double-scheduling GIF work.
      if (dup && !selectedVols.has(v)) return;
      lastToggleVol = v;
      lastToggleAtMs = now;

      var membershipChanged = !selectedVols.has(v);
      if (!membershipChanged) {
        // Volume already selected: re-anchor bold marker to the newly double-clicked point.
        beginVolAnimSelectionOverlay();
        if (fromPointNumber != null && fromPointNumber >= 0) {
          volDblClickPointNum[v] = fromPointNumber;
        }
        lastDblClickVol = v;
      } else {
        if (selectedVols.size >= maxSelectable()) {
          setAnimateStatus("Selection capped at the number of sketch volumes.", false, true);
          return;
        }
        beginVolAnimSelectionOverlay();
        selectedVols.add(v);
        lastDblClickVol = v;
        if (fromPointNumber != null && fromPointNumber >= 0) {
          volDblClickPointNum[v] = fromPointNumber;
        }
      }
      setAnimateStatus("", false, false);
      updateSelLabel();
      var hp = (fromPointNumber != null && fromPointNumber >= 0)
        ? new Promise(function(resolve) {
          requestAnimationFrame(function() {
            resolve(refreshSelectionHighlight());
          });
        })
        : refreshSelectionHighlight();

      // Only re-render ChimeraX GIF previews when the selected volume set changed.
      if (membershipChanged) {
        scheduleAutoGif("selection", hp);
      }
    }

    /**
     * Optional ``plotHighlightPromise``: when passed (volume selection paths), ``runAutoGif``
     * runs only after the scatter highlight + camera pin chain settles (still after the 600ms debounce).
     */
    function scheduleAutoGif(reason, plotHighlightPromise) {
      debouncedGifReason = reason || "selection";
      if (debounceGifTimer) clearTimeout(debounceGifTimer);
      var hp = plotHighlightPromise && typeof plotHighlightPromise.then === "function"
        ? plotHighlightPromise.catch(function() {})
        : Promise.resolve();
      debounceGifTimer = setTimeout(function() {
        debounceGifTimer = null;
        var r = debouncedGifReason;
        debouncedGifReason = "selection";
        hp.then(function() {
          runAutoGif(r);
        });
      }, 600);
    }

    function runAutoGif(reason) {
      reason = reason || "selection";
      var vols = Array.from(selectedVols).sort(function(a, b) { return a - b; });
      if (!vols.length) {
        landscapeAnimInFlight = false;
        safeSetPreviewAnimLoading(false);
        clearPreviewGrid();
        lastAnimToken = null;
        lastBatchMode = null;
        lastChimeraxViewMatrix = "";
        chimeraxViewMatrixUnavailable = false;
        syncViewRotationSummary();
        setAnimateStatus("", false, false);
        syncSaveGifButton();
        return;
      }
      var myGen = ++gifReqGen;
      landscapeAnimInFlight = true;
      lastChimeraxViewMatrix = "";
      chimeraxViewMatrixUnavailable = false;
      syncSaveGifButton();
      syncViewRotationSummary();
      var cpus = META && META.chimerax_cpus != null ? Number(META.chimerax_cpus) : 1;
      var cpuPhrase = cpus === 1 ? "1 CPU" : cpus + " CPUs";
      var keepPreviews = previewGrid && previewGrid.children.length > 0 && (
        reason === "color" || reason === "gif_mode" || reason === "gif_frames"
        || reason === "view_rotation"
      );
      var statusMsg;
      if (keepPreviews && reason === "color") {
        statusMsg = (
          "Re-rendering volume previews with the new plot colors — current GIFs keep playing until "
          + "the update finishes (ChimeraX, "
          + cpuPhrase + ")…"
        );
      } else if (keepPreviews && (reason === "gif_mode" || reason === "gif_frames")) {
        statusMsg = (
          "Re-rendering volume previews for the new animation settings — current GIFs keep playing until "
          + "the update finishes (ChimeraX, "
          + cpuPhrase + ")…"
        );
      } else if (keepPreviews && reason === "view_rotation") {
        statusMsg = (
          "Re-rendering volume previews for the new view angle — current GIFs keep playing until "
          + "the update finishes (ChimeraX, "
          + cpuPhrase + ")…"
        );
      } else {
        statusMsg = "Rendering GIFs with ChimeraX (using " + cpuPhrase + ")…";
      }
      safeSetPreviewAnimLoading(true);
      setAnimateStatus(statusMsg, false, false);
      if (!keepPreviews) {
        teardownPreviews();
        previewGrid.innerHTML = "";
      }
      var animPayload = {
        vol_indices: vols,
        mode: currentGifMode(),
        gif_frames: parseInt(gifFrames.value, 10) || 20,
        chimerax_cpus: cpus,
        color_mode: getColorMode(),
        view_rotations: currentViewRotationPayload(),
      };
      if (getPalette) {
        var pal = getPalette();
        if (pal) animPayload.palette = pal;
      }
      if (
        currentGifMode() === "cycle"
        && lastAnimToken
        && lastBatchMode === "rotate_each"
      ) {
        animPayload.reuse_rotate_keyframes_token = lastAnimToken;
      }
      postJson("/api/landscape_volpca/generate_animations", animPayload).then(function(j) {
        if (myGen !== gifReqGen) return;
        landscapeAnimInFlight = false;
        safeSetPreviewAnimLoading(false);
        lastAnimToken = j.token;
        lastBatchMode = j.batch_mode != null ? j.batch_mode : null;
        lastChimeraxViewMatrix = j.view_matrix != null ? String(j.view_matrix).trim() : "";
        if (!lastChimeraxViewMatrix && j.items && j.items.length && j.items[0].view_matrix != null) {
          lastChimeraxViewMatrix = String(j.items[0].view_matrix).trim();
        }
        chimeraxViewMatrixUnavailable = !lastChimeraxViewMatrix;
        syncViewRotationSummary();
        var ds = j.duration_s != null ? Number(j.duration_s) : NaN;
        var doneMsg = !isNaN(ds)
          ? ("Finished animating in " + Math.round(ds) + " s.")
          : "Finished animating.";
        setAnimateStatus(doneMsg, false, false);
        if (!j.items || !j.items.length) {
          syncSaveGifButton();
          return;
        }
        teardownPreviews();
        previewGrid.innerHTML = "";
        // For cycle_segments we only get back letters; invert the backend montage-label
        // mapping using the same sorted volume list we sent to the backend.
        var labelToVolByBackend = {};
        for (var li = 0; li < vols.length; li++) {
          labelToVolByBackend[montageLabelAt(li)] = vols[li];
        }
        j.items.forEach(function(it) {
          var fig = document.createElement("figure");
          var wrap = document.createElement("div");
          wrap.className = "volsketch-anim-wrap";
          var badge = document.createElement("span");
          badge.className = "volsketch-gif-overlay";
          var covBadge = document.createElement("span");
          covBadge.className = "volsketch-gif-overlay volsketch-gif-overlay--cov";
          var covVarEl = document.createElement("span");
          covVarEl.className = "volsketch-gif-cov-var";
          covVarEl.hidden = true;
          var covValEl = document.createElement("span");
          covValEl.className = "volsketch-gif-cov-val";
          covBadge.appendChild(covVarEl);
          covBadge.appendChild(covValEl);
          var img = document.createElement("img");
          img.alt = "";
          img.src = "data:image/gif;base64," + it.gif_b64;
          wrap.appendChild(img);
          wrap.appendChild(badge);
          wrap.appendChild(covBadge);
          fig.appendChild(wrap);
          var overlaySpec = it.preview_overlay;
          if (overlaySpec && overlaySpec.style === "static") {
            var v = it.vol != null ? parseInt(String(it.vol), 10) : NaN;
            if (!isNaN(v) && centroidVolsInSelection.has(v)) {
              overlaySpec = Object.assign({}, overlaySpec, {
                text: circledMontageLabelText(overlaySpec.text || ""),
              });
            }
          } else if (
            overlaySpec
            && overlaySpec.style === "cycle_segments"
            && Array.isArray(overlaySpec.segment_labels)
          ) {
            var segLabels = overlaySpec.segment_labels;
            var segLabelsCircled = segLabels.map(function(lbl) {
              var s = lbl == null ? "" : String(lbl);
              var vv = labelToVolByBackend[s];
              return vv != null && centroidVolsInSelection.has(vv)
                ? circledMontageLabelText(s)
                : s;
            });
            overlaySpec = Object.assign({}, overlaySpec, {
              segment_labels: segLabelsCircled,
            });
          }
          mountPreviewOverlay(wrap, badge, covBadge, covVarEl, covValEl, img, overlaySpec);
          wrap._l3dvaVol = it.vol != null ? parseInt(String(it.vol), 10) : NaN;
          wrap._l3dvaOverlaySpec = overlaySpec;
          wrap._l3dvaLabelToVol = labelToVolByBackend;
          previewGrid.appendChild(fig);
        });
        refreshPreviewOverlayLetterColors();
        syncPreviewGridLayout();
        if (j.rendered_vol_indices && j.rendered_vol_indices.length) {
          var synced = new Set();
          j.rendered_vol_indices.forEach(function(x) {
            var v = parseInt(x, 10);
            if (!isNaN(v)) synced.add(v);
          });
          selectedVols = synced;
          // Preserve bold/italic anchors for volumes that remain selected.
          // (Server may subsample render outputs; we only drop anchors for removed volumes.)
          if (lastDblClickVol != null && !selectedVols.has(lastDblClickVol)) {
            lastDblClickVol = null;
          }
          for (var k in volDblClickPointNum) {
            if (!Object.prototype.hasOwnProperty.call(volDblClickPointNum, k)) continue;
            var vv = parseInt(k, 10);
            if (!selectedVols.has(vv)) delete volDblClickPointNum[k];
          }
          updateSelLabel();
          refreshSelectionHighlight();
        }
        syncSaveGifButton();
      }).catch(function(e) {
        if (myGen !== gifReqGen) return;
        landscapeAnimInFlight = false;
        safeSetPreviewAnimLoading(false);
        chimeraxViewMatrixUnavailable = false;
        syncViewRotationSummary();
        setAnimateStatus(String(e), false, true);
        syncSaveGifButton();
      });
    }

    function ensurePlotEventsBound() {
      if (plotEventsBound) return;
      if (!gd || typeof gd.on !== "function") return;
      plotEventsBound = true;
      if (
        layoutHasVolAnim(gd)
        && typeof Plotly !== "undefined"
        && Plotly.relayout
        && typeof P3S.disableSceneSelectionRelayoutPatch === "function"
      ) {
        Plotly.relayout(gd, P3S.disableSceneSelectionRelayoutPatch()).catch(function() {});
      }
      function recordVolAnimPointerdownScene() {
        try {
          if (!layoutHasVolAnim(gd)) return;
          volAnimPointerdownSnap = P3S.snapshot(gd);
          if (sceneSnapHasCamera(volAnimPointerdownSnap) && typeof P3S.setPinnedSnap === "function") {
            P3S.setPinnedSnap(gd, volAnimPointerdownSnap);
          }
          volAnimPointerdownAt = Date.now();
        } catch (ePd) { /* ignore */ }
      }
      gd.addEventListener("pointerdown", recordVolAnimPointerdownScene, true);
      gd.on("plotly_doubleclick", function(ev) {
        // Plotly's own scatter3d double-click can clobber the scene camera + annotation text.
        // We restore the camera from the pointerdown snapshot; for cases where our
        // click-timer double-click detection doesn't fire (e.g. some points),
        // we optionally attempt the toggle here too.
        var s = mergeSnapsPreferPointerdown(
          volAnimPointerdownSnap,
          typeof P3S.getPinnedSnap === "function" ? P3S.getPinnedSnap(gd) : null
        );
        var pt = ev && ev.points && ev.points.length ? ev.points[0] : null;
        var shouldTryToggle = !!pt;
        var v = NaN;
        var fromPointNumber = pt ? pt.pointNumber : null;
        if (shouldTryToggle) {
          try {
            v = volFromPoint(gd, pt);
          } catch (eCd) {
            v = NaN;
          }
          shouldTryToggle = !isNaN(v);
        }

        function maybeToggleAfterRestore() {
          if (!shouldTryToggle) return;
          ignorePlotlyClickUntilMs = Date.now() + 900;
          beginVolAnimSelectionOverlay();
          toggleVol(v, fromPointNumber);
        }

        // If we don't have a snapshot yet, we can still try toggling.
        if (!s) {
          if (shouldTryToggle) requestAnimationFrame(maybeToggleAfterRestore);
          return;
        }

        var restoreFn = typeof P3S.restoreCameraOnly === "function"
          ? P3S.restoreCameraOnly
          : P3S.restore;
        restoreFn(gd, s);
        requestAnimationFrame(function() {
          restoreFn(gd, s);
          setTimeout(function() {
            restoreFn(gd, s);
            maybeToggleAfterRestore();
            volAnimPointerdownSnap = null;
            volAnimPointerdownAt = 0;
          }, 0);
        });
      });
      gd.on("plotly_click", function(ev) {
        if (!layoutHasVolAnim(gd)) return;
        if (!ev || !ev.points || !ev.points.length) return;
        if (Date.now() < ignorePlotlyClickUntilMs) return;
        var v = volFromPoint(gd, ev.points[0]);
        if (isNaN(v)) return;
        if (clickTimer !== null && clickLastVol === v) {
          clearTimeout(clickTimer);
          clickTimer = null;
          clickLastVol = null;
          ignorePlotlyClickUntilMs = Date.now() + 900;
          beginVolAnimSelectionOverlay();
          requestAnimationFrame(function() {
            toggleVol(v, ev.points[0].pointNumber);
          });
          return;
        }
        clickLastVol = v;
        if (clickTimer) clearTimeout(clickTimer);
        clickTimer = setTimeout(function() {
          clickTimer = null;
          clickLastVol = null;
        }, 320);
      });
    }

    function loadMeta(cb) {
      fetch("/api/landscape_volpca/meta")
        .then(function(r) { return r.json(); })
        .then(function(j) {
          if (!j.ok) {
            setAnimateStatus(j.error || "Meta failed", true);
            if (randomSelBtn) {
              randomSelBtn.disabled = true;
              randomSelBtn.textContent = "Choose volumes at random";
            }
            return;
          }
          META = j;
          syncRandomSelButton();
          updateSelLabel();
          if (cb) cb();
        })
        .catch(function(e) {
          setAnimateStatus(String(e), true);
          if (randomSelBtn) {
            randomSelBtn.disabled = true;
            randomSelBtn.textContent = "Choose volumes at random";
          }
        });
    }

    if (gifModeSwitch) {
      gifModeSwitch.addEventListener("change", function() {
        syncGifModeDescription();
        syncGifFramesVisibility();
        scheduleAutoGif("gif_mode");
      });
    }
    if (gifFrames) {
      gifFrames.addEventListener("change", function() {
        scheduleAutoGif("gif_frames");
      });
    }
    viewRotateBtns.forEach(function(btn) {
      btn.addEventListener("click", function() {
        var axis = String(btn.getAttribute("data-l3dva-view-axis") || "").toLowerCase();
        if (axis === "x" || axis === "y" || axis === "z") {
          applyViewRotation(axis);
        }
      });
    });

    var clearBtn = $("clear-sel");
    if (clearBtn) {
      clearBtn.addEventListener("click", function() {
        beginVolAnimSelectionOverlay();
        selectedVols.clear();
        volDblClickPointNum = {};
        lastDblClickVol = null;
        updateSelLabel();
        setAnimateStatus("", false, false);
        scheduleAutoGif("selection", refreshSelectionHighlight());
      });
    }

    if (randomSelBtn) {
      randomSelBtn.addEventListener("click", function() {
        var ids = allPlotVolIds();
        if (!ids.length) {
          setAnimateStatus("Wait for the scatter plot to finish loading.", false, true);
          return;
        }
        if (debounceGifTimer) {
          clearTimeout(debounceGifTimer);
          debounceGifTimer = null;
        }
        gifReqGen++;
        landscapeAnimInFlight = false;
        safeSetPreviewAnimLoading(false);
        clearPreviewGrid();
        lastAnimToken = null;
        lastBatchMode = null;
        setAnimateStatus("", false, false);
        syncSaveGifButton();
        function applyRandomPickAndHighlight() {
          var k = META && META.chimerax_cpus != null ? Number(META.chimerax_cpus) : 1;
          k = Math.max(1, Math.min(k, ids.length, maxSelectable()));
          beginVolAnimSelectionOverlay();
          selectedVols.clear();
          volDblClickPointNum = {};
          lastDblClickVol = null;
          shufflePick(ids, k).forEach(function(v) { selectedVols.add(v); });
          updateSelLabel();
          scheduleAutoGif("selection", refreshSelectionHighlight());
        }
        applyRandomPickAndHighlight();
      });
    }

    if (saveGifBtn) {
      saveGifBtn.addEventListener("click", function() {
        if (!lastAnimToken) {
          setAnimateStatus("Select volumes and wait for previews to finish rendering.", false, true);
          return;
        }
        var custom = window.prompt(
          "Save folder (leave empty for default: this run’s landscape k-means directory):",
          ""
        );
        if (custom === null) return;
        setAnimateStatus("Saving…", true, false);
        var body = { token: lastAnimToken };
        if (custom && custom.trim()) body.out_dir = custom.trim();
        postJson("/api/landscape_volpca/save_animations", body)
          .then(function(j) {
            setAnimateStatus("Saved: " + (j.paths || []).join(", "), false, false);
          })
          .catch(function(e) { setAnimateStatus(String(e), false, true); });
      });
    }

    loadMeta(function() {
      updateSelLabel();
      syncViewRotationSummary();
      syncSaveGifButton();
    });

    function reapplySelectionHighlightIfNeeded() {
      refreshPreviewOverlayLetterColors();
      return refreshSelectionHighlight();
    }

    return {
      onColorOrPaletteChanged: function() {
        refreshPreviewOverlayLetterColors();
        var pPlot = refreshSelectionHighlight();
        var chain = pPlot && typeof pPlot.then === "function"
          ? pPlot.catch(function() {})
          : Promise.resolve();
        chain.then(function() {
          scheduleAutoGif("color");
        });
      },
      afterPlotRedraw: function() {
        ensurePlotEventsBound();
        return reapplySelectionHighlightIfNeeded();
      },
      reapplySelectionHighlightIfNeeded: reapplySelectionHighlightIfNeeded,
    };
  }

  global.CryoLatent3dVolLandscapeAnim = { boot: boot };
})(typeof window !== "undefined" ? window : globalThis);
