/** Auto-fit dashboard button label typography to the button box. */
  // Auto-fit typography in all dashboard buttons:
  // - maximise font-size while keeping label inside the current button box
  // - leave an internal inset margin (implemented as padding):
  //   - 5% of button height (top/bottom)
  //   - 10% of button width (left/right)
  (function() {
    var FIT_MARGIN_RATIO = 0.05;
    var FIT_MARGIN_RATIO_X = FIT_MARGIN_RATIO * 2;
    var FIT_GLOBAL_FONT_SHRINK_RATIO = 0.13;
    // Derived from the particle-explorer cache-toolbar buttons:
    //   globalMinBtnHeight = round(max(refHeights) * 0.77)
    var GLOBAL_MIN_BTN_HEIGHT_FRAC = 0.77;
    var GLOBAL_MIN_BTN_SESSION_KEY = "cryoDashGlobalMinBtnHeightPx";
    var globalMinBtnHeightPx = null;
    var globalMinLoadedFromSession = false;
    var FIT_MIN_FONT_PX = 4;
    var FIT_MAX_FONT_PX = 240;
    var FIT_OVERFLOW_TOL_PX = 0.75;

    function parseFinitePx(v) {
      var n = parseFloat(v || "");
      return isFinite(n) && n > 0 ? n : null;
    }

    function loadGlobalMinBtnHeightFromSession() {
      if (globalMinLoadedFromSession) return globalMinBtnHeightPx;
      globalMinLoadedFromSession = true;
      try {
        var v = sessionStorage.getItem(GLOBAL_MIN_BTN_SESSION_KEY);
        globalMinBtnHeightPx = parseFinitePx(v);
      } catch (e) {}
      return globalMinBtnHeightPx;
    }

    function computeGlobalMinBtnHeightFromParticleExplorerRefs() {
      var refExpand = document.getElementById("btn-expand-cache");
      var refView = document.getElementById("btn-view-images");
      var refH = 0;
      if (refExpand) {
        var r1 = refExpand.getBoundingClientRect();
        if (r1 && isFinite(r1.height) && r1.height > 2) refH = Math.max(refH, r1.height);
      }
      if (refView) {
        var r2 = refView.getBoundingClientRect();
        if (r2 && isFinite(r2.height) && r2.height > 2) refH = Math.max(refH, r2.height);
      }
      if (!(refH > 2)) return null;
      return Math.round(refH * GLOBAL_MIN_BTN_HEIGHT_FRAC);
    }

    function ensureGlobalMinBtnHeightPx() {
      if (globalMinBtnHeightPx != null && isFinite(globalMinBtnHeightPx) && globalMinBtnHeightPx > 0) {
        return globalMinBtnHeightPx;
      }
      // Prefer cached value so the exact same height applies across interfaces.
      loadGlobalMinBtnHeightFromSession();
      if (globalMinBtnHeightPx != null && isFinite(globalMinBtnHeightPx) && globalMinBtnHeightPx > 0) {
        return globalMinBtnHeightPx;
      }
      var computed = computeGlobalMinBtnHeightFromParticleExplorerRefs();
      if (computed != null && isFinite(computed) && computed > 0) {
        globalMinBtnHeightPx = computed;
        try { sessionStorage.setItem(GLOBAL_MIN_BTN_SESSION_KEY, String(computed)); } catch (e2) {}
        return globalMinBtnHeightPx;
      }
      return null;
    }

    function isPlotlyButton(btn) {
      return !!(
        btn
        && btn.closest
        && (btn.closest(".js-plotly-plot") || btn.closest(".modebar"))
      );
    }

    function buttonHasVisibleText(btn) {
      var t = (btn && btn.textContent) ? String(btn.textContent) : "";
      // Treat NBSP as whitespace; collapse all whitespace so only real label chars remain.
      t = t.replace(/\u00A0/g, " ").replace(/\s+/g, "");
      return t.length > 0;
    }

    function applyDescendantFontInherit(btn) {
      if (!btn || btn.dataset.cryoFitButtonsInheritApplied === "1") return;
      var els = btn.querySelectorAll("*");
      for (var i = 0; i < els.length; i++) {
        var el = els[i];
        var tag = (el.tagName || "").toLowerCase();
        if (tag === "svg" || tag === "path" || tag === "g") continue;
        el.style.fontSize = "inherit";
      }
      btn.dataset.cryoFitButtonsInheritApplied = "1";
    }

    function getCurrentFontSizePx(btn) {
      var inline = parseFloat(btn.style.fontSize || "");
      if (isFinite(inline)) return inline;
      var cs = window.getComputedStyle ? window.getComputedStyle(btn) : null;
      return cs ? parseFloat(cs.fontSize) : NaN;
    }

    // Toggle-list "select" headers (palette / discrete collapse titles).
    function isToggleListSelectButton(btn) {
      if (!btn || !btn.classList) return false;
      return (
        btn.classList.contains("cryo-palette-select__title")
        || btn.classList.contains("pairplot-palette-select__title")
        || btn.classList.contains("cryo-cc-discrete-toggle-title")
      );
    }

    // Estimate how many "code-defined" lines a button label contains:
    // - `<br>` tags always count as a line break.
    // - literal newline characters only count when `white-space` preserves them
    //   (e.g. `pre-line` / `pre-wrap`), so template indentation doesn't inflate `l`.
    function estimateButtonCodeLineCount(btn) {
      var brCount = 0;
      if (btn && btn.getElementsByTagName) {
        brCount = btn.getElementsByTagName("br").length || 0;
      }

      var ws = "";
      try {
        ws = window.getComputedStyle ? window.getComputedStyle(btn).whiteSpace : "";
      } catch (e) {}

      var countNewlines = ws === "pre-line" || ws === "pre-wrap" || ws === "pre";

      var nlCount = 0;
      if (countNewlines) {
        var txt = (btn && btn.textContent) ? String(btn.textContent) : "";
        txt = txt.replace(/\r/g, "\n");
        for (var i = 0; i < txt.length; i++) {
          if (txt.charAt(i) !== "\n") continue;
          var prev = txt.charAt(i - 1);
          var next = txt.charAt(i + 1);
          if (prev && next && prev.trim() !== "" && next.trim() !== "") nlCount++;
        }
      }

      var l = brCount + nlCount + 1;
      if (!isFinite(l) || l < 1) l = 1;
      return l;
    }

    // Global height base (`h`) rules requested by the UI update:
    // - one line:  minHeight = 0.6h
    // - l lines:   minHeight = h + 0.7h(l-1)
    // - toggle headers: minHeight = 0.7h
    function computeMinHeightForButton(btn, h) {
      if (h == null || !isFinite(h) || h <= 0) return null;
      if (isToggleListSelectButton(btn)) return h * 0.7;

      var l = estimateButtonCodeLineCount(btn);
      if (l <= 1) return h * 0.6;
      return h + 0.7 * h * (l - 1);
    }

    function isVolumeSliderTickButton(btn) {
      return !!(btn && btn.classList && btn.classList.contains("cryo-vslice-volume-slider-tick"));
    }

    function isDiscreteLegendSizeButton(btn) {
      return !!(
        btn
        && btn.classList
        && (btn.classList.contains("cryo-cc-discrete-size-btn")
          || btn.classList.contains("cryo-traj-scatter-color-legend__size-btn"))
      );
    }

    function isTrajVolPickerExpandButton(btn) {
      return !!(btn && btn.classList && btn.classList.contains("cryo-traj-vol-picker-expand"));
    }

    function isTrajVolGroupButton(btn) {
      return !!(btn && btn.classList && btn.classList.contains("cryo-traj-vol-group-btn"));
    }

    function isFileBrowserBarButton(btn) {
      // Compact path bars use fixed CSS sizing; global fit freezes a tiny box and
      // overwrites font/padding via inline styles (e.g. covariate .pkl “Up”).
      return !!(btn && btn.closest && btn.closest(".cryo-file-browser-bar"));
    }

    function fitButtonLabel(btn) {
      if (!btn || btn.nodeType !== 1) return;
      if (btn.classList && btn.classList.contains("btn-copy")) return;
      if (isVolumeSliderTickButton(btn)) return;
      if (isDiscreteLegendSizeButton(btn)) return;
      if (isTrajVolPickerExpandButton(btn)) return;
      if (isTrajVolGroupButton(btn)) return;
      if (isFileBrowserBarButton(btn)) return;
      if (btn.dataset.cryoFitButtonsBusy === "1") return;
      if (btn.hidden) return;
      if (isPlotlyButton(btn)) return;
      // ChimeraX path modal: flex + [hidden] gives bogus rects; do not shrink these buttons.
      if (btn.closest && btn.closest("#cryo-chimerax-modal")) return;
      // Header ChimeraX: custom sizing; global fit overwrites width/height/font.
      if (btn.classList && btn.classList.contains("nav-chimerax-btn")) return;
      // Command builder (GitHub Pages): advanced-parameter submenu headers use fixed CSS
      // title/subtitle sizes; global button-fit would overwrite them when the card resizes.
      if (btn.classList && btn.classList.contains("cmd-group-card-trigger")) return;
      if (!buttonHasVisibleText(btn)) return;

      var isCacheButtonTextMatchTarget =
        btn.id === "btn-cache-selection-uncached"
        || btn.id === "btn-clear-image-cache";
      var isVolumeHeightMatchTarget =
        btn.id === "btn-volume-generate" || btn.id === "btn-volume-animate";

      // Apply the requested global min-height rules based on:
      // - `h` = ensureGlobalMinBtnHeightPx()
      // - `l` = code-defined label lines (br + template newlines where `white-space` preserves them)
      var h = ensureGlobalMinBtnHeightPx();
      if (h != null && isFinite(h) && h > 0 && window.getComputedStyle) {
        var desiredMinH = computeMinHeightForButton(btn, h);
        if (desiredMinH != null && isFinite(desiredMinH) && desiredMinH > 0) {
          var csMin = parseFloat(window.getComputedStyle(btn).minHeight || "");
          if (!isFinite(csMin)) csMin = parseFloat(btn.style.minHeight || "");
          if (!isFinite(csMin)) csMin = 0;
          if (csMin + 0.5 < desiredMinH) {
            btn.style.minHeight = desiredMinH + "px";
          }
        }
      }

      // For cache/volume buttons, also freeze `height` to the computed min-height
      // so the font-size probe runs on stable geometry.
      if ((isCacheButtonTextMatchTarget || isVolumeHeightMatchTarget) && window.getComputedStyle) {
        if (h != null && isFinite(h) && h > 0) {
          var desiredH2 = computeMinHeightForButton(btn, h);
          if (desiredH2 != null && isFinite(desiredH2) && desiredH2 > 0) {
            var curInlineH = parseFloat(btn.style.height || "");
            if (!isFinite(curInlineH) || Math.abs(curInlineH - desiredH2) > 0.5) {
              btn.style.height = desiredH2 + "px";
              btn.style.minHeight = desiredH2 + "px";
            }
          }
        }
      }

      var rect = btn.getBoundingClientRect();
      var h = rect.height;
      var w = rect.width;
      if (!(h > 2 && w > 2)) return;

      var lastH = parseFloat(btn.dataset.cryoFitButtonsLastH || "");
      var lastW = parseFloat(btn.dataset.cryoFitButtonsLastW || "");
      var lastFont = parseFloat(btn.dataset.cryoFitButtonsLastFont || "");
      var curFont = getCurrentFontSizePx(btn);
      if (
        !isCacheButtonTextMatchTarget
        &&
        isFinite(lastH)
        && isFinite(lastW)
        && isFinite(lastFont)
        && isFinite(curFont)
        && Math.abs(lastH - h) < 0.5
        && Math.abs(lastW - w) < 0.5
        && Math.abs(lastFont - curFont) < 0.5
      ) {
        return;
      }

      btn.dataset.cryoFitButtonsBusy = "1";

      // Keep/restore these so we don't permanently override layout constraints.
      var origHeight = btn.style.height;
      var origWidth = btn.style.width;
      var origBoxSizing = btn.style.boxSizing;

      try {
        applyDescendantFontInherit(btn);

        // Internal margin (inset) as requested:
        // - top/bottom inset = 5% of button height
        // - left/right inset = 10% of button width (doubled)
        var marginYPx = Math.max(0, h * FIT_MARGIN_RATIO);
        var marginXPx = Math.max(0, w * FIT_MARGIN_RATIO_X);

        // Freeze the outer box while we probe for the max font-size.
        btn.style.boxSizing = "border-box";
        btn.style.height = h + "px";
        btn.style.width = w + "px";

        // Implement inset as padding so label layout naturally wraps/centres inside.
        btn.style.paddingTop = marginYPx + "px";
        btn.style.paddingBottom = marginYPx + "px";
        btn.style.paddingLeft = marginXPx + "px";
        btn.style.paddingRight = marginXPx + "px";

        var maxContentH = Math.floor(h - 2 * marginYPx);
        if (!(maxContentH > 0)) maxContentH = 1;

        var lo = FIT_MIN_FONT_PX;
        var hi = Math.min(FIT_MAX_FONT_PX, Math.max(lo, maxContentH));
        var best = lo;

        // Special-case: clear-selection button should visually start around
        // 80% of the font-size of the “Selected: <x>/<y> particles” label.
        // We treat this as an initial probe inside the binary search.
        var startProbeBestPx = null; // internal “best” font-size (before global shrink)
        if (btn.id === "clear-explorer-selection") {
          var selCountEl = document.getElementById("sel-count");
          if (selCountEl && window.getComputedStyle) {
            var selCs = window.getComputedStyle(selCountEl);
            var selFontPx = parseFloat(selCs.fontSize || "");
            if (isFinite(selFontPx) && selFontPx > 0) {
              var targetAppliedPx = selFontPx * 0.8;
              startProbeBestPx = targetAppliedPx / (1 - FIT_GLOBAL_FONT_SHRINK_RATIO);
              startProbeBestPx = Math.round(startProbeBestPx);
              startProbeBestPx = Math.max(lo, Math.min(hi, startProbeBestPx));
              if (!isFinite(startProbeBestPx)) startProbeBestPx = null;
            }
          }
        }

        if (startProbeBestPx != null && startProbeBestPx !== lo) {
          btn.style.fontSize = startProbeBestPx + "px";
          var startFits =
            btn.scrollHeight <= btn.clientHeight + FIT_OVERFLOW_TOL_PX
            && btn.scrollWidth <= btn.clientWidth + FIT_OVERFLOW_TOL_PX;
          if (startFits) {
            best = startProbeBestPx;
            lo = startProbeBestPx + 1;
          } else {
            hi = startProbeBestPx - 1;
          }
        }

        while (lo <= hi) {
          var mid = (lo + hi) >> 1;
          btn.style.fontSize = mid + "px";

          var fits =
            btn.scrollHeight <= btn.clientHeight + FIT_OVERFLOW_TOL_PX
            && btn.scrollWidth <= btn.clientWidth + FIT_OVERFLOW_TOL_PX;

          if (fits) {
            best = mid;
            lo = mid + 1;
          } else {
            hi = mid - 1;
          }
        }

        // User requested extra safety/legibility: reduce label size slightly
        // after the "best-fit" probe.
        var appliedFontPx = best * (1 - FIT_GLOBAL_FONT_SHRINK_RATIO);
        // Trajectory anchor-picker buttons need a bit more text legibility.
        if (
          btn.id === "btn-anchor-manual-open"
          || btn.id === "btn-anchor-import-pkl"
          || btn.id === "btn-anchor-random"
        ) {
          appliedFontPx = appliedFontPx * 1.1;
        }
        // Ensure we never exceed the binary-search "best-fit" maximum.
        if (appliedFontPx > best) appliedFontPx = best;
        appliedFontPx = Math.round(appliedFontPx * 100) / 100;
        btn.style.fontSize = appliedFontPx + "px";

        // For these two particle-explorer cache buttons, match the clear
        // selection button's computed text size (if that reference exists
        // and the exact size still fits in this button).
        if (isCacheButtonTextMatchTarget) {
          var clearBtn = document.getElementById("clear-explorer-selection");
          if (clearBtn && clearBtn.dataset) {
            var refFont = parseFloat(clearBtn.dataset.cryoFitButtonsLastFont || "");
            if (isFinite(refFont) && refFont > 0) {
              btn.style.fontSize = refFont + "px";
              var refFits =
                btn.scrollHeight <= btn.clientHeight + FIT_OVERFLOW_TOL_PX
                && btn.scrollWidth <= btn.clientWidth + FIT_OVERFLOW_TOL_PX;
              if (refFits) {
                appliedFontPx = refFont;
              } else {
                btn.style.fontSize = appliedFontPx + "px";
              }
            }
          }
        }

        btn.dataset.cryoFitButtonsLastH = h.toFixed(1);
        btn.dataset.cryoFitButtonsLastW = w.toFixed(1);
        btn.dataset.cryoFitButtonsLastFont = String(appliedFontPx);
      } finally {
        btn.dataset.cryoFitButtonsBusy = "0";
        btn.style.height = origHeight;
        btn.style.width = origWidth;
        btn.style.boxSizing = origBoxSizing;
      }
    }

    var pending = new Set();
    var scheduled = false;

    function scheduleFit(btn) {
      if (!btn || btn.nodeType !== 1) return;
      // Avoid queueing plotly modebar stuff.
      if (isPlotlyButton(btn)) return;
      pending.add(btn);
      if (scheduled) return;
      scheduled = true;
      var raf = window.requestAnimationFrame || function(cb) { return setTimeout(cb, 16); };
      raf(function() {
        scheduled = false;
        // Ensure reference button is fitted first, so other buttons can
        // copy its computed font-size deterministically.
        var clearBtnEl = document.getElementById("clear-explorer-selection");
        if (clearBtnEl && pending.has(clearBtnEl)) {
          fitButtonLabel(clearBtnEl);
          pending.delete(clearBtnEl);
        }
        pending.forEach(function(b) { fitButtonLabel(b); });
        pending.clear();
      });
    }

    var ro = null;
    if (window.CryoResizeObserverUtils && window.CryoResizeObserverUtils.create) {
      ro = window.CryoResizeObserverUtils.create(function(entries) {
        for (var i = 0; i < entries.length; i++) {
          var t = entries[i].target;
          if (!t || t.dataset && t.dataset.cryoFitButtonsBusy === "1") continue;
          scheduleFit(t);
        }
      });
    } else if (window.ResizeObserver) {
      ro = new ResizeObserver(function(entries) {
        for (var i = 0; i < entries.length; i++) {
          var t = entries[i].target;
          if (!t || t.dataset && t.dataset.cryoFitButtonsBusy === "1") continue;
          scheduleFit(t);
        }
      });
    }

    function observeButton(btn) {
      if (!btn || !ro) return;
      if (btn.dataset.cryoFitButtonsObserved === "1") return;
      if (btn.closest && btn.closest("#cryo-chimerax-modal")) {
        btn.dataset.cryoFitButtonsObserved = "1";
        return;
      }
      if (btn.classList && btn.classList.contains("nav-chimerax-btn")) {
        btn.dataset.cryoFitButtonsObserved = "1";
        return;
      }
      if (btn.classList && btn.classList.contains("cmd-group-card-trigger")) {
        btn.dataset.cryoFitButtonsObserved = "1";
        return;
      }
      if (btn.classList && btn.classList.contains("btn-copy")) {
        btn.dataset.cryoFitButtonsObserved = "1";
        return;
      }
      if (isVolumeSliderTickButton(btn)) {
        btn.dataset.cryoFitButtonsObserved = "1";
        return;
      }
      if (isTrajVolPickerExpandButton(btn)) {
        btn.dataset.cryoFitButtonsObserved = "1";
        return;
      }
      if (isTrajVolGroupButton(btn)) {
        btn.dataset.cryoFitButtonsObserved = "1";
        return;
      }
      if (isFileBrowserBarButton(btn)) {
        btn.dataset.cryoFitButtonsObserved = "1";
        return;
      }
      ro.observe(btn);
      btn.dataset.cryoFitButtonsObserved = "1";
    }

    function fitExistingButtons() {
      var buttons = document.querySelectorAll("button");
      for (var i = 0; i < buttons.length; i++) {
        var btn = buttons[i];
        observeButton(btn);
        if (btn.closest && btn.closest("#cryo-chimerax-modal")) continue;
        if (btn.classList && btn.classList.contains("nav-chimerax-btn")) continue;
        if (btn.classList && btn.classList.contains("cmd-group-card-trigger")) continue;
        if (btn.classList && btn.classList.contains("btn-copy")) continue;
        if (isVolumeSliderTickButton(btn)) continue;
        if (isTrajVolPickerExpandButton(btn)) continue;
        if (isTrajVolGroupButton(btn)) continue;
        if (isFileBrowserBarButton(btn)) continue;
        if (buttonHasVisibleText(btn)) scheduleFit(btn);
      }
    }

    // Handle buttons created after initial page load.
    if (document.body && window.MutationObserver) {
      var mo = new MutationObserver(function(mutations) {
        for (var mi = 0; mi < mutations.length; mi++) {
          var m = mutations[mi];
          if (!m.addedNodes || m.addedNodes.length === 0) continue;
          for (var ni = 0; ni < m.addedNodes.length; ni++) {
            var node = m.addedNodes[ni];
            if (!node || node.nodeType !== 1) continue;

            if (node.tagName && String(node.tagName).toLowerCase() === "button") {
              observeButton(node);
              if (!(node.closest && node.closest("#cryo-chimerax-modal"))
                  && !(node.classList && node.classList.contains("nav-chimerax-btn"))
                  && !(node.classList && node.classList.contains("cmd-group-card-trigger"))
                  && !(node.classList && node.classList.contains("btn-copy"))
                  && !isVolumeSliderTickButton(node)
                  && !isTrajVolPickerExpandButton(node)
                  && !isTrajVolGroupButton(node)
                  && !isFileBrowserBarButton(node)
                  && buttonHasVisibleText(node)) {
                scheduleFit(node);
              }
            } else if (node.querySelectorAll) {
              var innerButtons = node.querySelectorAll("button");
              for (var bi = 0; bi < innerButtons.length; bi++) {
                var ib = innerButtons[bi];
                observeButton(ib);
                if (!(ib.closest && ib.closest("#cryo-chimerax-modal"))
                    && !(ib.classList && ib.classList.contains("nav-chimerax-btn"))
                    && !(ib.classList && ib.classList.contains("cmd-group-card-trigger"))
                    && !(ib.classList && ib.classList.contains("btn-copy"))
                    && !isVolumeSliderTickButton(ib)
                    && !isTrajVolPickerExpandButton(ib)
                    && !isTrajVolGroupButton(ib)
                    && !isFileBrowserBarButton(ib)
                    && buttonHasVisibleText(ib)) {
                  scheduleFit(ib);
                }
              }
            }
          }
        }
      });
      mo.observe(document.body, { childList: true, subtree: true });
    }

    // Initial pass after the browser lays out the page.
    var initialRaf = window.requestAnimationFrame || function(cb) { return setTimeout(cb, 16); };
    initialRaf(function() { fitExistingButtons(); });

    // Also catch size changes from user resizing (some layout shifts don't always fire ResizeObserver fast enough).
    var resizeT = null;
    window.addEventListener("resize", function() {
      if (resizeT) clearTimeout(resizeT);
      resizeT = setTimeout(fitExistingButtons, 120);
    });
  })();
