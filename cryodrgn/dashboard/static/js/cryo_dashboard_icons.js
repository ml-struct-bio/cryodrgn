/**
 * Shared dashboard UI icons (save floppy disk, etc.).
 */
(function (global) {
  "use strict";

  /* Outline 3.5" floppy — no outer frame; fills the icon slot via CSS width/height 100%. */
  var FLOPPY_SVG =
    "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 24 24\""
    + " aria-hidden=\"true\" focusable=\"false\" class=\"cryo-save-floppy-icon\""
    + " fill=\"none\" stroke=\"currentColor\" stroke-width=\"1.35\""
    + " stroke-linecap=\"round\" stroke-linejoin=\"round\">"
    + "<path d=\"M3.75 6.25h8.35V4.85h3.65l1.05 1.05v12.25H3.75V6.25z\"/>"
    + "<rect x=\"5.15\" y=\"7.2\" width=\"9.35\" height=\"2.75\" rx=\"0.3\"/>"
    + "<rect x=\"12.55\" y=\"7.6\" width=\"1.45\" height=\"1.75\" rx=\"0.15\"/>"
    + "<rect x=\"4.65\" y=\"11.65\" width=\"12.7\" height=\"6.55\" rx=\"0.3\"/>"
    + "<line x1=\"4.65\" y1=\"13.05\" x2=\"17.35\" y2=\"13.05\"/>"
    + "</svg>";

  var SAVE_BUTTON_SELECTORS = [
    "#save-indices-pkl",
    "#save-indices-custom",
    "#sel-fb-save-btn",
    "#save-pairplot-png",
    "#l3dva-save-gif",
    "#l3d-plot-gif-modal-save",
    "#volsketch-save-gif",
    "#fb-save-btn",
    "#btn-save-zpath-workdir",
    "#btn-save-zpath-as",
    "#btn-generate-volumes"
  ].join(",");

  function floppyDiskSvgInnerHTML() {
    return FLOPPY_SVG;
  }

  function decorateSaveButton(btn, opts) {
    if (!btn || btn.getAttribute("data-cryo-save-icon") === "1") {
      return btn;
    }
    opts = opts || {};
    btn.setAttribute("data-cryo-save-icon", "1");
    btn.classList.add("cryo-btn-with-save-icon");
    if (opts.iconOnly) {
      btn.classList.add("cryo-btn-with-save-icon--icon-only");
      btn.innerHTML =
        "<span class=\"cryo-btn-save-icon-slot\" aria-hidden=\"true\">"
        + floppyDiskSvgInnerHTML()
        + "</span>";
      return btn;
    }
    var slot = document.createElement("span");
    slot.className = "cryo-btn-save-icon-slot";
    slot.setAttribute("aria-hidden", "true");
    slot.innerHTML = floppyDiskSvgInnerHTML();
    btn.insertBefore(slot, btn.firstChild);
    return btn;
  }

  /** Update plain-text label on a decorated save button without dropping the icon. */
  function setSaveButtonLabel(btn, text) {
    if (!btn) return;
    decorateSaveButton(btn);
    var labelEl = btn.querySelector(".cryo-btn-save-label");
    if (!labelEl) {
      labelEl = document.createElement("span");
      labelEl.className = "cryo-btn-save-label";
      btn.appendChild(labelEl);
    }
    labelEl.textContent = text;
  }

  function decorateDashboardSaveButtons(root) {
    root = root || (typeof document !== "undefined" ? document : null);
    if (!root || !root.querySelectorAll) return;
    var nodes = root.querySelectorAll(SAVE_BUTTON_SELECTORS);
    for (var i = 0; i < nodes.length; i++) {
      decorateSaveButton(nodes[i]);
    }
  }

  global.CryoDashboardIcons = {
    floppyDiskSvgInnerHTML: floppyDiskSvgInnerHTML,
    decorateSaveButton: decorateSaveButton,
    setSaveButtonLabel: setSaveButtonLabel,
    decorateDashboardSaveButtons: decorateDashboardSaveButtons
  };

  if (typeof document !== "undefined") {
    var onReady = function () {
      decorateDashboardSaveButtons(document);
    };
    if (document.readyState === "loading") {
      document.addEventListener("DOMContentLoaded", onReady);
    } else {
      onReady();
    }
  }
})(typeof window !== "undefined" ? window : this);
