/**
 * Shared helpers for filling axis / colour &lt;select&gt; elements from covariate columns.
 */
(function (global) {
  "use strict";

  function displayLabel(displayMap, col) {
    if (displayMap && displayMap[col]) return displayMap[col];
    return col;
  }

  function fillSelect(sel, values, includeNone, displayMap) {
    if (!sel) return;
    sel.innerHTML = "";
    if (includeNone) {
      var none = document.createElement("option");
      none.value = "none";
      none.textContent = "None";
      sel.appendChild(none);
    }
    (values || []).forEach(function (c) {
      var o = document.createElement("option");
      o.value = c;
      o.textContent = displayLabel(displayMap, c);
      sel.appendChild(o);
    });
  }

  function fillAxisSelect(sel, values, val, displayMap) {
    if (!sel) return;
    sel.innerHTML = "";
    (values || []).forEach(function (c) {
      var o = document.createElement("option");
      o.value = c;
      o.textContent = displayLabel(displayMap, c);
      sel.appendChild(o);
    });
    if (values && values.indexOf(val) >= 0) {
      sel.value = val;
    } else if (values && values.length) {
      sel.value = values[0];
    }
  }

  function fillColorSelect(sel, values, displayMap) {
    fillSelect(sel, values, true, displayMap);
  }

  global.CryoCovariateSelects = {
    fillSelect: fillSelect,
    fillAxisSelect: fillAxisSelect,
    fillColorSelect: fillColorSelect,
  };
})(typeof window !== "undefined" ? window : this);
