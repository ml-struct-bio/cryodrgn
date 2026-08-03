/**
 * Shared UI for loading particle-indexed covariate numpy arrays from .pkl files.
 */
(function(global) {
  "use strict";

  var LOAD_OPTION_VALUE = "__load_covariate_pkl__";
  var LOAD_OPTION_LABEL = "Load covariate from .pkl\u2026";

  function byId(id) {
    return document.getElementById(id);
  }

  function appendLoadOptionToSelect(sel) {
    if (!sel) return;
    var existing = sel.querySelector('option[value="' + LOAD_OPTION_VALUE + '"]');
    if (existing) return;
    var o = document.createElement("option");
    o.value = LOAD_OPTION_VALUE;
    o.textContent = LOAD_OPTION_LABEL;
    sel.appendChild(o);
  }

  function CryoCovariatePklLoader(opts) {
    this.opts = opts || {};
    this.panel = byId("covariate-pkl-browser-panel");
    this.listEl = byId("covariate-pkl-list");
    this.pathEl = byId("covariate-pkl-path");
    this.statusEl = byId("covariate-pkl-status");
    this.currentDir = null;
    this.previousColorValue = "none";
    this._boundSelect = null;
    this._onSelectChange = null;
    this._wireModalControls();
  }

  CryoCovariatePklLoader.LOAD_OPTION_VALUE = LOAD_OPTION_VALUE;
  CryoCovariatePklLoader.LOAD_OPTION_LABEL = LOAD_OPTION_LABEL;

  CryoCovariatePklLoader.prototype._setStatus = function(msg, isError) {
    if (!this.statusEl) return;
    this.statusEl.textContent = msg || "";
    this.statusEl.style.color = isError ? "var(--accent,#c4703a)" : "";
  };

  CryoCovariatePklLoader.prototype._wireModalControls = function() {
    var self = this;
    var up = byId("covariate-pkl-up");
    if (up) {
      up.addEventListener("click", function() {
        if (!self.currentDir) return;
        self.loadDir(global.CryoFileBrowser
          ? CryoFileBrowser.parentDir(self.currentDir)
          : null);
      });
    }
    if (this.panel && global.CryoDashModal) {
      CryoDashModal.wire(this.panel, {
        closeAttr: "data-covariate-pkl-close",
        cancelId: "covariate-pkl-cancel",
        onClose: function() {
          self.currentDir = null;
          self._setStatus("", false);
        }
      });
    } else if (this.panel) {
      var cancel = byId("covariate-pkl-cancel");
      if (cancel) {
        cancel.addEventListener("click", function() { self.closeModal(); });
      }
      this.panel.querySelectorAll("[data-covariate-pkl-close]").forEach(function(el) {
        el.addEventListener("click", function() { self.closeModal(); });
      });
    }
  };

  CryoCovariatePklLoader.prototype.openModal = function() {
    if (!this.panel) return;
    this._setStatus("", false);
    if (global.CryoDashModal) {
      CryoDashModal.open(this.panel);
    } else {
      this.panel.hidden = false;
      this.panel.setAttribute("aria-hidden", "false");
      document.body.classList.add("cryo-explorer-save-modal-open");
    }
    this.loadDir(null);
  };

  CryoCovariatePklLoader.prototype.closeModal = function() {
    if (!this.panel) return;
    if (global.CryoDashModal) {
      CryoDashModal.close(this.panel);
    } else {
      this.panel.hidden = true;
      this.panel.setAttribute("aria-hidden", "true");
      document.body.classList.remove("cryo-explorer-save-modal-open");
    }
    this.currentDir = null;
    this._setStatus("", false);
  };

  CryoCovariatePklLoader.prototype.loadDir = function(dir) {
    var self = this;
    if (!this.listEl || !this.pathEl) return;
    if (global.CryoFileBrowser) {
      CryoFileBrowser.loadDir({
        listUrl: this.opts.listFilesUrl || "",
        dir: dir,
        kinds: "pkl",
        listEl: this.listEl,
        pathEl: this.pathEl,
        onDir: function(next) { self.loadDir(next); },
        onFile: function(path) { self.loadFile(path); },
        onLoaded: function(j) { self.currentDir = j.dir; }
      });
      return;
    }
    this.listEl.innerHTML = "<li class='cryo-file-browser-empty'>Loading\u2026</li>";
    var listUrl = this.opts.listFilesUrl || "";
    var q = "?kinds=pkl";
    if (dir) q += "&dir=" + encodeURIComponent(dir);
    fetch(listUrl + q)
      .then(function(r) { return r.json(); })
      .then(function(j) {
        if (!j.ok) {
          self.listEl.innerHTML =
            "<li class='cryo-file-browser-empty'>" + (j.error || "Error") + "</li>";
          return;
        }
        self.currentDir = j.dir;
        self.pathEl.textContent = j.dir;
        self.pathEl.title = j.dir;
        self.listEl.innerHTML = "";
        if (!j.entries || !j.entries.length) {
          self.listEl.innerHTML =
            "<li class='cryo-file-browser-empty'>No entries here</li>";
          return;
        }
        j.entries.forEach(function(ent) {
          var li = document.createElement("li");
          var icon = document.createElement("span");
          icon.className = "fb-icon";
          icon.textContent = ent.type === "dir" ? "\uD83D\uDCC1" : "\uD83D\uDCC4";
          var name = document.createElement("span");
          name.className = "fb-name";
          name.textContent = ent.name;
          li.appendChild(icon);
          li.appendChild(name);
          if (ent.type === "dir") {
            li.addEventListener("click", function() {
              self.loadDir(j.dir + "/" + ent.name);
            });
          } else {
            li.addEventListener("click", function() {
              self.loadFile(j.dir + "/" + ent.name);
            });
          }
          self.listEl.appendChild(li);
        });
      })
      .catch(function() {
        self.listEl.innerHTML =
          "<li class='cryo-file-browser-empty'>Could not list directory</li>";
      });
  };

  CryoCovariatePklLoader.prototype.loadFile = function(serverPath) {
    var self = this;
    var loadUrl = this.opts.loadCovariateUrl || "";
    this._setStatus("Loading covariate\u2026", false);
    fetch(loadUrl, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ path: serverPath })
    })
      .then(function(r) {
        return r.json().then(function(j) { return { ok: r.ok, j: j }; });
      })
      .then(function(res) {
        if (!res.ok || !res.j.ok) {
          var err = (res.j && res.j.error) || "Could not load covariate file.";
          self._setStatus(err, true);
          if (typeof self.opts.onError === "function") {
            self.opts.onError(err);
          }
          return;
        }
        self.closeModal();
        if (typeof self.opts.applyLoaded === "function") {
          self.opts.applyLoaded(res.j);
        }
        if (typeof self.opts.onLoaded === "function") {
          self.opts.onLoaded(res.j);
        }
      })
      .catch(function() {
        var err = "Covariate load request failed.";
        self._setStatus(err, true);
        if (typeof self.opts.onError === "function") {
          self.opts.onError(err);
        }
      });
  };

  CryoCovariatePklLoader.prototype.wireSelect = function(sel) {
    var self = this;
    if (!sel) return;
    appendLoadOptionToSelect(sel);
    this.previousColorValue = sel.value || "none";
    if (this._boundSelect === sel) return;
    if (this._boundSelect && this._onSelectChange) {
      this._boundSelect.removeEventListener("change", this._onSelectChange, true);
    }
    this._boundSelect = sel;
    this._onSelectChange = function(ev) {
      if (sel.value !== LOAD_OPTION_VALUE) {
        self.previousColorValue = sel.value;
        return;
      }
      var revert = self.previousColorValue || "none";
      sel.value = revert;
      if (ev && typeof ev.stopImmediatePropagation === "function") {
        ev.stopImmediatePropagation();
      }
      self.openModal();
    };
    sel.addEventListener("change", this._onSelectChange, true);
  };

  CryoCovariatePklLoader.prototype.wirePairGrid = function(container, getChecked, setChecked) {
    var self = this;
    if (!container) return;
    var btn = document.createElement("button");
    btn.type = "button";
    btn.className = "btn btn-secondary cryo-covariate-pkl-pair-btn";
    btn.textContent = LOAD_OPTION_LABEL;
    btn.style.marginTop = "0.45rem";
    btn.addEventListener("click", function() {
      self.previousColorValue = getChecked ? getChecked() : "none";
      self.openModal();
    });
    container.appendChild(btn);
    this._pairSetChecked = setChecked;
  };

  CryoCovariatePklLoader.prototype.appendPairGridRadios = function(
    container,
    columns,
    displayMap,
    primaryColumn
  ) {
    if (!container || !columns || !columns.length) return;
    var btn = container.querySelector(".cryo-covariate-pkl-pair-btn");
    var self = this;
    columns.forEach(function(col) {
      if (container.querySelector('input[name="color_cov"][value="' + col + '"]')) {
        return;
      }
      var label = document.createElement("label");
      label.className = "radio-line";
      var inp = document.createElement("input");
      inp.type = "radio";
      inp.name = "color_cov";
      inp.value = col;
      if (col === primaryColumn) inp.checked = true;
      var span = document.createElement("span");
      span.textContent = (displayMap && displayMap[col]) || col;
      label.appendChild(inp);
      label.appendChild(span);
      inp.addEventListener("change", function() {
        if (typeof self.opts.onPairColorChange === "function") {
          self.opts.onPairColorChange();
        }
      });
      if (btn) {
        container.insertBefore(label, btn);
      } else {
        container.appendChild(label);
      }
    });
  };

  CryoCovariatePklLoader.create = function(opts) {
    return new CryoCovariatePklLoader(opts);
  };

  global.CryoCovariatePklLoader = CryoCovariatePklLoader;
})(typeof window !== "undefined" ? window : this);
