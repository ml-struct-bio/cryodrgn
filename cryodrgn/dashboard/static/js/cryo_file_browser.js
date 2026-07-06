/**
 * Shared overlay modal + server file-browser listing for dashboard save/load dialogs.
 */
(function(global) {
  "use strict";

  var BODY_OPEN_CLASS = "cryo-explorer-save-modal-open";

  function parentDir(dir) {
    var parent = String(dir || "").replace(/[\\/]+$/, "");
    var idx = Math.max(parent.lastIndexOf("/"), parent.lastIndexOf("\\"));
    return idx > 0 ? parent.slice(0, idx) : parent;
  }

  function CryoDashModal() {}

  CryoDashModal.open = function(panel, opts) {
    opts = opts || {};
    if (!panel) return;
    panel.hidden = false;
    panel.setAttribute("aria-hidden", "false");
    document.body.classList.add(BODY_OPEN_CLASS);
    if (opts.focusEl) {
      requestAnimationFrame(function() {
        try {
          opts.focusEl.focus();
          if (typeof opts.focusEl.select === "function") {
            opts.focusEl.select();
          }
        } catch (eFocus) { /* ignore */ }
      });
    }
  };

  CryoDashModal.close = function(panel, opts) {
    opts = opts || {};
    if (!panel) return;
    panel.hidden = true;
    panel.setAttribute("aria-hidden", "true");
    document.body.classList.remove(BODY_OPEN_CLASS);
    var restore = opts.restoreFocus;
    if (restore && typeof restore.focus === "function") {
      try {
        restore.focus();
      } catch (eRestore) { /* ignore */ }
    }
  };

  /**
   * Wire backdrop / cancel / Escape close handlers on a modal panel.
   * opts: { closeAttr, cancelId, onClose, canClose, restoreFocusEl }
   */
  CryoDashModal.wire = function(panel, opts) {
    opts = opts || {};
    if (!panel) return;
    var closeAttr = opts.closeAttr || "data-modal-close";
    var onClose = typeof opts.onClose === "function" ? opts.onClose : function() {};
    var canClose = typeof opts.canClose === "function" ? opts.canClose : function() { return true; };

    function doClose() {
      if (!canClose()) return;
      CryoDashModal.close(panel, { restoreFocus: opts.restoreFocusEl });
      onClose();
    }

    panel.querySelectorAll("[" + closeAttr + "]").forEach(function(el) {
      el.addEventListener("click", doClose);
    });
    if (opts.cancelId) {
      var cancel = document.getElementById(opts.cancelId);
      if (cancel) {
        cancel.addEventListener("click", doClose);
      }
    }
    if (opts.escapeKey !== false) {
      document.addEventListener("keydown", function(ev) {
        if (ev.key !== "Escape" || panel.hidden) return;
        ev.preventDefault();
        doClose();
      });
    }
  };

  function CryoFileBrowser() {}

  /**
   * Populate a file-browser list from api_list_server_files.
   * opts: {
   *   listUrl, dir, kinds, listEl, pathEl,
   *   onDir, onFile, onLoaded, onError
   * }
   */
  CryoFileBrowser.loadDir = function(opts) {
    opts = opts || {};
    var listEl = opts.listEl;
    var pathEl = opts.pathEl;
    if (!listEl || !pathEl) return Promise.resolve(null);
    listEl.innerHTML = "<li class='cryo-file-browser-empty'>Loading\u2026</li>";
    var q = "";
    if (opts.kinds) {
      q += (q ? "&" : "?") + "kinds=" + encodeURIComponent(opts.kinds);
    }
    if (opts.dir) {
      q += (q ? "&" : "?") + "dir=" + encodeURIComponent(opts.dir);
    }
    var listUrl = opts.listUrl || "";
    return fetch(listUrl + q)
      .then(function(r) { return r.json(); })
      .then(function(j) {
        if (!j.ok) {
          listEl.innerHTML =
            "<li class='cryo-file-browser-empty'>" + (j.error || "Error") + "</li>";
          if (typeof opts.onError === "function") opts.onError(j.error || "Error");
          return null;
        }
        pathEl.textContent = j.dir;
        pathEl.title = j.dir;
        listEl.innerHTML = "";
        if (!j.entries || !j.entries.length) {
          listEl.innerHTML =
            "<li class='cryo-file-browser-empty'>No entries here</li>";
        } else {
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
                var next = j.dir + "/" + ent.name;
                if (typeof opts.onDir === "function") {
                  opts.onDir(next);
                } else {
                  CryoFileBrowser.loadDir(Object.assign({}, opts, { dir: next }));
                }
              });
            } else {
              li.addEventListener("click", function() {
                if (typeof opts.onFile === "function") {
                  opts.onFile(j.dir + "/" + ent.name, ent);
                }
              });
            }
            listEl.appendChild(li);
          });
        }
        if (typeof opts.onLoaded === "function") {
          opts.onLoaded(j);
        }
        return j;
      })
      .catch(function() {
        listEl.innerHTML =
          "<li class='cryo-file-browser-empty'>Could not list directory</li>";
        if (typeof opts.onError === "function") {
          opts.onError("Could not list directory");
        }
        return null;
      });
  };

  CryoFileBrowser.parentDir = parentDir;

  global.CryoDashModal = CryoDashModal;
  global.CryoFileBrowser = CryoFileBrowser;
})(typeof window !== "undefined" ? window : this);
