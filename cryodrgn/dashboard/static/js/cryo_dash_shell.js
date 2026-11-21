(function() {
  var cfg = window.CRYO_DASH_SHELL_CFG || {};

    if (window.cryoRefreshStatusTones) {
      window.cryoRefreshStatusTones(document);
      var statusToneObserver = new MutationObserver(function() {
        window.cryoRefreshStatusTones(document);
      });
      statusToneObserver.observe(document.body, {
        childList: true,
        subtree: true,
        characterData: true,
      });
    }

    /**
     * POST a switch of run folder or epoch and reload on success. On failure the page
     * stays on its data, so the menu goes back to its value and the reason is shown;
     * it used to keep the new value over the old data, silently.
     */
    function cryoPostSwitch(url, body, what, menu) {
      var previous = menu ? menu.getAttribute("data-cryo-current") : null;
      function fail(reason) {
        if (menu && previous !== null) menu.value = previous;
        console.error("[cryodrgn] switching " + what + " failed", reason);
        window.alert("Could not switch " + what + ": " + reason);
      }
      fetch(url, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(body)
      })
        .then(function(r) { return r.json().then(function(j) { return { ok: r.ok, j: j }; }); })
        .then(function(res) {
          if (res.ok && res.j.ok) {
            window.location.reload();
            return;
          }
          fail((res.j && res.j.error) || "the server refused it.");
        })
        .catch(function(err) {
          fail(String((err && err.message) || err || "the request failed."));
        });
    }

    function cryoPostWorkdir(wd, menu) {
      cryoPostSwitch(cfg.setWorkdirUrl, { workdir: wd ? wd : null }, "the run folder", menu);
    }

    var wdSel = document.getElementById("dash-workdir");
    if (wdSel) {
      wdSel.setAttribute("data-cryo-current", wdSel.value);
      wdSel.addEventListener("change", function() {
        cryoPostWorkdir(String(this.value || "").trim(), wdSel);
      });
    }

    var wdLaunchBtns = document.querySelectorAll(".cryo-dash-set-workdir[data-workdir]");
    for (var bi = 0; bi < wdLaunchBtns.length; bi++) {
      wdLaunchBtns[bi].addEventListener("click", function() {
        var w = String(this.getAttribute("data-workdir") || "").trim();
        if (!w) return;
        cryoPostWorkdir(w);
      });
    }

    (function() {
      var openBtn = document.getElementById("cryo-chimerax-setup-open");
      var modal = document.getElementById("cryo-chimerax-modal");
      if (!openBtn || !modal) return;
      var apiUrl = String(openBtn.getAttribute("data-api-url") || "").trim();
      var inp = document.getElementById("cryo-chimerax-path-input");
      var errEl = document.getElementById("cryo-chimerax-modal-err");
      var saveBtn = document.getElementById("cryo-chimerax-save");
      var clearBtn = document.getElementById("cryo-chimerax-clear");

      function setErr(s) {
        if (!errEl) return;
        errEl.textContent = s || "";
        if (window.cryoApplyStatusTone) window.cryoApplyStatusTone(errEl, s || "");
      }

      function openModal() {
        setErr("");
        modal.removeAttribute("hidden");
        modal.setAttribute("aria-hidden", "false");
        if (inp) {
          try {
            inp.focus();
            inp.select();
          } catch (e2) {}
        }
      }

      function closeModal() {
        modal.setAttribute("hidden", "hidden");
        modal.setAttribute("aria-hidden", "true");
        setErr("");
      }

      openBtn.addEventListener("click", openModal);

      modal.addEventListener("click", function(ev) {
        var t = ev.target;
        if (t && t.getAttribute && t.getAttribute("data-cryo-chimerax-close") === "1") {
          closeModal();
        }
      });

      document.addEventListener("keydown", function(ev) {
        if (ev.key === "Escape" && !modal.hasAttribute("hidden")) closeModal();
      });

      if (saveBtn) {
        saveBtn.addEventListener("click", function() {
          var p = inp ? String(inp.value || "").trim() : "";
          if (!p) {
            setErr("Enter the full path to the ChimeraX executable.");
            return;
          }
          setErr("");
          fetch(apiUrl, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ path: p })
          })
            .then(function(r) {
              return r.json().then(function(j) {
                return { ok: r.ok, j: j };
              });
            })
            .then(function(res) {
              if (res.ok && res.j && res.j.ok) {
                window.location.reload();
                return;
              }
              var msg = res.j && res.j.error ? String(res.j.error) : "Could not save path.";
              setErr(msg);
            })
            .catch(function() {
              setErr("Request failed.");
            });
        });
      }

      if (clearBtn) {
        clearBtn.addEventListener("click", function() {
          setErr("");
          fetch(apiUrl, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ path: null })
          })
            .then(function(r) {
              return r.json().then(function(j) {
                return { ok: r.ok, j: j };
              });
            })
            .then(function(res) {
              if (res.ok && res.j && res.j.ok) {
                window.location.reload();
                return;
              }
              var msg = res.j && res.j.error ? String(res.j.error) : "Could not clear.";
              setErr(msg);
            })
            .catch(function() {
              setErr("Request failed.");
            });
        });
      }
    })();

    var sel = document.getElementById("dash-epoch");
    if (sel && sel.options.length > 1) {
      sel.setAttribute("data-cryo-current", sel.value);
      sel.addEventListener("change", function() {
        var v = parseInt(this.value, 10);
        if (isNaN(v)) return;
        cryoPostSwitch(cfg.setEpochUrl, { epoch: v }, "the epoch", sel);
      });
    }
  })();
