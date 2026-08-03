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

    function cryoPostWorkdir(wd) {
      fetch(cfg.setWorkdirUrl, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ workdir: wd ? wd : null })
      })
        .then(function(r) { return r.json().then(function(j) { return { ok: r.ok, j: j }; }); })
        .then(function(res) {
          if (res.ok && res.j.ok) window.location.reload();
        })
        .catch(function() {});
    }

    var wdSel = document.getElementById("dash-workdir");
    if (wdSel) {
      wdSel.addEventListener("change", function() {
        cryoPostWorkdir(String(this.value || "").trim());
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
      sel.addEventListener("change", function() {
        var v = parseInt(this.value, 10);
        if (isNaN(v)) return;
        fetch(cfg.setEpochUrl, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ epoch: v })
        })
          .then(function(r) {
            return r.json().then(function(j) { return { ok: r.ok, j: j }; });
          })
          .then(function(res) {
            if (res.ok && res.j.ok) window.location.reload();
          })
          .catch(function() {});
      });
    }
  })();
