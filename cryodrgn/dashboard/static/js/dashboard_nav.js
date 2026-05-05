(function() {
  var boot = window.CRYO_DASH_NAV_BOOT || {};
  var apiSetWorkdir = boot.apiSetWorkdir || "";
  var apiSetEpoch = boot.apiSetEpoch || "";
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

  var wdSel = document.getElementById("dash-workdir");
  if (wdSel) {
    wdSel.addEventListener("change", function() {
      var wd = String(this.value || "").trim();
      fetch(apiSetWorkdir, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ workdir: wd ? wd : null }),
      })
        .then(function(r) { return r.json().then(function(j) { return { ok: r.ok, j: j }; }); })
        .then(function(res) {
          if (res.ok && res.j.ok) window.location.reload();
        })
        .catch(function() {});
    });
  }

  var sel = document.getElementById("dash-epoch");
  if (!sel || sel.options.length <= 1) return;
  sel.addEventListener("change", function() {
    var v = parseInt(this.value, 10);
    if (isNaN(v)) return;
    fetch(apiSetEpoch, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ epoch: v }),
    })
      .then(function(r) { return r.json().then(function(j) { return { ok: r.ok, j: j }; }); })
      .then(function(res) {
        if (res.ok && res.j.ok) window.location.reload();
      })
      .catch(function() {});
  });
})();
