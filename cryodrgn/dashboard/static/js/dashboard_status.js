    (function() {
      var ERROR_RE = /\b(error|failed|failure|invalid|unable|cannot|could not|not found|request failed|bad)\b/i;
      function asText(msg) {
        return String(msg == null ? "" : msg).trim();
      }
      window.cryoIsErrorMessage = function(msg) {
        var s = asText(msg);
        if (!s) return false;
        return ERROR_RE.test(s);
      };
      window.cryoApplyStatusTone = function(el, msg) {
        if (!el || !el.classList) return;
        var s = asText(msg);
        el.classList.toggle("cryo-status-error", window.cryoIsErrorMessage(s));
      };
      window.cryoSetStatus = function(el, msg) {
        if (!el) return;
        var s = asText(msg);
        el.textContent = s;
        window.cryoApplyStatusTone(el, s);
      };
      window.cryoRefreshStatusTones = function(scope) {
        var root = scope || document;
        var nodes = root.querySelectorAll("[id*='status'], .status, .cryo-dash-legend-note, .cmd-copy-status");
        for (var i = 0; i < nodes.length; i++) {
          var el = nodes[i];
          if (!el || !el.classList) continue;
          window.cryoApplyStatusTone(el, el.textContent || "");
        }
      };
    })();
