"""Access control: the dashboard answers only browsers that hold its launch token.

``cryodrgn dashboard`` listens on 127.0.0.1, which every user of a shared node can
reach, and several endpoints read or write files at paths given in the request
(covariate ``.pkl`` files are unpickled, z-paths, GIFs and volumes are saved where
asked). As in Jupyter, each launch therefore gets a random token: the command opens
and prints a link that carries it. The first request with the token answers with a
cookie and drops the token from the address bar; every later request must carry that
cookie. Scripts can send the token as the ``X-Cryodrgn-Token`` header or
as ``?token=`` on API calls instead.

``CRYODRGN_DASHBOARD_TOKEN`` fixes the token (for bookmarks); there is no way to turn
the check off from the command line.
Apps built by :func:`cryodrgn.dashboard.app.create_app` without a token (as the
tests do) answer everyone.
"""

from __future__ import annotations

import hashlib
import hmac
import os
import secrets
from urllib.parse import urlencode

from flask import Flask, Response, current_app, jsonify, redirect, request

TOKEN_ENV = "CRYODRGN_DASHBOARD_TOKEN"
TOKEN_HEADER = "X-Cryodrgn-Token"

# Package files only (scripts, styles, Plotly): nothing in them depends on the run.
_OPEN_ENDPOINTS = frozenset({"static", "bundled_plotly_js"})

_NEED_LINK_PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>cryoDRGN dashboard</title>
<style>body{font-family:Barlow,Roboto,sans-serif;max-width:40rem;margin:4rem auto;
padding:0 1rem;line-height:1.5;color:#1a2332}code{background:#f1ede6;padding:0 .2em}</style>
</head><body><h1>This dashboard needs its access link</h1>
<p>Open the link printed in the terminal where <code>cryodrgn dashboard</code> was
started. It ends in <code>?token=&hellip;</code> and changes with every launch.</p>
</body></html>
"""


def new_access_token() -> str:
    """The token for one launch: ``$CRYODRGN_DASHBOARD_TOKEN`` if set, else random."""
    fixed = os.environ.get(TOKEN_ENV, "").strip()

    return fixed or secrets.token_urlsafe(24)


def access_url(host: str, port: int, path: str, token: str) -> str:
    """The link that opens ``path`` with the token, for the browser and the terminal."""
    shown_host = "127.0.0.1" if host in ("0.0.0.0", "") else host
    if ":" in shown_host and not shown_host.startswith("["):
        shown_host = f"[{shown_host}]"
    sep = "&" if "?" in path else "?"

    return f"http://{shown_host}:{port}{path}{sep}{urlencode({'token': token})}"


def install_access_control(app: Flask, token: str) -> None:
    """Require ``token`` on every request to ``app`` but those for package files.

    Call before registering other ``before_request`` hooks, so that a refused request
    never loads the run. The cookie names derive from the token, so dashboards on two
    ports of one host keep their own cookies (browsers do not separate cookies by port).
    """
    tag = hashlib.sha256(token.encode()).hexdigest()[:12]
    app.config["DASHBOARD_ACCESS_TOKEN"] = token
    app.config["DASHBOARD_ACCESS_COOKIE"] = f"cryodrgn_dashboard_access_{tag}"
    app.config["SESSION_COOKIE_NAME"] = f"cryodrgn_dashboard_session_{tag}"
    app.config["SESSION_COOKIE_SAMESITE"] = "Lax"
    app.before_request(require_access_token)


def _same(given: str | None, token: str) -> bool:
    return bool(given) and hmac.compare_digest(given.encode(), token.encode())


def require_access_token() -> Response | None:
    """``before_request`` hook: let the request through, swap a URL token, or refuse."""
    token = current_app.config.get("DASHBOARD_ACCESS_TOKEN")

    if not token or request.endpoint in _OPEN_ENDPOINTS:
        return None

    cookie = current_app.config["DASHBOARD_ACCESS_COOKIE"]
    if _same(request.cookies.get(cookie), token) or _same(
        request.headers.get(TOKEN_HEADER), token
    ):
        return None

    if _same(request.args.get("token"), token):
        if request.method != "GET" or request.path.startswith("/api/"):
            return None

        # A page opened from the printed link: keep the token in a cookie instead of
        # the address bar, history and bookmarks.
        args = request.args.to_dict(flat=False)
        args.pop("token", None)
        target = request.path + (f"?{urlencode(args, doseq=True)}" if args else "")
        resp = redirect(target)
        resp.set_cookie(cookie, token, httponly=True, samesite="Strict")

        return resp

    if request.path.startswith("/api/") or (
        request.accept_mimetypes.best == "application/json"
    ):
        msg = (
            "This dashboard needs its access link: open the link printed in the "
            "terminal where `cryodrgn dashboard` was started."
        )

        return jsonify(error=msg), 403

    return Response(_NEED_LINK_PAGE, status=403, mimetype="text/html")
