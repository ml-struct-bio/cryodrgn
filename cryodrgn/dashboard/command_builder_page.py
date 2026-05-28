"""Build a static GitHub Pages bundle for the dashboard command builder.

Renders ``command_builder.html`` with a minimal Flask app (no torch / experiment
imports) so CI can publish ``https://<org>.github.io/<repo>/`` on each release tag.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
from pathlib import Path
from urllib.parse import quote

from flask import Flask, render_template

from cryodrgn.dashboard.command_builder_cli_help import (
    jinja_arg_display_name,
    load_command_module_docstrings,
)
from cryodrgn.dashboard.command_builder_data import (
    COMMAND_BUILDER_MANUSCRIPT_LABELS,
    COMMAND_BUILDER_MANUSCRIPT_URLS,
    COMMAND_BUILDER_REQUIRED_FIELD_TITLES,
    COMMAND_BUILDER_SCHEMA,
    arg_is_batch_size_denominated,
    arg_is_epoch_denominated,
    arg_is_num_epochs,
    arg_show_display_name,
    default_outdir_for_command,
)

_THIS_DIR = Path(__file__).resolve().parent
_TEMPLATE_DIR = _THIS_DIR / "templates"
_STATIC_DIR = _THIS_DIR / "static"

# Root-absolute paths in rendered HTML (Flask ``url_for`` / hard-coded ``/``).
_ROOT_PATH_RE = re.compile(r'(?P<attr>(?:href|src|action))="(?P<path>/[^"]*)"')
_PLOTLY_SCRIPT_RE = re.compile(
    r'\s*<script src="https://cdn\.plot\.ly/plotly[^"]*"[^>]*></script>\s*',
    re.IGNORECASE,
)


def _abbrev_middle(text: object, maxlen: int = 30) -> str:
    s = str(text or "").strip()
    if len(s) <= maxlen:
        return s
    keep = maxlen - 1
    left = keep // 2
    right = keep - left
    return s[:left] + "…" + s[-right:]


def _nav_interface_title(text: object) -> str:
    return _abbrev_middle(text, maxlen=36)


def _cryodrgn_version_context() -> dict[str, str]:
    ver = os.environ.get("CRYODRGN_VERSION", "").strip()
    if not ver:
        try:
            from cryodrgn import __version__ as ver  # type: ignore[attr-defined]
        except Exception:
            ver = "unknown"
    ver = str(ver).strip() or "unknown"
    short = ver.split("+", 1)[0].strip() if "+" in ver else ver
    return {"cryodrgn_version": ver, "cryodrgn_version_short": short}


# PEP 440 ``4.3.0a8`` → git tag ``4.3.0-a8`` (hyphen before pre-release letter).
_PEP440_PRERELEASE_TAG_RE = re.compile(
    r"^(\d+(?:\.\d+)*)(?<!-)([ab]\d+|rc\d+)(.*)$",
    re.IGNORECASE,
)


def _version_to_github_tree_tag(version: str) -> str:
    """Normalize a version string to the GitHub tag used for ``/tree/<tag>`` links."""
    tag = (version or "").strip().split("+", 1)[0].strip()
    if not tag:
        return tag
    m = _PEP440_PRERELEASE_TAG_RE.match(tag)
    if m:
        return f"{m.group(1)}-{m.group(2)}{m.group(3)}"
    return tag


def _github_repo_release_url(repo_url: str, version: str) -> str:
    """Link to the GitHub tree at the built cryoDRGN version tag (e.g. ``/tree/4.3.0-a8``)."""
    base = repo_url.rstrip("/")
    ver = (version or "").strip()
    if not ver or ver == "unknown":
        return base
    tag = _version_to_github_tree_tag(ver) or ver
    return f"{base}/tree/{quote(tag)}"


def _command_builder_template_kwargs() -> dict[str, object]:
    return {
        "default_particles": "",
        "default_ctf": "",
        "default_workdir": "",
        "default_epoch": "",
        "default_zdim": 8,
        "default_outdir_abinit": default_outdir_for_command("abinit"),
        "default_outdir_abinit_het_old": default_outdir_for_command("abinit_het_old"),
        "default_outdir_abinit_homo_old": default_outdir_for_command("abinit_homo_old"),
        "default_outdir_train_vae": default_outdir_for_command("train_vae"),
        "default_outdir_train_nn": default_outdir_for_command("train_nn"),
        "default_outdir_train_dec": default_outdir_for_command("train_dec"),
        "default_outdir_backproject_voxel": default_outdir_for_command(
            "backproject_voxel"
        ),
        "default_poses": "",
        "command_builder_schema": COMMAND_BUILDER_SCHEMA,
        "command_builder_required_field_titles": COMMAND_BUILDER_REQUIRED_FIELD_TITLES,
        "command_builder_command_docs": load_command_module_docstrings(),
        "command_builder_manuscript_urls": COMMAND_BUILDER_MANUSCRIPT_URLS,
        "command_builder_manuscript_labels": COMMAND_BUILDER_MANUSCRIPT_LABELS,
    }


def _command_builder_page_base_path(repo: str | None = None) -> str:
    """Project-site base path ``/<repo>/`` (trailing slash)."""
    if repo is None:
        repo = os.environ.get("GITHUB_REPOSITORY", "ml-struct-bio/cryodrgn")
    name = repo.split("/", 1)[-1].strip() or "cryodrgn"
    return f"/{name}/"


def _rewrite_root_paths(html: str, base_path: str) -> str:
    """Turn ``/…`` asset URLs into paths relative to ``<base href>``."""
    base = base_path if base_path.endswith("/") else f"{base_path}/"
    base_no_slash = base.rstrip("/")

    def repl(match: re.Match[str]) -> str:
        attr = match.group("attr")
        path = match.group("path")
        if path.startswith("//"):
            return match.group(0)
        if path == "/":
            return f'{attr}="{base}"'
        if path.startswith("/static/"):
            return f'{attr}="{path.removeprefix("/")}"'
        if path.startswith("/api/"):
            return f'{attr}="#"'
        if path in ("/command-builder", "/abinit-builder"):
            return f'{attr}="{base}"'
        if path.startswith(base_no_slash):
            rel = path[len(base_no_slash) :].lstrip("/")
            return f'{attr}="{rel or "."}"'
        return match.group(0)

    return _ROOT_PATH_RE.sub(repl, html)


def _adapt_html_for_command_builder_page(
    html: str, base_path: str, repo_url: str
) -> str:
    html = _PLOTLY_SCRIPT_RE.sub("\n", html)
    html = _rewrite_root_paths(html, base_path)
    # Brand link targets https://cryodrgn.cs.princeton.edu/ (set in the template).
    base = base_path if base_path.endswith("/") else f"{base_path}/"
    if "<base " not in html:
        html = html.replace("<head>", f'<head>\n  <base href="{base}"/>', 1)
    # ``base.html`` workdir/epoch handlers are dashboard-only;
    # disable on static hosting.
    html = html.replace('fetch("/api/set_workdir"', 'fetch("#"')
    html = html.replace('fetch("/api/set_epoch"', 'fetch("#"')
    return html


def render_command_builder_html(
    *,
    base_path: str | None = None,
    repo_url: str | None = None,
) -> str:
    """Return a single self-contained command-builder HTML document."""
    if base_path is None:
        base_path = _command_builder_page_base_path()
    if repo_url is None:
        repo = os.environ.get("GITHUB_REPOSITORY", "ml-struct-bio/cryodrgn")
        repo_url = f"https://github.com/{repo}"

    app = Flask(
        __name__,
        template_folder=str(_TEMPLATE_DIR),
        static_folder=str(_STATIC_DIR),
        static_url_path="/static",
    )
    app.jinja_env.filters["abbrev_middle"] = _abbrev_middle
    app.jinja_env.filters["nav_interface_title"] = _nav_interface_title
    app.jinja_env.filters["arg_is_epoch_denominated"] = arg_is_epoch_denominated
    app.jinja_env.filters[
        "arg_is_batch_size_denominated"
    ] = arg_is_batch_size_denominated
    app.jinja_env.filters["arg_is_num_epochs"] = arg_is_num_epochs
    app.jinja_env.filters["arg_display_name"] = jinja_arg_display_name
    app.jinja_env.filters["arg_show_display_name"] = arg_show_display_name
    app.config["COMMAND_BUILDER_ONLY"] = True

    def _noop_api() -> tuple[str, int]:
        return "", 404

    for rule, endpoint in (
        ("/api/set_workdir", "api_set_workdir"),
        ("/api/set_epoch", "api_set_epoch"),
    ):
        app.add_url_rule(rule, endpoint=endpoint, view_func=_noop_api, methods=["POST"])

    ver_ctx = _cryodrgn_version_context()
    ctx = {
        **_command_builder_template_kwargs(),
        **ver_ctx,
        "github_repo_url": repo_url,
        "github_repo_release_url": _github_repo_release_url(
            repo_url, ver_ctx["cryodrgn_version"]
        ),
        "exp_workdir": "",
        "exp_epoch": 0,
        "exp_kmeans": -1,
        "filter_plot_inds_default": "",
        "dashboard_epochs": [0],
        "cfg_model_type": "cryodrgn",
        "cfg_zdim": "",
        "cfg_train_command": "",
        "cfg_cmd_display_lines": ["No experiment loaded"],
        "command_builder_only": True,
        "discovered_workdirs": [],
        "discovered_workdir_options": [],
        "selected_workdir": "",
        "exp_epoch_output_stamp": None,
        "exp_epoch_output_stamp_title": None,
        "exp_run_log_cryodrgn_version": None,
        "exp_run_log_cryodrgn_version_short": None,
        "exp_run_log_cryodrgn_version_title": None,
        "command_builder_page_mode": True,
    }

    with app.test_request_context("/command-builder"):
        html = render_template("command_builder.html", **ctx)
    return _adapt_html_for_command_builder_page(html, base_path, repo_url)


def build_command_builder_page_site(
    out_dir: os.PathLike[str] | str,
    *,
    base_path: str | None = None,
    repo_url: str | None = None,
) -> Path:
    """Write ``index.html``, static assets, and ``.nojekyll`` into ``out_dir``."""
    root = Path(out_dir)
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    (root / ".nojekyll").touch()
    html = render_command_builder_html(base_path=base_path, repo_url=repo_url)
    (root / "index.html").write_text(html, encoding="utf-8")
    static_src = _STATIC_DIR
    if static_src.is_dir():
        shutil.copytree(static_src, root / "static")
    return root


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build static GitHub Pages files for the cryoDRGN command builder."
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("site"),
        help="output directory (default: site)",
    )
    parser.add_argument(
        "--base-path",
        default=None,
        help=(
            "GitHub Pages project base path "
            "(default: /$REPO_NAME/ from GITHUB_REPOSITORY)"
        ),
    )
    parser.add_argument(
        "--repo-url",
        default=None,
        help=(
            "repository URL for nav link "
            "(default: https://github.com/$GITHUB_REPOSITORY)"
        ),
    )
    args = parser.parse_args(argv)
    build_command_builder_page_site(
        args.output,
        base_path=args.base_path,
        repo_url=args.repo_url,
    )
    print(f"Wrote GitHub Pages site to {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
