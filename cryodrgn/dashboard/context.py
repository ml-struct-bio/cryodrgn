"""Workdir / epoch resolution, caches, and Jinja template context injectors.

This module owns the long-lived per-app state: which run and epoch the user is
viewing, which workdirs the command-builder page should list, and the context
dictionaries used by every rendered template. The HTTP views in :mod:`app`
read everything through these helpers so cache invalidation stays centralised.
"""

from __future__ import annotations

import os
import shlex
from pathlib import Path

import yaml
from flask import Flask, current_app, g, jsonify, redirect, request, session, url_for

from cryodrgn.dashboard.command_builder_data import (
    COMMAND_BUILDER_REQUIRED_FIELD_TITLES,
    COMMAND_BUILDER_SCHEMA,
)
from cryodrgn.dashboard.data import DashboardExperiment, list_z_epochs, load_experiment

# ---------------------------------------------------------------------------
# Module-level caches (invalidated whenever the active workdir/epoch changes).
# ---------------------------------------------------------------------------

EXP_CACHE: dict[tuple[str, int, int], DashboardExperiment] = {}
# (epoch, kmeans, xcol, ycol, selection_tuple_or_None) -> (rows, images_b64, elapsed_s).
PRELOAD_CACHE: dict[
    tuple[int, int, str, str, tuple[int, ...] | None],
    tuple[list[int], list[str], float],
] = {}
_EPOCHS_BY_WORKDIR_CACHE: dict[str, list[int]] = {}


def clear_experiment_caches() -> None:
    """Drop cached experiments / preloads / graph neighbors across the process."""
    from cryodrgn.dashboard.trajectory import _TRAJ_GRAPH_NEIGHBOR_CACHE

    EXP_CACHE.clear()
    PRELOAD_CACHE.clear()
    _TRAJ_GRAPH_NEIGHBOR_CACHE.clear()


# Endpoints that require a loaded experiment; every other endpoint is reachable
# in command-builder-only mode too.
EXP_REQUIRED_ENDPOINTS = frozenset(
    {
        "abinit_builder_redirect",
        "filter_page_redirect",
        "api_save_selection",
        "explorer",
        "api_explorer_volume_media",
        "api_scatter",
        "latent_3d_page",
        "api_scatter3d_z",
        "api_latent3d_preview_png",
        "api_preview_montage",
        "api_preload_images",
        "pairplot_page",
        "api_pairplot",
        "api_save_pairplot_png",
        "trajectory_creator_page",
        "api_trajectory_volumes",
        "api_trajectory_coords",
        "api_trajectory_save_volumes",
        "api_trajectory_save_zpath",
        "api_trajectory_import_anchors",
        "api_list_server_files",
        "api_trajectory_kmeans_centers",
        "api_trajectory_random_indices",
        "api_default_trajectory_endpoints",
        "landscape_volpca_page",
        "api_landscape_volpca_meta",
        "api_landscape_volpca_scatter",
        "api_landscape_volpca_generate_animations",
        "api_landscape_volpca_save_animations",
    }
)


# ---------------------------------------------------------------------------
# Workdir discovery (command-builder-only mode)
# ---------------------------------------------------------------------------


def _config_has_cryodrgn_cmd(config: object) -> bool:
    if not isinstance(config, dict):
        return False
    cmd = config.get("cmd")
    if isinstance(cmd, str):
        return "cryodrgn" in cmd.lower()
    if isinstance(cmd, list):
        return any("cryodrgn" in str(part).lower() for part in cmd)
    return False


def discover_cryodrgn_workdirs(cwd: str) -> list[str]:
    """Direct subfolders of ``cwd`` with ``config.yaml`` recording a cryodrgn command."""
    out: list[str] = []
    base = Path(cwd)
    if not base.is_dir():
        return out
    for child in sorted(base.iterdir(), key=lambda p: p.name.lower()):
        if not child.is_dir():
            continue
        cfg = child / "config.yaml"
        if not cfg.is_file():
            continue
        try:
            with cfg.open("r", encoding="utf-8") as fh:
                parsed = yaml.safe_load(fh)
        except Exception:
            continue
        if _config_has_cryodrgn_cmd(parsed):
            out.append(str(child.resolve()))
    return out


def _workdir_options(abs_paths: list[str], base_cwd: str) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for wd in abs_paths:
        try:
            label = os.path.relpath(wd, base_cwd)
        except Exception:
            label = wd
        out.append({"value": wd, "label": label})
    return out


# ---------------------------------------------------------------------------
# Active-workdir / epoch resolution
# ---------------------------------------------------------------------------


def _sync_discovery_session_boot() -> None:
    """Drop stale session workdir when restarting with CWD discovery.

    Only applies when the server was started without an outdir and at least one
    output folder was discovered — avoids auto-resuming a previous run from the
    cookie.
    """
    boot = current_app.config.get("DASHBOARD_DISCOVERY_BOOT_ID")
    if not boot:
        return
    if session.get("dashboard_discovery_boot") == boot:
        return
    session.pop("dashboard_workdir", None)
    session.pop("dashboard_epoch", None)
    session["dashboard_discovery_boot"] = boot


def active_workdir(app: Flask) -> str | None:
    default_wd = app.config.get("DASHBOARD_WORKDIR")
    candidates = set(app.config.get("DASHBOARD_DISCOVERED_WORKDIRS", []))
    if default_wd:
        candidates.add(default_wd)
    selected = session.get("dashboard_workdir")
    if selected and selected in candidates:
        return str(selected)
    return default_wd


def epochs_for_workdir(workdir: str) -> list[int]:
    cached = _EPOCHS_BY_WORKDIR_CACHE.get(workdir)
    if cached is not None:
        return cached
    epochs = list_z_epochs(workdir)
    _EPOCHS_BY_WORKDIR_CACHE[workdir] = epochs
    return epochs


def resolve_epoch(app: Flask) -> int:
    wd = active_workdir(app)
    if not wd:
        return 0
    epochs = epochs_for_workdir(wd)
    if not epochs:
        raise RuntimeError("No z.N.pkl epochs in workdir.")
    # New server process: drop stale session epoch so ``cryodrgn dashboard -e N`` applies.
    boot = app.config.get("DASHBOARD_SESSION_BOOT_ID")
    if boot and session.get("dashboard_session_boot") != boot:
        session["dashboard_session_boot"] = boot
        session.pop("dashboard_epoch", None)
    sess = session.get("dashboard_epoch")
    if sess is not None:
        try:
            ep = int(sess)
        except (TypeError, ValueError):
            ep = None
        else:
            if ep in epochs:
                return ep
            return max(epochs)
    start = app.config.get("DASHBOARD_START_EPOCH")
    if isinstance(start, int) and start in epochs:
        return int(start)
    return max(epochs)


def get_dashboard_exp(app: Flask) -> DashboardExperiment:
    wd = active_workdir(app)
    if not wd:
        raise RuntimeError("No output directory selected.")
    ep = resolve_epoch(app)
    km = int(app.config["DASHBOARD_KMEANS"])
    key = (wd, ep, km)
    if key not in EXP_CACHE:
        EXP_CACHE[key] = load_experiment(wd, epoch=ep, kmeans=km)
    return EXP_CACHE[key]


def bind_dashboard_exp() -> None:
    """``before_request`` hook: attach the loaded experiment to ``flask.g``."""
    _sync_discovery_session_boot()
    wd = active_workdir(current_app)
    if not wd:
        if request.endpoint in EXP_REQUIRED_ENDPOINTS:
            return redirect(url_for("index"), code=302)
        return
    g.dashboard_exp = get_dashboard_exp(current_app)


# ---------------------------------------------------------------------------
# Workdir / epoch switching routes (wired under ``/api/set_*``)
# ---------------------------------------------------------------------------


def _request_json_dict() -> dict:
    data = request.get_json(force=True, silent=True)
    return data if isinstance(data, dict) else {}


def api_set_epoch():
    wd = active_workdir(current_app)
    if not wd:
        return jsonify(error="Select an output folder first."), 400
    data = _request_json_dict()
    raw_epoch = data.get("epoch")
    if raw_epoch is None:
        return jsonify(error="Invalid epoch."), 400
    try:
        ep = int(raw_epoch)
    except (TypeError, ValueError):
        return jsonify(error="Invalid epoch."), 400
    if ep not in epochs_for_workdir(wd):
        return jsonify(error="Epoch not available for this output folder."), 400
    session["dashboard_epoch"] = ep
    clear_experiment_caches()
    return jsonify(ok=True, epoch=ep)


def api_set_workdir():
    data = _request_json_dict()
    raw = data.get("workdir")
    candidates = set(current_app.config.get("DASHBOARD_DISCOVERED_WORKDIRS", []))

    if raw is None or (isinstance(raw, str) and not raw.strip()):
        if not current_app.config.get("COMMAND_BUILDER_ONLY"):
            return jsonify(error="Cannot clear output folder in this mode."), 400
        if not candidates:
            return jsonify(error="No output folders available."), 400
        session.pop("dashboard_workdir", None)
        session.pop("dashboard_epoch", None)
        boot = current_app.config.get("DASHBOARD_DISCOVERY_BOOT_ID")
        if boot:
            session["dashboard_discovery_boot"] = boot
        clear_experiment_caches()
        return jsonify(ok=True, workdir=None)

    requested = str(raw).strip()
    if requested not in candidates:
        return jsonify(error="Invalid output folder."), 400
    epochs = epochs_for_workdir(requested)
    if not epochs:
        return jsonify(error="No analyzed epochs found in selected output folder."), 400
    session["dashboard_workdir"] = requested
    session["dashboard_epoch"] = max(epochs)
    boot = current_app.config.get("DASHBOARD_DISCOVERY_BOOT_ID")
    if boot:
        session["dashboard_discovery_boot"] = boot
    clear_experiment_caches()
    return jsonify(ok=True, workdir=requested, epoch=max(epochs))


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------


def abbrev_middle(text: object, maxlen: int = 30) -> str:
    """Shorten a string with a middle Unicode ellipsis when longer than ``maxlen``."""
    s = "" if text is None else str(text)
    if len(s) <= maxlen:
        return s
    if maxlen < 4:
        return s[:maxlen]
    ell = "\u2026"
    inner = maxlen - len(ell)
    left = inner // 2
    return s[:left] + ell + s[-(inner - left) :]


def _cmd_argv_for_nav_display(cmd_parts: list[str]) -> list[str]:
    """Drop the filesystem path to the cryodrgn entrypoint; show ``cryodrgn <args>``."""
    if not cmd_parts:
        return []
    parts = [str(x) for x in cmd_parts]
    if len(parts) >= 3 and parts[1] == "-m":
        mod = parts[2]
        if mod == "cryodrgn" or mod.startswith("cryodrgn."):
            return ["cryodrgn", *parts[3:]]
    base0 = os.path.basename(parts[0])
    if base0 == "cryodrgn" or base0.lower().startswith("cryodrgn."):
        return ["cryodrgn", *parts[1:]]
    if len(parts) >= 2:
        base1 = os.path.basename(parts[1])
        if base1 == "cryodrgn" or base1.lower().startswith("cryodrgn."):
            return ["cryodrgn", *parts[2:]]
    return parts


def _abbrev_middle_token(text: str, maxlen: int = 120) -> str:
    s = "" if text is None else str(text)
    if len(s) <= maxlen:
        return s
    if maxlen < 4:
        return s[:maxlen]
    ell = "\u2026"
    inner = maxlen - len(ell)
    left = inner // 2
    return s[:left] + ell + s[-(inner - left) :]


def _argv_four_command_lines(argv: list[str]) -> list[str]:
    """Format argv as at most four lines for the nav ribbon.

    Line 1 is ``argv[0:2]`` (so it shows e.g. ``cryodrgn abinit``). Args longer
    than 120 characters are abbreviated with a middle-ellipsis.
    """

    def _display_join(tokens: list[str]) -> str:
        return " ".join(_abbrev_middle_token(t) for t in tokens)

    if not argv:
        return []
    if len(argv) == 1:
        return [_abbrev_middle_token(argv[0])]
    if len(argv) == 2:
        return [_display_join(argv)]

    head_tokens = argv[0:2]
    head = _display_join(head_tokens)
    rest = argv[2:]

    def _can_break_after(token: str) -> bool:
        """Allow a line break only after argument values or ``key=value`` pairs."""
        if "=" in token:
            return True
        return not token.startswith("-")

    def chunk_weight(chunk: list[str]) -> int:
        # Use original token lengths (not abbreviated) for balancing.
        if not chunk:
            return 0
        return sum(len(x) for x in chunk) + max(0, len(chunk) - 1)

    if len(rest) == 1:
        return [head, _abbrev_middle_token(rest[0])]
    if len(rest) == 2:
        if _can_break_after(rest[0]):
            return [
                head,
                _abbrev_middle_token(rest[0]),
                _abbrev_middle_token(rest[1]),
            ]
        return [head, _display_join(rest)]

    # Brute-force two cut points for similar chunk "weight" (rest is small).
    n = len(rest)
    avg = chunk_weight(rest) / 3.0
    best_score: float | None = None
    best_i = 1
    best_j = n - 1
    for i in range(1, n - 1):
        if not _can_break_after(rest[i - 1]):
            continue
        for j in range(i + 1, n):
            if not _can_break_after(rest[j - 1]):
                continue
            c1, c2, c3 = rest[:i], rest[i:j], rest[j:]
            if not c1 or not c2 or not c3:
                continue
            w1 = chunk_weight(c1)
            w2 = chunk_weight(c2)
            w3 = chunk_weight(c3)
            score = (w1 - avg) ** 2 + (w2 - avg) ** 2 + (w3 - avg) ** 2
            if best_score is None or score < best_score:
                best_score = score
                best_i = i
                best_j = j
    if best_score is None:
        best_i, best_j = 1, n - 1

    return [
        head,
        _display_join(rest[:best_i]),
        _display_join(rest[best_i:best_j]),
        _display_join(rest[best_j:]),
    ]


def command_builder_template_kwargs(
    exp: DashboardExperiment | None,
) -> dict[str, object]:
    """Template variables for ``command_builder.html`` from experiment config."""
    if exp is None:
        return {
            "default_particles": "",
            "default_ctf": "",
            "default_zdim": 8,
            "default_outdir_abinit": "abinit_run",
            "default_outdir_train": "train_next",
            "default_poses": "",
            "command_builder_schema": COMMAND_BUILDER_SCHEMA,
            "command_builder_required_field_titles": COMMAND_BUILDER_REQUIRED_FIELD_TITLES,
        }

    cfg = exp.train_configs
    da = cfg.get("dataset_args", {}) or {}
    ma = cfg.get("model_args", {}) or {}
    raw_p = da.get("particles") or ""
    default_particles = raw_p if isinstance(raw_p, str) else str(raw_p)
    raw_c = da.get("ctf")
    default_ctf = raw_c if isinstance(raw_c, str) else (str(raw_c) if raw_c else "")
    try:
        default_zdim = int(ma.get("zdim", 8))
    except (TypeError, ValueError):
        default_zdim = 8
    raw_poses = da.get("poses")
    if isinstance(raw_poses, str) and raw_poses.strip():
        default_poses = raw_poses.strip()
    else:
        default_poses = os.path.join(exp.workdir, f"pose.{exp.epoch}.pkl")
    return {
        "default_particles": default_particles,
        "default_ctf": default_ctf,
        "default_zdim": default_zdim,
        "default_outdir_abinit": os.path.join(exp.workdir, "abinit_run"),
        "default_outdir_train": os.path.join(exp.workdir, "train_next"),
        "default_poses": default_poses,
        "command_builder_schema": COMMAND_BUILDER_SCHEMA,
        "command_builder_required_field_titles": COMMAND_BUILDER_REQUIRED_FIELD_TITLES,
    }


# ---------------------------------------------------------------------------
# Jinja context processors
# ---------------------------------------------------------------------------


def inject_meta() -> dict:
    """Template context with a loaded experiment (default mode)."""
    e: DashboardExperiment = g.dashboard_exp
    discovered_workdirs = current_app.config.get("DASHBOARD_DISCOVERED_WORKDIRS", [])
    discovery_cwd = current_app.config.get("DASHBOARD_DISCOVERY_CWD", os.getcwd())
    epochs = epochs_for_workdir(e.workdir)
    cfg = e.train_configs
    cmd_list = cfg.get("cmd", [])
    zdim = cfg.get("model_args", {}).get("zdim", "?")
    raw_parts = [str(x) for x in cmd_list] if isinstance(cmd_list, list) else []
    cmd_parts = _cmd_argv_for_nav_display(raw_parts)
    if len(cmd_parts) > 1:
        model_type = cmd_parts[1]
    elif len(cmd_list) > 1:
        model_type = cmd_list[1]
    else:
        model_type = "unknown"
    cfg_train_command = shlex.join(cmd_parts) if cmd_parts else ""
    cfg_cmd_display_lines = (
        _argv_four_command_lines(cmd_parts) if cmd_parts else [f"{model_type} z{zdim}"]
    )
    return {
        "exp_workdir": e.workdir,
        "exp_epoch": e.epoch,
        "exp_kmeans": e.kmeans_folder_id,
        "filter_plot_inds_default": current_app.config["FILTER_PLOT_INDS"] or "",
        "dashboard_epochs": epochs if epochs else [e.epoch],
        "cfg_model_type": model_type,
        "cfg_zdim": zdim,
        "cfg_train_command": cfg_train_command,
        "cfg_cmd_display_lines": cfg_cmd_display_lines,
        "command_builder_only": False,
        "discovered_workdirs": discovered_workdirs,
        "discovered_workdir_options": _workdir_options(
            discovered_workdirs, discovery_cwd
        ),
        "selected_workdir": e.workdir,
    }


def inject_meta_command_builder_only() -> dict:
    """Template context for the command-builder-only launch mode.

    Falls back to :func:`inject_meta` once a workdir has been selected.
    """
    active_wd = active_workdir(current_app)
    if active_wd and hasattr(g, "dashboard_exp"):
        return inject_meta()
    discovered_workdirs = current_app.config.get("DASHBOARD_DISCOVERED_WORKDIRS", [])
    discovery_cwd = current_app.config.get("DASHBOARD_DISCOVERY_CWD", os.getcwd())
    epochs = epochs_for_workdir(active_wd) if active_wd else []
    exp_epoch = resolve_epoch(current_app) if active_wd and epochs else 0
    return {
        "exp_workdir": "",
        "exp_epoch": exp_epoch,
        "exp_kmeans": -1,
        "filter_plot_inds_default": "",
        "dashboard_epochs": epochs or [0],
        "cfg_model_type": "cryodrgn",
        "cfg_zdim": "",
        "cfg_train_command": "",
        "cfg_cmd_display_lines": ["No experiment loaded"],
        "command_builder_only": not bool(active_wd),
        "discovered_workdirs": discovered_workdirs,
        "discovered_workdir_options": _workdir_options(
            discovered_workdirs, discovery_cwd
        ),
        "selected_workdir": active_wd or "",
    }
