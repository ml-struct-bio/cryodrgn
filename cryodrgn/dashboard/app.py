"""Flask app for the cryoDRGN analysis dashboard."""

from __future__ import annotations

import base64
import io
import logging
import os
import pickle
import shlex
import time

import matplotlib.pyplot as plt
import numpy as np
from flask import (
    Flask,
    Response,
    current_app,
    g,
    jsonify,
    redirect,
    render_template,
    request,
    session,
    url_for,
)

from cryodrgn import utils
from cryodrgn.dashboard.data import (
    DashboardExperiment,
    list_z_epochs,
    load_experiment,
    particle_image_array,
)
from cryodrgn.dashboard.command_builder_data import (
    COMMAND_BUILDER_REQUIRED_FIELD_TITLES,
    COMMAND_BUILDER_SCHEMA,
)
from cryodrgn.dashboard.explorer_volumes import (
    DEFAULT_CHIMERAX_PARALLEL,
    DEFAULT_GIF_FRAMES,
    explorer_volumes_eligible,
    generate_montage_volume_pngs,
    generate_trajectory_volume_pngs,
    volume_cell_gif_from_cache,
)
from cryodrgn.dashboard.mpl_style import ezlab_matplotlib_rc
from cryodrgn.dashboard.plots import (
    normalize_continuous_palette,
    pair_grid_png,
    pair_grid_skeleton_placeholder_layout,
    scatter3d_z_json,
    scatter_json,
)

logger = logging.getLogger(__name__)

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_TEMPLATE_DIR = os.path.join(_THIS_DIR, "templates")
_STATIC_DIR = os.path.join(_THIS_DIR, "static")


def _encode_particle_batch(
    mrcfile: str,
    datadir: str | None,
    global_indices: list[int],
    max_px: int,
) -> list[str]:
    """Load and encode a batch of particle images as base64 JPEG strings.

    Designed to run in a ``ProcessPoolExecutor`` worker — all imports
    are local so the function is self-contained after a fork.
    """
    import base64 as _b64
    import io as _io

    import numpy as _np
    from PIL import Image as PILImage

    from cryodrgn.dataset import ImageDataset

    ds = ImageDataset(mrcfile=mrcfile, lazy=True, datadir=datadir)
    out: list[str] = []
    for gidx in global_indices:
        raw = ds.src.images(gidx, as_numpy=True)
        if raw.ndim == 3:
            raw = raw[0]
        arr = _np.asarray(raw, dtype=_np.float32)
        lo, hi = _np.percentile(arr, (2, 98))
        u8 = (_np.clip((arr - lo) / (hi - lo + 1e-9), 0, 1) * 255).astype(_np.uint8)
        pil = PILImage.fromarray(u8, mode="L")
        if max(pil.size) > max_px:
            pil = pil.resize((max_px, max_px), PILImage.LANCZOS)
        buf = _io.BytesIO()
        pil.save(buf, format="JPEG", quality=85)
        out.append(_b64.b64encode(buf.getvalue()).decode("ascii"))
    return out


# (epoch, kmeans) -> loaded experiment (invalidated when epoch changes in session).
_EXP_CACHE: dict[tuple[int, int], DashboardExperiment] = {}
# (epoch, kmeans, xcol, ycol, selection_tuple_or_None) -> (rows, images_b64, elapsed_s).
_PRELOAD_CACHE: dict[
    tuple[int, int, str, str, tuple[int, ...] | None],
    tuple[list[int], list[str], float],
] = {}


def _filter_ui_scatter_max_points() -> int:
    """Max points when ``filter_ui=1`` on ``/api/scatter`` (env ``CRYODRGN_DASHBOARD_FILTER_MAX_POINTS``)."""
    try:
        v = int(os.environ.get("CRYODRGN_DASHBOARD_FILTER_MAX_POINTS", "500000"))
    except ValueError:
        v = 500_000
    return max(50_000, min(v, 2_000_000))


def _default_xy_cols(cols: list[str]) -> tuple[str, str]:
    """Pick sensible default X/Y axes (UMAP if available, else first two)."""
    x = "UMAP1" if "UMAP1" in cols else cols[0]
    y = "UMAP2" if "UMAP2" in cols else cols[min(1, len(cols) - 1)]
    return x, y


def _trajectory_default_xy_cols(cols: list[str]) -> tuple[str, str]:
    """Trajectory creator defaults to PC1 vs PC2 when present, else UMAP / first columns."""
    if "PC1" in cols and "PC2" in cols:
        return "PC1", "PC2"
    return _default_xy_cols(cols)


def _covariate_display_name(name: str) -> str:
    """Human-friendly covariate names in dashboard selectors."""
    if name == "labels":
        return "k-means labels"
    return name


def _montage_bytes(exp: DashboardExperiment, row_indices: list[int]) -> bytes:
    rows = [int(r) for r in row_indices[:25]]
    with ezlab_matplotlib_rc():
        if not rows:
            fig, ax = plt.subplots(figsize=(4, 4))
            ax.text(
                0.5, 0.5, "Hover points to preview particles", ha="center", va="center"
            )
            ax.axis("off")
            buf = io.BytesIO()
            fig.savefig(buf, format="png", bbox_inches="tight")
            plt.close(fig)
            buf.seek(0)
            return buf.getvalue()

        imgs = [particle_image_array(exp, r) for r in rows]
        n = len(imgs)
        ncol = int(np.ceil(n**0.5))
        nrow = int(np.ceil(n / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(1.8 * ncol, 1.8 * nrow))
        axes_flat = np.atleast_1d(axes).ravel()
        for ax in axes_flat[n:]:
            ax.axis("off")
        for i, img in enumerate(imgs):
            ax = axes_flat[i]
            lo, hi = np.percentile(img, (2, 98))
            disp = np.clip((img - lo) / (hi - lo + 1e-9), 0, 1)
            ax.imshow(disp, cmap="gray")
            ax.set_title(f"idx {int(exp.plot_df.iloc[rows[i]]['index'])}", fontsize=8)
            ax.axis("off")
        plt.tight_layout()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight", dpi=120)
        plt.close(fig)
        buf.seek(0)
        return buf.getvalue()


def _particle_thumbnail_b64_from_row(
    exp: DashboardExperiment, row_index: int, max_side: int = 72
) -> str | None:
    """Small greyscale JPEG (base64) for trajectory NN inset previews."""
    if not exp.can_preview_particles:
        return None
    try:
        from PIL import Image

        img = particle_image_array(exp, int(row_index))
        lo, hi = np.percentile(img, (2, 98))
        disp = np.clip((img - lo) / (hi - lo + 1e-9), 0, 1)
        u8 = (disp * 255).astype(np.uint8)
        pil = Image.fromarray(u8, mode="L")
        w, h = pil.size
        s = max(w, h)
        if s > max_side:
            scale = max_side / float(s)
            pil = pil.resize(
                (max(1, int(w * scale)), max(1, int(h * scale))),
                Image.LANCZOS,
            )
        buf = io.BytesIO()
        pil.save(buf, format="JPEG", quality=88)
        return base64.standard_b64encode(buf.getvalue()).decode("ascii")
    except Exception:
        logger.debug("particle thumbnail encode failed", exc_info=True)
        return None


def _plot_row_particle_index(exp: DashboardExperiment, row_index: int) -> int:
    ri = int(row_index)
    if "index" in exp.plot_df.columns:
        return int(exp.plot_df.iloc[ri]["index"])
    return ri


def _round_direct_mode_traj_xy(traj_xy: np.ndarray) -> np.ndarray:
    """Snap axis coordinates to a short decimal representation (2 or 3 places per axis)."""
    out = np.asarray(traj_xy, dtype=np.float64).copy()
    for j in range(out.shape[1]):
        col = out[:, j]
        finite = col[np.isfinite(col)]
        if finite.size == 0:
            continue
        mx = float(np.nanmax(np.abs(finite)))
        decimals = 2 if mx >= 100.0 else 3
        out[:, j] = np.round(col, decimals)
    return out


def _plot_df_rows_for_dataset_indices(
    exp: DashboardExperiment, dataset_indices: np.ndarray
) -> list[int]:
    """Map saved dataset indices (values in ``indices.pkl``) to ``plot_df`` row positions."""
    want_arr = np.asarray(dataset_indices).ravel().astype(int, copy=False)
    if want_arr.size == 0:
        return []
    ai = np.asarray(exp.all_indices)
    mask = np.isin(ai, want_arr, assume_unique=False)
    return np.nonzero(mask)[0].astype(int).tolist()


def _load_plot_df_rows_from_plot_inds_file(
    exp: DashboardExperiment, plot_inds_path: str | None
) -> list[int]:
    """Map ``--plot-inds`` pickle dataset indices to ``plot_df`` rows; [] on missing path or error."""
    path = (plot_inds_path or "").strip()
    if not path or not os.path.isfile(path):
        return []
    try:
        pre = utils.load_pkl(path)
        return _plot_df_rows_for_dataset_indices(exp, np.asarray(pre))
    except Exception as err:
        logger.warning("Could not load plot-inds from %s: %s", path, err)
        return []


def _has_umap_columns(exp: DashboardExperiment) -> bool:
    return (
        exp.umap is not None
        and "UMAP1" in exp.plot_df.columns
        and "UMAP2" in exp.plot_df.columns
    )


def _has_pc_columns(exp: DashboardExperiment) -> bool:
    return "PC1" in exp.plot_df.columns and "PC2" in exp.plot_df.columns


def _parse_preselect_rows_param(raw: str | None) -> tuple[list[int] | None, str | None]:
    """Parse ``preselect_rows`` query (comma-separated ints). Returns (rows, error_message)."""
    s = (raw or "").strip()
    if not s:
        return None, None
    try:
        return [int(p) for p in s.split(",") if p.strip()][:5000], None
    except ValueError as exc:
        return None, f"invalid preselect_rows: {exc}"


def _command_builder_template_kwargs(exp: DashboardExperiment) -> dict[str, object]:
    """Template variables for ``command_builder.html`` from experiment config."""
    cfg = exp.train_configs
    da = cfg.get("dataset_args", {}) or {}
    ma = cfg.get("model_args", {}) or {}
    raw_p = da.get("particles") or ""
    default_particles = raw_p if isinstance(raw_p, str) else str(raw_p)
    raw_c = da.get("ctf")
    default_ctf = raw_c if isinstance(raw_c, str) else (str(raw_c) if raw_c else "")
    zd = ma.get("zdim", 8)
    try:
        default_zdim = int(zd)
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


def _redirect(endpoint: str):
    return redirect(url_for(endpoint), code=302)


def _resolve_epoch(app: Flask) -> int:
    epochs: list[int] = app.config["DASHBOARD_EPOCHS"]
    if not epochs:
        raise RuntimeError("No z.N.pkl epochs in workdir.")
    sess = session.get("dashboard_epoch")
    if sess is None:
        return int(app.config["DASHBOARD_START_EPOCH"])
    try:
        ep = int(sess)
    except (TypeError, ValueError):
        return int(app.config["DASHBOARD_START_EPOCH"])
    if ep not in epochs:
        return max(epochs)
    return ep


def _get_dashboard_exp(app: Flask) -> DashboardExperiment:
    ep = _resolve_epoch(app)
    km = int(app.config["DASHBOARD_KMEANS"])
    key = (ep, km)
    if key not in _EXP_CACHE:
        _EXP_CACHE[key] = load_experiment(
            app.config["DASHBOARD_WORKDIR"],
            epoch=ep,
            kmeans=km,
        )
    return _EXP_CACHE[key]


def _stratified_xy_row_indices(
    coords: np.ndarray,
    rng: np.random.Generator,
    total_k: int,
) -> set[int]:
    """Local row indices 0..m-1: random + high mean dist to refs + high NN dist (cf. preload)."""
    from scipy.spatial import cKDTree
    from scipy.spatial.distance import cdist as _cdist

    m = int(coords.shape[0])
    if m == 0:
        return set()
    total_k = min(m, total_k)
    k_random = total_k * 2 // 3
    k_mean = total_k // 6
    random_rows = set(
        rng.choice(m, size=min(k_random, m), replace=False).tolist(),
    )

    n_ref = min(m, 500)
    ref_coords = coords[rng.choice(m, size=n_ref, replace=False)]
    avg_dists = np.zeros(m)
    batch = 20_000
    for start in range(0, m, batch):
        avg_dists[start : start + batch] = _cdist(
            coords[start : start + batch],
            ref_coords,
        ).mean(axis=1)

    mean_rows: set[int] = set()
    for idx in np.argsort(-avg_dists):
        if len(mean_rows) >= k_mean:
            break
        r = int(idx)
        if r not in random_rows:
            mean_rows.add(r)

    tree = cKDTree(coords)
    nn_dists = tree.query(coords, k=2)[0][:, 1]

    k_nn = max(0, total_k - k_random - len(mean_rows))
    excluded = random_rows | mean_rows
    nn_rows: set[int] = set()
    for idx in np.argsort(-nn_dists):
        if len(nn_rows) >= k_nn:
            break
        r = int(idx)
        if r not in excluded:
            nn_rows.add(r)

    return random_rows | mean_rows | nn_rows


def _sample_plot_df_rows_for_preload(
    exp: DashboardExperiment,
    xcol: str,
    ycol: str,
    *,
    restrict_to_rows: list[int] | None = None,
) -> tuple[list[int], list[int]]:
    """Pick up to 5k plot_df rows for thumbnail preload (random + spread-out in XY).

    If ``restrict_to_rows`` is set, only those plot_df indices are eligible.
    """
    n = len(exp.plot_df)
    coords_full = exp.plot_df[[xcol, ycol]].values.astype(np.float64)
    rng = np.random.default_rng(42)

    if restrict_to_rows is not None:
        sub_idx = sorted({int(r) for r in restrict_to_rows if 0 <= int(r) < n})
        if not sub_idx:
            return [], []
        coords = coords_full[np.array(sub_idx, dtype=int)]
        inv_map = sub_idx
    else:
        coords = coords_full
        inv_map = list(range(n))

    m = len(coords)
    total_k = min(m, 5000)
    local_sel = _stratified_xy_row_indices(coords, rng, total_k)
    plot_rows_set = {inv_map[i] for i in local_sel}
    pairs = sorted(
        ((r, int(exp.all_indices[r])) for r in plot_rows_set),
        key=lambda p: p[1],
    )
    plot_rows = [p[0] for p in pairs]
    global_indices = [p[1] for p in pairs]
    return plot_rows, global_indices


def _preload_cache_time_estimate_bounds(cpus: int) -> tuple[int, int]:
    """Rough wall-clock range (seconds) for caching up to 5k 96px JPEG thumbs.

    Scales inversely with ``cpus`` (preload worker count), anchored at ~20–30s
    with 4 cores (empirical). I/O and fixed overhead keep a modest floor at high
    core counts.
    """
    cpus = max(1, int(cpus))
    low = max(8, int(round(20 * 4 / cpus)))
    high_parallel = int(round(30 * 4 / cpus))
    high = max(max(18, low + 6), high_parallel)
    return low, high


def _format_preload_cache_time_hint(cpus: int) -> str:
    lo, hi = _preload_cache_time_estimate_bounds(cpus)
    cpus = max(1, int(cpus))
    if lo == hi:
        span = f"~{lo}s"
    else:
        span = f"~{lo}\u2013{hi}s"
    cword = "core" if cpus == 1 else "cores"
    return f"{span} for {cpus} {cword}"


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
    right = inner - left
    return s[:left] + ell + s[-right:]


def _bind_dashboard_exp() -> None:
    g.dashboard_exp = _get_dashboard_exp(current_app)


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


def _argv_four_command_lines(argv: list[str]) -> list[str]:
    """Format argv as at most four lines.

    Line 1 is argv[0:2] (so it shows e.g. ``cryodrgn abinit``).
    Args longer than 120 characters are abbreviated with an ellipsis-in-the-middle.
    """

    def _abbrev_middle_token(text: str, maxlen: int = 120) -> str:
        s = "" if text is None else str(text)
        if len(s) <= maxlen:
            return s
        if maxlen < 4:
            return s[:maxlen]
        ell = "…"
        inner = maxlen - len(ell)
        left = inner // 2
        right = inner - left
        return s[:left] + ell + s[-right:]

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
        """Allow a line break only after argument values or key=value pairs."""
        if "=" in token:
            return True
        # Argument "values" are tokens not starting with "-" / "--".
        return not token.startswith("-")

    def chunk_weight(chunk: list[str]) -> int:
        # Keep line-splitting based on original token lengths (not abbreviated display),
        # so we still get "similar amounts of space".
        if not chunk:
            return 0
        return sum(len(x) for x in chunk) + max(0, len(chunk) - 1)

    # If the tail is tiny, keep it compact (fewer than 4 lines).
    if len(rest) == 1:
        return [head, _abbrev_middle_token(rest[0])]
    if len(rest) == 2:
        # Only split into 3 lines if the first tail token is a valid break point.
        if _can_break_after(rest[0]):
            return [
                head,
                _abbrev_middle_token(rest[0]),
                _abbrev_middle_token(rest[1]),
            ]
        # Otherwise keep both tail tokens on the same line.
        return [head, _display_join(rest)]

    # Split rest into three chunks (lines 2-4) with similar token "weight".
    # Brute force two cut points; `rest` is small (typically CLI args), so O(n^2) is fine.
    n = len(rest)
    avg = chunk_weight(rest) / 3.0
    best_score: float | None = None
    best_i = 1
    best_j = n - 1

    # Cut points i and j are boundaries:
    # - boundary at i breaks after rest[i-1] (so rest[i-1] must be breakable)
    # - boundary at j breaks after rest[j-1] (so rest[j-1] must be breakable)
    for i in range(1, n - 1):
        if not _can_break_after(rest[i - 1]):
            continue
        for j in range(i + 1, n):
            if not _can_break_after(rest[j - 1]):
                continue
            c1 = rest[:i]
            c2 = rest[i:j]
            c3 = rest[j:]
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

    # If no feasible cut points exist (e.g. tail is all flags), fall back to the
    # unconstrained split so we still get a readable display.
    if best_score is None:
        best_score = 0.0
        best_i = 1
        best_j = n - 1

    mid1 = _display_join(rest[:best_i])
    mid2 = _display_join(rest[best_i:best_j])
    last = _display_join(rest[best_j:])
    return [head, mid1, mid2, last]


def inject_meta():
    e: DashboardExperiment = g.dashboard_exp
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
        "dashboard_epochs": current_app.config["DASHBOARD_EPOCHS"],
        "cfg_model_type": model_type,
        "cfg_zdim": zdim,
        "cfg_train_command": cfg_train_command,
        "cfg_cmd_display_lines": cfg_cmd_display_lines,
    }


def api_set_epoch():
    data = request.get_json(force=True, silent=True) or {}
    raw_epoch = data.get("epoch")
    if raw_epoch is None:
        return jsonify(error="Invalid epoch."), 400
    try:
        ep = int(raw_epoch)
    except (TypeError, ValueError):
        return jsonify(error="Invalid epoch."), 400
    if ep not in current_app.config["DASHBOARD_EPOCHS"]:
        return jsonify(error="Epoch not available for this output folder."), 400
    session["dashboard_epoch"] = ep
    _EXP_CACHE.clear()
    _PRELOAD_CACHE.clear()
    return jsonify(ok=True, epoch=ep)


def index():
    e: DashboardExperiment = g.dashboard_exp
    return render_template(
        "index.html",
        can_images=e.can_preview_particles,
        zdim=int(e.z.shape[1]),
        show_trajectory_creator=explorer_volumes_eligible(e),
    )


def command_builder_page():
    e: DashboardExperiment = g.dashboard_exp
    return render_template(
        "command_builder.html",
        **_command_builder_template_kwargs(e),
    )


def abinit_builder_redirect():
    return _redirect("command_builder_page")


def filter_page_redirect():
    return _redirect("explorer")


def api_save_selection():
    e: DashboardExperiment = g.dashboard_exp
    data = request.get_json(force=True, silent=True) or {}
    rows_raw = data.get("rows")
    if not isinstance(rows_raw, list) or len(rows_raw) == 0:
        return jsonify(error="No particles selected."), 400
    try:
        rows = sorted({int(r) for r in rows_raw})
    except (TypeError, ValueError):
        return jsonify(error="Invalid row list."), 400
    n_df = len(e.plot_df)
    if any(r < 0 or r >= n_df for r in rows):
        return jsonify(error="Row index out of range."), 400
    selected_ds = np.asarray(e.all_indices[rows], dtype=int)
    force = bool(data.get("force"))
    save_inverse = bool(data.get("save_inverse"))
    basename = (data.get("basename") or "indices").strip() or "indices"
    basename = os.path.basename(basename)
    if force:
        out_base = os.path.join(e.workdir, "indices")
    else:
        sel_dir = (data.get("sel_dir") or "").strip()
        dir_abs = os.path.abspath(sel_dir) if sel_dir else os.path.abspath(e.workdir)
        if not os.path.isdir(dir_abs):
            return jsonify(error=f"Directory does not exist: {dir_abs}"), 400
        out_base = os.path.join(dir_abs, basename)
    path_main = out_base + ".pkl"
    path_inv = out_base + "_inverse.pkl"
    parent = os.path.dirname(path_main)
    if parent:
        os.makedirs(parent, exist_ok=True)
    inverse_ds = np.setdiff1d(np.asarray(e.all_indices, dtype=int), selected_ds)
    try:
        with open(path_main, "wb") as fh:
            pickle.dump(selected_ds, fh)
        if save_inverse:
            with open(path_inv, "wb") as fh:
                pickle.dump(inverse_ds, fh)
    except OSError as err:
        logger.exception("save selection failed")
        return jsonify(error=str(err)), 500
    payload = {
        "ok": True,
        "path": path_main,
        "n_selected": int(selected_ds.size),
    }
    if save_inverse:
        payload["inverse_path"] = path_inv
    return jsonify(payload)


def explorer():
    e: DashboardExperiment = g.dashboard_exp
    if not e.can_preview_particles:
        return (
            render_template(
                "no_images.html",
                reason="Tilt-series data does not include single-particle thumbnails in this view.",
            ),
            200,
        )
    cols = e.numeric_columns
    dx, dy = _default_xy_cols(cols)
    initial_rows = _load_plot_df_rows_from_plot_inds_file(
        e,
        current_app.config.get("FILTER_PLOT_INDS"),
    )
    pc = int(current_app.config["PRELOAD_CPUS"])
    return render_template(
        "scatter_explorer.html",
        numeric_cols=cols,
        covariate_display_map={c: _covariate_display_name(c) for c in cols},
        default_x=dx,
        default_y=dy,
        initial_rows=initial_rows,
        total_particles=int(len(e.all_indices)),
        workdir=e.workdir,
        preload_cache_time_hint=_format_preload_cache_time_hint(pc),
        show_volume_explorer=explorer_volumes_eligible(e),
    )


def api_explorer_volume_media():
    e: DashboardExperiment = g.dashboard_exp
    if not explorer_volumes_eligible(e):
        return (
            jsonify(
                error="Volume explorer needs a CUDA GPU, single-particle data, "
                "and weights for the current epoch.",
            ),
            400,
        )
    data = request.get_json(force=True, silent=True) or {}
    rows_raw = data.get("rows")
    if not isinstance(rows_raw, list) or len(rows_raw) == 0:
        return jsonify(error="No montage rows supplied."), 400
    try:
        rows = [int(r) for r in rows_raw]
    except (TypeError, ValueError):
        return jsonify(error="Invalid row indices."), 400
    if len(rows) > 64:
        return jsonify(error="At most 64 montage cells are supported."), 400
    mode = str(data.get("mode") or "static").strip().lower()
    try:
        if mode == "static":
            blobs, cache_token = generate_montage_volume_pngs(e, rows)
            fmt = "png"
            b64s = [base64.standard_b64encode(b).decode("ascii") for b in blobs]
            return jsonify(
                ok=True,
                format=fmt,
                images=b64s,
                rows=rows,
                volume_cache_id=cache_token,
            )
        if mode in ("animate", "gif"):
            cache_id = data.get("volume_cache_id")
            if not cache_id or not isinstance(cache_id, str):
                return (
                    jsonify(
                        error="Animate requires volume_cache_id from Generate volumes.",
                    ),
                    400,
                )
            raw_ci = data.get("cell_index")
            if raw_ci is None:
                return jsonify(error="cell_index is required for animate."), 400
            try:
                cell_index = int(raw_ci)
            except (TypeError, ValueError):
                return jsonify(error="cell_index must be an integer."), 400
            gf = int(data.get("gif_frames", DEFAULT_GIF_FRAMES))
            cc = int(data.get("chimerax_cpus", DEFAULT_CHIMERAX_PARALLEL))
            gif_bytes = volume_cell_gif_from_cache(
                cache_id,
                cell_index,
                rows_expected=tuple(rows),
                gif_frames=gf,
                chimerax_cpus=cc,
            )
            b64_one = base64.standard_b64encode(gif_bytes).decode("ascii")
            return jsonify(
                ok=True,
                format="gif",
                image=b64_one,
                cell_index=cell_index,
                rows=rows,
            )
        return jsonify(error='mode must be "static" or "animate".'), 400
    except EnvironmentError as err:
        return jsonify(error=str(err), need_chimerax=True), 503
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("explorer volume media failed")
        return jsonify(error=str(err)), 500


def api_scatter():
    e: DashboardExperiment = g.dashboard_exp
    xcol = request.args.get("x", e.numeric_columns[0])
    ycol = request.args.get("y", e.numeric_columns[0])
    ccol = request.args.get("color") or "none"
    filter_ui = request.args.get("filter_ui") == "1"
    full = request.args.get("full") == "1"
    if xcol not in e.plot_df.columns or ycol not in e.plot_df.columns:
        return jsonify(error="bad axis column"), 400
    if ccol != "none" and ccol not in e.plot_df.columns:
        return jsonify(error="bad color column"), 400
    if filter_ui:
        max_pts = _filter_ui_scatter_max_points()
    elif full:
        max_pts = None
    else:
        max_pts = 200_000
    preselect_rows, pre_err = _parse_preselect_rows_param(
        request.args.get("preselect_rows"),
    )
    if pre_err:
        return jsonify(error=pre_err), 400
    marker_size = 4.0
    raw_ms = request.args.get("marker_size")
    if raw_ms:
        try:
            marker_size = max(0.5, min(float(raw_ms), 20))
        except ValueError:
            pass
    try:
        js = scatter_json(
            e,
            xcol,
            ycol,
            None if ccol == "none" else ccol,
            max_points=max_pts,
            preselect_plot_df_rows=preselect_rows,
            use_webgl=False,
            marker_size=marker_size,
            continuous_palette=request.args.get("palette"),
        )
    except Exception as err:
        logger.exception("scatter plot failed")
        return jsonify(error=str(err)), 500
    return Response(js, mimetype="application/json")


def latent_3d_page():
    e: DashboardExperiment = g.dashboard_exp
    zdim = int(e.z.shape[1])
    if zdim < 3:
        return (
            render_template(
                "pair_grid_need_more_cols.html",
                kind="z3",
                n=zdim,
            ),
            200,
        )
    z_cols = [f"z{i}" for i in range(zdim)]
    cols = e.numeric_columns
    return render_template(
        "latent_3d.html",
        z_cols=z_cols,
        numeric_cols=cols,
        covariate_display_map={c: _covariate_display_name(c) for c in cols},
        default_x="z0",
        default_y="z1",
        default_z="z2",
    )


def api_scatter3d_z():
    e: DashboardExperiment = g.dashboard_exp
    xcol = request.args.get("x", "z0")
    ycol = request.args.get("y", "z1")
    zcol = request.args.get("z", "z2")
    ccol = request.args.get("color") or "none"
    if ccol != "none" and ccol not in e.plot_df.columns:
        return jsonify(error="bad color column"), 400
    try:
        js = scatter3d_z_json(
            e,
            xcol,
            ycol,
            zcol,
            None if ccol == "none" else ccol,
            continuous_palette=request.args.get("palette"),
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("3-D Latent Space Visualizer failed")
        return jsonify(error=str(err)), 500
    return Response(js, mimetype="application/json")


def api_preview_montage():
    e: DashboardExperiment = g.dashboard_exp
    raw = request.args.get("rows", "")
    if not raw.strip():
        data = _montage_bytes(e, [])
    else:
        parts = [p.strip() for p in raw.split(",") if p.strip()]
        try:
            idxs = [int(p) for p in parts]
        except ValueError:
            return jsonify(error="rows must be integers"), 400
        data = _montage_bytes(e, idxs)
    return Response(data, mimetype="image/png")


def api_preload_images():
    """Return a subsample of particle images as base64 JPEGs for the viewer.

    Selection: ~2/3 random, ~1/6 high mean distance to reference points,
    ~1/6 high nearest-neighbour distance (see ``_sample_plot_df_rows_for_preload``).

    Use **POST** with a JSON body when ``selected_rows`` is large (lasso selections);
    query strings hit proxy/server URI length limits (414).
    """
    e: DashboardExperiment = g.dashboard_exp
    if not e.can_preview_particles:
        return jsonify(rows=[], images=[], elapsed=0)

    cols = e.plot_df.columns
    restrict_list: list[int] | None = None

    if request.method == "POST":
        data = request.get_json(force=True, silent=True) or {}
        xcol = str(data.get("x") or "")
        ycol = str(data.get("y") or "")
        raw_list = data.get("selected_rows")
        if raw_list is None:
            restrict_list = None
        elif not isinstance(raw_list, list):
            return jsonify(error="selected_rows must be a JSON array of integers."), 400
        else:
            try:
                restrict_list = [int(r) for r in raw_list]
            except (TypeError, ValueError):
                return (
                    jsonify(error="selected_rows must be a JSON array of integers."),
                    400,
                )
            if not restrict_list:
                restrict_list = None
    else:
        xcol = request.args.get("x", "")
        ycol = request.args.get("y", "")
        raw_sel = (request.args.get("selected_rows") or "").strip()
        if raw_sel:
            try:
                restrict_list = [
                    int(x.strip()) for x in raw_sel.split(",") if x.strip()
                ]
            except ValueError:
                return (
                    jsonify(error="selected_rows must be comma-separated integers."),
                    400,
                )
            if not restrict_list:
                restrict_list = None

    if xcol not in cols:
        xcol = str(cols[0])
    if ycol not in cols:
        ycol = str(cols[min(1, len(cols) - 1)])
    sel_key: tuple[int, ...] | None = (
        tuple(sorted(set(restrict_list))) if restrict_list else None
    )

    key = (e.epoch, e.kmeans_folder_id, xcol, ycol, sel_key)
    if key in _PRELOAD_CACHE:
        rows, imgs, elapsed = _PRELOAD_CACHE[key]
        return jsonify(rows=rows, images=imgs, elapsed=elapsed)

    t0 = time.monotonic()
    cpus = int(current_app.config.get("PRELOAD_CPUS") or 4)
    rows, global_indices = _sample_plot_df_rows_for_preload(
        e, xcol, ycol, restrict_to_rows=restrict_list
    )
    if not global_indices:
        return jsonify(rows=[], images=[], elapsed=0.0)

    if cpus > 1 and len(global_indices) > cpus:
        from concurrent.futures import ProcessPoolExecutor

        chunk_sz = -(-len(global_indices) // cpus)
        chunks = [
            global_indices[i : i + chunk_sz]
            for i in range(0, len(global_indices), chunk_sz)
        ]
        with ProcessPoolExecutor(max_workers=len(chunks)) as pool:
            futures = [
                pool.submit(
                    _encode_particle_batch,
                    e.particles_path,
                    e.datadir,
                    ch,
                    96,
                )
                for ch in chunks
            ]
            imgs: list[str] = []
            for f in futures:
                imgs.extend(f.result())
    else:
        imgs = _encode_particle_batch(
            e.particles_path,
            e.datadir,
            global_indices,
            96,
        )

    elapsed = round(time.monotonic() - t0, 1)
    _PRELOAD_CACHE[key] = (rows, imgs, elapsed)
    return jsonify(rows=rows, images=imgs, elapsed=elapsed)


def pairplot_page():
    e: DashboardExperiment = g.dashboard_exp
    zdim = int(e.z.shape[1])
    if zdim < 2:
        return (
            render_template(
                "pair_grid_need_more_cols.html",
                kind="zdim",
                n=zdim,
            ),
            200,
        )
    has_pc = _has_pc_columns(e)
    has_umap = _has_umap_columns(e)
    if not has_pc:
        return (
            render_template(
                "pair_grid_need_more_cols.html",
                kind="pca",
                n=zdim,
            ),
            200,
        )
    z_names = {f"z{i}" for i in range(zdim)}
    color_choices = [c for c in e.numeric_columns if c not in z_names]
    if not color_choices:
        return (
            render_template(
                "pair_grid_need_more_cols.html",
                kind="numeric",
                n=0,
            ),
            200,
        )
    default_color = next(
        (c for c in ("labels", "znorm", "UMAP1", "PC1") if c in color_choices),
        color_choices[0],
    )
    skeleton_cells = pair_grid_skeleton_placeholder_layout(zdim)
    return render_template(
        "pair_grid.html",
        color_choices=color_choices,
        covariate_display_map={c: _covariate_display_name(c) for c in color_choices},
        default_color=default_color,
        has_umap=has_umap,
        zdim=zdim,
        skeleton_placeholder_cells=skeleton_cells,
        pairplot_save_default_name="zdim_pairplot.png",
        pairplot_save_default_dir=os.path.join(e.workdir, f"analyze.{e.epoch}"),
    )


def api_pairplot():
    e: DashboardExperiment = g.dashboard_exp
    payload = request.get_json(force=True, silent=True) or {}
    color_col = payload.get("color_col") or payload.get("lower_color_col")
    if not color_col or not isinstance(color_col, str):
        return jsonify(error="Choose a color covariate."), 400
    if color_col not in e.plot_df.columns:
        return jsonify(error="Invalid color column."), 400
    if color_col in {f"z{i}" for i in range(int(e.z.shape[1]))}:
        return (
            jsonify(error="Latent z columns cannot be used as the color covariate."),
            400,
        )
    raw_diag = payload.get("diagonal_emb")
    if raw_diag is None or (isinstance(raw_diag, str) and raw_diag.strip() == ""):
        diagonal_emb = "umap" if _has_umap_columns(e) else "pc"
    else:
        diagonal_emb = str(raw_diag).lower()
    upper_style = (payload.get("upper_style") or "scatter").lower()
    raw_palette = payload.get("palette")
    pair_palette = normalize_continuous_palette(
        str(raw_palette) if raw_palette is not None else None,
    )
    if diagonal_emb not in ("pc", "umap"):
        return jsonify(error="diagonal_emb must be pc or umap."), 400
    if upper_style not in ("scatter", "hex"):
        return jsonify(error="upper_style must be scatter or hex."), 400
    if diagonal_emb == "umap" and not _has_umap_columns(e):
        return jsonify(error="UMAP is not available for this run."), 400
    if diagonal_emb == "pc" and not _has_pc_columns(e):
        return jsonify(error="PCA components are not available."), 400
    try:
        png, cells = pair_grid_png(
            e,
            lower_color_col=color_col,
            diagonal_emb=diagonal_emb,
            upper_style=upper_style,
            continuous_palette=pair_palette,
        )
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("pair grid failed")
        return jsonify(error=str(err)), 500
    b64 = base64.standard_b64encode(png).decode("ascii")
    return jsonify(png_b64=b64, cells=cells)


def api_save_pairplot_png():
    e: DashboardExperiment = g.dashboard_exp
    payload = request.get_json(force=True, silent=True) or {}
    color_col = payload.get("color_col") or payload.get("lower_color_col")
    if not color_col or not isinstance(color_col, str):
        return jsonify(error="Choose a color covariate."), 400
    if color_col not in e.plot_df.columns:
        return jsonify(error="Invalid color column."), 400
    if color_col in {f"z{i}" for i in range(int(e.z.shape[1]))}:
        return (
            jsonify(error="Latent z columns cannot be used as the color covariate."),
            400,
        )
    raw_diag = payload.get("diagonal_emb")
    if raw_diag is None or (isinstance(raw_diag, str) and raw_diag.strip() == ""):
        diagonal_emb = "umap" if _has_umap_columns(e) else "pc"
    else:
        diagonal_emb = str(raw_diag).lower()
    upper_style = (payload.get("upper_style") or "scatter").lower()
    if diagonal_emb not in ("pc", "umap"):
        return jsonify(error="diagonal_emb must be pc or umap."), 400
    if upper_style not in ("scatter", "hex"):
        return jsonify(error="upper_style must be scatter or hex."), 400
    if diagonal_emb == "umap" and not _has_umap_columns(e):
        return jsonify(error="UMAP is not available for this run."), 400
    if diagonal_emb == "pc" and not _has_pc_columns(e):
        return jsonify(error="PCA components are not available."), 400

    raw_palette_save = payload.get("palette")
    pair_palette_save = normalize_continuous_palette(
        str(raw_palette_save) if raw_palette_save is not None else None,
    )

    raw_name = str(payload.get("filename") or "zdim_pairplot.png").strip()
    filename = os.path.basename(raw_name) or "zdim_pairplot.png"
    if not filename.lower().endswith(".png"):
        filename += ".png"
    out_dir = os.path.join(e.workdir, f"analyze.{e.epoch}")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, filename)
    try:
        png, _ = pair_grid_png(
            e,
            lower_color_col=color_col,
            diagonal_emb=diagonal_emb,
            upper_style=upper_style,
            continuous_palette=pair_palette_save,
        )
        with open(out_path, "wb") as fh:
            fh.write(png)
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except OSError as err:
        logger.exception("save pairplot failed")
        return jsonify(error=str(err)), 500
    except Exception as err:
        logger.exception("save pairplot failed")
        return jsonify(error=str(err)), 500
    return jsonify(ok=True, path=out_path, filename=filename)


def api_default_trajectory_endpoints():
    """Return start/end in plot space along the long axis of the point cloud (2D PCA).

    The segment passes through the centroid and spans from the minimum to the maximum
    projection of points onto the first principal direction — a line through the
    middle of the mass with endpoints at opposite extents of the occupied region.
    """
    e: DashboardExperiment = g.dashboard_exp
    if not explorer_volumes_eligible(e):
        return (
            jsonify(
                error="Trajectory creator needs a CUDA GPU, single-particle data, "
                "and weights for the current epoch.",
            ),
            400,
        )
    xcol = request.args.get("x", "")
    ycol = request.args.get("y", "")
    if xcol not in e.plot_df.columns or ycol not in e.plot_df.columns:
        return jsonify(error="bad axis column"), 400

    sub = e.plot_df[[xcol, ycol]].dropna()
    if len(sub) < 2:
        return jsonify(error="not enough finite points for default trajectory"), 400
    xy = sub.values.astype(np.float64)
    mu = xy.mean(axis=0)
    xc = xy - mu
    span = xc.max(axis=0) - xc.min(axis=0)
    if not np.any(np.isfinite(span)) or float(np.nanmax(span)) < 1e-15:
        return jsonify(
            ok=True,
            start=mu.tolist(),
            end=(mu + np.array([1e-6, 0.0])).tolist(),
        )

    if float(span[0]) >= float(span[1]):
        v = np.array([1.0, 0.0])
    else:
        v = np.array([0.0, 1.0])

    try:
        _u, _s, vt = np.linalg.svd(xc, full_matrices=False)
        if vt.shape[0] >= 1:
            cand = vt[0].astype(np.float64)
            nrm = float(np.linalg.norm(cand))
            if nrm > 1e-15:
                v = cand / nrm
    except np.linalg.LinAlgError:
        pass

    t = xc @ v
    t_min = float(np.min(t))
    t_max = float(np.max(t))
    if np.isclose(t_min, t_max):
        bump = max(float(np.nanmax(span)), 1.0) * 0.05
        start = (mu + (t_min - bump) * v).tolist()
        end = (mu + (t_max + bump) * v).tolist()
    else:
        start = (mu + t_min * v).tolist()
        end = (mu + t_max * v).tolist()
    return jsonify(ok=True, start=start, end=end)


def trajectory_creator_page():
    e: DashboardExperiment = g.dashboard_exp
    if not explorer_volumes_eligible(e):
        return (
            render_template(
                "no_images.html",
                reason="Trajectory creator needs a CUDA GPU and model weights for the current epoch.",
            ),
            200,
        )
    cols = e.numeric_columns
    dx, dy = _trajectory_default_xy_cols(cols)
    return render_template(
        "trajectory_creator.html",
        numeric_cols=cols,
        covariate_display_map={c: _covariate_display_name(c) for c in cols},
        default_x=dx,
        default_y=dy,
        zdim=int(e.z.shape[1]),
    )


def api_trajectory_volumes():
    e: DashboardExperiment = g.dashboard_exp
    if not explorer_volumes_eligible(e):
        return (
            jsonify(
                error="Trajectory creator needs a CUDA GPU, single-particle data, "
                "and weights for the current epoch.",
            ),
            400,
        )
    data = request.get_json(force=True, silent=True) or {}
    mode = str(data.get("mode", "direct")).strip().lower()
    if mode not in ("direct", "nearest"):
        return jsonify(error='mode must be "direct" or "nearest".'), 400

    start_xy = data.get("start")
    end_xy = data.get("end")
    if (
        not isinstance(start_xy, list)
        or len(start_xy) != 2
        or not isinstance(end_xy, list)
        or len(end_xy) != 2
    ):
        return jsonify(error="start and end must be [x, y] coordinate pairs."), 400
    try:
        sx0, sy0 = float(start_xy[0]), float(start_xy[1])
        ex0, ey0 = float(end_xy[0]), float(end_xy[1])
    except (TypeError, ValueError):
        return jsonify(error="start/end coordinates must be numeric."), 400

    xcol = str(data.get("x", ""))
    ycol = str(data.get("y", ""))
    if xcol not in e.plot_df.columns or ycol not in e.plot_df.columns:
        return jsonify(error="bad axis column"), 400

    raw_n = data.get("n_points", 4)
    try:
        n_points = max(2, min(int(raw_n), 20))
    except (TypeError, ValueError):
        n_points = 4

    coords = e.plot_df[[xcol, ycol]].values.astype(np.float64)

    try:
        if mode == "direct":
            start_pt = np.array([sx0, sy0])
            end_pt = np.array([ex0, ey0])
            dists_start = np.sum((coords - start_pt) ** 2, axis=1)
            dists_end = np.sum((coords - end_pt) ** 2, axis=1)
            start_row = int(np.argmin(dists_start))
            end_row = int(np.argmin(dists_end))
            z_start = e.z[start_row]
            z_end = e.z[end_row]
            t = np.linspace(0.0, 1.0, n_points, dtype=np.float64)
            z_traj = np.outer(1.0 - t, z_start) + np.outer(t, z_end)
            traj_rows = None
            traj_xy = np.outer(1.0 - t, start_pt) + np.outer(t, end_pt)
            traj_xy = _round_direct_mode_traj_xy(traj_xy)
        else:
            t = np.linspace(0.0, 1.0, n_points, dtype=np.float64)
            line_xy = np.outer(1.0 - t, np.array([sx0, sy0])) + np.outer(
                t, np.array([ex0, ey0])
            )
            traj_rows = []
            for pt in line_xy:
                dists = np.sum((coords - pt) ** 2, axis=1)
                traj_rows.append(int(np.argmin(dists)))
            z_traj = e.z[np.array(traj_rows)]
            traj_xy = coords[np.array(traj_rows)]

        blobs, cache_token = generate_trajectory_volume_pngs(e, z_traj)
        b64s = [base64.standard_b64encode(b).decode("ascii") for b in blobs]
        payload: dict = {
            "ok": True,
            "images": b64s,
            "n_points": n_points,
            "mode": mode,
            "traj_rows": traj_rows,
            "traj_xy": traj_xy.tolist(),
            "xcol": xcol,
            "ycol": ycol,
            "volume_cache_id": cache_token,
        }
        if mode == "nearest" and traj_rows is not None:
            payload["traj_particle_indices"] = [
                _plot_row_particle_index(e, r) for r in traj_rows
            ]
            payload["particle_thumbs"] = [
                _particle_thumbnail_b64_from_row(e, int(r)) for r in traj_rows
            ]
        return jsonify(payload)
    except EnvironmentError as err:
        return jsonify(error=str(err), need_chimerax=True), 503
    except ValueError as err:
        return jsonify(error=str(err)), 400
    except Exception as err:
        logger.exception("trajectory volume generation failed")
        return jsonify(error=str(err)), 500


def create_app(
    workdir: str,
    epoch: int = -1,
    kmeans: int = -1,
    filter_plot_inds: str | None = None,
    cpus: int = 4,
) -> Flask:
    app = Flask(
        __name__,
        template_folder=_TEMPLATE_DIR,
        static_folder=_STATIC_DIR,
        static_url_path="/static",
    )
    app.jinja_env.filters["abbrev_middle"] = abbrev_middle
    app.config["PRELOAD_CPUS"] = max(1, cpus)
    app.secret_key = os.environ.get(
        "CRYODRGN_DASHBOARD_SECRET", "cryodrgn-dashboard-dev-key"
    )
    workdir = os.path.abspath(workdir)
    epochs = list_z_epochs(workdir)
    if not epochs:
        raise ValueError(
            f"No analyzed epochs under {workdir!r} — need z.N.pkl and analyze.N/ "
            "(run `cryodrgn analyze` first)."
        )
    start_epoch = epoch if epoch != -1 else max(epochs)
    if start_epoch not in epochs:
        start_epoch = max(epochs)
    app.config["DASHBOARD_WORKDIR"] = workdir
    app.config["DASHBOARD_KMEANS"] = kmeans
    app.config["DASHBOARD_EPOCHS"] = epochs
    app.config["DASHBOARD_START_EPOCH"] = start_epoch
    app.config["FILTER_PLOT_INDS"] = filter_plot_inds

    app.before_request(_bind_dashboard_exp)
    app.context_processor(inject_meta)
    app.add_url_rule("/api/set_epoch", view_func=api_set_epoch, methods=["POST"])
    app.add_url_rule("/", view_func=index, methods=["GET"])
    app.add_url_rule(
        "/command-builder", view_func=command_builder_page, methods=["GET"]
    )
    app.add_url_rule(
        "/abinit-builder", view_func=abinit_builder_redirect, methods=["GET"]
    )
    app.add_url_rule("/filter", view_func=filter_page_redirect, methods=["GET"])
    app.add_url_rule(
        "/api/save_selection", view_func=api_save_selection, methods=["POST"]
    )
    app.add_url_rule("/explorer", view_func=explorer, methods=["GET"])
    app.add_url_rule(
        "/api/explorer_volume_media",
        view_func=api_explorer_volume_media,
        methods=["POST"],
    )
    app.add_url_rule("/api/scatter", view_func=api_scatter, methods=["GET"])
    app.add_url_rule("/latent-3d", view_func=latent_3d_page, methods=["GET"])
    app.add_url_rule("/api/scatter3d_z", view_func=api_scatter3d_z, methods=["GET"])
    app.add_url_rule(
        "/api/preview_montage", view_func=api_preview_montage, methods=["GET"]
    )
    app.add_url_rule(
        "/api/preload_images", view_func=api_preload_images, methods=["GET", "POST"]
    )
    app.add_url_rule("/pairplot", view_func=pairplot_page, methods=["GET"])
    app.add_url_rule("/api/pairplot", view_func=api_pairplot, methods=["POST"])
    app.add_url_rule(
        "/api/save_pairplot_png", view_func=api_save_pairplot_png, methods=["POST"]
    )
    app.add_url_rule("/trajectory", view_func=trajectory_creator_page, methods=["GET"])
    app.add_url_rule(
        "/api/trajectory_volumes",
        view_func=api_trajectory_volumes,
        methods=["POST"],
    )
    app.add_url_rule(
        "/api/default_trajectory_endpoints",
        view_func=api_default_trajectory_endpoints,
        methods=["GET"],
    )

    return app


def run_server(
    workdir: str,
    epoch: int = -1,
    kmeans: int = -1,
    plot_inds: str | None = None,
    host: str = "127.0.0.1",
    port: int = 5050,
    debug: bool = False,
    cpus: int = 4,
) -> None:
    app = create_app(
        workdir=workdir,
        epoch=epoch,
        kmeans=kmeans,
        filter_plot_inds=plot_inds,
        cpus=cpus,
    )

    app.run(host=host, port=port, debug=debug, threaded=True)
