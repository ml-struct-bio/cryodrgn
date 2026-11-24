"""Load heterogeneous-reconstruction outputs for dashboard UIs.

Mirrors the ``filter`` command for the dashboard pages.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from functools import cached_property
from typing import Any

import numpy as np
import pandas as pd
import yaml
from scipy.spatial import transform

from cryodrgn import analysis, utils
from cryodrgn.dataset import ImageDataset, TiltSeriesData


def list_z_epochs(workdir: str) -> list[int]:
    """Sorted epochs with both ``z.{{epoch}}.pkl`` and ``analyze.{{epoch}}/``
    under ``workdir``."""
    if not os.path.isdir(workdir):
        return []
    out: list[int] = []
    with os.scandir(workdir) as it:
        for entry in it:
            if not entry.is_file():
                continue
            m = re.fullmatch(r"z\.([0-9]+)\.pkl", entry.name)
            if not m:
                continue
            ep = int(m.group(1))
            if os.path.isdir(os.path.join(workdir, f"analyze.{ep}")):
                out.append(ep)
    return sorted(out)


def discover_analyzed_workdirs(parent: str) -> list[str]:
    """Direct subfolders of ``parent`` with at least one analyzed epoch.

    Each candidate must have a ``z.N.pkl`` next to an ``analyze.N/`` folder.
    """
    if not os.path.isdir(parent):
        return []
    out: list[str] = []
    try:
        names = sorted(os.listdir(parent), key=str.lower)
    except OSError:
        return []
    for name in names:
        path = os.path.join(parent, name)
        if os.path.isdir(path) and list_z_epochs(path):
            out.append(os.path.abspath(path))
    return out


TILT_NO_THUMBNAILS_MSG = (
    "Tilt-series data does not include single-particle thumbnails in this view."
)
TILT_ABINIT_NO_IMAGES_MSG = (
    "This is a tilt-series ab initio run: each point is a particle made of several "
    "tilt images, and the dashboard does not yet map particles to their images, so "
    "particle thumbnails are turned off for it."
)
TILT_ABINIT_NO_EXPORT_MSG = (
    "Selections cannot be saved for tilt-series ab initio runs yet: their --ind "
    "files list tilt images, and the dashboard does not yet convert the selected "
    "particles into tilt images."
)
TILT_ABINIT_NO_PLOT_INDS_MSG = (
    "index files of tilt-series ab initio runs may list tilt images, and the "
    "dashboard does not yet convert tilt images into particles."
)


def is_tilt_series_abinit(train_configs: dict[str, Any]) -> bool:
    """True for a tilt-series (subtomogram averaging) ``cryodrgn abinit`` run.

    These runs do not set ``model_args.encode_mode``; they are marked by
    ``training.subtomogram_averaging``. Their ``z.N.pkl`` has one row per
    particle, while the particle stack, the CTF table and ``pose.N.pkl`` have
    one row per tilt image. Older configs store that tilt-image filter in
    ``ind``. Current abinit stores particle indices in ``ind`` and any
    tilt-image filter in ``tilt_ind``.
    """
    training = train_configs.get("training")
    return isinstance(training, dict) and bool(training.get("subtomogram_averaging"))


@dataclass
class DashboardExperiment:
    """Everything the dashboard needs for scatter / pair-grid / image preview."""

    workdir: str
    epoch: int
    kmeans_folder_id: int
    plot_df: pd.DataFrame
    all_indices: np.ndarray
    train_configs: dict[str, Any]
    particles_path: str
    datadir: str
    enc_mode: str
    kmeans_labels: np.ndarray
    umap: np.ndarray | None
    pc: np.ndarray
    z: np.ndarray
    # Tilt-series ab initio run: rows are particles, images are per tilt.
    tilt_abinit: bool = False
    _image_source: object | None = field(
        default=None, init=False, repr=False, compare=False
    )
    user_covariate_columns: list[str] = field(
        default_factory=list, init=False, repr=False, compare=False
    )

    @cached_property
    def numeric_columns(self) -> list[str]:
        """Numeric ``plot_df`` columns except ``index``.

        Cached per instance; ``plot_df`` is fixed after load.
        """
        return [
            c
            for c in self.plot_df.select_dtypes(include=[np.number]).columns
            if c != "index"
        ]

    @cached_property
    def color_covariate_columns(self) -> list[str]:
        """Columns available for scatter colour (numeric plus user-loaded labels)."""
        out = list(self.numeric_columns)
        for c in self.user_covariate_columns:
            if c in self.plot_df.columns and c not in out:
                out.append(c)
        return out

    @property
    def can_preview_particles(self) -> bool:
        return self.enc_mode != "tilt"

    @property
    def particle_images_unavailable(self) -> str | None:
        """Why this run has no image per plotted row (``None`` when it has one).

        For a tilt-series ab initio run, row ``i`` of the table is particle ``i``
        but image ``i`` of the stack is a tilt image of some other particle, so
        showing it would be wrong rather than merely unhelpful.
        """
        if self.enc_mode == "tilt":
            return TILT_NO_THUMBNAILS_MSG
        if self.tilt_abinit:
            return TILT_ABINIT_NO_IMAGES_MSG
        return None

    @property
    def selection_export_unavailable(self) -> str | None:
        """Why selected rows cannot be written as a dataset-index file (or ``None``)."""
        if self.tilt_abinit:
            return TILT_ABINIT_NO_EXPORT_MSG
        return None

    @property
    def analyze_dir(self) -> str:
        return os.path.join(self.workdir, f"analyze.{self.epoch}")

    @property
    def kmeans_dir(self) -> str:
        return os.path.join(self.analyze_dir, f"kmeans{self.kmeans_folder_id}")

    @property
    def weights_path(self) -> str:
        return os.path.join(self.workdir, f"weights.{self.epoch}.pkl")

    @property
    def z_pkl_path(self) -> str:
        return os.path.join(self.workdir, f"z.{self.epoch}.pkl")

    @property
    def landscape_dir(self) -> str:
        return os.path.join(self.workdir, f"landscape.{self.epoch}")

    @property
    def landscape_full_dir(self) -> str:
        return os.path.join(self.workdir, f"landscape_full.{self.epoch}")

    def particle_image_source(self):
        """Lazily opened particle stack shared by explorer thumbnail preload."""
        if self._image_source is None:
            from cryodrgn.source import ImageSource

            self._image_source = ImageSource.from_file(
                self.particles_path, lazy=True, datadir=self.datadir or ""
            )
        return self._image_source


def load_experiment(
    workdir: str,
    epoch: int = -1,
    kmeans: int = -1,
) -> DashboardExperiment:
    """Build a :class:`DashboardExperiment` from a training output folder.

    When ``epoch`` is ``-1``, uses the latest epoch that has both ``z.N.pkl``
    and ``analyze.N/``. When ``kmeans`` is ``-1``, uses the first ``kmeans*``
    directory found under that analysis folder.
    """
    train_configs_file = os.path.join(workdir, "config.yaml")
    if not os.path.isfile(train_configs_file):
        raise ValueError("Missing config.yaml file in given output folder!")

    has_z_pkl = any(re.fullmatch(r"z\.[0-9]+\.pkl", fl) for fl in os.listdir(workdir))
    if not has_z_pkl:
        raise NotImplementedError(
            "Dashboard scatter views need heterogeneous outputs (z.N.pkl)."
        )

    with open(train_configs_file) as f:
        train_configs = yaml.safe_load(f)

    if epoch == -1:
        analyzed = list_z_epochs(workdir)
        if not analyzed:
            raise ValueError(
                "No epochs with analysis outputs — run `cryodrgn analyze` for at "
                f"least one epoch (need {workdir}/analyze.N/ next to z.N.pkl)."
            )
        epoch = max(analyzed)

    anlzdir = os.path.join(workdir, f"analyze.{epoch}")
    if not os.path.isdir(anlzdir):
        raise ValueError(
            f"No analysis available for epoch {epoch} — "
            f"first run `cryodrgn analyze {workdir} {epoch}`"
        )

    z = utils.load_pkl(os.path.join(workdir, f"z.{epoch}.pkl"))

    ds_args = train_configs["dataset_args"]
    pose_pkl = ds_args.get("poses") or os.path.join(
        workdir,
        f"pose.{epoch}.pkl",
    )
    # ``abinit --no-trans`` saves the rotations alone instead of a (rot, trans) pair.
    poses = utils.load_pkl(pose_pkl)
    if isinstance(poses, (tuple, list)):
        rot = poses[0]
        trans = poses[1] if len(poses) > 1 else None
    else:
        rot, trans = poses, None

    ctf_path = ds_args["ctf"]
    ctf_params = utils.load_pkl(ctf_path) if ctf_path is not None else None

    enc_mode = train_configs["model_args"].get("encode_mode", "autodec")

    ind_raw = ds_args["ind"]
    if isinstance(ind_raw, int):
        indices = slice(ind_raw)
    elif ind_raw is not None:
        indices = utils.load_pkl(ind_raw)
    else:
        indices = None

    imgs_fl = ds_args["particles"]
    datadir = ds_args.get("datadir") or ""
    tilt_abinit = is_tilt_series_abinit(train_configs)

    if tilt_abinit:
        # Rows of z are particles. Older runs stored tilt-image indices in
        # ``ind``. Current abinit stores particle indices there and puts a
        # tilt-image list in ``tilt_ind`` when one was passed.
        all_indices = np.arange(z.shape[0])
    elif ctf_params is not None and enc_mode != "tilt":
        all_indices = np.arange(ctf_params.shape[0])
    elif enc_mode == "tilt":
        pt, _tp = TiltSeriesData.parse_particle_tilt(imgs_fl)
        all_indices = np.arange(len(pt))
    else:
        all_indices = np.arange(len(ImageDataset(mrcfile=imgs_fl, lazy=True)))

    # Per-tilt tables (CTF, poses) are sliced by tilt-image indices. A new
    # abinit config records those separately from particle ``ind``.
    table_indices = indices
    if tilt_abinit and "tilt_ind" in ds_args:
        tilt_raw = ds_args.get("tilt_ind")
        if isinstance(tilt_raw, int):
            table_indices = slice(tilt_raw)
        elif tilt_raw is not None:
            table_indices = utils.load_pkl(tilt_raw)
        else:
            table_indices = None
    if table_indices is not None:
        ctf_params = ctf_params[table_indices, :] if ctf_params is not None else None
        if not tilt_abinit:
            all_indices = all_indices[table_indices]
        if "poses" in ds_args:
            rot = rot[table_indices, :, :]
            trans = trans[table_indices, :] if trans is not None else None

    pc, _ = analysis.run_pca(z)
    umap_path = os.path.join(anlzdir, "umap.pkl")
    umap = utils.load_pkl(umap_path) if os.path.isfile(umap_path) else None

    if kmeans == -1:
        kmeans_dirs = [
            d
            for d in os.listdir(anlzdir)
            if os.path.isdir(os.path.join(anlzdir, d))
            and re.match(r"^kmeans[0-9]+$", d)
        ]
        if not kmeans_dirs:
            raise RuntimeError("No k-means outputs found under the analysis folder.")
        kmeans_dir = os.path.join(anlzdir, kmeans_dirs[0])
    else:
        kmeans_dir = os.path.join(anlzdir, f"kmeans{kmeans}")
        if not os.path.isdir(kmeans_dir):
            raise ValueError(f"No k-means results for k={kmeans}.")

    km_match = re.search(r"kmeans(\d+)$", os.path.basename(kmeans_dir))
    kmeans_folder_id = int(km_match.group(1)) if km_match else kmeans

    kmeans_lbls = utils.load_pkl(os.path.join(kmeans_dir, "labels.pkl"))
    znorm = np.sum(z**2, axis=1) ** 0.5

    if rot.shape[0] == z.shape[0]:
        plot_df = analysis.load_dataframe(
            z=z,
            pc=pc,
            euler=transform.Rotation.from_matrix(rot).as_euler("zyz", degrees=True),
            trans=trans,
            labels=kmeans_lbls,
            umap=umap,
            znorm=znorm,
        )
        if ctf_params is not None:
            plot_df["df1"] = ctf_params[:, 2]
            plot_df["df2"] = ctf_params[:, 3]
            plot_df["dfang"] = ctf_params[:, 4]
            plot_df["phase"] = ctf_params[:, 8]
    else:
        plot_df = analysis.load_dataframe(
            z=z, pc=pc, labels=kmeans_lbls, umap=umap, znorm=znorm
        )

    return DashboardExperiment(
        workdir=workdir,
        epoch=epoch,
        kmeans_folder_id=kmeans_folder_id,
        plot_df=plot_df,
        all_indices=all_indices,
        train_configs=train_configs,
        particles_path=imgs_fl,
        datadir=datadir,
        enc_mode=enc_mode,
        kmeans_labels=np.asarray(kmeans_lbls),
        umap=umap,
        pc=pc,
        z=z,
        tilt_abinit=tilt_abinit,
    )


def particle_image_array(exp: DashboardExperiment, row_index: int) -> np.ndarray:
    """Return a single 2D real-space image for dataframe row ``row_index``."""
    unavailable = exp.particle_images_unavailable
    if unavailable:
        raise RuntimeError(unavailable)
    g = int(exp.all_indices[int(row_index)])
    img = exp.particle_image_source().images(g, as_numpy=True)
    if img.ndim == 3:
        img = img[0]
    return np.asarray(img, dtype=np.float32)
