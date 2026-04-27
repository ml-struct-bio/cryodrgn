"""Load heterogeneous-reconstruction outputs for dashboard UIs (mirrors `filter` command)."""

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
    """Return sorted epochs with both ``z.{{epoch}}.pkl`` and ``analyze.{{epoch}}/`` under ``workdir``."""
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
    _image_dataset: ImageDataset | None = field(
        default=None, init=False, repr=False, compare=False
    )

    @cached_property
    def numeric_columns(self) -> list[str]:
        """Numeric ``plot_df`` columns except ``index`` (cached per instance; ``plot_df`` is fixed after load)."""
        return [
            c
            for c in self.plot_df.select_dtypes(include=[np.number]).columns
            if c != "index"
        ]

    @property
    def can_preview_particles(self) -> bool:
        return self.enc_mode != "tilt"


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
                "No epochs with analysis outputs — run `cryodrgn analyze` for at least one epoch "
                f"(need {workdir}/analyze.N/ next to z.N.pkl)."
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
    rot, trans = utils.load_pkl(pose_pkl)

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

    if ctf_params is not None and enc_mode != "tilt":
        all_indices = np.arange(ctf_params.shape[0])
    elif enc_mode == "tilt":
        pt, _tp = TiltSeriesData.parse_particle_tilt(imgs_fl)
        all_indices = np.arange(len(pt))
    else:
        all_indices = np.arange(len(ImageDataset(mrcfile=imgs_fl, lazy=True)))

    if indices is not None:
        ctf_params = ctf_params[indices, :] if ctf_params is not None else None
        all_indices = all_indices[indices]
        if "poses" in ds_args:
            rot = rot[indices, :, :]
            trans = trans[indices, :]

    pc, _ = analysis.run_pca(z)
    umap = utils.load_pkl(os.path.join(anlzdir, "umap.pkl"))

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
    )


def particle_image_array(exp: DashboardExperiment, row_index: int) -> np.ndarray:
    """Return a single 2D real-space image for dataframe row ``row_index``."""
    if not exp.can_preview_particles:
        raise RuntimeError("Particle previews are not supported for tilt-series data.")
    g = int(exp.all_indices[int(row_index)])
    if exp._image_dataset is None:
        exp._image_dataset = ImageDataset(
            mrcfile=exp.particles_path, lazy=True, datadir=exp.datadir
        )
    img = exp._image_dataset.src.images(g, as_numpy=True)
    if img.ndim == 3:
        img = img[0]
    return np.asarray(img, dtype=np.float32)
