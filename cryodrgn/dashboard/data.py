"""Load heterogeneous-reconstruction outputs for dashboard UIs (mirrors `filter` command)."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Any, Optional

import numpy as np
import pandas as pd
import yaml
from scipy.spatial import transform

from cryodrgn import analysis, utils
from cryodrgn.dataset import ImageDataset, TiltSeriesData


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
    umap: Optional[np.ndarray]
    pc: np.ndarray
    z: np.ndarray

    @property
    def numeric_columns(self) -> list[str]:
        cols = [
            c
            for c in self.plot_df.select_dtypes(include=[np.number]).columns
            if c != "index"
        ]
        return cols

    @property
    def can_preview_particles(self) -> bool:
        return self.enc_mode != "tilt"


def load_experiment(
    workdir: str,
    epoch: int = -1,
    kmeans: int = -1,
) -> DashboardExperiment:
    train_configs_file = os.path.join(workdir, "config.yaml")
    if not os.path.exists(train_configs_file):
        raise ValueError("Missing config.yaml file in given output folder!")

    conf_fls = [fl for fl in os.listdir(workdir) if re.fullmatch(r"z\.[0-9]+\.pkl", fl)]
    if not conf_fls:
        raise NotImplementedError(
            "Dashboard filtering views need heterogeneous outputs (z.N.pkl)."
        )

    with open(train_configs_file, "r") as f:
        train_configs = yaml.safe_load(f)

    if epoch == -1:
        epoch = max(int(x.split(".")[1]) for x in conf_fls)

    anlzdir = os.path.join(workdir, f"analyze.{epoch}")
    if not os.path.isdir(anlzdir):
        raise ValueError(
            f"No analysis available for epoch {epoch} — "
            f"first run `cryodrgn analyze {workdir} {epoch}`"
        )

    z = utils.load_pkl(os.path.join(workdir, f"z.{epoch}.pkl"))

    if "poses" in train_configs["dataset_args"]:
        pose_pkl = train_configs["dataset_args"]["poses"]
    else:
        pose_pkl = os.path.join(workdir, f"pose.{epoch}.pkl")

    rot, trans = utils.load_pkl(pose_pkl)
    if train_configs["dataset_args"]["ctf"] is not None:
        ctf_params = utils.load_pkl(train_configs["dataset_args"]["ctf"])
    else:
        ctf_params = None

    if "encode_mode" in train_configs["model_args"]:
        enc_mode = train_configs["model_args"]["encode_mode"]
    else:
        enc_mode = "autodec"

    if isinstance(train_configs["dataset_args"]["ind"], int):
        indices = slice(train_configs["dataset_args"]["ind"])
    elif train_configs["dataset_args"]["ind"] is not None:
        indices = utils.load_pkl(train_configs["dataset_args"]["ind"])
    else:
        indices = None

    imgs_fl = train_configs["dataset_args"]["particles"]
    datadir = train_configs["dataset_args"].get("datadir") or ""

    if ctf_params is not None and enc_mode != "tilt":
        all_indices = np.array(range(ctf_params.shape[0]))
    elif enc_mode == "tilt":
        pt, tp = TiltSeriesData.parse_particle_tilt(imgs_fl)
        all_indices = np.array(range(len(pt)))
    else:
        all_indices = np.array(range(len(ImageDataset(mrcfile=imgs_fl, lazy=True))))

    if indices is not None:
        ctf_params = ctf_params[indices, :] if ctf_params is not None else None
        all_indices = all_indices[indices]
        if "poses" in train_configs["dataset_args"]:
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
        if len(kmeans_dirs) == 0:
            raise RuntimeError("No k-means outputs found under the analysis folder.")
        kmeans_dir = os.path.join(anlzdir, kmeans_dirs[0])
    else:
        kmeans_dir = os.path.join(anlzdir, f"kmeans{kmeans}")
        if not os.path.exists(kmeans_dir):
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
    ds = ImageDataset(mrcfile=exp.particles_path, lazy=True, datadir=exp.datadir)
    img = ds.src.images(g, as_numpy=True)
    if img.ndim == 3:
        img = img[0]
    return np.asarray(img, dtype=np.float32)
