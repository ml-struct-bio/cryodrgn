"""Utilities shared across all types of models."""

import os
import logging
import yaml
from typing import Union, Optional
import torch
from torch import nn
from cryodrgn.lattice import Lattice
from cryodrgn.models.variational_autoencoder import HetOnlyVAE
from cryodrgn.models.amortized_inference import HyperVolume
from cryodrgn.trainers.hps_trainer import HierarchicalPoseSearchTrainer
from cryodrgn.trainers.amortinf_trainer import AmortizedInferenceTrainer

logger = logging.getLogger(__name__)
trainer_classes = {
    "hps": HierarchicalPoseSearchTrainer,
    "amort": AmortizedInferenceTrainer,
}


def load_from_config(
    config_file: str,
) -> Union[HierarchicalPoseSearchTrainer, AmortizedInferenceTrainer]:
    """Retrieves all configurations that have been saved to file."""
    if os.path.exists(config_file):
        with open(config_file, "r") as f:
            configs = yaml.safe_load(f)
    else:
        raise FileNotFoundError(
            f"Cannot find training configurations file `{config_file}` "
            f"— has this model been trained yet?"
        )

    cfg_dict = {
        sub_k: sub_v
        for k, v in configs.items()
        if isinstance(v, dict)
        for sub_k, sub_v in v.items()
    }
    cfg_dict.update({k: v for k, v in configs.items() if not isinstance(v, dict)})

    trainer_class = trainer_classes[cfg_dict["model"]]
    cfg_dict = {
        k: v for k, v in cfg_dict.items() if k in set(trainer_class.parameters())
    }

    return trainer_class(cfg_dict)


def load_model(
    cfg: dict[str, dict], weights=Union[str, None], device=Optional[str]
) -> tuple[torch.nn.Module, Lattice, int]:
    """Instantiate a volume model from a config.yaml

    Inputs:
        config (str, dict): Path to config.yaml or loaded config.yaml
        weights (str): Path to weights.pkl
        device: torch.device object

    Returns:
        nn.Module instance, Lattice instance
    """

    # cryodrgn v3 and v4
    if "dataset_args" in cfg:
        logger.info("loading a cryoDRGN v<=3 model...")

        c = cfg["lattice_args"]
        lat = Lattice(c["D"], extent=c["extent"], device=device)
        c = cfg["model_args"]

        if c["enc_mask"] > 0:
            enc_mask = lat.get_circular_mask(c["enc_mask"])
            in_dim = int(enc_mask.sum())
        else:
            assert c["enc_mask"] == -1
            enc_mask = None
            in_dim = lat.D**2

        activation = {"relu": nn.ReLU, "leaky_relu": nn.LeakyReLU}[c["activation"]]
        model = HetOnlyVAE(
            lat,
            c["qlayers"],
            c["qdim"],
            c["players"],
            c["pdim"],
            in_dim,
            c["z_dim"],
            encode_mode=c["encode_mode"],
            enc_mask=enc_mask,
            enc_type=c["pe_type"],
            enc_dim=c["pe_dim"],
            domain=c["domain"],
            activation=activation,
            feat_sigma=c["feat_sigma"],
            tilt_params=c.get("tilt_params", {}),
        )

        if weights is not None:
            ckpt = torch.load(weights, map_location=device)
            model.load_state_dict(ckpt["model_state_dict"])
        if device is not None:
            model.to(device)

        radius_mask = None

    else:
        logger.info("loading a DRGNai model...")

        checkpoint = torch.load(weights, map_location=device)
        hypervolume_params = checkpoint["hypervolume_params"]
        model = HyperVolume(**hypervolume_params)
        model.load_state_dict(checkpoint["hypervolume_state_dict"])
        model.to(device)

        lat = Lattice(
            checkpoint["hypervolume_params"]["resolution"], extent=0.5, device=device
        )

        radius_mask = (
            checkpoint["output_mask_radius"]
            if "output_mask_radius" in checkpoint
            else None
        )

    return model, lat, radius_mask
