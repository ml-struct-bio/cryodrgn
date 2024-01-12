"""Utilities shared across all types of models."""

import logging
from typing import Union, Optional
import torch
from torch import nn
import cryodrgn.config
from cryodrgn.lattice import Lattice
from cryodrgn.models.variational_autoencoder import HetOnlyVAE
from cryodrgn.models.amortized_inference import HyperVolume

logger = logging.getLogger(__name__)


def load_model(
    config: Union[str, dict], weights=Union[str, None], device=Optional[str]
) -> tuple[torch.nn.Module, Lattice, int]:
    """Instantiate a volume model from a config.yaml

    Inputs:
        config (str, dict): Path to config.yaml or loaded config.yaml
        weights (str): Path to weights.pkl
        device: torch.device object

    Returns:
        nn.Module instance, Lattice instance
    """
    cfg = cryodrgn.config.load(config)

    if "cmd" in cfg:
        logger.info("loading a cryoDRGN2 model...")

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
            c["zdim"],
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
