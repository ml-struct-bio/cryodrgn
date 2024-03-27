"""Utilities shared across all types of models."""

import logging
from typing import Union, Optional
import torch
from torch import nn
from cryodrgn.lattice import Lattice
from cryodrgn.models.variational_autoencoder import HetOnlyVAE
from cryodrgn.models.amortized_inference import DRGNai, HyperVolume

logger = logging.getLogger(__name__)


def load_model(
    cfg: dict[str, dict], weights=Union[str, None], device=Optional[str]
) -> tuple[Union[HetOnlyVAE, DRGNai], Lattice, int]:
    """Instantiate a volume model from a config.yaml

    Inputs:
        config (str, dict): Path to config.yaml or loaded config.yaml
        weights (str): Path to weights.pkl
        device: torch.device object

    Returns:
        nn.Module instance, Lattice instance
    """

    # cryodrgn v3 and v4
    if "model_args" in cfg and cfg["model_args"]["model"] == "hps":
        lattice_args = cfg["lattice_args"]
        model_args = cfg["model_args"]

        lattice = Lattice(
            lattice_args["D"], extent=lattice_args["extent"], device=device
        )

        logger.info("loading a cryoDRGN v<=3 model...")

        # TODO: merge with Trainer.mask_dimensions?
        if model_args["enc_mask"] > 0:
            enc_mask = lattice.get_circular_mask(model_args["enc_mask"])
            in_dim = int(enc_mask.sum())
        else:
            assert model_args["enc_mask"] == -1
            enc_mask = None
            in_dim = lattice.D**2

        # TODO: pull activation from a global dictionary?
        model = HetOnlyVAE(
            lattice,
            model_args["qlayers"],
            model_args["qdim"],
            model_args["players"],
            model_args["pdim"],
            in_dim,
            model_args["z_dim"],
            encode_mode=model_args["encode_mode"],
            enc_mask=enc_mask,
            enc_type=model_args["pe_type"],
            enc_dim=model_args["pe_dim"],
            domain=model_args["domain"],
            activation={"relu": nn.ReLU, "leaky_relu": nn.LeakyReLU}[
                model_args["activation"]
            ],
            feat_sigma=model_args["feat_sigma"],
            tilt_params=model_args.get("tilt_params", {}),
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

        lattice = Lattice(
            checkpoint["hypervolume_params"]["resolution"], extent=0.5, device=device
        )

        radius_mask = (
            checkpoint["output_mask_radius"]
            if "output_mask_radius" in checkpoint
            else None
        )

    return model, lattice, radius_mask
