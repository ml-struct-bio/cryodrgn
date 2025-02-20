"""Utilities shared across all types of models."""

from typing import Any, Optional, Union
import torch
from torch import nn
from cryodrgn.trainers.reconstruction import (
    ReconstructionModelTrainer,
    ReconstructionModelConfigurations,
)
from cryodrgn.trainers.amortinf_trainer import (
    AmortizedInferenceTrainer,
    AmortizedInferenceConfigurations,
)
from cryodrgn.trainers.hps_trainer import (
    HierarchicalPoseSearchTrainer,
    HierarchicalPoseSearchConfigurations,
)
from cryodrgn.models.amortized_inference import DRGNai
from cryodrgn.models.variational_autoencoder import HetOnlyVAE
from cryodrgn.models.neural_nets import get_decoder
from cryodrgn.lattice import Lattice
from cryodrgn.masking import CircularMask, FrequencyMarchingMask


def update_configs(cfg: dict[str, Any]) -> dict[str, Any]:
    if "model_args" in cfg:
        if "qlayers" in cfg["model_args"]:
            cfg["model_args"]["enc_layers"] = cfg["model_args"]["qlayers"]
            del cfg["model_args"]["qlayers"]
        if "players" in cfg["model_args"]:
            cfg["model_args"]["dec_layers"] = cfg["model_args"]["players"]
            del cfg["model_args"]["players"]
        if "qdim" in cfg["model_args"]:
            cfg["model_args"]["enc_dim"] = cfg["model_args"]["qdim"]
            del cfg["model_args"]["qdim"]
        if "pdim" in cfg["model_args"]:
            cfg["model_args"]["dec_dim"] = cfg["model_args"]["pdim"]
            del cfg["model_args"]["pdim"]
        if "domain" in cfg["model_args"]:
            cfg["model_args"]["volume_domain"] = cfg["model_args"]["domain"]
            del cfg["model_args"]["domain"]

    return cfg


def get_model_trainer(
    cfg: dict[str, Any], add_cfgs: Optional[list[str]] = None
) -> ReconstructionModelTrainer:
    cfg = update_configs(cfg)

    if "model" in cfg:
        model = cfg["model"]
    elif "model_args" in cfg and "model" in cfg["model_args"]:
        model = cfg["model_args"]["model"]
    else:
        model = "hps"
        cfg["model"] = "hps"

    if model == "amort":
        trainer_cls = AmortizedInferenceTrainer
    elif model == "hps":
        trainer_cls = HierarchicalPoseSearchTrainer
    else:
        raise ValueError(f"Unrecognized model `{model}` specified in config!")

    if add_cfgs:
        cfg.update(trainer_cls.config_cls.parse_cfg_keys(add_cfgs))

    return trainer_cls.load_from_config(cfg)


def get_model_configurations(
    cfg: dict[str, Any], add_cfgs: Optional[list[str]] = None
) -> ReconstructionModelConfigurations:
    cfg = ReconstructionModelConfigurations.parse_config(update_configs(cfg))

    if "model" not in cfg:
        configs_cls = HierarchicalPoseSearchConfigurations
    elif cfg["model"] == "amort":
        configs_cls = AmortizedInferenceConfigurations
    elif cfg["model"] == "hps":
        configs_cls = HierarchicalPoseSearchConfigurations
    else:
        raise ValueError(
            f"Model unrecognized by cryoDRGN: `{cfg['model']}` specified in config!"
        )

    if add_cfgs:
        cfg.update(configs_cls.parse_cfg_keys(add_cfgs))

    return configs_cls(**cfg)


def get_model(
    cfg: dict[str, Any],
    add_cfgs: Optional[list[str]] = None,
    weights=None,
    device=None,
) -> Union[HetOnlyVAE, DRGNai]:
    configs = get_model_configurations(cfg, add_cfgs)
    lattice = Lattice(
        cfg["lattice_args"]["D"], extent=cfg["lattice_args"]["extent"], device=device
    )

    if isinstance(configs, AmortizedInferenceConfigurations):
        if configs.output_mask == "circ":
            radius = configs.max_freq or lattice.D // 2
            output_mask = CircularMask(lattice, radius)

        elif configs.output_mask == "frequency_marching":
            output_mask = FrequencyMarchingMask(
                lattice,
                radius=configs.l_start_fm,
                radius_max=lattice.D // 2,
                add_one_every=configs.add_one_frequency_every,
            )
        else:
            raise NotImplementedError

        if "particle_count" in cfg["dataset_args"]:
            particle_count = cfg["dataset_args"]["particle_count"]
            if "image_count" not in cfg["dataset_args"]:
                image_count = particle_count
            else:
                image_count = cfg["dataset_args"]["image_count"]
        else:
            trainer = get_model_trainer(cfg, add_cfgs)
            particle_count, image_count = trainer.particle_count, trainer.image_count

        model = DRGNai(
            lattice=lattice,
            output_mask=output_mask,
            n_particles_dataset=particle_count,
            n_tilts_dataset=image_count,
            cnn_params=cfg["model_args"]["cnn_params"],
            conf_regressor_params=cfg["model_args"]["conf_regressor_params"],
            hypervolume_params=cfg["model_args"]["hypervolume_params"],
            resolution_encoder=configs.resolution_encoder,
            no_trans=configs.no_trans,
            use_gt_poses=configs.pose_estimation == "fixed",
            use_gt_trans=configs.use_gt_trans,
            will_use_point_estimates=False,
            ps_params=cfg["model_args"]["ps_params"],
            verbose_time=configs.verbose_time,
            pretrain_with_gt_poses=configs.pretrain_with_gt_poses,
            n_tilts_pose_search=configs.n_tilts_pose_search,
        )

    elif configs.model == "hps":
        activation = {"relu": nn.ReLU, "leaky_relu": nn.LeakyReLU}[configs.activation]
        if configs.z_dim > 0:
            if (
                cfg["model_args"]["enc_mask"] is not None
                and cfg["model_args"]["enc_mask"] > 0
            ):
                enc_mask = lattice.get_circular_mask(cfg["model_args"]["enc_mask"])
                in_dim = int(enc_mask.sum())
            else:
                enc_mask = None
                in_dim = lattice.D**2

            model = HetOnlyVAE(
                lattice=lattice,
                qlayers=cfg["model_args"]["enc_layers"],
                qdim=cfg["model_args"]["enc_dim"],
                players=cfg["model_args"]["dec_layers"],
                pdim=cfg["model_args"]["dec_dim"],
                in_dim=in_dim,
                z_dim=configs.z_dim,
                encode_mode=cfg["model_args"]["encode_mode"],
                enc_mask=enc_mask,
                enc_type=configs.pe_type,
                enc_dim=configs.pe_dim,
                domain=configs.volume_domain,
                activation=activation,
                feat_sigma=cfg["model_args"]["feat_sigma"],
                tilt_params=cfg["model_args"].get("tilt_params", {}),
            )
        else:
            model = get_decoder(
                in_dim=3,
                D=lattice.D,
                layers=cfg["model_args"]["qlayers"],
                dim=cfg["model_args"]["qlayers"],
                domain=configs.volume_domain,
                enc_type=configs.pe_type,
                enc_dim=configs.pe_dim,
                activation=activation,
                feat_sigma=configs.feat_sigma,
            )
        if weights is not None:
            ckpt = torch.load(weights, device=device)
            model.load_state_dict(ckpt["model_state_dict"])
        if device is not None:
            model.to(device)
    else:
        raise ValueError(f"Unrecognized model `{configs.model}` specified in config!")

    return model
