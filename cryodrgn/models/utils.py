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


def undeprecate_configs(cfg: dict[str, Any]) -> dict[str, Any]:
    """Replace parameter labels found from previous cryoDRGN versions."""

    model_args = cfg["model_args"] if "model_args" in cfg else cfg
    if "qlayers" in model_args:
        model_args["enc_layers"] = model_args["qlayers"]
        del model_args["qlayers"]
    if "players" in model_args:
        model_args["dec_layers"] = model_args["players"]
        del model_args["players"]
    if "qdim" in model_args:
        model_args["enc_dim"] = model_args["qdim"]
        del model_args["qdim"]
    if "pdim" in model_args:
        model_args["dec_dim"] = model_args["pdim"]
        del model_args["pdim"]
    if "domain" in model_args:
        model_args["volume_domain"] = model_args["domain"]
        del model_args["domain"]
    if "zdim" in model_args:
        model_args["z_dim"] = model_args["zdim"]
        del model_args["zdim"]
    if "ntilts" in model_args:
        model_args["n_tilts"] = model_args["ntilts"]
        del model_args["ntilts"]

    dataset_args = cfg["dataset_args"] if "dataset_args" in cfg else cfg
    if "keepreal" in dataset_args:
        dataset_args["use_real"] = dataset_args["keepreal"]
        del dataset_args["keepreal"]

    lattice_args = cfg["lattice_args"] if "lattice_args" in cfg else cfg
    if "extent" in lattice_args:
        lattice_args["l_extent"] = lattice_args["extent"]
        del lattice_args["extent"]

    return cfg


def get_model_trainer(
    cfg: dict[str, Any], outdir: str, add_cfgs: Optional[list[str]] = None
) -> ReconstructionModelTrainer:
    cfg = undeprecate_configs(cfg)
    model_args = cfg["model_args"] if "model_args" in cfg else cfg
    model = model_args["model"] if "model" in model_args else "cryodrgn"
    if "model" not in model_args:
        model_args["model"] = "cryodrgn"

    if model == "cryodrgn-ai":
        trainer_cls = AmortizedInferenceTrainer
    elif model == "cryodrgn":
        trainer_cls = HierarchicalPoseSearchTrainer
    else:
        raise ValueError(f"Unrecognized model `{model}` specified in config!")

    if add_cfgs:
        cfg.update(trainer_cls.config_cls.parse_cfg_keys(add_cfgs))

    return trainer_cls.load_from_config(cfg, outdir)


def get_model_configurations(
    cfg: dict[str, Any], add_cfgs: Optional[list[str]] = None
) -> ReconstructionModelConfigurations:
    cfg = undeprecate_configs(cfg)
    model_args = cfg["model_args"] if "model_args" in cfg else cfg

    if "model" not in model_args:
        configs_cls = AmortizedInferenceConfigurations
    elif model_args["model"] == "cryodrgn-ai":
        configs_cls = AmortizedInferenceConfigurations
    elif model_args["model"] == "cryodrgn":
        configs_cls = HierarchicalPoseSearchConfigurations
    else:
        raise ValueError(
            f"Model unknown by cryoDRGN `{model_args['model']}` specified in config!"
        )

    cfg = configs_cls.parse_config(cfg)
    if add_cfgs:
        cfg.update(configs_cls.parse_cfg_keys(add_cfgs))

    return configs_cls(**cfg)


# TODO: redundancy with `make_reconstruction_model` from trainers?
def get_model(
    cfg: dict[str, Any],
    outdir: str,
    add_cfgs: Optional[list[str]] = None,
    weights=None,
    device=None,
) -> Union[HetOnlyVAE, DRGNai]:
    configs = get_model_configurations(cfg, add_cfgs)
    lattice = Lattice(
        cfg["lattice_args"]["D"], extent=cfg["lattice_args"]["l_extent"], device=device
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

        if "n_particles" in cfg["dataset_args"]:
            particle_count = cfg["dataset_args"]["n_particles"]
            if "n_images" not in cfg["dataset_args"]:
                image_count = particle_count
            else:
                image_count = cfg["dataset_args"]["n_images"]
        else:
            trainer = get_model_trainer(cfg, outdir, add_cfgs)
            particle_count, image_count = trainer.particle_count, trainer.image_count

        if cfg["model_args"]["hypervolume_params"]["pe_dim"] is None:
            cfg["model_args"]["hypervolume_params"]["pe_dim"] = lattice.D // 2

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
            ps_params=cfg["model_args"]["ps_params"]
            if "ps_params" in cfg["model_args"]
            else None,
            verbose_time=configs.verbose_time,
            pretrain_with_gt_poses=configs.pretrain_with_gt_poses,
            n_tilts_pose_search=configs.n_tilts_pose_search,
        )

    elif isinstance(configs, HierarchicalPoseSearchConfigurations):
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

            if cfg["model_args"]["encode_mode"] == "tilt":
                tilt_params = dict(
                    t_emb_dim=cfg["model_args"]["t_emb_dim"],
                    tdim=cfg["model_args"]["tdim"],
                    tlayers=cfg["model_args"]["tlayers"],
                    ntilts=cfg["model_args"]["n_tilts"],
                )
            else:
                tilt_params = dict()

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
                tilt_params=tilt_params,
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
