"""Utilities shared across all types of models."""

from typing import Any, Optional
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
