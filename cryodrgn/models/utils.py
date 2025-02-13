"""Utilities shared across all types of models."""

from typing import Any
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


def get_model_trainer(cfg: dict[str, Any]) -> ReconstructionModelTrainer:
    if "model" in cfg:
        model = cfg["model"]
    elif "model_args" in cfg and "model" in cfg["model_args"]:
        model = cfg["model_args"]["model"]
    else:
        model = "hps"

    if model == "amort":
        trainer_cls = AmortizedInferenceTrainer
    elif model == "hps":
        trainer_cls = HierarchicalPoseSearchTrainer
    else:
        raise ValueError(f"Unrecognized model `{model}` specified in config!")

    return trainer_cls.load_from_config(cfg)


def get_model_configurations(cfg: dict[str, Any]) -> ReconstructionModelConfigurations:
    cfg = ReconstructionModelConfigurations.parse_config(cfg)

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
    cfg = {
        k: configs_cls.fields_dict()[k].type(v) if v is not None else None
        for k, v in cfg.items()
    }

    return configs_cls(**cfg)
