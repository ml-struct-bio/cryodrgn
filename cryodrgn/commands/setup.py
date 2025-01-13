"""Utilities for creating experiment output folders and configuration files."""

import os
import argparse
import yaml
from typing import Optional, Union
from cryodrgn.trainers import ModelConfigurations
from cryodrgn.utils import load_yaml


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("config_file", help="experiment config file (.yaml)")

    parser.add_argument(
        "--model",
        "-m",
        default="amort",
        choices=["amort", "hps"],
        help="which generation of cryoDRGN learning models to apply",
    )

    parser.add_argument("--dataset", help="which dataset to run the experiment on")
    parser.add_argument(
        "--particles", help="path to the picked particles (.mrcs/.star /.txt)"
    )
    parser.add_argument("--ctf", help="path to the CTF parameters (.pkl)")
    parser.add_argument("--poses", help="path to the poses (.pkl)")

    parser.add_argument(
        "--capture-setup",
        default="spa",
        choices=["spa", "et"],
        help="`spa` for single-particle imaging (default) "
        "or `et` for electron tomography",
    )
    parser.add_argument(
        "--reconstruction-type",
        default="homo",
        choices=["het", "homo"],
        help="homogeneous (default) or heterogeneous reconstruction?",
    )
    parser.add_argument(
        "--pose-estimation",
        default="abinit",
        choices=["abinit", "refine", "fixed"],
        help="`abinit` for no initialization (default), `refine` to refine "
        "ground truth poses by gradient descent or `fixed` to use ground "
        "truth poses",
    )

    parser.add_argument(
        "--cfgs",
        "-c",
        nargs="+",
        help="additional configuration parameters to pass to the model "
        "in the form of 'CFG_KEY1=CFG_VAL1' 'CFG_KEY2=CFG_VAL2' ... ",
    )


class SetupHelper:
    def __init__(self, config_file: str, update_existing: bool = True) -> None:

        self.configs_file = config_file
        self.update_existing = update_existing
        config_dir = os.path.dirname(self.configs_file)
        if config_dir:
            os.makedirs(config_dir, exist_ok=True)

        if os.path.exists(self.configs_file):
            self.old_configs = load_yaml(self.configs_file)
        else:
            self.old_configs = dict()

    def create_configs(
        self,
        model: Optional[str] = None,
        dataset: Optional[str] = None,
        particles: Optional[str] = None,
        ctf: Optional[str] = None,
        poses: Optional[str] = None,
        additional_cfgs: Optional[list[str]] = None,
        **cfg_args,
    ) -> dict:
        configs: dict[str, Union[str, dict, None]] = self.old_configs.copy()

        if model:
            configs["model"] = model
        elif "model" not in configs:
            configs["model"] = "amort"

        if dataset:
            configs["dataset"] = dataset
        if particles:
            configs["particles"] = particles
        if ctf:
            configs["ctf"] = ctf
        if poses:
            configs["poses"] = poses
        if additional_cfgs:
            configs = {**configs, **ModelConfigurations.parse_cfg_keys(additional_cfgs)}

        configs = {**configs, **cfg_args}

        # turn anything that looks like a relative path into an absolute path
        for k in list(configs):
            if isinstance(configs[k], str) and not os.path.isabs(configs[k]):
                new_path = os.path.abspath(configs[k])
                if os.path.exists(new_path):
                    configs[k] = new_path

        paths_file = os.environ.get("DRGNAI_DATASETS")
        if paths_file:
            with open(paths_file, "r") as f:
                data_paths = yaml.safe_load(f)
        else:
            data_paths = None

        # handling different ways of specifying the input data, starting with a
        # file containing the data files
        if "dataset" in configs:
            if os.path.exists(configs["dataset"]):
                paths = load_yaml(configs["dataset"])

                # resolve paths relative to the dataset file if they look relative
                for k in list(paths):
                    if paths[k] and not os.path.isabs(paths[k]):
                        paths[k] = os.path.abspath(
                            os.path.join(configs["dataset"], paths[k])
                        )

                del configs["dataset"]

            elif data_paths and configs["dataset"] not in data_paths:
                raise ValueError(
                    f"Given dataset {configs['dataset']} is not a "
                    "label in the list of known datasets!"
                )

            elif data_paths is None:
                raise ValueError(
                    "To specify datasets using a label, first specify"
                    "a .yaml catalogue of datasets using the "
                    "environment variable $DRGNAI_DATASETS!"
                )

            # you can also give the dataset as a label in the global dataset list
            else:
                paths = data_paths[configs["dataset"]]

        # one can also specify the dataset files themselves in the config file
        elif "particles" in configs and "ctf" in configs:
            paths = {"particles": configs["particles"], "ctf": configs["ctf"]}

            if "pose" in configs and configs["pose"]:
                paths["pose"] = configs["pose"]

            if "dataset" in configs:
                del configs["dataset"]

        else:
            raise ValueError(
                "Must specify either a dataset label stored in "
                f"{paths_file} or the paths to a particles and "
                "ctf settings file!"
            )

        configs = {**configs, **paths}

        if self.update_existing:
            with open(self.configs_file, "w") as f:
                yaml.dump(configs, f, sort_keys=False)

        return configs


def main(args: argparse.Namespace) -> None:
    setup_helper = SetupHelper(args.config_file)
    setup_helper.create_configs(
        args.model,
        args.dataset,
        args.particles,
        args.ctf,
        args.poses,
        args.cfgs,
        capture_setup=args.capture_setup,
        reconstruction_type=args.reconstruction_type,
        pose_estimation=args.pose_estimation,
    )
