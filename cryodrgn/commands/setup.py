"""Utilities for creating experiment output folders and configuration files."""

import os
import shutil
import argparse
import yaml
from typing import Optional, Union
from cryodrgn.utils import load_yaml


def add_args(parser):
    parser.add_argument("outdir", help="experiment output location")

    parser.add_argument(
        "--remove-existing", action="store_true", help="remove existing output folder"
    )

    parser.add_argument(
        "--model",
        "-m",
        default="v4",
        choices=["v2", "v3", "v4"],
        help="which generation of cryoDRGN learning models to apply",
    )

    parser.add_argument("--dataset", help="which dataset to run the experiment on")
    parser.add_argument(
        "--particles", help="path to the picked particles (.mrcs/.star /.txt)"
    )
    parser.add_argument("--ctf", help="path to the CTF parameters (.pkl)")
    parser.add_argument("--pose", help="path to the poses (.pkl)")

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
        "--conf-estimation",
        default="autodecoder",
        choices=["encoder", "autodecoder", "refine"],
        help="conformation estimation mode for heterogenous reconstruction: "
        "`autodecoder` (default), "
        "`encoder` or `refine` to refine conformations by gradient "
        "descent (you must then define initial_conf)",
    )


class SetupHelper:
    def __init__(
        self, outdir: str, remove_existing: bool = False, update_existing: bool = True
    ) -> None:

        self.outdir = outdir
        self.configs_file = os.path.join(outdir, "configs.yaml")
        self.update_existing = update_existing

        if remove_existing:
            shutil.rmtree(self.outdir)

        os.makedirs(self.outdir, exist_ok=True)
        if os.path.exists(self.configs_file):
            self.old_configs = load_yaml(self.configs_file)
        else:
            self.old_configs = dict()

        if "quick_config" not in self.old_configs:
            self.old_configs["quick_config"] = dict()

    def create_configs(
        self,
        model: Optional[str] = None,
        dataset: Optional[str] = None,
        particles: Optional[str] = None,
        ctf: Optional[str] = None,
        pose: Optional[str] = None,
        capture_setup: Optional[str] = None,
        reconstruction_type: Optional[str] = None,
        pose_estimation: Optional[str] = None,
        conf_estimation: Optional[str] = None,
    ) -> dict:
        configs: dict[str, Union[str, dict, None]] = self.old_configs.copy()

        if model:
            configs["model"] = model
        else:
            configs["model"] = "amort"

        if dataset:
            configs["dataset"] = dataset
        if particles:
            configs["particles"] = particles
        if ctf:
            configs["ctf"] = ctf
        if pose:
            configs["pose"] = pose

        if capture_setup:
            configs["quick_config"]["capture_setup"] = capture_setup
        if reconstruction_type:
            configs["quick_config"]["reconstruction_type"] = reconstruction_type
        if pose_estimation:
            configs["quick_config"]["pose_estimation"] = pose_estimation
        if conf_estimation:
            configs["quick_config"]["conf_estimation"] = conf_estimation

        # turn anything that looks like a relative path into an absolute path
        for k in list(configs):
            if isinstance(configs[k], str) and not os.path.isabs(configs[k]):
                new_path = os.path.abspath(os.path.join(self.outdir, configs[k]))

                if os.path.exists(new_path):
                    configs[k] = new_path

        if "reconstruction_type" in configs:
            if configs["quick_config"]["reconstruction_type"] == "homo":
                configs["quick_config"]["conf_estimation"] = None

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

        # parameters only the model code needs to use, not the user
        configs["outdir"] = os.path.join(self.outdir, "out")

        return configs


def main(args):
    setup_helper = SetupHelper(args.outdir, args.remove_existing)

    setup_helper.create_configs(
        args.model,
        args.dataset,
        args.particles,
        args.ctf,
        args.pose,
        args.capture_setup,
        args.reconstruction_type,
        args.pose_estimation,
        args.conf_estimation,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
