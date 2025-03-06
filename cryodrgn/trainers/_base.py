"""Base classes for model training engines and their parameter configurations."""

import argparse
import sys
import difflib
import pandas as pd
from abc import ABC
from dataclasses import dataclass, fields, Field, MISSING, asdict
from typing import Any
from typing_extensions import Self
import yaml
import logging
import numpy as np


@dataclass
class BaseConfigurations(ABC):
    """The abstract base data class for config parameter sets used by cryoDRGN engines.

    This class' variables constitute the core parameters used by all cryoDRGN
    configuration parameter sets. These configurations are used by various engines
    that train cryoDRGN models; the base class for these engines is defined below as
    `BaseTrainer`. Python data classes inherit their parents' data fields, so all
    of this abstract class' children configuration classes contain these parameters
    in addition to the ones they themselves define.

    Note that unlike regular data classes these config classes must define defaults for
    all their parameters to ensure that default engine behaviour is explicitly stated,
    with an AssertionError being thrown upon initialization otherwise.

    Arguments
    ---------
    verbose:        An integer specifiying the verbosity level for this engine, with
                    the default value of 0 generally specifying no/minimum verbosity.
    seed:           A non-negative integer used to fix the stochasticity of the random
                    number generators used by this engine for reproducibility.
                    The default is to not fix stochasticity and thus use a different
                    random seed upon each run of the engine.
    test_installation:  Only perform a smoke test that this module has been installed
                        correctly and exit immediately without running anything if this
                        boolean value is set to `True`.
                        Default is not to run this test.
    """

    # A parameter belongs to this configuration set if and only if it has a type and a
    # default value defined here, note that children classes inherit these parameters
    verbose: bool = False
    seed: int = None
    test_installation: bool = False

    def __post_init__(self) -> None:
        """Parsing given configuration parameter values and checking their validity."""
        for this_field in fields(self):
            assert this_field.default is not MISSING, (
                f"`{self.__class__.__name__}` class has no default value defined "
                f"for parameter `{this_field.name}`!"
            )

        if self.test_installation:
            print("Installation was successful!")
            sys.exit()
        elif self.test_installation is not False:
            raise ValueError(
                f"Given `test_installation` value `{self.test_installation}` "
                f"cannot be interpreted as a boolean!"
            )

        if not isinstance(self.verbose, int) or self.verbose < 0:
            raise ValueError(
                f"Given verbosity `{self.verbose}` is not a positive integer!"
            )

        if self.seed is None:
            self.seed = np.random.randint(0, 10000)
        if not isinstance(self.seed, int):
            raise ValueError(
                "Configuration `seed` must be given as an integer, "
                f"given `{self.seed}` instead!"
            )

    def __iter__(self):
        return iter(asdict(self).items())

    def __str__(self):
        return "\n".join([f"{par}{str(val):>20}" for par, val in self])

    def __contains__(self, val) -> bool:
        return val in {k for k, _ in self}

    @property
    def file_dict(self) -> dict[str, Any]:
        return {"seed": self.seed}

    def write(self, fl: str) -> None:
        """Saving configurations to file using the original order."""

        with open(fl, "w") as f:
            yaml.dump(self.file_dict, f, default_flow_style=False, sort_keys=False)

    @classmethod
    def fields_dict(cls) -> dict[str, Field]:
        """Convenience method to get all fields indexed by field label."""
        return {fld.name: fld for fld in fields(cls)}

    @classmethod
    def parse_config(cls, configs: dict[str, Any]) -> dict[str, Any]:
        """Retrieves all configurations that have been saved to file."""
        cfg = {
            k.split("___")[-1]: v[0]
            for k, v in pd.json_normalize(configs, sep="___").items()
        }
        cfg = {k: v for k, v in cfg.items() if k in cls.fields_dict()}
        cfg = {
            k: cls.fields_dict()[k].type(v) if v is not None else None
            for k, v in cfg.items()
        }

        return cfg

    @classmethod
    def parse_cfg_keys(cls, cfg_keys: list[str]) -> dict[str, Any]:
        """Retrieve the parameter values given in a list of --cfgs command line entries.

        This method parses the parameters given by a user via a `--cfgs` flag defined
        for commands such as `drgnai setup` to provide an arbitrary set of
        configuration parameters through the command line interface.

        """
        cfgs = dict()

        for cfg_str in cfg_keys:
            if cfg_str.count("=") != 1:
                raise ValueError(
                    "--cfgs entries must have exactly one equals sign "
                    "and be in the form 'CFG_KEY=CFG_VAL'!"
                )
            cfg_key, cfg_val = cfg_str.split("=")

            if cfg_val is None or cfg_val == "None":
                cfgs[cfg_key] = None

            else:
                for fld in fields(cls):
                    if cfg_key == fld.name:
                        if fld.type is str:
                            cfgs[cfg_key] = str(cfg_val)
                        else:
                            cfgs[cfg_key] = fld.type(eval(cfg_val))

                        # accounting for parameters like `ind` which can be paths
                        # to files as well as integers
                        if isinstance(cfgs[cfg_key], str) and cfgs[cfg_key].isnumeric():
                            cfgs[cfg_key] = int(cfgs[cfg_key])

                        break

                else:
                    close_keys = difflib.get_close_matches(
                        cfg_key, [fld.name for fld in fields(cls)]
                    )

                    if close_keys:
                        close_str = f"\nDid you mean one of:\n{', '.join(close_keys)}"
                    else:
                        close_str = ""

                    raise ValueError(
                        f"--cfgs parameter `{cfg_key}` is not a "
                        f"valid configuration parameter!{close_str}"
                    )

        return cfgs


class BaseTrainer(ABC):
    """Abstract base class for training engines used by cryoDRGN.

    Arguments
    ---------
    configs (dict)      The raw configuration parameters for this engine.
                        Will be parsed using the engine's configuration class

    Attributes
    ----------
    config_cls          The configuration class that will be used to parse parameter
                        sets for this engine class. Children implementations of this
                        class will thus use children of `BaseConfigurations` here.
    outdir:             Path to where output produced by the engine will be saved.
    label:              String used to refer to this engine for e.g. logging messages.

    configs (BaseConfigurations)    The parsed parameter configs for this engine.
    logger (Logger)     Logging utility used to create info and warning messages.
    """

    config_cls = BaseConfigurations
    label = "cDRGN training"

    def __init__(self, configs: dict[str, Any], outdir: str) -> None:
        self.configs = self.config_cls(**configs)
        self.outdir = outdir
        np.random.seed(self.configs.seed)
        self.logger = logging.getLogger(self.label)

    @classmethod
    def defaults(cls) -> dict[str, Any]:
        """The user-set parameters governing the behaviour of this model."""
        return {fld.name: fld.default for fld in cls.config_cls.fields()}

    @classmethod
    def parameters(cls) -> list[str]:
        """The user-set parameters governing the behaviour of this model."""
        return [fld.name for fld in cls.config_cls.fields()]

    @classmethod
    def load_from_config(cls, configs: dict[str, Any], outdir: str) -> Self:
        """Retrieves all configurations that have been saved to file."""

        return cls(cls.config_cls.parse_config(configs), outdir)

    @classmethod
    def parse_args(cls, args: argparse.Namespace, outdir: str) -> Self:
        """Utility for initializing using a namespace as opposed to a dictionary."""

        return cls(
            {
                par: getattr(args, par) if hasattr(args, par) else cls.defaults()[par]
                for par in cls.parameters()
            },
            outdir,
        )
