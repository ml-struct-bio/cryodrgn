"""Base classes for model training engines and their parameter configurations."""

import os
import argparse
import sys
import difflib
import inspect
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

    This class also defines special behaviour for the `quick_config` class variable,
    which is not treated as a data field and instead defines a set of shortcuts used as
    values for the data field parameters listed as its keys. These shortcuts each define
    a list of fields and values that are used as the new defaults when the shortcut is
    used, but can still be overridden by values specified by the user.

    Note that unlike regular data classes these config classes must define defaults for
    all their parameters to ensure that default engine behaviour is explicitly stated,
    with an AssertionError being thrown upon initialization otherwise.

    Arguments
    ---------
    verbose:        An integer specifiying the verbosity level for this engine, with
                    the default value of 0 generally specifying no/minimum verbosity.

    outdir:         Path to where output produced by the engine will be saved.

    seed:           A non-negative integer used to fix the stochasticity of the random
                    number generators used by this engine for reproducibility.
                    The default is to not fix stochasticity and thus use a different
                    random seed upon each run of the engine.

    test_installation:  Only perform a smoke test that this module has been installed
                        correctly and exit immediately without running anything if this
                        boolean value is set to `True`.
                        Default is not to run this test.

    Attributes
    ----------
    quick_config:   A dictionary with keys consisting of special `quick_config` shortcut
                    parameters; each value is a dictionary of non-quick_config
                    parameter keys and shortcut values that are used when the
                    corresponding quick configuration parameter value is used.
    """

    # This class variable is not a dataclass field and is instead used to define shortcut
    # labels to set values for a number of other fields
    quick_config = dict()

    # A parameter belongs to this configuration set if and only if it has a type and a
    # default value defined here, note that children classes inherit these parameters
    verbose: int = 0
    outdir: str = os.getcwd()
    seed: int = None
    test_installation: bool = False

    def __init__(self, **config_args: dict[str, Any]) -> None:
        """Setting given config values as attributes; saving values given by user."""
        self.given_configs = config_args

        # for configuration values not given by the user we use the defined defaults
        for this_field in self.fields():
            assert this_field.default is not MISSING, (
                f"`{self.__class__.__name__}` class has no default value defined "
                f"for parameter `{this_field.name}`!"
            )

        # set values specified explicitly by the user as attributes of this class
        for k, v in self.given_configs.items():
            setattr(self, k, v)

        self.__post_init__()

    def __post_init__(self) -> None:
        """Parsing given configuration parameter values and checking their validity."""
        self.outdir = os.path.abspath(self.outdir)

        for quick_cfg_k, quick_cfg_dict in self.quick_config.items():
            assert quick_cfg_k in self, (
                f"Configuration class `{self.__class__.__name__}` has a `quick_config` "
                f"entry `{quick_cfg_k}` that is not a valid configuration parameter!"
            )
            for quick_cfg_label, quick_label_dict in quick_cfg_dict.items():
                for quick_cfg_param, quick_cfg_val in quick_label_dict.items():
                    assert quick_cfg_param in self, (
                        f"Configuration class `{self.__class__.__name__}` has a "
                        f"`quick_config` entry `{quick_cfg_label}` under "
                        f"`{quick_cfg_k}` with a value for `{quick_cfg_param}` which "
                        f"is not a valid configuration parameter!"
                    )

            if quick_cfg_k in self.given_configs:
                quick_cfg_val = getattr(self, quick_cfg_k)
                if quick_cfg_val is not None:
                    if quick_cfg_val not in self.quick_config[quick_cfg_k]:
                        raise ValueError(
                            f"Given value `{quick_cfg_val}` is not a valid entry "
                            f"for quick config shortcut parameter `{quick_cfg_k}`!"
                        )

                    # We only use the `quick_config` value if the parameter is not
                    # also being set explicitly by the user
                    for param_k, param_val in self.quick_config[quick_cfg_k][
                        quick_cfg_val
                    ].items():
                        if param_k not in self.given_configs:
                            setattr(self, param_k, param_val)

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

    def write(self, fl: str) -> None:
        """Saving configurations to file using the original order."""

        with open(fl, "w") as f:
            yaml.dump(asdict(self), f, default_flow_style=False, sort_keys=False)

    @classmethod
    def fields(cls) -> list[Field]:
        """Returning all fields defined for this class without needing an instance.

        The default Python dataclass `fields` method does not have a counterpart for
        classes, which we need in cases like `parse_cfg_keys()` which we want to call
        without using an instance of the data class!

        """
        members = inspect.getmembers(cls)
        return list(
            list(filter(lambda x: x[0] == "__dataclass_fields__", members))[0][
                1
            ].values()
        )

    @property
    def fields_dict(self) -> dict[str, Field]:
        return {fld.name: fld for fld in fields(self)}

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
                for fld in cls.fields():
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
                        cfg_key, [fld.name for fld in cls.fields()]
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
    configs (dict):     The raw configuration parameters for this engine.
                        Will be parsed by the

    Attributes
    ----------
    config_cls:         The configuration class that will be used to parse parameter
                        sets for this engine class. Children implementations of this
                        class will thus use children of `BaseConfigurations` here.
    label:              String used to refer to this engine for e.g. logging messages.

    configs (BaseConfigurations):    The parsed parameter configs for this engine.
    outdir  (str):      The path where output produced by the engine will be saved.
    logger  (Logger):   Logging utility used to create info and warning messages.
    """

    config_cls = BaseConfigurations
    label = "cDRGN training"

    @classmethod
    def parse_args(cls, args: argparse.Namespace) -> Self:
        """Utility for initializing using a namespace as opposed to a dictionary."""
        return cls(
            {
                par: (
                    getattr(args, par)
                    if hasattr(args, par)
                    else cls.config_cls.defaults[par]
                )
                for par in tuple(cls.config_cls.defaults)
            }
        )

    def __init__(self, configs: dict[str, Any]) -> None:
        if "load" in configs and configs["load"] == "latest":
            configs = self.get_latest_configs()

        self.configs = self.config_cls(**configs)
        self.outdir = self.configs.outdir
        np.random.seed(self.configs.seed)
        self.logger = logging.getLogger(self.label)

    @classmethod
    def defaults(cls) -> dict[str, Any]:
        """The user-set parameters governing the behaviour of this model."""
        return {fld.name: fld.default for fld in fields(cls.config_cls)}

    @classmethod
    def parameters(cls) -> list[str]:
        """The user-set parameters governing the behaviour of this model."""
        return [fld.name for fld in fields(cls.config_cls)]
