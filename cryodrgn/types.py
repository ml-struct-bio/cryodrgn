"""Helper module to define cryoDRGN-specific types for static type-checking.
Usage of this module is recommended for all type annotations instead of direct
imports of types, to avoid circular import issues.

It is recommended that client code use `import cryodrgn.types` instead of reaching-in for types.
"""
from typing import TYPE_CHECKING, Union, List
import numpy as np


if TYPE_CHECKING:  # Set if type-checking
    # Avoid importing any cryodrgn-specific submodules outside this block.
    from cryodrgn.mrc import LazyImage

    ImageArray = Union[np.ndarray, List[LazyImage]]

else:
    """
    Any names declared above should also be declared here.
    This is to cover us in cases where client codes reaches in to get a hold of specific types.
    """
    ImageArray = None
