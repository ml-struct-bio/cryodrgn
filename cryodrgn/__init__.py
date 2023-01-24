import os

# The _version.py file is managed by setuptools-scm
#   and is not in version control.
try:
    from cryodrgn._version import version as __version__  # type: ignore
except ModuleNotFoundError:
    # We're likely running as a source package without installation
    __version__ = "src"

_ROOT = os.path.abspath(os.path.dirname(__file__))

# Temporary - till everything works regardless of whether the following is True/False
USE_NEW_DATASET_API = True
