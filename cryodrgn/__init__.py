import os

# The _version.py file is managed by setuptools-scm
#   and is not in version control.
from ._version import version as __version__

_ROOT = os.path.abspath(os.path.dirname(__file__))
