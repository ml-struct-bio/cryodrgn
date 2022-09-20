import os

# The _version.py file is managed by setuptools-scm
#   and is not in version control.
try:
    from ._version import version as __version__
except ModuleNotFoundError:
    # We're likely running as a source package without installation
    __version__ = 'src'

_ROOT = os.path.abspath(os.path.dirname(__file__))
