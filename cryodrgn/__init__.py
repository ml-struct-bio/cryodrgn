import os
import logging.config


# The _version.py file is managed by setuptools-scm
#   and is not in version control.
try:
    from cryodrgn._version import version as __version__  # type: ignore
except ModuleNotFoundError:
    # We're likely running as a source package without installation
    __version__ = "src"

_ROOT = os.path.abspath(os.path.dirname(__file__))

logging.config.dictConfig(
    {
        "version": 1,
        "formatters": {
            "standard": {"format": "%(asctime)s %(message)s (%(filename)s:%(lineno)d)"}
        },
        "handlers": {
            "default": {
                "level": "NOTSET",
                "formatter": "standard",
                "class": "logging.StreamHandler",
                "stream": "ext://sys.stdout",
            }
        },
        "loggers": {"": {"handlers": ["default"], "level": "INFO"}},
    }
)

# Temporary - till everything works regardless of whether the following is True/False
USE_NEW_DATASET_API = False
# Use preallocated arrays with extra dimension in x/y to reduce memory footprint during fft symmetrizing operations?
PREALLOCATED = True

# train_vae /media/vineetb/t5/10076_128/particles.128.txt --poses /home/vineetb/cryodrgn/cryodrgn_empiar/empiar10076/inputs/poses.pkl --ctf /home/vineetb/cryodrgn/cryodrgn_empiar/empiar10076/inputs/ctf.pkl -o benchmark0 --zdim 16 --enc-dim 32 --dec-dim 32 -n 1 --max-threads 8 --no-amp --num-workers-per-gpu 1 --log-interval 1 --batch-size 32 --lazy
# works
#   old/non-prealloc (0:47:26.111794 / 0:58:05.464952)
#   old/prealloc     (0:42:42.554901 / 0:55:22.231400)
#   new/non-prealloc (0:47:35.429347 / 0:58:46.271654)
# doesn't work
