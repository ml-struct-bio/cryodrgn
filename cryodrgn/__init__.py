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
            "standard": {
                "format": "(%(levelname)s) (%(filename)s) (%(asctime)s) %(message)s",
                "datefmt": "%d-%b-%y %H:%M:%S",
            }
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

# train_vae /media/vineetb/t5/10076_128/particles.128.txt --poses /home/vineetb/cryodrgn/cryodrgn_empiar/empiar10076/inputs/poses.pkl --ctf /home/vineetb/cryodrgn/cryodrgn_empiar/empiar10076/inputs/ctf.pkl -o benchmark0 --zdim 16 --enc-dim 32 --dec-dim 32 -n 1 --max-threads 8 --no-amp --num-workers-per-gpu 1 --log-interval 1 --batch-size 32 --lazy
# works
#   old/lazy/non-prealloc  (400 it/s)
#   old/lazy/prealloc      (397 it/s)
#   old/eager/non-prealloc (397 it/s)
#   old/eager/prealloc     (397 it/s)
#   new/lazy/non-prealloc  (405 it/s)
#   new/lazy/prealloc      (405 it/s)
#   new/eager/non-prealloc (400 it/s)
#   new/eager/prealloc     (400 it/s)
