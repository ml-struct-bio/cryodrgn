import logging


logger = logging.getLogger(__name__)


def numeric_object(which):
    if which == "cupy":
        from .cupy import Cupy as NumericClass
    elif which == "numpy":
        from .numpy import Numpy as NumericClass
    elif which == "torch":
        from .torch import Torch as NumericClass
    else:
        raise RuntimeError(f"Invalid selection for numeric module: {which}")
    return NumericClass()


xp = numeric_object("torch")


def fft_object(which):
    if which == "pyfftw":
        from .pyfftw_fft import PyfftwFFT as FFTClass
    elif which == "cupy":
        from .cupy_fft import CupyFFT as FFTClass
    elif which == "scipy":
        from .scipy_fft import ScipyFFT as FFTClass
    elif which == "numpy":
        from .numpy_fft import NumpyFFT as FFTClass
    elif which == "torch":
        from .torch_fft import TorchFFT as FFTClass
    else:
        raise RuntimeError(f"Invalid selection for fft class: {which}")
    return FFTClass()


fft = fft_object("torch")


