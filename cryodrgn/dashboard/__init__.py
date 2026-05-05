"""Web dashboard for cryoDRGN interactive analyses.

Application factory: :func:`cryodrgn.dashboard.app.create_app`. HTML templates live
under ``cryodrgn/dashboard/templates/``; shared CSS/JS under ``static/``. The
Flask layer stays thin: experiment state and caches are in
:mod:`cryodrgn.dashboard.context`; numerical and I/O work in sibling modules
(``plots``, ``preload``, ``particle_explorer``, etc.).
"""

# Agg must be selected before any submodule imports matplotlib.pyplot.
import matplotlib as _mpl

_mpl.use("Agg")
del _mpl
