"""Matplotlib rc settings aligned with https://ezlab.princeton.edu/ (Barlow / Roboto)."""

from __future__ import annotations

from contextlib import contextmanager

import matplotlib.pyplot as plt

# Same ordering as ezlab main.css (Barlow first, Roboto second); local fallbacks for Agg.
EZLAB_SANS_SERIF = [
    "Barlow",
    "Roboto",
    "DejaVu Sans",
    "Helvetica",
    "Arial",
    "sans-serif",
]


@contextmanager
def ezlab_matplotlib_rc():
    with plt.rc_context(
        {
            "font.family": "sans-serif",
            "font.sans-serif": EZLAB_SANS_SERIF,
        }
    ):
        yield
