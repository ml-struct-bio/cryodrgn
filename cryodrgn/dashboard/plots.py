"""Plotly figures for the dashboard; pair grid uses Matplotlib + Seaborn (PNG).

Implementation is split across sibling modules (``plots_*``) to keep files small and
avoid import cycles; this module re-exports the public API.
"""

from __future__ import annotations

from cryodrgn.dashboard.palette_config import (
    DEFAULT_DASHBOARD_CONTINUOUS_PALETTE,
    mpl_cmap_for_palette,
    normalize_continuous_palette,
)
from cryodrgn.dashboard.plots_color_covariate import (
    _continuous_series_stats,
    _labels_colors_and_legend_items,
    _lower_legend_entry_label,
    covariate_legend_context_payload,
    plot_df_color_filter_mask,
)
from cryodrgn.dashboard.plots_figure_utils import (
    _DASHBOARD_CREAM,
    _PLOTLY_FONT,
    _plotly_to_json,
    plot_df_row_indices_for_explorer_scatter,
)
from cryodrgn.dashboard.plots_pair_grid import (
    pair_grid_figure_aspect_ratio,
    pair_grid_margin_fractions_for_js,
    pair_grid_png,
    pair_grid_skeleton_placeholder_layout,
)
from cryodrgn.dashboard.plots_scatter import (
    scatter3d_z_json,
    scatter3d_z_preview_png,
    scatter_json,
)

__all__ = [
    "DEFAULT_DASHBOARD_CONTINUOUS_PALETTE",
    "mpl_cmap_for_palette",
    "normalize_continuous_palette",
    "_continuous_series_stats",
    "_labels_colors_and_legend_items",
    "_lower_legend_entry_label",
    "covariate_legend_context_payload",
    "plot_df_color_filter_mask",
    "_DASHBOARD_CREAM",
    "_PLOTLY_FONT",
    "_plotly_to_json",
    "plot_df_row_indices_for_explorer_scatter",
    "pair_grid_figure_aspect_ratio",
    "pair_grid_margin_fractions_for_js",
    "pair_grid_png",
    "pair_grid_skeleton_placeholder_layout",
    "scatter3d_z_json",
    "scatter3d_z_preview_png",
    "scatter_json",
]
