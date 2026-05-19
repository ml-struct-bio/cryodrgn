"""Shared dataframe column name constants for dashboard views."""

from __future__ import annotations

# Row index in ``DashboardExperiment.plot_df`` for each sampled point (hover ``row``).
VOL_LANDSCAPE_3D_PLOT_DF_ROW = "_dashboard_plot_df_row"
# Nearest k-means sketch ``vol_NNN.mrc`` index (L2 in shared volume-PCA dims) for GIF animations.
VOL_LANDSCAPE_NEAREST_SKETCH_VOL = "_dashboard_nearest_sketch_vol"
# ``1`` when this particle row is the ``centers_ind.txt`` anchor for its nearest sketch volume.
VOL_LANDSCAPE_IS_SKETCH_CENTROID = "_dashboard_sketch_centroid_point"
