"""Central cache ownership for dashboard experiment and preload state.

The process-wide :data:`EXPERIMENT_STORE` singleton holds loaded
:class:`~cryodrgn.dashboard.data.DashboardExperiment` instances, explorer
thumbnail preload batches, per-workdir epoch lists, and trajectory k-NN graph
arrays.  Invalidation is coordinated through :meth:`ExperimentStore.clear_all`
and :meth:`ExperimentStore.clear_preloads_for_experiment`.
"""

from __future__ import annotations
from cryodrgn.dashboard.data import DashboardExperiment, list_z_epochs, load_experiment

ExperimentKey = tuple[str, int, int]
PreloadKey = tuple[int, int, str, str, tuple[int, ...] | None]
PreloadValue = tuple[list[int], list[str], float]
TrajGraphKey = tuple[str, int, int, int]
TrajGraphValue = tuple["np.ndarray", "np.ndarray", "_csr_matrix"]


class ExperimentStore:
    """Owns long-lived dashboard caches for a single server process."""

    def __init__(self) -> None:
        self.experiments: dict[ExperimentKey, DashboardExperiment] = {}
        self.preloads: dict[PreloadKey, PreloadValue] = {}
        self.traj_graph_neighbors: dict[TrajGraphKey, TrajGraphValue] = {}
        self._epochs_by_workdir: dict[str, list[int]] = {}

    def get_experiment(
        self,
        workdir: str,
        epoch: int,
        kmeans: int,
    ) -> DashboardExperiment:
        """Return a cached experiment, loading from disk on first access."""
        key: ExperimentKey = (workdir, int(epoch), int(kmeans))
        if key not in self.experiments:
            self.experiments[key] = load_experiment(
                workdir, epoch=int(epoch), kmeans=int(kmeans)
            )
        return self.experiments[key]

    def epochs_for_workdir(self, workdir: str) -> list[int]:
        """Epochs with both ``z.N.pkl`` and ``analyze.N/`` (memoised per workdir)."""
        cached = self._epochs_by_workdir.get(workdir)
        if cached is not None:
            return cached
        epochs = list_z_epochs(workdir)
        self._epochs_by_workdir[workdir] = epochs
        return epochs

    def clear_epochs_for_workdir(self, workdir: str) -> None:
        """Drop the memoised epoch list for one workdir."""
        self._epochs_by_workdir.pop(workdir, None)

    def clear_preloads_for_experiment(self, exp: DashboardExperiment) -> int:
        """Drop explorer thumbnail preload entries for this epoch and k-means id.

        Preload keys are ``(epoch, kmeans_folder_id, xcol, ycol, selection)``.
        """
        prefix = (int(exp.epoch), int(exp.kmeans_folder_id))
        to_drop = [k for k in list(self.preloads.keys()) if k[:2] == prefix]
        for key in to_drop:
            del self.preloads[key]
        return len(to_drop)

    def clear_all(self) -> None:
        """Drop experiments, preloads, and trajectory graph neighbours."""
        self.experiments.clear()
        self.preloads.clear()
        self.traj_graph_neighbors.clear()


EXPERIMENT_STORE = ExperimentStore()
