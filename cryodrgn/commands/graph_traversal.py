"""Construct the shortest path along a nearest neighbor graph in the latent z-space.

Example usage
-------------
# Find the path between kmeans cluster centers from cryodrgn analyze
$ cryodrgn graph_traversal my_workdir/z.50.pkl \
                           --anchors my_workdir/analyze.50/centers_ind.txt \
                           -o graph_traversal/z-path.txt \
                           --outind graph_traversal/z-path.ind.txt

See also
--------
`cryodrgn eval_vol` for generating volumes from the path (See the cryodrgn docs for more info)
`cryodrgn direct_traversal` for direct interpolation between points

"""
import argparse
import os
import pickle
import logging
import numpy as np
import pandas as pd
from typing import List, Tuple
from heapq import heappop, heappush
import torch
from cryodrgn.commands.direct_traversal import parse_anchors
from datetime import datetime as dt

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    """The command-line arguments for use with `cryodrgn graph_traversal`."""

    parser.add_argument("zfile", help="Input .pkl file containing z-embeddings")
    parser.add_argument(
        "--anchors",
        required=True,
        nargs="+",
        help="List of anchor point indices in the given file, either given "
        "directly as integers, or in a .txt file(s).",
    )
    parser.add_argument(
        "--max-neighbors",
        type=int,
        default=10,
        help="Maximum number of nearest neighbors to consider when constructing "
        "the nearest neighbor graph. This limits the number of connections per node, "
        "improving computational efficiency. (default: %(default)s)",
    )
    parser.add_argument(
        "--avg-neighbors",
        type=float,
        default=5,
        help="Average number of neighbors to aim for each point in the graph. "
        "This parameter adjusts the maximum distance to ensure a balanced graph "
        "density. (default: %(default)s)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=1000,
        help="Number of data points to process at once when computing nearest "
        "neighbors. Larger batch sizes may improve speed but require more memory. "
        "(default: %(default)s)",
    )

    parser.add_argument(
        "--outtxt",
        "-o",
        type=os.path.abspath,
        default="z-path.txt",
        metavar="Z-PATH.TXT",
        help="Output .txt file for path z-values (default: %(default)s)",
    )
    parser.add_argument(
        "--outind",
        type=os.path.abspath,
        default="z-path-indices.txt",
        metavar="IND-PATH.TXT",
        help="Output .txt file for path indices (default: %(default)s)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print path values, indices to screen",
    )


class GraphLatentTraversor:
    """An engine for finding a path between images in a z-latent-space embedding."""

    def __init__(self, edges: List[Tuple[int, int, float]]) -> None:
        """
        Everything after here is derived from (weights, actions, probs)
        for computational efficiency

        Arguments
        ---------
            edges: A list of tuples (src, dest, distance)

        """
        # FIXME: nodes and goal nodes should be the same
        self.nodes = set([x[0] for x in edges] + [x[1] for x in edges])
        self.edges = {x: set() for x in self.nodes}
        self.edge_length = dict()

        for s, d, L in edges:
            assert isinstance(s, int) and isinstance(d, int) and isinstance(L, float)
            self.edges[s].add(d)
            self.edge_length[(s, d)] = L

    def find_path(self, src: int, dest: int) -> Tuple[List[int], float]:
        """Find the shortest path between two nodes in the graph."""

        visited = set()
        unvisited = list()
        distances = dict()
        predecessors = dict()
        distances[src] = 0
        heappush(unvisited, (0, src))

        while unvisited:
            # Visit the neighbors
            dist, v = heappop(unvisited)
            if v in visited or v not in self.edges:
                continue

            visited.add(v)
            if v == dest:
                path = list()
                pred = v

                while pred is not None:
                    path.append(pred)
                    pred = predecessors.get(pred, None)

                return path[::-1], dist

            for neighbor in self.edges[v]:
                if neighbor not in visited:
                    new_distance = dist + self.edge_length[(v, neighbor)]

                    if new_distance < distances.get(neighbor, float("inf")):
                        distances[neighbor] = new_distance
                        heappush(unvisited, (new_distance, neighbor))
                        predecessors[neighbor] = v

        # Couldn't find a path
        return None, None


def main(args: argparse.Namespace) -> None:
    """Find the shortest path in a z-space neighbor graph (see `add_args` above)."""
    zdata_np = pickle.load(open(args.zfile, "rb"))
    zdata = torch.from_numpy(zdata_np)

    use_cuda = torch.cuda.is_available()
    print(f"Use cuda {use_cuda}")
    device = torch.device("cuda" if use_cuda else "cpu")
    zdata = zdata.to(device)
    N, D = zdata.shape

    anchors = parse_anchors(args.anchors, zdata, args.zfile)
    n2 = (zdata * zdata).sum(-1, keepdim=True)
    B = min(args.batch_size, N)
    max_neighbors = min(args.max_neighbors, B)
    avg_neighbors = min(args.avg_neighbors, B)
    ndist = torch.empty(N, max_neighbors, device=device)
    neighbors = torch.empty(N, max_neighbors, dtype=torch.long, device=device)

    for i in range(0, N, B):
        # (a-b)^2 = a^2 + b^2 - 2ab
        print(f"Working on images {i}-{min(N, i+B)}")
        batch_dist = n2[i : i + B] + n2.t() - 2 * torch.mm(zdata[i : i + B], zdata.t())
        ndist[i : i + B], neighbors[i : i + B] = batch_dist.topk(
            max_neighbors, dim=-1, largest=False
        )

    assert ndist.min() >= -1e-3, ndist.min()

    # convert d^2 to d
    ndist = ndist.clamp(min=0).pow(0.5)
    if args.avg_neighbors:
        total_neighbors = int(N * avg_neighbors)
        max_dist = ndist.view(-1).topk(total_neighbors, largest=False)[0][-1]
    else:
        max_dist = None
    print(
        f"Max dist between neighbors: {max_dist:.4g}  "
        f"(to enforce average of {avg_neighbors} neighbors)"
    )

    if max_dist is not None:
        max_dist = max_dist.to("cpu")
    neighbors = neighbors.to("cpu")
    ndist = ndist.to("cpu")
    edges = []
    for i in range(neighbors.shape[0]):
        for j in range(neighbors.shape[1]):
            if max_dist is None or ndist[i, j] < max_dist:
                edges.append((int(i), int(neighbors[i, j]), float(ndist[i, j])))

    graph = GraphLatentTraversor(edges)
    full_path = list()
    zdata_df = pd.DataFrame()

    t1 = dt.now()

    for i in range(len(anchors) - 1):
        anchor_str = f"{anchors[i]}->{anchors[i + 1]}"
        src, dest = anchors[i], anchors[i + 1]
        path, total_distance = graph.find_path(src, dest)
        path_zdata = zdata[path].cpu().numpy()
        dists = ((path_zdata[1:, :] - path_zdata[0:-1, :]) ** 2).sum(axis=1) ** 0.5

        if path is not None:
            if full_path and full_path[-1] == path[0]:
                path = path[1:]
            else:
                dists = [0] + dists.tolist()

            new_df = pd.DataFrame(
                zdata_np[path],
                index=path,
                columns=[f"z{i + 1}" for i in range(D)],
            )
            new_df["dist"] = dists
            zdata_df = pd.concat([zdata_df, new_df])
            full_path += path

            euc_dist = ((path_zdata[0] - path_zdata[-1]) ** 2).sum() ** 0.5
            print(f"Total path distance {anchor_str}: {total_distance:.4g}")
            print(f" Euclidean distance {anchor_str}: {euc_dist:.4g}")
        else:
            logger.warning(f"Could not find a {anchor_str} path!")
            full_path = None
            break

    if full_path:
        if args.verbose:
            zdata_df.index = [
                f"{i}(a)" if i in anchors else str(i) for i in zdata_df.index
            ]
            zdata_df.index.name = "ind"
            print_data = zdata_df.round(3).to_csv(sep="\t")
            logger.info(f"Found shortest nearest-neighbor path:\n{print_data}")

        if args.outind:
            if not os.path.exists(os.path.dirname(args.outind)):
                os.makedirs(os.path.dirname(args.outind))
            logger.info(
                f"Saving path indices relative to {args.zfile} to {args.outind}"
            )
            np.savetxt(args.outind, full_path, fmt="%d")

        if args.outtxt:
            if not os.path.exists(os.path.dirname(args.outtxt)):
                os.makedirs(os.path.dirname(args.outtxt))
            logger.info(f"Saving path z-values to {args.outtxt}")
            np.savetxt(args.outtxt, zdata_np[full_path])

        t2 = dt.now()
        elapsed_time = t2 - t1
        logger.info(f"Graph traversal completed in {elapsed_time}")

    elif len(anchors) > 2:
        logger.warning(f"Could not find a path between {anchors[0]} and {anchors[-1]}!")
