"""Construct the shortest path along a nearest neighbor graph in the latent z-space.

Example usages
--------------
# find the graph path between the first two points in the file, print to screen
$ cryodrgn graph_traversal zvals.pkl --anchors 0 1

# find the graph path connecting the first three points in order, save to file
$ cryodrgn graph_traversal zvals.pkl --anchors 0 1 2 -o

# connect the points whose indices are listed in a file; save only path indices
$ cryodrgn graph_traversal zvals.pkl --anchors path-anchors.txt --outind path-ind.txt

"""
import argparse
import os
import pickle
import logging
import numpy as np
import pandas as pd
from heapq import heappop, heappush
import torch
from cryodrgn.commands.direct_traversal import parse_anchors

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("zfile", help="Input .pkl file containing z-embeddings")
    parser.add_argument(
        "--anchors",
        required=True,
        nargs="+",
        help="List of anchor point indices in the given file, either given "
        "directly as integers, or in a .txt file(s).",
    )

    parser.add_argument("--max-neighbors", type=int, default=10)
    parser.add_argument("--avg-neighbors", type=float, default=5)
    parser.add_argument("--batch-size", type=int, default=1000)
    parser.add_argument("--max-images", type=int, default=None)

    parser.add_argument(
        "--outtxt",
        "-o",
        type=os.path.abspath,
        nargs="?",
        const="z-path.txt",
        metavar="Z-PATH.TXT",
        help="output .txt file for path z-values; "
        "choose name automatically if flag given with no name",
    )
    parser.add_argument(
        "--outind",
        type=os.path.abspath,
        nargs="?",
        const="z-path-indices.txt",
        metavar="IND-PATH.TXT",
        help="Output .txt file for path indices",
    )


class GraphLatentTraversor(object):
    def __init__(self, edges):  # edges is a list of tuples (src, dest, distance)
        # everything after here is derived from (weights, actions, probs)
        # for computational efficiency

        # FIXME: nodes and goal nodes should be the same
        self.nodes = set([x[0] for x in edges] + [x[1] for x in edges])
        self.edges = {x: set() for x in self.nodes}
        self.edge_length = {}
        for s, d, L in edges:
            assert type(s) == int and type(d) == int and type(L) == float
            self.edges[s].add(d)
            self.edge_length[(s, d)] = L

    def find_path(self, src, dest):
        visited = set()
        unvisited = []
        distances = {}
        predecessors = {}

        distances[src] = 0
        heappush(unvisited, (0, src))

        while unvisited:
            # visit the neighbors
            dist, v = heappop(unvisited)
            if v in visited or v not in self.edges:
                continue
            visited.add(v)
            if v == dest:
                # We build the shortest path and display it
                path = []
                pred = v
                while pred is not None:
                    path.append(pred)
                    pred = predecessors.get(pred, None)
                return path[::-1], dist

            neighbors = list(self.edges[v])

            for idx, neighbor in enumerate(neighbors):
                if neighbor not in visited:
                    new_distance = dist + self.edge_length[(v, neighbor)]
                    if new_distance < distances.get(neighbor, float("inf")):
                        distances[neighbor] = new_distance
                        heappush(unvisited, (new_distance, neighbor))
                        predecessors[neighbor] = v

        # couldn't find a path
        return None, None


def main(args):
    data_np = pickle.load(open(args.zfile, "rb"))
    data = torch.from_numpy(data_np)

    if args.max_images is not None:
        data = data[: args.max_images]

    use_cuda = torch.cuda.is_available()
    print(f"Use cuda {use_cuda}")
    device = torch.device("cuda" if use_cuda else "cpu")
    data = data.to(device)
    N, D = data.shape

    anchors = parse_anchors(args.anchors, data, args.zfile)
    n2 = (data * data).sum(-1, keepdim=True)
    B = args.batch_size
    ndist = torch.empty(data.shape[0], args.max_neighbors, device=device)
    neighbors = torch.empty(
        data.shape[0], args.max_neighbors, dtype=torch.long, device=device
    )
    for i in range(0, data.shape[0], B):
        # (a-b)^2 = a^2 + b^2 - 2ab
        print(f"Working on images {i}-{min(data.shape[0], i+B)}")
        batch_dist = n2[i : i + B] + n2.t() - 2 * torch.mm(data[i : i + B], data.t())
        ndist[i : i + B], neighbors[i : i + B] = batch_dist.topk(
            args.max_neighbors, dim=-1, largest=False
        )

    assert ndist.min() >= -1e-3, ndist.min()

    # convert d^2 to d
    ndist = ndist.clamp(min=0).pow(0.5)
    if args.avg_neighbors:
        total_neighbors = int(N * args.avg_neighbors)
        max_dist = ndist.view(-1).topk(total_neighbors, largest=False)[0][-1]
    else:
        max_dist = None
    print(
        f"Max dist between neighbors: {max_dist:.4g}  "
        f"(to enforce average of {args.avg_neighbors} neighbors)"
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
    full_path = []
    data_df = pd.DataFrame()

    for i in range(len(anchors) - 1):
        anchor_str = f"{anchors[i]}->{anchors[i + 1]}"
        src, dest = anchors[i], anchors[i + 1]
        path, total_distance = graph.find_path(src, dest)
        dd = data[path].cpu().numpy()
        dists = ((dd[1:, :] - dd[0:-1, :]) ** 2).sum(axis=1) ** 0.5

        if path is not None:
            if full_path and full_path[-1] == path[0]:
                path = path[1:]
            else:
                dists = [0] + dists.tolist()

            new_df = pd.DataFrame(
                data_np[path],
                index=path,
                columns=[f"z{i + 1}" for i in range(D)],
            )

            new_df["dist"] = dists
            data_df = pd.concat([data_df, new_df])
            full_path += path

            euc_dist = ((dd[0] - dd[-1]) ** 2).sum() ** 0.5
            print(f"Total path distance {anchor_str}: {total_distance:.4g}")
            print(f" Euclidean distance {anchor_str}: {euc_dist:.4g}")
        else:
            print(f"Could not find a {anchor_str} path!")

    data_df.index = [f"{i}(a)" if i in anchors else str(i) for i in data_df.index]
    data_df.index.name = "ind"

    if args.outind:
        if not os.path.exists(os.path.dirname(args.outind)):
            os.makedirs(os.path.dirname(args.outind))
        logger.info(f"Saving path indices relative to {args.zfile} to {args.outind}!")
        np.savetxt(args.outind, full_path, fmt="%d")
    elif args.outtxt:
        logger.info(f"Found shortest nearest-neighbor path with indices:\n{full_path}")

    if args.outtxt:
        if not os.path.exists(os.path.dirname(args.outtxt)):
            os.makedirs(os.path.dirname(args.outtxt))
        logger.info(f"Saving path z-values to {args.outtxt}!")
        np.savetxt(args.outtxt, data_np[full_path])
    else:
        if args.outind:
            logger.info(f"Found shortest nearest-neighbor path:\n{data_np[full_path]}")
        else:
            print_data = data_df.round(3).to_csv(sep="\t")
            logger.info(f"Found shortest nearest-neighbor path:\n{print_data}")
