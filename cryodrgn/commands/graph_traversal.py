'''
Find shortest path along nearest neighbor graph
'''
import torch
import argparse
import pickle
import numpy as np
import os

from heapq import heappush, heappop

def add_args(parser):
    parser.add_argument('data', help='Input z.pkl embeddings')
    parser.add_argument('--anchors', type=int, nargs='+', required=True, help='Index of anchor points')
    parser.add_argument('--max-neighbors', type=int, default=10)
    parser.add_argument('--avg-neighbors', type=float, default=5)
    parser.add_argument('--batch-size', type=int, default=1000)
    parser.add_argument('--max-images', type=int, default=None)
    parser.add_argument('-o', metavar='PATH.TXT', type=os.path.abspath, required=True, help='Output .txt file for path indices')
    parser.add_argument('--out-z', metavar='Z.PATH.TXT', type=os.path.abspath, required=True, help='Output .txt file for path z-values')
    return parser

class Graph(object):

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
                    if new_distance < distances.get(neighbor, float('inf')):
                        distances[neighbor] = new_distance
                        heappush(unvisited, (new_distance, neighbor))
                        predecessors[neighbor] = v

        # couldn't find a path
        return None, None



def main(args):
    data_np = pickle.load(open(args.data, 'rb'))
    data = torch.from_numpy(data_np)

    if args.max_images is not None:
        data = data[:args.max_images]

    use_cuda = torch.cuda.is_available()
    print(f'Use cuda {use_cuda}')
    device = torch.device('cuda' if use_cuda else 'cpu')
    data = data.to(device)

    N, D = data.shape
    for i in args.anchors:
        assert i < N
    assert len(args.anchors) >= 2

    n2 = (data*data).sum(-1, keepdim=True)
    B = args.batch_size
    ndist = torch.empty(data.shape[0], args.max_neighbors, device=device)
    neighbors = torch.empty(data.shape[0], args.max_neighbors, dtype=torch.long, device=device)
    for i in range(0, data.shape[0], B):
        # (a-b)^2 = a^2 + b^2 - 2ab
        print(f"Working on images {i}-{i+B}")
        batch_dist = n2[i:i+B] + n2.t() - 2 * torch.mm(data[i:i+B], data.t())
        ndist[i:i+B], neighbors[i:i+B] = batch_dist.topk(args.max_neighbors, dim=-1, largest=False)

    assert ndist.min() >= -1e-3, ndist.min()

    # convert d^2 to d
    ndist = ndist.clamp(min=0).pow(0.5)
    if args.avg_neighbors:
        total_neighbors = int(N * args.avg_neighbors)
        max_dist = ndist.view(-1).topk(total_neighbors, largest=False)[0][-1]
    else:
        max_dist = None
    print(f"Max dist between neighbors: {max_dist}  (to enforce average of {args.avg_neighbors} neighbors)")

    max_dist = max_dist.to("cpu")
    neighbors = neighbors.to("cpu")
    ndist = ndist.to("cpu")
    edges = []
    for i in range(neighbors.shape[0]):
        for j in range(neighbors.shape[1]):
            if max_dist is None or ndist[i, j] < max_dist:
                edges.append((int(i), int(neighbors[i, j]), float(ndist[i, j])))

    graph = Graph(edges)
    full_path = []
    for i in range(len(args.anchors)-1):
        src, dest = args.anchors[i], args.anchors[i+1]
        path, total_distance = graph.find_path(src, dest)
        dd = data[path].cpu().numpy()
        dists = ((dd[1:,:] - dd[0:-1,:])**2).sum(axis=1)**.5
        
        if path is not None:
            if full_path and full_path[-1] == path[0]:
                full_path.extend(path[1:])
            else:
                full_path.extend(path)

        print()
        if path is not None:
            print("Path:")
            for id in path:
                print(id)
            print()
            print(f"Total distance: {total_distance}")
            print()
            print('Neighbor distance:')
            for d in dists:
                print(d)
            print()
            print('Euclidean distance: {}'.format(((dd[0]-dd[-1])**2).sum()**.5))
        else:
            print("Could not find path!")
    
    if not os.path.exists(os.path.dirname(args.o)):
        os.makedirs(os.path.dirname(args.o))
    print(args.o)
    np.savetxt(args.o,full_path,fmt='%d')
    print(args.out_z)
    np.savetxt(args.out_z, data_np[full_path])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
