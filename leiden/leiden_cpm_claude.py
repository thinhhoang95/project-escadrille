import random
from collections import defaultdict

class Graph:
    def __init__(self):
        self.adj = defaultdict(set)
        self.weights = defaultdict(float)
        self.node_weights = defaultdict(float)
        self.total_weight = 0.0

    def add_edge(self, u, v, weight=1.0):
        self.adj[u].add(v)
        self.adj[v].add(u)
        self.weights[(u, v)] = weight
        self.weights[(v, u)] = weight
        self.total_weight += weight

    def nodes(self):
        return list(self.adj.keys())

    def neighbors(self, node):
        return self.adj[node]

class Partition:
    def __init__(self, graph):
        self.graph = graph
        self.membership = {node: node for node in graph.nodes()}
        self.csize = {node: 1 for node in graph.nodes()}
        self.total_weight_in_comm = defaultdict(float)
        self.total_weight_from_comm = defaultdict(float)

    def move_node(self, node, new_comm):
        old_comm = self.membership[node]
        if old_comm != new_comm:
            self.membership[node] = new_comm
            self.csize[old_comm] -= 1
            self.csize[new_comm] += 1
            if self.csize[old_comm] == 0:
                del self.csize[old_comm]

            for neighbor in self.graph.neighbors(node):
                weight = self.graph.weights[(node, neighbor)]
                old_neigh_comm = self.membership[neighbor]
                
                self.total_weight_in_comm[old_comm] -= weight
                self.total_weight_from_comm[old_comm] -= weight
                
                if old_neigh_comm == new_comm:
                    self.total_weight_in_comm[new_comm] += weight
                else:
                    self.total_weight_from_comm[new_comm] += weight

                if old_neigh_comm == old_comm:
                    self.total_weight_in_comm[old_comm] -= weight
                    self.total_weight_from_comm[old_comm] += weight

    def quality(self, resolution):
        q = 0.0
        for comm in self.csize:
            q += self.total_weight_in_comm[comm] - resolution * self.csize[comm] * (self.csize[comm] - 1) / 2
        return q

class Optimiser:
    def __init__(self, resolution=1.0):
        self.resolution = resolution

    def move_nodes(self, partition):
        improved = False
        nodes = list(partition.graph.nodes())
        random.shuffle(nodes)

        for node in nodes:
            old_comm = partition.membership[node]
            best_comm = old_comm
            best_quality = 0.0

            # Consider moving to neighboring communities
            considered_comms = set()
            for neighbor in partition.graph.neighbors(node):
                neigh_comm = partition.membership[neighbor]
                if neigh_comm not in considered_comms:
                    quality = self.calculate_delta_quality(partition, node, neigh_comm)
                    if quality > best_quality:
                        best_quality = quality
                        best_comm = neigh_comm
                    considered_comms.add(neigh_comm)

            if best_comm != old_comm:
                partition.move_node(node, best_comm)
                improved = True

        return improved

    def calculate_delta_quality(self, partition, node, new_comm):
        old_comm = partition.membership[node]
        delta_q = 0.0

        # Calculate the change in quality
        for neighbor in partition.graph.neighbors(node):
            neigh_comm = partition.membership[neighbor]
            weight = partition.graph.weights[(node, neighbor)]
            if neigh_comm == new_comm:
                delta_q += weight
            if neigh_comm == old_comm:
                delta_q -= weight

        delta_q -= self.resolution * (partition.csize[new_comm] - partition.csize[old_comm] + 1)

        return delta_q

def leiden(graph, resolution=1.0, max_iter=10):
    partition = Partition(graph)
    optimiser = Optimiser(resolution)

    for _ in range(max_iter):
        improved = optimiser.move_nodes(partition)
        if not improved:
            break

        # Aggregate the graph
        new_graph = Graph()
        comm_to_node = {}
        for node, comm in partition.membership.items():
            if comm not in comm_to_node:
                comm_to_node[comm] = len(comm_to_node)
            new_node = comm_to_node[comm]
            for neighbor in graph.neighbors(node):
                neigh_comm = partition.membership[neighbor]
                new_neigh = comm_to_node[neigh_comm]
                if new_neigh != new_node:
                    weight = graph.weights[(node, neighbor)]
                    new_graph.add_edge(new_node, new_neigh, weight)

        graph = new_graph
        partition = Partition(graph)

    return partition.membership

# Example usage
if __name__ == "__main__":
    g = Graph()
    edges = [
        (1, 2), (1, 3), (2, 3),
        (4, 5), (5, 6), (4, 6),
        (3, 4),
        (7, 8), (8, 9), (7, 9),
        (9, 4)
    ]
    for u, v in edges:
        g.add_edge(u, v)

    communities = leiden(g, resolution=0.5)
    print("Final communities:", communities)