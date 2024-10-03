from collections import defaultdict, deque
import copy

class Graph:
    def __init__(self):
        self.adj = defaultdict(set)
        self.weights = defaultdict(int)  # For weighted edges
        self.total_weight = 0

    def add_edge(self, u, v, weight=1):
        if v not in self.adj[u]:
            self.adj[u].add(v)
            self.adj[v].add(u)
            self.weights[(u, v)] = weight
            self.weights[(v, u)] = weight
            self.total_weight += weight
        else:
            # If edge already exists, update the weight
            self.weights[(u, v)] += weight
            self.weights[(v, u)] += weight
            self.total_weight += weight

    def nodes(self):
        return list(self.adj.keys())

    def neighbors(self, u):
        return self.adj[u]

def leiden(graph, gamma=1.0, resolution_parameter=None, seed=None):
    """
    Leiden algorithm implementation with CPM and recursive mapping.
    :param graph: Graph object
    :param gamma: CPM resolution parameter
    :return: list of communities (each community is a set of original nodes)
    """
    if resolution_parameter is not None:
        gamma = resolution_parameter

    # Initialize each node to its own community
    partition = {node: node for node in graph.nodes()}
    # Initialize community hierarchy: maps community ID to original nodes
    community_hierarchy = {node: {node} for node in graph.nodes()}
    
    # Initialize community metrics: m_c and n_c
    community_m = {node: 0 for node in graph.nodes()}  # Internal edges
    community_n = {node: 1 for node in graph.nodes()}  # Number of nodes

    # Compute initial internal edges for each community
    for node in graph.nodes():
        for neighbor in graph.neighbors(node):
            if partition[neighbor] == partition[node]:
                community_m[partition[node]] += graph.weights[(node, neighbor)]
    # Since each internal edge is counted twice in undirected graphs
    for comm in community_m:
        community_m[comm] = community_m[comm] // 2

    improvement = True

    while improvement:
        improvement = False
        # Local moving phase
        nodes = graph.nodes()
        for node in nodes:
            current_comm = partition[node]
            neighbor_comms = defaultdict(int)

            # Collect the weights to neighboring communities
            for neighbor in graph.neighbors(node):
                neighbor_comm = partition[neighbor]
                neighbor_comms[neighbor_comm] += graph.weights[(node, neighbor)]

            # Compute the best community to move to
            best_comm = current_comm
            best_gain = 0
            total_weight = graph.total_weight

            for comm, edge_weight in neighbor_comms.items():
                if comm == current_comm:
                    continue
                # Compute Delta CPM using the correct formula
                # ΔCPM = (2 * k_i_in_new_comm / m) - gamma * (n_new_comm + n_current_comm -1) / m
                k_i_in_new_comm = edge_weight
                k_i_in_current_comm = 0
                # Since node is currently in current_comm, we need to find k_i_in_current_comm
                # which is the number of edges node has within current_comm
                k_i_in_current_comm = neighbor_comms.get(current_comm, 0)

                n_current_comm = community_n[current_comm]
                n_new_comm = community_n[comm]

                delta_cpm = (2 * k_i_in_new_comm) / graph.total_weight - gamma * (n_new_comm + n_current_comm - 1) / graph.total_weight

                if delta_cpm > best_gain:
                    best_gain = delta_cpm
                    best_comm = comm

            if best_comm != current_comm and best_gain > 0:
                # Move node to best_comm
                partition[node] = best_comm
                improvement = True

                # Update community metrics
                # Decrease metrics of current_comm
                community_m[current_comm] -= k_i_in_current_comm
                community_n[current_comm] -= 1

                # Increase metrics of best_comm
                community_m[best_comm] += k_i_in_new_comm
                community_n[best_comm] += 1

        if improvement:
            # Refinement phase to ensure communities are well-connected
            partition = refine(graph, partition)

            # Recompute community metrics after refinement
            communities = defaultdict(list)
            for node, comm in partition.items():
                communities[comm].append(node)

            # Reset community metrics
            community_m = {}
            community_n = {}
            for comm, nodes_in_comm in communities.items():
                community_n[comm] = len(nodes_in_comm)
                internal_edges = 0
                for node in nodes_in_comm:
                    for neighbor in graph.neighbors(node):
                        if partition[neighbor] == comm:
                            internal_edges += graph.weights[(node, neighbor)]
                # Each internal edge is counted twice
                community_m[comm] = internal_edges // 2

            # Aggregation phase
            # Assign unique community IDs
            new_comm_ids = {comm: idx for idx, comm in enumerate(communities.keys())}

            # Update community hierarchy: map new communities to original nodes
            new_community_hierarchy = {}
            for comm, nodes_in_comm in communities.items():
                new_comm_id = new_comm_ids[comm]
                # Union of all original nodes in the constituent communities
                aggregated_nodes = set()
                for node in nodes_in_comm:
                    aggregated_nodes.update(community_hierarchy[node])
                new_community_hierarchy[new_comm_id] = aggregated_nodes

            # Build a new aggregated graph
            new_graph = Graph()
            # To track inter-community edge weights
            inter_comm_edges = defaultdict(int)

            for comm, nodes_in_comm in communities.items():
                for node in nodes_in_comm:
                    for neighbor in graph.neighbors(node):
                        neighbor_comm = partition[neighbor]
                        if neighbor_comm != comm:
                            # To avoid double-counting, ensure consistent ordering
                            sorted_comm = tuple(sorted((new_comm_ids[comm], new_comm_ids[neighbor_comm])))
                            inter_comm_edges[sorted_comm] += graph.weights[(node, neighbor)]

            # Add edges to the new graph based on inter-community edge weights
            for (comm1, comm2), weight in inter_comm_edges.items():
                new_graph.add_edge(comm1, comm2, weight)

            # Update total_weight for the new graph
            new_graph.total_weight = sum(new_graph.weights.values()) // 2  # Since undirected

            # Update community hierarchy
            community_hierarchy = new_community_hierarchy

            # Prepare for next iteration
            graph = new_graph
            # Reinitialize partition: each new community is its own community
            partition = {node: node for node in graph.nodes()}

            # Initialize community metrics for the new graph
            community_m = {node: 0 for node in graph.nodes()}  # Internal edges
            community_n = {node: 1 for node in graph.nodes()}  # Number of nodes

            for node in graph.nodes():
                for neighbor in graph.neighbors(node):
                    if partition[neighbor] == partition[node]:
                        community_m[partition[node]] += graph.weights[(node, neighbor)]
            # Since each internal edge is counted twice in undirected graphs
            for comm in community_m:
                community_m[comm] = community_m[comm] // 2

    # Final partition mapping back to original nodes
    final_communities = defaultdict(set)
    for node, comm in partition.items():
        # Each node here is an aggregated community; map to original nodes
        original_nodes = community_hierarchy[node]
        final_communities[comm].update(original_nodes)

    # Convert sets to lists for readability
    final_communities = [sorted(list(members)) for members in final_communities.values()]
    return final_communities

def refine(graph, partition):
    """
    Refinement step to ensure communities are well-connected.
    Splits communities if they are not internally connected.
    :param graph: Graph object
    :param partition: current partition mapping
    :return: refined partition mapping
    """
    communities = defaultdict(list)
    for node, comm in partition.items():
        communities[comm].append(node)

    refined_partition = {}
    new_comm_id = max(partition.values()) + 1  # Start new community IDs from here

    for comm, nodes in communities.items():
        # Check connectivity using BFS
        subgraph = build_subgraph(graph, nodes)
        components = get_connected_components(subgraph, nodes)
        if len(components) == 1:
            # Community is connected
            for node in components[0]:
                refined_partition[node] = comm
        else:
            # Split into separate communities
            for component in components:
                for node in component:
                    refined_partition[node] = new_comm_id
                new_comm_id += 1

    return refined_partition

def build_subgraph(graph, nodes):
    """
    Builds a subgraph containing only the specified nodes.
    :param graph: Graph object
    :param nodes: list of nodes to include in the subgraph
    :return: subgraph as a Graph object
    """
    sub = Graph()
    for u in nodes:
        for v in graph.neighbors(u):
            if v in nodes:
                sub.add_edge(u, v, graph.weights[(u, v)])
    return sub

def get_connected_components(subgraph, nodes):
    """
    Returns the list of connected components in the subgraph.
    :param subgraph: Graph object
    :param nodes: list of nodes in the subgraph
    :return: list of lists, each inner list is a connected component
    """
    visited = set()
    components = []

    for node in nodes:
        if node not in visited:
            component = []
            queue = deque()
            queue.append(node)
            visited.add(node)
            while queue:
                current = queue.popleft()
                component.append(current)
                for neighbor in subgraph.neighbors(current):
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
            components.append(component)
    return components

# Example Usage
if __name__ == "__main__":
    # Create a sample graph
    g = Graph()
    edges = [
        (1, 2), (1, 3), (2, 3),
        (4, 5), (5, 6), (4, 6),
        (3, 4),  # Connect the two communities
        (7, 8), (8, 9), (7, 9),
        (9, 4)   # Connect to the main community
    ]
    for u, v in edges:
        g.add_edge(u, v)

    # Run Leiden algorithm with default gamma
    print("Running Leiden algorithm with gamma=1.0")
    communities = leiden(g, gamma=1.0)
    print("Detected communities:")
    for i, comm in enumerate(communities):
        print(f"Community {i + 1}: {sorted(comm)}")

    # Run Leiden algorithm with higher gamma for finer communities
    print("\nRunning Leiden algorithm with gamma=1.5")
    communities = leiden(g, gamma=1.5)
    print("Detected communities:")
    for i, comm in enumerate(communities):
        print(f"Community {i + 1}: {sorted(comm)}")
