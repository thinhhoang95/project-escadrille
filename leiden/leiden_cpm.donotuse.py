from collections import defaultdict, deque
import copy

class Graph:
    def __init__(self):
        self.adj = defaultdict(set)
        self.weights = defaultdict(int)  # For weighted edges, if needed
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
    Leiden algorithm implementation with CPM.
    :param graph: Graph object
    :param gamma: CPM resolution parameter
    :return: list of communities
    """
    if resolution_parameter is not None:
        gamma = resolution_parameter

    # Initialize each node to its own community
    partition = {node: node for node in graph.nodes()}
     # Initialize community hierarchy: maps community ID to original nodes
    community_hierarchy = {node: {node} for node in graph.nodes()}
    # Set improvement to True to enter the while loop
    improvement = True

    while improvement:
        print('--------------------------------')
        improvement = False
        print('*** Before local moving phase (node: community assignment): ', partition)
        # Local moving phase
        nodes = graph.nodes()
        for node in nodes:
            print('Considering node: ', node)
            current_comm = partition[node]
            neighbor_comms = defaultdict(int)

            # Collect the weights to neighboring communities
            for neighbor in graph.neighbors(node):
                neighbor_comm = partition[neighbor]
                neighbor_comms[neighbor_comm] += graph.weights[(node, neighbor)]
                
            print('Neighbor communities (community: weight): ', neighbor_comms)

            # Compute the best community to move to
            best_comm = current_comm
            best_gain = 0
            total_weight = graph.total_weight

            for comm, edge_weight in neighbor_comms.items():
                print('Considering moving node ', node, ' to community ', comm, ' with edge weight ', edge_weight)
                if comm == current_comm:
                    continue
                # Calculate the gain in CPM
                # ΔCPM = (edge_weight / total_weight) - gamma * (n_c_new / (2 * total_weight))
                # For accurate CPM, we need to track community sizes and internal edges
                # For simplicity, we'll use a heuristic similar to modularity
                delta = (edge_weight / total_weight) - gamma
                print('Edge weight: ', edge_weight, ', total weight: ', total_weight)
                print('Delta: ', delta, ', best gain: ', best_gain)
                if delta > best_gain:
                    print('Change could be accepted. Might move node ', node, ' to community ', comm, ' with gain ', delta)
                    best_gain = delta
                    best_comm = comm
                else:
                    print('Change could NOT be accepted. Keep node ', node, ' in community ', current_comm)
                
            if best_comm != current_comm:
                print('Finally, node ', node, ' is assigned to community ', best_comm)
                partition[node] = best_comm
                improvement = True

        # Refinement phase to ensure communities are well-connected
        print('*** After local moving phase, before refinement (node: community assignment): ', partition)
        partition = refine(graph, partition)
        print('*** After refinement, before aggregation (node: community assignment): ', partition)

        # ======================================================
        # Aggregation phase
        # ======================================================
        # Check if any improvement was made during local moving and refinement
        # If yes, proceed to aggregation
        if improvement:
            # Build communities based on the current partition
            communities = defaultdict(list)
            for node, comm in partition.items():
                communities[comm].append(node)

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
                print('Adding edge between ', comm1, ' and ', comm2, ' with weight ', weight)
                new_graph.add_edge(comm1, comm2, weight)

            # Update total_weight for the new graph
            new_graph.total_weight = sum(new_graph.weights.values()) // 4  # Since undirected

            # Update community hierarchy
            community_hierarchy = new_community_hierarchy

            # Prepare for next iteration
            graph = new_graph
            # Reinitialize partition: each new community is its own community
            partition = {node: node for node in graph.nodes()}
            print('*** After aggregation (node: community assignment): ', partition)
            print('*** Community hierarchy updated ***')
            print(community_hierarchy)
        # ======================================================
        # End of Aggregation phase
        # ======================================================
        
    # Final partition mapping back to original nodes
    final_communities = defaultdict(set)
    if not partition:
        final_communities = community_hierarchy
    else:
        for node, comm in partition.items():
            # Each node here is an aggregated community; map to original nodes
            original_nodes = community_hierarchy[node]
            final_communities[comm].update(original_nodes)

    # Convert sets to lists for readability
    final_communities = [list(members) for members in final_communities.values()]
    return final_communities

def refine(graph, partition):
    """
    Refinement step to ensure communities are well-connected.
    Splits communities if they are not internally connected.
    """
    if not partition:
        return {}

    communities = defaultdict(list)
    for node, comm in partition.items():
        communities[comm].append(node)

    refined_partition = {}
    new_comm_id = max(partition.values(), default=0) + 1  # Start new community IDs from here

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

    # Run Leiden algorithm
    communities = leiden(g, gamma=0.09)

    # Print communities
    print("Detected communities:")
    for i, comm in enumerate(communities):
        print(f"Community {i + 1}: {comm}")
