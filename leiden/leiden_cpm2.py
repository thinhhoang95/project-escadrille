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
    
def compute_community_metrics(graph, partition):
    community_m = defaultdict(int)  # Internal edges {community: internal edges}
    community_n = defaultdict(int)  # Number of nodes {community: number of nodes}

    for node in graph.nodes():
        community = partition[node]
        community_n[community] += 1
        for neighbor in graph.neighbors(node):
            if partition[neighbor] == community:
                community_m[community] += graph.weights[(node, neighbor)]

    # Since each internal edge is counted twice in undirected graphs
    for comm in community_m:
        community_m[comm] //= 2

    return dict(community_m), dict(community_n)
    
def update_community_metrics(community_m, community_n, move_buffer, verbose=False):
    if verbose:
        print('\033[94mUpdating community metrics...\033[0m')
    # move_buffer: (node, current_comm, best_comm, k_i_in_current_comm, k_i_in_new_comm, n_current_comm, n_new_comm)
    for (node, current_comm, best_comm, k_i_in_current_comm, k_i_in_new_comm, n_current_comm, n_new_comm) in move_buffer:
        if verbose:
            print('')
            print('Moving node ', node, ' from community ', current_comm, ' to community ', best_comm)
            print('Before moving: ')
            print('Community internal edges: ', community_m)
            print('Community number of nodes: ', community_n)
        # Decrease metrics of current_comm
        n_moved = n_current_comm
        community_m[current_comm] -= k_i_in_current_comm
        community_n[current_comm] -= n_moved

        # Increase metrics of best_comm
        community_m[best_comm] += k_i_in_new_comm
        community_n[best_comm] += n_moved
        
        if verbose:
            print('After moving: ')
            print('Community internal edges after moving (community: internal edges): ', community_m)
            print('Community number of nodes after moving (community: number of nodes): ', community_n)
        
    

def leiden(graph, gamma=1.0, resolution_parameter=None, seed=None, verbose=False):
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
    # Move buffer: for local moving phase, we move nodes from current community to neighboring communities
    # we write the changes into the buffer, so we can update the community metrics later
    move_buffer = [] # node, current community, new community
    # Set improvement to True to enter the while loop
    improvement = True
    
    
    
    
    # Initialize community metrics: m_c and n_c
    community_m = {node: 0 for node in graph.nodes()}  # Internal edges {community: internal edges}
    community_n = {node: 1 for node in graph.nodes()}  # Number of nodes {community: number of nodes}
    
    

    # Compute initial internal edges for each community
    for node in graph.nodes():
        for neighbor in graph.neighbors(node):
            if partition[neighbor] == partition[node]:
                community_m[partition[node]] += graph.weights[(node, neighbor)]
    # Since each internal edge is counted twice in undirected graphs
    for comm in community_m:
        community_m[comm] = community_m[comm] // 2
        
        
        
    level = 0
    while improvement:
        print('\033[94m--------------------------------\033[0m')
        print('\033[94mLEVEL ', level + 1, '\033[0m')
        level += 1
        print('\033[94m--------------------------------\033[0m')
        
        print('Community internal edges: ', community_m)
        print('Community number of nodes: ', community_n)
        
        improvement = False
        if verbose:
            print('*** Before local moving phase (node: community assignment): ', partition)
            print('The expected community assignment is: every node is its own community, e.g., 1:1, 2:2, 3:3, etc.')
        # ======================================================
        # Local moving phase
        # ======================================================
        nodes = graph.nodes()
        for node in nodes:
            if verbose:
                print('Considering node: ', node)
            current_comm = partition[node]
            neighbor_comms = defaultdict(int)

            # Collect the weights to neighboring communities
            for neighbor in graph.neighbors(node):
                neighbor_comm = partition[neighbor]
                neighbor_comms[neighbor_comm] += graph.weights[(node, neighbor)]
                


            for comm, edge_weight in neighbor_comms.items():
                if verbose:
                    print('Connectedness to community ', comm, ': ', edge_weight)
            
            # Compute the best community to move to
            best_comm = current_comm
            best_gain = 0

            for comm, edge_weight in neighbor_comms.items():
                if verbose:
                    print('Considering moving node ', node, ' to community ', comm, ' with edge weight ', edge_weight)
                if comm == current_comm:
                    if verbose:
                        print('Node ', node, ' is already in community ', comm, '. Skip.')
                    continue
                # Compute Delta CPM using the correct formula
                # ΔCPM = (2 * k_i_in_new_comm / m) - gamma * (n_new_comm + n_current_comm -1) / m
                k_i_in_new_comm = edge_weight 
                # k_i_in_current_comm = 0
                # Since node is currently in current_comm, we need to find k_i_in_current_comm
                # which is the number of edges node has within current_comm
                # k_i_in_current_comm = neighbor_comms.get(current_comm, 0)
                k_i_in_current_comm = community_m[current_comm]

                n_current_comm = community_n[current_comm]
                n_new_comm = community_n[comm]


                if verbose:
                    print('*')
                    print('Total weight: ', graph.total_weight)
                    print(f'k_i_in_new_comm (connectedness to neighbor community {comm} / total edges from node {node} to community {comm}): ', k_i_in_new_comm)
                    print(f'k_i_in_current_comm (connectedness to current community {current_comm} / total edges from node {node} to other nodes in current community {current_comm}): ', k_i_in_current_comm)
                    print(f'n_current_comm (number of nodes in current community {current_comm}): ', n_current_comm)
                    print(f'n_new_comm (number of nodes in new community {comm} before moving): ', n_new_comm)
                    # print('Delta CPM 1: ', (2 * k_i_in_new_comm) / graph.total_weight)
                    # print('Delta CPM 2: ', gamma * (n_new_comm + n_current_comm - 1) / graph.total_weight)
                    print(f'Delta CPM 1: {(k_i_in_new_comm - k_i_in_current_comm) / graph.total_weight}')
                    print(f'Delta CPM 2: {gamma * (n_new_comm - n_current_comm + 1) / graph.total_weight}')
                    # print(f'Delta CPM 2: {n_new_comm * n_current_comm / graph.total_weight}')
                    print('*')
                
                
                # delta_cpm = (2 * k_i_in_new_comm) / graph.total_weight - gamma * (n_new_comm - n_current_comm + 1) / graph.total_weight
                delta_cpm = (k_i_in_new_comm - k_i_in_current_comm) / graph.total_weight - gamma * (n_new_comm - n_current_comm + 1) / graph.total_weight
                
                if verbose:
                    print('Delta CPM: ', delta_cpm, ', best gain: ', best_gain)
                
                
                if delta_cpm > best_gain:
                    best_gain = delta_cpm
                    best_comm = comm
                else:
                    if verbose:
                        print('Change could NOT be accepted. Keep node ', node, ' in community ', current_comm)
                
            if best_comm != current_comm and best_gain > 0:
                # Move node to best_comm
                partition[node] = best_comm
                if verbose:
                    print('\033[91mNode ', node, ' is moved to community ', best_comm, '\033[0m')
                improvement = True
                
        
        # Check if any improvement was made during local moving and refinement
        # If yes, proceed to aggregation
        if improvement:
            # ======================================================
            # Refinement phase
            # ======================================================
            # Refinement phase to ensure communities are well-connected
            print('')
            print('')
            if verbose:
                print('*** After local moving phase, before refinement (node: community assignment): ', partition)
            partition = refine(graph, partition)
            if verbose:
                print('*** After refinement, before aggregation (node: community assignment): ', partition)
            
            
            # ======================================================
            # Aggregation phase
            # ======================================================
            # Build communities based on the current partition
            communities = defaultdict(list)
            for node, comm in partition.items():
                communities[comm].append(node)
                
            # Reset community metrics
            community_m = {} # recall: community_m is the sum of internal edges for each community
            community_n = {} # recall: community_n is the number of nodes for each community
            for comm, nodes_in_comm in communities.items():
                community_n[comm] = len(nodes_in_comm)
                internal_edges = 0
                for node in nodes_in_comm:
                    for neighbor in graph.neighbors(node):
                        if partition[neighbor] == comm:
                            internal_edges += graph.weights[(node, neighbor)]
                # Each internal edge is counted twice
                community_m[comm] = internal_edges // 2


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
            inter_comm_edges = defaultdict(int) # {community pair: edge weight}
            intra_comm_edges = defaultdict(int) # {community: edge weight}
            
            # Compute inter-community edges
            for comm, nodes_in_comm in communities.items():
                for node in nodes_in_comm:
                    for neighbor in graph.neighbors(node):
                        neighbor_comm = partition[neighbor]
                        if neighbor_comm != comm:
                            # To avoid double-counting, ensure consistent ordering
                            sorted_comm = tuple(sorted((new_comm_ids[comm], new_comm_ids[neighbor_comm])))
                            inter_comm_edges[sorted_comm] += graph.weights[(node, neighbor)]
                        else:
                            intra_comm_edges[comm] += graph.weights[(node, neighbor)]

            # Add edges to the new graph based on inter-community edge weights
            for (comm1, comm2), weight in inter_comm_edges.items():
                if verbose:
                    print('Adding edge between ', comm1, ' and ', comm2, ' with weight ', weight)
                new_graph.add_edge(comm1, comm2, weight // 2)

            # Update total_weight for the new graph
            # new_graph.total_weight = sum(new_graph.weights.values()) // 2  # Since undirected
            new_graph.total_weight = graph.total_weight # total weight is the same as the original graph 
            
            # Update community hierarchy
            community_hierarchy = new_community_hierarchy

            # Prepare for next iteration
            graph = new_graph
            
            # Reinitialize partition: each new community is its own community
            partition = {node: node for node in graph.nodes()}
            
            # Reindexing community metrics
            community_m_copy = copy.deepcopy(community_m)
            community_n_copy = copy.deepcopy(community_n)
            community_m = {node: 0 for node in graph.nodes()}  # Internal edges {community: internal edges}
            community_n = {node: 1 for node in graph.nodes()}  # Number of nodes {community: number of nodes}
            for comm in community_m_copy: # comm is the old community ID, e.g., 2, 9, 6, 8; new_comm_ids[comm] is the corresponding new community ID, e.g., 0, 1, 2, 3
                community_m[new_comm_ids[comm]] = community_m_copy[comm]
                community_n[new_comm_ids[comm]] = community_n_copy[comm]
            
            if verbose:
                print('*** After aggregation (node: community assignment): ', partition)
                print('Expect: every node is its own community, e.g., 1:1, 2:2, 3:3, etc.')
                print('Community internal edges: ', community_m)
                print('Community number of nodes: ', community_n)
                print('\033[92m*** Community hierarchy updated ***\033[0m')
                print(community_hierarchy)
                print('Current graph: ', graph.adj)
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

    
def check_if_graph_is_connected(graph):
    """
    Checks if the given graph is connected.
    
    :param graph: Graph object
    :return: Boolean indicating whether the graph is connected
    """
    if not graph.nodes():
        return True  # An empty graph is considered connected

    start_node = next(iter(graph.nodes()))
    visited = set()
    queue = deque([start_node])

    while queue:
        node = queue.popleft()
        if node not in visited:
            visited.add(node)
            queue.extend(neighbor for neighbor in graph.neighbors(node) if neighbor not in visited)

    return len(visited) == len(graph.nodes())

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
        
        
    # Check if the graph is connected
    graph_is_connected = check_if_graph_is_connected(g)
    if not graph_is_connected:
        raise ValueError('The graph is not connected. Please check the input graph.')

    # Run Leiden algorithm
    communities = leiden(g, gamma=0.5, verbose=True) # the higher the gamma, the more communities will be detected
    # the lower the gamma, the fewer communities will be detected (more changes will be accepted)

    # Print communities
    print("Detected communities:")
    for i, comm in enumerate(communities):
        print(f"Community {i + 1}: {comm}")
