from collections import defaultdict, deque
import copy

class Graph:
    def __init__(self):
        self.adj = defaultdict(set) # adjacency list {node: set(neighbors)}
        self.weights = defaultdict(int)  # For weighted edges, if needed
        self.total_weight = 0 # total weight of edges
        # supernode features
        self.internal_edges = defaultdict(int) # internal edges {community: internal edges}
        self.num_nodes = defaultdict(int) # number of nodes {community: number of nodes}

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
            
    def update_supernode_features(self, u, v, internal_edges_u, internal_edges_v, num_nodes_u, num_nodes_v):
        # Update supernode features
        self.internal_edges[u] = internal_edges_u
        self.internal_edges[v] = internal_edges_v
        self.num_nodes[u] = num_nodes_u
        self.num_nodes[v] = num_nodes_v
        
    def nodes(self):
        return list(self.adj.keys())

    def neighbors(self, u):
        return self.adj[u]
    
    def get_supernode_features(self, u):
        return self.internal_edges[u], self.num_nodes[u]
    
    def get_supernode_features_all(self):
        return self.internal_edges, self.num_nodes
    
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
    
    
def leiden(graph, gamma=1.0, resolution_parameter=None, seed=None, verbose=False):
    graph_ori = copy.deepcopy(graph) # store the original graph
    
    if resolution_parameter is not None:
        gamma = resolution_parameter

    # Initialize each node to its own community
    partition = {node: node for node in graph.nodes()}
    # Initialize community hierarchy: maps community ID to original nodes
    community_hierarchy = {node: {node} for node in graph.nodes()}
    # Set improvement to True to enter the while loop
    improvement = True
        
    level = 0
    
    while improvement:
        if verbose:
            print('\033[94m--------------------------------\033[0m')
            print('\033[94mLEVEL ', level + 1, '\033[0m')
            level += 1
            print('\033[94m--------------------------------\033[0m')
        else:
            level += 1
        
        improvement = False
        if verbose:
            print('*** Before local moving phase (node: community assignment): ', partition)
            print('    (The expected community assignment is: every node is its own community, e.g., 1:1, 2:2, 3:3, etc.)')
        
        
        # ======================================================
        # Local moving phase
        # ======================================================
        
        nodes = graph.nodes()
        for node in nodes:
            if verbose:
                print('Considering supernode: ', node)
            current_comm = partition[node]
            neighbor_comms = defaultdict(int)

            # Collect the weights to neighboring communities (on this graph)
            for neighbor in graph.neighbors(node):
                neighbor_comm = partition[neighbor]
                neighbor_comms[neighbor_comm] += graph.weights[(node, neighbor)]
                
            # neighbor_comms collects the weights to neighboring communities from this node (on this graph)
            # and yes, the graph is a supernode graph so the weights are the sum of weights of edges between nodes in the supernode
            if verbose:
                print(f'Neighbor supernodes of supernode {node} is {graph.neighbors(node)}, and their communities are {[partition[neighbor] for neighbor in graph.neighbors(node)]}')
            for comm, edge_weight in neighbor_comms.items():
                if verbose:
                    print('Connectedness to community ', comm, ': ', edge_weight)
                    
            
        
            # Compute the best community to move to
            best_comm = current_comm
            best_gain = 0
            
            # Get the supernode features of the current node
            community_m, community_n = graph.get_supernode_features_all() # community_m is the total of internal edge weights
            # community_n is the number of nodes in the supernode

            for comm, edge_weight in neighbor_comms.items():
                if verbose:
                    print('Considering moving supernode ', node, ' to community ', comm, ' with edge weight ', edge_weight)
                if comm == current_comm:
                    if verbose:
                        print('Supernode ', node, ' is already in community ', comm, '. Skip.')
                    continue
                
                # Delta CPM computation
                # 1. Number of edges in the current community
                # equals the sum of internal edges within the supernode
                # plus the number of edges from the current supernode to other supernodes in the current community
                
                
                
                # Compute k_i_in_current_comm (number of edges in the current community)
                # Internal edges within supernode i
                e_ii = graph.internal_edges[node]

                # Edges from supernode i to other supernodes in its current community (excluding itself)
                k_i_in_current_comm = e_ii
                for neighbor in graph.neighbors(node):
                    if partition[neighbor] == current_comm and neighbor != node:
                        k_i_in_current_comm += graph.weights[(node, neighbor)]

                
                # Compute k_i_in_new_comm
                k_i_in_new_comm = neighbor_comms[comm]
                
                # Compute n_i, n_current_comm, n_new_comm
                n_i = graph.num_nodes[node]
                n_current_comm = community_n[current_comm]
                n_new_comm = community_n[comm]


                if verbose:
                    print('*')
                    print('Total weight: ', graph.total_weight)
                    print(f'k_i_in_new_comm (connectedness to neighbor community {comm} / total edges from supernode {node} to community {comm}): ', k_i_in_new_comm)
                    print(f'k_i_in_current_comm (connectedness to current community {current_comm} / total edges from supernode {node} to other supernodes in current community {current_comm}): ', k_i_in_current_comm)
                    print(f'n_current_comm (number of nodes in current community {current_comm}): ', n_current_comm)
                    print(f'n_new_comm (number of nodes in new community {comm} before moving): ', n_new_comm)
                    print('*')
                
                
                n_i = graph.num_nodes[node]
                n_current_comm = community_n[current_comm]
                n_new_comm = community_n[comm]

                delta_e = k_i_in_new_comm - k_i_in_current_comm

                delta_r = gamma * (
                    ( (n_new_comm + n_i)*(n_new_comm + n_i - 1) - n_new_comm*(n_new_comm - 1)
                    + (n_current_comm - n_i)*(n_current_comm - n_i - 1) - n_current_comm*(n_current_comm - 1)
                    )
                ) / 2
                
                if verbose:
                    print('Delta E: ', delta_e)
                    print('Delta R: ', delta_r)

                delta_cpm = (delta_e - delta_r) / graph.total_weight

                
                if verbose:
                    print('Delta CPM: ', delta_cpm)
                    
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
        
        # ======================================================
        # End of local moving phase
        # ======================================================
        if improvement:
            # ======================================================
            # Refinement phase
            # ======================================================
            # Refinement phase to ensure communities are well-connected
            
            if verbose:
                print('')
                print('')
                print('*** After local moving phase, before refinement assignment (node: community assignment): ', partition)
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
            comm_nodes = defaultdict(list) # {community: list of nodes}
            
            # Compute inter community edges
            for comm, nodes_in_comm in communities.items():
                for node in nodes_in_comm:
                    for neighbor in graph.neighbors(node):
                        neighbor_comm = partition[neighbor]
                        if neighbor_comm != comm:
                            # To avoid double-counting, ensure consistent ordering
                            sorted_comm = tuple(sorted((new_comm_ids[comm], new_comm_ids[neighbor_comm])))
                            inter_comm_edges[sorted_comm] += graph.weights[(node, neighbor)]
                            
            # Compute community nodes
            for comm, nodes_in_comm in communities.items():
                count = 0
                for node in nodes_in_comm: # intrinsic number of nodes in the supernode
                    count += graph.num_nodes[node]
                comm_nodes[new_comm_ids[comm]] = count
                
            # Compute internal edges
            for comm, nodes_in_comm in communities.items(): # example: {2: {1, 3}, 3: {2}, 10: {4, 5, 6}, 11: {8}, 8: {9, 7}}
                for node in nodes_in_comm: # example: nodes_in_comm = {1, 3}, node = 1
                    intra_comm_edges[comm] += graph.internal_edges[node] # intrinsic edges inside the supernode
                    for neighbor in graph.neighbors(node):
                        if neighbor in nodes_in_comm:
                            # print('Updating comm ', comm)
                            intra_comm_edges[comm] += graph.weights[(node, neighbor)] / 2 # edges between supernodes that belong to the same community
                            
            # Reindexing intra_comm_edges
            intra_comm_edges_new = {new_comm_ids[comm]: intra_comm_edges[comm] for comm in intra_comm_edges}
            intra_comm_edges = copy.deepcopy(intra_comm_edges_new)
            
            # Update supernode features
            for comm, nodes_in_comm in communities.items():
                new_graph.update_supernode_features(new_comm_ids[comm], new_comm_ids[comm], intra_comm_edges[new_comm_ids[comm]], intra_comm_edges[new_comm_ids[comm]], comm_nodes[new_comm_ids[comm]], comm_nodes[new_comm_ids[comm]])

            # Add edges to the new graph based on inter-community edge weights
            for (comm1, comm2), weight in inter_comm_edges.items():
                if verbose:
                    print('Adding edge between ', comm1, ' and ', comm2, ' with weight ', weight)
                new_graph.add_edge(comm1, comm2, weight // 2)

            # Update total_weight for the new graph
            new_graph.total_weight = sum(new_graph.weights.values()) // 2  # Since undirected
            # new_graph.total_weight = graph_ori.total_weight # total weight is the same as the original graph 
            
            # Verify the total weight of the new graph by summing up the internal edges of all supernodes and the total weight of the new graph
            total_weight_new_graph = sum(new_graph.internal_edges.values()) + new_graph.total_weight
            new_graph.total_weight = total_weight_new_graph
            
            if verbose:
                print('Total weight of the new graph: ', total_weight_new_graph)
            
            # Update community hierarchy
            community_hierarchy = new_community_hierarchy

            # Prepare for next iteration
            graph = copy.deepcopy(new_graph)
            
            # Reinitialize partition: each new community is its own community
            partition = {node: node for node in graph.nodes()}
            
            if verbose:
                print('*** After aggregation (node: community assignment): ', partition)
                print('    Expect: every node is its own community, e.g., 1:1, 2:2, 3:3, etc.')
                print('\033[92m*** Community hierarchy updated ***\033[0m')
                print(community_hierarchy)
                print('Current graph: ', graph.adj)
                print('New graph internal edges: ', graph.internal_edges)
                print('New graph number of nodes: ', graph.num_nodes)
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

# Example Usage
if __name__ == "__main__":
    # Create a sample graph
    
    # g = Graph()
    # edges = [
    #     (1, 2), (1, 3), (2, 3),
    #     (4, 5), (5, 6), (4, 6),
    #     (3, 4),  # Connect the two communities
    #     (7, 8), (8, 9), (7, 9),
    #     (9, 4)   # Connect to the main community
    # ]
    # for u, v in edges:
    #     g.add_edge(u, v)
    #     g.update_supernode_features(u, v, 0, 0, 1, 1)
    
    from graph_gen import get_graph
    g = get_graph()
        
    # Run Leiden algorithm
    communities = leiden(g, gamma=0.5, verbose=True) # the higher the gamma, the more communities will be detected
    # the lower the gamma, the fewer communities will be detected (more changes will be accepted)

    print('Final communities: ', communities)