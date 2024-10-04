import networkx as nx
from collections import defaultdict

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
    
def convert_networkx_to_custom_graph(nx_graph):
    """
    Converts a NetworkX graph to a custom Graph class instance.

    Parameters:
    - nx_graph (networkx.Graph): The NetworkX graph to convert.

    Returns:
    - Graph: An instance of the custom Graph class with the same structure.
    """
    custom_graph = Graph()

    # Add all edges to the custom graph
    for u, v, data in nx_graph.edges(data=True):
        # If the NetworkX graph has weights, use them; otherwise, default to 1
        weight = data.get('weight', 1)
        custom_graph.add_edge(u, v, weight)
        custom_graph.update_supernode_features(u, v, 0, 0, 1, 1)

    # Optionally, ensure all nodes are added (even isolated nodes)
    for node in nx_graph.nodes():
        if node not in custom_graph.adj:
            custom_graph.adj[node]  # This will initialize the set for the node

    return custom_graph

def convert_custom_graph_to_networkx(custom_graph):
    """
    Converts a custom Graph class instance to a NetworkX graph.

    Parameters:
    - custom_graph (Graph): The custom Graph instance to convert.

    Returns:
    - networkx.Graph: A NetworkX graph with the same structure and weights.
    """
    nx_graph = nx.Graph()

    # Add all nodes and edges to the NetworkX graph
    for u in custom_graph.nodes():
        nx_graph.add_node(u)
        for v in custom_graph.neighbors(u):
            if u < v:  # To avoid adding the same edge twice
                weight = custom_graph.weights[(u, v)]
                nx_graph.add_edge(u, v, weight=weight)

    # Add supernode features as node attributes
    internal_edges, num_nodes = custom_graph.get_supernode_features_all()
    nx.set_node_attributes(nx_graph, internal_edges, 'internal_edges')
    nx.set_node_attributes(nx_graph, num_nodes, 'num_nodes')

    return nx_graph

import matplotlib.pyplot as plt

def visualize_custom_graph(g, communities = None): # g: custom graph
    # Visualize graph g with networkx 
    g_nx = convert_custom_graph_to_networkx(g)
    # Visualize the graph
    pos = nx.spring_layout(g_nx)
    plt.figure(figsize=(6, 4))
    
    # Assign colors to communities
    if communities:
        node_map = {node_id_in_graph: j for j, node_id_in_graph in enumerate(g.nodes())} # reindexing the nodes in the graph to zero-indexed
        color_map = plt.cm.get_cmap('viridis')
        node_colors = {}
        for i, community in enumerate(communities):
            for node in community:
                node_colors[node_map[node]] = color_map(i / len(communities)) # i: community index
        colors = [node_colors.get(node_map[node], 'lightgrey') for node in g_nx.nodes()]
    else:
        colors = 'lightblue'  # Default color if no communities are provided
    
    nx.draw(g_nx, pos, with_labels=True, node_color=colors, 
            node_size=500, font_size=10, font_weight='bold')
    
    plt.title("Stochastic Block Model Graph")
    plt.axis('off')
    plt.tight_layout()
    plt.show()
