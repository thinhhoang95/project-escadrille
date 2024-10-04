import networkx as nx
import matplotlib.pyplot as plt
import random

def get_sbm_graph():
    # Set random seed for reproducibility
    random.seed(42)

    # Define the number of communities and their sizes
    num_communities = 3
    community_sizes = [5, 5, 5]  # Four communities with different sizes

    # Define the probability of edges within and between communities
    p_intra = 1  # Probability of edges within the same community
    p_inter = 0.08  # Probability of edges between different communities

    # Create the Stochastic Block Model
    probs = [[p_intra if i == j else p_inter for j in range(num_communities)] for i in range(num_communities)]
    G = nx.stochastic_block_model(community_sizes, probs, seed=42)

    # Ensure no unconnected nodes
    isolated_nodes = list(nx.isolates(G))
    while isolated_nodes:
        for node in isolated_nodes:
            # Connect to a random non-isolated node
            connected_nodes = [n for n in G.nodes() if n not in isolated_nodes]
            if connected_nodes:
                target = random.choice(connected_nodes)
                G.add_edge(node, target)
        isolated_nodes = list(nx.isolates(G))
        
    return G

# Function to convert NetworkX graph to custom Graph class
from graph_cot import Graph, convert_networkx_to_custom_graph

def get_graph():
    G = get_sbm_graph()
    g = convert_networkx_to_custom_graph(G)
    return g