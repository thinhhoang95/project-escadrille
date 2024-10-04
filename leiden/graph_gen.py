import networkx as nx
import matplotlib.pyplot as plt
import random

def get_sbm_graph(num_communities = 4, community_sizes = [10, 8, 6, 10], p_intra = 1, p_inter = 0.03):
    # Set random seed for reproducibility
    random.seed(42)

    # Define the number of communities and their sizes
    # num_communities = 4
    # community_sizes = [10, 8, 6, 10]  # Four communities with different sizes

    # Define the probability of edges within and between communities
    # p_intra = 1  # Probability of edges within the same community
    # p_inter = 0.03  # Probability of edges between different communities

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

def get_graph(num_communities = 4, community_sizes = [10, 8, 6, 10], p_intra = 1, p_inter = 0.03):
    G = get_sbm_graph(num_communities, community_sizes, p_intra, p_inter)
    g = convert_networkx_to_custom_graph(G)
    return g

from graph_cot import convert_custom_graph_to_networkx

if __name__ == "__main__":
    g = get_graph()
    # Visualize graph g with networkx 
    g_nx = convert_custom_graph_to_networkx(g)
    # Visualize the graph
    pos = nx.spring_layout(g_nx)
    plt.figure(figsize=(12, 8))
    
    # Assign colors to communities
    communities = nx.algorithms.community.label_propagation_communities(g_nx)
    color_map = plt.cm.get_cmap('viridis')
    colors = [color_map(i / len(communities)) for i, com in enumerate(communities) for _ in com]
    
    nx.draw(g_nx, pos, with_labels=True, node_color=colors, 
            node_size=500, font_size=10, font_weight='bold')
    
    plt.title("Stochastic Block Model Graph")
    plt.axis('off')
    plt.tight_layout()
    plt.show()