# centrality.py
import networkx as nx

def compute_all_centralities(G):
    """
    Compute common centrality metrics on a NetworkX graph G.
    Returns a dictionary of centrality results keyed by metric name.
    """
    centralities = {}

    try:
        centralities["degree"] = nx.degree_centrality(G)
        centralities["betweenness"] = nx.betweenness_centrality(G)
        centralities["closeness"] = nx.closeness_centrality(G)
        centralities["eigenvector"] = nx.eigenvector_centrality(G, max_iter=1000)
    except nx.NetworkXException as e:
        print(f"Centrality computation error: {e}")

    return centralities

def get_top_n_nodes(centrality_scores, n=5):
    """
    Return top-N nodes from a centrality score dictionary.
    """
    sorted_nodes = sorted(centrality_scores.items(), key=lambda x: x[1], reverse=True)
    return sorted_nodes[:n]

def normalize_scores(score_dict):
    """
    Normalize values in a centrality dictionary to [0, 1] range for visualization.
    """
    if not score_dict:
        return {}
    max_val = max(score_dict.values())
    if max_val == 0:
        return {k: 0 for k in score_dict}
    return {k: v / max_val for k, v in score_dict.items()}

def color_and_scale_nodes(cmd, top_nodes, color="red", base_radius=0.3, scale=1.5):
    """
    Color and scale top nodes in PyMOL based on normalized centrality scores.
    Assumes node IDs match residue indices.
    """
    for node, score in top_nodes:
        scaled_radius = base_radius + (score * scale)
        selection = f"name CA and resi {node}"
        cmd.set("sphere_scale", scaled_radius, selection)
        cmd.color(color, selection)
