from itertools import combinations

import numpy as np
import igraph

from modular_painter.parent_selection import get_breakpoints, find_recombinations


def cluster_phages(graphs, gamma=0.75, feature="breakpoints"):
    bk = get_breakpoints(graphs)

    recomb_graph = igraph.Graph()
    recomb_graph.add_vertices(list(graphs.keys()))

    if feature == "breakpoint":
        edges = bk.groupby(["bk_id", "parents"]).ref.agg(list)
    else:
        rc = find_recombinations(bk)
        edges = rc.groupby(["bk_ids", "parents"]).ref.agg(list)

    for refs in edges.values:
        for (r1, r2) in combinations(refs, 2):
            recomb_graph.add_edge(r1, r2)

    communities = recomb_graph.community_leiden(
        objective_function="CPM",
        resolution_parameter=gamma,
        n_iterations=-1)

    vertices = np.array(recomb_graph.vs['name'])
    communities = sorted([vertices[idx] for idx in communities], key=lambda x: -len(x))

    return communities
    
