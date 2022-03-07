from itertools import combinations

import numpy as np
import igraph



def cluster_breakpoints(bk, gamma=0.75):
    graph = igraph.Graph()
    graph.add_vertices(bk.target.unique())

    bk_groups = bk.groupby('bin_id').agg(list)

    for group, data in bk_groups.iterrows():
        for (t1, t2) in combinations(data.target, 2):
            graph.add_edge(t1, t2, parent=data.parents[0])

    communities = graph.community_leiden(
        objective_function="CPM",
        resolution_parameter=gamma,
        n_iterations=-1)

    vertices = np.array(graph.vs['name'])
    communities = [vertices[idx] for idx in communities]

    return communities

