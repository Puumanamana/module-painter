from itertools import combinations

import numpy as np
import igraph

from module_painter.parent_selection import get_breakpoints, find_recombinations


def cluster_phages(graphs, gamma=0.75, feature="breakpoints", outdir=None):
    bk = get_breakpoints(graphs)
    rc = find_recombinations(bk)

    # save recombinations to file
    rc_summary = rc[["ref", "parents"]].value_counts().unstack(fill_value=0)
    rc_summary.to_csv(f"{outdir}/recombinations_summary.csv")

    # Make feature graph before clustering 
    recomb_graph = igraph.Graph()
    recomb_graph.add_vertices([graph["ref"] for graph in graphs])

    if feature == "breakpoint":
        edges = bk.groupby(["bk_id", "parents"]).ref.agg(list)
    else:
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

    # Save clusters to file
    with open(f"{outdir}/clusters.csv", "w") as writer:
        writer.write("seq_id,cluster\n")
        for i, community in enumerate(communities):
            for c in community:
                writer.write(f"{c},{i}\n")

    return communities
    
