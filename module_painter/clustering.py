import logging
from itertools import combinations

import numpy as np
import pandas as pd
import igraph

from module_painter.parent_selection import get_breakpoints, find_recombinations

logger = logging.getLogger("module-painter")

def get_links(breakpoints, feature="recombination", outdir=None):
    if "breakpoint" in feature:
        links = breakpoints.groupby("bk_id").agg(dict(sacc=list, parents="first"))
    else:
        links = find_recombinations(breakpoints)
        links = links.groupby("bk_ids").agg(dict(sacc=list, parents="first"))
        links.index = pd.MultiIndex.from_tuples(links.index)
        links.index.names = ["bk_id_1", "bk_id_2"]

    if outdir is not None:
        links_str = links.copy()
        links_str.sacc = links_str.sacc.apply("-".join)
        # save recombinations to file
        links_str.to_csv(f"{outdir}/links.csv")

    return links

def cluster_phages(links, outdir=None, method="connected_components", **kwargs):
    links.sacc = links.sacc.apply(lambda cl: list(combinations(cl, 2)))
    links = links.explode("sacc").dropna()
    
    recomb_graph = igraph.Graph()
    recomb_graph.add_vertices(links.sacc.explode().unique())
    recomb_graph.add_edges(links.sacc.tolist())

    if method == "connected_components":
        communities = recomb_graph.components()
    elif method == "leiden":
        communities = recomb_graph.community_leiden(n_iterations=-1, **kwargs)
    else:
        raise("Clustering method not yet suported")

    vertices = np.array(recomb_graph.vs['name'])
    communities = sorted([vertices[c] for c in communities], key=lambda x: -len(x))

    # Save clusters to file
    with open(f"{outdir}/clusters.csv", "w") as writer:
        writer.write("seq_id,cluster\n")
        for i, community in enumerate(communities):
            for c in community:
                writer.write(f"{c},{i}\n")

    # Log the clusters
    for i, c in enumerate(communities):
        logger.info(f"Cluster #{i}: {','.join(c)}")
                
    return communities
    
