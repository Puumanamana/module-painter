import logging
import re
from itertools import combinations

import numpy as np
import pandas as pd
import igraph

from module_painter.parent_selection import get_breakpoints, find_recombinations
from module_painter.display import get_cmap

logger = logging.getLogger("module-painter")

def summarize_recombinations(breakpoints, outdir):
    rc = find_recombinations(breakpoints)
    rc = (rc.groupby(["bk_ids", "parents"]).agg(dict(sacc=list))
          .droplevel("bk_ids")
          .sort_values(by="sacc", key=lambda x: x.str.len(), ascending=False))
    rc.sacc = rc.sacc.apply("-".join)
    rc.to_csv(f"{outdir}/recombinations.csv")

    return rc

def get_links(df, feature="recombination", outdir=None):
    
    keys = ["bk_id", "parents"]

    if "recombination" in feature:
        keys[0] = "bk_ids"
        df = find_recombinations(df)

    links = (df.groupby(keys).sacc.agg(lambda x: list(combinations(x, 2)))
             .explode().dropna()
             .droplevel(keys[0]).reset_index()
             .set_index("sacc").parents)

    if not links.empty:
        links.index = pd.MultiIndex.from_tuples(links.index)
        links.sort_index(inplace=True)

    return links

def cluster_phages(links, outdir=None, method="leiden", group_pattern=None, resolution=0.2):

    if links.empty: return []

    recomb_graph = igraph.Graph()
    recomb_graph.add_vertices(list({x for tup in links.index for x in tup}))
    edge_weights = links.index.value_counts()
    recomb_graph.add_edges(edge_weights.index, attributes=dict(weight=edge_weights.values))

    vnames = recomb_graph.vs["name"]
    aes = dict(vertex_label=vnames, edge_width=recomb_graph.es["weight"])
               # vertex_size=30, vertex_label_size=20)
    
    if group_pattern is not None:
        groups = [re.findall(group_pattern, x)[0] for x in vnames]
        aes["vertex_label"] = [re.sub(group_pattern, "", vname) for vname in vnames]
        cmap = get_cmap(groups)
        aes["vertex_color"] = [cmap[g] for g in groups]
    try:
        igraph.plot(recomb_graph, target=f"{outdir}/recombination_graph.pdf", **aes)
    except AttributeError:
        print("Skipping graph plot: neither packages 'pycairo' or 'cairocffi' are installed")

    if method == "connected_components":
        communities = recomb_graph.components()
    elif method == "leiden":
        communities = recomb_graph.community_leiden(n_iterations=-1, resolution_parameter=resolution, weights="weight",
                                                    objective_function="CPM")
    elif method == "infomap":
        communities = recomb_graph.community_infomap(edge_weights="weight")
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
    
