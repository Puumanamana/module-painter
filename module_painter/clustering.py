import logging
import re
from itertools import combinations

import numpy as np
import pandas as pd
import igraph
from bokeh.palettes import Colorblind, Turbo256

from module_painter.parent_selection import get_breakpoints, find_recombinations

logger = logging.getLogger("module-painter")

def get_colors(groups):
    categories = set(groups)
    n = len(categories)
    if n > 256:
        print("Too many groups to display. Skipping coloring")
        return []
    elif n <= 8:
        colors = Colorblind[max(3, n)]
    else:
        colors = Turbo256

    cmap = dict(zip(categories, colors))
        
    return [cmap[g] for g in groups]
    

def get_links(breakpoints, feature="recombination", outdir=None):
    if "breakpoint" in feature:
        links = breakpoints.groupby("bk_id").agg(dict(sacc=list, parents="first"))
    else:
        links = find_recombinations(breakpoints)
        links = links.groupby("bk_ids").agg(dict(sacc=list, parents="first"))

    if outdir is not None:
        links_str = links.copy().sort_values(by="sacc", key=lambda x: x.str.len(), ascending=False)
        links_str.sacc = links_str.sacc.apply("-".join)
        # save recombinations to file
        links_str.to_csv(f"{outdir}/recombinations.csv")

    return links

def cluster_phages(links, outdir=None, method="leiden", group_pattern=None, resolution=0.2):

    links = links.copy()[links.sacc.apply(len) > 1]
    
    if links.empty: return []

    links.sacc = links.sacc.apply(lambda cl: list(combinations(cl, 2)))
    links = links.explode("sacc")
    
    recomb_graph = igraph.Graph()
    recomb_graph.add_vertices(links.sacc.explode().unique())
    recomb_graph.add_edges(links.sacc.tolist())

    vnames = recomb_graph.vs["name"]
    aes = dict(vertex_label=vnames)
    
    if group_pattern is not None:
        groups = [re.findall(group_pattern, x)[0] for x in vnames]
        aes["vertex_label"] = [re.sub(group_pattern, "", vname) for vname in vnames]
        colors = get_colors(groups)
        if colors:
            aes["vertex_color"] = colors
    try:
        igraph.plot(recomb_graph, target=f"{outdir}/recombination_graph.pdf", **aes)
    except AttributeError:
        print("Skipping graph plot: neither packages 'pycairo' or 'cairocffi' are installed")

    if method == "connected_components":
        communities = recomb_graph.components()
    elif method == "leiden":
        communities = recomb_graph.community_leiden(n_iterations=-1, resolution_parameter=resolution)
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
    
