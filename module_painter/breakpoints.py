from itertools import product, combinations
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import igraph
from sklearn.cluster import AgglomerativeClustering

from module_painter.util import subset_fasta, build_homology_graph


def set_breakpoint_ids(graphs, fasta, min_overlap=50, max_dist=100, outdir=None, threads=1):
    """
    Extract breakpoints in all genomes and maps them together
    """
    # Fetch and format breakpoint data
    breakpoint_data = pd.concat([
        pd.DataFrame({attr: graph.es[attr] for attr in graph.es.attributes()})
        for graph in graphs
    ])
    breakpoint_data["eid"] = [eid for graph in graphs for eid in graph.es.indices]
    breakpoint_data = breakpoint_data.groupby(["sacc", "sstart", "send"]).eid.agg(
        lambda x: ";".join(map(str, x))
    )

    #===== Assign bin ids ====#
    bk_fasta = subset_fasta(fasta, breakpoint_data, min_len=min_overlap, outprefix=f"{outdir}/breakpoints")
    homology_graph = build_homology_graph(bk_fasta, min_id=0.95, min_size_ratio=0.5, verbose=0, threads=threads)

    #==== Extract connected components ====#
    bins = {}
    vertices_array = np.array(homology_graph.vs['name'])
    for i, component in enumerate(homology_graph.components()):
        for vertex in component:
            (sacc, _, _, eids) = vertices_array[vertex].split("^")
            for eid in eids.split(";"):
                bins[(sacc, int(eid))] = i

    bins = pd.Series(bins, name="bin").rename_axis(index=["sacc","eid"]).sort_index()

    for graph in graphs:
        if graph["sacc"] in bins.index.get_level_values("sacc"):
            graph.es["bk_id"] = bins.loc[graph["sacc"]].to_numpy()

def map_missing_parents(genomes, coverages, min_id=0.97, threads=1, outdir=None):
    """
    Check if the missing data is approximately the same in different genomes
    """
    # Intervals with fillers
    nocovs_df = pd.DataFrame([
        [coverage.sacc, arc.sstart, arc.send]
        for coverage in coverages
        for arc in coverage.arcs if arc.qacc == {"NA"}
    ], columns=["sacc", "sstart", "send"])

    if len(nocovs_df) == 0:
        return

    nocovs_df["uid"] = np.arange(len(nocovs_df))
    nocovs_df = nocovs_df.set_index(["sacc", "sstart", "send"]).uid.sort_index()

    # Save to fasta for minimap alignment
    nocov_fasta = subset_fasta(genomes, nocovs_df, outprefix=f"{outdir}/missing_data")

    # Build homology graph between NA segments
    homology_graph = build_homology_graph(
        nocov_fasta, min_id=min_id, min_size_ratio=0.3,
        verbose=0, threads=threads
    )

    # Extract connected components
    bin_mapping = {}
    vertices_array = np.array(homology_graph.vs['name'])
    for i, component in enumerate(homology_graph.components()):
        for vertex in component:
            (name, sstart, send, _) = vertices_array[vertex].split("^")
            bin_mapping[(name, int(sstart), int(send))] = i

    # Update NAs with corresponding bin id
    for coverage in coverages:
        found = defaultdict(int) # in case repeats happen
        for arc in coverage.arcs:
            name = (coverage.sacc, arc.sstart, arc.send)
            if name in bin_mapping:
                bin_id = bin_mapping[name]
                suffix = ""
                if bin_id in found:
                    suffix = f"-repeat{found[bin_id]}"
                arc.qacc = {f"NA-{bin_id}{suffix}"}
                found[bin_id] += 1
