from itertools import product, combinations
from collections import Counter

import numpy as np
import pandas as pd
import igraph
from sklearn.cluster import AgglomerativeClustering

from util import subset_fasta, build_homology_graph


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
    breakpoint_data = breakpoint_data.groupby(["ref", "start", "end"]).eid.agg(lambda x: ";".join(map(str, x)))

    #===== Assign bin ids ====#
    bk_fasta = subset_fasta(fasta, breakpoint_data, min_len=min_overlap, outprefix=f"{outdir}/breakpoints")
    homology_graph = build_homology_graph(bk_fasta, min_id=0.95, verbose=0, threads=threads)

    #==== Extract connected components ====#
    bins = {}
    vertices_array = np.array(homology_graph.vs['name'])
    for i, component in enumerate(homology_graph.components()):
        for vertex in component:
            (ref, _, _, eids) = vertices_array[vertex].split("|")
            for eid in eids.split(";"):
                bins[(ref, eid)] = i
    
    bins = pd.Series(bins, name="bin").rename_axis(index=["ref","eid"]).sort_index()

    for graph in graphs:
        ref = graph.vs["ref"][0]
        graph.es["bk_id"] = bins.loc[ref].to_numpy()

def handle_missing_data(genomes, coverages, min_id=0.9, threads=1, outdir=None):
    """
    Check if the missing data is approximately the same in different genomes
    """

    # Intervals with fillers
    nocovs_df = pd.DataFrame([
        [coverage.ref, arc.start, arc.end]
        for coverage in coverages
        for arc in coverage.arcs if arc.meta == {"NA"}
    ], columns=["ref", "start", "end"])

    nocovs_df["uid"] = np.arange(len(nocovs_df))
    nocovs_df = nocovs_df.set_index(["ref", "start", "end"]).uid.sort_index()

    # Save to fasta for minimap alignment
    nocov_fasta = subset_fasta(genomes, nocovs_df, outprefix=f"{outdir}/missing_data")

    # Build homology graph between NA segments
    homology_graph = build_homology_graph(nocov_fasta, min_id=min_id,
                                          verbose=0, threads=threads)

    # Extract connected components
    bin_mapping = {}
    vertices_array = np.array(homology_graph.vs['name'])
    for i, component in enumerate(homology_graph.components()):
        for vertex in component:
            (name, start, end, _) = vertices_array[vertex].split("|")
            bin_mapping[(name, int(start), int(end))] = i

    # Update NAs with corresponding bin id
    for coverage in coverages:
        for arc in coverage.arcs:
            name = (coverage.ref, arc.start, arc.end)
            if name in bin_mapping:
                bin_id = bin_mapping[name]
                arc.meta = {f"NA-{bin_id}"}
