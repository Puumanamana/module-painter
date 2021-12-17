from itertools import product, combinations
from collections import Counter

import numpy as np
import pandas as pd
import igraph

from util import subset_fasta, get_minimap_graph


def group_breakpoints(paintings, fasta, min_len=10, min_id=0.8, min_cov=0.8,
                      verbose=False, threads=1, min_dist=100):

    # Find breakpoints in all genomes
    breaks = (pd.concat([painting.get_junctions() for painting in paintings])
              .reset_index(drop=True).reset_index())
    breaks['bin_id'] = -1

    # Handle differently ponctual and interval break points
    # ponctual case: based on common parents and position on genome
    breaks_ponct = breaks[breaks.end-breaks.start <= min_len]
    breaks.loc[breaks_ponct.index, 'bin_id'] = bin_punctual_bk(
        breaks_ponct, min_dist, min_len
    )
    # interval case: based on sequence similarity and coverage
    breaks_itv = breaks[breaks.end-breaks.start > min_len]
    breaks_itv_grouped = breaks_itv.groupby('target').agg(list)
    bk_fasta = subset_fasta(fasta, breaks_itv_grouped, min_len=min_len)
    bins_itv = get_minimap_graph(bk_fasta, min_id=min_id, min_cov=min_cov,
                                          verbose=0, threads=threads)
    breaks.loc[bins_itv.index, 'bin_id'] = bins_itv + breaks.bin_id.max() + 1

    # duplicate breakpoint for all possible combination of parent1/parent2
    breaks.parents = [list(product(x.split('/'), y.split('/')))
                      for (x,y) in breaks.parents.str.split(' <-> ')]
    breaks = breaks.drop(columns='index').explode('parents')

    # Remove combination with the same parent
    breaks = breaks[breaks.parents.map(set).map(len) == 2].reset_index(drop=True)

    return breaks

def bin_punctual_bk(breaks, min_dist, min_len):
    """
    Find groups of breakpoints with the same parents such that
    the distance between 2 successive breakpoints is smaller than min_dist bp

    Returns bin assignment of each breakpoints
    """

    bk_punct = breaks.groupby('parents').agg(dict(start=sorted, index=list))

    bin_assignment = pd.Series(-1, index=breaks.index)
    bin_id = 0
    for _, (starts, indices) in bk_punct.iterrows():

        if len(starts) == 1:
            bin_assignment[indices[0]] = bin_id
            bin_id += 1
            continue

        # Increment bin_id if the distance between 2 successive breakpoints
        # is greater than min_dist bp
        bin_assignment[indices] = bin_id
        bin_assignment[indices[1:]] += np.cumsum(np.diff(starts) > min_dist)

        bin_id = bin_assignment[indices[-1]] + 1

    return bin_assignment


def handle_missing_data(genomes, paintings, min_id=0.8, min_cov=0.8, threads=1):
    """
    Check if the missing data is approximately the same in different genomes
    """

    # Intervals with fillers
    nocovs = pd.DataFrame([[painting.target, arc.start, arc.end]
                           for painting in paintings
                           for arc in painting.arcs
                           if arc.meta == "NA"],
                          columns=['target', 'start', 'end'])

    nocovs_groups = nocovs.reset_index().groupby('target').agg(list)
    # Save to fasta for minimap alignment
    nocov_fasta = subset_fasta(genomes, nocovs_groups)

    # Build graph and extract connected bin_ids
    bin_ids = get_minimap_graph(nocov_fasta, min_id=min_id, min_cov=min_cov,
                                verbose=0, threads=threads)
    nocovs['bin_id'] = bin_ids.astype(str).str.replace(r'^(\d)', r'NoCov-\1', regex=True).sort_index()
    nocovs = nocovs.set_index(['target', 'start', 'end']).bin_id

    for painting in paintings:
        for arc in painting.arcs:
            name = (painting.target, arc.start, arc.end)
            if name in nocovs.index:
                arc.meta = nocovs.loc[name]
