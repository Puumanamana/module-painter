from itertools import product, combinations
from collections import Counter

import numpy as np
import pandas as pd
import igraph

from util import subset_fasta, get_minimap_graph


def group_breakpoints(paintings, fasta, min_len=10, min_id=0.8, min_cov=0.8,
                      verbose=False, threads=1, min_dist=100):

    # Find breakpoints in all genomes
    breaks = (pd.concat([painting.get_junctions() for painting in paintings.values()])
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

def count_recombinations(bk):
    rc = (
        bk
        .groupby(['par_set', 'target']).bin_id
        .agg(lambda x: len(x)//2)
        .sum(level='par_set')
    )

    return rc

def remove_alt_parents(bk, choice):
    bin_ids = set(bk.loc[bk.par_set == choice, 'bin_id'])
    bk = (bk.set_index(['target', 'i1', 'i2']).sort_index()
          .assign(idt=np.arange(bk.shape[0])))

    to_rm = []

    for (keys, e) in bk.iterrows():
        (target, i1, i2) = keys

        if (
                e.par_set == choice # this is a correct entry
                or e.bin_id not in bin_ids # the bin_id is not affected by the change
                # this target does not have `choice` as a parent for this `bin_id`
                or all(bk[bk.bin_id == e.bin_id].loc[target, 'par_set'] != choice)
        ):
            continue

        # parents in the correct order
        ref = bk[(bk.bin_id == e.bin_id) & (bk.par_set == choice)].parents[0]
        
        to_rm.append(e.idt)

        n_itv = bk.loc[target].index.get_level_values(level=1).max()
        # fix the previous interval
        if i1 > 0:
            prev_i = (target, i1-1, i1)
        else:
            prev_i = (target, n_itv, i1)

        prev_to_rm = bk.loc[prev_i, 'idt'][
            bk.loc[prev_i, 'parents'].str[1] != ref[0]
        ].tolist()
        to_rm += prev_to_rm
            
        # fix the next interval
        if i2 < n_itv:
            next_i = (target, i2, i2+1)
        else:
            next_i = (target, i2, 0)

        next_to_rm = bk.loc[next_i, 'idt'][
            bk.loc[next_i, 'parents'].str[0] != ref[1]
        ].tolist()
        to_rm += next_to_rm

    bk = bk.reset_index().set_index('idt').drop(to_rm)

    return bk

def select_parents_by_sparsity(bk):
    """
    Minimize number of different parent pairs
    """

    ## WARNING: choice of breakpoints affects the next: A B/C A: if I choose AB then I choose BA after + careful order of breakpoints
    multicov = True
    selected = {}
    # Sorted list of parents, so that we count (A,B) + (A,B) the same way as (A,B) + (B, A)
    bk['par_set'] = bk.parents.apply(lambda x: tuple(sorted(x)))

    while multicov > 0:
        # Get parents with the most recombinations
        rc = count_recombinations(bk[~bk.par_set.isin(selected)])
        rc_max = rc.max()
        parent_pairs = rc.index[rc==rc_max]

        # Solve equalities by selecting the parents that cover the most phages overall
        p_max = parent_pairs[0] # parent pair with the max freq
        s_max = 0 # score of the best parent

        if len(parent_pairs) > 1:
            parent_freq = Counter([p for pair in bk.parents for p in pair])

            for (p1, p2) in parent_pairs:
                score = ( # penalty on NoCov
                    parent_freq[p1] + ('NoCov' in p1) * -100 +
                    parent_freq[p2] + ('NoCov' in p2) * -100
                )
                if score > s_max:
                    s_max = score
                    p_max = (p1, p2)

        bk = remove_alt_parents(bk, p_max)

        multicov = bk[['bin_id', 'target']].duplicated().sum()

        selected[p_max] = rc_max

    return bk.drop(columns='par_set')

def cluster_breakpoints(bk):
    graph = igraph.Graph()
    graph.add_vertices(bk.target.unique())

    bk_groups = bk.groupby('bin_id').agg(list)

    for group, data in bk_groups.iterrows():
        for (t1, t2) in combinations(data.target, 2):
            graph.add_edge(t1, t2, parent=data.parents[0])

    communities = graph.community_leiden(
        objective_function="CPM",
        resolution_parameter=0.5,
        n_iterations=-1)

    vertices = np.array(graph.vs['name'])
    communities = [vertices[idx] for idx in communities]

    return communities


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
    nocovs = pd.DataFrame([[target, arc.start, arc.end]
                           for (target, painting) in paintings.items()
                           for arc in painting.data
                           if 'NoCov' in arc.data.index],
                          columns=['target', 'start', 'end'])

    nocovs_groups = nocovs.reset_index().groupby('target').agg(list)
    # Save to fasta for minimap alignment
    nocov_fasta = subset_fasta(genomes, nocovs_groups)

    # Build graph and extract connected bin_ids
    bin_ids = get_minimap_graph(nocov_fasta, min_id=min_id, min_cov=min_cov,
                                verbose=0, threads=threads)
    nocovs['bin_id'] = bin_ids.astype(str).str.replace(r'^(\d)', r'NoCov-\1', regex=True).sort_index()
    nocovs = nocovs.set_index(['target', 'start', 'end']).bin_id

    for target, painting in paintings.items():
        for arc in painting.data:
            name = (target, arc.start, arc.end)
            if name in nocovs.index:
                arc.data.rename({'NoCov': nocovs.loc[name]}, inplace=True)
                # arc.data.index = pd.Index([nocovs.loc[name]] * len(arc.data), name='source')
