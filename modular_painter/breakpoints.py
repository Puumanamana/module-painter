from itertools import product, combinations
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
    breaks = breaks[breaks.parents.map(set).map(len) == 2]

    return breaks


def find_recombinations(bk):
    return (
        bk
        .groupby(['parents', 'target']).bin_id
        .agg(lambda x: len(x)//2)
        .sum(level='parents')
    )

def select_parents_2(bk):

    parent_stats = (
        bk.groupby(['parents', 'target']).bin_id
        .agg(bin_id=tuple, prevalence=len)
    )
    while True:
        n_recomb = find_recombinations(bk)
        best_pair = n_recomb.index[n_recomb==n_recomb.max()]

        import ipdb;ipdb.set_trace()

def select_parents(bk):

    select_parents_2(bk)
    # parent selection: use the parent pair with the highest prevalence
    scores = (bk.groupby(['parents', 'bin_id']).target.agg(tuple)
              .reset_index()
              .groupby(['parents', 'target']).bin_id
              .agg(bin_id=tuple, prevalence=len))

    scores['depth'] = (scores.index
                       .get_level_values('target')
                       .map(lambda x: 1+len(x)))

    scores['known'] = (2-scores.index
                       .get_level_values('parents')
                       .map(lambda x: ''.join(x).count('NoCov')))

    scores.sort_values(by=['depth', 'prevalence', 'known'], ascending=False, inplace=True)

    rc = find_recombinations(bk).sort_values(ascending=False)
    import ipdb;ipdb.set_trace()
    # Make a parent pair choice for all breakpoints
    choice = (scores.reset_index().explode('bin_id').explode('target')
              .groupby(['bin_id', 'target'], as_index=False).parents.first()
              .drop_duplicates()
              .set_index(['bin_id', 'target']))

    bk.parents = choice.reindex(index=list(zip(bk.bin_id, bk.target))).to_numpy()

    return bk.drop_duplicates()

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
                arc.data.index = [nocovs.loc[name]] * len(arc.data)
