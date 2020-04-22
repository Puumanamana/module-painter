from itertools import product
import numpy as np
import pandas as pd

from util import subset_fasta, get_minimap_graph

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

    # Build graph and extract connected components
    components = get_minimap_graph(nocov_fasta, min_id=min_id, min_cov=min_cov, verbose=0, threads=threads)
    nocovs['components'] = components.astype(str).str.replace('^(\d)', r'NoCov-\1').sort_index()
    nocovs = nocovs.set_index(['target', 'start', 'end']).components

    for target, painting in paintings.items():
        for arc in painting.data:
            name = (target, arc.start, arc.end)
            if name in nocovs.index:
                arc.data.index = [nocovs.loc[name]] * len(arc.data)

def get_punctual_bk(breaks, min_dist, min_len):

    bk_punct = (breaks[breaks.end-breaks.start <= min_len]
                .reset_index()
                .groupby('parents').agg({'start': sorted, 'index': list}))

    components = np.zeros(len(breaks), dtype=int) - 1
    component_id = 0
    for parents, (intervals, indices) in bk_punct.iterrows():

        if len(intervals) == 1:
            components[indices[0]] = component_id
            component_id += 1
            continue

        groups = np.cumsum(np.diff(intervals) > min_dist)
        components[indices[0]] = component_id
        components[indices[1:]] = component_id + groups

        component_id = groups[-1] + component_id + 1

    components = pd.Series(dict(enumerate(components)))

    return components

def group_breakpoints(paintings, fasta, min_len=10, min_id=0.8, min_cov=0.8, verbose=False, threads=1, min_dist=100):

    breaks = (pd.concat([painting.get_junctions() for painting in paintings.values()])
              .reset_index(drop=True))

    components = get_punctual_bk(breaks, min_dist, min_len)

    breaks_grouped = breaks.reset_index().groupby(['target'])[["index", "start", "end"]].agg(list)
    
    bk_fasta = subset_fasta(fasta, breaks_grouped, min_len=min_len) 
    components_interv = get_minimap_graph(bk_fasta, min_id=min_id, min_cov=min_cov, verbose=0, threads=threads)
    components.loc[components_interv.index] = components_interv + components.max() + 1

    breaks['component'] = components

    breaks.parents = [list(product(x.split('/'), y.split('/')))
                      for (x,y) in breaks.parents.str.split(' <-> ')]

    breaks = breaks.explode('parents')
    breaks_grouped = breaks.groupby(['parents', 'component']).target

    scores = (breaks_grouped
              .agg(lambda x: '-'.join(x))
              .reset_index()
              .groupby(['parents', 'target']).component
              .agg(components=list, prevalence=len))

    scores['depth'] = pd.MultiIndex.get_level_values(scores.index, 'target').map(lambda x: 1+x.count('-'))

    scores.sort_values(by=['depth', 'prevalence'], ascending=False, inplace=True)

    # Then:
    # 1. Select rows until completely covered
    # 2. Resolve equalities for given target T between parents {(Xi, Yi)} with
    #    a) most frequent pair
    #    b) most frequent parent
    #    c) least frequent parent
    #    d) random choice
    import ipdb;ipdb.set_trace()
    return breaks

def cluster_breakpoints(breakpoints):
    import ipdb;ipdb.set_trace()
    return
