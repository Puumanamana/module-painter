from collections import Counter
from itertools import combinations

import numpy as np
import igraph


def count_recombinations(bk):
    rc = (
        bk
        .groupby(['par_set', 'target']).bin_id
        .agg(lambda x: len(x)//2)
        .sum(level='par_set')
    )

    return rc

def remove_alt_parents(bk, choice):
    '''
    Only keep `choice` in `bk` where multiple parents are possible
    '''
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
    - Iteratively choose parents with the most total recombinations
    - Solve equalities by choosing the parents that cover the most phages
    """

    multicov = True
    selected = {}
    # Sorted list of parents so that we count (A,B) + (A,B) the same way as (A,B) + (B, A)
    bk['par_set'] = bk.parents.apply(lambda x: tuple(sorted(x)))

    while multicov > 0:
        # Get parents with the most recombinations
        rc = count_recombinations(bk[~bk.par_set.isin(selected)])
        rc_max = rc.max()
        parent_pairs = rc.index[rc==rc_max]

        # Solve equalities by selecting the parents that cover the most phages overall
        (parents_max, score_max) = (parent_pairs[0], 0)

        if len(parent_pairs) > 1:
            parent_freq = Counter([p for pair in bk.parents for p in pair])

            for (p1, p2) in parent_pairs:
                score = ( # penalty on NoCov
                    parent_freq[p1] + ('NoCov' in p1) * -100 +
                    parent_freq[p2] + ('NoCov' in p2) * -100
                )
                if score > score_max:
                    score_max = score # score of the best parent
                    parents_max = (p1, p2) # parent pair with the max freq

        bk = remove_alt_parents(bk, parents_max)

        multicov = bk[['bin_id', 'target', 'i1', 'i2']].duplicated().sum()

        selected[parents_max] = rc_max

    return bk.drop(columns='par_set')

def cluster_breakpoints(bk, gamma=0.75):
    graph = igraph.Graph()
    graph.add_vertices(bk.target.unique())

    bk_groups = bk.groupby('bin_id').agg(list)

    for group, data in bk_groups.iterrows():
        for (t1, t2) in combinations(data.target, 2):
            graph.add_edge(t1, t2, parent=data.parents[0])

    communities = graph.community_leiden(
        objective_function="CPM",
        resolution_parameter=gamma,
        n_iterations=-1)

    vertices = np.array(graph.vs['name'])
    communities = [vertices[idx] for idx in communities]

    return communities

