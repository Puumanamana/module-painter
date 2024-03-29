import random
from itertools import combinations

import pandas as pd
import numpy as np

random.seed(42)


def iter_prev_next(l):
    return zip(l, l[1:] + [l[0]])

def is_interfering(parents, freqs):
    """
    Interferences happen when we pick a breakpoint:
    - for which the parent pair appear an odd number of times
    on the reference
    - there are other parent pairs at the same location that 
    appear an even number of times
    """
    if len(parents) == 1:
        return False
    (ref, _) = parents.name
    freqs_for_bk = freqs.loc[ref].loc[parents].values
    is_rc = freqs_for_bk > 1
    is_odd = freqs_for_bk % 2 == 1

    return any(~is_odd) & is_odd & is_rc

def get_breakpoints(graphs):
    breakpoints = []

    for graph in graphs:
        for e in graph.es:
            breakpoints.append([graph["sacc"], e["bk_id"], e["parents"], e["pos"]])

    breakpoints = pd.DataFrame(breakpoints, columns=["sacc", "bk_id", "parents", "pos"])
    breakpoints["mult"] = breakpoints.duplicated(subset=["sacc", "pos"], keep=False)

    # Conditions for interference
    # the parent pair appears >1 time on the same reference with another parent pair with even number of breakpoints
    # 1) compute freqs = Counter[(sacc, parent_pair) for ...]
    # 2) groupby(["sacc", "bk_id"]).parents -> agg(lambda x: freqs[x]))
    freqs_per_sacc = breakpoints[["sacc", "parents"]].value_counts()
    breakpoints["interfere"] = breakpoints.groupby(["sacc", "bk_id"]).parents.transform(is_interfering, freqs_per_sacc)

    return breakpoints

def find_recombinations(bk):
    """
    Identify each possible recombination
    """
    def get_sorted_recombs(x):
        # avoid repeated breakpoints on the same child (due to repeats mostly)
        return [tuple(sorted(xi)) for xi in combinations(set(x), 2)]

    # Get all possible recombinations
    rc = bk.groupby(["sacc", "parents"]).bk_id.agg(get_sorted_recombs)
    # Reformat
    rc = rc.explode().dropna().reset_index().rename(columns=dict(bk_id="bk_ids"))

    mult = bk.loc[bk.mult].set_index(["sacc", "bk_id"]).sort_index().index
    rc["mult"] = [any((sacc, bk_id) in mult for bk_id in bk_ids)
                  for (sacc, bk_ids) in rc[['sacc', 'bk_ids']].values]

    interfere = bk.loc[bk.interfere].set_index(["sacc", "bk_id"]).sort_index().index
    rc["interference"] = [sum((sacc, bk_id) in interfere for bk_id in bk_ids)
                         for (sacc, bk_ids) in rc[['sacc', 'bk_ids']].values]
    return rc

def filter_recombinations(rc, bk):
    # Compute recombinations metrics
    scores = rc[rc.mult].groupby(["parents", "bk_ids"]).agg(dict(
        sacc=lambda x: len(set(x)),
        interference=sum
    )).rename(columns=dict(sacc="rc_prev"))

    # Filter based on: 1) prevalence 2) minimum interference
    scores = scores[scores.rc_prev==scores.rc_prev.max()]
    scores = scores[scores.interference==scores.interference.min()]

    # 3) on each breakpoint
    bk_abund = bk.groupby(["parents", "bk_id"]).sacc.agg(len).to_dict()
    scores["bk_abund"] = [sum(bk_abund[(parents, bk_id)] for bk_id in bk_ids)
                          for (parents, bk_ids) in scores.index]
    # Filter based on breakpoints
    # We pick the smallest bk abundance to not disrupt other potential recombinations
    scores = scores[scores.bk_abund==scores.bk_abund.max()]

    return scores

def filter_breakpoints(bk):
    scores = bk[bk.mult].groupby(["parents", "bk_id"]).sacc.agg(
        bk_prev=lambda x: len(set(x))
    )
    # Filter based on prevalence
    scores = scores[scores.bk_prev==scores.bk_prev.max()]

    return scores

def filter_parents(scores, bk):
    # Compute parent occurrence
    par_abund = bk.parents.str.split("/").explode().value_counts().to_dict()
    scores["par_abund"] = [sum(par_abund[parent] for parent in parents.split("/"))
                           for (parents, _) in scores.index]
    # Filter based on parent occurrence
    scores = scores[scores.par_abund==scores.par_abund.max()]

    return scores

def select_by_recombinations(graphs):
    """
    Minimize number of different parent pairs
    - Iteratively choose parents with the most total recombinations
    - Solve equalities by choosing the parents that cover the most phages
    """
    prev = ""
    while True:
        bk = get_breakpoints(graphs)
        rc = find_recombinations(bk)

        if not rc.mult.any():
            return

        # Scoring and filtering
        scores = filter_recombinations(rc, bk)
        scores = filter_parents(scores, bk)
        scores = scores.sample(1) # at random if equalities remain

        # Best recomb
        (parents, bk_ids) = scores.index[0]

        # Check for infinite loop
        if prev == scores.index[0]:
            print("Infinite loop...")
            print(rc[rc.bk_ids==prev[1]])
            import ipdb;ipdb.set_trace()
        prev = scores.index[0]
        
        # Remove from graph
        for graph in graphs:
            remove_alternative_breakpoints(graph, set(bk_ids), parents)
        
        # print(f"Selection {parents} {set(bk_ids)}. Remaining: {rc.mult.sum()}")

def select_by_breakpoints(graphs):
    prev = ""
    while True:
        bk = get_breakpoints(graphs)

        if not bk.mult.any():
            return

        # Scoring and filtering
        scores = filter_breakpoints(bk)
        scores = filter_parents(scores, bk)
        scores = scores.sample(1) # at random if equalities remain

        # Best recomb
        (parents, bk_id) = scores.index[0]

        # check for infinite loop
        if prev == scores.index[0]:
            print("Infinite loop...")
            import ipdb;ipdb.set_trace()
        prev = scores.index[0]
        
        # Remove from graph
        for graph in graphs:
            remove_alternative_breakpoints(graph, {bk_id}, parents)
        
        # print(f"Selection {parents} {bk_id}. Remaining: {bk.mult.sum()}")
            
def remove_alternative_breakpoints(graph, bk_, parents):
    if not graph.es: # no breakpoints
        return

    edges_with_bk = graph.es.select(bk_id_in=bk_)
    # make sure graph has both breakpoints and parents
    # > happens when the same bk_id happens on the same child
    # multiple times (most likely due to repeats)
    if not len(edges_with_bk.select(parents=parents)) >= len(bk_):
        return

    parents = set(parents.split("/"))

    # Remove from graph
    graph.delete_vertices(
        vi for e in edges_with_bk
        for vi in e.tuple
        if not graph.vs["parent"][vi] in parents
    )
