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

    for ref, graph in graphs.items():
        for e in graph.es:
            breakpoints.append([ref, e["bk_id"], e["parents"]])

    breakpoints = pd.DataFrame(breakpoints, columns=["ref", "bk_id", "parents"])
    breakpoints["mult"] = breakpoints.duplicated(subset=["ref", "bk_id"], keep=False)

    # Conditions for interference
    # the bk appears >1 time on the same reference with another parent pair with even number of breakpoints
    # 1) compute freqs = Counter[(ref, parent_pair) for ...]
    # 2) groupby(["ref", "bk_id"]).parents -> agg(lambda x: freqs[x]))
    freqs_per_ref = breakpoints[["ref", "parents"]].value_counts()
    breakpoints["interfere"] = breakpoints.groupby(["ref", "bk_id"]).parents.transform(is_interfering, freqs_per_ref)

    return breakpoints

def find_recombinations(bk):
    """
    Identify each possible recombination
    """
    def get_sorted_recombs(x):
        return [tuple(sorted(xi)) for xi in combinations(x, 2)]

    # Get all possible recombinations
    rc = bk.groupby(["ref", "parents"]).bk_id.agg(get_sorted_recombs)
    # Reformat
    rc = rc.explode().dropna().reset_index().rename(columns=dict(bk_id="bk_ids"))

    mult = bk.loc[bk.mult].set_index(["ref", "bk_id"]).sort_index().index
    rc["mult"] = [any((ref, bk_id) in mult for bk_id in bk_ids)
                  for (ref, bk_ids) in rc[['ref', 'bk_ids']].values]

    interfere = bk.loc[bk.interfere].set_index(["ref", "bk_id"]).sort_index().index
    rc["interference"] = [sum((ref, bk_id) in interfere for bk_id in bk_ids)
                         for (ref, bk_ids) in rc[['ref', 'bk_ids']].values]
    return rc

def filter_recombinations(rc, bk):
    # Compute recombinations metrics
    scores = rc[rc.mult].groupby(["parents", "bk_ids"]).agg(dict(
        ref=lambda x: len(set(x)),
        interference=sum
    )).rename(columns=dict(ref="rc_prev"))

    # Filter based on: 1) prevalence 2) minimum interference
    scores = scores[scores.rc_prev==scores.rc_prev.max()]
    scores = scores[scores.interference==scores.interference.min()]

    # 3) on each breakpoint
    bk_abund = bk.groupby(["parents", "bk_id"]).ref.agg(len).to_dict()
    scores["bk_abund"] = [sum(bk_abund[(parents, bk_id)] for bk_id in bk_ids)
                          for (parents, bk_ids) in scores.index]
    # Filter based on breakpoints
    # We pick the smallest bk abundance to not disrupt other potential recombinations
    scores = scores[scores.bk_abund==scores.bk_abund.max()]

    return scores

def filter_breakpoints(bk):
    scores = bk[bk.mult].groupby(["parents", "bk_id"]).ref.agg(
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
            import ipdb;ipdb.set_trace()
        prev = scores.index[0]
        
        # Remove from graph
        for ref, graph in graphs.items():
            remove_alternative_breakpoints(graph, set(bk_ids), parents)
        
        # print(f"Selection {parents} {set(bk_ids)}. Remaining: {rc.mult.sum()}")

def select_by_breakpoints(graphs):
    """
    """
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
        for ref, graph in graphs.items():
            remove_alternative_breakpoints(graph, {bk_id}, parents)
        
        # print(f"Selection {parents} {bk_id}. Remaining: {bk.mult.sum()}")
            
def remove_alternative_breakpoints(graph, bk_, parents):
    if not graph.es: # no breakpoints
        return

    # make sure graph has all breakpoints in bk_
    if not len(graph.es.select(parents=parents, bk_id_in=bk_)) == len(bk_):
        return

    parents = set(parents.split("/"))

    # Remove from graph
    graph.delete_vertices(
        vi for e in graph.es.select(bk_id_in=bk_)
        for vi in e.tuple
        if not graph.vs["parent"][vi] in parents
    )
