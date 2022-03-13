import random
from itertools import combinations, groupby
import pandas as pd
import numpy as np

random.seed(42)

def get_breakpoints(graphs):
    breakpoints = []

    for ref, graph in graphs.items():
        for e in graph.es:
            breakpoints.append([e.index, ref, e["bk_id"], e["parents"]])

    breakpoints = pd.DataFrame(breakpoints, columns=["eid", "ref", "bk_id", "parents"])
    breakpoints["mult"] = breakpoints.duplicated(subset=["ref", "bk_id"], keep=False)

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
    
    return rc

def filter_recombinations(rc, bk):
    # Compute recombinations metrics
    scores = rc[rc.mult].groupby(["parents", "bk_ids"]).ref.agg(
        rc_abund=len,
        rc_prev=lambda x: len(set(x))
    )
    # Filter based on: 1) prevalence 2) abundance
    scores = scores[scores.rc_prev==scores.rc_prev.max()]
    scores = scores[scores.rc_abund==scores.rc_abund.max()]

    # 3) on each breakpoint
    bk_abund = bk.groupby(["parents", "bk_id"]).ref.agg(len).to_dict()
    scores["bk_abund"] = [sum(bk_abund[(parents, bk_id)] for bk_id in bk_ids)
                          for (parents, bk_ids) in scores.index]
    # Filter based on breakpoints
    # We pick the smallest bk abundance to not disrupt other potential recombinations
    scores = scores[scores.bk_abund==scores.bk_abund.min()]

    return scores

def filter_breakpoints(bk):
    scores = bk[bk.mult].groupby(["parents", "bk_id"]).ref.agg(
        bk_abund=len,
        bk_prev=lambda x: len(set(x))
    )
    # Filter based on abundance
    scores = scores[scores.bk_abund==scores.bk_abund.max()]

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

        # Remove from graph
        for ref, graph in graphs.items():
            remove_alternative_breakpoints(graph, set(bk_ids), parents)
        
        print(f"Selection {parents} {set(bk_ids)}. Remaining: {rc.mult.sum()}")

def select_by_breakpoints(graphs):
    """
    """
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

        # Remove from graph
        for ref, graph in graphs.items():
            remove_alternative_breakpoints(graph, {bk_id}, parents)
        
        print(f"Selection {parents} {bk_id}. Remaining: {bk.mult.sum()}")
            
def remove_alternative_breakpoints(graph, bk_, parents):
    if len(graph.vs) < 2: # no breakpoints
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
    
def select_parents(graphs):

    select_by_recombinations(graphs)
    select_by_breakpoints(graphs)

    import ipdb;ipdb.set_trace()
