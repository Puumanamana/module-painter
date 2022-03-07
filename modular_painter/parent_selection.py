from itertools import combinations
import pandas as pd



def get_bk_info(graphs):
    bin_info = []

    for graph in graphs:
        for e in graph.es:
            parents = "/".join(sorted(graph.vs["parent"][i] for i in e.tuple))
            bin_info.append([e['ref'], e.index, e['bk_id'], parents])

    bin_info = pd.DataFrame(bin_info, columns=["ref", "eid", "bk_id", "parents"])
    bin_info.index = bin_info.ref + "-" + bin_info.eid.astype(str)

    return bin_info

def recombination_prevalence_score(bk):
    """
    Compute prevalence of each possible recombination at each position
    """
    def get_sorted_recombs(x):
        return [tuple(sorted(xi)) for xi in combinations(x, 2)]
    rc_prev = bk.groupby(["ref", "parents"]).bk_id.agg(get_sorted_recombs).explode().dropna().reset_index()
    rc_prev["rc_score"] = rc_prev.groupby("bk_id").ref.transform(lambda x: len(set(x))) # recombination prevalence
    rc_prev = rc_prev.explode("bk_id").drop_duplicates().set_index(["ref", "parents", "bk_id"])
    indices = list(zip(bk.bk_id, bk.parents, bk.bk_id))

    return rc_prev.reindex(index=indices).fillna(0).astype(int).values

def recombination_abundance_score(bk):
    """
    Compute #total of recombinations for each parent
    """
    rc_abund = bk.groupby(["ref", "parents"]).bk_id.agg(lambda x: len(x)//2).sum(level="parents")

    return rc_abund.loc[bk.parents].values

def breakpoint_prevalence_score(bk):
    """
    Compute prevalence of each breakpoint
    """
    bk_prev = bk.groupby(["bk_id", "parents"]).ref.agg(len)
    indices = list(zip(bk.bk_id, bk.parents))
    return bk_prev.loc[indices].values

def parents_abundance_score(bk):
    """
    Compute sum of abundance of each parent in breakpoint (across the whole dataset)
    """
    abund_scores = (
        bk.assign(parents=bk.parents.str.split("/"))
        .explode("parents")
        .groupby(["ref", "parents"]).eid.agg(len)
        .to_dict()
    )
    scores = [sum(abund_scores[ref, p] for p in parents.split("/"))
              for (ref, parents) in bk[["ref", "parents"]].values]
    return scores

def score_parents(bk):
    # Compute metrics for parent selection
    # 1) subset bk_ids such that
    #    a) mult
    #    b) if only one choice with recomb, pick it
    #    c) if 0 or more than 2, then we think
    # Maybe:
    # 1) identify all potential recombinations everywhere
    # 2) rank them in terms of #links created
    # 3) going from best to worst, resolve mult involving the recombination #i in ranked order.
    # Q: how do we update this list so that we can use edges ids? --> recombination info needs to list the edge ids
    # we need a dataframe for each recombination with: ref, parents, bk_id pair, edge_id pair
    # 1) we score each recombination by grouping by bk_id and looking at the length (=prevalence). The scores need to be recomputed each time we update the edges
    # 2) then we use a parent score: how many times does the parent pair occur overall in all recombinations
    # 3) then we use breakpoint score: how often does this breakpoint occur regardless of any recombinations
    # 4) then a phage score: how many times does the phage cover the genome where it appears

    # 1) Recombination scores
    scores = bk[["ref", "parents", "bk_id", "eid"]].assign(
        rc_score=recombination_prevalence_score(bk),
        pp_score=recombination_abundance_score(bk),
        bk_score=breakpoint_prevalence_score(bk),
        parent_score=parents_abundance_score(bk)
    )

    def choose(df, cols=None):
        (ref, bk_id) = df.iloc[0][["ref", "bk_id"]]
        i = 0
        while len(df) > 1 and i < len(cols):
            col = cols[i]
            best = df[col].max()
            df = df[df[col]==best]
            i += 1
        else:
            df = df.sample(1)
            print(f"Random choice for ref={ref} and bk_id={bk_id}")
        return df.iloc[0].drop(["ref", "bk_id"])

    # Issue with choose: we delete breakpoints again, without looking at the other, matching breakpoint.
    # They both need to be handled simulatenously: groupby(bk_av, bk_ap) -> choose
    
    scores = scores.groupby(["ref", "bk_id"]).apply(
        choose,
        cols=[f"{col}_score" for col in ["rc", "pp", "bk", "parent"]]
    )
    
    return scores

def select_parents(graphs):
    """
    Minimize number of different parent pairs
    - Iteratively choose parents with the most total recombinations
    - Solve equalities by choosing the parents that cover the most phages
    """

    bk_info = get_bk_info(graphs)

    scores = score_parents(bk_info)

    eids_kept = scores.groupby("ref").eid.agg(set).to_dict()

    for graph in graphs:
        ref = graph.vs[0]['ref']
        edges_to_rm = set(graph.es.indices) - eids_kept[ref]
        import ipdb;ipdb.set_trace()
        graph.delete_edges(edges_to_rm)
    
    # # Sorted list of parents so that we count (A,B) + (A,B) the same way as (A,B) + (B, A)
    # bk['par_set'] = bk.parents.apply(lambda x: tuple(sorted(x)))

    # prev = "" # to remove later, checks for infinite loops that can happen in rare cases
    # while True:
    #     # Update bin info with remaining edges
    #     remaining = [f"{e['ref']}-{e.index}" for g in graphs for e in g.es]
    #     bk_info = bk_info.loc[remaining]

    #     # Identify multiple edges coming from a vertex -- stop if none
    #     bk_info["mult"] = bk_info.groupby(["bk_id", "ref"]).parents.transform(
    #         lambda x: len(set(x)) > 1
    #     )

    #     if not bk_info.mult.any():
    #         break

    #     # Count, for each parent pair,
    #     # groupby(parent_pair, bin_id) --> how many genomes does this breakpoint affect
    #     # then we sum the top 2 breakpoints for each parent --> score
    #     parents_max = score_parents(bk_info)
        

        # # Get breakpoints with multiple options
        # bk["mult"] = bk[["bin_id", "i1", "i2", "target"]].duplicated(keep=False)

    #     if not bk.mult.any():
    #         break

    #     # Get parents with the most recombinations (and with multiple bk issues)
    #     rc = count_recombinations(bk)
    #     parent_pairs = rc.index[rc==rc.max()]

    #     # Solve equalities by selecting the parents that cover the most phages overall
    #     (parents_max, score_max) = (parent_pairs[0], 0)

    #     if len(parent_pairs) > 1:
    #         parent_freq = Counter([p for pair in bk.parents for p in pair])

    #         for (p1, p2) in parent_pairs:
    #             score = ( # penalty on NA
    #                 parent_freq[p1] + ('NA' in p1) * -100 +
    #                 parent_freq[p2] + ('NA' in p2) * -100
    #             )
    #             if score > score_max:
    #                 score_max = score # score of the best parent
    #                 parents_max = (p1, p2) # parent pair with the max freq

    #     print(f'selecting {parents_max}')
    #     # remove alternative parents where parents_max is present
    #     bin_ids_with_parent = set(bk.bin_id[(bk.par_set == parents_max) & bk.mult])
    #     is_bin_id = bk.bin_id.isin(bin_ids_with_parent) & bk.mult
        
    #     entries_to_rm = bk.index[bk.bin_id.isin(bin_ids_with_parent) &
    #                              bk.mult &
    #                              (bk.par_set != parents_max)].to_list()
    #     # remove pairs for the preceding breakpoint
    #     for bid in bin_ids_with_parent:
    #         prev_indices = bk# []
        
    #     # bk.drop(entries_to_rm, inplace=True)

    #     # to remove later, checks for infinite loops that can happen in rare cases
    #     if prev == parents_max:
    #         print(f"Infinite loop with parents: {prev}")
    #         print(bk[bk.bin_id.isin(bin_ids_with_parent)])
    #         import ipdb;ipdb.set_trace()
    #     prev = parents_max

    # return bk.drop(columns=['par_set', 'mult'])

