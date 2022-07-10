import numpy as np
from Bio import SeqIO, Align
from sklearn.cluster import AgglomerativeClustering

from module_painter.util import wrapping_substr


def refine_alignments(alns, fastas, min_pident=None, arc_eq_diffs=None, min_module_size=None):

    seq_data = {seq.id: seq.seq for fasta in fastas for seq in SeqIO.parse(fasta, "fasta")}
    # parent-wise simplification
    alns = (
        alns[alns.pident >= min_pident].copy()
        .groupby("sacc", as_index=False).apply(sync_boundaries, "sstart", arc_eq_diffs)
        .groupby("sacc", as_index=False).apply(sync_boundaries, "send", arc_eq_diffs)
        # .pipe(clean_boundaries)
        # .sort_values(by=["sstart", "send"])
        # .groupby(["sacc", "qacc"], sort=False, as_index=False).apply(remove_embedded)
        .groupby(["sacc", "qacc"], sort=False, as_index=False).apply(reverse_if_opposite)
        .groupby(["sacc", "qacc"], sort=False, as_index=False).apply(check_unaligned, seq_data, min_module_size)
    )
    # # global simplification
    # alns = (
    #     alns[~alns.merge_with_next]
    #     .groupby(["sacc", "sstart", "send"]).agg(dict(qacc=set, slen="first")) # merge equal intervals
    #     .reset_index()
    #     .sort_values(by=["sstart", "send"])
    #     .groupby("sacc", sort=False, as_index=False).apply(remove_embedded)
    # )

    return alns

def sync_boundaries(intervals, column, max_dist):
    '''
    Change intervals boundaries to make them similar
    '''
    if len(intervals) < 2:
        return

    model = AgglomerativeClustering(
        linkage='complete', affinity='l1',
        distance_threshold=max_dist, n_clusters=None
    ).fit(intervals[[column]])

    intervals[column] = intervals.groupby(model.labels_)[column].transform(min if "start" in column else max)

    return intervals

def clean_boundaries(intervals):
    full_intervals = intervals.send - intervals.sstart >= intervals.slen
    intervals.loc[full_intervals, "sstart"] = 0
    intervals.loc[full_intervals, "send"] = intervals[full_intervals].slen - 1
    return intervals

def remove_embedded(intervals):
    """
    Assumes intervals are (start, end) sorted
    - Case 1: [e(i), s(i)] C [s(i+1), e(i+1)]
    Checking the next one is enough
    s(i) = s(i+1) (since e(i) <= e(i+1))
    - Case 2: [e(i), s(i)] C [s(j), e(j)], j<i
    e(i) <= max(e(j), j<i) (and s(i) < s(j, j<i)
    The looping condition is handled separately
    """
    size = intervals.slen.iloc[0]
    max_end = np.maximum.accumulate(intervals.send).shift(1)
    embedded_in_prev = intervals.send <= max_end
    embedded_in_next = intervals.sstart == intervals.sstart.shift(-1)
    
    # looping condition
    loop_embedded = intervals.send < intervals.send.max() - size

    embedded = embedded_in_prev | embedded_in_next | loop_embedded

    return intervals[~embedded]

def check_unaligned(hits, seq_data, min_gap=30, max_gap=1e3, min_size_ratio=0.8, min_seq_id=0.7):
    if len(hits) < 2:
        return hits.assign(merge_with_next=False)
    # to align gap sequences later
    aligner = Align.PairwiseAligner(mode="global")

    # get gap boundaries
    (sacc, qacc) = hits[["sacc", "qacc"]].iloc[0]

    sseq = seq_data[sacc]
    qseq = seq_data[qacc]

    hits["next_q"] = get_next_hit(hits.qstart, len(qseq))
    hits["next_s"] = get_next_hit(hits.sstart, len(sseq))

    gap_cols = ["qend", "next_q", "send", "next_s"]

    scores = np.zeros(len(hits))

    # Align gaps on reference and query and fill gap
    for i, (gap_start_q, gap_end_q, gap_start_s, gap_end_s) in enumerate(hits[gap_cols].values):
        gaps = [gap_end_q-gap_start_q, gap_end_s-gap_start_s]
        gap_ratio = gaps[0] / gaps[1] # what matters is that the query is longer than the reference

        if gaps[1] <= min_gap:
            scores[i] = 1.1

        # if the gaps have the same size of reference and query, we check the identity
        elif gap_ratio > min_size_ratio and max(gaps) < max_gap:
            gap_q = wrapping_substr(qseq, gap_start_q, gap_end_q)
            gap_s = wrapping_substr(sseq, gap_start_s, gap_end_s)

            aln = aligner.align(gap_q, gap_s)
            scores[i] = aln.score / len(gap_s)

    hits["gap_score"] = scores
    hits["merge_with_next"] = scores > min_seq_id

    hits = merge_hits(hits)

    return hits

def get_next_hit(starts, size):
    """
    Get start position of the next hit
    The last start value wraps around (=starts[0]+size)
    """
    next_start = np.roll(starts.values, -1)
    next_start[-1] += size
    return next_start

def merge_hits(hits):
    indices = hits.index

    for i, merge in enumerate(hits.merge_with_next):
        if merge:
            current_idx = indices[i]
            next_idx = indices[(i+1) % len(hits)]

            if (i+1) == len(indices): # wrapping condition
                if hits.merge_with_next.all():
                    next_idx = hits.index[0]
                else:
                    next_idx = hits.index[~hits.merge_with_next][0]

            # left-extend next interval's start
            for which in ["q", "s"]:
                field = f"{which}start"
                hits.loc[next_idx, field] = hits.loc[current_idx, field]

                if (i+1) == len(indices): # wrapping condition
                    hits.loc[next_idx, f"{which}end"] += hits.loc[current_idx, f"{which}len"]

    # special case when everything is merged
    if hits.merge_with_next.all(): # we keep the first one
        hits.loc[indices[0], "merge_with_next"] = False

    return hits

def reverse_if_opposite(hits):
    most_likely_orient = hits.groupby("sstrand").nident.sum().idxmax()
    hits = hits[hits.sstrand == most_likely_orient]
    
    if most_likely_orient == '-':
        # Reverse the coordinates
        (hits.qstart, hits.qend) = (hits.qlen-hits.qend-1, hits.qlen-hits.qstart-1)

    return hits

def fill_gaps(interval_df, max_dist=None, cols=["qacc", "slen", "sstart", "send"]):
    """
    - Removes embedded intervals
    - Extend consecutive intervals if dist <= max_dist
    intervals: (id, start, end) list or np.array
    max_dist: int
    """
    if len(interval_df) <= 2:
        return interval_df[cols]

    intervals = interval_df[cols].to_numpy()
    
    results = []
    (prev, *intervals) = intervals

    while len(intervals) > 0:
        (cur, *intervals) = intervals
        (s1, e1) = prev
        (s2, e2) = cur

        if e1 >= e2: # cur C prev => keep prev
            pass
        elif s1 == s2 and e1 < e2: # prev C cur => keep cur
            prev = cur
        else:
            dist = s2 - e1
            if dist > max_dist:
                results.append(prev)
                prev = cur
            else: # combine prev and cur
                prev = (s1, e2)

    results.append(prev)

    results = pd.DataFrame(results, columns=cols)

    return results

 
