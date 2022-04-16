"""
Important
---------
> If a list of intervals L is ordered on start and its intervals have 
> no containment relation, then it is also sorted on stop.

1) Sticky/Synchronized boundaries
2) Extend interval from the same parent
3) Remove embedded
"""

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering


def apply_rules(alignments, min_itv_dist=15, min_module_size=30):
    # ==== Sticky boundaries ==== #
    for attr in ["sstart", "send"]:
        alignments = sync_boundaries(alignments, attr, min_itv_dist)
    alignments.sort_values(by=["qacc", "sstart", "send"], inplace=True)

    # ==== Simplify coverage of each parent separately ==== #
    # Remove embedded modules from the same parent
    # important: sort=False
    embedded = alignments.groupby("qacc", sort=False, as_index=False).apply(is_embedded)
    alignments = alignments[~embedded.values]
    # Merge short-distance consecutive modules
    # At this point, alignments for the same parent are increasing
    # for both start and end (i.e. no embedding)
    alignments = (alignments.groupby("qacc", sort=False, as_index=False)
                  .apply(fuse_close_modules, min_module_size)
                  .sort_values(by=["sstart", "send"])
                  .reset_index(drop=True))
    # # ==== All parents together ==== #
    # Remove embedded
    embedded = is_embedded(alignments)
    alignments = alignments[~embedded.values]
    # Extend interval to touch or add NA
    alignments = fill_gaps(alignments, min_module_size).sort_values(by=["sstart", "send"])

    print(alignments)
    import ipdb;ipdb.set_trace()

    return alignments

def fill_gaps(intervals, min_module_size):
    """
    Assumes:
    - intervals are (start, end) sorted
    - no embedding
    - same parents are already fused
    """
    size = intervals.slen.iloc[0]
    next_start = np.roll(intervals.sstart, -1)
    next_start[-1] += size

    fillers = []
    for i, (start, end) in enumerate(zip(intervals.send, next_start)):
        gap = end-start
        if gap <= 0:
            continue
        elif gap > min_module_size:
            fillers.append(["NA", start, end, size])
        else:
            if i+1 == len(intervals):
                intervals.iloc[0, 1] = 0
                intervals.iloc[-1, 2] = size-1
            else:
                intervals.iloc[i, 2] += gap//2
                intervals.iloc[i+1, 1] += (gap-gap//2)

    fillers = pd.DataFrame(fillers, columns=["qacc", "sstart", "send", "slen"])
    intervals = pd.concat([intervals, fillers])

    return intervals

def is_embedded(intervals):
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

    return embedded_in_prev | embedded_in_next | loop_embedded

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

def compute_gaps(intervals):
    prev_end = np.roll(intervals.send, 1)
    prev_end[0] -= intervals.slen.iloc[0]

    return intervals.sstart - prev_end
    
def merge_modules(intervals_df, max_dist=None):
    """
    - Removes embedded intervals
    - Extend consecutive intervals if dist <= max_dist
    - Align with NW if dist > max_dist
    intervals_df: (id, qstart, qend, qlen, sstart, send, slen) list or array
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
    results["slen"] = interval_df.slen.iloc[0]
    results = try_circularize(results, max_dist)

    return results

