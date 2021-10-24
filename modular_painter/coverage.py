from itertools import groupby

import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import pdist, squareform

from modular_painter.arc import Arc
from modular_painter.coverage_util import get_filler_arc, get_successor


class Coverage:

    def __init__(self, *arcs, target=None):
        self.arcs = list(arcs)
        self.target = target

        if arcs:
            self.size = self.arcs[0].size

    def __len__(self):
        return sum(1 for arc in self.arcs if not arc.flagged)

    def get_last_arc(self):
        max_end = max(arc.end for arc in self.arcs if not arc.flagged)
        arc_max_end = [arc for arc in self.arcs if not arc.flagged and arc.end == max_end]

        min_start = min(arc.start for arc in arc_max_end)
        arc_max_end_min_start = next(arc for arc in arc_max_end if arc.start == min_start)

        return arc_max_end_min_start

    def iter_arcs(self):
        arcs = [arc for arc in self.arcs if not arc.flagged]
        last_arc = self.get_last_arc()
        return iter([last_arc] + arcs)
    
    def iter_arcs_parent(self, parent):
        arcs = [arc for arc in self.arcs if not arc.flagged and parent in arc.meta]
        return iter([arcs[-1]] + arcs)

    def __repr__(self, subset=None):
        display = [f"Target: {self.target} (L={self.size})"]
        for i, arc in enumerate(self.arcs):
            if arc.flagged:
                continue
            if subset is not None:
                if not arc.meta.intersection(subset):
                    continue
            display.append(f"{i:3} - {arc.start:>9}{arc.end:>9}: {arc.meta}")
        return "\n".join(display)

    def get_all_parents(self):
        parents = set()
        for arc in self.arcs:
            if not arc.flagged:
                parents |= arc.meta
        return parents
    
    def sort(self):
        self.arcs = sorted(self.arcs, key=lambda arc: (arc.start, arc.end))

    @classmethod
    def from_coverages(cls, *covs):
        size = covs[0].size
        target = covs[0].target
        assert all(cov.size==size for cov in covs)
        assert all(cov.target==target for cov in covs)

        arcs = [arc for cov in covs for arc in cov.arcs]

        return cls(*arcs, target=target)

    def to_pandas(self):
        n = len(self)
        data = dict(
            start=np.zeros(n),
            end=np.zeros(n),
            parent=np.zeros(n),
            flag=np.zeros(n)
        )

        for i, arc in enumerate(self.arcs):
            data["start"][i] = arc.start
            data["end"][i] = arc.end
            data["parent"][i] = arc.meta
            data["flag"][i] = arc.flagged

        df = pd.DataFrame(data).assign(target=self.target)

        return df

    def is_covered(self):
        arcs = self.iter_arcs()
        current_end = next(arcs).end - self.size

        for i, arc in enumerate(arcs):
            if arc.start > current_end: # There is a hole
                print(f"{arc.start} < {current_end}")
                return False

            if arc.end > current_end: # There is a hole
                current_end = arc.end

        return True

    def sync_boundaries(self, max_dist, which='start'):
        '''
        Change intervals boundaries to make them similar
        '''
        arcs = list(self.iter_arcs())
        boundary = np.array([getattr(arc, which) for arc in arcs])

        model = AgglomerativeClustering(
            linkage='complete',
            distance_threshold=max_dist,
            affinity='l1',
            n_clusters=None
        ).fit(boundary.reshape(-1, 1))

        cluster_info = pd.DataFrame({
            "indices": np.arange(len(arcs)),
            "cluster": model.labels_,
            which: boundary,
        })

        agg_fn = dict(start=min, end=max)
        agg_fn = {"indices": list, which: agg_fn[which]}
        
        cluster_info = (
            cluster_info
            .groupby("cluster").filter(lambda x: len(x) > 1)
            .groupby("cluster").agg(agg_fn)
        )

        for (indices, bound) in cluster_info.values:
            if len(indices) < 2:
                continue
            for idx in indices:
                setattr(arcs[idx], which, bound)

        self.sort()

    def extend_arcs_per_parent(self, max_extend_dist):
        """
        Extract all consecutive arcs for each parent
        - if embedded, flag the smaller one
        - if the distance is shorter than max_extend_dist, merge them
        """
        for parent in self.get_all_parents():
            arcs = self.iter_arcs_parent(parent)
            last_arc = next(arcs)
            for arc in arcs:
                if last_arc.is_embedded(arc):
                    last_arc.flag()
                else:
                    last_arc.try_fuse_with(arc, max_extend_dist)
                last_arc = arc
        self.sort()

    def merge_equal_intervals(self):
        arcs = self.iter_arcs()
        prev_arc = next(arcs)
        for arc in arcs:
            prev_arc.try_merge_with(arc)
            prev_arc = arc

    def fill_gaps(self, max_extend_dist):
        """
        Fill gaps (if any) between arcs
        """
        arcs = self.iter_arcs()
        prev_arc = next(arcs)
        fillers = []

        for arc in arcs:
            dist = prev_arc.dist_to_next(arc)

            if not arc.is_embedded(prev_arc):
                if dist > max_extend_dist:
                    # add a filler arc
                    filler = get_filler_arc(prev_arc, arc, dist)
                    fillers.append(filler)
                else:
                    prev_arc.try_extend_end_with(arc, max_extend_dist, dist=dist)

                prev_arc = arc
        self.arcs += fillers
        self.sort()

    def simplify_embedded(self):
        arcs = self.iter_arcs()
        prev_arc = next(arcs)

        for arc in arcs:
            if arc.is_embedded(prev_arc):
                arc.flag()
            else:
                if prev_arc.is_embedded(arc):
                    prev_arc.flag()
                prev_arc = arc

    def get_minimal_coverage(self):
        """
        Lee and lee 1983 algorithm
        On a circle-cover minimization problem
        Conditions:
        - self is fully covered (no gaps)
        - No embedded intervals
        """

        arcs = [arc for arc in self.arcs if not arc.flagged]
        n_interv = len(self)

        if n_interv == 1:
            return self

        successors = [get_successor(i, arcs) for (i, arc) in enumerate(arcs)]

        # Order of the B_i
        t_ = np.zeros(n_interv, dtype=int)

        (i, B_1) = (1, 0)

        S = [arcs[B_1]]
        start0, current_end = arcs[B_1].bounds()

        for arc in arcs:
            arc.marked = False

        arcs[B_1].marked = True # Mark the 1st arc
        t_[B_1] = 1
        B_i = B_1

        while self.size-current_end+(start0-1) > 0:
            B_i = successors[B_i]
            S.append(arcs[B_i])
            arcs[B_i].marked = True

            current_end = arcs[B_i].end

            t_[B_i] = i + 1
            i += 1

        # At this point we have an initial cover
        k = i

        is_opt = False

        while not is_opt:
            i += 1
            B_i = successors[B_i]

            S.append(arcs[B_i])

            if(arcs[B_i].marked or
               ((i % (k-1) == 1) # A total of k disjoint zones has been found
                and arcs[B_i].dist_to_next(arcs[B_1]) > 0)
            ):
                is_opt = True
                if i == t_[B_i] + (k-1): # Current arc Bi is marked and is the same as B[i-k]]
                    indices = range(i-(k-1), i)
                else:
                    indices = range(k)
            else: # B_i is in S or the number of disjoint zones is k, size of the minimal cover=k
                arcs[B_i].marked = True
                t_[B_i] = i

        for arc in self.arcs:
            arc.flag()

        for i in indices:
            S[i].unflag()

        self.sort()
