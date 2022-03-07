from itertools import groupby

import igraph
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import pdist, squareform

from modular_painter.arc import Arc
from modular_painter.coverage_util import get_filler_arc, get_successor


class Coverage:

    def __init__(self, *arcs, ref=None):
        self.arcs = list(arcs)
        self.ref = ref

        if arcs:
            self.size = self.arcs[0].size
        self.sort()

    def __len__(self):
        return sum(1 for arc in self.arcs if not arc.flagged)

    def get_last_arc(self):
        max_end = max(arc.end for arc in self.arcs if not arc.flagged)
        arc_max_end = [arc for arc in self.arcs if not arc.flagged and arc.end == max_end]

        min_start = min(arc.start for arc in arc_max_end)
        arc_max_end_min_start = next(arc for arc in arc_max_end if arc.start == min_start)

        return arc_max_end_min_start

    def iter_arcs(self, wrap=False):
        arcs = [arc for arc in self.arcs if not arc.flagged]
        if wrap:
            last_arc = self.get_last_arc()
            arcs = [last_arc] + arcs
        return iter(arcs)
    
    def iter_arcs_for_parent(self, parent):
        """
        Iter over arc from {parent}. Assumes self.arcs is (start, end) sorted.
        """
        arcs = [arc for arc in self.arcs if not arc.flagged and parent in arc.meta]
        return iter([arcs[-1]] + arcs)

    def __repr__(self, subset=None):
        if isinstance(subset, str):
            subset = set(subset)
        display = [f"Ref: {self.ref} (L={self.size})"]
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
        ref = covs[0].ref
        assert all(cov.size==size for cov in covs)
        assert all(cov.ref==ref for cov in covs)

        arcs = [arc for cov in covs for arc in cov.arcs]

        return cls(*arcs, ref=ref)

    @classmethod
    def from_pandas(cls, df, size):
        ref = df.sacc.iloc[0]
        assert all(df.sacc == ref)

        arcs = [Arc(row.sstart, row.send, size, {row.qacc}) for _, row in df.iterrows()]

        return cls(*arcs, ref=ref)

    def to_pandas(self):
        n = len(self)
        data = dict(
            start=np.zeros(n),
            end=np.zeros(n),
            parent=np.empty(n, dtype="<U32"),
            flag=np.zeros(n)
        )

        for i, arc in enumerate(self.arcs):
            data["start"][i] = arc.start
            data["end"][i] = arc.end
            data["parent"][i] = list(arc.meta)[0]
            data["flag"][i] = arc.flagged

        df = pd.DataFrame(data).assign(ref=self.ref)

        return df

    def remove_flagged(self):
        self.arcs = [arc for arc in self.arcs if not arc.flagged]

    def is_covered(self):
        """
        Assumes self.arcs is (start, end) sorted
        """
        arcs = self.iter_arcs(wrap=True)
        current_arc = next(arcs)

        for i, arc in enumerate(arcs):
            dist = current_arc.dist_to_next(arc)
            
            if dist > 0: # There is a hole
                print(f"{current_arc.end} < {arc.start} (dist={dist})")
                return False

            if arc.end > current_arc.end:
                current_arc = arc

        return True

    def sync_boundaries(self, max_dist, attrs):
        '''
        Change intervals boundaries to make them similar
        '''
        arcs = list(self.iter_arcs(wrap=False))
        boundaries = np.array([[getattr(arc, attr) for attr in attrs] for arc in arcs])

        # Wrapping around condition
        wrap = boundaries[:, -1] > self.size
        boundaries[wrap, -1] -= self.size

        model = AgglomerativeClustering(
            linkage='complete',
            distance_threshold=max_dist*len(attrs),
            affinity='l1',
            n_clusters=None
        ).fit(boundaries)

        cluster_info = pd.DataFrame({
            "indices": np.arange(len(arcs)),
            "cluster": model.labels_
        })
        for i, attr in enumerate(attrs):
            cluster_info[attr] = boundaries[:, i]

        # Compute each cluster's boundary: min of starts, max of ends
        agg_fn = dict(start=min, end=max)
        agg_fn = {attr: agg_fn[attr] for attr in attrs}
        agg_fn["indices"] = list
        
        cluster_info = (
            cluster_info
            .groupby("cluster").filter(lambda x: len(x) > 1)
            .groupby("cluster").agg(agg_fn)
        )

        # Sync boundaries of arcs in each cluster
        for v in cluster_info.values:
            bounds = v[:-1]
            indices = v[-1]
            if len(indices) < 2:
                continue
            for idx in indices:
                arc = arcs[idx]
                
                for (attr, bound) in zip(attrs, bounds):
                    # A bit more complicated than it should be because
                    # we need to handle the wrapping condition
                    current = getattr(arc, attr)
                    diff = bound - (current % self.size)
                    if diff != 0:
                        setattr(arc, attr, current+diff)
        self.sort()

    def extend_arcs_per_parent(self, max_extend_dist):
        """
        Extract all consecutive arcs for each parent
        - if embedded, flag the smaller one
        - if the distance is shorter than max_extend_dist, merge them
        """
        for parent in self.get_all_parents():
            arcs = self.iter_arcs_for_parent(parent)
            last_arc = next(arcs)
            for arc in arcs:
                if arc.flagged:
                    continue
                if last_arc.is_embedded(arc, strict=False):
                    last_arc.flag()
                else:
                    last_arc.try_fuse_with(arc, max_extend_dist)
                last_arc = arc
        self.sort()

    def merge_equal_intervals(self):
        if len(self) < 2:
            return
        arcs = self.iter_arcs(wrap=True)
        prev_arc = next(arcs)
        for arc in arcs:
            prev_arc.try_merge_with(arc)
            prev_arc = arc

    def fill_gaps(self, max_extend_dist):
        """
        Fill gaps (if any) between arcs
        """
        arcs = self.iter_arcs(wrap=True)
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
        arcs = self.iter_arcs(wrap=True)
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
        
    def get_junctions(self):
        '''
        Assumes perfect coverage
        '''

        if len(self) == 1:
            return pd.DataFrame([])

        data = np.zeros(len(self), dtype=[
            ('parents', '<U128'),
            # ('parents_set', '<U128'),
            ('start', 'uint32'),
            ('end', 'uint32'),
            ('i1', 'uint32'),
            ('i2', 'uint32')
        ])

        prev_i = len(self) - 1
        prev = self.arcs[prev_i]

        for i in range(len(self)):
            curr = self.arcs[i]

            (p1, p2) = ('/'.join(sorted(prev.meta)),
                        '/'.join(sorted(curr.meta)))

            data[i] = (
                "{} <-> {}".format(p1, p2),
                curr.start,
                prev.end % self.size,
                prev_i,
                i
            )
            prev_i = i
            prev = curr

        df = pd.DataFrame(data)
        df['ref'] = self.ref

        return df

    def get_overlap_graph(self, min_overlap=10):
        g = igraph.Graph()

        # add vertices
        for i, arc in enumerate(self.iter_arcs(wrap=False)):
            for meta in arc.meta:
                uid = f"{self.ref}.{meta}.{arc.start}-{arc.end}"
                g.add_vertex(name=uid, parent=meta, order=i, ref=self.ref,
                             start=arc.start, end=arc.end)

        # add edges
        arcs = self.iter_arcs(wrap=True)
        prev_arc = next(arcs)

        for arc in arcs:
            for meta1 in prev_arc.meta:
                uid1 = f"{self.ref}.{meta1}.{prev_arc.start}-{prev_arc.end}"
                start = arc.start
                end = prev_arc.end % self.size # wrapping condition. This should be the only arc with end >= size
                
                for meta2 in arc.meta:
                    uid2 = f"{self.ref}.{meta2}.{arc.start}-{arc.end}"
                    g.add_edge(uid1, uid2, ref=self.ref, start=start, end=end)
            prev_arc = arc

        return g
