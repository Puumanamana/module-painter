from itertools import groupby

import igraph
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from Bio import Align

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

    def __repr__(self, subset=None):
        if isinstance(subset, str):
            subset = set(subset)
        display = [f"Ref: {self.ref} (L={self.size})"]
        for i, arc in enumerate(self.arcs):
            if subset is not None:
                if not arc.meta.intersection(subset):
                    continue
            display.append(f"{i:3} - {arc.start:>9}{arc.end:>9}: {arc.meta}")
        return "\n".join(display)

    def sort(self):
        self.arcs = sorted(self.arcs, key=lambda arc: (arc.start, arc.end))

    def remove_flagged(self):
        self.arcs = [arc for arc in self.arcs if not arc.flagged]

    def show(self):
        parents = []
        for i, arc in enumerate(self.arcs):
            parents.append("|".join(sorted(arc.meta)))# + f"({arc.end-arc.start:,})")
        parents = "->".join(f"{p:^5}" for p in parents)
        return f"(Reference: {self.ref})  >>>  {parents}"

    def iter_arcs(self, wrap=False, parent=None):
        """
        Important: for wrap=True, self.arcs needs to be (start, end) sorted
        """
        if parent is not None:
            arcs = [arc for arc in self.arcs if not arc.flagged and parent in arc.meta]
        else:
            arcs = [arc for arc in self.arcs if not arc.flagged]

        if wrap and len(arcs) > 1:
            max_end = max(arc.end for arc in arcs)
            last_arc = next(arc for arc in arcs if arc.end == max_end)            
            arcs = [last_arc] + arcs

        return iter(arcs)
    
    def get_all_parents(self):
        parents = set()
        for arc in self.iter_arcs(wrap=False):
            parents |= arc.meta
        return parents
    
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

        arcs = [Arc(row.sstart, row.send, size, meta={row.qacc}, qstart=row.qstart, qend=row.qend)
                for _, row in df.iterrows()]

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

            if not arc.is_embedded(current_arc):
                current_arc = arc

        return True

    def sync_boundaries(self, attr, max_dist):
        """
        Change intervals boundaries to make them similar
        """
        if len(self) < 2:
            return

        features = pd.Series([getattr(arc, attr) for arc in self.iter_arcs(wrap=False)])

        model = AgglomerativeClustering(
            linkage="complete", affinity="l1",
            distance_threshold=max_dist, n_clusters=None
        ).fit(features.values[:, None])

        data = features.groupby(model.labels_).transform(min if "start" in attr else max)
        data.groupby(level=0).filter(lambda x: len(x) > 1)
        
        # Set new boundaries
        for i, v in data.iteritems():
            setattr(self.arcs[i], attr, v)
        self.sort()

    def fuse_close_modules(self, seq_data, max_dist, min_pident, min_size_ratio=0.8):
        """
        - Removed embedded modules from the same parent (keep largest)
        - Fuse consecutive modules if distance <= min_module_size
        - if dist > min_module_size, check pident with N-W
        Arcs need to be (start, end) sorted
        """
        aligner = Align.PairwiseAligner(mode="global", open_gap_score=-0.1)

        for parent in sorted(self.get_all_parents()):
            arcs = self.iter_arcs(wrap=True, parent=parent)
            prev = next(arcs)
            i = 0

            while True:
                try:
                    cur = next(arcs)
                except StopIteration:
                    break
                #===== 1) Check for embedded intervals ====#
                # a) cur C prev => keep prev
                if prev.start <= cur.start and prev.end >= cur.end:
                    cur.flag()
                # b) prev C cur => keep cur
                elif prev.start == cur.start and prev.end < cur.end:
                    prev.flag()
                    prev = cur
                #==== 2) The intervals are not touching
                else:
                    # Check if they are close => we fuse
                    exitcode = prev.try_fuse_with(cur, max_dist)
                    # If not, we align with NW and check the identity
                    if exitcode == 0:
                        # Get the unaligned sequence of ref and query
                        qseq = seq_data[parent]
                        sseq = seq_data[self.ref]
                        if i == 0: # wrap around
                            qseq = qseq[prev.qend:] + qseq[:cur.qstart]
                            sseq = sseq[prev.end:] + sseq[:cur.start]
                        else: # general case
                            # need for modulo if the last is merged
                            qseq = qseq[(prev.qend % len(qseq)):cur.qstart]
                            sseq = sseq[(prev.end % len(sseq)):cur.start]

                        # Align if computationally feasible
                        if len(sseq) < 1e3 and len(qseq) / len(sseq) > min_size_ratio:
                            aln = aligner.align(qseq, sseq)
                            score = aln.score / len(sseq)
                            
                            if score > min_pident:
                                prev.try_fuse_with(cur, force=True)
                    prev = cur
                i += 1
        self.remove_flagged()
        self.sort()

    def merge_equal_intervals(self):
        if len(self) < 2:
            return
        arcs = self.iter_arcs(wrap=True)
        prev_arc = next(arcs)
        for arc in arcs:
            prev_arc.try_merge_with(arc)
            prev_arc = arc

        self.remove_flagged()        

    def fill_gaps(self, max_extend_dist):
        """
        Fill gaps (if any) between arcs
        """
        arcs = self.iter_arcs(wrap=True)
        prev_arc = next(arcs)
        fillers = []

        for arc in arcs:
            if not prev_arc.is_embedded(arc):
                dist = prev_arc.dist_to_next(arc)
                if dist <= 0:
                    pass
                elif dist > max_extend_dist:
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
            if arc.is_embedded(prev_arc, strict=True):
                arc.flag()
            else:
                if prev_arc.is_embedded(arc, strict=True):
                    prev_arc.flag()
                prev_arc = arc
        self.remove_flagged()

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

        if not self.is_covered():
            print(self)
            raise ValueError(f"{self.ref} is not covered. Aborting")

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

        self.remove_flagged()
        self.sort()

    def get_overlap_graph(self, min_overlap=10):
        g = igraph.Graph(directed=True)
        g["ref"] = self.ref
        g["size"] = self.size

        # add vertices
        for i, arc in enumerate(self.iter_arcs(wrap=False)):
            for meta in arc.meta:
                uid = f"{self.ref}.{meta}.{arc.start}-{arc.end}"
                g.add_vertex(name=uid, parent=meta, order=i, ref=self.ref,
                             size=self.size, start=arc.start, end=arc.end)

        if len(self) == 1:
            return g

        # add edges
        arcs = self.iter_arcs(wrap=True)
        prev_arc = next(arcs)

        for i, arc in enumerate(arcs):
            for meta1 in prev_arc.meta:
                uid1 = f"{self.ref}.{meta1}.{prev_arc.start}-{prev_arc.end}"
                start = arc.start
                end = prev_arc.end % self.size # wrapping condition. This should be the only arc with end >= size
                for meta2 in arc.meta:
                    uid2 = f"{self.ref}.{meta2}.{arc.start}-{arc.end}"
                    parents = "/".join(sorted([meta1, meta2]))
                    positions = ((i-1) % len(self), i)
                    g.add_edge(uid1, uid2, ref=self.ref, pos=positions,
                               start=start, end=end, parents=parents)
            prev_arc = arc

        return g
