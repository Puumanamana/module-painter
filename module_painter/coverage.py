from pathlib import Path
import logging
from collections import Counter
from itertools import groupby, chain, product, combinations

from igraph import Graph
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering

from module_painter.arc import Arc
from module_painter.coverage_util import get_filler_arc, get_successor


logger = logging.getLogger("module-painter")

class Coverage:
    def ignore_if_singleton(func):
        def wrapper(*args, **kwargs):
            if len(args[0]) < 2:
                return
            return func(*args, **kwargs)
        return wrapper

    def __init__(self, *arcs, sacc=None):
        self.arcs = list(arcs)
        self.sacc = sacc

        if arcs:
            self.slen = self.arcs[0].slen

        self.sort()

    def __len__(self):
        return sum(1 for arc in self.arcs if not arc.flagged)

    def __repr__(self, subset=None):
        if isinstance(subset, str):
            subset = set(subset)
        display = [f"id: {self.sacc} (L={self.slen})"]
        for i, arc in enumerate(self.arcs):
            if subset is not None:
                if not arc.qacc.intersection(subset):
                    continue
            display.append(f"{i:3} - {arc.sstart:>9}{arc.send:>9}: {arc.qacc}")
        return "\n".join(display)

    def sort(self):
        self.arcs = sorted(self.arcs, key=lambda arc: (arc.sstart, arc.send))

    def remove_flagged(self):
        self.arcs = [arc for arc in self.arcs if not arc.flagged]

    def show(self):
        parents = []
        for i, arc in enumerate(self.arcs):
            parents.append("|".join(sorted(arc.qacc)))# + f"({arc.send-arc.sstart:,})")
        parents = "->".join(f"{p:^5}" for p in parents)
        return f"(ID: {self.sacc})  >>>  {parents}"

    def iter(self, wrap=False, parent=None):
        """
        Important: for wrap=True, self.arcs needs to be (sstart, send) sorted
        """
        arcs = [arc for arc in self.arcs if not arc.flagged]
    
        if parent is not None:
            arcs = [arc for arc in arcs if parent in arc.qacc]

        if wrap:
            max_end = max(arc.send for arc in arcs)
            last_arc = next(arc for arc in arcs if arc.send == max_end)
            arcs = [last_arc] + arcs

        return iter(arcs)

    def get_all_parents(self):
        parents = set()
        for arc in self.iter(wrap=False):
            parents |= arc.qacc
        return parents
    
    @classmethod
    def from_coverages(cls, *covs):
        slen = covs[0].slen
        sacc = covs[0].sacc
        assert all(cov.slen==slen for cov in covs)
        assert all(cov.sacc==sacc for cov in covs)

        arcs = [arc for cov in covs for arc in cov.arcs]

        return cls(*arcs, sacc=sacc)

    @classmethod
    def from_pandas(cls, df):
        sacc = df.sacc.iloc[0]
        slen = df["slen"].iloc[0]
        assert all(df.sacc == sacc)

        arcs = [Arc(row.sstart, row.send, slen, qacc={row.qacc}, qstart=row.qstart, qend=row.qend)
                for _, row in df.iterrows()]

        return cls(*arcs, sacc=sacc)

    @classmethod
    def from_csv(cls, fpath):
        df = pd.read_csv(fpath)
        return cls.from_pandas(df)
        
    def to_pandas(self):
        n = len(self)
        data = {attr: np.zeros(n, dtype=int) for attr in ["sstart", "send", "qstart", "qend", "flag"]}
        data["qacc"] = np.empty(n, dtype="<U32")

        for i, arc in enumerate(self.arcs):
            data["sstart"][i] = arc.sstart
            data["send"][i] = arc.send
            data["qacc"][i] = list(arc.qacc)[0]
            data["flag"][i] = arc.flagged

        df = pd.DataFrame(data).assign(sacc=self.sacc, slen=self.slen)

        return df
    
    def to_csv(self, fpath):
        Path(fpath).parent.mkdir(exist_ok=True)
        df = self.to_pandas()
        df.to_csv(fpath)

    def is_covered(self):
        """
        Assumes self.arcs is (sstart, send) sorted
        """
        arcs = self.iter(wrap=True)
        current_arc = next(arcs)

        for arc in arcs:
            dist = current_arc.dist_to_next(arc)
            
            if dist > 0: # There is a hole
                print(f"{current_arc.send} < {arc.sstart} (dist={dist})")
                return False

            if not arc.is_embedded(current_arc):
                current_arc = arc

        return True

    @ignore_if_singleton
    def get_breakpoints(self, skip_na=True):
        arcs = list(self.iter(wrap=True))
        breakpoints = zip(arcs[:-1], arcs[1:])
        
        if skip_na:
            breakpoints = filter(lambda bk: all("@" not in arc.qacc for arc in bk), breakpoints)
        
        return breakpoints

    @ignore_if_singleton
    def iter_recombinations(self, skip_na=True):
        bk = self.get_breakpoints(skip_na=skip_na)

        for bk1, bk2 in combinations(bk, 2):
            inter = [pi.qacc.intersection(qi.qacc) for pi in bk1 for qi in bk2]

            if inter[1] and inter[2]:
                yield inter[1], inter[2]
            elif inter[0] and inter[-1]:
                yield inter[0], inter[-1]

    def count_rc(self, skip_na=True):
        if len(self) < 2:
            return 0
        return sum(1 for _ in self.iter_recombinations(skip_na=skip_na))

    @ignore_if_singleton
    def sync_boundaries(self, attr, max_dist):
        """
        Change intervals boundaries to make them similar
        """
        features = pd.Series([getattr(arc, attr) for arc in self.iter(wrap=False)])

        model = AgglomerativeClustering(
            linkage="complete", affinity="l1",
            distance_threshold=max_dist, n_clusters=None
        ).fit(features.values[:, None])

        data = features.groupby(model.labels_).transform(min if "sstart" in attr else max)
        data = data.groupby(model.labels_).filter(lambda x: len(x) > 1)

        logger.debug(f"Changing {attr} boundary for {data.shape[0]} arcs")
        
        # Set new boundaries
        for i, v in data.iteritems():
            setattr(self.arcs[i], attr, v)
            self.arcs[i].fix_boundaries()
        self.sort()

    @ignore_if_singleton
    def fuse_close_modules(self, seq_data, **kw):
        """
        - Removed embedded modules from the same parent (keep largest)
        - Fuse consecutive modules if distance <= min_module_size
        - if dist > min_module_size, check pident with N-W
        Arcs need to be (sstart, send) sorted
        """
        for parent in sorted(self.get_all_parents()):
            # First, we remove any embedded modules
            self.simplify_embedded(parent=parent)
            # From this point onwards, intervals are sorted on both sstart and send
            iterator = self.iter(wrap=False, parent=parent)
            prev = next(iterator)

            (qseq, sseq) = (seq_data.get(parent), seq_data.get(self.sacc))
            while True:
                try:
                    cur = next(iterator)
                except StopIteration:
                    break
                prev.try_fuse_with(cur, qseq=qseq, sseq=sseq, **kw)
                prev = cur
            
            self.remove_flagged()

            # check if we can circularize
            iterator = self.iter(wrap=True, parent=parent)
            last = next(iterator)
            first = next(iterator)
            last.try_fuse_with(first, qseq=qseq, sseq=sseq, **kw)
            
        for arc in self.arcs:
            arc.fix_boundaries()

        self.remove_flagged()
        self.sort()

    @ignore_if_singleton
    def merge_equal_intervals(self):
        arcs = self.iter(wrap=False)
        prev_arc = next(arcs)
        for arc in arcs:
            prev_arc.try_merge_with(arc)
            prev_arc = arc

        self.remove_flagged()        

    @ignore_if_singleton
    def fill_gaps(self, max_extend_dist):
        """
        Fill gaps (if any) between arcs
        Assumes no embedding between 2 consecutive arcs
        """
        arcs = self.iter(wrap=True)
        prev_arc = next(arcs)
        fillers = []

        for arc in arcs:
            is_extended = prev_arc.try_extend_end_with(arc, max_extend_dist)
            if not is_extended: # add a filler arc
                filler = get_filler_arc(prev_arc, arc)
                fillers.append(filler)
            prev_arc = arc

        self.arcs += fillers
        self.sort()
    
    @ignore_if_singleton
    def simplify_embedded(self, parent=None):
        arcs = self.iter(wrap=True, parent=parent)
        prev_arc = next(arcs)

        for arc in arcs:
            if arc.is_embedded(prev_arc, strict=True):
                # logger.debug(f"Discarding {arc.qacc} (embedded in {prev_arc.qacc})")
                arc.flag()
            else:
                if prev_arc.is_embedded(arc, strict=True):
                    # logger.debug(f"Discarding {prev_arc.qacc} (embedded in {arc.qacc})")
                    prev_arc.flag()
                prev_arc = arc
        self.remove_flagged()

    @ignore_if_singleton        
    def get_minimal_coverage(self):
        """
        Lee and lee 1983 algorithm
        On a circle-cover minimization problem
        Conditions:
        - self is fully covered (no gaps)
        - No embedded intervals
        """

        arcs = [arc for arc in self.arcs if not arc.flagged]

        if not self.is_covered():
            print(self)
            logger.error(f"{self.sacc} is not covered. Aborting")
            exit(1)

        successors = [get_successor(i, arcs) for (i, arc) in enumerate(arcs)]

        # Order of the B_i
        t_ = np.zeros(len(self), dtype=int)

        (i, B_1) = (1, 0)

        S = [arcs[B_1]]
        start0, current_end = arcs[B_1].bounds()

        for arc in arcs:
            arc.marked = False

        arcs[B_1].marked = True # Mark the 1st arc
        t_[B_1] = 1
        B_i = B_1

        while self.slen-current_end+(start0-1) > 0:
            B_i = successors[B_i]
            S.append(arcs[B_i])
            arcs[B_i].marked = True

            current_end = arcs[B_i].send

            t_[B_i] = i + 1
            i += 1
            if i > 1e6:
                print("infinite loop...")
                print(self)
                import ipdb;ipdb.set_trace()

        # At this point we have an initial cover
        k = i

        is_opt = False

        while not is_opt:
            i += 1
            if i > 1e6:
                print("infinite loop...")
                print(self)
                import ipdb;ipdb.set_trace()
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
        g = Graph(directed=True)
        g["sacc"] = self.sacc
        g["slen"] = self.slen

        # add vertices
        for i, arc in enumerate(self.iter(wrap=False)):
            for qacc in arc.qacc:
                uid = f"{self.sacc}.{qacc}.{arc.sstart}-{arc.send}"
                g.add_vertex(name=uid, parent=qacc, order=i, sacc=self.sacc,
                             slen=self.slen, sstart=arc.sstart, send=arc.send)

        if len(self) == 1:
            return g

        # add edges
        arcs = self.iter(wrap=True)
        prev_arc = next(arcs)

        for i, arc in enumerate(arcs):
            for qacc1 in prev_arc.qacc:
                uid1 = f"{self.sacc}.{qacc1}.{prev_arc.sstart}-{prev_arc.send}"
                sstart = arc.sstart
                send = prev_arc.send
                # wrapping condition. This should be the only arc with send >= slen-1
                if prev_arc.send == self.slen-1:
                    send = 0
                elif prev_arc.send > self.slen:
                    send -= self.slen
                for qacc2 in arc.qacc:
                    uid2 = f"{self.sacc}.{qacc2}.{arc.sstart}-{arc.send}"
                    parents = "/".join(sorted([qacc1, qacc2]))
                    positions = ((i-1) % len(self), i)
                    g.add_edge(uid1, uid2, sacc=self.sacc, pos=positions,
                               sstart=sstart, send=send, parents=parents)
            prev_arc = arc

        return g

    def apply_rules(self, arc_eq_diffs, min_module_size, seq_data, min_nw_id, skip_nw=False):
        self.sync_boundaries("sstart", arc_eq_diffs)
        self.sync_boundaries("send", arc_eq_diffs)
        self.fuse_close_modules(seq_data, max_dist=min_module_size,
                                min_size_ratio=0.8, min_nw_id=min_nw_id, skip_nw=skip_nw)
        self.merge_equal_intervals()
        self.simplify_embedded()
        self.fill_gaps(min_module_size)
        self.get_minimal_coverage()

        logger.info(f"{self.sacc}: {self.count_rc()} recombinations found")
