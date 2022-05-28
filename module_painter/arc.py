from Bio import Align

class Arc:
    def __init__(self, start, end, size, meta="X", **attrs):
        self.size = size
        self.meta = meta
        self.flagged = False
        self.start = start
        self.end = end
        
        for key, val in attrs.items():
            setattr(self, key, val)

        self.fix_boundaries()        

    def fix_boundaries(self):
        seq_len = self.end-self.start+1
        if seq_len <= 0:
            raise ValueError(self.__repr__())
        if seq_len >= self.size:
            self.start = 0
            self.end = self.size - 1

    def __repr__(self):
        return "{}: ({}, {}), size={} (flagged: {})".format(
            self.meta,
            self.start,
            self.end,
            self.size,
            self.flagged
        )

    def __eq__(self, other):
        return (
            self.start == other.start and
            self.end == other.end and
            self.size == other.size and
            self.meta == other.meta
        )

    def __len__(self):
        return self.end - self.start + 1

    def flag(self):
        self.flagged = True

    def unflag(self):
        self.flagged = False

    def bounds(self):
        return (self.start, self.end)

    def split_at_end(self):
        if self.end < self.size:
            assert ValueError("Arc does not wrap around, cannot split at end")
        arc_1 = Arc(self.start, self.size-1, self.size, meta=self.meta)
        arc_2 = Arc(0, self.end % self.size, self.size, meta=self.meta)

        return (arc_1, arc_2)

    def is_embedded(self, other, strict=False):
        """
        is self is embedded in other?
        """
        # Special case: equal intervals
        if strict and self.bounds() == other.bounds():
            return False
        # self loops around
        if self.end >= self.size:
            # split interval in 2
            (arc1, arc2) = self.split_at_end()
            return arc1.is_embedded(other) and arc2.is_embedded(other)
        # other loops around
        if other.end >= other.size:
            (arc1, arc2) = other.split_at_end()
            return self.is_embedded(arc1) or self.is_embedded(arc2)
        # General case: no looping around
        return other.start <= self.start <= self.end <= other.end

    def dist_to_next(self, other):
        """
        Distance from self to other (note: assymmetric)
        - Assumes no embedding relationship (undefined behavior)
        - If other starts before self, then we wrap around
        - If the intervals just intersect, returns the overlap 
        (negative distance)
        """
        # case 1: other precedes self => we wrap around
        if other.start <= self.start:
            return (self.size-self.end) + other.start - 1
        # case 2: both interval extend past size
        if other.end >= other.size and self.end >= self.size:
            return other.start - self.end - 1
        # case 3: General case
        return other.start - (self.end % self.size) - 1


    def try_extend_end_with(self, other, max_dist, dist=None):
        """
        1. Check distance between self and other, proceed if 0 < d < max_dist
        2. Extend evenly both intervals to touch. Intervals should not 
           share parents (should have already been fused)
        """
        if dist is None:
            dist = self.dist_to_next(other)
        if dist <= 0:
            return True # already touching
        if dist > max_dist:
            return False # too far

        if self.start <= other.start: # General case: 1--1...2--2
            left_extend = dist // 2
        else: # Loop: ..2---2...1--1..
            # We prioritize extending towards 0
            left_extend = min(other.start, dist)
        right_extend = dist - left_extend
        self.end += right_extend
        other.start -= left_extend

        return True

    def try_merge_with(self, other):
        """
        Merge self with other if they have the same start and end
        The new meta is the union of both meta
        Flag "other"
        """
        if self.bounds() == other.bounds():
            self.flag()
            other.meta |= self.meta

    def check_similarity(self, other, qseq=None, sseq=None, min_size_ratio=None, max_len=1e3):
        aligner = Align.PairwiseAligner(mode="global", open_gap_score=-0.1)
        
        if self.start >= other.start: # go to end and come back
            qseq = qseq[self.qend:] + qseq[:other.qstart+1]
            sseq = sseq[self.end:] + sseq[:other.start+1]
        else: # general case
            # need for modulo if self wraps around
            qseq = qseq[(self.qend % len(qseq)):other.qstart+1]
            sseq = sseq[(self.end % len(sseq)):other.start+1]

        # Align if computationally feasible
        if len(sseq) < 1e3 and len(qseq) / len(sseq) > min_size_ratio:
            aln = aligner.align(qseq, sseq)
            score = aln.score / len(sseq)
            return score
        return -1
        
    def try_fuse_with(self, other, max_dist=None, nw=True, min_nw_id=None, **aln_kw):
        dist = self.dist_to_next(other)

        if dist > max_dist and not nw:
            return
        if dist > max_dist:
            nw_score = self.check_similarity(other, **aln_kw)
            if nw_score < min_nw_id:
                return

        shared_meta = self.meta.intersection(other.meta)
        if not shared_meta:
            return

        # we fuse
        self.flag()
        other.unflag()
        other.meta = shared_meta

        if self.start < other.start: # General case: 1--1...2--2
            other.start = self.start
        else: # Loop: ..2---2...1--1..
            other.start = self.start
            other.end += self.size
        
