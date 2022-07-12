import logging
from Bio import Align


logger = logging.getLogger("module-painter")

class Arc:
    def __init__(self, sstart, send, slen, qacc="X", **attrs):
        self.slen = slen
        self.qacc = qacc
        self.flagged = False
        self.sstart = sstart
        self.send = send
        
        for key, val in attrs.items():
            setattr(self, key, val)

        self.fix_boundaries()        

    def fix_boundaries(self):
        seq_len = self.send-self.sstart+1
        if seq_len <= 0:
            logger.error("send < sstart for arc: {self.qacc}")
            logger.error(self.__repr__())
            exit(1)
        if seq_len >= self.slen:
            self.sstart = 0
            self.send = self.slen - 1

    def __repr__(self):
        return "{}: ({}, {}), slen={} (flagged: {})".format(
            self.qacc,
            self.sstart,
            self.send,
            self.slen,
            self.flagged
        )

    def __eq__(self, other):
        return (
            self.sstart == other.sstart and
            self.send == other.send and
            self.slen == other.slen and
            self.qacc == other.qacc
        )

    def __len__(self):
        return self.send - self.sstart + 1

    def flag(self):
        self.flagged = True

    def unflag(self):
        self.flagged = False

    def bounds(self):
        return (self.sstart, self.send)

    def split_at_end(self):
        if self.send < self.slen:
            logger.error(f"Arc {self.qacc} does not wrap around, cannot split at send")
            exit(1)
        arc_1 = Arc(self.sstart, self.slen-1, self.slen, qacc=self.qacc)
        arc_2 = Arc(0, self.send % self.slen, self.slen, qacc=self.qacc)

        return (arc_1, arc_2)

    def is_embedded(self, other, strict=False):
        """
        is self is embedded in other?
        """
        # Special case: equal intervals
        if strict and self.bounds() == other.bounds():
            return False
        # self loops around
        if self.send >= self.slen:
            # split interval in 2
            (arc1, arc2) = self.split_at_end()
            return arc1.is_embedded(other) and arc2.is_embedded(other)
        # other loops around
        if other.send >= other.slen:
            (arc1, arc2) = other.split_at_end()
            return self.is_embedded(arc1) or self.is_embedded(arc2)
        # General case: no looping around
        return other.sstart <= self.sstart <= self.send <= other.send

    def dist_to_next(self, other):
        """
        Distance from self to other (note: assymmetric)
        - Assumes no embedding relationship (undefined behavior)
        - If other starts before self, then we wrap around
        - If the intervals just intersect, returns the overlap 
        (negative distance)
        """
        # case 1: other precedes self => we wrap around
        if other.sstart <= self.sstart:
            return (self.slen-self.send) + other.sstart - 1
        # case 2: both interval extend past slen
        if other.send >= other.slen and self.send >= self.slen:
            return other.sstart - self.send - 1
        # case 3: General case
        return other.sstart - (self.send % self.slen) - 1


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

        # logger.debug("Extending: {self.qacc}:{self.sstart:,}-{self.send:,} -> {other.qacc}:{other.sstart:,}-{other.send:,}")

        if self.sstart <= other.sstart: # General case: 1--1...2--2
            left_extend = dist // 2
        else: # Loop: ..2---2...1--1..
            # We prioritize extending towards 0
            left_extend = min(other.sstart, dist)
        right_extend = dist - left_extend
        self.send += right_extend
        other.sstart -= left_extend

        return True

    def try_merge_with(self, other):
        """
        Merge self with other if they have the same sstart and send
        The new qacc is the union of both qacc
        Flag "other"
        """
        if self.bounds() == other.bounds():
            self.flag()
            other.qacc |= self.qacc

    def check_similarity(self, other, qseq=None, sseq=None, min_size_ratio=None):
        aligner = Align.PairwiseAligner(mode="global", open_gap_score=-0.1)
        
        if self.sstart >= other.sstart: # go to end and come back
            qseq = qseq[self.qend:] + qseq[:other.qstart+1]
            sseq = sseq[self.send:] + sseq[:other.sstart+1]
        else: # general case
            # need for modulo if self wraps around
            qseq = qseq[(self.qend % len(qseq)):other.qstart+1]
            sseq = sseq[(self.send % len(sseq)):other.sstart+1]

        # Align if computationally feasible
        if len(sseq) < 1e3 and len(qseq) / len(sseq) > min_size_ratio:
            aln = aligner.align(qseq, sseq)
            score = aln.score / len(sseq)
            return score
        return -1
        
    def try_fuse_with(self, other, max_dist=None, min_nw_id=None, skip_nw=False, **aln_kw):
        """
        Fuse self with other into one arc with the qacc shared by both arcs
        Flag self, unflag other
        """
        dist = self.dist_to_next(other)

        if dist > max_dist:
            if skip_nw:
                return
            nw_score = self.check_similarity(other, **aln_kw)
            if nw_score < min_nw_id:
                return

        shared_qacc = self.qacc.intersection(other.qacc)
        if not shared_qacc:
            return

        # we fuse
        self.flag()
        other.unflag()
        other.qacc = shared_qacc

        # logger.debug(f"Fusing: {self.qacc}:{self.sstart:,}-{self.send:,} -> {other.qacc}:{other.sstart:,}-{other.send:,}")
        
        if self.sstart < other.sstart: # General case: 1--1...2--2
            other.sstart = self.sstart
        else: # Loop: ..2---2...1--1..
            other.sstart = self.sstart
            other.send += self.slen
        
