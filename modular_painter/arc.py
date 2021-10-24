class Arc:
    def __init__(self, start, end, size, meta="X"):
        self.size = size
        self.start = start
        self.end = end
        self.meta = meta
        self.flagged = False

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
            self.size == other.size
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
        arc_1 = Arc(self.start, self.size-1, self.size)
        arc_2 = Arc(0, self.end % self.size, self.size)

        return (arc_1, arc_2)

    def is_embedded(self, other, strict=True):
        """
        is self is embedded in other?
        """
        # Special case: equal intervals
        if strict and self.bounds() == other.bounds():
            return False
        # General case: no looping around or both looping around
        if other.start <= self.start <= self.end <= other.end:
            return True
        # "other" loops around
        if other.end >= other.size and 0 <= self.start <= self.end <= other.end-self.size:
            return True
        # if only "self" loops around, it's not embedded
        return False

    def dist_to_next(self, other):
        """
        assymmetric distance between two intervals.
        """
        # case 1: both intervals wrap around
        if self.end >= self.size and other.end >= self.size:
            return 0
        # case 2: wrap around
        if other.start < self.start:
            return max(0, other.start + (self.size-self.end))
        # general case
        return max(0, other.start - self.end)

    def try_merge_with(self, other):
        """
        Merge self with other if they have the same start and end
        The new meta is the union of both meta
        Flag "other"
        """
        if self.bounds() == other.bounds():
            self.flag()
            other.meta |= self.meta

    def try_fuse_with(self, other, max_dist, dist=None):
        """
        Extend self with other if they share parents
        """
        shared_meta = self.meta.intersection(other.meta)

        if not shared_meta:
            return
        if dist is None:
            dist = self.dist_to_next(other)
        if dist > max_dist:
            return
        self.flag()
        other.unflag()
        other.meta = shared_meta

        if self.start < other.start: # General case: 1--1...2--2
            other.start = self.start
        else: # Loop: ..2---2...1--1..
            other.start = self.start
            other.end += self.size

    def try_extend_end_with(self, other, max_dist, dist=None):
        """
        1. Check distance between self and other, proceed if < max_dist
        2. Extend evenly both intervals to touch. Intervals should not 
           share parents (should have already been fused)
        """
        if dist is None:
            dist = self.dist_to_next(other)

        if dist > max_dist:
            return

        if self.start <= other.start: # General case: 1--1...2--2
            left_extend = dist // 2            
        else: # Loop: ..2---2...1--1..
            # We prioritize extending towards 0
            left_extend = min(other.start, dist)
        right_extend = dist - left_extend
        self.end += right_extend
        other.start -= left_extend
