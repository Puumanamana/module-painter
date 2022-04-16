class Arc:
    def __init__(self, start, end, size, meta="X", **attrs):
        self.size = size
        self.meta = meta
        self.flagged = False
        self.set_bounds(start, end)

        for key, val in attrs.items():
            setattr(self, key, val)
        

    def set_bounds(self, start, end):
        if end <= start:
            raise ValueError(f"{self.meta}: {start} >= {end}")
        if end-start >= self.size:
            self.start = 0
            self.end = self.size - 1
        else:
            self.start = start
            self.end = end

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
        - If other starts before self, then we wrap around
        - If the intervals intersect, returns the overlap 
        (negative distance)
        """
        # case 1: self wraps around, not other
        if other.end < self.size <= self.end:
            return other.start - (self.end % self.size) - 1
        # case 2: other precedes self, but self does not wrap
        if other.start < self.start < self.end:
            return self.size - self.end + other.start - 1
        # general case
        return other.start - self.end - 1

    def try_merge_with(self, other):
        """
        Merge self with other if they have the same start and end
        The new meta is the union of both meta
        Flag "other"
        """
        if self.bounds() == other.bounds():
            self.flag()
            other.meta |= self.meta

    def try_fuse_with(self, other, max_dist=0, dist=None, force=False):
        """
        Extend self with other if they share parents
        Return 1 if fused else 0
        """
        shared_meta = self.meta.intersection(other.meta)

        if not shared_meta:
            return 0
        if dist is None:
            dist = self.dist_to_next(other)
        if not force and dist > max_dist:
            return 0
        self.flag()
        other.unflag()
        other.meta = shared_meta

        if self.start < other.start: # General case: 1--1...2--2
            other.start = self.start
        else: # Loop: ..2---2...1--1..
            other.start = self.start
            other.end += self.size
        return 1

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
