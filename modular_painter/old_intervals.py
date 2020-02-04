import sys
from bisect import bisect

import pandas as pd

class Arc:

    def __init__(self, names, start, end, max_len, identities):
        self.start = start
        self.end = end
        self.mod = max_len

        if not isinstance(names, list):
            names = [names]
        if not isinstance(identities, list):
            identities = [identities]
        self.identities = pd.Series(identities, index=names)

        if end < start:
            self.end += max_len

    def __repr__(self):
        return "{}: ({}, {}), max: {}, id: {}".format(
            '/'.join(self.identities.tolist()),
            self.start,
            self.end,
            self.mod)

    def __eq__(self, other):
        if len(self.identities) != len(other.identities):
            return False

        return (
            self.identities.to_dict() == other.identities.to_dict()
            and (self.start == other.start)
            and (self.end == other.end)
            and (self.mod == other.mod)
        )

    def len(self):
        return self.end - self.start

    def intersect(self, other):
        if self.mod != other.mod:
            sys.exit("Cannot compare, different max length {} vs {}".format(self.mod, other.mod))

        if ((other.start <= self.start <= other.end)
            or (other.start <= self.end % self.mod <= other.end)
            or (self.start <= other.end % other.mod <= self.end)):
            return True

        return False

class Coverage:

    def __init__(self, target, *arcs):
        self.target = target
        self.arcs = sorted(arcs, key=lambda x: x.start)
        self.starts = [arc.start for arc in self.arcs]
        self.mod = arcs[0].mod

    @classmethod
    def from_pandas(cls, data, target, max_len=-1):
        if max_len < 0:
            max_len = data.loc[:, ['tstart', 'tend']].max().max()

        cov = cls(target, *[
            Arc(row.source, row.tstart, row.tend, max_len, row.identities)
            for _, row in data.iterrows()
        ])

        return cov

    def __getitem__(self, i):
        return self.arcs[i]

    def __len__(self):
        return len(self.arcs)

    def __repr__(self):
        return 'Target: {}\n'.format(self.target) +'\n'.join(arc.__str__() for arc in self.arcs)
    
    def to_pandas(self):
        data = [[arc.identities.index.tolist(),
                 arc.start,
                 arc.end,
                 arc.identities.tolist()] for arc in self.arcs]
        df = pd.DataFrame(data, columns=['source', 'tstart', 'tend', 'identities'])

        return df

    def simplify(self):
        '''
        - Remove any interval that is strictly embedded in another
        - Assumes no equal interval boundaries
        '''

        self.arcs = sorted(self.arcs, key=lambda x: (x.start, x.end))

        cum_max = self.arcs.end[0]

        for i, arc in enumerate(self.arcs):
            if i < len(self):
                # Keep only the last repeated entry
                if arc.start == self[i+1].start:
                    self.arcs.remove(arc)
            
            
            
        self.starts = [arc.start for arc in arcs]
        return
        
    def insert(self, other):
        '''
        Does not check if the arc already exists
        '''

        if self.mod != other.mod:
            sys.exit('Cannot insert: the 2 intervals have different max_len ({}, {})'
                     .format(self.mod, other.mod))

        ins_pos = bisect(self.starts, other.start)
        self.arcs.insert(ins_pos, other)
        self.starts.insert(ins_pos, other.start)

    def get_successor(self, i):
        offset = 0
        end_i = self[i].end
        
        not_found = True
        j = i

        while not_found:
            j += 1
        
            if j >= len(self):
                j = 0
                offset = self.mod
            
            if self[j].start + offset > end_i:
                return (j - 1) % len(self)

    def set_all_successors(self):
        self.successors = [self.get_sucessor(i) for i in range(len(self))]

    def fill_holes(self, max_merge_dist):
        i_max = 0
        current_end = self.arcs[0].end

        fillers = []
        for i, arc in enumerate(self.arcs[1:], 1):
            dist = arc.start - current_end

            if 0 < dist <= max_merge_dist:
                # extend the previous arc
                self[i-1].end += dist
            elif dist > max_merge_dist:
                filler = Arc(['NoCov'], current_end, arc.start, self.mod, 0)
                fillers.append(filler)

            if arc.end >= current_end:
                # update the max
                i_max = i
                current_end = arc.end

        # Looping condition
        loop_dist = (self[0].start + self.mod) - current_end

        if loop_dist > max_merge_dist:
            filler = Arc(['NoCov'], current_end, arc.start, self.mod, 0)
            fillers.append(filler)

        if 0 < loop_dist <= max_merge_dist:
            self[i_max].end += loop_dist
            
        # Add the fillers
        for filler in fillers:
            self.insert(filler)

    def is_covered(self):
        current_end = self.arcs[0].end

        for arc in self.arcs:
            if arc.start > current_end: # There is a hole
                return False
            if arc.end > current_end:
                current_end = arc.end

        # Looping condition
        if current_end % self.mod >= self.arcs[0].start:
            return True
        return False
