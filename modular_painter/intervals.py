from itertools import groupby, compress
import sys
import pandas as pd
import numpy as np

def get_sucessor(data, i, mod):
    '''
    - Data must not include any duplicated intervals
    - The while interval must be covered
    '''

    offset = 0
    end_i = data[i, 1]

    not_found = True
    j = i

    while not_found:
        j += 1
        
        if j >= len(data):
            j = 0
            offset = mod
            
        if data[j, 0] + offset > end_i + 1:
            return (j - 1) % len(data)


def compare_arc_pair(arc_prev, arc_next, max_extend_dist):
    '''
    Check if arcs are overlapping or close enough and merge them if they have common sources
    - arc_prev and arc_next need to be sorted in increasing start order
    - The 2 arcs cannot have the exact same bounds
    '''

    close = (arc_next.start - arc_prev.end) % arc_next.mod <= max_extend_dist
    overlap = arc_prev.intersect(arc_next)

    if arc_next.end <= arc_prev.end and overlap:
        return ('discard', arc_prev)

    if arc_prev.start == arc_next.start and arc_prev.end <= arc_next.end:
        return ('discard_prev', arc_next)

    if close or overlap:
        common = np.intersect1d(arc_prev.data.index, arc_next.data.index)
        if common.size > 0:
            merger_data = arc_prev.data.loc[common].reset_index()
            # Set the end to the next end.
            merger_data['end'] = arc_next.end
            return ('merge_extend', Arc(merger_data, arc_prev.target, arc_prev.mod))

        elif close and not overlap:
            return ('extend', arc_next)
        
        return ('do_nothing', arc_next)
            
    # The two arcs are too far away. We make up a fake arc
    filler = pd.DataFrame([['NoCov', arc_prev.end+1, arc_next.start-1, 0]],
                          columns=['source', 'start', 'end', 'identity'])
    return ('fill', Arc(filler, arc_prev.target, arc_prev.mod))

class Arc:

    def __init__(self, data, target, modulo):
        self.mod = modulo
        self.target = target
        self.start = data.start[0]
        self.end = data.end[0]
        self.fix_bounds()

        self.data = data.set_index('source')
        self.marked = False # For lee and lee algorithm        
        
    def __repr__(self):
        return "{}: ({}, {}), modulo={}".format(
            '/'.join(self.data.index),
            *self.bounds(), self.mod)

    def __eq__(self, other):
        return self.data.to_dict() == other.data.to_dict()

    def __len__(self):
        return self.end - self.start + 1

    def set(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def bounds(self):
        return [self.start, self.end]

    def fix_bounds(self):
        if self.start >= self.mod:
            import ipdb;ipdb.set_trace()
            sys.exit("{} >= {}: this should never happen".format(self.start, self.mod))
        if self.end < self.start:
            self.end += self.mod

    def intersect(self, other):
        '''
        - 2 arcs intersect if they have no gap inbetween
        Thus, [0, 1] and [2,4] intersect since we work with discrete intervals
        '''
        if self.target != other.target or self.mod != other.mod:
            sys.exit("Cannot compare, different target {} vs {}".format(self.target, other.target))

        # 3 cases: 1) easy, 2) self spans 0, 3) other spans 0
        # |---|    
        #        0
        #           |-----|
        start1, end1 = self.bounds()
        start2, end2 = other.bounds()
        if self.end >= self.mod:
            start1, end1 = (start1-self.mod, end1-self.mod)
        if other.end >= other.mod:
            start2, end2 = (start2-other.mod, end2-other.mod)

        return (
            start1 - 1 <= start2 <= end1 + 1
            or start1 - 1 <= end2 <= end1 + 1
            or (end1 + 1) % self.mod == start2
            or (end2 + 1) % self.mod == start1
        )

    def mark(self):
        self.marked = True
    
class Coverage:

    def __init__(self, *arcs):
        if arcs:
            self.target = arcs[0].target
            self.mod = arcs[0].mod
            self.data = sorted(arcs, key=lambda x: (x.bounds()))

        self.dirty = False

    @classmethod
    def from_pandas(cls, data, target, modulo=-1):
        '''
        start/end as columns in data
        '''
        if modulo < 0:
            modulo = data[['start', 'end']].max().max()
            
        # data.loc[data.tend < data.tstart, 'tend'] += modulo
        df = data.set_index(['start', 'end']).sort_index()

        arcs = [Arc(df.loc[idx].reset_index(), target, modulo) for idx in df.index.unique()]

        return cls(*arcs)

    def to_pandas(self):
        result = []

        for arc in self.data:
            df = arc.data.reset_index()
            df['start'] = arc.start
            df['end'] = arc.end
            result.append(df)

        return pd.concat(result)

    def __len__(self):
        return len(self.data)

    def __repr__(self):

        locations = ['{:3} - {:>8.1%}{:>8.1%}: {}'.format(i, arc.start/self.mod, arc.end/self.mod, '/'.join(arc.data.index))
                     for i, arc in enumerate(self.data)]
        message = "Target: {} (L={})\n{}".format(self.target, self.mod, '\n'.join(locations))
        return message

    def __getitem__(self, i):
        if self.dirty:
            self.data = sorted(self.data, key=lambda x: x.bounds())
            self.dirty = False
        return self.data[i]

    def get_intervals(self):
        return np.array([self[i].bounds() for i in range(len(self))])

    def is_covered(self):
        current_end = self[0].end

        for i in range(len(self)):
            (start, end) = self[i].bounds()
            if start > current_end + 1: # There is a hole
                print("{} > {}".format(start, current_end))
                return False
            if end > current_end:
                current_end = end

        # Looping condition
        return self[-1].intersect(self[0])

    def merge_close_intervals(self, arc_eq_diffs):
        '''
        Group intervals that have about the same boundaries 
        '''

        intervals = self.get_intervals()
        dists = np.abs(intervals - np.roll(intervals, 1, axis=0)).sum(axis=1)
        # is_far has an incrementing index that stay constant when dist <= arc_eq_diffs
        is_far = np.cumsum(dists > arc_eq_diffs)

        final_indices = np.ones(len(self)).astype(bool)

        for _, indices in groupby(range(len(self)), key=lambda i: is_far[i]):
            indices = list(indices)

            if len(indices) > 1:
                start = intervals[indices, 0].min()
                end = intervals[indices, 1].max()

                i0 = indices.pop()
                self[i0].set(start=start, end=end, data=pd.concat([self[j].data for j in indices]))

                final_indices[indices] = False

        self.data = list(compress(self.data, final_indices))
        
    def fill_or_extend(self, max_extend_dist):
        '''
        ISSUE: because we delete arcs, we might still have consecutive arcs with the same label
        FIX: 
           1 - Loop to merge consecutive element
           2 - Sort by (start, end) position --> groupby "start" key--> last element
        '''
        
        i_max = 0
        i = 1
        to_discard = set()

        while True:
            if i < len(self):
                arc = self[i]
            else:
                break

            # Decide whether we should merge or fill gap between the 2 arcs or discard an arc
            outcome, new_arc = compare_arc_pair(self[i_max], arc, max_extend_dist)

            if outcome == 'discard_prev':
                to_discard.add(i_max % len(self))

            elif outcome in {'discard', 'merge_extend'}:
                # 2 cases that require to delete the current arc:
                # - merge_extend: 2 consecutive arcs with the same source --> group
                # - discard: the arc is embedded in the previous
                if i % len(self) != i_max:
                    to_discard.add(i % len(self))

                if outcome == 'merge_extend':
                    self.data[i_max] = new_arc

            elif outcome == 'extend':
                # small gap between prev and current: we extend prev to current.start
                self.data[i_max].end = arc.start - 1
                self.data[i_max].fix_bounds()

            elif outcome == 'fill':
                # Fake arc to fill the gap between the 2 neighbors
                self.data.insert(i, new_arc)
                i += 1 # we increment one more because we inserted

            if new_arc.end > self[i_max].end:
                i_max = i
                
            i += 1

        if not self[0].intersect(self[i_max]):
            start0, end0 = self[0].bounds()
            startf, endf = self[i_max].bounds()
    
            if (start0 + self.mod) - endf <= max_extend_dist:
                self[i_max].end = self.mod + (start0 - 1)
            else:
                filler = pd.DataFrame([['NoCov', endf+1, start0-1 % self.mod, 0]],
                                      columns=['source', 'start', 'end', 'identity'])
                self.data.append(Arc(filler, self.target, self.mod))

        self.data = [arc for (i, arc) in enumerate(self.data) if i not in to_discard]

        if not self.is_covered():
            print('Something went wrong. Coverage is missing some parts')
            import ipdb;ipdb.set_trace()

    def set_all_successors(self):
        arc_boundaries = self.get_intervals()
        self.successors = [get_sucessor(arc_boundaries, i, self.mod) for i in range(len(self))]

    def get_minimal_coverage(self):
        '''
        Lee and lee 1983 algorithm
        On a circle-cover minimization problem
        '''

        n_interv = len(self)

        if n_interv == 1:
            return self

        self.set_all_successors()

        # Order of the B_i
        t_ = np.zeros(n_interv, dtype=int)

        (i, B_1) = (1, 0)

        S = [self[B_1]]
        start0, current_end = self[B_1].bounds()

        self[B_1].mark() # Mark the 1st arc
        t_[B_1] = 1

        B_i = B_1

        while current_end < self.mod-1 or (current_end+1) % self.mod < start0:
            B_i = self.successors[B_i]
            S.append(self[B_i])
            self[B_i].mark()

            current_end = self[B_i].end

            t_[B_i] = i + 1
            i += 1

        # At this point we have an initial cover
        k = i

        is_opt = False

        while not is_opt:
            i += 1
            B_i = self.successors[B_i]

            # S.insert(self.arc(B_i), inorder=True)
            S.append(self[B_i])

            if(self[B_i].marked or
               ((i % (k-1) == 1) # A total of k disjoint zones has been found
                and not self[B_i].intersect(self[B_1]))
            ):
                is_opt = True
                if i == t_[B_i] + (k-1): # Current arc Bi is marked and is the same as B[i-k]]
                    indices = range(i-(k-1), i)
                else:
                    indices = range(k)
            else: # B_i is in S or the number of disjoint zones is k, size of the minimal cover=k
                self[B_i].mark()
                t_[B_i] = i

        S = Coverage(*[S[i] for i in indices])

        return S
        
    def get_junctions(self):
        '''
        Assumes perfect coverage
        '''
        data = np.zeros(len(self), dtype=[('source1', '<U64'), ('source2', '<U64'),
                                          ('start', 'uint32'), ('end', 'uint32')])
        prev = self[len(self)-1]
        
        for i in range(len(self)):
            curr = self[i]
            
            data[i] = ('/'.join(prev.data.index),
                       '/'.join(curr.data.index),
                       curr.start,
                       (prev.end + 1) % self.mod)
            prev = curr

        return data

