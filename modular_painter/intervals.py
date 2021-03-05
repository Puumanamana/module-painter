from itertools import compress
import sys
import pandas as pd
import numpy as np

def get_successor(data, i, mod):
    '''
    - Data must not include any duplicated intervals
    - The while interval must be covered
    '''

    offset = 0
    end_i = data.end[i]

    not_found = True
    j = i

    while not_found:
        j += 1

        if j >= len(data):
            j = 0
            offset = mod

        if data.start[j] + offset > end_i + 1:
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

    if overlap:
        return ('do_nothing', arc_next)

    if close:
        return ('extend', arc_next)

    # The gap is too big --> we fill it with a fake arc
    filler = pd.DataFrame([['NoCov', 0]], columns=['source', 'identity'])

    return ('fill', Arc(arc_prev.end+1, arc_next.start-1, arc_prev.mod, target=arc_prev.target, data=filler))

class Arc:

    def __init__(self, start, end, modulo, target=None, data=None):
        self.mod = modulo
        self.target = target
        self.start = start
        self.end = end
        self.fix_bounds()

        if data is not None:
            self.data = data.set_index('source')
        self.marked = False # For lee and lee algorithm

    @classmethod
    def from_pandas(cls, data, modulo, target=None):
        '''
        start/end as columns in data
        '''
        return cls(data.start[0], data.end[0], modulo, target=target, data=data.drop(['start', 'end'], axis=1))

    def __repr__(self):
        return "{}: ({}, {}), modulo={}".format(
            '/'.join(self.data.index),
            *self.bounds(), self.mod)

    def __eq__(self, other):
        return self.data.to_dict() == other.data.to_dict()

    def __len__(self):
        return self.end - self.start + 1

    def __getitem__(self, key):
        if key in self.__dict__:
            return getattr(self, key)
        else:
            return self.data[key]

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
        2 arcs intersect if they have no gap inbetween
        Note: [0, 1] and [2,4] intersect since we work with discrete intervals
        '''
        if self.target != other.target or self.mod != other.mod:
            sys.exit("Cannot compare, different targets {} vs {}".format(self.target, other.target))

        if self.end >= self.mod:
            arc1_1 = Arc(self.start, self.mod-1, self.mod, self.target)
            arc1_2 = Arc(0, self.end % self.mod, self.mod, self.target)

            return arc1_1.intersect(other) or arc1_2.intersect(other)

        if other.end >= self.mod:
            arc2_1 = Arc(other.start, other.mod-1, other.mod, self.target)
            arc2_2 = Arc(0, other.end % other.mod, other.mod, self.target)

            return self.intersect(arc2_1) or self.intersect(arc2_2)

        if self.start > other.start:
            return other.intersect(self)

        # we only compare arcs that do not cross 0 and where "self" starts before "other"
        return (
            self.start <= other.start <= self.end + 1
            or (other.end + 1) % other.mod == self.start
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

        df = data.set_index(['start', 'end']).sort_index()

        arcs = [Arc(*bounds, modulo, data=df.loc[bounds], target=target) for bounds in df.index.unique()]

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

        locations = ['{:3} - {:>9}{:>9}: {}'.format(i, arc.start, arc.end, '/'.join(arc.data.index))
                     for i, arc in enumerate(self.data)]
        message = "Target: {} (L={})\n{}".format(self.target, self.mod, '\n'.join(locations))
        return message

    def __getitem__(self, i):
        if self.dirty:
            self.data = sorted(self.data, key=lambda x: x.bounds())
            self.dirty = False
        return self.data[i]

    def get(self, *keys):
        values = [[arc[k] for k in keys] for arc in self.data]

        return pd.DataFrame(values, columns=keys)

    def change_parent(self, start, end, p1, p2):

        bk = Arc(start, end, self.data[0].mod, target=self.data[0].target)

        for (arc1, arc2) in zip(self.data, self.data[1:]+[self.data[0]]):
            if arc1.intersect(bk):
                if p1 in arc1.data.index and p2 in arc2.data.index:
                    arc1.data = arc1.data.loc[[p1]]
                    arc2.data = arc2.data.loc[[p2]]
                    break
                elif p2 in arc1.data.index and p1 in arc2.data.index:
                    arc1.data = arc1.data.loc[[p2]]  
                    arc2.data = arc2.data.loc[[p1]]
                    break
        else:
            import ipdb;ipdb.set_trace()

        if bk.target=='X':
            print(start, end, p1, p2)
            import ipdb;ipdb.set_trace()
                
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

    def extend_close_intervals(self, arc_eq_diffs):
        '''
        "Sticky boundaries": extend start/end that are close to another interval
        '''

        transf = {'start': min, 'end': max}

        for key, func in transf.items():
            feature = self.get(key)[key].sort_values()

            same_feature = (feature - np.roll(feature, 1) > arc_eq_diffs).cumsum().rename('group')
            grouped = (same_feature.reset_index()
                       .groupby('group').filter(lambda x: len(x)>1)
                       .groupby('group')['index'].agg(list))

            for indices in grouped.values:
                new = func(feature[indices])
                for i in indices:
                    setattr(self.data[i], key, new)

    def merge_equal_intervals(self):
        '''
        Group intervals that have the same boundaries
        '''

        intervals = self.get('start', 'end').reset_index()

        # is_far has an incrementing index that stay constant
        # when both boundaries_dist <= arc_eq_diffs
        grouped = intervals.groupby(['start', 'end'])['index'].agg(list)

        final_indices = np.ones(len(intervals)).astype(bool)

        for indices in grouped.values:
            i0 = indices.pop()
            self[i0].data = pd.concat([self[j].data for j in indices+[i0]])
            # Remove all indices-i0
            final_indices[indices] = False

        self.data = list(compress(self.data, final_indices))

    def fill_or_extend(self, max_extend_dist):
        '''
        '''

        (i_max, i) = (0, 0)
        to_discard = set()

        while True:
            i += 1
            if i < len(self):
                arc = self[i]
            else:
                break

            # Decide whether we should merge or fill gap between the 2 arcs or discard an arc
            outcome, new_arc = compare_arc_pair(self[i_max], arc, max_extend_dist)

            if outcome == 'discard_prev':
                to_discard.add(i_max % len(self))

            elif outcome == 'discard':
                if i % len(self) != i_max:
                    to_discard.add(i % len(self))

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

        if not self[0].intersect(self[i_max]):
            start0, end0 = self[0].bounds()
            startf, endf = self[i_max].bounds()

            if (start0 + self.mod) - endf <= max_extend_dist:
                self[i_max].end = self.mod + (start0 - 1)
            else:
                filler = pd.DataFrame([['NoCov', 0]], columns=['source', 'identity'])
                self.data.append(Arc(endf+1, start0-1 % self.mod, self.mod, data=filler, target=self.target))

        self.data = [arc for (i, arc) in enumerate(self.data) if i not in to_discard]

        if not self.is_covered():
            print('Something went wrong. Coverage is missing some parts')
            import ipdb;ipdb.set_trace()

    def set_all_successors(self):
        arc_boundaries = self.get('start', 'end')
        self.successors = [get_successor(arc_boundaries, i, self.mod) for i in range(len(self))]

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

    def circularize(self):
        common_sources = np.intersect1d(self[0].data.index, self[-1].data.index)
        if len(common_sources) > 0 and len(self) > 1:
            self[-1].end = self.mod + self[0].end
            self[-1].data = self[0].data.loc[common_sources]
            self.data.pop(0)

    def get_junctions(self):
        '''
        Assumes perfect coverage
        '''

        if len(self) == 1:
            return pd.DataFrame([])

        data = np.zeros(len(self), dtype=[('parents', '<U128'),
                                          ('start', 'uint32'),
                                          ('end', 'uint32')])
        prev = self[len(self)-1]

        for i in range(len(self)):
            curr = self[i]

            (p1, p2) = sorted(('/'.join(sorted(prev.data.index)),
                               '/'.join(sorted(curr.data.index))))

            data[i] = ("{} <-> {}".format(p1, p2),
                       curr.start,
                       (prev.end+1) % self.mod)
            prev = curr

        df = pd.DataFrame(data)
        df['target'] = self.target

        return df
