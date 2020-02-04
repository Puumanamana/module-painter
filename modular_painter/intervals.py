import sys
import pandas as pd
import numpy as np

def get_sucessor(data, i, mod):
    '''
    Data must not include any duplicated intervals
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
            
        if data[j, 0] + offset > end_i:
            return (j - 1) % len(data)

class Arc:

    def __init__(self, data, target, modulo):
        self.mod = modulo
        self.target = target
        
        self.start = data.start[0]

        self.end = data.end[0]
        if self.end < self.start:
            self.end += modulo
        
        
        self.data = data[['source', 'identity']]
        self.marked = False # For lee and lee algorithm

    def __repr__(self):
        return "{}: ({}, {}), modulo={}".format(
            '/'.join(self.data.source.unique().tolist()),
            *self.bounds(), self.mod)

    def __eq__(self, other):
        return self.data.to_dict() == other.data.to_dict()

    def __len__(self):
        return self.end - self.start + 1

    def bounds(self):
        return [self.start, self.end]

    def intersect(self, other):
        if self.target != other.target:
            sys.exit("Cannot compare, different target {} vs {}".format(self.target, other.target))

        if ((other.start <= self.start <= other.end)
            or (other.start <= self.end % self.mod <= other.end)
            or (self.start <= other.end % other.mod <= self.end)):
            return True

        return False

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
            df = arc.data
            df['source'] = arc.data.source
            df['start'] = arc.start
            df['end'] = arc.end
            result.append(df)

        return pd.concat(result)

    def __len__(self):
        return len(self.data)

    def __repr__(self):
        locations = ['{}-{}'.format(*bound) for bound in self.get_intervals()]
        message = "Target: {}\n{}".format(self.target, '\n'.join(locations))
        return message

    def __getitem__(self, i):
        if self.dirty:
            self.data = sorted(self.data, key=lambda x: x.bounds())
        return self.data[i]

    def get_intervals(self):
        return np.array([self[i].bounds() for i in range(len(self))])

    def simplify(self):
        '''
        - Remove any interval that is strictly embedded in another
        - Assumes no equal interval boundaries
        '''


        # Remove any duplicate start and keep the last one (larger interval)
        intervals = pd.DataFrame(self.get_intervals(), columns=['start', 'end'])
        intervals = intervals[~intervals.start.duplicated(keep='last')]

        # At this point, no pair of interval share the same start point
        # Reorder the intervals to put the ones that span the "0" point first
        intervals[intervals.end > self.mod] -= self.mod
        intervals.sort_values(by=['start', 'end'], inplace=True)

        # Discard an intervals if its end is not striclty increasing
        cum_max = intervals.end.shift(1).cummax().fillna(intervals.start.min())
        intervals = intervals[intervals.end > cum_max]
        interv_to_keep = {(start, end) if start > 0
                          else (start + self.mod, end + self.mod)
                          for (start, end) in intervals.to_numpy()}

        self.data = [arc for arc in self.data if tuple(arc.bounds()) in interv_to_keep]

    def insert(self, other):
        '''
        Insert an arc to the coverage
        - Does not sort the array by start/end each time
        '''

        if self.target != other.target:
            sys.exit('Cannot insert: different targets ({}, {})'
                     .format(self.target, other.target))

        self.data.append(other)
        self.dirty = True

    def set_all_successors(self):
        arc_boundaries = self.get_intervals()
        self.successors = [get_sucessor(arc_boundaries, i, self.mod) for i in range(len(self))]

    def merge_close_intervals(self, arc_eq_diffs):
        '''
        Group intervals that have about the same boundaries 
        '''

        intervals = pd.DataFrame(self.get_intervals(), columns=['start', 'end'])

        close = np.abs((intervals.shift(1)-intervals)).sum(axis=1) <= arc_eq_diffs

        # groups of similar bounds
        groups = np.zeros(len(close), dtype=int)
        last_id = 0

        for i, same in enumerate(close.to_numpy()[1:], 1):
            if not same:
                last_id += 1
            groups[i] = last_id

        intervals.set_index(groups, inplace=True)
        # Set the new bound by extending the previous arc
        intervals = intervals[intervals.index.duplicated(keep=False)]
        group_boundaries = intervals.groupby(level=0).agg({'start': min, 'end': max})

        # Update the arcs
        for i in range(len(self)):
            if groups[i] in group_boundaries.index:
                self[i].start = group_boundaries.loc[groups[i], 'start']
                self[i].end = group_boundaries.loc[groups[i], 'end']

        self.dirty = True

    def fill_missing_data(self, max_merge_dist):
        i_max = 0
        intervals = self.get_intervals()
        
        current_end = self[0].end

        fillers = []
        filler_cols = ['source', 'start', 'end', 'identity']

        for i, (start_i, end_i) in enumerate(intervals[1:], 1):
            dist = start_i - current_end

            if 0 < dist <= max_merge_dist:
                # extend the previous arc
                self.data[i-1].end += dist
            elif dist > max_merge_dist:
                filler = Arc(
                    pd.DataFrame([['NoCov', current_end, start_i, 0]], columns=filler_cols),
                    self.target, self.mod)
                fillers.append(filler)

            if end_i >= current_end:
                # update the max
                i_max = i
                current_end = end_i

        # Looping condition
        loop_dist = (intervals[0, 0] + self.mod) - current_end

        if loop_dist > max_merge_dist:
            filler = Arc(pd.DataFrame([['NoCov', current_end, intervals[0, 0], 0]], columns=filler_cols),
                         self.target, self.mod)
            fillers.append(filler)

        if 0 < loop_dist <= max_merge_dist:
            self[i_max].end += loop_dist
            
        # Add the fillers
        for filler in fillers:
            self.insert(filler)

        self.dirty = True

        if not self.is_covered():
            print('Something went wrong. Coverage is missing some parts')
            import ipdb;ipdb.set_trace()

    def is_covered(self):
        current_end = self.data[0].end

        for i in range(len(self)):
            (start, end) = self[i].bounds()
            if start > current_end: # There is a hole
                return False
            if end > current_end:
                current_end = end

        # Looping condition: the largest end needs to overlap the smallest start
        if (
                (current_end % self.mod >= self[0].start)
                and (current_end >= self.mod)
        ):
            return True
        return False

    def get_break_points(self):
        '''
        Assume we already have an optimal coverage
        --> a location on the circle is covered by 2 arcs tops
        '''

        return

    def get_minimal_coverage(self):
        '''
        Lee and lee 1983 algorithm
        On a circle-cover minimization problem
        '''

        n_interv = len(self)

        self.set_all_successors()

        # Order of the B_i
        t_ = np.zeros(n_interv, dtype=int)

        (i, B_1) = (1, 0)

        S = [self[B_1]]
        start0, current_end = self[B_1].bounds()

        self[B_1].mark() # Mark the 1st arc
        t_[B_1] = 1

        B_i = B_1

        while (current_end < self.mod) or (current_end % self.mod < start0):

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
        
