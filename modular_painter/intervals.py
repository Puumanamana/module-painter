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

    def __init__(self, data, target, modulo, check_bounds=True):
        self.mod = modulo
        self.target = target

        if check_bounds:
            data.loc[data.tend < data.tstart, 'tend'] += modulo
        self.data = data
        self.marked = False # For lee and lee algorithm

    def __repr__(self):
        return "{}: ({}, {}), modulo={}".format(
            '/'.join(self.data.source.unique().tolist()),
            self.start(), self.end(),
            self.mod)

    def __eq__(self, other):
        return self.data.to_dict() == other.data.to_dict()

    def __len__(self):
        return self.end() - self.start() + 1

    def start(self):
        return self.data.tstart.iloc[0]

    def end(self):
        return self.data.tend.iloc[0]

    def intersect(self, other):
        if self.target != other.target:
            sys.exit("Cannot compare, different target {} vs {}".format(self.target, other.target))

        if ((other.start() <= self.start() <= other.end())
            or (other.start() <= self.end() % self.mod <= other.end())
            or (self.start() <= other.end() % other.mod <= self.end())):
            return True

        return False

    def mark(self):
        self.marked = True
    
class Coverage:

    def __init__(self, *arcs):
        if arcs:
            self.target = arcs[0].target
            self.mod = arcs[0].mod
            self.data = pd.concat([arc.data for arc in arcs])
        self.dirty = True # Keep track of whether the location are sorted
        self._intervals = None

    @property
    def intervals(self):
        if self.dirty:
            self._intervals = (
                self.data[['tstart', 'tend']]
                .drop_duplicates()
                .sort_values(by=['tstart', 'tend'])
                .values
            )
            self.dirty = False
        return self._intervals

    @classmethod
    def from_pandas(cls, data, target, modulo=-1):
        cov = cls()
        if modulo < 0:
            modulo = data[['tstart', 'tend']].max().max()

        cov.target = target
        cov.mod = modulo
        cov.data = data

        cov.data = data.copy()
        cov.data.loc[data.tend < data.tstart, 'tend'] += modulo

        return cov

    def __len__(self):
        return len(self.intervals)

    def __repr__(self):
        locations = ['{}-{}'.format(*loc) for loc in self.intervals]
        message = "Target: {}\n{}".format(self.target, '\n'.join(locations))
        return message

    def arc(self, i):
        # !!!! set_index(start,end) !!!
        start, end = self.intervals[i]
        interv = self.data[
            (self.data.tstart == start)
            & (self.data.tend == end)].copy()
        
        return Arc(interv, self.target, self.mod)

    def filter_interv(self, intervals):
        to_keep = [i for i, x in enumerate(self.data[['tstart', 'tend']].to_numpy())
                   if tuple(x) in intervals]

        self.data = self.data.iloc[to_keep]
        self.dirty = True

    def simplify(self):
        '''
        - Remove any interval that is strictly embedded in another
        - Assumes no equal interval boundaries
        '''

        # Remove any duplicate start and keep the last one (larger interval)
        intervals = pd.DataFrame(self.intervals, columns=['tstart', 'tend'])
        intervals = intervals[~intervals.tstart.duplicated(keep='last')]

        # At this point, no pair of interval share the same start point
        # Reorder the intervals to put the ones that span the "0" point first
        intervals[intervals.tend > self.mod] -= self.mod
        intervals.sort_values(by=['tstart', 'tend'], inplace=True)

        # Discard an intervals if its end is not striclty increasing
        intervals = intervals[intervals.tend > intervals.tend.shift(1).cummax().fillna(intervals.tstart.min())]
        interv_to_keep = {(start, end) if start > 0
                          else (start + self.mod, end + self.mod)
                          for (start, end) in intervals.values}

        self.filter_interv(interv_to_keep)

    def insert(self, other, inorder=False):
        '''
        Insert an arc to the coverage
        - Does not sort the array by start/end each time
        '''

        if self.target != other.target:
            sys.exit('Cannot insert: different targets ({}, {})'
                     .format(self.target, other.target))

        self.data = pd.concat([self.data, other.data])

        if not inorder:
            self.dirty = True
        else:
            if self.intervals is None:
                _ = self.intervals
            self._intervals = np.vstack([
                self._intervals,
                other.data[['tstart', 'tend']].drop_duplicates().to_numpy()
            ])

    def set_all_successors(self):
        self.successors = [get_sucessor(self.intervals, i, self.mod) for i in range(len(self))]

    def merge_close_intervals(self, arc_eq_diffs):
        '''
        Group intervals that have about the same boundaries 
        '''

        intervals = pd.DataFrame(self.intervals, columns=['tstart', 'tend'])
        close = np.abs((intervals.shift(1)-intervals)).sum(axis=1) <= arc_eq_diffs
        
        groups = np.zeros(len(close), dtype=int)
        last_id = 0

        for i, same in enumerate(close.to_numpy()[1:], 1):
            if not same:
                last_id += 1
            groups[i] = last_id

        intervals.set_index(groups, inplace=True)
        intervals = intervals[intervals.index.duplicated(keep=False)]
        intervals = intervals.groupby(['tstart', 'tend']).agg({'tstart': max, 'tend': min})

        
        for idx, row in zip(self.data.index, self.data.to_numpy()):
            position = (row[0], row[1])
            if position in intervals.index:
                entry = intervals.loc[position].tolist() + [row.identity-arc_eq_diffs]
                self.data.loc[idx, ['tstart', 'tend', 'identity']] = entry

        self.dirty = True

    def fill_missing_data(self, max_merge_dist):
        i_max = 0
        current_end = self.intervals[0, 1]

        fillers = []
        filler_cols = ['source', 'tstart', 'tend', 'identity']

        for i, (start_i, end_i) in enumerate(self.intervals[1:], 1):
            dist = start_i - current_end

            if 0 < dist <= max_merge_dist:
                # extend the previous arc
                prev_start, prev_end = self.intervals[i-1]
                self.data.loc[(self.data.tstart == prev_start)
                              & (self.data.tend == prev_end), 'tend'] += dist
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
        loop_dist = (self.intervals[0, 0] + self.mod) - current_end

        if loop_dist > max_merge_dist:
            filler = Arc(pd.DataFrame([['NoCov', current_end, self.intervals[0, 0], 0]], columns=filler_cols),
                         self.target, self.mod)
            fillers.append(filler)

        if 0 < loop_dist <= max_merge_dist:
            start_max, end_max = self.intervals[i_max]
            self.data.loc[(self.data.tstart == start_max)
                          & (self.data.tend == end_max), 'tend'] += loop_dist
            self.dirty = True

        # Add the fillers
        for filler in fillers:
            self.insert(filler)

        if not self.is_covered():
            print('Something went wrong. Coverage is missing some parts')
            import ipdb;ipdb.set_trace()

    def is_covered(self):
        current_end = self.intervals[0, 0]

        for (start, end) in self.intervals:
            if start > current_end: # There is a hole
                return False
            if end > current_end:
                current_end = end

        # Looping condition: the largest end needs to overlap the smallest start
        if (
                (current_end % self.mod >= self.intervals[0, 0])
                and (current_end >= self.mod)
        ):
            return True
        return False

    def get_break_points(self):
        '''
        Assume we already have an optimal coverage
        --> a location on the circle is covered by 2 arcs tops
        '''

        # tmp = self.data.groupby(['tstart', 'tend'], as_index=False).agg(list)

        # breaks = np.vstack([
        #     np.roll(tmp[['tend', 'source']]),
        #     tmp[['tstart', 'source']].to_numpy()
        # ]).T

        # breakpoints = np.vstack([
        #     np.roll(self.intervals[:, 0], -1),
        #     self.intervals[:, 1],
        # ]).T

        # Looping condition
        # loops = breakpoints[breakpoints[:, 1] - breakpoints[:, 0] > self.mod]
        # loops[:, 1] = loops[:, 1] % self.mod
        return

    def get_minimal_coverage(self):
        '''
        Lee and lee 1983 algorithm
        On a circle-cover minimization problem
        '''

        n_interv = len(self)

        tmp = self.data.set_index(['tstart', 'tend']).sort_index()
        arc_list = [Arc(tmp.loc[s, e].reset_index(), self.target, self.mod, check_bounds=False)
                    for (s, e) in tmp.index.unique()]
        self.set_all_successors()

        # Order of the B_i
        t_ = np.zeros(n_interv, dtype=int)

        (i, B_1) = (1, 0)

        S = [arc_list[B_1]]
        start0, current_end = self.intervals[B_1]

        arc_list[B_1].mark() # Mark the 1st arc
        t_[B_1] = 1

        B_i = B_1

        while (current_end < self.mod) or (current_end % self.mod < start0):

            B_i = self.successors[B_i]
            S.append(arc_list[B_i])
            arc_list[B_i].mark()

            current_end = arc_list[B_i].end()

            t_[B_i] = i + 1
            i += 1
        
        # S = Coverage(*S)

        # print(Coverage(*S))
        # At this point we have an initial cover
        k = i

        is_opt = False

        while not is_opt:
            i += 1
            B_i = self.successors[B_i]

            # S.insert(self.arc(B_i), inorder=True)
            S.append(arc_list[B_i])

            if(arc_list[B_i].marked or
               ((i % (k-1) == 1) # A total of k disjoint zones has been found
                and not arc_list[B_i].intersect(arc_list[B_1]))
            ):
                is_opt = True
                if i == t_[B_i] + (k-1): # Current arc Bi is marked and is the same as B[i-k]]
                    indices = range(i-(k-1), i)
                else:
                    indices = range(k)
            else: # B_i is in S or the number of disjoint zones is k, size of the minimal cover=k
                arc_list[B_i].mark()
                t_[B_i] = i

        S = Coverage(*[S[i] for i in indices])

        return S
        
