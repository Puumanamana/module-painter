from itertools import chain
import pandas as pd
import numpy as np
from Bio import SeqIO

from display import donut_display
    
def rotate_genome(genome):
    '''
    Possible approaches: 
    - BLAST on terminase?
    - Align based on one module present in all individuals in the species group
    '''

    return str(genome)

def trim_results(painting):
    '''
    Filter the whole painting by removing any alignment that is strictly embedded in another
    - Starts with a dataframe sorted by start position
    - Select unique start position entries with the largest size
    '''

    painting_dedup = (painting
        .groupby(['tstart', 'tend'], as_index=False)
        .agg(list)
    )

    # Remove any duplicate start and keep the last one (larger interval)
    # At this point, painting_filt has strictly start
    painting_filt = painting_dedup[~painting_dedup.tstart.duplicated(keep='last')]

    painting_filt = painting_filt[
        (painting_filt.tend > painting_filt.tend.shift(1).cummax())
        | (painting_filt.tend == painting_filt.tend.min())
    ]

    return painting_filt

def fill_missing_data(painting, genome_length, max_merge_dist=3):
    '''
    Make up missing pieces for the lee_and_lee algorithm
    '''

    fill = []
    prev_end = 1

    for idx in painting.index:
        dist = painting.tstart[idx] - prev_end

        if 0 < dist <= max_merge_dist:
            painting.loc[idx, 'tstart'] -= dist
        elif dist > max_merge_dist:
            fill.append([prev_end, painting.tstart[idx], ['NC'], [1]])
        prev_end = painting.tend[idx]

    fill = pd.DataFrame(fill, columns=['tstart', 'tend', 'source', 'identity'])

    painting = pd.concat([painting, fill]).sort_values(by=['tstart', 'tend'])
    painting.index = range(len(painting))

    return painting

def get_successors(i, intervals):

    starts = intervals.tstart.copy()
    end_i = intervals.loc[i, 'tend']

    not_found = True
    j = i

    while not_found:
        j += 1
        
        if j >= len(intervals):
            j = 0
            starts += intervals.tend.max()
            
        if starts[j] > end_i:
            return (j - 1) % len(intervals)

def lee_and_lee(painting):
    '''
    Run Lee and Lee algorithm to compute the optimal cover and the alternative parents
    {painting} is a pandas DataFrame with 4 columns: 
      - [tstart, tend] are the HSP boundaries
      - source is the name of the HSP source ("NC" corresponds to no coverage)
      - identity is the pct_identity for this alignement
    '''

    successors = [get_successors(i, painting) for i in painting.index]

    # print(pd.Series(list(zip(painting.index, successors))))

    # fig, ax = plt.subplots()
    # for i, (start, end) in enumerate(painting[['tstart', 'tend']].values):
    #     ax.hlines(i, start, end, color='r', linewidth=3)
    # plt.show()

    # To change since genomes are circular
    genome_len = painting.tend.max()
    
    n_interv = len(painting)

    # All arcs are unmarked at first
    r_ = np.arange(n_interv)
    t_ = np.zeros(n_interv)

    i = 1
    B_i = 0
    S = [B_i]

    # Mark the first arc
    r_[0] += n_interv
    t_[0] = 1

    covered = False

    while not covered:
        i += 1

        B_i = successors[B_i]
        S.append(B_i)
        r_[B_i] += n_interv
        t_[B_i] = i

        # To change since the genome is circular
        covered = (painting.tend[B_i] - painting.tstart[S[0]] + 1) >= genome_len

    k = i
    # At this point we have an initial cover of size k

    done = False
    while not done:
        i += 1
        B_i = successors[B_i]

        inter_B1 = (painting.tstart[S[0]] < painting.tend[B_i] < painting.tend[S[0]]) \
            or (painting.tstart[S[0]] < painting.tstart[B_i] < painting.tend[S[0]])

        if((r_[B_i] > n_interv) # r_[B_i] > n_interv means B_i has already been marked
           or (i % (k-1) == 1) and not inter_B1): # A total of k disjoint zones has been found
            done = True

            if i == t_[B_i] + (k-1): # Current arc Bi is marked and is the same as B[i-(k+1)]]
                S = S[i-(k-1):i]
        else: # B_i is in S or the number of disjoint zones is k, the size of the minimal cover is k
            r_[B_i] = r_[B_i] + n_interv
            t_[B_i] = i

    return painting.loc[S]
    
def paint_all(fastas, species, max_merge_dist=None):
    '''
    Paint all species groups in population
    '''

    genome_lengths = dict(chain(*[
        [(sequence.id, len(sequence.seq))
         for sequence in SeqIO.parse(fasta, 'fasta')]
        for fasta in fastas
    ]))

    result = {}
    for target in species.target.unique():
        print('processing {}'.format(target))
        subset = species.loc[species.target == target].drop('target', axis=1)
        painting_trimmed = trim_results(subset)
        painting_full = fill_missing_data(painting_trimmed, genome_lengths[target], max_merge_dist=max_merge_dist)
        optimal_coverage = lee_and_lee(painting_full)

        result[target] = optimal_coverage

        donut_display(optimal_coverage.copy(), '/tmp/cedric/modular_painting_tests/painting_{}.html'.format(target), circular=False)
        
    return result
