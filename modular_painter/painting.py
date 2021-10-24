import numpy as np
from Bio import SeqIO
from sklearn.neighbors import BallTree
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import pdist, squareform

from coverage import Coverage


def paint_genomes(fastas, species, arc_eq_diffs=None, min_module_size=None):
    '''
    Paint all species groups in population
    For each phage, do:
      1) Sticky boundaries
      2) Merge intervals with the same boundaries
      3) Circularize coverage
      4) First simplification of the coverage using embedding rule: 
         - if I C J, we remove J
         - if I = J, we combine I and J
      5) Extend consecutive intervals if they're close, fill holes with putative parents
      6) Compute minimum coverage with Lee & Lee
    '''

    genome_lengths = dict([(sequence.id, len(sequence.seq)) for fasta in fastas
                           for sequence in SeqIO.parse(fasta, 'fasta')])

    result = {}
    for target, itv in species.items():
        print(f'===== processing {target} (L={genome_lengths[target]}) =====')
        itv = coalesce_interval_boundaries(itv, arc_eq_diffs)
        coverage = Coverage.from_pandas(itv, target)
        coverage.fill_gaps(min_module_size)
        coverage.remove_embedded()
        result[target] = coverage.get_minimal_coverage()

    return result

def coalesce_interval_boundaries(itv, arc_eq_diffs):
    '''
    Cluster intervals using agglomerative clustering, such that:
    - the distance between 2 intervals (s1, t1) and (s2, t2) is
       dist( (s1, t1), (s2, t2) ) = max(|s2-s1|, |t2-t1|)
    - distance between 2 clusters is the max of all pair of data points (complete linkage)
    Then merges intervals' parents if they share the same boundaries
    '''
    dist_fn = lambda x,y: np.max(np.abs(x-y))
    dist_mat = squareform(pdist(itv[['start', 'end']], metric=dist_fn))

    model = AgglomerativeClustering(
        linkage='complete',
        distance_threshold=arc_eq_diffs,
        affinity='precomputed',
        n_clusters=None
    ).fit(dist_mat)

    itv.start = itv.groupby(model.labels_).start.transform(min)
    itv.end = itv.groupby(model.labels_).end.transform(max)

    # Combine intervals with the same boundaries
    itv = itv.groupby(['start', 'end']).parent.agg(set).sort_index()

    return itv

        
