from itertools import chain
from Bio import SeqIO
import pandas as pd

from intervals import Coverage
from modular_painter.display import donut_display
    
def paint_all(fastas, species, arc_eq_diffs=None, min_module_size=None):
    '''
    Paint all species groups in population
    '''

    genome_lengths = dict(chain(*[
        [(sequence.id, len(sequence.seq))
         for sequence in SeqIO.parse(fasta, 'fasta')]
        for fasta in fastas
    ]))

    result = {}
    for target in species:
        print('processing {}'.format(target))
        coverage_t = Coverage.from_pandas(species[target], target, genome_lengths[target])
        coverage_t.merge_close_intervals(arc_eq_diffs)
        coverage_t.fill_or_extend(min_module_size)
        optimal_coverage = coverage_t.get_minimal_coverage()
        result[target] = optimal_coverage

        print(optimal_coverage)
        print(pd.DataFrame(optimal_coverage.get_junctions()).head(10))

        donut_display(optimal_coverage, '/tmp/cedric/modular_painting_tests/painting_{}.html'
                      .format(target), circular=False)
        
    return result
