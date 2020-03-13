from Bio import SeqIO
import pandas as pd

from intervals import Coverage
from modular_painter.display import display_single_genome
    
def paint_all(fastas, species, arc_eq_diffs=None, min_module_size=None):
    '''
    Paint all species groups in population
    '''

    genome_lengths = dict([(sequence.id, len(sequence.seq)) for fasta in fastas
                           for sequence in SeqIO.parse(fasta, 'fasta')])

    result = {}
    for target in species:

        print('processing {}'.format(target))
        coverage_t = Coverage.from_pandas(species[target], target, genome_lengths[target])
        coverage_t.extend_close_intervals(arc_eq_diffs)
        coverage_t.merge_equal_intervals()
        coverage_t.fill_or_extend(min_module_size)
        coverage_t.circularize()

        result[target] = coverage_t.get_minimal_coverage()

        # display_single_genome(result[target], '/tmp/cedric/modular_painting_tests/painting_{}.html'
        #                       .format(target), circular=False)
        
    return result
