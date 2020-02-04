from itertools import chain
from Bio import SeqIO

from intervals import Coverage
from modular_painter.display import donut_display
    
def rotate_genome(genome):
    '''
    Possible approaches: 
    - BLAST on terminase?
    - Align based on one module present in all individuals in the species group
    '''

    return str(genome)
    
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
    for target in species.target.unique():
        print('processing {}'.format(target))
        coverage_t = Coverage.from_pandas(
            species.loc[species.target == target].drop('target', axis=1),
            target,
            genome_lengths[target]
        )

        coverage_t.merge_close_intervals(arc_eq_diffs)
        coverage_t.simplify()
        coverage_t.fill_missing_data(min_module_size)
        optimal_coverage = coverage_t.get_minimal_coverage()
        result[target] = optimal_coverage
        optimal_coverage.get_break_points()

        donut_display(optimal_coverage, '/tmp/cedric/modular_painting_tests/painting_{}.html'.format(target), circular=False)
        
    return result
