from Bio import SeqIO

def rotate_genome(sequence):
    '''
    Possible approaches: 
    - BLAST on terminase?
    - Align based on one module present in all individuals in the species group
    '''

def paint_sequence(p, q):
    '''
    Paint phage q with phage p:
    - Simple blast?
    - Result shape = array/dataframe with 5 columns: [s1, e1, s2, e2, identity]
      sorted by {1: s2, 2: e2}
    '''

def clean_results(painting):
    '''
    - Filter the whole painting by removing any alignment that is strictly embedded in another
    - Fill holes in painting
    '''

    # Remove any duplicate start that does not have the same end, keep the last one (larger interval)
    no_dup_start = ~painting.s2.duplicated(keep='last') | painting.duplicated(keep=False)

    painting_filt = painting[no_dup_start]
    # At this point, painting_filt has strictly start values, with equality only if the end is also equal.
    
    painting_prev = painting_filt.shift(1)
    cum_max = painting_prev.e2.cummax()

    to_keep = (
        (painting_filt.e2 > cum_max) # end position is strictly bigger compared to before (and the start can't be the same after)
        | (painting_filt.e2 == cum_max) # or it's a duplicate maximal entry (that we keep)
        & ((painting_filt == painting_prev).sum(axis=1) == 2)
    )

    to_keep.iloc[0] = True

    return painting_filt[to_keep]

def lee_and_lee(painting):
    '''
    Run Lee and Lee algorithm to compute the optimal cover and the alternative parents
    '''

def paint_all(fasta, species):
    '''
    Paint all species groups in population
    '''

    genomes = {sequence.id: rotate_genome(sequence.seq)
               for sequence in SeqIO.parse(fasta, 'fasta')}

    paintings = {
        subject: {
            painter.id: paint_sequence(genomes[painter], genomes[subject])
            for painter in species[subject.id]
        }
        for subject in species
    }

    result = {}
    for subject, painting in paintings.items():
        painting_trimmed = clean_results(painting)
        optimal_coverage = lee_and_lee(painting_trimmed)

        result[subject] = optimal_coverage

    return optimal_coverage
