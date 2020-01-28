from tempfile import mkdtemp
from collections import defaultdict

from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from Bio.Blast import NCBIXML

import matplotlib.pyplot as plt

## TO DO: only create database once for each sample

def display_alignment(record, selected, display=False, output=None):
    qlen = record.query_length
    n_hits = len(record.alignments)

    fig, ax = plt.subplots()
    ax.minorticks_on()
    ax.xaxis.grid(True, which='both')
    ax.yaxis.grid(False, which='minor')

    for i, hit in enumerate(record.alignments):
        ax.hlines(n_hits-i, 0, qlen, color='k', linewidth=.2)

        for hsp in hit.hsps:
            (start, end) = (hsp.query_start, hsp.query_end)
            color = 'g' if hit.hit_def in selected else 'r'
            ax.hlines(n_hits-i, start, end, color=color, linewidth=4, label=hit.hit_def)

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), fontsize=8, ncol=2,
               bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(right=0.7)

    plt.xlim([0, qlen])

    plt.title("Blast hits of {}".format(record.query))

    if display:
        plt.show()

    if output is not None:
        plt.savefig("{}/{}.pdf".format(output, record.query), transparent=True)

def blast_sequences(query, db):
    if True:
        return '/tmp/cedric/modular_painting_tests/{}_on_{}'.format(query.stem, db.stem)
    
    tmpdir = mkdtemp()

    db_path = '{}/{}'.format(tmpdir, db.stem)
    blast_path = '{}/{}_on_{}'.format(tmpdir, query.stem, db.stem)
    
    make_db = NcbimakeblastdbCommandline(dbtype="nucl", input_file=str(db), out=db_path)
    print(make_db)
    make_db()
    
    blastn_on_db = NcbiblastnCommandline(query=str(query), db=db_path, out=blast_path, outfmt=5)
    print(blastn_on_db)
    blastn_on_db()

    return blast_path

def filter_blast_results(blast_file, output, min_id=0.9, min_cov=0.5):

    species = defaultdict(lambda: [])

    for record in NCBIXML.parse(open(blast_file)):
        # 1 record is one blast result from the query sequences
        if record.alignments:            
            for hit in record.alignments:
                # 1 alignment is one hit (with potentially multiple HSPs)

                matches = sum(hsp.identities for hsp in hit.hsps)
                covered = sum(hsp.align_length for hsp in hit.hsps)

                pct_identity = matches / covered
                pct_coverage = covered / record.query_length

                # print("{} on {}: Coverage={}, Identity={}".format(record.query, hit.hit_def, pct_coverage, pct_identity))
                if pct_coverage > 1:
                    boundaries = [(hsp.sbjct_start, hsp.sbjct_end) for hsp in hit.hsps]
                    print('\n'.join(map(str, boundaries)))
                    import ipdb;ipdb.set_trace()

                if pct_identity > min_id and pct_coverage > min_cov:
                    species[record.query].append(hit.hit_def)
        if record.query in species:
            if len(species[record.query]) > 2:
                display_alignment(record, species[record.query], output=output)
    return species

def summarize_results(species):
    distribution = [len(sp) for sp in species.values()]

    print(distribution)

def check_modules_order():
    '''
    Make sure the order of the modules is the same in all the individuals in the species group
    '''

    print("Module order not implemented yet!")
    print("NB: coverage can be > 1 (overlapping HSPs) --> BLAST with non maximal coverage?")

def split_communities(query, db, output=None, min_id=None, min_cov=None):

    blast_file = blast_sequences(query, db)
    species = filter_blast_results(blast_file, output, min_id=min_id, min_cov=min_cov)
    summarize_results(species)

    check_modules_order()

    return species


if __name__ == '__main__':

    from pathlib import Path

    db = Path('../test/population-B_117m_simple.fasta')
    query = Path('../test/population-C_250m_simple.fasta')
    output = Path('/tmp/cedric/modular_painting_tests')
    
    split_communities(query, db, output, 0.9, 0.5)
    
    
