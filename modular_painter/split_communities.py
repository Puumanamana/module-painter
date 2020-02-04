from pathlib import Path
from tempfile import mkdtemp

import pandas as pd
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from Bio.Blast import NCBIXML

import matplotlib.pyplot as plt

## TO DO: only create database once for each sample

def display_alignment(record, selected, output=None):
    qlen = record.query_length
    n_hits = len(record.alignments)

    fig, ax = plt.subplots()
    ax.minorticks_on()
    ax.xaxis.grid(True, which='both')
    ax.yaxis.grid(False, which='minor')

    for i, hit in enumerate(record.alignments):
        hit_id = hit.hit_def.split(' ')[0]

        ax.hlines(n_hits-i, 0, qlen, color='k', linewidth=.2)

        for hsp in hit.hsps:
            (start, end) = (hsp.query_start, hsp.query_end)
            color = 'g' if hit_id in selected else 'r'
            ax.hlines(n_hits-i, start, end, color=color, linewidth=4, label=hit_id)

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), fontsize=8, ncol=2,
               bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(right=0.7)

    plt.xlim([0, qlen])

    plt.title("Blast hits of {}".format(record.query))

    # if output is not None:
    #     plt.savefig("{}/{}.pdf".format(output, record.query), transparent=True)

def blast_sequences(query, db):
    output = Path('/tmp/cedric/modular_painting_tests/{}_on_{}'.format(query.stem, db.stem))
    
    if output.is_file():
        print('No need to do BLAST, found previous results.')
        return output
    
    # tmpdir = mkdtemp()
    tmpdir = '/tmp/cedric/modular_painting_tests/'

    db_path = '{}/{}'.format(tmpdir, db.stem)
    blast_path = '{}/{}_on_{}'.format(tmpdir, query.stem, db.stem)
    
    make_db = NcbimakeblastdbCommandline(dbtype="nucl", input_file=str(db), out=db_path)
    print(make_db)
    make_db()
    
    blastn_on_db = NcbiblastnCommandline(query=str(query), db=db_path, out=blast_path, outfmt=5, gapopen=10, gapextend=5)
    print(blastn_on_db)
    blastn_on_db()

    return blast_path

def filter_blast_results(blast_file, output, min_id=0.8, min_cov=0.5, min_hsp_len=30, show=False):

    species = []

    for record in NCBIXML.parse(open(blast_file)):
        # 1 record is one blast result from the query sequences
        if record.alignments:
            record_id = record.query.split(' ')[0]
            selected = []

            for hit in record.alignments:
                # 1 alignment is one hit (with potentially multiple HSPs)

                matches = sum(hsp.identities for hsp in hit.hsps)
                covered = sum(hsp.align_length for hsp in hit.hsps)

                pct_identity = matches / covered
                pct_coverage = covered / record.query_length

                if pct_identity > min_id and pct_coverage > min_cov:
                    hit_id = hit.hit_def.split(' ')[0]

                    entry = [[record_id, hit_id, hsp.query_start-1, hsp.query_end-1, hsp.identities]
                             for hsp in hit.hsps if abs(hsp.query_end-hsp.query_start) > min_hsp_len]
                    species += entry
                    selected.append(hit_id)

        if (len(selected) > 10) and show:
            display_alignment(record, selected, output=output)
            plt.show()

    return pd.DataFrame(species, columns=['target', 'source', 'tstart', 'tend', 'identity'])

def summarize_results(species):
    print(species.groupby('target').source.agg(lambda x: len(set(x))))

def split_communities(query, db, output=None, show=False, **blast_filter_prms):

    blast_file = blast_sequences(query, db)
    species = filter_blast_results(blast_file, output, show=show, **blast_filter_prms)
    summarize_results(species)

    return species

if __name__ == '__main__':

    query = Path('../test/population-B_117m_simple.fasta')
    db = Path('../test/population-C_250m_simple.fasta')
    output = Path('/tmp/cedric/modular_painting_tests')
    
    split_communities(query, db, output, 0.9, 0.5)
    
    
