from itertools import compress
from pathlib import Path
from tempfile import mkdtemp

import pandas as pd
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from Bio.Blast import NCBIXML

import matplotlib.pyplot as plt

from display import get_palette

## TO DO: only create database once for each sample



def split_communities(query, db, output=None, show=False, **blast_filter_prms):

    blast_file = blast_sequences(query, db)
    species = filter_blast_results(blast_file, output, show=show, **blast_filter_prms)

    return species

def blast_sequences(query, db, cache=True):
    output = Path('/tmp/cedric/modular_painting_tests/{}_on_{}'.format(query.stem, db.stem))

    if output.is_file() and cache:
        print('No need to do BLAST, found previous results.')
        return output

    # tmpdir = mkdtemp()
    tmpdir = '/tmp/cedric/modular_painting_tests/'

    db_path = '{}/{}'.format(tmpdir, db.stem)
    blast_path = '{}/{}_on_{}'.format(tmpdir, query.stem, db.stem)

    make_db = NcbimakeblastdbCommandline(dbtype="nucl", input_file=str(db), out=db_path)
    print(make_db)
    make_db()

    blastn_on_db = NcbiblastnCommandline(query=str(query), db=db_path, out=blast_path,
                                         outfmt=5, gapopen=5, gapextend=2)
    print(blastn_on_db)
    blastn_on_db()

    return blast_path

def filter_blast_results(blast_file, output,
                         min_id=0.8,
                         min_cov=0.5,
                         min_module_size=30,
                         show=False):

    species = {}
    selected = []

    for record in NCBIXML.parse(open(blast_file)):
        # 1 record is one blast result from the query sequences
        if record.alignments:
            record_id = record.query.split(' ')[0]
            species[record_id] = []

            for hit in record.alignments:
                # 1 alignment is one hit (with potentially multiple HSPs)

                matches = sum(hsp.identities for hsp in hit.hsps)
                covered = sum(hsp.align_length for hsp in hit.hsps)

                pct_identity = matches / covered
                pct_coverage = covered / record.query_length

                if pct_identity > min_id and pct_coverage > min_cov:
                    hit_id = hit.hit_def.split(' ')[0]
                    hsps = merge_contiguous_hsps(hit.hsps, record.query_length, min_module_size)

                    entry = [[hit_id, hsp.query_start, hsp.query_end, hsp.identities]
                             for hsp in hsps]
                    species[record_id] += entry
                    # species += entry
                    selected.append(hit_id)

            species[record_id] = pd.DataFrame(species[record_id], columns=['source', 'start', 'end', 'identity'])

        if (len(selected) > 10) and show:
            display_alignment(record, selected, output=output)
            plt.show()
    return species

def merge_contiguous_hsps(hsps_in, query_len, min_module_size):
    hsps = sorted(hsps_in, key=lambda x: (x.query_start, x.query_end))
    i_max = 0

    to_keep = [True] * len(hsps)

    for i, hsp in enumerate(hsps[1:], 1):
        hsp.query_start -= 1
        hsp.query_end -= 1

        qmax = hsps[i_max].query_end

        if hsp.query_start <= qmax+1+min_module_size:
            to_keep[i] = False

            if hsp.query_end > qmax:
                hsps[i_max].query_end = hsp.query_end
                hsps[i_max].identities += hsp.query_end - qmax
        else:
            i_max = i

    return compress(hsps, to_keep)


def display_alignment(record, selected, output=None):
    qlen = record.query_length
    n_hits = len(record.alignments)

    fig, ax = plt.subplots()
    ax.minorticks_on()
    ax.xaxis.grid(True, which='both')
    ax.yaxis.grid(False, which='minor')

    palette = get_palette(len(record.alignments))
    for i, hit in enumerate(record.alignments):
        hit_id = hit.hit_def.split()[0]

        ax.hlines(n_hits-i, 0, qlen, color='k', linewidth=.2)

        for hsp in hit.hsps:
            (start, end) = (hsp.query_start, hsp.query_end)
            # color = 'g' if hit_id in selected else 'r'
            color = palette[i]
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
