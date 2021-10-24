from collections import defaultdict
from pathlib import Path
from tempfile import mkdtemp

import pandas as pd
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO

import matplotlib.pyplot as plt

from modular_painter.display import get_palette
from modular_painter.arc import Arc
from modular_painter.coverage import Coverage
from modular_painter.coverage_util import merge_contiguous_blast_hits

## TO DO: only create database once for each sample

def get_raw_coverage(query, db, output=None, **blast_filter_prms):

    blast_file = blast_sequences(query, db)

    genome_lengths = {sequence.id: len(sequence.seq) for sequence in SeqIO.parse(db, 'fasta')}

    coverages = format_blast_results(blast_file, genome_lengths, **blast_filter_prms)

    return coverages

def blast_sequences(query, db, cache=True):
    output = Path('/tmp/cedric/modular_painting_tests/{}_on_{}'.format(query.stem, db.stem))

    if output.is_file() and cache:
        print('No need to do BLAST, found previous results.')
        return output

    # tmpdir = mkdtemp()
    tmpdir = '/tmp/cedric/modular_painting_tests/'

    db_path = Path(tmpdir, db.stem)
    blast_path = Path(tmpdir, f'{query.stem}_on_{db.stem}')

    make_db = NcbimakeblastdbCommandline(dbtype="nucl", input_file=str(db), out=db_path)
    print(make_db)
    make_db()

    blastn_on_db = NcbiblastnCommandline(query=str(query), db=db_path, out=blast_path,
                                         outfmt=5, gapopen=5, gapextend=2, strand='plus')
    print(blastn_on_db)
    blastn_on_db()

    return blast_path

def format_blast_results(blast_file, genome_lengths,
                         min_id=0.8, min_cov=0.5, extend_min_dist=0):

    coverages = defaultdict(lambda: [])

    for record in NCBIXML.parse(open(blast_file)):
        # 1 record is one blast result from the query sequences
        if record.alignments:
            record_id = record.query.split(' ')[0]

            for hit in record.alignments:
                # 1 alignment is one parent (query, with potentially multiple HSPs)
                matches = sum(hsp.identities for hsp in hit.hsps)
                covered = sum(hsp.align_length for hsp in hit.hsps)

                pct_identity = matches / covered
                pct_coverage = covered / record.query_length

                if pct_identity > min_id and pct_coverage > min_cov:
                    hit_id = hit.hit_def.split(' ')[0]

                    hsps = Coverage(*
                        (Arc(hsp.sbjct_start, hsp.sbjct_end, genome_lengths[hit_id], {record_id})
                         for hsp in hit.hsps),
                        target=hit_id
                    )
                    coverages[hit_id].append(hsps)

    coverages = [Coverage.from_coverages(*covs) for (g, covs) in coverages.items()]

    return coverages
