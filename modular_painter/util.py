from pathlib import Path
from tempfile import mkdtemp

import pandas as pd
import numpy as np
import igraph
import mappy

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def subset_fasta(filename, subset_df, output=None, min_len=0):

    sequences = []

    with open(filename, 'r') as handle:
        for (title, seq) in SimpleFastaParser(handle):
            v_id = title.split()[0]
            if v_id not in subset_df.index:
                continue

            indices = map(str, subset_df.loc[v_id, 'index'])
            starts = subset_df.loc[v_id, 'start']
            ends = subset_df.loc[v_id, 'end']

            for (idx, start, end) in zip(indices, starts, ends):
                if end-start < min_len:
                    continue

                bp_id = "{}_{}-{}".format(v_id, start, end)
                bp_sequence = SeqRecord(Seq(seq[start-1:end]), id=idx, description=bp_id)
                sequences.append(bp_sequence)

    if not output:
        # tmpdir = mkdtemp()
        tmpdir = "/tmp/cedric/modular_painting_tests"
        output = "{}/subset_{}".format(tmpdir, Path(filename).name)

    SeqIO.write(sequences, output, 'fasta')

    return output


def get_minimap_graph(fasta, min_id=0.8, min_cov=0.5, verbose=0, threads=1):
    db = mappy.Aligner(fasta, n_threads=threads, k=16)  # build index
    if not db: raise Exception("ERROR: failed to load/build index")

    (vertices, edges) = (set(), set())

    for name, seq, _ in mappy.fastx_read(fasta):
        vertices.add(name)
        for hit in db.map(seq): # traverse alignments
            if name == hit.ctg:
                continue

            if verbose > 0:
                print("{} vs {} ({}/{} matches): ref {}-{}, query {}-{}".format(
                    name, hit.ctg, hit.mlen, hit.blen, hit.r_st, hit.r_en, hit.q_st, hit.q_en))

            if (hit.mlen/hit.blen > min_id and
                any(hit.blen/l > min_cov for l in [hit.ctg_len, len(seq)])):
                edge = sorted((name, hit.ctg))
                edges.add(tuple(edge))
                vertices |= {name, hit.ctg}

    graph = igraph.Graph()
    graph.add_vertices(list(vertices))
    graph.add_edges(list(edges))

    edges_array = np.array(graph.vs['name']).astype(int)

    components = dict((edges_array[edge], i) for (i, component) in enumerate(graph.components()) for edge in component)

    return pd.Series(components)
