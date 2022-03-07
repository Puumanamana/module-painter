from pathlib import Path

import pandas as pd
import numpy as np
import igraph
import mappy

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def subset_fasta(filename, subset_df, min_len=0, outprefix="/tmp"):
    sequences = []

    # Make output directory if it doesn't exist
    output = Path(outprefix, Path(filename).name)
    output.parent.mkdir(exist_ok=True, parents=True)

    refs = set(subset_df.index.get_level_values("ref"))

    with open(filename, 'r') as reader, open(output, 'w') as writer:
        for (title, seq) in SimpleFastaParser(reader):
            seq_id = title.split()[0]
            if seq_id not in refs:
                continue

            for (start, end) in subset_df.loc[seq_id].index:
                metadata = str(subset_df.loc[(seq_id, start, end)])

                seq_len = end-start
                extension = 0
                if seq_len < min_len: # artificially extend segment
                    extension = (min_len-seq_len + 1) // 2

                subseq = wrapping_substr(seq, start-extension, end+extension+1)
                writer.write(f">{seq_id}|{start}|{end}|{metadata}\n{subseq}\n")

    return output

def build_homology_graph(fasta, min_id=0.9, verbose=0, threads=1):
    db = mappy.Aligner(str(fasta), n_threads=threads)  # build index
    if not db: raise Exception("ERROR: failed to load/build index")

    (vertices, edges) = (set(), set())

    weights = {}
    
    print(f"Building minimap2 graph for {fasta}")
    for name, seq, _ in mappy.fastx_read(str(fasta)):
        vertices.add(name)

        for hit in db.map(seq): # traverse alignments
            if name == hit.ctg:
                continue

            pident = hit.mlen/min([hit.ctg_len, len(seq)])

            if verbose > 0:
                print("{} vs {} (id={:.1%}, mapq={})".format(
                    name, hit.ctg, pident, hit.mapq,
                    hit.r_st, hit.r_en, hit.q_st, hit.q_en))

            if pident > min_id:
                edge = tuple(sorted((name, hit.ctg)))
                edges.add(edge)
                vertices |= {name, hit.ctg}

    graph = igraph.Graph()
    graph.add_vertices(sorted(vertices))
    graph.add_edges(sorted(edges))

    return graph

def wrapping_substr(s, i, j):
    if 0 <= i < j <= len(s):
        return s[i:j]
    if i < 0 < j <= len(s):
        return s[i:] + s[:j]
    if 0 < i < len(s) < j:
        return s[i:] + s[:(j % len(s))]
    # if 0 < j < len(s) < i:
    #     return s[(i % len(s)):j]
        
    raise ValueError(f"String subset not defined for start={i} and end={j} (L={len(s)})")
