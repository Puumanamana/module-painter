from pathlib import Path

import pandas as pd
import numpy as np
from igraph import Graph
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

    saccs = set(subset_df.index.get_level_values("sacc"))

    with open(filename, 'r') as reader, open(output, 'w') as writer:
        for (title, seq) in SimpleFastaParser(reader):
            seq_id = title.split()[0]
            if seq_id not in saccs:
                continue

            for (sstart, send) in subset_df.loc[seq_id].index:
                metadata = str(subset_df.loc[(seq_id, sstart, send)])

                seq_len = send-sstart
                extension = 0
                if seq_len < min_len: # artificially extend segment
                    extension = (min_len-seq_len + 1) // 2

                subseq = wrapping_substr(seq, sstart-extension, send+extension+1)
                writer.write(f">{seq_id}^{sstart}^{send}^{metadata}\n{subseq}\n")

    return output

def build_homology_graph(fasta, min_id=0.9, min_size_ratio=0, verbose=0, threads=1):
    db = mappy.Aligner(str(fasta), n_threads=threads)  # build index
    if not db: raise Exception("ERROR: failed to load/build index")

    vertices = set()
    edges = {}

    weights = {}

    if verbose:
        print(f"Building minimap2 graph for {fasta}")
    for name, seq, _ in mappy.fastx_read(str(fasta)):
        vertices.add(name)

        for hit in db.map(seq): # traverse alignments
            if name == hit.ctg:
                continue

            lengths = [hit.ctg_len, len(seq)]
            (min_len, max_len) = min(lengths), max(lengths)
            pident = hit.mlen/hit.blen

            if verbose > 0:
                print("{} vs {} (id={:.1%}, mapq={})".format(
                    name, hit.ctg, pident, hit.mapq,
                    hit.r_st, hit.r_en, hit.q_st, hit.q_en))

            if pident > min_id and min_len/max_len > min_size_ratio:
                edge = tuple(sorted((name, hit.ctg)))
                edges[edge] = dict(pident=pident, matches=hit.mlen, length=hit.blen, qlen=len(seq), rlen=hit.ctg_len, qseq=seq)
                vertices |= {name, hit.ctg}

    graph = Graph()
    graph.add_vertices(sorted(vertices))
    graph.add_edges(edges.keys())
    if edges:
        for attr in edges[edge]:
            graph.es[attr] = [v[attr] for v in edges.values()]
    graph.es["name"] = ["<->".join(graph.vs["name"][vi] for vi in e.tuple) for e in graph.es]

    return graph

def wrapping_substr(s, i, j):
    if 0 <= i < j <= len(s):
        return s[i:j]
    if i < 0 < j <= len(s):
        return s[i:] + s[:j]
    if 0 < i < len(s) < j:
        return s[i:] + s[:(j % len(s))]
    if 0 < j < len(s) < i:
        return s[(i % len(s)):j]
    raise ValueError(f"String subset not defined for sstart={i} and send={j} (L={len(s)})")

def concatenate_fasta(*fa, outdir="./", resume=False, rename=False, prefix="seq", min_length=1000):
    if not rename and len(fa) == 1:
        return fa[0]

    outprefix = "-".join(sorted(Path(f).stem for f in fa))
    if len(outprefix) > 20:
        outprefix = outprefix[:20] + "_trunc"
    output = Path(outdir, outprefix).with_suffix(".fasta")

    if resume and output.is_file():
        return output

    count = 0
    with open(output, "w") as writer:
        for f in fa:
            with open(f) as reader:
                for (title, seq) in SimpleFastaParser(reader):
                    if len(seq) > min_length:
                        if rename:
                            title = f"{prefix}_{count}"
                        writer.write(f">{title}\n{seq}\n")
                        count += 1
    return output
