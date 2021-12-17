from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
import mappy

import matplotlib.pyplot as plt

from modular_painter.display import get_palette
from modular_painter.arc import Arc
from modular_painter.coverage import Coverage

## TO DO: only create database once for each sample

def minimap(query, db, threads=10):
    db = mappy.Aligner(
        str(db),
        best_n=20,
        # preset="splice",
        scoring=[2,4,8,6],
        n_threads=threads
    )  # build index
    if not db: raise Exception("ERROR: failed to load/build index")

    arcs = defaultdict(lambda: [])
    for qname, qseq, _ in mappy.fastx_read(str(query)):
        for hit in db.map(qseq):
            if hit.mapq > 0:
                arcs[hit.ctg].append(
                    Arc(hit.r_st, hit.r_en, hit.ctg_len, {qname})
                )

    coverages = [Coverage(*arc_list, target=g) for (g, arc_list) in arcs.items()]

    for cov in coverages:
        if cov.target in "RSWXd":
            print(cov)
    import ipdb;ipdb.set_trace()
    return coverages
