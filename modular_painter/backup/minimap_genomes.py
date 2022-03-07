from collections import defaultdict
from pathlib import Path
import io
import subprocess as sp

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import Align
from Bio.SeqIO.FastaIO import as_fasta

import mappy

from modular_painter.coverage import Coverage


### minimap2 parameters
# default: -s 80 -z400,200 -N  5 -U  10,1e6
# asm20  : -s200 -z200,200 -N 50 -U  50,500    -r100k -k19 --rmq -g10k -A1 -B4 -O6,26 -E2,1
# synteny: -s100 -z200,200 -N 50 -U 100,1e6 -c
# here: synteny and #-r {min_module_size}

PAF_HEADER = ["q_acc", "q_len", "q_start", "q_end", "orient",
              "t_acc", "t_len", "t_start", "t_end",
              "matches", "length", "mapq"]

def run_minimap(parents, child, seq_data,
                k=15, r=50, s=100, z=100, N=50, U=100, c=True, threads=10):
    """
    Minimap2 parameters used here:
    -k INT	
      Minimizer k-mer length [15]
    -s INT	
      Minimal peak DP alignment score to output [40]. 
      The peak score is computed from the final CIGAR. It is the score of the max scoring segment 
      in the alignment and may be different from the total alignment score.
    -z INT1[,INT2]
      Truncate an alignment if the running alignment score drops too quickly along the diagonal 
      of the DP matrix (diagonal X-drop, or Z-drop) [400,200]. 
      If the drop of score is above INT2, minimap2 will reverse complement the query in the related
      region and align again to test small inversions. Minimap2 truncates alignment if there is an
      inversion or the drop of score is greater than INT1. 
      Decrease INT2 to find small inversions at the cost of performance and false positives. 
      Increase INT1 to improves the contiguity of alignment at the cost of poor alignment in the 
      middle.
    -r NUM1[,NUM2]
      Bandwidth for chaining and base alignment [500,20k]. 
      NUM1 is used for initial chaining and alignment extension; 
      NUM2 for RMQ-based re-chaining and closing gaps in alignments.
    -G NUM	
      Maximum gap on the reference (effective with -xsplice/--splice). 
      This option also changes the chaining and alignment band width to NUM. 
      Increasing this option slows down spliced alignment. [200k]
    -n INT 
      Discard chains consisting of <INT number of minimizers [3]
    -N INT	
      Output at most INT secondary alignments [5]. This option has no effect when -X is applied.
    -U INT1[,INT2]
      Lower and upper bounds of k-mer occurrences [10,1000000]. The final k-mer occurrence 
      threshold is max{INT1, min{INT2, -f}}. This option prevents excessively small or large -f 
      estimated from the input reference. Available since r1034 and deprecating --min-occ-floor 
      in earlier versions of minimap2.
    -t INT	
      Number of threads [3]. Minimap2 uses at most three threads when indexing target sequences, 
      and uses up to INT+1 threads when mapping (the extra thread is for I/O, which is frequently 
      idle and takes little CPU time).
    """
    cmd = ["minimap2", "-", parents,
           "--no-long-join",
           "-k", k, 
           "-s", s, "-z", f"{z},{z}",
           "-N", N,
           "-M", 0.01,
           # "-r", r,
           "-n", 2,
           "-U", U,
           "-t", threads] + ["-c"]*c

    aln_str = sp.check_output(
        map(str, cmd),
        input=as_fasta(child),
        text=True,
        stderr=sp.DEVNULL
    )
    aln = pd.read_csv(io.StringIO(aln_str), header=None, sep="\t",
                      usecols=range(len(PAF_HEADER)), names=PAF_HEADER)
    aln['pident'] = aln.matches / aln.length

    return aln


def check_gaps(aln, seq_data, child):
    aln = (aln.sort_values(by=["q_acc", "t_start", "t_end"])
           .groupby("q_acc")
           .apply(check_gaps, qseqs=seq_data, tseq=child.seq)
           .reset_index(drop=True))

    aln = Coverage.from_pandas(aln, size=len(child))

    return aln

def run_mappy(query, ref, threads=10, N=50, preset="asm20", **kwargs):
    db = mappy.Aligner(str(ref), best_n=N, preset=preset, n_threads=threads)
    if not db: raise Exception("ERROR: failed to load/build index")

    hits = []
    for qname, qseq, _ in mappy.fastx_read(str(query)):
        for hit in db.map(qseq):
            hits.append([
                qname, len(qseq), hit.q_st, hit.q_en, hit.strand,
                hit.ctg, hit.ctg_len, hit.r_st, hit.r_en,
                hit.mlen, hit.blen, hit.mapq
            ])

    hits = pd.DataFrame(hits, columns=PAF_HEADER)

    return hits

def check_gaps(hits, qseqs, tseq, min_gap=10, max_gap=1e4, min_size_ratio=0.9, min_seq_id=0.75):
    # to align gap sequences later
    aligner = Align.PairwiseAligner(
        mode="global",
    )

    # get gap boundaries
    qseq = qseqs[hits.name]

    hits["next_q"] = get_next_hit(hits.q_start, len(qseq))
    hits["next_t"] = get_next_hit(hits.t_start, len(tseq))

    gap_cols = ["q_end", "next_q", "t_end", "next_t"]

    scores = np.zeros(len(hits))
    
    for i, (gap_start_q, gap_end_q, gap_start_t, gap_end_t) in enumerate(hits[gap_cols].values):
        gaps = [gap_end_q-gap_start_q, gap_end_t-gap_start_t]
        gap_ratio = min(gaps) / max(gaps)

        # if the gaps have the same size of reference and query, we check the identity
        if gap_ratio > min_size_ratio and max(gaps) < max_gap:
            if all(g < 0 for g in gaps):
                scores[i] = 1.1
            else:
                gap_q = qseq[gap_start_q:gap_end_q]
                gap_t = tseq[gap_start_t:gap_end_t]

                if gap_end_q > len(qseq):
                    gap_q += qseq[:(gap_end_q % len(qseq))]
                if gap_end_t > len(tseq):
                    gap_t += tseq[:(gap_end_t % len(tseq))]

                aln = aligner.align(gap_q, gap_t)

                scores[i] = aln.score / max(len(gap_q), len(gap_t))

    hits["gap_score"] = scores
    hits["merge_with_next"] = scores > min_seq_id

    hits = merge_hits(hits)

    return hits

def get_next_hit(starts, size):
    """
    Get start position of the next hit
    The last start value wraps around (=starts[0]+size)
    """
    next_start = np.roll(starts.values, -1)
    next_start[-1] += size
    return next_start

def merge_hits(hits):
    indices = hits.index
    to_remove = []
    
    for i, merge in enumerate(hits.merge_with_next):
        if merge:
            current_idx = indices[i]
            next_idx = indices[(i+1) % len(indices)]
            to_remove.append(next_idx)

            hits.loc[current_idx, "t_end"] = hits.loc[next_idx, "t_end"]
            hits.loc[current_idx, "q_end"] = hits.loc[next_idx, "q_end"]

            if (i+1) == len(indices):
                hits.loc[current_idx, "t_end"] += hits.loc[next_idx, "t_len"]
                hits.loc[current_idx, "q_end"] += hits.loc[next_idx, "q_len"]                

    hits = hits.drop(to_remove)

    return hits
                  
