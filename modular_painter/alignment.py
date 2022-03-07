from pathlib import Path
import io
import subprocess as sp

import numpy as np
import pandas as pd

from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from Bio import SeqIO
from Bio import Align
from Bio.SeqIO.FastaIO import as_fasta

from modular_painter.util import wrapping_substr


ALN_HEADER = ["qacc", "qlen", "qstart", "qend", "sstrand",
              "sacc", "slen", "sstart", "send",
              "nident", "length"]
PAF_HEADER = ALN_HEADER + ["mapq"]
BLAST_HEADER = ALN_HEADER + ["bitscore", "pident"]


def align(parents, children, algo='blastn', output=None, min_pident=.9, **aln_params):
    print(f"Initial alignment with {algo}")
    if "blast" in algo.lower():
        aln = run_blastn(parents, children, output, **aln_params)
    else:
        aln = run_minimap(parents, children, **aln_params)

    # Retrieve raw sequences
    seq_data = {seq.id: seq.seq for fname in [parents, children] for seq in SeqIO.parse(fname, "fasta")}

    print("Refining alignment gaps with Needleman-Wunsch")
    # Align gap and try to extend if possible
    aln = (aln[aln.pident >=min_pident]
           .sort_values(by=["sacc", "qacc", "sstart", "send"])
           .groupby(["sacc", "qacc"], sort=False)
           .apply(align_gaps, seq_data=seq_data)
           .reset_index(drop=True))
    aln = aln[~aln.merge_with_next]

    return aln

def run_blastn(parents, children, output=None, **blast_filter_prms):

    # Set output folder
    blast_dir = Path(output, 'blast')
    blast_dir.mkdir(exist_ok=True, parents=True)

    # Make blast database
    db_prefix = Path(blast_dir, children.stem)
    if not Path(db_prefix.with_suffix('.nhr')).is_file():
        make_db = NcbimakeblastdbCommandline(dbtype="nucl", input_file=children, out=db_prefix)
        print(make_db)
        make_db()

    # Run blast
    blast_file = Path(blast_dir, f'{parents.stem}_on_{children.stem}.tsv')
    if not blast_file.is_file():
        fmt = " ".join(BLAST_HEADER)
        blastn_on_db = NcbiblastnCommandline(
            query=parents, db=db_prefix, out=blast_file,
            gapopen=5, gapextend=2, strand='plus',
            outfmt=f'6 {fmt}')
        print(blastn_on_db)
        blastn_on_db()

    aln = pd.read_csv(blast_file, sep="\t", header=None, names=BLAST_HEADER)
    aln.pident = aln.pident/100.
    aln.sstrand = aln.sstrand.replace(dict(plus='+', minus='-'))

    # Fix start/end to be 0-indexed
    for attr in ["qstart", "qend", "sstart", "send"]:
        aln[attr] -= 1

    return aln

def run_minimap(parents, children,
                k=15, r=150, s=100, z=100, N=50, U=100, c=True, threads=10):
    """
    Minimap2 parameters used here:
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
    -U INT1[,INT2]
      Lower and upper bounds of k-mer occurrences [10,1000000]. The final k-mer occurrence 
      threshold is max{INT1, min{INT2, -f}}. This option prevents excessively small or large -f 
      estimated from the input reference. Available since r1034 and deprecating --min-occ-floor 
      in earlier versions of minimap2.
    """
    cmd = ["minimap2", "-", parents,
           "--no-long-join",
           "-k", k, 
           "-s", s, "-z", f"{z},{z}",
           "-N", N,
           "-M", 0.01,
           "-r", r,
           "-n", 2,
           "-U", U,
           "-t", threads] + ["-c"]*c

    alns = []
    
    for child in SeqIO.parse(children, "fasta"):
        aln_str = sp.check_output(
            map(str, cmd),
            input=as_fasta(child),
            text=True,
            stderr=sp.DEVNULL
        )

        aln = pd.read_csv(io.StringIO(aln_str), header=None, sep="\t",
                          usecols=range(len(PAF_HEADER)), names=PAF_HEADER)
        aln['pident'] = aln.nident / aln.length

        alns.append(aln[aln.mapq>30])

    alns = pd.concat(alns)
    
    # Fix end to be 0-indexed
    for attr in ["qend", "send"]:
        alns[attr] -= 1

    return alns

def align_gaps(hits, seq_data, min_gap=10, max_gap=1e3, min_size_ratio=0.8, min_seq_id=0.75):
    # to align gap sequences later
    aligner = Align.PairwiseAligner(
        mode="global",
        open_gap_score=-0.1,
        # open_gap_score=-1,
        # extend_gap_score=-0.1,
        # target_open_gap_score=-1,
        # target_extend_gap_score=-0.2,
        # query_open_gap_score=-1,
        # query_extend_gap_score=0.5,
    )

    # get gap boundaries
    (sacc, qacc) = hits.name

    sseq = seq_data[sacc]
    qseq = seq_data[qacc]

    hits["next_q"] = get_next_hit(hits.qstart, len(qseq))
    hits["next_s"] = get_next_hit(hits.sstart, len(sseq))

    gap_cols = ["qend", "next_q", "send", "next_s"]

    scores = np.zeros(len(hits))

    # Align gaps on reference and query and fill gap
    for i, (gap_start_q, gap_end_q, gap_start_s, gap_end_s) in enumerate(hits[gap_cols].values):
        gaps = [gap_end_q-gap_start_q, gap_end_s-gap_start_s]
        gap_ratio = gaps[0] / gaps[1] # what matters is that the query is longer than the reference

        if gaps[0] >= 0 and -min_gap <= gaps[1] <= min_gap: # insertion on the parent, but not the child
            scores[i] = 1.1

        # if the gaps have the same size of reference and query, we check the identity
        elif gap_ratio > min_size_ratio and max(gaps) < max_gap:
            if all(g < 0 for g in gaps):
                continue
                # scores[i] = 1.1
            else:
                gap_q = wrapping_substr(qseq, gap_start_q, gap_end_q)
                gap_s = wrapping_substr(sseq, gap_start_s, gap_end_s)

                aln = aligner.align(gap_q, gap_s)

                scores[i] = aln.score / len(gap_s)

    hits["gap_score"] = scores
    hits["merge_with_next"] = scores > min_seq_id

    # if sacc == "a" and qacc in {"B", "G"}:
    #     print(sacc, qacc)
    #     print(hits[["qstart", "qend", "qlen", "sstart", "send", "slen", "pident", "gap_score", "mapq"]])
    #     (gap_start_q, gap_end_q, gap_start_s, gap_end_s) = hits[gap_cols].values[1]
    #     gap_ratio = gaps[0] / gaps[1]
    #     print(gap_ratio)
    #     gap_q = wrapping_substr(qseq, gap_start_q, gap_end_q)
    #     gap_s = wrapping_substr(sseq, gap_start_s, gap_end_s)

    #     aln = aligner.align(gap_q, gap_s)
    #     # open("test.afa", "w").write(str(aln[0]))
    #     print(aln.score/len(gap_s))
    #     import ipdb;ipdb.set_trace()
        
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
    
    for i, merge in enumerate(hits.merge_with_next):
        if merge:
            current_idx = indices[i]
            next_idx = indices[(i+1) % len(indices)]

            for which in ["q", "s"]:
                field = f"{which}start"
                hits.loc[next_idx, field] = hits.loc[current_idx, field]

                if (i+1) == len(indices): # wrapping condition
                    hits.loc[next_idx, f"{which}end"] += hits.loc[current_idx, f"{which}len"]

    # special case when everything is merged
    if hits.merge_with_next.all(): # we keep the first one
        hits.loc[indices[0], "merge_with_next"] = False

    return hits
                  
