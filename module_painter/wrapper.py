import io
import logging
import shutil
import subprocess as sp
from pathlib import Path

import pandas as pd
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from Bio.SeqIO.FastaIO import as_fasta


ALN_HEADER = ["qacc", "qlen", "qstart", "qend", "sstrand",
              "sacc", "slen", "sstart", "send",
              "nident", "length"]
PAF_HEADER = ALN_HEADER + ["mapq"]
BLAST_HEADER = ALN_HEADER + ["bitscore", "pident"]


logger = logging.getLogger('module-painter')

def check_if_exists(func):
    def wrapper(*args, **kwargs):
        binary = func.__name__
        if not shutil.which(binary):
            raise FileNotFoundError(f"Did you install {binary}?")
        res = func(*args, **kwargs)
        return res
    return wrapper

@check_if_exists
def minimap2(query, ref, **kwargs):
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
    cmd = ["minimap2"]
    for (k, v) in kwargs.items():
        cmd.append(f"-{k}" if len(k) == 1 else f"--{k}")
        if not isinstance(v, bool):
            cmd.append(str(v))
            
    if isinstance(ref, Path):
        cmd += [str(ref), str(query)]
        in_stream = {}
    else:
        cmd += ["-", str(query)]
        in_stream = {"input": as_fasta(ref)}


    logger.info(f"Running: {' '.join(cmd)}")
    aln_str = sp.check_output(
        cmd,
        text=True,
        stderr=sp.DEVNULL,
        **in_stream
    )

    alns = pd.read_csv(io.StringIO(aln_str), header=None, sep="\t",
                      usecols=range(len(PAF_HEADER)), names=PAF_HEADER)
    alns['pident'] = alns.nident / alns.length

    # Fix end to be 0-indexed
    for attr in ["qend", "send"]:
        alns[attr] -= 1

    logger.info(f"Found {alns.shape[0]:,} alignments with minimap2")
    logger.info(f"Aligned: #children: {alns.sacc.nunique()}, #parents: {alns.qacc.nunique()}")

    return alns

@check_if_exists
def blastn(query, ref, outdir=None, **blast_opt):

    # Set output folder
    blast_dir = Path(outdir, 'blast')
    blast_dir.mkdir(exist_ok=True, parents=True)

    # Make blast database
    db_prefix = Path(blast_dir, ref.stem)
    if not Path(db_prefix.with_suffix('.nhr')).is_file():
        make_db = NcbimakeblastdbCommandline(dbtype="nucl", input_file=ref, out=db_prefix)
        logger.info(make_db)
        make_db()

    # Run blast
    blast_file = Path(blast_dir, f'{query.stem}_on_{ref.stem}.tsv')

    fmt = " ".join(BLAST_HEADER)
    blastn_on_db = NcbiblastnCommandline(
        query=query, db=db_prefix, out=blast_file, outfmt=f'6 {fmt}', **blast_opt
    )
    logger.info(blastn_on_db)
    blastn_on_db()

    alns = pd.read_csv(blast_file, sep="\t", header=None, names=BLAST_HEADER)
    alns.pident = alns.pident/100.
    alns.sstrand = alns.sstrand.replace(dict(plus='+', minus='-'))

    # swap sstart and send if strand="-" so that sstart < send
    reverse = alns.sstrand == "-"
    sends = alns.loc[reverse, "send"]
    alns.loc[reverse, "send"] = alns.loc[reverse, "sstart"]
    alns.loc[reverse, "sstart"] = sends

    # Fix start/end to be 0-indexed
    for attr in ["qstart", "qend", "sstart", "send"]:
        alns[attr] -= 1

    logger.info(f"Found {alns.shape[0]:,} alignments with blastn")
    logger.info(f"#Aligned: children: {alns.sacc.nunique()}, #parents: {alns.qacc.nunique()}")

    return alns

def filter_alignments(alns, min_id=0.9):
    # Keep hits with pident > min_id
    logger.info(f"pident >= {min_id:.1%} -> {sum(alns.pident>=min_id):,} alignments remaining")
    alns = alns[alns.pident >= min_id]

    # Stop if no alignments left
    if alns.empty:
        logger.error(f"All alignments were discarded. Try decreasing --min-id.")
        return alns

    # For each query, keep hits with sstrand of the best hit
    sstrands = alns.groupby(["sacc", "qacc"]).apply(lambda df: df.set_index("sstrand").nident.idxmax()).to_dict()
    alns = alns.groupby(["sacc", "qacc"], group_keys=False).apply(lambda df: df[df.sstrand==sstrands[df.name]])
    logger.debug(f"Strand filtering -> {alns.shape[0]:,} alignments remaining")

    # Discard children if only one parent is found
    alns = alns.groupby("sacc").filter(lambda x: x.qacc.nunique() > 1)
    logger.info(f"{alns.sacc.nunique():,} children remaining (>2 parents found)")

    return alns
