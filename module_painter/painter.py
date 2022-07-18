import logging
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from igraph import Graph

from module_painter.io import find_outputs
from module_painter.util import concatenate_fasta
from module_painter.coverage import Coverage
from module_painter.wrapper import blastn, minimap2, filter_alignments
from module_painter.breakpoints import map_missing_parents, set_breakpoint_ids
from module_painter.parent_selection import select_by_recombinations, select_by_breakpoints, get_breakpoints
from module_painter.clustering import get_links, summarize_recombinations

TRUTH = dict(
    R="MAJAM",
    S="(JM)FJA(JM)",
    W="JFJAJ",
    X="(JM)F(AJ)MA(JM)",
    d="AMJA",
    Z="B-B",
    a="B-B",
    T="LB-B--L(BCEFHIK)-B--(BDG)-(BDGHK)-(ACDEFGHIJKL)B(AJM)L",
    U="LB-B-BL(BCEFHIK)-B--(BDG)-(BDGHK)-(ACDEFGHIJKL)B(AJM)L",
    V="LB-B--BL(BCEFHIK)-B--B(AJM)L",
    Y="(DGK)(ACFHJM)LAL(AJM)B(CFH)(BDGK)A(DGK)(AJM)LGB(DGK)",
    b="A(EI)-BLAL(AJ)(IL)A",
    c="L(EI)L(BEI)-(ABCFHJM)LBL(BDGK)(EI)(CFH)BLBIHBL"
)

ALN_PARAMS = dict(
    minimap2={"z": 50, "U": 100,
              "no-long-join": True, "c": True, "P": True},
    blastn={"gapopen": 5, "gapextend": 2}
)

def paint(parents=None, children=None, outdir=None, resume=False, rename=False, aligner="minimap2",
          min_length=0, min_module_size=0, min_id=0.9, arc_eq_diffs=0, min_nw_id=0.8, skip_nw=False,
          clustering_feature="breakpoint", clustering_gamma=0.2, plot_fmt="html",
          threads=20, **kwargs):

    logger = logging.getLogger("module-painter")

    (folders, outputs) = find_outputs(outdir, aligner)

    populations = [
        concatenate_fasta(*parents, outdir=outdir, min_length=min_length,
                          resume=resume, rename=rename, prefix="p"),
        concatenate_fasta(*children, outdir=outdir, min_length=min_length,
                          resume=resume, rename=rename, prefix="c")
    ]
    seq_data = {seq.id: seq.seq for fasta in populations for seq in SeqIO.parse(fasta, "fasta")}

    # Raw alignment
    if not outputs["raw_aln"].is_file() or not resume:        
        if "blast" in aligner:
            alns = blastn(*populations, outdir, **ALN_PARAMS["blastn"])
        else:
            alns = minimap2(*populations, r=min_module_size, t=threads,
                            **ALN_PARAMS["minimap2"])
        alns = filter_alignments(alns, min_id=min_id)
        alns.to_csv(outputs["raw_aln"])
    else:
        logger.info(f"Loading {aligner} alignment (resume)")
        alns = pd.read_csv(outputs["raw_aln"], index_col=0)

    coverages = []
    for child in alns.sacc.unique():
        cov_file = Path(folders["raw_cov"], f"{child}.csv")
        # If cov_file doesn't exist but next step started, then len(cov) < 2
        next_step_started = any(outputs["@mapped"])
        if resume and (cov_file.is_file() or next_step_started):
            logger.info(f"Loading {child} (resume)")
            if cov_file.is_file(): 
                coverages.append(Coverage.from_csv(cov_file))
            continue

        coverage = Coverage.from_pandas(alns[alns.sacc==child])
        coverage.apply_rules(arc_eq_diffs, min_module_size, seq_data, min_nw_id, skip_nw=skip_nw)
        coverage.to_csv(cov_file)

        if len(coverage) < 2: logger.debug(f"Discarding {child} (only one parent left)")
        else: coverages.append(coverage)

    if not coverages:
        logger.warning(f"All coverages were filtered out. Exiting")
        return dict()

    if resume and all(cov.sacc in outputs["@mapped"] for cov in coverages):
        logger.info(f"Loading coverages with missing parents mapped (resume)")
        coverages = [Coverage.from_csv(out) for out in outputs["@mapped"].values()]
    else:
        logger.info(f"Mapping missing parents (n_children={len(coverages)})")
        map_missing_parents(populations[1], coverages, outdir=folders["@mapped"], threads=threads)

    if resume and all(cov.sacc in outputs["graphs"] for cov in coverages):
        logger.info(f"Loading overlap graphs (resume)")
        overlap_graphs = [Graph.Read_Pickle(f) for f in outputs["graphs"].values()]
    else:
        overlap_graphs = [cov.get_overlap_graph(min_overlap=50) for cov in coverages]
        logger.info("Mapping breakpoints")
        set_breakpoint_ids(overlap_graphs, populations[1],
                           outdir=folders["graphs"], threads=threads)
        logger.info("Parent selection: by recombinations")
        select_by_recombinations(overlap_graphs)
        logger.info("Parent selection: by breakpoint")
        select_by_breakpoints(overlap_graphs)

        for graph in overlap_graphs:
            graph.write_pickle(Path(folders["graphs"], f"{graph['sacc']}.pickle"))

    breakpoints = get_breakpoints(overlap_graphs)
    rc = summarize_recombinations(breakpoints, outdir)
    links = get_links(breakpoints, feature=clustering_feature)

    return dict(
        graphs=overlap_graphs,
        breakpoints=breakpoints,
        rc=rc,
        links=links
    )
