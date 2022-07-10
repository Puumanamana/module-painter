import logging
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from igraph import Graph

from module_painter.util import concatenate_fasta
from module_painter.coverage import Coverage
from module_painter.wrapper import blastn, minimap2, filter_alignments
from module_painter.breakpoints import map_missing_parents, set_breakpoint_ids
from module_painter.parent_selection import summarize_breakpoints, select_by_recombinations, select_by_breakpoints
from module_painter.clustering import cluster_phages
from module_painter.display import display_genomes


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
    minimap2={"z": "50,50", "N": 100, "U": 100,
              "no-long-join": True, "c": True, "P": True},
    blastn={"gapopen": 5, "gapextend": 2}
)

def paint(parents=None, children=None, outdir=None, resume=False, rename=False, aligner="minimap2",
          min_length=0, min_module_size=0, min_id=0.9, arc_eq_diffs=0, min_nw_id=0.8, skip_nw=False,
          clustering_feature="breakpoint", clustering_gamma=0.2, plot_fmt="html",
          threads=20, **kwargs):

    logger = logging.getLogger("module-painter")

    populations = [
        concatenate_fasta(*parents, outdir=outdir, min_length=min_length,
                          resume=resume, rename=rename, prefix="p"),
        concatenate_fasta(*children, outdir=outdir, min_length=min_length,
                          resume=resume, rename=rename, prefix="c")
    ]
    seq_data = {seq.id: seq.seq for fasta in populations for seq in SeqIO.parse(fasta, "fasta")}

    # Raw alignment
    aln_file = Path(outdir, f"{aligner}_{populations[0].stem}_on_{populations[1].stem}.csv")
    if resume and aln_file.is_file():
        alns = pd.read_csv(aln_file, index_col=0)
    else:
        if "blast" in aligner:
            alns = blastn(*populations, outdir, **ALN_PARAMS["blastn"])
        else:
            alns = minimap2(*populations, r=min_module_size, t=threads,
                            **ALN_PARAMS["minimap2"])
        alns.to_csv(aln_file)

    alns = filter_alignments(alns, min_id=min_id)
        
    logger.info("Applying coverage rules")
    coverages = []
    for child in alns.sacc.unique():
        cov_file = Path(f"{outdir}/simplified_coverage/{child}.csv")
        if resume and cov_file.is_file():
            logger.info(f"{child} already processed. Skipping")
            coverage = Coverage.from_csv(cov_file)
        else:
            coverage = Coverage.from_pandas(alns[alns.sacc==child])
            coverage.apply_rules(arc_eq_diffs, min_module_size, seq_data, min_nw_id, skip_nw=skip_nw)
            coverage.write_csv(f"{outdir}/simplified_coverage")

        if len(coverage) < 2:
            logger.debug(f"Discarding {coverage.sacc} (only one parent left)")
        else:
            coverages.append(coverage)

    if not coverages: return

    logger.info(f"Mapping missing parents (n_children={len(coverages)})")
    map_missing_parents(populations[1], coverages, outdir=outdir, threads=threads)

    graph_dir = Path(outdir, "overlap_graphs")
    graph_dir.mkdir(exist_ok=True)
    graph_paths = {cov.sacc: Path(graph_dir, f"{cov.sacc}.pickle") for cov in coverages}

    if resume and all(f.is_file() for f in graph_paths.values()):
        overlap_graphs = [Graph.Read_Pickle(f) for f in graph_paths.values()]
    else:
        overlap_graphs = [cov.get_overlap_graph(min_overlap=50) for cov in coverages]
        logger.info("Mapping breakpoints")
        set_breakpoint_ids(overlap_graphs, populations[1], outdir=outdir, threads=threads)

        logger.info("Parent selection: by recombinations")
        select_by_recombinations(overlap_graphs)
        logger.info("Parent selection: by breakpoint")
        select_by_breakpoints(overlap_graphs)

        summarize_breakpoints(overlap_graphs)

        for graph in overlap_graphs:
            graph.write_pickle(graph_paths[graph["sacc"]])
    
    logger.info(f"Clustering (feature: {clustering_feature})")
    clusters = cluster_phages(overlap_graphs, gamma=clustering_gamma, feature=clustering_feature, outdir=outdir)
    clusters = [c for c in clusters if len(c) > 1]

    if not clusters:
        logger.warning("No species cluster identified.")
        return

    for i, c in enumerate(sorted(clusters, key=lambda x: -len(x))):
        logger.info(f"Cluster #{i}: {','.join(c)}")

    logger.info(f"Interactive plot in {outdir}")
    display_genomes(overlap_graphs, clusters=clusters, norm=True, outdir=outdir, fmt=plot_fmt)

