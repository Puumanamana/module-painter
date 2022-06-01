from pathlib import Path
import pandas as pd
from Bio import SeqIO
from igraph import Graph

from module_painter.io import parse_args, setup_logger
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
    minimap2={"z": "50,50", "N": 50, "U": 100,
              "no-long-join": True, "c": True, "P": True},
    blastn={"gapopen": 5, "gapextend": 2}
)

def main():
    
    args = parse_args()
    logger = setup_logger('module-painter', Path(args.outdir, 'log.txt'))
    logger.info("Starting module-painter")
    logger.info(f"Parent files: {[str(x) for x in args.parents]}")
    logger.info(f"Children files: {[str(x) for x in args.children]}")

    populations = [
        concatenate_fasta(*args.parents, outdir=args.outdir, min_length=args.min_length,
                          resume=args.resume, rename=args.rename, prefix="p"),
        concatenate_fasta(*args.children, outdir=args.outdir, min_length=args.min_length,
                          resume=args.resume, rename=args.rename, prefix="c")
    ]

    # Raw alignment
    aln_file = Path(args.outdir, f"{args.aligner}_{populations[0].stem}_on_{populations[1].stem}.csv")
    if args.resume and aln_file.is_file():
        alns = pd.read_csv(aln_file, index_col=0)
    else:
        if "blast" in args.aligner:
            alns = blastn(*populations, args.outdir, **ALN_PARAMS["blastn"])
        else:
            alns = minimap2(*populations, r=args.min_module_size, t=args.threads,
                            **ALN_PARAMS["minimap2"])
        alns.to_csv(aln_file)

    alns = filter_alignments(alns, min_id=args.min_id)

    if alns.empty:
        return
        
    seq_data = {seq.id: seq.seq for fasta in populations for seq in SeqIO.parse(fasta, "fasta")}

    logger.info("Applying coverage rules for all children")
    coverages = []
    for child in alns.sacc.unique():
        aln = alns[alns.sacc==child]
        if aln.qacc.nunique() == 1:
            logger.warning(f"Too few parents cover {child}. Skipping")
            continue
        logger.debug(f"Processing {child}")
        if args.use_ground_truth:
            aln = alns[alns.qacc.isin(set(TRUTH[child]))]

        coverage = Coverage.from_pandas(aln, size=aln.slen.iloc[0])
        # sticky boundaries
        coverage.sync_boundaries("start", args.arc_eq_diffs)
        coverage.sync_boundaries("end", args.arc_eq_diffs)
        # Fuse intervals from the same parents if they are close
        coverage.fuse_close_modules(seq_data,
                                    max_dist=args.min_module_size,
                                    min_size_ratio=0.8,
                                    min_nw_id=args.min_nw_id)
        # simplify coverage by grouping equal intervals
        coverage.merge_equal_intervals()
        # Remove embedded intervals
        coverage.simplify_embedded()
        # Fill all gaps
        coverage.fill_gaps(args.min_module_size)

        logger.debug(coverage.show())

        if len(coverage) < 2:
            logger.debug(f"Discarding {coverage.ref} (only one parent left)")
            continue
        # Lee and Lee
        coverage.get_minimal_coverage()
        coverages.append(coverage)

    logger.info(f"{len(coverages)} children remaining.")
        
    if not coverages:
        return

    tmp = {c.ref: c for c in coverages}

    logger.info("Mapping missing parents")
    map_missing_parents(populations[1], coverages, outdir=args.outdir, threads=1)

    graph_dir = Path(args.outdir, "overlap_graphs")
    graph_dir.mkdir(exist_ok=True)
    graph_paths = {cov.ref: Path(graph_dir, f"{cov.ref}.pickle") for cov in coverages}

    if args.resume and all(f.is_file() for f in graph_paths.values()):
        overlap_graphs = [Graph.Read_Pickle(f) for f in graph_paths.values()]
    else:
        overlap_graphs = [cov.get_overlap_graph(min_overlap=50) for cov in coverages]
        logger.info("Mapping breakpoints")
        set_breakpoint_ids(overlap_graphs, populations[1], outdir=args.outdir, threads=1)

        logger.info("Parent selection: by recombinations")
        select_by_recombinations(overlap_graphs)
        logger.info("Parent selection: by breakpoint")
        select_by_breakpoints(overlap_graphs)

        summarize_breakpoints(overlap_graphs)

        for graph in overlap_graphs:
            graph.write_pickle(graph_paths[graph["ref"]])
    
    logger.info(f"Clustering (feature: {args.clustering_feature})")
    clusters = cluster_phages(overlap_graphs, gamma=args.clustering_gamma, feature=args.clustering_feature, outdir=args.outdir)
    clusters = [c for c in clusters if len(c) > 1]

    if not clusters:
        logger.warning("No species cluster identified.")
        return

    logger.info("\n".join(",".join(c) for c in sorted(clusters, key=lambda x: -len(x))))

    logger.info(f"Interactive plot in {args.outdir}")
    display_genomes(overlap_graphs, clusters=clusters, norm=True, outdir=args.outdir)

if __name__ == '__main__':
    main()
