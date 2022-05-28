import logging
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from igraph import Graph

from module_painter.parser import parse_args
from module_painter.util import concatenate_fasta
from module_painter.coverage import Coverage
from module_painter.wrapper import blastn, minimap2
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
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(levelname)s:%(message)s',
                        datefmt='%H:%M:%S')

    args = parse_args()

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
        logging.info(f"Initial alignment with {args.aligner}")
        if "blast" in args.aligner:
            alns = blastn(*populations, args.outdir, **ALN_PARAMS["blastn"])
        else:
            alns = minimap2(*populations, r=args.min_module_size, t=args.threads,
                            **ALN_PARAMS["minimap2"])
        alns.to_csv(aln_file)

    logging.debug(f"{sum(alns.pident>=args.min_id):,} alignments remaining (pident>={args.min_id:.1%})")
    alns = alns[alns.pident >= args.min_id]
    # For each query, keep hits with sstrand of the best hit
    sstrands = alns.groupby(["sacc", "qacc"]).apply(lambda df: df.set_index("sstrand").nident.idxmax()).to_dict()
    alns = alns.groupby(["sacc", "qacc"], group_keys=False).apply(lambda df: df[df.sstrand==sstrands[df.name]])

    if alns.empty:
        logging.error(f"No alignments found with more than 2 parents between parents and children")
        return

    n_children = alns.sacc.nunique()
    alns = alns.groupby("sacc").filter(lambda x: x.qacc.nunique() > 1)
    logging.info(f"{n_children-alns.sacc.nunique():,} children discarded (<2 parents found)")

    seq_data = {seq.id: seq.seq for fasta in populations for seq in SeqIO.parse(fasta, "fasta")}

    logging.info("Applying coverage rules for all children")
    coverages = []
    for child in alns.sacc.unique():
        aln = alns[alns.sacc==child]
        if aln.qacc.nunique() == 1:
            logging.warning(f"Too few parents cover {child}. Skipping")
            continue
        logging.debug(f"Processing {child}")
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
        if len(coverage) < 2:
            continue
        # Lee and Lee
        coverage.get_minimal_coverage()
        coverages.append(coverage)
    
    if not coverages:
        logging.warning("No shared species between parents and children. Aborting.")
        return

    tmp = {c.ref: c for c in coverages}

    logging.info("Mapping missing parents")
    map_missing_parents(populations[1], coverages, outdir=args.outdir, threads=1)

    for c in coverages:
        # print(c)
        logging.debug(c.show())
        # logging.debug(TRUTH.get(c.ref))

    logging.info(f"Remaining: {[c.ref for c in coverages]}")
    
    graph_dir = Path(args.outdir, "overlap_graphs")
    graph_dir.mkdir(exist_ok=True)
    graph_paths = {cov.ref: Path(graph_dir, f"{cov.ref}.pickle") for cov in coverages}

    if args.resume and all(f.is_file() for f in graph_paths.values()):
        overlap_graphs = [Graph.Read_Pickle(f) for f in graph_paths.values()]
    else:
        overlap_graphs = [cov.get_overlap_graph(min_overlap=50) for cov in coverages]
        logging.info("Mapping breakpoints")
        set_breakpoint_ids(overlap_graphs, populations[1], outdir=args.outdir, threads=1)

        logging.info("Parent selection: by recombinations")
        select_by_recombinations(overlap_graphs)
        logging.info("Parent selection: by breakpoint")
        select_by_breakpoints(overlap_graphs)

        summarize_breakpoints(overlap_graphs)

        for graph in overlap_graphs:
            graph.write_pickle(graph_paths[graph["ref"]])
    
    logging.info(f"Clustering (feature: {args.clustering_feature})")
    clusters = cluster_phages(overlap_graphs, gamma=args.clustering_gamma, feature=args.clustering_feature)
    clusters = [c for c in clusters if len(c) > 1]

    if not clusters:
        logging.warning("No species cluster identified.")
        return

    logging.info("\n".join(",".join(c) for c in sorted(clusters, key=lambda x: -len(x))))

    logging.info(f"Interactive plot in {args.outdir}")
    display_genomes(overlap_graphs, clusters=clusters, norm=True, outdir=args.outdir)

if __name__ == '__main__':
    main()
