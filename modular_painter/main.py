import logging
from pathlib import Path
import pandas as pd
from Bio import SeqIO

from modular_painter.parser import parse_args
from modular_painter.util import concatenate_fasta
from modular_painter.coverage import Coverage
from modular_painter.wrapper import blastn, minimap2
from modular_painter.breakpoints import map_missing_parents, set_breakpoint_ids
from modular_painter.parent_selection import select_by_recombinations, select_by_breakpoints
from modular_painter.clustering import cluster_phages
from modular_painter.display import display_genomes


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
    minimap2={"s": 100, "z": "100,100", "N": 50,
              "no-long-join": True, "c": True, "P": True},
    blastn={"gapopen": 5, "gapextend": 2}
)

def main():
    logging.basicConfig(level=logging.DEBUG,
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

    logging.info(f"{sum(alns.pident<args.min_id):,} alignments discarded (pident<{args.min_id:.1%})")
    alns = alns[alns.pident >= args.min_id]
    n_children = alns.sacc.nunique()
    alns = alns.groupby("sacc").filter(lambda x: x.qacc.nunique() > 1)
    logging.info(f"{n_children-alns.sacc.nunique():,} children discarded (< 2 parents found)")

    if alns.empty:
        logging.error(f"No alignments found with {args.aligner} between parents and children")
        return

    seq_data = {seq.id: seq.seq for fasta in populations for seq in SeqIO.parse(fasta, "fasta")}

    coverages = []
    for child in alns.sacc.unique():
        aln = alns[alns.sacc==child]
        if aln.qacc.nunique() == 1:
            logging.warning(f"Too few parents cover {child}. Skipping")
            continue
        logging.info(f"Processing {child}")
        if args.use_ground_truth:
            aln = alns[alns.qacc.isin(set(TRUTH[child]))]

        coverage = Coverage.from_pandas(aln, size=aln.slen.iloc[0])
        # sticky boundaries
        coverage.sync_boundaries("start", args.arc_eq_diffs)
        coverage.sync_boundaries("end", args.arc_eq_diffs)
        # Fuse intervals from the same parents if they are close
        coverage.fuse_close_modules(seq_data, args.min_module_size, args.min_nw_id)
        # simplify coverage by grouping equal intervals
        coverage.merge_equal_intervals()
        # Remove embedded intervals
        coverage.simplify_embedded()
        # Fill all gaps
        coverage.fill_gaps(args.min_module_size)
        # Lee and Lee
        coverage.get_minimal_coverage()

        if len(coverage) > 1:
            coverages.append(coverage)

    if not coverages:
        logging.warning("No shared species between parents and children. Aborting.")
        return

    logging.info("Mapping missing parents")
    map_missing_parents(populations[1], coverages, outdir=args.outdir, threads=1)

    for c in coverages:
        print(c)
        logging.debug(c.show())
        logging.debug(TRUTH.get(c.ref))

    overlap_graphs = {cov.ref: cov.get_overlap_graph(min_overlap=50) for cov in coverages}
    logging.info("Mapping breakpoints")
    set_breakpoint_ids(overlap_graphs, populations[1], outdir=args.outdir, threads=1)

    logging.info("Parent selection: by recombinations")
    select_by_recombinations(overlap_graphs)
    logging.info("Parent selection: by breakpoint")
    select_by_breakpoints(overlap_graphs)

    logging.info(f"Clustering (feature: {args.clustering_feature})")
    clusters = cluster_phages(overlap_graphs, gamma=0.5, feature=args.clustering_feature)
    clusters = [c for c in clusters if len(c) > 1]

    if not clusters:
        logging.warning("No species cluster identified.")
        return

    logging.debug([list(c) for c in sorted(clusters, key=lambda x: -len(x))])

    logging.info("Interactive plot")
    display_genomes(overlap_graphs, clusters=clusters, norm=True)

if __name__ == '__main__':
    main()
