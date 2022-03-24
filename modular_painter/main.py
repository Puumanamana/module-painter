from modular_painter.parser import parse_args
from modular_painter.coverage import Coverage
from modular_painter.alignment import align
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

def main():

    args = parse_args()
    
    # Raw alignment
    print(f"Initial alignment with {args.aln}")
    alns = align(args.parents, args.children, algo=args.aln, output=args.output)

    coverages = []
    for child in "WSRXdZaTUVYbc":
        print(f"Simplification ({child})")
        aln = alns[alns.sacc==child]
        if args.use_ground_truth:
            aln = alns[alns.qacc.isin(set(TRUTH[child]))]
        
        coverage = Coverage.from_pandas(aln, size=aln.slen.iloc[0])
        # sticky boundaries
        coverage.sync_boundaries(args.arc_eq_diffs, ["start"])
        coverage.sync_boundaries(args.arc_eq_diffs, ["end"])
        # simplify coverage by grouping equal intervals
        coverage.merge_equal_intervals()
        # Fill all gaps
        coverage.fill_gaps(args.min_module_size)
        # Remove embedded intervals
        coverage.simplify_embedded()
        # Lee and Lee
        coverage.get_minimal_coverage()
        # Cleanup
        coverage.remove_flagged()
        coverages.append(coverage)

    print("\nMapping missing parents")
    map_missing_parents(args.children, coverages, outdir=args.output, threads=1)

    for c in coverages:
        print(c)
        print(TRUTH.get(c.ref))
        # import ipdb;ipdb.set_trace()

    overlap_graphs = {cov.ref: cov.get_overlap_graph(min_overlap=50) for cov in coverages}
    print("Mapping breakpoints")
    set_breakpoint_ids(overlap_graphs, args.children, outdir=args.output, threads=1)

    print("Parent selection: by recombinations")
    select_by_recombinations(overlap_graphs)
    print("Parent selection: by breakpoint")
    select_by_breakpoints(overlap_graphs)

    print(f"Clustering (feature: {args.clustering_feature})")
    clusters = cluster_phages(overlap_graphs, gamma=0.5, feature=args.clustering_feature)

    print("Interactive plot")
    display_genomes(overlap_graphs, clusters=clusters, norm=True)

if __name__ == '__main__':
    main()
