from modular_painter.parser import parse_args
from modular_painter.coverage import Coverage
from modular_painter.alignment import align
from modular_painter.breakpoints import handle_missing_data, set_breakpoint_ids
from modular_painter.parent_selection import select_parents
from modular_painter.clustering import cluster_breakpoints
from modular_painter.display import display_genomes


TRUTH = dict(
    R="MAJA",
    S="(JM)FJA",
    W="JFJA",
    X="(JM)F(AJ)MA",
    d="AMJ",
    Z="B-B",
    a="B-B",
    T="LB-B--L(BCEFHIK)-B--(BDG)-(BDGHK)-(ABCDEFGHIJKL)B(AJM)L",
    U="LB-B-BL(BCEFHIK)-B--(BDG)-(BDGHK)-(ABCDEFGHIJKL)B(AJM)L",
    V="LB-B--BL(BCEFHIK)-B--B(AJM)L"
)

def main():

    args = parse_args()
    
    # Raw alignment
    alns = align(args.parents, args.children, algo=args.aln, output=args.output)

    coverages = []
    # for child in "RSWXdZa": # for testing purposes, change later
    for child in "RSWXdZaTUVYbc":
        print(f"======== Processing: {child} ========")
        # convert to coverage object
        aln = alns[alns.sacc==child]
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

    handle_missing_data(args.children, coverages, outdir=args.output, threads=1)

    for c in coverages:
        print(c)
        print(TRUTH.get(c.ref))

    overlap_graphs = {cov.ref: cov.get_overlap_graph(min_overlap=50)
                      for cov in coverages}
    set_breakpoint_ids(overlap_graphs, args.children, outdir=args.output, threads=1)
    
    bk_unique = select_parents(overlap_graphs)

    coverages = {cov.ref: cov for cov in coverages}

    # For clustering and display, make parents unique in coverage
    for (parents, _, _ , i1, i2, ref, _) in bk_unique.to_numpy():
        (left, right) = (coverages[ref].arcs[i1], coverages[ref].arcs[i2])

        if not parents[0] in left.meta:
            (left, right) = (right, left)
            if not parents[0] in left.meta:
                print(coverages[ref])
                print(bk_unique[bk_unique.ref==ref])
                print(parents)
                print(left, right)
                import ipdb;ipdb.set_trace()
        left.meta = {parents[0]}
        right.meta = {parents[1]}

    clusters = cluster_breakpoints(bk_unique, gamma=0.75)
    
    display_genomes(coverages, clusters=clusters, norm=True)
    

if __name__ == '__main__':
    main()
