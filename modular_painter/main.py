from modular_painter.parser import parse_args
from modular_painter.blast_genomes import get_raw_coverage
# from modular_painter.minimap_genomes import minimap

from modular_painter.breakpoints import handle_missing_data, group_breakpoints
from modular_painter.clustering import select_parents_by_sparsity, cluster_breakpoints
from modular_painter.display import display_genomes


def main():

    args = parse_args()

    coverages = get_raw_coverage(
        *args.fasta,
        output=args.output,
        extend_min_dist=args.min_module_size,
        min_id=args.min_id,
        min_cov=args.min_cov
    )
    # At this point, for each coverage, all parents are distant of at least {min_module_size}
    # and should be looping around if the start and end are the same
    # This also takes care of potential embedding relationship between hits, which helps for later
    
    for coverage in coverages:
        print(f"======== Processing: {coverage.target} ========")
        # merge arcs from the same parent if they are closer than {min_module_size} and
        # remove embedded arcs (from the same parent)
        coverage.extend_arcs_per_parent(args.min_module_size)
        # sticky boundaries
        coverage.sync_boundaries(args.arc_eq_diffs, which="start")
        coverage.sync_boundaries(args.arc_eq_diffs, which="end")
        # Boundaries have changed a bit, we check if we have embeddings or fusing to do
        coverage.extend_arcs_per_parent(args.min_module_size)
        # simplify coverage by grouping equal intervals
        coverage.merge_equal_intervals()
        # fill all gaps
        coverage.fill_gaps(args.min_module_size)
        # Remove embedded intervals
        coverage.simplify_embedded()
        # Lee and Lee
        coverage.get_minimal_coverage()
        # Cleanup
        coverage.remove_flagged()

        if coverage.target in "RSWXdZa":
            print(coverage)

    handle_missing_data(args.fasta[1], coverages, min_id=0.8,
                        min_cov=0.8, threads=1)

    breakpoints = group_breakpoints(coverages, args.fasta[1],
                                    min_len=10, min_id=0.8, min_cov=0.8,
                                    threads=1)
    bk_unique = select_parents_by_sparsity(breakpoints)

    # For clustering and display, make parents unique in coverage
    for (target, i1, i2, parents, _, _, _) in bk_unique.to_numpy():
        (l, r) = (paintings[target].data[i1], paintings[target].data[i2])

        if not parents[0] in l.data.index:
            (l, r) = (r, l)

        l.data = l.data.loc[[parents[0]]]
        r.data = r.data.loc[[parents[1]]]

    print(*list(paintings.values()))

    clusters = cluster_breakpoints(bk_unique, gamma=0.75)

    display_genomes(coverages, clusters=clusters, norm=True)
    

if __name__ == '__main__':
    main()
