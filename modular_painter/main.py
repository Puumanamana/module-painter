from modular_painter.parser import parse_args
from modular_painter.split_communities import split_communities
from modular_painter.painting import paint_all
from modular_painter.breakpoints import (
    handle_missing_data,
    group_breakpoints,
    select_parents_by_sparsity,
    cluster_breakpoints
    )
from modular_painter.display import display_genomes

def main():

    args = parse_args()

    species = split_communities(*args.fasta,
                                output=args.output,
                                min_id=args.min_id,
                                min_cov=args.min_cov,
                                min_module_size=args.min_module_size,
                                show=args.show)

    paintings = paint_all(args.fasta, species,
                          min_module_size=args.min_module_size,
                          arc_eq_diffs=args.arc_eq_diffs)

    handle_missing_data(args.fasta[0], paintings, min_id=0.8,
                        min_cov=0.8, threads=1)

    breakpoints = group_breakpoints(paintings, args.fasta[0],
                                    min_len=10, min_id=0.8, min_cov=0.8,
                                    threads=1)
    bk_unique = select_parents_by_sparsity(breakpoints)

    # For display, make parents unique in coverage
    for (target, i1, i2, parents, _, _, _) in bk_unique.to_numpy():
        (l, r) = (paintings[target].data[i1], paintings[target].data[i2])

        if not parents[0] in l.data.index:
            (l, r) = (r, l)

        l.data = l.data.loc[[parents[0]]]
        r.data = r.data.loc[[parents[1]]]

    print(*list(paintings.values()))
    display_genomes(paintings)
    import ipdb;ipdb.set_trace()
    clusters = cluster_breakpoints(bk_unique)

if __name__ == '__main__':
    main()
