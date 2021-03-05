from modular_painter.parser import parse_args
from modular_painter.split_communities import split_communities
from modular_painter.painting import paint_all
from modular_painter.breakpoints import (
    handle_missing_data,
    group_breakpoints,
    select_parents,
    cluster_breakpoints
    )
from modular_painter.display import display_linear_genome

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
    bk_unique = select_parents(breakpoints)

    # For display, make parents unique in coverage
    for (parents, start, end, target, _) in bk_unique.to_numpy():
        paintings[target].change_parent(start, end, *parents)

    clusters = cluster_breakpoints(bk_unique)

    for cluster in clusters:
        display_linear_genome([painting[c] for c in cluster])

if __name__ == '__main__':
    main()
