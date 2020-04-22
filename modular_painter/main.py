from modular_painter.parser import parse_args
from modular_painter.split_communities import split_communities
from modular_painter.painting import paint_all
from modular_painter.breakpoints import handle_missing_data, group_breakpoints, cluster_breakpoints

def main():
    '''
    '''

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

    handle_missing_data(args.fasta[0], paintings, min_id=0.8, min_cov=0.8, threads=1)
    breakpoints = group_breakpoints(paintings, args.fasta[0],
                                    min_len=10, min_id=0.8, min_cov=0.8, threads=1)
    cluster_breakpoints(breakpoints)

if __name__ == '__main__':
    main()
