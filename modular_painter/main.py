from parser import parse_args
from split_communities import split_communities
from painting import paint_all
from clustering import get_subspecies

def main():
    '''
    '''

    args = parse_args()

    species = split_communities(*args.fasta, output=args.output,
                                min_id=args.min_id, min_cov=args.min_cov, min_hsp_len=args.min_hsp_len,
                                show=args.show)

    paintings = paint_all(args.fasta, species, max_merge_dist=args.max_merge_dist)

    get_subspecies(paintings)

if __name__ == '__main__':
    main()
