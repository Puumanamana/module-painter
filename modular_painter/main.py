from parser import parse_args
from split_communities import split_communities
from painting import paint_all

def main():
    '''
    '''

    args = parse_args()

    species = split_communities(*args.fasta, output=args.output,
                                min_id=args.min_id, min_cov=args.min_cov,
                                show=args.show)

    paintings = paint_all(args.fasta, species)

if __name__ == '__main__':
    main()
