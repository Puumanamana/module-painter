from parser import parse_args
from split_communities import split_communities

def main():
    '''
    '''

    args = parse_args()

    split_communities(*args.fasta, output=args.output,
                      min_id=args.min_id, min_cov=args.min_cov)

if __name__ == '__main__':
    main()
