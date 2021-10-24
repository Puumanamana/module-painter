from pathlib import Path
import argparse


ROOT_DIR = Path(__file__).resolve().parent.parent
TEST_FA = Path(ROOT_DIR, 'tests').glob('*.fasta')


def parse_args():
    '''
    Command line parser for ModularPainter
    '''

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--fasta', type=str, nargs='+', help='List of path to viral samples (fasta format)', default=list(TEST_FA))
    parser.add_argument('--output', type=str, help='Output folder', default='/tmp/cedric/modular_painting_tests')    
    parser.add_argument('--min-id', type=float, default=0.8, help='Minimum sequence identity for species definition')
    parser.add_argument('--min-cov', type=float, default=0.5, help='Minimum sequence coverage for species definition')
    parser.add_argument('--min-module-size', type=int, default=50, help='Minimum size of a module/HSP')
    parser.add_argument('--arc-eq-diffs', type=int, default=10, help='Maximum distance between embedded modules to consider them identical.')
    parser.add_argument('--show', action='store_true', default=False)

    args = parser.parse_args()

    args.fasta = [Path(fasta) for fasta in args.fasta]

    return args

