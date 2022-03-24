from pathlib import Path
import argparse


ROOT_DIR = Path(__file__).resolve().parent.parent
TEST_DIR = Path(ROOT_DIR, "tests")

def parse_args():
    '''
    Command line parser for ModularPainter
    '''

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--parents', type=str, default=Path(TEST_DIR, "fromagerie_1.fasta"),
                        help='Parent sequences (fasta)')
    parser.add_argument('-c', '--children', type=str, default=Path(TEST_DIR, "fromagerie_3.fasta"),
                        help='Children sequences (fasta)')
    parser.add_argument('--aln', type=str, default='minimap2', choices=['blast', 'minimap2'],
                        help='Alignment tool')    
    parser.add_argument('--output', type=str, default='/tmp/cedric/modular_painting_tests',
                        help='Output folder')
    parser.add_argument('--min-id', type=float, default=0.8,
                        help='Minimum sequence identity coverage')
    parser.add_argument('--min-module-size', type=int, default=20,
                        help='Minimum size of a module/HSP')
    parser.add_argument('--arc-eq-diffs', type=int, default=15,
                        help='Maximum distance between modules boundaries to consider them identical.')
    parser.add_argument('--clustering-feature', type=str, default='breakpoint', choices=['breakpoint', 'recombination'],
                        help='Feature to use to cluster phages')
    parser.add_argument('--use-ground-truth', action="store_true",
                        help='argument to remove later after the paper is published')
    

    args = parser.parse_args()

    args.parents = Path(args.parents)
    args.children = Path(args.children)

    return args

