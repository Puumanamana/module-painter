from pathlib import Path
import argparse


ROOT_DIR = Path(__file__).resolve().parent.parent

def parse_args():
    '''
    Command line parser for ModularPainter
    '''

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--parents', type=str, help='Parent sequences (fasta)',
                        default=Path(ROOT_DIR, "tests", "fromagerie_1.fasta"))
    parser.add_argument('-c', '--children', type=str, help='Children sequences (fasta)',
                        default=Path(ROOT_DIR, "tests", "fromagerie_3.fasta"))
    parser.add_argument('--aln', type=str, help='Alignment tool', default='minimap2', choices=['blast', 'minimap2'])    
    parser.add_argument('--output', type=str, help='Output folder', default='/tmp/cedric/modular_painting_tests')    
    parser.add_argument('--min-id', type=float, default=0.8, help='Minimum sequence identity coverage')
    parser.add_argument('--min-module-size', type=int, default=20, help='Minimum size of a module/HSP')
    parser.add_argument('--arc-eq-diffs', type=int, default=20, help='Maximum distance between modules boundaries to consider them identical.')

    args = parser.parse_args()

    args.parents = Path(args.parents)
    args.children = Path(args.children)

    return args

