from pathlib import Path
import argparse

def parse_args():
    '''
    Command line parser for ModularPainter
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str, nargs='+', help='List of path to viral samples (fasta format)', default=['../test/population-B_117m_simple.fasta', '../test/population-C_250m_simple.fasta'])
    parser.add_argument('--output', type=str, help='Output folder', default='/tmp/cedric/modular_painting_tests')    
    parser.add_argument('--min-id', type=float, default=0.8, help='Minimum sequence identity for species definition')
    parser.add_argument('--min-cov', type=float, default=0.5, help='Minimum sequence coverage for species definition')
    parser.add_argument('--min-hsp-len', type=int, default=30, help='Minimum HSP length in BLAST alignments')
    parser.add_argument('--max-merge-dist', type=int, default=3, help='Maximum distance between modules to consider them contiguous.')        
    parser.add_argument('--show', action='store_true', default=False)

    args = parser.parse_args()

    args.fasta = [Path(fasta) for fasta in args.fasta]

    return args

