from pathlib import Path
import argparse
import re


ROOT_DIR = Path(__file__).resolve().parent.parent
TEST_DIR = Path(ROOT_DIR, "tests")

def parse_args():
    """
    Command line parser for ModularPainter
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-d", "--dataset", type=str, default="fromageries",
                        help="Test dataset choice")
    parser.add_argument("-p", "--populations", type=str, nargs="+",
                        help="All populations (fasta)")
    parser.add_argument("-c", "--children", type=str, nargs="+",
                        help="Children sequences within the population (list of glob pattern)")
    parser.add_argument("--aligner", type=str, default="blastn", choices=["blastn", "minimap2"],
                        help="Alignment tool")    
    parser.add_argument("--outdir", type=str, default="/tmp/cedric/modular_painting/output",
                        help="Output folder")
    parser.add_argument("--min-length", type=int, default=5000,
                        help="Minimum contig length")
    parser.add_argument("--min-id", type=float, default=0.9,
                        help="Minimum sequence identity")
    parser.add_argument("--min-nw-id", type=float, default=0.8,
                        help="Minimum sequence identity for gap closing")
    parser.add_argument("--min-module-size", type=int, default=40,
                        help="Minimum size of a module/HSP")
    parser.add_argument("--arc-eq-diffs", type=int, default=15,
                        help="Maximum distance between modules boundaries to consider them identical.")
    parser.add_argument("--clustering-feature", type=str, default="breakpoint", choices=["breakpoint", "recombination"],
                        help="Feature to use to cluster phages")
    parser.add_argument("--clustering-gamma", type=float, default=0.5,
                        help="Cluster density")
    parser.add_argument("--use-ground-truth", action="store_true",
                        help="argument to remove later after the paper is published")
    parser.add_argument("--resume", action="store_true",
                        help="Resume analysis if files already exist in output folder")
    parser.add_argument("--rename", action="store_true",
                        help="Rename contigs")
    parser.add_argument("--threads", type=int, default=20,
                        help="Number of threads for alignment")
    args = parser.parse_args()

    args.outdir = Path(args.outdir)
    args.outdir.mkdir(exist_ok=True)

    if not (args.populations or args.children):
        args.populations = sorted(Path(TEST_DIR, args.dataset).glob("*.fasta"))

        if args.dataset == "fromageries" and "sim" not in args.dataset:
            args.children = ["fromagerie_3"]
        elif "sim" in args.dataset:
            args.children = ["children"]
        elif "delong" in args.dataset:
            args.rename = True
            args.children = ["D025"]
        else:
            raise ValueError(f"Unknown dataset {args.dataset}")
    else:
        args.populations = [Path(p) for p in args.populations]

    args.children = [f for f in args.populations
                     if any(re.match(".*"+ri.strip("*")+".*", f.name)
                            for ri in args.children)]

    setattr(args, "parents", [f for f in args.populations if f not in args.children])

    args.parents = [Path(p) for p in args.parents]
    args.children = [Path(c) for c in args.children]

    return args

