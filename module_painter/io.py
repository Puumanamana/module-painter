from pathlib import Path
import argparse
import re
import logging
import psutil


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

    if not args.children:
        raise ValueError(f"No children found")

    setattr(args, "parents", [f for f in args.populations if f not in args.children])

    args.parents = [Path(p) for p in args.parents]
    args.children = [Path(c) for c in args.children]

    return args

def setup_logger(name, log_file, level=logging.INFO):
    """
    Setup logging if not set, or return logger if already exists
    Args:
        name (str): name of logger
        log_file (str): path to save logs
        level (int): log level for stderr
    Returns:
        logging.Logger
    """

    # Create the Logger
    logger = logging.getLogger(name)
    logger.addFilter(MemoryTracer())

    logger.setLevel(logging.DEBUG)

    if not any(isinstance(hdl, logging.FileHandler) for hdl in logger.handlers):
        if logger.hasHandlers():
            logger.handlers.clear()
        logger.propagate = False

        # Create a Formatter for formatting the log messages
        formatter = logging.Formatter(
            '{asctime} (Mem:{mem}) <{name}> {levelname}: {message}',
            '%H:%M:%S',
            style="{"
        )

        # Create the Handler for logging data to a file
        Path(log_file.parent).mkdir(exist_ok=True)
        logger_handler = logging.FileHandler(str(log_file), mode="w+")
        logger_handler.setLevel(logging.DEBUG)
        logger_handler.setFormatter(formatter)
        logger.addHandler(logger_handler)

        # Create the Handler for logging data to console.
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.getLevelName(level))
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    return logger

class MemoryTracer(logging.Filter):
    """
    To track memory usage at different steps
    Used memory is computed as the PSS of the main program
    + the sum of the PSS of its children for linux, and RSS
    + the sum of the USS of its children for MacOS
    """

    def filter(self, record):
        process = psutil.Process()
        mem = process.memory_full_info()

        if hasattr(mem, 'pss'):
            mem = mem.pss
        else: # No PSS info for MacOS
            mem = mem.rss

        for child in process.children(recursive=True):
            try:
                mem += child.memory_full_info().pss
            except AttributeError: # No PSS info for MacOS
                mem += child.memory_full_info().uss
            except psutil.NoSuchProcess:
                pass
            except psutil.AccessDenied:
                pass

        record.mem = f'{mem/2**30:>5.1f} GB'

        return True
