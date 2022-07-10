from pathlib import Path
import argparse
import re
import logging
import psutil


ROOT_DIR = Path(__file__).resolve().parent.parent
TEST_DIR = Path(ROOT_DIR, "tests/data")

def parse_args():
    """
    Command line parser for painting and simulation
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(help="sub-command help", dest="cmd")

    # Main parser
    main_parser = subparsers.add_parser("run", help="Paint set of child phages with parents",
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    main_parser.add_argument("-d", "--dataset", type=str, default="fromageries",
                        help="Test dataset choice")
    main_parser.add_argument("-p", "--populations", type=str, nargs="+",
                        help="All populations (fasta)")
    main_parser.add_argument("-c", "--children", type=str, nargs="+",
                        help="Children sequences within the population (list of glob pattern)")
    main_parser.add_argument("--exclude", type=str, nargs="+",
                        help="Files to exclude from the population (list of glob pattern)")
    main_parser.add_argument("--aligner", type=str, default="blastn", choices=["blastn", "minimap2"],
                        help="Alignment tool")    
    main_parser.add_argument("--outdir", type=str, default="/tmp/cedric/modular_painting/output",
                        help="Output folder")
    main_parser.add_argument("--min-length", type=int, default=5000,
                        help="Minimum contig length")
    main_parser.add_argument("--min-id", type=float, default=0.9,
                        help="Minimum sequence identity")
    main_parser.add_argument("--min-nw-id", type=float, default=0.8,
                        help="Minimum sequence identity for gap closing")
    main_parser.add_argument("--skip-nw", action="store_true",
                        help="Skip NW alignment for coverage refinement")
    main_parser.add_argument("--min-module-size", type=int, default=40,
                        help="Minimum size of a module/HSP")
    main_parser.add_argument("--arc-eq-diffs", type=int, default=10,
                        help="Maximum distance between modules boundaries to consider them identical.")
    main_parser.add_argument("--clustering-feature", type=str, default="breakpoint", choices=["breakpoint", "recombination"],
                        help="Feature to use to cluster phages")
    main_parser.add_argument("--clustering-gamma", type=float, default=0.2,
                        help="Cluster density")
    main_parser.add_argument("--resume", action="store_true",
                        help="Resume analysis if files already exist in output folder")
    main_parser.add_argument("--rename", action="store_true",
                        help="Rename contigs")
    main_parser.add_argument("--plot-fmt", type=str, default="pdf", choices=["html", "svg"],
                        help="Figure format")
    main_parser.add_argument("--threads", type=int, default=20,
                        help="Number of threads for alignment")

    # Simulation
    sim_parser = subparsers.add_parser("simulate", help="Simulate recombinant populations",
                                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sim_parser.add_argument("-m", "--n-modules", type=int, default=30)        
    sim_parser.add_argument("-n", "--n-forefathers", type=int, default=10)        
    sim_parser.add_argument("-k", "--n-subpopulations", type=int, default=3)
    sim_parser.add_argument("-r", "--n-rc", type=int, default=20)
    sim_parser.add_argument("-o", "--outdir", type=str, default=TEST_DIR)
    sim_parser.add_argument("--module-size-range", type=int, nargs=2, default=(200, 500))
    sim_parser.add_argument("--num-variants-range", type=int, nargs=2, default=(5, 10))        

    args = parser.parse_args()

    if args.cmd == "run":
        args = setup_populations(args)

    # Prepare outputs
    args.outdir = Path(args.outdir)
    args.outdir.mkdir(exist_ok=True, parents=True)

    return args

def setup_populations(args):
    if not (args.populations or args.children):
        args.populations = sorted(Path(TEST_DIR, args.dataset).glob("*.fasta"))

        if args.dataset == "fromageries":
            args.children = ["fromagerie_3"]
        elif "delong" in args.dataset:
            args.children = ["D117"]
        elif "sim" in args.dataset:
            args.children = ["children"]
        else:
            raise ValueError(f"Unknown dataset {args.dataset}")
    else:
        args.populations = [Path(p) for p in args.populations]

    if args.exclude:
        args.populations = [p for p in args.populations if not
                            any(re.match(".*"+e.strip("*")+".*", p.name) for e in args.exclude)]
        
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
            "{asctime} (Mem:{mem}) <{name}> {levelname}: {message}",
            "%H:%M:%S",
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

        if hasattr(mem, "pss"):
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

        record.mem = f"{mem/2**30:>5.1f} GB"

        return True
