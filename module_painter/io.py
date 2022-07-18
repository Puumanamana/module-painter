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
    main_parser.add_argument("--threads", type=int, default=10,
                        help="Number of threads for alignment")
    main_parser.add_argument("--resume", action="store_true",
                        help="Resume analysis if files already exist in output folder")
    main_parser.add_argument("--rename", action="store_true",
                        help="Rename contigs")
    main_parser.add_argument("--min-length", type=int, default=5000,
                        help="Minimum contig length")
    io_args = main_parser.add_argument_group('I/O')
    io_args.add_argument("-d", "--dataset", default="fromageries",
                        help="Test dataset choice")
    io_args.add_argument("-p", "--populations", nargs="+",
                        help="All populations (fasta)")
    io_args.add_argument("-c", "--children", nargs="*",
                        help="Children sequences within the population (list of glob pattern)")
    io_args.add_argument("--exclude", nargs="+",
                        help="Files to exclude from the population (list of glob pattern)")
    io_args.add_argument("--outdir", default="/tmp/cedric/modular_painting/output",
                        help="Output folder")
    aln_args = main_parser.add_argument_group('Alignment')
    aln_args.add_argument("--aligner", default="minimap2", choices=["blastn", "minimap2"],
                        help="Alignment tool")    
    aln_args.add_argument("--min-id", type=float, default=0.9,
                        help="Minimum sequence identity")
    aln_args.add_argument("--min-nw-id", type=float, default=0.9,
                        help="Minimum sequence identity for gap closing")
    aln_args.add_argument("--skip-nw", action="store_true",
                        help="Skip NW alignment for coverage refinement")
    aln_args.add_argument("--min-module-size", type=int, default=100,
                        help="Minimum size of a module/HSP")
    aln_args.add_argument("--arc-eq-diffs", type=int, default=15,
                        help="Maximum distance between modules boundaries to consider them identical.")
    cluster_args = main_parser.add_argument_group('Clustering')
    cluster_args.add_argument("--clustering-feature", default="recombination", choices=["breakpoint", "recombination"],
                        help="Feature to use to cluster phages")
    cluster_args.add_argument("--clustering-method", default="leiden", choices=["connected_components", "leiden", "infomap"],
                             help="Phage clustering method")
    cluster_args.add_argument("--resolution", type=float, default=0.8,
                              help="Cluster density (CPM threshold for community detection for Leiden method)")
    cluster_args.add_argument("--rotate-parent", action="store_true",
                        help="Cluster datasets by rotating parent set")
    plot_args = main_parser.add_argument_group('Plotting')
    plot_args.add_argument("--plot-coverages", action="store_true",
                        help="Display coverage for phages")
    plot_args.add_argument("--plot-fmt", default="html",
                        help="Figure format")
    plot_args.add_argument("--group-pattern", default=None,
                        help="Regex pattern for coloring groups for recombination graph")

    # Simulation
    sim_parser = subparsers.add_parser("simulate", help="Simulate recombinant populations",
                                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    genetic_args = sim_parser.add_argument_group('Genetic pool')
    genetic_args.add_argument("--n-modules", type=int, default=30)
    genetic_args.add_argument("--n-variants-range", type=int, nargs=2, default=(2, 5))
    genetic_args.add_argument("--module-size-range", type=int, nargs=2, default=(200, 500))

    population_args = sim_parser.add_argument_group('Populations')
    population_args.add_argument("--n-subpop", type=int, default=3)
    population_args.add_argument("--subpop-size-lam", type=int, default=3)
    population_args.add_argument("--n-gen-range", type=int, nargs=2, default=(1, 6))
    population_args.add_argument("--rate-range", type=float, nargs=2, default=(0.1, 0.3))
    sim_parser.add_argument("--outdir", default=TEST_DIR)

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

    if not args.rotate_parent:
        args.children = [f for f in args.populations
                         if any(re.match(".*"+ri.strip("*")+".*", f.name)
                                for ri in args.children)]

        if not args.children:
            raise ValueError(f"No children found")

        setattr(args, "parents", [f for f in args.populations if f not in args.children])

        args.parents = [Path(p) for p in args.parents]
        args.children = [Path(c) for c in args.children]

    return args

def find_outputs(outdir, aligner="minimap2"):
    interm_dir = Path(outdir, "interm")
    interm_dir.mkdir(exist_ok=True)
    
    folders = {
        key: Path(outdir, "interm", name) for (key, name) in [
            ("raw_cov", "simplified_coverage"),
            ("@mapped", "missing_data"),
            ("graphs", "overlap_graphs")
        ]
    }
    outputs = dict(raw_aln=Path(interm_dir, aligner).with_suffix(".csv"))

    for (key, folder) in folders.items():
        outputs[key] = {}
        folder.mkdir(exist_ok=True)
        
        ext = "pickle" if "graph" in key else "csv"
        for fpath in Path(folder).glob(f"*.{ext}"):
            outputs[key][fpath.stem] = fpath

    return (folders, outputs)

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
