from pathlib import Path
import pandas as pd

from module_painter.io import parse_args, setup_logger
from module_painter.painter import paint
from module_painter.simulation import simulate
from module_painter.clustering import cluster_phages
from module_painter.display import display_phages_hv, display_phages_plt

def main():
    args = parse_args()

    logger = setup_logger('module-painter', Path(args.outdir, 'log.txt'))
    logger.info(f"Starting module-painter {args.cmd}")

    #========================#
    #====== CLI parsing =====#
    #========================#
    logger.debug("==== Parameters ====")
    for k,v in vars(args).items():
        logger.debug(f"{k}={v}")

    #========================#
    #======= Simulate =======#
    #========================#
    if args.cmd == "simulate":
        simulate(**vars(args))
        return
    
    #========================#
    #========== Run =========#
    #========================#
    paintings = []
    if args.rotate_parent:
        logger.info("Rotating parent mode")
        outdir = args.outdir
        for parent in args.populations:
            logger.info(f"==> Setting {parent} as parent population")
            args.parents = [parent]
            args.children = [p for p in args.populations if p!=parent]
            args.outdir = Path(outdir, f"parent-{parent.stem}")
            args.outdir.mkdir(exist_ok=True)
            paintings.append(paint(**vars(args)))

        args.outdir = outdir

        # Aggregate results for all runs
        ## Recombinations
        rc = pd.concat([p["rc"] for p in paintings])
        rc.to_csv(f"{args.outdir}/recombinations.csv")        
        
        ## Shared breakpoints/recombinations
        all_links = [p["links"] for p in paintings if not p["links"].empty]
        if not all_links:
            return
        links = pd.concat(all_links)
        links.to_csv(f"{args.outdir}/links.csv")
    else:
        painting = paint(**vars(args))

        if not painting: return
        
        paintings.append(painting)
        links = paintings[0]["links"]
        links.to_csv(f"{args.outdir}/links.csv")

    #========================#
    #======= Clustering =====#
    #========================#
    clusters = cluster_phages(
        links, outdir=args.outdir,
        method=args.clustering_method,
        group_pattern=args.group_pattern,
        resolution=args.resolution
    )
    
    clusters = [x for x in clusters if len(x) > 1]

    if not clusters:
        logger.warning("No species cluster identified.")
        return

    #========================#
    #======== Display =======#
    #========================#
    if not args.rotate_parent and args.plot_coverages:
        logger.info(f"Interactive plot in {args.outdir}")

        if args.plot_fmt in {"html", "svg"}:
            display_phages_hv(paintings[0]["graphs"], clusters=clusters,
                              norm=True, outdir=args.outdir, fmt=args.plot_fmt)
        else:
            display_phages_plt(paintings[0]["graphs"], clusters=clusters,
                               norm=True, outdir=args.outdir, fmt=args.plot_fmt)
    
