from pathlib import Path
import pandas as pd

from module_painter.io import parse_args, setup_logger
from module_painter.painter import paint
from module_painter.simulation import simulate
from module_painter.clustering import cluster_phages
from module_painter.display import display_phages

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
            logger.info(f"Setting {parent} as parent population")
            args.parents = [parent]
            args.children = [p for p in args.populations if p!=parent]
            args.outdir = Path(outdir, f"parent-{parent.stem}")
            args.outdir.mkdir(exist_ok=True)
            paintings.append(paint(**vars(args)))
        args.outdir = outdir
        links = pd.concat([p["links"] for p in paintings])
        # save to file
        links_str = links.copy().sort_values(by="sacc", key=lambda x: x.str.len(), ascending=False)
        links_str.sacc = links_str.sacc.apply("-".join)
        links_str.to_csv(f"{args.outdir}/links-combined.csv")
    else:
        painting = paint(**vars(args))

        if not painting: return
        
        paintings.append(painting)
        links = paintings[0]["links"]

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
        display_phages(paintings[0]["graphs"], clusters=clusters,
                        norm=True, outdir=args.outdir, fmt=args.plot_fmt)
    
