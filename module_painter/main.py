from pathlib import Path

from module_painter.io import parse_args, setup_logger
from module_painter.painter import paint
from module_painter.simulation import simulate

def main():
    args = parse_args()

    logger = setup_logger('module-painter', Path(args.outdir, 'log.txt'))
    
    if args.cmd == "run":
        paint(**vars(args))
    elif args.cmd == "simulate":
        simulate(**vars(args))
    else:
        raise ValuError(f"Unknown cmd {cmd}")
