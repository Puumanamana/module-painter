import argparse
from copy import deepcopy
from pathlib import Path

import numpy as np
from Bio import SeqIO


np.random.seed(123)

TEST_DIR = Path(Path(__file__).resolve().parent.parent, "tests", "sim")
NUCL = list("ACGT")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--id', type=int, default=0)
    parser.add_argument('--sim-dir', type=str, default=TEST_DIR)
    parser.add_argument('--module-size-range', type=int, nargs=2, default=(200, 500))
    parser.add_argument('--num-variants-range', type=int, nargs=2, default=(5, 10))        
    parser.add_argument('--n-modules', type=int, default=30)        
    parser.add_argument('--n-forefathers', type=int, default=10)        
    parser.add_argument('--n-subpopulations', type=int, default=3)
    parser.add_argument('--n-rc', type=int, default=20)
    args = parser.parse_args()

    args.sim_dir = Path(args.sim_dir).with_suffix(f".{args.id}")
    args.sim_dir.mkdir(exist_ok=True)

    return args

class color:
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\u001b[0m'

    def get_color(i):
        return f"\u001b[1;38;5;{i}m"

    def display(x, n=255, bold=False, underlined=False):
        prefix = color.get_color(n)
        if bold:
            prefix += color.BOLD
        if underlined:
            prefix += color.UNDERLINE

        return f"{prefix}{x}{color.END}"

def show_rc(parent, pos, n=1):
    parent_str = [str(m) for m in parent]
    parent_str[pos] = color.display(parent_str[pos], n=n, bold=True, underlined=True)
    print("|".join(parent_str))

def generate_forefathers(n_variants, n=10):
    forefathers = []
    for i in range(n):
        seq = [np.random.randint(ni) for ni in n_variants]
        forefathers.append(tuple(seq))
    return forefathers

def partition_population(individuals, n_partitions=2):
    # np.random.shuffle(individuals)
    return np.array_split(individuals, n_partitions)

def recombine_population(population, n_recombinations, return_rc=False):
    rcs = []
    for i in range(n_recombinations):
        print(f"======== Iteration {i:,} =======")
        (i1, i2) = np.random.choice(len(population), 2, replace=False)
        (p1, p2) = (population[i1], population[i2])
        pos = recombine(p1, p2)
        rcs.append((i1, i2, pos))

    if return_rc:
        return (population, rcs)
    return population

def recombine(s1, s2, show=True):
    pos = np.random.randint(len(s1))
    if show:
        for s in [s1, s2]:
            show_rc(s, pos)
        print("-"*30)
    (s2[pos], s1[pos]) = (s1[pos], s2[pos])
    if show:
        for s in [s1, s2]:
            show_rc(s, pos)
    return pos

def generate_module(min_size, max_size, n_variants):
    mean_size = np.random.randint(min_size, max_size)
    sizes = np.random.poisson(mean_size, size=n_variants)
    variants = ["".join(np.random.choice(NUCL, size=vi)) for vi in sizes]

    return tuple(variants)

if __name__ == '__main__':
    args = parse_args()

    n_variants = [np.random.randint(*args.num_variants_range)
                  for _ in range(args.n_modules)]
    forefathers = generate_forefathers(n_variants, n=args.n_forefathers)
    forefathers = partition_population(forefathers, args.n_subpopulations)
    children = [recombine_population(deepcopy(subpop), args.n_rc) for subpop in forefathers]

    forefathers = {f"F{i}.{k}": forefather
                   for k, subpop in enumerate(forefathers)
                   for i, forefather in enumerate(subpop)}
    children = {f"C{i}.{k}": child
                for k, subpop in enumerate(children)
                for i, child in enumerate(subpop)}
    
    modules = [generate_module(*args.module_size_range, ni) for ni in n_variants]
    junctions = [generate_module(*args.module_size_range, 1)[0] for _ in n_variants]

    with open(f"{args.sim_dir}/forefathers.fasta", "w") as writer:
        for (seq_id, forefather) in forefathers.items():
            meta = "-".join(map(str, forefather))
            seq = ""
            for pos, (variant, junction) in enumerate(zip(forefather, junctions)):
                seq += (junction + modules[pos][variant])
            writer.write(f">{seq_id} {meta}\n{seq}\n")

    with open(f"{args.sim_dir}/children.fasta", "w") as writer:
        for (seq_id, child) in children.items():
            meta = "-".join(map(str, child))
            seq = ""
            for pos, (variant, junction) in enumerate(zip(child, junctions)):
                seq += (junction + modules[pos][variant])
            writer.write(f">{seq_id} {meta}\n{seq}\n")
