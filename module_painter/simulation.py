import logging
from math import ceil
import random
from copy import deepcopy
from pathlib import Path

import numpy as np
from Bio import SeqIO


logger = logging.getLogger("module-painter")

NUCL = list("ACGT")

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
    logger.debug("|".join(parent_str))

def generate_parents(n_variants, n=10):
    parents = []
    for i in range(n):
        seq = [random.randint(0, ni-1) for ni in n_variants]
        parents.append(tuple(seq))
    return parents

def partition_population(individuals, n_partitions=2):
    # np.random.shuffle(individuals)
    return np.array_split(individuals, n_partitions)

def expand(population, lam=3):
    n = np.random.poisson(lam, size=len(population))
    return [seq for (ni, seq) in zip(n, population) for _ in range(ni)]

def recombine_population(population, rate_range=(0.2, 1), n_gen=1, return_rc=True):
    rate = random.uniform(*rate_range)
    n_recombinations = ceil(rate*len(population)*n_gen)

    rcs = []
    for i in range(n_recombinations):
        logger.debug(f"======== Iteration {i:,} =======")
        (i1, i2) = np.random.choice(len(population), 2, replace=False)
        (p1, p2) = (population[i1], population[i2])
        pos = recombine(p1, p2)
        rcs.append((i1, i2, pos))

    if return_rc:
        return (population, rcs)
    return population

def recombine(s1, s2, show=True):
    pos = random.randint(0, len(s1)-1)
    if show:
        for s in [s1, s2]:
            show_rc(s, pos)
        logger.debug("-"*30)
    (s2[pos], s1[pos]) = (s1[pos], s2[pos])
    if show:
        for s in [s1, s2]:
            show_rc(s, pos)
    return pos

def generate_module(min_size, max_size, n_variants):
    mean_size = random.randint(min_size, max_size)
    sizes = np.random.poisson(mean_size, size=n_variants)
    variants = ["".join(np.random.choice(NUCL, size=vi)) for vi in sizes]

    return tuple(variants)

def simulate(n_modules=None, n_parents=None,
             n_subpopulations=None, module_size_range=None, outdir=None,
             n_variants_range=None, n_gen_range=None, rate_range=None,
             **kwargs):

    n_gen = random.randint(*n_gen_range)
    n_variants = [random.randint(*n_variants_range) for _ in range(n_modules)]

    logger.info(f"#Parents: {n_parents}")
    logger.info(f"#Subpopulations: {n_subpopulations}")
    logger.info(f"#Variants per module: {','.join(map(str, n_variants))} (range: {'-'.join(map(str, n_variants_range))})")
    logger.info(f"Module size range: {'-'.join(map(str, module_size_range))}")    
    logger.info(f"#Generations: {n_gen} (range: {'-'.join(map(str, n_gen_range))})")    
    logger.info(f"Recombination rate range: {'-'.join(map(str, rate_range))}")    
    
    parents = generate_parents(n_variants, n=n_parents)
    parents = partition_population(parents, n_subpopulations)
    parents = [expand(subpop) for subpop in parents]

    children = []
    n_rc = []
    for subpop in parents:
        (recombined, rc) = recombine_population(deepcopy(subpop), rate_range=rate_range, n_gen=n_gen)
        children.append(recombined)
        n_rc.append(len(rc))

    pop_sizes = [len(subpop) for subpop in parents]
    with open(f"{outdir}/summary.csv", "w") as handle:
        handle.write(",".join(["N", "n_gen", "n_rc"]) + "\n")
        vals = map(str, [sum(pop_sizes), n_gen, sum(n_rc)])
        handle.write(",".join(vals) + "\n")
    
    parents = {f"P{k}.{i}": parent
                   for k, subpop in enumerate(parents)
                   for i, parent in enumerate(subpop)}

    children = {f"C{k}.{i}": child
                for k, subpop in enumerate(children)
                for i, child in enumerate(subpop)}
    
    modules = [generate_module(*module_size_range, ni) for ni in n_variants]
    junctions = [generate_module(*module_size_range, 1)[0] for _ in n_variants]

    with open(f"{outdir}/parents.fasta", "w") as writer:
        for (seq_id, parent) in parents.items():
            meta = "-".join(map(str, parent))
            seq = ""
            for pos, (variant, junction) in enumerate(zip(parent, junctions)):
                seq += (junction + modules[pos][variant])
            writer.write(f">{seq_id} {meta}\n{seq}\n")

    with open(f"{outdir}/children.fasta", "w") as writer:
        for (seq_id, child) in children.items():
            meta = "-".join(map(str, child))
            seq = ""
            for pos, (variant, junction) in enumerate(zip(child, junctions)):
                seq += (junction + modules[pos][variant])
            writer.write(f">{seq_id} {meta}\n{seq}\n")

    # Some more logging info
    logger.debug("====Parents====")
    cmap = {seq_id: i for (i, seq_id) in enumerate(parents)}

    for (seq_id, parent) in parents.items():
        f_str =  f"{seq_id}: " + " | ".join(map(str, parent))
        logger.debug(color.display(f_str, n=cmap[seq_id]))

    logger.debug("====Children====")
    for (seq_id, child) in children.items():
        cluster = seq_id.split(".")[0][1:]
        parent_cluster = [it for it in parents.items()
                              if it[0][1:].startswith(cluster)]
        c_str = []
        for j, variant in enumerate(child):
            parent = next(it[0] for it in parent_cluster if it[1][j] == variant)
            c_str.append(color.display(variant, n=cmap[parent]))
        logger.debug(f"{seq_id}: " + " | ".join(c_str))
            
