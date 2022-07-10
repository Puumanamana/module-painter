import logging
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
        logger.debug(f"======== Iteration {i:,} =======")
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
        logger.debug("-"*30)
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

def simulate(num_variants_range=None, n_modules=None, n_forefathers=None, n_rc=None,
             n_subpopulations=None, module_size_range=None, outdir=None, **kwargs):

    logger.info(f"#Forefathers: {n_forefathers}")
    logger.info(f"#Subpopulations: {n_subpopulations}")
    logger.info(f"#Variants per module: {'-'.join(map(str, num_variants_range))}")
    logger.info(f"Module size range: {'-'.join(map(str, module_size_range))}")    
    
    n_variants = [np.random.randint(*num_variants_range)
                  for _ in range(n_modules)]
    forefathers = generate_forefathers(n_variants, n=n_forefathers)
    forefathers = partition_population(forefathers, n_subpopulations)
    children = [recombine_population(deepcopy(subpop), n_rc) for subpop in forefathers]

    forefathers = {f"F{k}.{i}": forefather
                   for k, subpop in enumerate(forefathers)
                   for i, forefather in enumerate(subpop)}

    children = {f"C{k}.{i}": child
                for k, subpop in enumerate(children)
                for i, child in enumerate(subpop)}

    modules = [generate_module(*module_size_range, ni) for ni in n_variants]
    junctions = [generate_module(*module_size_range, 1)[0] for _ in n_variants]

    with open(f"{outdir}/forefathers.fasta", "w") as writer:
        for (seq_id, forefather) in forefathers.items():
            meta = "-".join(map(str, forefather))
            seq = ""
            for pos, (variant, junction) in enumerate(zip(forefather, junctions)):
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
    logger.debug("====Forefathers====")
    cmap = {seq_id: i for (i, seq_id) in enumerate(forefathers)}

    for (seq_id, forefather) in forefathers.items():
        f_str =  f"{seq_id}: " + " | ".join(map(str, forefather))
        logger.debug(color.display(f_str, n=cmap[seq_id]))

    logger.debug("====Children====")
    for (seq_id, child) in children.items():
        cluster = seq_id.split(".")[0][1:]
        forefather_cluster = [it for it in forefathers.items()
                              if it[0][1:].startswith(cluster)]
        c_str = []
        for j, variant in enumerate(child):
            parent = next(it[0] for it in forefather_cluster if it[1][j] == variant)
            c_str.append(color.display(variant, n=cmap[parent]))
        logger.debug(f"{seq_id}: " + " | ".join(c_str))
            
