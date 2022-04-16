from copy import deepcopy

import numpy as np
from Bio import SeqIO

np.random.seed(123)

NUCL = list("ACGT")
(MIN_MODULE_SIZE, MAX_MODULE_SIZE) = (50, 1000)
(MIN_VARIANTS, MAX_VARIANTS) = (2, 5)


def generate_forefathers(n_variants, n=10):
    forefathers = []
    for i in range(n):
        seq = [np.random.randint(ni) for ni in n_variants]
        forefathers.append(seq)
    return forefathers

def generate_children(population, n_recombinations):
    for _ in range(n_recombinations):
        (i1, i2) = np.random.choice(len(population), 2, replace=False)
        recombine(population[i1], population[i2])
    return population

def recombine(s1, s2):
    pos = np.random.randint(len(s1))
    (s2[pos], s1[pos]) = (s1[pos], s2[pos])

if __name__ == '__main__':
    outdir = "/tmp/cedric/modular_painting_tests/sim"
    poisson_l = 300
    n_modules = 10
    n_forefathers = 10
    n_recombinations = 30

    n_variants = [np.random.randint(MIN_VARIANTS, MAX_VARIANTS) for _ in range(n_modules)]
    forefathers = generate_forefathers(n_variants, n=n_forefathers)
    children = generate_children(deepcopy(forefathers), n_recombinations)

    modules = []
    junctions = []
    for i, ni in enumerate(n_variants):
        mi = np.random.randint(MIN_MODULE_SIZE, MAX_MODULE_SIZE)
        variant_sizes = np.random.poisson(mi, size=ni)
        variants = ["".join(np.random.choice(NUCL, size=vi)) for vi in variant_sizes]
        modules.append(variants)

        junction_length = np.random.randint(MIN_MODULE_SIZE, MAX_MODULE_SIZE)
        junction = "".join(np.random.choice(NUCL, size=junction_length))
        junctions.append(junction)

    with open(f"{outdir}/forefathers.fasta", "w") as writer:
        for i, variants in enumerate(forefathers):
            seqid = chr(65+i) + " " + "-".join(map(str, variants))
            seq = ""
            for pos, (variant, junction) in enumerate(zip(variants, junctions)):
                seq += (junction + modules[pos][variant])
            writer.write(f">{seqid}\n{seq}\n")

    with open(f"{outdir}/children.fasta", "w") as writer:
        for i, variants in enumerate(children):
            seqid = chr(97+i) + " " + "-".join(map(str, variants))
            seq = ""
            for pos, (variant, junction) in enumerate(zip(variants, junctions)):
                seq += (junction + modules[pos][variant])
            writer.write(f">{seqid}\n{seq}\n")
            
