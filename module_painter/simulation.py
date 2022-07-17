import logging
from functools import total_ordering
from itertools import combinations
from math import ceil
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
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

@total_ordering
class Phage:
    def __init__(self, variants, idx=None):
        assert len(variants) > 1
        self.idx = idx
        self.variants = variants
        self.history = set()
        self.abund = 1
        
    def __repr__(self):
        return f"{self.idx}. " + "|".join(map(str, self.variants))

    def __len__(self):
        return len(self.variants)

    def __hash__(self):
        return hash(tuple(self.variants))

    def __lt__(self, other):
        return tuple(self.variants) < tuple(other.variants)

    def __eq__(self, other):
        return tuple(self.variants) == tuple(other.variants)
    
    def as_rc(self, pos1, pos2):
        variants_str = [str(v) for v in self.variants]
        for pos in range(pos1, pos2):
            variants_str[pos] = color.display(variants_str[pos], n=1, bold=True, underlined=True)
        return "|".join(variants_str)

    @classmethod
    def random(cls, n_variants, idx=None):
        variants = [np.random.randint(0, ni) for ni in n_variants]
        phage = Phage(variants, idx=idx)
        phage.evolve()
        return phage

    @classmethod
    def from_recombination(cls, phages, positions, idx=None):
        np.random.shuffle(phages)

        variants = phages[0].variants.copy()
        for pos in range(*positions):
            variants[pos] = phages[1].variants[pos]

        phage = Phage(variants, idx=idx)
        rc_info = (phages[0].idx, phages[1].idx, positions[0], positions[1])
        phage.history = phages[0].history | phages[1].history
        phage.history.add(rc_info)

        return phage
    
    def evolve(self, lam=3):
        self.abund = 1 + np.random.poisson(lam)

    def as_fasta(self, modules, junctions):
        meta = "-".join(map(str, self.variants))
        seq = "".join(module[variant]+junction[0] for (module, variant, junction) in
                      zip(modules, self.variants, junctions))

        return f">{self.idx} {meta}\n{seq}\n"

class PhagePopulation:
    def __init__(self, *phages, rc_rate=0.5, name=""):
        self.name = name
        self.phages = set(phages)
        self.rc_rate = rc_rate
        self.rc_count = 0
        self.n_gen = 0

    def __repr__(self):
        return "\n".join(p.__repr__() for p in self.phages)

    def __len__(self):
        return len(self.phages)

    def __getitem__(self, key):
        try:
            return next(phage for phage in self.phages if phage.idx == key)
        except StopIteration:
            print(f"Phage {key} not found in population")
    
    def add(self, phage):
        phage.idx = f"{self.name}.{phage.idx}"
        self.phages.add(phage)
    
    def summarize(self):
        return dict(
            name=self.name, pop_size=len(self),
            n_rc=self.rc_count, rc_rate=self.rc_rate, n_gen=self.n_gen
        )

    @classmethod
    def random(self, lam, n_variants, rc_rate_range=None, name=""):
        n = max(2, np.random.poisson(lam))
        rc_rate = np.random.uniform(*rc_rate_range)
        
        population = PhagePopulation(name=name, rc_rate=rc_rate)
        for i in range(n):
            phage = Phage.random(n_variants, idx=i)
            population.add(phage)
        return population

    @classmethod
    def from_parents(cls, parents, n_generations, name=""):
        population = deepcopy(parents)
        population.name = name

        for i in range(n_generations):
            population.next_generation()

        population.phages = population.phages - parents.phages
        population.n_gen = n_generations

        return population
        
    def next_generation(self):
        n_rc = ceil(self.rc_rate*len(self))
        logger.debug(f"Next generation: {n_rc} recombinations")
        for i in range(n_rc):
            self.random_recombination()

    def random_recombination(self, verbose=False):
        indices = np.random.choice(len(self), 2, replace=False)
        phages = [phage for i, phage in enumerate(self.phages) if i in indices]
        pos = sorted(np.random.choice(range(len(phages[0])), 2))

        child = Phage.from_recombination(phages, pos, idx=len(self))

        if child in self.phages:
            return

        child.evolve() # not included in model for now
        self.add(child)

        if verbose:
            for phage in phages:
                logger.debug(phage.as_rc(*pos))
            logger.info("-"*30)
            logger.info(child.as_rc(*pos))

        self.rc_count += 1

    def shared_evolution(self, return_dict=False):
        shared = {}
        for (phage1, phage2) in combinations(self.phages, 2):
            shared[(phage1, phage2)] = len(
                phage1.history & phage2.history
            )

        if return_dict:
            return shared
        
        shared = pd.Series(shared)
        shared.index = shared.index.map(lambda p: (p[0].idx, p[1].idx))

        return shared

    def subpopulations(self):
        kinship = self.shared_evolution(return_dict=True)
        
        neighbors = {p: [] for p in self.phages}
        for (p1, p2), x in kinship.items():
            if x == 0:
                continue
            neighbors[p1].append(p2)
            neighbors[p2].append(p1)

        visited = set()
        
        for p in self.phages:
            if p in visited:
                continue
            component = []
            remaining = [p]
            while remaining:
                q = remaining.pop()
                if q in visited:
                    continue
                visited.add(q)
                remaining += neighbors[q]
                component.append(q)

            yield component
            
    def to_fasta(self, modules, junctions, output=None, mode="a"):
        with open(output, mode) as writer:
            for (i, phage) in enumerate(self.phages):
                entry = phage.as_fasta(modules, junctions)
                writer.write(entry)

def generate_module(min_size, max_size, n_variants):
    mean_size = np.random.randint(min_size, max_size)
    sizes = np.random.poisson(mean_size, size=n_variants)
    variants = ["".join(np.random.choice(NUCL, size=vi)) for vi in sizes]

    return tuple(variants)

def simulate(n_modules=None, n_variants_range=None, module_size_range=None,
             n_subpop=None, subpop_size_lam=None,
             n_gen_range=None, rate_range=None,
             outdir=None, **kwargs):

    n_variants = [np.random.randint(*n_variants_range) for _ in range(n_modules)]
    n_gen = np.random.randint(*n_gen_range)

    logger.info(f"#Modules: {n_modules}")
    logger.info("#Variant per module ~ U({}, {})".format(*n_variants_range))
    logger.info("Module size (bp) ~ U({}, {})".format(*module_size_range))
    logger.info(f"#Subpopulations: {n_subpop}")
    logger.info(f"#Subpopulation size lambda ~ Poisson({subpop_size_lam})")
    logger.info("#Generations: {} (~ U({}, {}))".format(n_gen, *n_gen_range))
    logger.info("Recombination rate ~ U({}, {})".format(*rate_range))    

    parents = [
        PhagePopulation.random(
            lam=subpop_size_lam,
            n_variants=n_variants,
            rc_rate_range=rate_range,
            name=f"P{k}"
        )
        for k in range(n_subpop)
    ]

    children = [PhagePopulation.from_parents(subpop, n_gen, name=f"C{k}")
                for (k, subpop) in enumerate(parents)]

    modules = [generate_module(*module_size_range, ni) for ni in n_variants]
    junctions = [generate_module(*module_size_range, 1) for _ in n_variants]
    
    # Write populations to fasta file
    for (name, data) in [("parents", parents), ("children", children)]:
        outfile = Path(outdir, f"{name}.fasta")
        outfile.unlink(missing_ok=True)
        for k, subpop in enumerate(data):
            subpop.to_fasta(modules, junctions, output=outfile)

    # Save some metadata
    clusters = {}
    k = -1
    for i, pop in enumerate(children):
        for subpop in pop.subpopulations():
            k += 1
            for p in subpop:
                clusters[p.idx] = k

    clusters = pd.Series(clusters, name="subpop_id")
    clusters.index.name = "seq_id"
    clusters.to_csv(f"{outdir}/subpopulations.csv")

    summary = pd.DataFrame([pop.summarize() for pop in children])
    summary["subpop_size_lambda"] = subpop_size_lam
    summary.to_csv(f"{outdir}/summary.csv", index=False)

    shared = pd.concat([
        subpop.shared_evolution()
        for subpop in children
    ])
    shared.to_csv(f"{outdir}/shared_evolution.csv", header=None)
