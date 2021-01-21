# Summary of the method

## Inputs

Two fasta files measuring bacteriophages populations in two different conditions. There must be common species between them. One fasta is the source population (S) and the other one is the target population (T). The goal is to "paint" the target population with the source population.

## Identification of bacterial species

1. We blast S on T 
2. For each genome in T, we retrieve all hits with min\_identity=0.8 and min\_coverage=0.5

Note: We merge supplementary alignments if the distance between them is lower than min\_module\_size=40bp

We store the hits (potential parents of the genomes in S) in a dictionary.

## Painting

This step find the optimal coverage/painting of each genome in S with their parents.

1. Extend close intervals: extend the start/end of an interval if it's less than 100bp away from another one.
2. Merge intervals with the same boundaries
3. Fill missing coverage with putative genomes or extend existing coverage up to module_size=40bp
4. Circularize genome
5. Extract minimal coverage (Lee and Lee algorithm)


