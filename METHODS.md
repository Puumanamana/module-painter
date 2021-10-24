# Summary of the method

## Inputs

Two fasta files measuring bacteriophages populations in two different conditions. There must be common species between them. One fasta is the source population (S) and the other one is the target population (T). The goal is to "paint" the target population with the source population.

## Identification of bacterial species

1. We blast S on T 
2. For each genome in T, we retrieve all hits with `min_identity=0.8` and `min_coverage=0.5`

Note: We merge supplementary alignments if the distance between them is lower than `min_module_size`. This also makes coverage circular if the parent spans both the start and end of the genome.

We store the hits (potential parents of the genomes in S) in a dictionary.

## Painting simplification

## Criteria
i) Make the intervals comparable between parents. In other words, 2 different intervals that have "almost" the same boundaries should be considered as equal
ii) There should not be any gaps in coverage. If there is, we either:
  - add a putative parent to fill the gap if it's large
  - extend both parents to touch if it's small
iii) There should not be any interval embedded in another one. If this happens, we keep the larger one.
	
Once these conditions are met, we can use the algorithm from Lee & Lee to find the optimal coverage
The order in which we satisfy those condition matter, so we need to be careful when we do it.

Finally, the circularity of theses interval needs to be handled carefully at each step since the reasoning is a bit different for those.

### Approach

This step finds the optimal coverage/painting of each genome in S with their parents.

1. Clustering of intervals from fuzzy interval boundaries:
   a) We cluster the start of all intervals using Agglomerative clustering:
	   - Each interval starts in its own cluster, and clusters are successively merged together
	   - Complete linkage: We minimize the maximum distance between intervals of pairs of clusters
	   - Stop criteria: When any two clusters c1, c2 cannot be merged without having `dist(c1, c2) <  min_module_size`
   b) Repeat step a with the ends of each intervals
   
2. Extension of close parent intervals. For each parent:
   - Interate throught the sorted intervals
   - Merge (interval1, interval2) if `dist(interval1,interval2)<=min_module_size`

3. Equal interval simplification
   - Merge any interval (from different parents) if they have the exact same boundaries

4. Missing coverage: If there is a gap G in coverage between intervals I and J, then:
   - if `G > min_module_size`, we add a putative parent
   - if `0 < G <= min_module_size`: we extend J to touch I.
   We could do it at the same time as 2), but doing it this way simplifies the approach when iterating though consecutive arcs.
   
5. Remove embedded intervals and keep the largest one
   
6. Extract minimal coverage (Lee and Lee algorithm)

## Breakpoints extraction

This steps extracts and maps breakpoints between phages

1. Blast each breakpoint against the others
2. Build a graph where the vertices are the breakpoints and we draw an edge if the blast between one breakpoints covers the other with 80% coverage and 80% identity and is at least 10bp long

## Breakpoint clustering

1. Parent selection: If multiple choice of parents are possible for a given breakpoints, we use the following procedure to pick one:
   While there are breakpoints with multiple parents:
       1. Find the parent pair with the most recombinations in the population
       2. Solve equalities by selecting the parent pair that covers the most phages
       3. Whenever this parent pair has multiple choices in a breakpoints, remove all other choices
2. Clustering: Leiden agglomerative clustering. We build a graph with phages as vertices. We add an edge between 2 phages is they have a breakpoint in common. 




