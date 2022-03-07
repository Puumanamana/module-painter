import igraph
import numpy as np

def compute_arc_dist(i1, i2, size):
    """
    assymmetric distance between two intervals.
    """
    (x1, y1) = i1
    (x2, y2) = i2
    
    # case 1: both intervals wrap around
    if y1 >= size and y2 >= size:
        return 0
    # case 2: i1 wraps around
    if x2 < x1:
        return max(0, x2-(y1-size))
    # general case
    return max(0, x2-y1)

class OverlapGraph(igraph.Graph):
    def __init__(self, arcs, size=None, target=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.target = target
        self.genome_size = size

        for arc in arcs:
            self.add_arc(arc)

        self.compute_all_overlaps()

    def __len__(self):
        return len(self.vs.indices)
        
    def add_arc(self, arc):
        self.add_vertex(
            start=arc.start, end=arc.end,
            meta=arc.meta
        )

    def compute_all_overlaps(self):
        """
        Go iterate over all vertices in sorted order
        - 2 embedded loops to compute all overlaps
        - inner loop stops when there is no overlap with next
        """
        indices = np.array(self.vs.indices)
        boundaries = np.array([
            self.vs["start"],
            self.vs["end"]
        ]).T
        
        order = np.argsort(boundaries[:, 0] + 1.j*boundaries[:, 1])
        indices = indices[order]
        boundaries = boundaries[order, :]

        for i, itv1 in enumerate(boundaries):
            overlap = True
            k = i

            while overlap:
                k = (k+1) % len(self)
                itv2 = boundaries[k, :]
                dist = compute_arc_dist(itv1, itv2, self.genome_size)
                overlap = (dist <= 0) and (k != i)

                if overlap:
                    self.add_edge(indices[i], indices[k])

    def __repr__(self, subset=None):
        display = [f"Target: {self.target} (L={self.genome_size})"]
        for i, v in enumerate(self.vs):
            if subset is not None:
                if not v.meta.intersection(subset):
                    continue
            display.append(f"{i:3} - {v['start']:>9}{v['end']:>9}: {v['meta']}")
        return "\n".join(display)
                
