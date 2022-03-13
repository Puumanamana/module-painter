import igraph
import numpy as np


class OverlapGraph(igraph.Graph):
    def __init__(self, coverage, size=None, target=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.target = target
        self.genome_size = size

        for arc in arcs:
            self.add_arc(arc)

        self.compute_all_overlaps()

    def __repr__(self):
        display = [f"Target: {self.target} (L={self.genome_size})"]
        for i, v in enumerate(self.vs):
            display.append(f"{i:3} - {v['start']:>9}{v['end']:>9}: {v['parent']}")
        return "\n".join(display)

    @classmethod
    def from_coverage(cls, cov):
        g = igraph.Graph()

        # add vertices
        for i, arc in enumerate(cov.iter_arcs(wrap=False)):
            for meta in arc.meta:
                uid = f"{cov.ref}.{meta}.{arc.start}-{arc.end}"
                g.add_vertex(name=uid, parent=meta, order=i, ref=cov.ref,
                             start=arc.start, end=arc.end)

        # add edges
        arcs = cov.iter_arcs(wrap=True)
        prev_arc = next(arcs)

        for arc in arcs:
            for meta1 in prev_arc.meta:
                uid1 = f"{cov.ref}.{meta1}.{prev_arc.start}-{prev_arc.end}"
                start = arc.start
                end = prev_arc.end % cov.size # wrapping condition. This should be the only arc with end >= size

                for meta2 in arc.meta:
                    uid2 = f"{cov.ref}.{meta2}.{arc.start}-{arc.end}"
                    parents = "/".join(sorted([meta1, meta2]))
                    g.add_edge(uid1, uid2, ref=cov.ref, start=start, end=end, parents=parents)
            prev_arc = arc

        return g
        
