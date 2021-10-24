from modular_painter.arc import Arc

def get_filler_arc(arc1, arc2, dist):
    """
    Get filler arc between arc1 and arc2
    """
    size = arc1.size

    # Rare case (when looping): the filler fits at the beginning
    if arc2.start < arc1.start and arc2.start >= dist:
        return Arc(arc2.start-dist, arc2.start, size, meta="NA")
    # General case
    return Arc(arc1.end, arc1.end+dist, size, meta="NA")

def get_successor(i, arcs):
    """
    Get position of successor of {i] in {arcs}
    """
    j = i
    while True:
        j = (j+1) % len(arcs)

        if arcs[i].dist_to_next(arcs[j]) > 0:
            break
        if j == i:
            break
    return (j-1) % len(arcs)
