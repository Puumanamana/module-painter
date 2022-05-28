from module_painter.arc import Arc


def get_filler_arc(arc1, arc2):
    """
    Get filler arc between arc1 and arc2
    """
    # handle arc1 extending past size
    if arc1.end >= arc1.size:
        (arc1a, arc1b) = arc1.split_at_end()
        return get_filler_arc(arc1b, arc2)
    if arc1.end == arc1.size-1: # filler starts at 0
        return Arc(0, arc2.start-1, arc1.size, {"NA"})
    # after this point, arc1.end <= size - 2
    # first, handle case where arc2 precedes arc1 (wrap around)
    if arc1.start > arc2.start: 
        # add filler after arc1
        return Arc(arc1.end+1, arc1.size+arc2.start-1, arc1.size, {"NA"})
    # General case: arc1 precedes arc2
    return Arc(arc1.end+1, arc2.start-1, arc1.size, {"NA"})

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
