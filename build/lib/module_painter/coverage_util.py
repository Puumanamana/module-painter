from module_painter.arc import Arc


def get_filler_arc(arc1, arc2):
    """
    Get filler arc between arc1 and arc2
    """
    # handle arc1 extending past slen
    if arc1.send >= arc1.slen:
        (arc1a, arc1b) = arc1.split_at_end()
        return get_filler_arc(arc1b, arc2)
    if arc1.send == arc1.slen-1: # filler starts at 0
        return Arc(0, arc2.sstart-1, arc1.slen, {"NA"})
    # after this point, arc1.send <= slen - 2
    # first, handle case where arc2 precedes arc1 (wrap around)
    if arc1.sstart > arc2.sstart: 
        # add filler after arc1
        return Arc(arc1.send+1, arc1.slen+arc2.sstart-1, arc1.slen, {"NA"})
    # General case: arc1 precedes arc2
    return Arc(arc1.send+1, arc2.sstart-1, arc1.slen, {"NA"})

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
