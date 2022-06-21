from module_painter.arc import Arc

def test_len():
    arc = Arc(90, 110, 100)
    assert len(arc) == 21

def test_eq():
    arc1 = Arc(90, 110, 100, {"A", "C"})
    arc2 = Arc(90, 110, 100, {"C", "A"})

    assert arc1 == arc2

def test_neq():
    arc1 = Arc(90, 110, 100, {"A", "C"})
    arc2 = Arc(90, 110, 100, {"B"})

    assert arc1 != arc2

def test_flag():
    arc = Arc(90, 110, 100, {"A", "C"})
    arc.flag()
    assert arc.flagged

def test_flag():
    arc = Arc(90, 110, 100, {"A", "C"})
    arc.flag()
    arc.unflag()
    assert not arc.flagged

def test_bounds():
    arc = Arc(90, 110, 100, {"A"})

    assert arc.bounds() == (90, 110)
    
def test_split_at_end():
    arc = Arc(90, 110, 100)
    arc1, arc2 = arc.split_at_end()

    assert (arc1.sstart, arc1.send) == (90, 99)
    assert (arc2.sstart, arc2.send) == (0, 10)

def test_is_embedded_no_wrap():
    arc1 = Arc(10, 110, 100)
    arc2 = Arc(30, 90, 100)

    assert arc2.is_embedded(arc1)
    assert not arc1.is_embedded(arc2)

def test_is_embedded_wrap():
    arc1 = Arc(10, 130, 100)
    arc2 = Arc(0, 20, 100)

    assert arc2.is_embedded(arc1)
    assert not arc1.is_embedded(arc2)

def test_is_embedded_wrap_full():
    arc1 = Arc(29, 130, 100)
    arc2 = Arc(10, 20, 100)

    assert arc2.is_embedded(arc1)
    
def test_pos_dist_to_next():
    arc1 = Arc(10, 20, 50)
    arc2 = Arc(30, 40, 50)

    assert arc1.dist_to_next(arc2) == 9
    assert arc2.dist_to_next(arc1) == 19

def test_overlap_dist_to_next():
    arc1 = Arc(10, 20, 50)
    arc2 = Arc(15, 40, 50)

    assert arc1.dist_to_next(arc2) == -6

def test_self_wrap_dist_to_next():
    arc1 = Arc(45, 55, 50)
    arc2 = Arc(10, 40, 50)

    assert arc1.dist_to_next(arc2) == 4

def test_other_wrap_dist_to_next():
    arc1 = Arc(10, 40, 50)
    arc2 = Arc(45, 55, 50)

    assert arc1.dist_to_next(arc2) == 4

def test_both_wrap_dist_to_next():
    arc1 = Arc(30, 55, 50)
    arc2 = Arc(45, 60, 50)

    assert arc1.dist_to_next(arc2) == -11

def test_both_wrap_overlap_dist_to_next():
    arc1 = Arc(30, 55, 50)
    arc2 = Arc(45, 85, 50) # 45 -> 35

    assert arc1.dist_to_next(arc2) == -11
    
def test_try_merge_equal():
    arc1 = Arc(10, 20, 50, {"A"})
    arc2 = Arc(10, 20, 50, {"B"})
    arc1.try_merge_with(arc2)

    assert arc1.flagged and not arc2.flagged
    assert arc2.qacc == {"A", "B"}

def test_try_merge_different():
    arc1 = Arc(10, 20, 50)
    arc2 = Arc(15, 20, 50)
    arc1.try_merge_with(arc2)

    assert not arc1.flagged and not arc2.flagged

def test_try_fuse_overlap():
    arc1 = Arc(10, 20, 50, {"A"})
    arc2 = Arc(15, 30, 50, {"A", "B"})
    arc1.try_fuse_with(arc2, max_dist=5, skip_nw=True)
    
    assert arc1.flagged and not arc2.flagged
    assert arc2.qacc == {"A"}
    assert arc2.bounds() == (10, 30)

def test_try_fuse_close():
    arc1 = Arc(10, 20, 50, {"A"})
    arc2 = Arc(26, 30, 50, {"A", "B"})
    arc1.try_fuse_with(arc2, max_dist=10, skip_nw=True)
    
    assert arc1.flagged and not arc2.flagged
    assert arc2.qacc == {"A"}
    assert arc2.bounds() == (10, 30)
    
def test_try_fuse_wrap():
    arc1 = Arc(30, 55, 50, {"A"})
    arc2 = Arc(7, 15, 50, {"A", "B"})
    arc1.try_fuse_with(arc2, max_dist=10, skip_nw=True)
    
    assert arc1.flagged and not arc2.flagged
    assert arc2.qacc == {"A"}
    assert arc2.bounds() == (30, 65)
    
def test_try_fuse_far():
    arc1 = Arc(10, 20, 50, {"A"})
    arc2 = Arc(26, 30, 50, {"A", "B"})
    arc1.try_fuse_with(arc2, max_dist=2, skip_nw=True)
    
    assert not arc1.flagged and not arc2.flagged

def test_try_extend_close():
    arc1 = Arc(10, 20, 50, {"A"})
    arc2 = Arc(26, 30, 50, {"B"})
    arc1.try_extend_end_with(arc2, max_dist=10)
    
    assert not arc1.flagged and not arc2.flagged
    assert arc1.bounds() == (10, 23)
    assert arc2.bounds() == (24, 30)
    
def test_try_extend_close_and_loop():
    arc1 = Arc(30, 55, 50, {"A"})
    arc2 = Arc(7, 15, 50, {"B"})
    arc1.try_extend_end_with(arc2, max_dist=10)
    
    assert arc1.bounds() == (30, 55)
    assert arc2.bounds() == (6, 15)

def test_try_extend_far():
    arc1 = Arc(40, 55, 50, {"A"})
    arc2 = Arc(25, 30, 50, {"B"})
    arc1.try_extend_end_with(arc2, max_dist=5)
    
    assert arc1.bounds() == (40, 55)
    assert arc2.bounds() == (25, 30)
