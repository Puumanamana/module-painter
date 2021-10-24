import pandas as pd
from modular_painter.arc import Arc
from modular_painter.coverage import Coverage
from modular_painter.coverage_util import get_successor


def test_len():
    arcs = [Arc(0, 12, 100, {"A"}), Arc(0, 10, 100, {"B"})]
    coverage = Coverage(*arcs)
    assert len(coverage) == len(arcs)

def test_to_pandas():
    pass

def test_from_coverages():
    pass

def test_sort():
    cov = Coverage(
        Arc(0, 12, 100, {"A"}),
        Arc(0, 10, 100, {"B"}),
        Arc(8, 30, 100, {"C"}),
        Arc(11, 60, 100, {"D"}),
        Arc(60, 111, 100, {"E"}),
        Arc(65, 80, 100, {"F"})
    )
    cov.sort()
    names = "".join(par for arc in cov.arcs for par in arc.meta)
    assert names == "BACDEF"

def test_get_all_parents():
    cov = Coverage(
        Arc(0, 1, 10, {"A", "B"}),
        Arc(1, 2, 10, {"C", "B"}),
        Arc(0, 5, 10, {"D", "E"})        
    )
    assert cov.get_all_parents() == set("ABCDE")

def test_is_covered():
    cov = Coverage(
        Arc(0, 12, 100),
        Arc(11, 102, 100)
    )
    assert cov.is_covered()
 
def test_not_covered():
    cov = Coverage(
        Arc(0, 12, 100),
        Arc(11, 80, 100)
    )
    assert not cov.is_covered()

def test_sync_boundaries():
    cov = Coverage(
        Arc(0, 10, 100),
        Arc(1, 12, 100),
        Arc(2, 5, 100),
        Arc(12, 45, 100),
        Arc(11, 40, 100),        
        Arc(16, 40, 100),        
    )

    cov.sync_boundaries(3, which="start")

    assert cov.arcs[1].start == 0
    assert cov.arcs[2].start == 0
    assert cov.arcs[3].start == 11

def test_extend_arcs_per_parent():
    cov = Coverage(
        Arc(0, 2, 20, {"A"}),
        Arc(1, 4, 20, {"A"}),
        Arc(3, 7, 20, {"A"}),
        Arc(0, 5, 20, {"B"}),
        Arc(4, 7, 20, {"B"}),
        Arc(0, 2, 20, {"C"}),
        Arc(7, 10, 20, {"C"}),
    )
    cov.sort()
    cov.extend_arcs_per_parent(2)

    assert sum(1 for arc in cov.arcs if arc.flagged) == 3
    assert any(arc.bounds() == (0, 7) for arc in cov.arcs if arc.meta == {"A"})
    assert any(arc.bounds() == (0, 7) for arc in cov.arcs if arc.meta == {"B"})

def test_merge_equal_intervals():
    cov = Coverage(
        Arc(0, 3, 10, {"A"}),
        Arc(0, 3, 10, {"B"}),
        Arc(0, 5, 10, {"C"})
    )
    cov.merge_equal_intervals()
    assert cov.arcs[0].flagged and not cov.arcs[1].flagged
    assert cov.arcs[1].meta == {"A", "B"}

def test_fill_gaps():
    cov = Coverage(
        Arc(0, 2, 20, {"A"}),
        Arc(4, 6, 20, {"B"}),
        Arc(13, 18, 20, {"C"}),
    )
    cov.sort()
    cov.fill_gaps(2)

    na_arc = [arc for arc in cov.arcs if arc.meta == "NA"]
    assert len(na_arc) == 1
    assert na_arc[0].bounds() == (6, 13)
    assert cov.arcs[0].bounds() == (0, 3)
    assert cov.arcs[1].bounds() == (3, 6)
    assert cov.arcs[-1].bounds() == (13, 20)

def test_simplify_embedded():
    cov = Coverage(
        Arc(0, 3, 10, {"A"}),
        Arc(1, 2, 10, {"B"}),
        Arc(5, 10, 10, {"C"}),
    )
    cov.sort()
    cov.simplify_embedded()
    assert cov.arcs[1].flagged

def test_simplify_embedded_loop():
    cov = Coverage(
        Arc(0, 2, 10, {"A"}),
        Arc(1, 4, 10, {"B"}),
        Arc(1, 10, 10, {"C"}),
    )
    cov.sort()
    cov.simplify_embedded()
    assert cov.arcs[1].flagged

def test_get_successors():
    cov = Coverage(
        Arc(0, 3, 10),
        Arc(1, 4, 10),
        Arc(2, 5, 10),
        Arc(4, 10, 10)
    )
    successors = [get_successor(i, cov.arcs) for (i, arc) in enumerate(cov.arcs)]

    assert successors == [2, 3, 3, 0]

def test_lee_and_lee_case_1():
    cov = Coverage(*
        (Arc(*boundaries, 100) for boundaries in [
            (0,22),(21,40),(39,60),(59,80),(79,110),
            (9,30),(29,55),(54,70),(69,90)
        ])
    )
    cov.sort()
    opt_cov = cov.get_minimal_coverage()

    assert len(opt_cov) == 5
    assert opt_cov.is_covered()
    assert opt_cov.arcs == [
        Arc(*boundaries, 100) for boundaries in
        [(0,22), (21,40), (39,60), (59,80), (79,110)]
    ]
    
def test_lee_and_lee_case_2():
    cov = Coverage(*
        (Arc(*boundaries, 100) for boundaries in [
            (0,20),(19,40),(39,60),(59,80),(79,111),
            (10,25),(24,50),(49,70),(69,110)
        ])
    )
    cov.sort()
    opt_cov = cov.get_minimal_coverage()

    assert len(opt_cov) == 4
    assert opt_cov.is_covered()
    assert opt_cov.arcs == [
        Arc(*boundaries, 100) for boundaries in
        [(10,25), (24,50), (49,70), (69,110)]
    ]

def test_lee_and_lee_case_3():
    cov = Coverage(*
        (Arc(*boundaries, 100) for boundaries in [
            (0,22),(21,40),(39,60),(59,80),(79,110),
            (9,30),(29,55),(54,70),(69,105),(4,27),
            (27,48)
        ])
    )
    cov.sort()
    opt_cov = cov.get_minimal_coverage()

    assert len(opt_cov) == 5
    assert opt_cov.is_covered()
    assert opt_cov.arcs == [
        Arc(*boundaries, 100) for boundaries in
        [(0,22), (21,40), (39,60), (59,80), (79,110)]
    ]
    
