import pandas as pd
from modular_painter.intervals import Arc, Coverage

def arc_data(x=1, y=20, n=2):
    intervals = [
        ['V{}'.format(i), x, y, 1] for i in range(n)
    ]

    data = pd.DataFrame(intervals, columns=['source', 'start', 'end', 'identity'])
    return data

def coverage_data():
    intervals = [
        ['V1', 0, 50, 40],
        ['V2', 40, 60, 15],
        ['V3', 40, 60, 18],        
        ['V1', 55, 98, 20],
        ['V3', 90, 30, 20],
        ['V4', 90, 30, 25],        
        ['V4', 20, 54, 30],
        ['V3', 50, 95, 40]
    ]
     
    data = pd.DataFrame(intervals, columns=['source', 'start', 'end', 'identity'])

    return data

def test_arc_low_high():
    arc = Arc(arc_data(1, 20, 1), 'target', 100)
    assert arc.mod == 100
    assert arc.target == 'target'
    assert len(arc) == 20

def test_arc_high_low():
    arc = Arc(arc_data(90, 20), 'target', 100)
    assert arc.start == 90
    assert arc.end == 120
    assert len(arc) == 31

def test_arc_intersect():
    arc1 = Arc(arc_data(1, 10), 'target', 100)
    arc2 = Arc(arc_data(5, 20, 2), 'target', 100)
    arc3 = Arc(arc_data(90, 4), 'target', 100)

    assert arc1.intersect(arc2)
    assert arc1.intersect(arc3)
    assert arc3.intersect(arc1)
    assert not arc3.intersect(arc2)
    assert not arc2.intersect(arc3)

def test_coverage_init():
    arc1 = Arc(arc_data(1, 10, 1), 'target', 100)
    arc2 = Arc(arc_data(9, 80), 'target', 100)

    cov = Coverage(arc1, arc2)

    assert not cov.dirty
    assert cov.mod == arc1.mod
    assert cov.target == arc2.target

def test_coverage_from_pandas():
    data = coverage_data()

    cov1 = Coverage.from_pandas(data, 'TARGET', modulo=100)
    cov2 = Coverage.from_pandas(data, 'TARGET')

    pos = cov1.get_intervals()

    assert cov1.mod == 100 and cov2.mod == 98
    assert pos[0, 1] == 50
    assert pos[5, 0] == 90 and pos[5, 1] == 130
    assert cov2.target == 'TARGET'

def test_coverage_insert():
    cov = Coverage.from_pandas(coverage_data(), 'target', modulo=100)
    new = Arc(arc_data(25, 55), 'target', 100)

    l1 = len(cov)

    cov.insert(new)

    assert cov.dirty
    assert len(cov.data) == l1 + 1

def test_coverage_successor():
    cov = Coverage.from_pandas(coverage_data(), 'target', modulo=100)
    cov.set_all_successors()

    assert cov.successors[0] == 3
    assert cov.successors[5] == 1
    assert cov.successors[1] == 3

def test_coverage_fill():
    cov = Coverage(
        Arc(arc_data(1, 10), 'target', 100),
        Arc(arc_data(30, 50), 'target', 100),        
        Arc(arc_data(53, 70), 'target', 100),
        Arc(arc_data(94, 99), 'target', 100)
    )
    cov.fill_missing_data(max_merge_dist=3)

    assert set(cov[1].bounds()) == {10, 30}
    assert set(cov[4].bounds()) == {70, 94}
    assert sum(sum(x.data.source=='NoCov') for x in cov) == 2
    assert cov[2].end == 53
    assert cov[5].end == 101

def test_coverage_is_covered():
    cov = Coverage(
        Arc(arc_data(1, 10), 'target', 100),
        Arc(arc_data(30, 90), 'target', 100),
        Arc(arc_data(94, 99), 'target', 100),
    )

    not_cov = cov.is_covered()
    
    cov.fill_missing_data(3)
    
    assert not not_cov
    assert cov.is_covered()

def test_merge_close_intervals():
    cov = Coverage(
        Arc(arc_data(2, 9), 'target', 100),
        Arc(arc_data(3, 10), 'target', 100),
        Arc(arc_data(1, 12), 'target', 100),        
    )

    cov.merge_close_intervals(3)
    assert cov[1].end == 10    
    assert cov[2].start == 2
    assert set(cov[0].bounds()) == {1, 12}
    
def test_simplify():
    cov = Coverage(
        Arc(arc_data(1, 10), 'target', 100),
        Arc(arc_data(3, 10), 'target', 100),
        Arc(arc_data(3, 20), 'target', 100),        
        Arc(arc_data(30, 90), 'target', 100),
        Arc(arc_data(40, 70), 'target', 100),
        Arc(arc_data(95, 9), 'target', 100),        
        Arc(arc_data(99, 15), 'target', 100),
    )

    cov.simplify()
    assert len(cov) == 4
    
def test_lee_and_lee_case_1():
    cov = Coverage(
        Arc(arc_data(0, 22), 'target', 100),
        Arc(arc_data(21, 40), 'target', 100),        
        Arc(arc_data(39, 60), 'target', 100),
        Arc(arc_data(59, 80), 'target', 100),
        Arc(arc_data(79, 10), 'target', 100),        
        Arc(arc_data(9, 30), 'target', 100),
        Arc(arc_data(29, 55), 'target', 100),
        Arc(arc_data(54, 70), 'target', 100),
        Arc(arc_data(69, 90), 'target', 100),
    )
    
    opt_cov = cov.get_minimal_coverage()

    print(opt_cov)
    assert len(opt_cov) == 5
    assert opt_cov.is_covered()


def test_lee_and_lee_case_2():
    cov = Coverage(
        Arc(arc_data(0, 20), 'target', 100),
        Arc(arc_data(19, 40), 'target', 100),        
        Arc(arc_data(39, 60), 'target', 100),
        Arc(arc_data(59, 80), 'target', 100),
        Arc(arc_data(79, 11), 'target', 100),        
        Arc(arc_data(10, 25), 'target', 100),
        Arc(arc_data(24, 50), 'target', 100),
        Arc(arc_data(49, 70), 'target', 100),
        Arc(arc_data(69, 10), 'target', 100),
    )

    opt_cov = cov.get_minimal_coverage()

    print(opt_cov)
    assert len(opt_cov) == 4
    assert opt_cov.is_covered()

def test_lee_and_lee_case_3():
    cov = Coverage(
        Arc(arc_data(0, 22), 'target', 100),
        Arc(arc_data(21, 40), 'target', 100),        
        Arc(arc_data(39, 60), 'target', 100),
        Arc(arc_data(59, 80), 'target', 100),
        Arc(arc_data(79, 10), 'target', 100),        
        Arc(arc_data(9, 30), 'target', 100),
        Arc(arc_data(29, 55), 'target', 100),
        Arc(arc_data(54, 70), 'target', 100),
        Arc(arc_data(69, 5), 'target', 100),
        Arc(arc_data(4, 28), 'target', 100),
        Arc(arc_data(27, 48), 'target', 100),
    )
    
    opt_cov = cov.get_minimal_coverage()

    print(opt_cov)
    assert len(opt_cov) == 5
    assert opt_cov.is_covered()

    
if __name__ == '__main__':

    # test_merge_close_intervals()
    # test_coverage_is_covered()
    
    # test_coverage_from_pandas()
    # test_coverage_successor()
    # test_coverage_fill()
    test_simplify()
    # test_coverage_insert()
    # test_lee_and_lee_case_1()
