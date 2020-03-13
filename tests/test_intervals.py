import pandas as pd
from modular_painter.intervals import Arc, Coverage

def arc_data(x=1, y=20, n=2, labels=None):
    if labels is None:
        labels = ['V{}'.format(i) for i in range(n)]
    elif isinstance(labels, str):
        labels = [labels]

    intervals = list(zip(labels, [x]*n, [y]*n, [1]*n))

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
    arc = Arc.from_pandas(arc_data(1, 20, 1), 100, 'target')
    assert arc.mod == 100
    assert arc.target == 'target'
    assert len(arc) == 20

def test_arc_high_low():
    arc = Arc.from_pandas(arc_data(90, 20), 100, 'target')
    assert arc.start == 90
    assert arc.end == 120
    assert len(arc) == 31

def test_arc_intersect():
    arc1 = Arc.from_pandas(arc_data(1, 10), 100, 'target')
    arc2 = Arc.from_pandas(arc_data(5, 20, 2), 100, 'target')
    arc3 = Arc.from_pandas(arc_data(90, 4), 100, 'target')
    arc4 = Arc.from_pandas(arc_data(12, 20, 3), 100, 'target')

    assert arc1.intersect(arc2)
    assert arc1.intersect(arc3)
    assert arc3.intersect(arc1)
    assert arc3.intersect(arc2)
    assert arc2.intersect(arc3)
    assert not arc1.intersect(arc4)

def test_coverage_init():
    arc1 = Arc.from_pandas(arc_data(1, 10, 1), 100, 'target')
    arc2 = Arc.from_pandas(arc_data(9, 80), 100, 'target')

    cov = Coverage(arc1, arc2)

    assert not cov.dirty
    assert cov.mod == arc1.mod
    assert cov.target == arc2.target

def test_coverage_from_pandas():
    data = coverage_data()

    cov1 = Coverage.from_pandas(data, 'TARGET', modulo=100)
    cov2 = Coverage.from_pandas(data, 'TARGET')

    pos = cov1.get('start', 'end')

    assert cov1.mod == 100 and cov2.mod == 98
    assert pos.end[0] == 50
    assert pos.values[5, 0] == 90 and pos.values[5, 1] == 130
    assert cov2.target == 'TARGET'

def test_coverage_successor():
    cov = Coverage.from_pandas(coverage_data(), 'target', modulo=100)
    cov.set_all_successors()

    assert cov.successors[0] == 3
    assert cov.successors[5] == 1
    assert cov.successors[1] == 4

def test_merge_close_intervals():
    cov = Coverage(
        Arc.from_pandas(arc_data(2, 9), 100, 'target'),
        Arc.from_pandas(arc_data(3, 10), 100, 'target'),
        Arc.from_pandas(arc_data(1, 12), 100, 'target'),        
    )

    cov.extend_close_intervals(1)
    cov.merge_equal_intervals()

    assert set(cov[0].bounds()) == {1, 12}
    assert set(cov[1].bounds()) == {1, 10}

def test_coverage_fill_or_extend():
    cov = Coverage(
        Arc.from_pandas(arc_data(1, 10, 1, 'V2'), 100, 'target'),
        Arc.from_pandas(arc_data(30, 50, 1, 'V1'), 100, 'target'),       
        Arc.from_pandas(arc_data(53, 80, 1, 'V2'), 100, 'target'),
        Arc.from_pandas(arc_data(72, 80, 2, ['V4', 'V2']), 100, 'target'),
        Arc.from_pandas(arc_data(83, 85, 1, 'V3'), 100, 'target'),
        Arc.from_pandas(arc_data(96, 99), 100, 'target'),
        Arc.from_pandas(arc_data(96, 98), 100, 'target'),
    )

    print(cov)
    cov.fill_or_extend(10)

    assert set(cov[0].bounds()) == {1, 10}
    assert set(cov[1].bounds()) == {11, 29}
    assert sum(sum(x.data.index.str.startswith('NoCov')) for x in cov) == 2
    assert set(cov[3].bounds()) == {53, 82}
    assert set(cov[-1].bounds()) == {96, 100}

def test_coverage_is_covered():
    cov = Coverage(
        Arc.from_pandas(arc_data(1, 10), 100, 'target'),
        Arc.from_pandas(arc_data(30, 90), 100, 'target'),
        Arc.from_pandas(arc_data(94, 99), 100, 'target'),
    )

    not_cov = cov.is_covered()
    
    cov.fill_or_extend(20)
    
    assert not not_cov
    assert cov.is_covered()

def test_lee_and_lee_case_1():
    cov = Coverage(
        Arc.from_pandas(arc_data(0, 22), 100, 'target'),
        Arc.from_pandas(arc_data(21, 40), 100, 'target'),        
        Arc.from_pandas(arc_data(39, 60), 100, 'target'),
        Arc.from_pandas(arc_data(59, 80), 100, 'target'),
        Arc.from_pandas(arc_data(79, 10), 100, 'target'),        
        Arc.from_pandas(arc_data(9, 30), 100, 'target'),
        Arc.from_pandas(arc_data(29, 55), 100, 'target'),
        Arc.from_pandas(arc_data(54, 70), 100, 'target'),
        Arc.from_pandas(arc_data(69, 90), 100, 'target'),
    )
    
    opt_cov = cov.get_minimal_coverage()

    print(opt_cov)
    assert len(opt_cov) == 5
    assert opt_cov.is_covered()


def test_lee_and_lee_case_2():
    cov = Coverage(
        Arc.from_pandas(arc_data(0, 20), 100, 'target'),
        Arc.from_pandas(arc_data(19, 40), 100, 'target'),        
        Arc.from_pandas(arc_data(39, 60), 100, 'target'),
        Arc.from_pandas(arc_data(59, 80), 100, 'target'),
        Arc.from_pandas(arc_data(79, 11), 100, 'target'),        
        Arc.from_pandas(arc_data(10, 25), 100, 'target'),
        Arc.from_pandas(arc_data(24, 50), 100, 'target'),
        Arc.from_pandas(arc_data(49, 70), 100, 'target'),
        Arc.from_pandas(arc_data(69, 10), 100, 'target'),
    )

    opt_cov = cov.get_minimal_coverage()

    print(opt_cov)
    assert len(opt_cov) == 4
    assert opt_cov.is_covered()

def test_lee_and_lee_case_3():
    cov = Coverage(
        Arc.from_pandas(arc_data(0, 22), 100, 'target'),
        Arc.from_pandas(arc_data(21, 40), 100, 'target'),        
        Arc.from_pandas(arc_data(39, 60), 100, 'target'),
        Arc.from_pandas(arc_data(59, 80), 100, 'target'),
        Arc.from_pandas(arc_data(79, 10), 100, 'target'),        
        Arc.from_pandas(arc_data(9, 30), 100, 'target'),
        Arc.from_pandas(arc_data(29, 55), 100, 'target'),
        Arc.from_pandas(arc_data(54, 70), 100, 'target'),
        Arc.from_pandas(arc_data(69, 5), 100, 'target'),
        Arc.from_pandas(arc_data(4, 27), 100, 'target'),
        Arc.from_pandas(arc_data(27, 48), 100, 'target'),
    )

    
    opt_cov = cov.get_minimal_coverage()

    print(opt_cov)
    assert len(opt_cov) == 5
    assert opt_cov.is_covered()

    
if __name__ == '__main__':

    # test_arc_intersect()
    # test_coverage_is_covered()
    # test_coverage_from_pandas()
    # test_coverage_successor()
    test_merge_close_intervals()
    # test_coverage_fill_or_extend()
    # test_lee_and_lee_case_1()
