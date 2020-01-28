import numpy as np
import pandas as pd

from modular_painter.painting import clean_results

def test_clean():
    painting = pd.DataFrame(np.array(
        [[3, 9],
         [3, 9],
         [0, 5],
         [2, 7],
         [2, 12],
         [2, 12],
         [2, 12],
         [3, 12],         
         [5, 8],
         [7, 8],
         [7, 9]]), columns=['s2','e2'])
    painting.sort_values(by=['s2','e2'], inplace=True)

    assert np.sum(clean_results(painting).values != np.array([[0, 5], [2, 12], [2, 12], [2, 12]])) == 0

if __name__ == '__main__':
    test_clean()
