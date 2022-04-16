import pandas as pd

from modular_painter.rules import (sync_boundaries, is_embedded, circularize, merge_if_close)



def test_sync_start():
    df = pd.DataFrame(dict(
        sstart=[0, 1, 10, 20, 22, 26, 30]
    ))
    df = sync_boundaries(df, "sstart", 3)
    assert df.sstart.to_list() == [0, 0, 10, 20, 20, 26, 30]

def test_sync_end():
    df = pd.DataFrame(dict(
        send=[0, 1, 10, 20, 22, 26, 30]
    ))
    df = sync_boundaries(df, "send", 3)
    assert df.send.to_list() == [1, 1, 10, 22, 22, 26, 30]
    
