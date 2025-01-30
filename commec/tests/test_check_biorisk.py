import pandas as pd
import pytest


from commec.screeners.check_biorisk import _get_nt_len_from_query

def test_get_nt_len_from_query():
    df = pd.DataFrame({
        'query name': ['query1', 'query2', 'query1'],
        'qlen': [100, 200, 197]
    })
    
    class MockQueryResult:
        def __init__(self, length):
            self.query_length = length
            
    class MockScreenResult:
        def get_query(self, name):
            return MockQueryResult(150) if name == 'query1' else None
    
    lengths = _get_nt_len_from_query(df, MockScreenResult())
    assert lengths.tolist() == [150, 600, 150]