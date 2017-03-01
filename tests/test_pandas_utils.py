import pytest
import pandas as pd
from lcdblib.pandas import utils

@pytest.fixture(scope='session')
def sample_table():
    metadata = {
        'sample': ['one', 'two'],
        'tissue': ['ovary', 'testis']
        }
    return pd.DataFrame(metadata)

def test_cartesian_df(sample_table):
    df2 = pd.DataFrame({'num': [100, 200]})
    result = utils.cartesian_product(sample_table, df2)

    # Compare a slice
    sf = result.iloc[0, :].sort_index()
    sf.name = ''
    test_sf = pd.Series({'sample': 'one', 'tissue': 'ovary', 'num': 100}, name='').sort_index()
    assert sf.equals(test_sf)
    assert result.shape == (4, 3)

def test_cartesian_sf(sample_table):
    sf2 = pd.Series([100, 200], name='num')
    result = utils.cartesian_product(sample_table, sf2)

    # Compare a slice
    sf = result.iloc[0, :].sort_index()
    sf.name = ''
    test_sf = pd.Series({'sample': 'one', 'tissue': 'ovary', 'num': 100}, name='').sort_index()
    assert sf.equals(test_sf)
    assert result.shape == (4, 3)

def test_cartesian_dict(sample_table):
    df2 = {'num': [100, 200]}
    result = utils.cartesian_product(sample_table, df2)

    # Compare a slice
    sf = result.iloc[0, :].sort_index()
    sf.name = ''
    test_sf = pd.Series({'sample': 'one', 'tissue': 'ovary', 'num': 100}, name='').sort_index()
    assert sf.equals(test_sf)
    assert result.shape == (4, 3)
