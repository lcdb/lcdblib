import pandas as pd
from itertools import product

def tidy_dataframe(df, column, sep='|'):
    """
    Given a dataframe with a delimiter-separated string, create a new dataframe
    with separate rows for each value in each corresponding string.

    E.g.::

               gene  score
        0     g1|g2      1
        1        g3      5
        2  g4|g5|g6      9

    becomes:

        tidy_dataframe(df, 'gene', sep='|')
        #   gene  score
        # 0   g1      1
        # 0   g2      1
        # 1   g3      5
        # 2   g4      9
        # 2   g5      9
        # 2   g6      9
    """
    s = df[column].str.split(sep, expand=True).stack()
    i = s.index.get_level_values(0)
    df2 = df.loc[i].copy()
    df2[column] = s.values
    return df2

def cartesian_product(df1, df2):
    """ Calculates the carteisan product and returns expanded DataFrame.

    Given a pandas.DataFrame:

    | sample | tissue |
    |--------|--------|
    | one    | ovary  |
    | two    | testis |

    and some set of values `{'num': [100, 200]}` build:

    | sample | tissue | num |
    |--------|--------|-----|
    | one    | ovary  | 100 |
    | one    | ovary  | 200 |
    | two    | testis | 100 |
    | two    | testis | 200 |

    Parameters
    ----------
    df1: pandas.DataFrame
        A DataFrame that you want to expand.
    df2: dict of array-like | pandas.DataFrame | pandas.Series
        The set of values that you want to expand df1 by.

    """
    if isinstance(df2, dict):
        df2 = pd.DataFrame(df2)
    elif isinstance(df2, pd.Series):
        df2 = df2.to_frame()

    rows = product(df1.iterrows(), df2.iterrows())
    df = pd.DataFrame(left.append(right) for (_, left), (_, right) in rows)

    return df.reset_index(drop=True)
