import pandas as pd

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
