"""
Module for easily running and reporting Fishers exact tests.
"""

from scipy.stats import fisher_exact


def print_2x2_table(table, row_labels=['class 1', 'not'],
                    col_labels=['class 2', 'not'], fmt="%d"):
    """
    Prints a table used for Fisher's exact test. Adds row, column, and grand
    totals.

    Parameters
    ----------
    table : list or tuple, length=4
        The four cells of a 2x2 table: [r1c1, r1c2, r2c1, r2c2]

    row_labels : list
        A length-2 list of row names
    col_labels : list
        A length-2 list of column names
    """
    grand = sum(table)

    # Separate table into components and get row/col sums
    t11, t12, t21, t22 = table

    # Row sums, col sums, and grand total
    r1 = t11 + t12
    r2 = t21 + t22
    c1 = t11 + t21
    c2 = t12 + t22

    # Re-cast everything as the appropriate format
    t11, t12, t21, t22, c1, c2, r1, r2, grand = [
        fmt % i for i in [t11, t12, t21, t22, c1, c2, r1, r2, grand]]

    # Construct rows and columns the long way...
    rows = [
        [""] + col_labels + ['total'],
        [row_labels[0], t11, t12, r1],
        [row_labels[1], t21, t22, r2],
        ['total', c1, c2, grand],
    ]

    cols = [
        [row[0] for row in rows],
        [col_labels[0], t11, t21, c1],
        [col_labels[1], t12, t22, c2],
        ['total', r1, r2, grand],
    ]

    # Get max column width for each column; need this for nice justification
    widths = []
    for col in cols:
        widths.append(max(len(i) for i in col))

    # ReST-formatted header
    sep = ['=' * i for i in widths]

    # Construct the table one row at a time with nice justification
    s = []
    s.append(' '.join(sep))
    s.append(' '.join(i.ljust(j) for i, j in zip(rows[0], widths)))
    s.append(' '.join(sep))
    for row in rows[1:]:
        s.append(' '.join(i.ljust(j) for i, j in zip(row, widths)))
    s.append(' '.join(sep) + '\n')
    return '\n'.join([i.rstrip() for i in s])


def print_row_perc_table(table, row_labels=['class 1', 'not'],
                         col_labels=['class 2', 'not']):
    """
    Returns a string form of the table showing row percentages.

    Parameters
    ----------
    table : list or tuple, length=4
        The four cells of a 2x2 table: [r1c1, r1c2, r2c1, r2c2]

    row_labels : list
        A length-2 list of row names
    col_labels : list
        A length-2 list of column names
    """
    r1c1, r1c2, r2c1, r2c2 = map(float, table)
    row1 = r1c1 + r1c2
    row2 = r2c1 + r2c2

    blocks = [
        (r1c1, row1),
        (r1c2, row1),
        (r2c1, row2),
        (r2c2, row2)]

    new_table = []

    for cell, row in blocks:
        try:
            x = cell / row
        except ZeroDivisionError:
            x = 0
        new_table.append(x)

    s = print_2x2_table(new_table, row_labels, col_labels, fmt="%.2f")
    s = s.splitlines(True)
    del s[5]
    return '\n'.join([i.rstrip() for i in s])


def print_col_perc_table(table, row_labels=['class 1', 'not'],
                         col_labels=['class 2', 'not']):
    """
    Returns a string form of the table showing column percentages.

    Parameters
    ----------
    table : list or tuple, length=4
        The four cells of a 2x2 table: [r1c1, r1c2, r2c1, r2c2]

    row_labels : list
        A length-2 list of row names

    col_labels : list
        A length-2 list of column names
    """
    r1c1, r1c2, r2c1, r2c2 = map(float, table)
    col1 = r1c1 + r2c1
    col2 = r1c2 + r2c2

    blocks = [
        (r1c1, col1),
        (r1c2, col2),
        (r2c1, col1),
        (r2c2, col2)]

    new_table = []

    for cell, row in blocks:
        try:
            x = cell / row
        except ZeroDivisionError:
            x = 0
        new_table.append(x)

    s = print_2x2_table(new_table, row_labels, col_labels, fmt="%.2f")
    s = s.splitlines(False)
    last_space = s[0].rindex(" ")
    new_s = [i[:last_space] for i in s]
    return '\n'.join([i.rstrip() for i in new_s])


def fisher(table):
    """
    Perform a Fisher's exact test.

    Parameters
    ----------
    table : list or tuple, length=4
        The four cells of a 2x2 table: [r1c1, r1c2, r2c1, r2c2]

    Returns
    -------
    Tuple of (odds ratio, two-sided pvalue)

    """
    return fisher_exact(
        [
            [table[0], table[1]],
            [table[2], table[3]],
        ]
    )


def table_from_bool(ind1, ind2):
    """
    Given two boolean arrays, return the 2x2 contingency table

    ind1, ind2 : array-like
        Arrays of the same length
    """
    return [
            sum(ind1 & ind2),
            sum(ind1 & ~ind2),
            sum(~ind1 & ind2),
            sum(~ind1 & ~ind2),
        ]


def fisher_tables(table, row_names=['class 1', 'not'],
                  col_names=['class 2', 'not'], title=None):
    """
    Print all tables along with the FET results
    """
    s = []
    if title is not None:
        s.append(title)
        s.append('-' * len(title))
        s.append('')
    s.append(print_2x2_table(table, row_names, col_names))
    s.append('')
    s.append(print_row_perc_table(table, row_names, col_names))
    s.append('')
    s.append(print_col_perc_table(table, row_names, col_names))
    oddsratio, pval = fisher(table)
    s.append('odds ratio: {}'.format(oddsratio))
    s.append('2-sided pval: {}'.format(pval))
    return '\n'.join(s)
