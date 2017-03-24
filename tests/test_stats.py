from textwrap import dedent
import numpy as np
from lcdblib.stats import fisher

table = [12, 5, 29, 2]


def fix(s):
    # Remove common leading whitespace and blank lines at top or bottom
    ss = [i.rstrip() for i in s.splitlines(False)]
    if len(ss[0]) == 0:
        ss = ss[1:]
    if len(ss[-1]) == 0:
        ss = ss[:-1]
    return dedent("\n".join(ss))


def test_2x2():
    s = fisher.print_2x2_table(
        table,
        row_labels=['Selected', 'Not selected'],
        col_labels=['Having the property', 'Not having the property']
    )
    assert s == fix("""
    ============ =================== ======================= =====
                 Having the property Not having the property total
    ============ =================== ======================= =====
    Selected     12                  5                       17
    Not selected 29                  2                       31
    total        41                  7                       48
    ============ =================== ======================= =====
    """)


def test_table_from_bool():
    ind1 = np.ones(sum(table)) == 0
    ind2 = np.ones(sum(table)) == 0
    ind1[:table[0]] = True
    ind2[:table[0]] = True
    ind1[table[0]:table[0] + table[1]] = True
    ind2[table[0] + table[1]:table[0] + table[1] + table[2]] = True
    assert fisher.table_from_bool(ind1, ind2) == table


def test_fisher_tables():
    s = fisher.fisher_tables(table, row_names=['has property', 'does not'],
                      col_names=['selected', 'not'], title='testing')
    assert s == fix("""
    testing
    -------

    ============ ======== === =====
                 selected not total
    ============ ======== === =====
    has property 12       5   17
    does not     29       2   31
    total        41       7   48
    ============ ======== === =====

    ============ ======== ==== =====
                 selected not  total
    ============ ======== ==== =====
    has property 0.71     0.29 1.00
    does not     0.94     0.06 1.00
    ============ ======== ==== =====

    ============ ======== ====
                 selected not
    ============ ======== ====
    has property 0.29     0.71
    does not     0.71     0.29
    total        1.00     1.00
    ============ ======== ====
    odds ratio: 0.165517
    2-sided pval: 0.0802686
    """)
