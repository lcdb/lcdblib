#!/usr/bin/env python
""" Set of basic plotting functions. """
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def corrfunc(x, y, type='spearman', **kwargs):
    """ Adds text to the current axes with either the sparman or pearson r.

    Parameters
    ----------
    x: array-like
        Vector of values.

    y: array-like
        Vector of values.

    type: str
        Which type of correlation coefficient to use. Either "spearman" or
        "pearson". Uses scipy.stats module.

    """
    from scipy.stats import spearmanr, pearsonr

    if type == 'spearman':
        corr = spearmanr
    elif type == 'pearson':
        corr = pearsonr

    ax = plt.gca()
    ax.text(0.1, .9, "r={:0.4}".format(corr(x, y)[0]), transform=ax.transAxes, **kwargs)


def maPlot(df, x, y, **kwargs):
    """ Creates a MA plot.

    Parameters
    ----------
    df: pandas.DataFrame
        DataFrame of values to compare.
        `
    x: str
        First column name for comparison.

    y: str
        Second column name for comparison.

    """

    x_dat = df[x]
    y_dat = df[y]
    difference = np.array(np.log2(x_dat / y_dat))
    average = np.array(np.log10((x_dat + y_dat)/2))

    # Convert nan to 0
    difference = np.nan_to_num(difference)
    average = np.nan_to_num(average)
    sns.regplot(average, difference, fig_reg=True, ci=False, **kwargs)
    ax.get_gca()
    ax.set_xlabel('log10(({x} + {y}) / 2)'.format(x=x, y=y))
    ax.set_ylabel('log2({x} / {y})'.format(x=x, y=y))
    ax.axhline(0, ls='--', color='r')

    
def lowerTriangle(df, func, func_kw={}, pairgrid_kw={}, **kwargs):
    """ Create a PairGrid lower triangle panel. 

    Parameters
    ----------
    df: pandas.DataFrame
        DataFrame containing the columns you want to plot in a PairGrid.

    func: function or list of functions
        The function you want to plot in a PairGrid
    
    Returns
    -------
    seaborn.PairGrid

    """
    p = sns.PairGrid(df, **pairgrid_kw)
    if isinstance(func, list):
        for f in zip(func, func_kw):
            kw = func_kw.get(str(f), {})
            p.map_lower(f, **kw)
    else:
        p.map_lower(func, **func_kw)

    for i, j in zip(*np.triu_indices_from(p.axes, 0)):
        p.axes[i, j].set_visible(False)

    return p

if __name__ == '__main__':
    pass
