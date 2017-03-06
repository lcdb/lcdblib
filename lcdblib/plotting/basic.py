#!/usr/bin/env python
""" Set of basic plotting functions. """
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from seaborn import utils


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
    ax.text(0.1, .9, "r={:0.4}".format(corr(x, y)[0]), transform=ax.transAxes,
            **kwargs)


def maPlot(x, y, data=None, title=None, log=False, **kwargs):
    """ Creates a MA plot.

    Parameters
    ----------
    x, y: string, series, or vector array
        Input variables. If strings, these should correspond with column names
        in ``data``. When pandas objects are used, axes will be labeled with
        the series name.
    data : DataFrame
        Tidy ("long-form") dataframe where each column is a variable and each
        row is an observation.
    title: str
        Title to add to the plot
    log: bool
        If true the log2 of the Geometric mean and ratio will be used instead
        of the mean and difference.
    **kwargs: dict
        Values to pass to seaborn.regplot

    Returns
    -------
    matplotlib.axes.Axes: A matplotlib axes.

    Example
    -------
    Given that 'x' and 'y' are column names in df:
        >>> maPlot('x', 'y', data=df, log=True, fig_reg=True,
        ... line_kws={'color': 'red'})

    Given that 'x' and 'y' are pandas.Series
        >>> maPlot(x, y, log=True, fig_reg=True, line_kws={'color': 'red'})


    """
    if isinstance(x, str) and isinstance(y, str) and data is not None:
        dat = data[[x, y]].copy()
        xlabel = x
        ylabel = y
    elif isinstance(x, pd.Series) and isinstance(y, pd.Series):
        dat = pd.DataFrame(pd.concat([x, y], axis=1))
        xlabel = x.name
        ylabel = y.name
    else:
        raise ValueError("x and y must be a column name in data or "
                         "pandas.Series")

    dat.reset_index(drop=True, inplace=True)

    if log:
        dat = dat[~(dat == 0).any(axis=1)]

        # Calculate differences and mean
        dat['difference'] = np.log2(dat.iloc[:, 0] / dat.iloc[:, 1])
        dat['average'] = np.log2(np.sqrt(dat.iloc[:, 0] * dat.iloc[:, 1]))
        diff = 'log2({x} / {y})'.format(x=xlabel, y=ylabel)
        avg = 'log2(sqrt({x} * {y}))'.format(x=xlabel, y=ylabel)
    else:
        dat['difference'] = dat.iloc[:, 0] - dat.iloc[:, 1]
        dat['average'] = (dat.iloc[:, 0] + dat.iloc[:, 1]) / 2
        diff = '{x} - {y}'.format(x=xlabel, y=ylabel)
        avg = '({x} + {y}) / 2'.format(x=xlabel, y=ylabel)

    # Plot
    if 'fit_reg' not in kwargs:
        kwargs['fit_reg'] = False
    if 'ci' not in kwargs:
        kwargs['ci'] = False

    ax = sns.regplot('average', 'difference', data=dat, **kwargs)
    ax.set_xlabel(avg)
    ax.set_ylabel(diff)
    ax.axhline(0, ls='--', color='r')

    if title is not None:
        ax.set_title(title)

    return ax


class PairGrid(sns.PairGrid):
    def __init__(self, data, hue=None, hue_order=None, palette=None,
                 hue_kws=None, vars=None, x_vars=None, y_vars=None,
                 diag_sharey=True, size=2.5, aspect=1,
                 despine=True, dropna=True, subplots_kws={}):
        """A slight modification of seaborn.PairGrid.

        PairGrid required that the axes of the plots be shared. I contacted the
        author about allowing the axes to be separate, but he did not like that
        idea. I have subclassed seaborn.PairGrid and just modified the
        __init__() to allow separate axes.

        Initialize the plot figure and PairGrid object.

        Parameters
        ----------
        data : DataFrame
            Tidy (long-form) dataframe where each column is a variable and
            each row is an observation.
        hue : string (variable name), optional
            Variable in ``data`` to map plot aspects to different colors.
        hue_order : list of strings
            Order for the levels of the hue variable in the palette
        palette : dict or seaborn color palette
            Set of colors for mapping the ``hue`` variable. If a dict, keys
            should be values  in the ``hue`` variable.
        hue_kws : dictionary of param -> list of values mapping
            Other keyword arguments to insert into the plotting call to let
            other plot attributes vary across levels of the hue variable (e.g.
            the markers in a scatterplot).
        vars : list of variable names, optional
            Variables within ``data`` to use, otherwise use every column with
            a numeric datatype.
        {x, y}_vars : lists of variable names, optional
            Variables within ``data`` to use separately for the rows and
            columns of the figure; i.e. to make a non-square plot.
        size : scalar, optional
            Height (in inches) of each facet.
        aspect : scalar, optional
            Aspect * size gives the width (in inches) of each facet.
        despine : boolean, optional
            Remove the top and right spines from the plots.
        dropna : boolean, optional
            Drop missing values from the data before plotting.

        See Also
        --------
        pairplot : Easily drawing common uses of :class:`PairGrid`.
        FacetGrid : Subplot grid for plotting conditional relationships.

        Examples
        --------

        Draw a scatterplot for each pairwise relationship:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> import seaborn as sns; sns.set()
            >>> iris = sns.load_dataset("iris")
            >>> g = sns.PairGrid(iris)
            >>> g = g.map(plt.scatter)

        Show a univariate distribution on the diagonal:

        .. plot::
            :context: close-figs

            >>> g = sns.PairGrid(iris)
            >>> g = g.map_diag(plt.hist)
            >>> g = g.map_offdiag(plt.scatter)

        (It's not actually necessary to catch the return value every time,
        as it is the same object, but it makes it easier to deal with the
        doctests).

        Color the points using a categorical variable:

        .. plot::
            :context: close-figs

            >>> g = sns.PairGrid(iris, hue="species")
            >>> g = g.map_diag(plt.hist)
            >>> g = g.map_offdiag(plt.scatter)
            >>> g = g.add_legend()

        Use a different style to show multiple histograms:

        .. plot::
            :context: close-figs

            >>> g = sns.PairGrid(iris, hue="species")
            >>> g = g.map_diag(plt.hist, histtype="step", linewidth=3)
            >>> g = g.map_offdiag(plt.scatter)
            >>> g = g.add_legend()

        Plot a subset of variables

        .. plot::
            :context: close-figs

            >>> g = sns.PairGrid(iris, vars=["sepal_length", "sepal_width"])
            >>> g = g.map(plt.scatter)

        Pass additional keyword arguments to the functions

        .. plot::
            :context: close-figs

            >>> g = sns.PairGrid(iris)
            >>> g = g.map_diag(plt.hist, edgecolor="w")
            >>> g = g.map_offdiag(plt.scatter, edgecolor="w", s=40)

        Use different variables for the rows and columns:

        .. plot::
            :context: close-figs

            >>> g = sns.PairGrid(iris,
            ...                  x_vars=["sepal_length", "sepal_width"],
            ...                  y_vars=["petal_length", "petal_width"])
            >>> g = g.map(plt.scatter)

        Use different functions on the upper and lower triangles:

        .. plot::
            :context: close-figs

            >>> g = sns.PairGrid(iris)
            >>> g = g.map_upper(plt.scatter)
            >>> g = g.map_lower(sns.kdeplot, cmap="Blues_d")
            >>> g = g.map_diag(sns.kdeplot, lw=3, legend=False)

        Use different colors and markers for each categorical level:

        .. plot::
            :context: close-figs

            >>> g = sns.PairGrid(iris, hue="species", palette="Set2",
            ...                  hue_kws={"marker": ["o", "s", "D"]})
            >>> g = g.map(plt.scatter, linewidths=1, edgecolor="w", s=40)
            >>> g = g.add_legend()

        """

        # Sort out the variables that define the grid
        if vars is not None:
            x_vars = list(vars)
            y_vars = list(vars)
        elif (x_vars is not None) or (y_vars is not None):
            if (x_vars is None) or (y_vars is None):
                raise ValueError("Must specify `x_vars` and `y_vars`")
        else:
            numeric_cols = self._find_numeric_cols(data)
            x_vars = numeric_cols
            y_vars = numeric_cols

        if np.isscalar(x_vars):
            x_vars = [x_vars]
        if np.isscalar(y_vars):
            y_vars = [y_vars]

        self.x_vars = list(x_vars)
        self.y_vars = list(y_vars)
        self.square_grid = self.x_vars == self.y_vars

        # Create the figure and the array of subplots
        figsize = len(x_vars) * size * aspect, len(y_vars) * size

        kwargs = {'figsize': figsize, 'sharex': "col",
                  'sharey': "row", 'squeeze': False}
        kwargs.update(subplots_kws)
        fig, axes = plt.subplots(len(y_vars), len(x_vars), **kwargs)

        self.fig = fig
        self.axes = axes
        self.data = data

        # Save what we are going to do with the diagonal
        self.diag_sharey = diag_sharey
        self.diag_axes = None

        # Label the axes
        self._add_axis_labels()

        # Sort out the hue variable
        self._hue_var = hue
        if hue is None:
            self.hue_names = ["_nolegend_"]
            self.hue_vals = pd.Series(["_nolegend_"] * len(data),
                                      index=data.index)
        else:
            hue_names = utils.categorical_order(data[hue], hue_order)
            if dropna:
                # Filter NA from the list of unique hue names
                hue_names = list(filter(pd.notnull, hue_names))
            self.hue_names = hue_names
            self.hue_vals = data[hue]

        # Additional dict of kwarg -> list of values for mapping the hue var
        self.hue_kws = hue_kws if hue_kws is not None else {}

        self.palette = self._get_palette(data, hue, hue_order, palette)
        self._legend_data = {}

        # Make the plot look nice
        if despine:
            utils.despine(fig=fig)
        fig.tight_layout()


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
