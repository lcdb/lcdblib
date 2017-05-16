from textwrap import dedent
import numpy as np
import pandas
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.transforms import blended_transform_factory
from matplotlib.ticker import MaxNLocator
from matplotlib.collections import EventCollection
import gffutils
import pybedtools
from pybedtools import featurefuncs

from .. utils import utils
from . import colormap_adjust

_base_doc = """%s
The underlying pandas.DataFrame is always available with the `data`
attribute.

Any attributes not explicitly in this class will be looked for in the
underlying pandas.DataFrame.

Parameters
----------
data : string or pandas.DataFrame
    If string, assumes it's a filename and calls
    pandas.read_table(data, **import_kwargs).

db : string or gffutils.FeatureDB
    Optional database that can be used to generate features

import_kwargs : dict
    These arguments will be passed to pandas.read_table() if `data` is
    a filename.
"""


class ResultsTable(object):
    __doc__ = _base_doc % dedent(
        """
        Wrapper around a pandas.DataFrame that adds additional functionality.
        """)

    def __init__(self, data, db=None, import_kwargs=None):
        if isinstance(data, str):
            import_kwargs = import_kwargs or {}
            data = pandas.read_table(data, **import_kwargs)
        if not isinstance(data, pandas.DataFrame):
            raise ValueError("`data` is not a pandas.DataFrame")
        self.data = data

        self._kwargs = dict(db=db, import_kwargs=import_kwargs)
        self.attach_db(db)
        self._cached_features = None

    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self.data, attr)

    def __getitem__(self, attr):
        if isinstance(attr, str):
            return self.data.__getitem__(attr)
        else:
            return self.__class__(self.data.__getitem__(attr), **self._kwargs)

    def update(self, dataframe):
        """
        Updates the current data with a new dataframe.

        This extra step is required to get around the fancy pandas.DataFrame
        indexing (like .ix, .iloc, etc).
        """
        return self.__class__(dataframe, **self._kwargs)

    def copy(self):
        data = self.data.copy(deep=True)
        return self.__class__(data, db=self.db, import_kwargs=self._kwargs)

    def __repr__(self):
        s = []
        s.append("<%s instance, wrapping the following:"
                 % self.__class__.__name__)
        s.append('')
        s.extend('\t' + i for i in repr(self.data).splitlines(False))
        s.append('>')
        return '\n'.join(s)

    def attach_db(self, db):
        """
        Attach a gffutils.FeatureDB for access to features.

        Useful if you want to attach a db after this instance has already been
        created.

        Parameters
        ----------
        db : gffutils.FeatureDB
        """
        if db is not None:
            if isinstance(db, str):
                db = gffutils.FeatureDB(db)
            if not isinstance(db, gffutils.FeatureDB):
                raise ValueError(
                    "`db` must be a filename or a gffutils.FeatureDB")
        self._kwargs['db'] = db
        self.db = db

    def features(self, ignore_unknown=False):
        """
        Generator of features.

        If a gffutils.FeatureDB is attached, returns a pybedtools.Interval for
        every feature in the dataframe's index.

        Parameters
        ----------
        ignore_unknown : bool
            If True, silently ignores features that are not found in the db.
        """
        if not self.db:
            raise ValueError("Please attach a gffutils.FeatureDB")
        for i in self.data.index:
            try:
                yield gffutils.helpers.asinterval(self.db[i])
            except gffutils.FeatureNotFoundError:
                if ignore_unknown:
                    continue
                else:
                    raise gffutils.FeatureNotFoundError('%s not found' % i)

    def reindex_to(self, x, attribute=None):
        """
        Returns a copy that only has rows corresponding to feature names in x.


        Parameters
        ----------
        x : str or pybedtools.BedTool
            If str, then assume it's a filename. BED, GFF, GTF, or VCF where
            the "Name" field (that is, the value returned by feature['Name'])
            or any arbitrary attribute

        attribute : str or int or None
            If `x` is GFF or GTF format, and `attribute` is str, then attribute
            containing the name of the feature to use. If `x` format is BED and
            `attribute` is str, then use getattr on the interval (e.g., 'name'
            or 'score'). If `attribute` is int, then use that column. If None,
            then use the "name" attribute of the Interval, which falls back to
            one of "gene_id", "Name", "transcript_id" for GFF/GTF.
        """
        if (
            attribute is not None and (
                x.file_type == 'gff' or
                isinstance(attribute, int)
            )
        ):
            names = [i[attribute] for i in x]

        else:
            if attribute is None:
                attribute = 'name'
            names = [getattr(i, attribute) for i in x]

        new = self.copy()
        new.data = new.data.reindex(names)
        return new

    def align_with(self, other):
        """
        Align the dataframe's index with another.
        """
        return self.__class__(self.data.reindex_like(other), **self._kwargs)

    def __len__(self):
        return len(self.data)

    def scatter(self, x, y, xfunc=None, yfunc=None, xscale=None, yscale=None,
                xlab=None, ylab=None, genes_to_highlight=None,
                marginal_histograms=False,
                general_kwargs=dict(color="k", alpha=0.2, picker=True),
                general_hist_kwargs=None, offset_kwargs={}, label_kwargs=None,
                ax=None, one_to_one=None, callback=None, hist_size=0.3,
                hist_pad=0.0, nan_offset=0.015, pos_offset=0.99,
                linelength=0.01, neg_offset=0.005):
        """
        Do-it-all method for making annotated scatterplots.

        Includes rugplots for NaN/Inf/-Inf, a default callback that prints the
        entries of the underlying dataframe when a point is clicked, point
        labeling options, marginal histograms for multiple subsets, arbitrary
        styling of points for arbitrary subsets, and more.

        Parameters
        ----------

        x, y : array-like
            Variables to plot.  Must be names in self.data's DataFrame.

        xfunc, yfunc : callable
            Functions to apply to `xvar` and `yvar` respectively. If xlab or
            ylab is not set separately, the function name will be used along
            with the column name to label the corresponding axis. This lets you
            play around with transformation functions (e.g., np.log, np.log1p)
            without having to add the corresponding column to the underlying
            dataframe.

        xlab, ylab : string
            Labels for x and y axes; default is to use function names for
            `xfunc` and `yfunc` and variable names `xvar` and `yvar`, e.g.,
            "log2(baseMeanA)"

        ax : None or Axes object
            If not None then plot on the provided Axes, otherwise create a new
            figure and axes.

        general_kwargs : dict
            Kwargs for matplotlib.scatter; specifies how all points look. Note
            that if you override this, you should include at least
            `picker=True` so that the callback function will work.

        genes_to_highlight : list of 2-tuples or 3-tuples
            Provides lots of control to colors.  It is a list of (`ind`,
            `kwargs`) tuples, where each `ind` specifies genes to plot with
            `kwargs`. `ind` is anything that can be used with `DataFrame.ix`.

            For example::

                [
                    (
                        x.log2FoldChange < 0,
                        dict(color='b', label='downregulated')
                    ),
                    (
                        x.log2FoldChange > 0,
                        dict(color='r', label='upregulated')
                    ),
                ]

            Each dictionary updates a copy of `general_kwargs`. If
            `genes_to_highlight` has a "name" kwarg, this must be a list that't
            the same length as `ind`.  It will be used to label the genes in
            `ind` using `label_kwargs`.

            For example::

                [
                    (
                        ['ENSG001', 'ENSG002'],
                        dict(color='r', name=['geneA', 'geneB'])
                    ),
                ]

            The tuples can also be 3-tuples of (ind, scatter_kwargs,
            hist_kwargs). The first two items act as above, and the third can
            be used to control histogram kwargs if `marginal_hists` is True.

            Note that, unless overridden, the color and alpha of the histograms
            will be inherited from the scatter kwargs, which in turn are
            inherited from general_kwargs. So the following will modifiy
            marginal histograms to have 100 bins::

                [
                    (
                        x.log2FoldChange < 0,
                        dict(color='b', label='down'),
                        dict(bins=100)
                    )
                ]

        callback : callable
            Function to call upon clicking a point. Must accept a single
            argument which is the gene ID. The function can do whatever it
            wants with it, but probably will want to access the underlying
            dataframe. Default is to print the corresponding row from the
            underlying dataframe.

        one_to_one : None or dict
            If not None, a dictionary of matplotlib.plot kwargs that will be
            used to plot a 1:1 line, e.g., dict(color='r', linestyle=':').

        label_kwargs : dict
            Kwargs for labeled genes, e.g., dict=(style='italic').  Will only
            be used if an entry in `genes_to_highlight` has a `name` key.

        offset_kwargs : dict
            Kwargs to be passed to matplotlib.transforms.offset_copy, used for
            adjusting the positioning of gene labels in relation to the actual
            point.

        xlab_prefix, ylab_prefix : str
            Optional label prefix that will be added to the beginning of `xlab`
            and/or `ylab`.

        marginal_histograms : bool
            If True, for each subset in `genes_to_highlight`, add marginal
            histograms along x and y axes subject to the various histogram
            controls below.

        hist_size : float
            Size of marginal histograms

        hist_pad : float
            Spacing between marginal histograms

        nan_offset, pos_offset, neg_offset : float
            Offset, in units of "fraction of axes" for the NaN, +inf, and -inf
            "rug plots"

        linelength : float
            Line length for the rug plots

        """
        _x = self.data[x]
        _y = self.data[y]

        # Construct defaults---------------------------------------------------
        def identity(x):
            return x.copy()

        # Axis label setup
        if xlab is None:
            xlab = x

            if xfunc is not None:
                xlab = "%s(%s)" % (xfunc.__name__, str(x))
            else:
                xlab = "%s" % (str(x))

        if ylab is None:
            ylab = y
            if yfunc is not None:
                ylab = "%s(%s)" % (yfunc.__name__, str(y))
            else:
                ylab = "%s" % (str(y))

        if xfunc is None:
            xfunc = identity

        if yfunc is None:
            yfunc = identity

        if general_kwargs is None:
            general_kwargs = {}

        if general_hist_kwargs is None:
            general_hist_kwargs = {}

        if genes_to_highlight is None:
            genes_to_highlight = []

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        if label_kwargs is None:
            label_kwargs = dict(
                horizontalalignment='right',
                verticalalignment='center',
                style='italic',
                bbox=dict(facecolor='w', edgecolor='None', alpha=0.5)
            )

        # Clean data ---------------------------------------------------------
        xi = xfunc(_x)
        yi = yfunc(_y)

        # handle inf, -inf, and NaN so we can build the rug plots appropriately
        x_is_pos_inf = np.isinf(xi) & (xi > 0)
        x_is_neg_inf = np.isinf(xi) & (xi < 0)
        x_is_nan = np.isnan(xi)
        y_is_pos_inf = np.isinf(yi) & (yi > 0)
        y_is_neg_inf = np.isinf(yi) & (yi < 0)
        y_is_nan = np.isnan(yi)

        # Indexes for valid values
        x_valid = ~(x_is_pos_inf | x_is_neg_inf | x_is_nan)
        y_valid = ~(y_is_pos_inf | y_is_neg_inf | y_is_nan)

        # global min/max
        gmin = max(xi[x_valid].min(), yi[y_valid].min())
        gmax = min(xi[x_valid].max(), yi[y_valid].max())

        # Now we remove any genes from in allind (which will be plotted using
        # `general_kwargs`) that will be plotted by something in
        # genes_to_highlight.  This avoids double-plotting.
        allind = np.zeros_like(xi) == 0
        for block in genes_to_highlight:
            ind = block[0]
            allind[ind] = False

        # Copy over the color and alpha from general_kwargs (the scatter
        # points) to general_hist_kwargs (marginal histograms) if they're not
        # otherwise specified
        general_hist_kwargs = utils.updatecopy(
            orig=general_hist_kwargs, update_with=general_kwargs,
            keys=['color', 'alpha'])

        # Put the non-highlighted genes at the beginning of _genes_to_highlight
        # list so we can just iterate over one list. If genes_to_highlight
        # styles don't otherwise override z-order, this means the "background"
        # gets plotted first and therefore under the highlighted points.
        _genes_to_highlight = genes_to_highlight[:]
        _genes_to_highlight.insert(
            0,
            (allind, general_kwargs, general_hist_kwargs)
        )

        # Set up the object that will handle the marginal histograms and
        # plotting.
        self.marginal = MarginalHistScatter(
            ax, hist_size=hist_size, pad=hist_pad)

        # Set up kwargs for x and y rug plots
        rug_x_kwargs = dict(
            linelength=linelength,
            transform=blended_transform_factory(ax.transData, ax.transAxes),
        )
        rug_y_kwargs = dict(
            linelength=linelength,
            transform=blended_transform_factory(ax.transAxes, ax.transData),
            orientation='vertical',
        )

        # EventCollection objects need a color as a 3-tuple, so set up
        # a converter here.
        color_converter = matplotlib.colors.ColorConverter().to_rgb

        # Plot the one-to-one line, if kwargs were specified. If zorder is not
        # overridden, points will overlay it.
        if one_to_one:
            ax.plot([gmin, gmax],
                    [gmin, gmax],
                    **one_to_one)

        # Plot 'em all, and label if specified

        # In order to avoid calling the callback function multiple times when
        # we have overlapping genes to highlight (e.g., a gene that is both
        # upregulated AND has a peak), here we keep track of everything that's
        # been added so far.
        self._seen = np.ones_like(xi) == 0
        for block in _genes_to_highlight:
            ind = block[0]
            kwargs = block[1]

            if len(block) == 3:
                hist_kwargs = block[2]
            else:
                hist_kwargs = {}

            names = kwargs.pop('names', None)

            _marginal_histograms = (
                kwargs.pop('marginal_histograms', False) or
                marginal_histograms)

            updated_kwargs = utils.updatecopy(
                orig=kwargs, update_with=general_kwargs)

            updated_hist_kwargs = utils.updatecopy(
                orig=hist_kwargs, update_with=general_hist_kwargs)

            updated_hist_kwargs = utils.updatecopy(
                orig=updated_hist_kwargs, update_with=kwargs,
                keys=['color', 'alpha'], override=True)

            xhist_kwargs = updated_kwargs.pop('xhist_kwargs', None)
            yhist_kwargs = updated_kwargs.pop('yhist_kwargs', None)

            # MarginalHistScatter does most of the actual plotting work.
            scatter_ind = ind & x_valid & y_valid
            self.marginal.append(
                xi[scatter_ind],
                yi[scatter_ind],
                scatter_kwargs=dict(**updated_kwargs),
                hist_kwargs=updated_hist_kwargs,
                xhist_kwargs=xhist_kwargs,
                yhist_kwargs=yhist_kwargs,
                marginal_histograms=_marginal_histograms,
            )

            # Callback functions have access to the Collection object that was
            # picked as well as the index into the collection. Since ultimately
            # we want to access a row of the dataframe, to the collection we
            # attach the dataframe itself and the index into that dataframe
            # represented by the scatter points in the collection.
            coll = self.marginal.scatter_ax.collections[-1]
            coll.df = self.data
            coll.ind = scatter_ind

            color = color_converter(updated_kwargs['color'])
            rug_x_kwargs['color'] = color
            rug_y_kwargs['color'] = color

            # Note: if both x and y are not valid, then they will not be on the
            # plot at all (since there's no obvious place to put them).
            items = [
                # top rug, y is +inf and x is valid
                (xi, ind & x_valid & y_is_pos_inf, pos_offset, rug_x_kwargs),

                # one of the bottom rugs, where y is NaN
                (xi, ind & x_valid & y_is_nan, nan_offset, rug_x_kwargs),

                # bottom rug, y is -inf
                (xi, ind & x_valid & y_is_neg_inf, neg_offset, rug_x_kwargs),

                # right rug, x is +inf
                (yi, ind & y_valid & x_is_pos_inf, pos_offset, rug_y_kwargs),

                # one of the left rugs; x is NaN
                (yi, ind & y_valid & x_is_nan, nan_offset, rug_y_kwargs),

                # left rug, x is -inf
                (yi, ind & y_valid & x_is_neg_inf, neg_offset, rug_y_kwargs),
            ]
            for values, index, offset, kwargs in items:

                # provide np.array rather than pandas.Series
                coll = EventCollection(
                    values[index].values, lineoffset=offset, **kwargs)
                coll.df = self.data
                coll.ind = index
                ax.add_collection(coll)

            # Plot the text names if configured
            if names:
                transOffset = matplotlib.transforms.offset_copy(
                    ax.transData, fig=ax.figure, **offset_kwargs)

                for xii, yii, name in zip(xi[ind], yi[ind], names):
                    ax.text(xii,
                            yii,
                            name,
                            transform=transOffset,
                            **label_kwargs)

        # register callback
        if callback is None:
            callback = self._default_callback

        # This trick is so that we can say, in the docstring, that a callback
        # function must accept an index name (i.e. gene ID). Much easier than
        # saying "the callback function accepts an event, and the event's
        # artist has an `ind` and `df` attribute".
        def wrapped_callback(event):
            for _id in self._id_callback(event):
                callback(_id)

        # Connect the callback.
        ax.figure.canvas.mpl_connect('pick_event', wrapped_callback)

        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.axis('tight')

        return ax

    def radviz(self, column_names, transforms=dict(), **kwargs):
        """
        Radviz plot.

        Useful for exploratory visualization, a radviz plot can show
        multivariate data in 2D.  Conceptually, the variables (here, specified
        in `column_names`) are distributed evenly around the unit circle.  Then
        each point (here, each row in the dataframe) is attached to each
        variable by a spring, where the stiffness of the spring is proportional
        to the value of corresponding variable.  The final position of a point
        represents the equilibrium position with all springs pulling on it.

        In practice, each variable is normalized to 0-1 (by subtracting the
        mean and dividing by the range).

        This is a very exploratory plot.  The order of `column_names` will
        affect the results, so it's best to try a couple different orderings.
        For other caveats, see [1].

        Additional kwargs are passed to self.scatter, so subsetting, callbacks,
        and other configuration can be performed using options for that method
        (e.g., `genes_to_highlight` is particularly useful).

        Parameters
        ----------
        column_names : list
            Which columns of the dataframe to consider.  The columns provided
            should only include numeric data, and they should not contain any
            NaN, inf, or -inf values.

        transforms : dict
            Dictionary mapping column names to transformations that will be
            applied just for the radviz plot.  For example, np.log1p is
            a useful function. If a column name is not in this dictionary, it
            will be used as-is.

        ax : matplotlib.Axes
            If not None, then plot the radviz on this axes.  If None, then
            a new figure will be created.

        kwargs : dict
            Additional arguments are passed to self.scatter.  Note that not all
            possible kwargs for self.scatter are necessarily useful for
            a radviz plot (for example, margninal histograms would not be
            meaningful).

        Notes
        -----
        This method adds two new variables to self.data: "radviz_x" and
        "radviz_y".  It then calls the self.scatter method, using these new
        variables.

        The data transformation was adapted from the
        pandas.tools.plotting.radviz function.

        References
        ----------
        [1]  Hoffman,P.E. et al. (1997) DNA visual and analytic data mining. In
             the Proceedings of the IEEE Visualization. Phoenix, AZ, pp.
             437-441.
        [2] http://www.agocg.ac.uk/reports/visual/casestud/brunsdon/radviz.htm
        [3] http://pandas.pydata.org/pandas-docs/stable/visualization.html\
                #radviz
        """
        # make a copy of data
        x = self.data[column_names].copy()

        for k, v in transforms.items():
            x[k] = v(x[k])

        def normalize(series):
            mn = min(series)
            mx = max(series)
            return (series - mn) / (mx - mn)

        df = x.apply(normalize)

        to_plot = []

        n = len(column_names)

        s = np.array([(np.cos(t), np.sin(t))
                      for t in [2.0 * np.pi * (i / float(n))
                                for i in range(n)]])

        for i in range(len(x)):
            row = df.irow(i).values
            row_ = np.repeat(np.expand_dims(row, axis=1), 2, axis=1)
            to_plot.append((s * row_).sum(axis=0) / row.sum())

        x_, y_ = zip(*to_plot)
        self.data['radviz_x'] = x_
        self.data['radviz_y'] = y_

        ax = self.scatter('radviz_x', 'radviz_y', **kwargs)

        ax.add_patch(patches.Circle((0.0, 0.0), radius=1.0, facecolor='none'))
        for xy, name in zip(s, column_names):
            ax.add_patch(patches.Circle(xy, radius=0.025, facecolor='gray'))
            if xy[0] < 0.0 and xy[1] < 0.0:
                ax.text(xy[0] - 0.025, xy[1] - 0.025, name,
                        ha='right', va='top', size='small')
            elif xy[0] < 0.0 and xy[1] >= 0.0:
                ax.text(xy[0] - 0.025, xy[1] + 0.025, name,
                        ha='right', va='bottom', size='small')
            elif xy[0] >= 0.0 and xy[1] < 0.0:
                ax.text(xy[0] + 0.025, xy[1] - 0.025, name,
                        ha='left', va='top', size='small')
            elif xy[0] >= 0.0 and xy[1] >= 0.0:
                ax.text(xy[0] + 0.025, xy[1] + 0.025, name,
                        ha='left', va='bottom', size='small')

        ax.axis('equal')
        return ax

    def _id_callback(self, event):
        # event.ind is the index into event's x and y data.
        #
        # event.artist.ind is the index of the entire artist into the original
        # dataframe.
        subset_df = event.artist.df.ix[event.artist.ind]
        for i in event.ind:
            _id = subset_df.index[i]
            yield _id

    def _default_callback(self, i):
        print(self.data.ix[i])


class DifferentialExpressionResults(ResultsTable):

    __doc__ = _base_doc % dedent("""
    A ResultsTable subclass for working with differential expression results.

    Adds methods for up/down regulation, ma_plot, and sets class variables for
    which columns should be considered for pval, log fold change, and mean
    values. This class acts as a parent for subclasses like DESeqResults,
    EdgeRResults, and others/
    """)
    pval_column = 'padj'
    lfc_column = 'log2FoldChange'
    mean_column = 'baseMean'

    def __init__(self, data, db=None, header_check=True, **kwargs):
        import_kwargs = kwargs.pop('import_kwargs', {})
        if header_check and isinstance(data, str):
            comment_char = import_kwargs.get('comment', '#')
            import_kwargs['comment'] = comment_char
        import_kwargs['na_values'] = ['nan']
        import_kwargs['index_col'] = import_kwargs.pop('index_col', 0)
        super(DifferentialExpressionResults, self).__init__(
            data=data, db=db, import_kwargs=import_kwargs, **kwargs)

    def changed(self, alpha=0.05, lfc=0, idx=True):
        """
        Changed features.

        Helper function to get where the pval is <= alpha and the
        absolute value log2foldchange is >= lfc.

        Parameters
        ----------
        alpha : float

        lfc : float

        idx : bool
            If True, a boolean index will be returned.  If False, a new object
            will be returned that has been subsetted.
        """
        ind = (
            (self.data[self.pval_column] <= alpha) &
            (np.abs(self.data[self.lfc_column]) >= lfc)
        )

        if idx:
            return ind
        return self[ind]

    def unchanged(self, alpha=0.05, lfc=0, idx=True):
        """
        Unchanged features.

        Helper function to get where the pval is > alpha and the
        absolute value of the log2foldchange is < lfc.

        Parameters
        ----------
        alpha : float

        lfc : float

        idx : bool
            If True, a boolean index will be returned.  If False, a new object
            will be returned that has been subsetted.
        """
        ind = ~self.changed(alpha, lfc, idx)
        if idx:
            return ind
        return self[ind]

    def upregulated(self, alpha=0.05, lfc=0, idx=True):
        """
        Upregulated features.

        Helper function to get where the pval is <= alpha and the
        log2foldchange is >= lfc.

        Parameters
        ----------
        alpha : float

        lfc : float

        idx : bool
            If True, a boolean index will be returned.  If False, a new object
            will be returned that has been subsetted.
        """
        ind = (
            (self.data[self.pval_column] <= alpha) &
            (self.data[self.lfc_column] >= lfc)
        )
        if idx:
            return ind
        return self[ind]

    def downregulated(self, alpha=0.05, lfc=0, idx=True):
        """
        Downregulated features.

        Helper function to get where the pval is <= alpha and the
        log2foldchange is <= lfc.

        Parameters
        ----------
        alpha : float

        lfc : float

        idx : bool
            If True, a boolean index will be returned.  If False, a new object
            will be returned that has been subsetted.
        """
        ind = (
            (self.data[self.pval_column] <= alpha) &
            (self.data[self.lfc_column] <= lfc)
        )
        if idx:
            return ind
        return self[ind]

    def ma_plot(self, alpha, up_kwargs=None, dn_kwargs=None,
                zero_line=None, **kwargs):
        """
        MA plot.

        Plots the average read count across treatments (x-axis) vs the log2
        fold change (y-axis).

        Additional kwargs are passed to self.scatter (useful ones might include
        `genes_to_highlight`)

        Parameters
        ----------
        alpha : float
            Features with values <= `alpha` will be highlighted in the plot.

        up_kwargs, dn_kwargs : None or dict
            Kwargs passed to matplotlib's scatter(), used for styling up/down
            regulated features (defined by `alpha` and `col`)

        zero_line : None or dict
            Kwargs passed to matplotlib.axhline(0).

        """
        genes_to_highlight = kwargs.pop('genes_to_highlight', [])
        genes_to_highlight.append(
            (self.upregulated(alpha),
             up_kwargs or dict(color='r')))
        genes_to_highlight.append(
            (self.downregulated(alpha),
             dn_kwargs or dict(color='b')))
        if zero_line is None:
            zero_line = {}
        x = self.mean_column
        y = self.lfc_column

        if 'xfunc' not in kwargs:
            kwargs['xfunc'] = np.log
        ax = self.scatter(
            x=x,
            y=y,
            genes_to_highlight=genes_to_highlight,
            **kwargs)
        if zero_line:
            ax.axhline(0, **zero_line)
        return ax


class EdgeRResults(DifferentialExpressionResults):
    __doc__ = _base_doc % dedent(
        """
        Class for working with results from edgeR.

        Just like a DifferentialExpressionResults object, but sets the
        pval_column, lfc_column, and mean_column to the names used in edgeR's
        output.
        """)
    pval_column = 'FDR'
    lfc_column = 'logFC'
    mean_column = 'logCPM'


class DESeqResults(DifferentialExpressionResults):
    __doc__ = _base_doc % dedent(
        """
        Class for working with results from DESeq.

        Just like a DifferentialExpressionResults object, but sets the
        pval_column, lfc_column, and mean_column to the names used in DESeq
        (v1) output.
        """)

    def colormapped_bedfile(self, genome, cmap=None):
        """
        Create a BED file with padj encoded as color

        Features will be colored according to adjusted pval (phred
        transformed).  Downregulated features have the sign flipped.

        Parameters
        ----------
        cmap : matplotlib colormap
            Default is matplotlib.cm.RdBu_r

        Notes
        -----
        Requires a FeatureDB to be attached.
        """
        if self.db is None:
            raise ValueError("FeatureDB required")
        db = gffutils.FeatureDB(self.db)

        def scored_feature_generator(d):
            for i in range(len(d)):
                try:
                    feature = db[d.ix[i]]
                except gffutils.FeatureNotFoundError:
                    raise gffutils.FeatureNotFoundError(d.ix[i])
                score = -10 * np.log10(d.padj[i])
                lfc = d.log2FoldChange[i]
                if np.isnan(lfc):
                    score = 0
                if lfc < 0:
                    score *= -1
                feature.score = str(score)
                feature = featurefuncs.extend_fields(
                    featurefuncs.gff2bed(
                        gffutils.helpers.asinterval(feature)), 9)
                fields = feature.fields[:]
                fields[6] = fields[1]
                fields[7] = fields[2]
                fields.append(str(d.padj[i]))
                fields.append(str(d.pval[i]))
                fields.append('%.3f' % d.log2FoldChange[i])
                fields.append('%.3f' % d.baseMeanB[i])
                fields.append('%.3f' % d.baseMeanB[i])
                yield pybedtools.create_interval_from_list(fields)

        x = pybedtools.BedTool(scored_feature_generator(self)).saveas()
        norm = x.colormap_normalize()
        if cmap is None:
            cmap = matplotlib.cm.RdBu_r
        cmap = colormap_adjust.cmap_center_point_adjust(
            cmap, [norm.vmin, norm.vmax], 0)

        def score_zeroer(f):
            f.score = '0'
            return f
        return x.each(featurefuncs.add_color, cmap=cmap, norm=norm)\
                .sort()\
                .each(score_zeroer)\
                .truncate_to_chrom(genome)\
                .saveas()

    def autosql_file(self):
        """
        Generate the autosql for DESeq results (to create bigBed)

        Returns a temp filename containing the autosql defining the extra
        fields.

        This for creating bigBed files from BED files created by
        colormapped_bed.  When a user clicks on a feature, the DESeq results
        will be reported.
        """
        fn = pybedtools.BedTool._tmp()

        AUTOSQL = dedent(
            """
            table example
            "output from DESeq"
            (
            string  chrom;  "chromosome"
            uint chromStart; "start coord"
            uint chromEnd; "stop coord"
            string name; "name of feature"
            uint score; "always zero"
            char[1] strand; "+ or - for strand"
            uint    thickStart; "Coding region start"
            uint    thickEnd;  "Coding region end"
            uint reserved; "color according to score"
            string padj; "DESeq adjusted p value"
            string pval; "DESeq raw p value"
            string logfoldchange; "DESeq log2 fold change"
            string basemeana; "DESeq baseMeanA"
            string basemeanb; "DESeq baseMeanB"
        )
        """)

        fout = open(fn, 'w')
        fout.write(AUTOSQL)
        fout.close()
        return fn


class DESeq2Results(DESeqResults):
    __doc__ = _base_doc % dedent(
        """
        Class for working with results from DESeq2.

        Just like a DifferentialExpressionResults object, but sets the
        pval_column, lfc_column, and mean_column to the names used in edgeR's
        output.
        """)
    pval_column = 'padj'
    lfc_column = 'log2FoldChange'
    mean_column = 'baseMean'


class LazyDict(object):
    def __init__(self, fn_dict, index_file=None, index_from=None, extra=None,
                 cls=DESeqResults):
        """
        Dictionary-like object that lazily-loads ResultsTable objects.

        Primary use-case for this is for organizing a large number of
        differential expression results that you'll be cross-comparing and 1)
        you don't want to immediately load EVERYTHING into a dataframe, and 2)
        you want to ensure the indexes are aligned across dataframes.

        Only upon accessing a value from the LazyDict will the data be loaded
        into a DataFrame.

        Parameters
        ----------
        fn_dict : dict
            Keys of `fn_dict` will be the keys of this LazyDict object.  Values
            should be filenames which will be loaded into ResultsTable object
            upon access for the first time.

        index_file : str
            Path to a file that contains one ID per line.  This file is used to
            ensure all ResultsTable objects are aligned to the same index. If
            you don't want to provide this, then set it to None and consider
            the `index_from` kwarg.

        index_from : str
            Key of the dataframe whose index will be used to align other
            dataframes. If this is provided, this dataframe will have to be
            loaded before any others, so the first access will trigger two
            dataframe loads.

        cls : ResultsTable class or subclass
            Each filename in `fn_dict` will be converted using this class.

        """
        self.fn_dict = fn_dict

        # this acts as the cache
        self._dict = {}

        if index_file is not None and index_from is not None:
            raise ValueError(
                "Only one of `index_file` or `index_from` should be provided."
            )
        self.index_file = index_file
        self.index_from = index_from

        if index_file:
            self.index = [i.strip() for i in open(index_file)]
        else:
            self.index = None
        self._cls = cls

    def _load(self, key):
        return self._cls(self.fn_dict[key])

    def __getitem__(self, key):
        if self.index is None:
            obj = self._load(self.index_from)
            self._dict[self.index_from] = obj
            self.index = obj.index

        if key not in self._dict:
            obj = self._load(key)
            obj.data = obj.data.ix[self.index]
            self._dict[key] = obj

        return self._dict[key]

    def __repr__(self):
        s = "<%s> with possible keys\n:%s\n" \
            % (self.__class__.__name__, self.fn_dict.keys())
        s += "and existing keys:\n"
        s += repr(self._dict)
        return s

    def keys(self):
        return self.fn_dict.keys()

    def values(self):
        return [self._dict[key] for key in self.keys()]

    def items(self):
        return list((key, self._dict[key]) for key in self.keys())


class MarginalHistScatter(object):
    def __init__(self, ax, hist_size=0.6, pad=0.05):
        """
        Class to enable incremental appending of scatterplots, each of which
        generate additional marginal histograms.

        The `append` method does most of the work.
        """
        self.scatter_ax = ax
        self.fig = ax.figure

        self.divider = make_axes_locatable(self.scatter_ax)

        self.top_hists = []
        self.right_hists = []
        self.hist_size = hist_size
        self.pad = pad
        self.xfirst_ax = None
        self.yfirst_ax = None

        # will hold histogram data
        self.hxs = []
        self.hys = []

    @property
    def xmax(self):
        return self.scatter_ax.dataLim.xmax

    @property
    def ymax(self):
        return self.scatter_ax.dataLim.ymax

    @property
    def xmin(self):
        return self.scatter_ax.dataLim.xmin

    @property
    def ymin(self):
        return self.scatter_ax.dataLim.ymin

    @property
    def limits(self):
        return (self.xmin, self.xmax, self.ymin, self.ymax)

    def append(self, x, y, scatter_kwargs, hist_kwargs=None, xhist_kwargs=None,
               yhist_kwargs=None, num_ticks=3, labels=None, hist_share=False,
               marginal_histograms=True):
        """
        Adds a new scatter to self.scatter_ax as well as marginal histograms
        for the same data, borrowing addtional room from the axes.

        Parameters
        ----------

        x, y : array-like
            Data to be plotted

        scatter_kwargs : dict
            Keyword arguments that are passed directly to scatter().

        hist_kwargs : dict
            Keyword arguments that are passed directly to hist(), for both the
            top and side histograms.

        xhist_kwargs, yhist_kwargs : dict
            Additional, margin-specific kwargs for the x or y histograms
            respectively.  These are used to update `hist_kwargs`

        num_ticks : int
            How many tick marks to use in each histogram's y-axis

        labels : array-like
            Optional NumPy array of labels that will be set on the collection
            so that they can be accessed by a callback function.

        hist_share : bool
            If True, then all histograms will share the same frequency axes.
            Useful for showing relative heights if you don't want to use the
            hist_kwarg `normed=True`

        marginal_histograms : bool
            Set to False in order to disable marginal histograms and just use
            as a normal scatterplot.
        """
        scatter_kwargs = scatter_kwargs or {}
        hist_kwargs = hist_kwargs or {}
        xhist_kwargs = xhist_kwargs or {}
        yhist_kwargs = yhist_kwargs or {}
        yhist_kwargs.update(dict(orientation='horizontal'))

        # Plot the scatter
        coll = self.scatter_ax.scatter(x, y, **scatter_kwargs)
        coll.labels = labels

        if not marginal_histograms:
            # that was easy!
            return

        xhk = utils.updatecopy(hist_kwargs, xhist_kwargs)
        yhk = utils.updatecopy(hist_kwargs, yhist_kwargs)

        axhistx = self.divider.append_axes(
            'top', size=self.hist_size,
            pad=self.pad, sharex=self.scatter_ax, sharey=self.xfirst_ax)

        axhisty = self.divider.append_axes(
            'right', size=self.hist_size,
            pad=self.pad, sharey=self.scatter_ax, sharex=self.yfirst_ax)

        axhistx.yaxis.set_major_locator(
            MaxNLocator(nbins=num_ticks, prune='both'))

        axhisty.xaxis.set_major_locator(
            MaxNLocator(nbins=num_ticks, prune='both'))

        if not self.xfirst_ax and hist_share:
            self.xfirst_ax = axhistx

        if not self.yfirst_ax and hist_share:
            self.yfirst_ax = axhisty

        # Keep track of which axes are which, because looking into fig.axes
        # list will get awkward....
        self.top_hists.append(axhistx)
        self.right_hists.append(axhisty)

        # Scatter will deal with NaN, but hist will not.  So clean the data
        # here.
        hx = x[np.isfinite(x)]
        hy = x[np.isfinite(y)]

        self.hxs.append(hx)
        self.hys.append(hy)

        # Only plot hists if there's valid data
        if len(hx) > 0:
            if len(hx) == 1:
                _xhk = utils.updatecopy(
                    orig=xhk, update_with=dict(bins=[hx[0], hx[0]]),
                    keys=['bins'])
                axhistx.hist(hx, **_xhk)
            else:
                axhistx.hist(hx, **xhk)
        if len(hy) > 0:
            if len(hy) == 1:
                _yhk = utils.updatecopy(
                    orig=yhk, update_with=dict(bins=[hy[0], hy[0]]),
                    keys=['bins'])
                axhisty.hist(hy, **_yhk)
            else:
                axhisty.hist(hy, **yhk)

        # Turn off unnecessary labels -- for these, use the scatter's axes
        # labels
        for txt in axhisty.get_yticklabels() + axhistx.get_xticklabels():
            txt.set_visible(False)

        for txt in axhisty.get_xticklabels():
            txt.set_rotation(-90)

    def add_legends(self, xhists=True, yhists=False, scatter=True, **kwargs):
        """
        Add legends to axes.
        """
        axs = []
        if xhists:
            axs.extend(self.hxs)
        if yhists:
            axs.extend(self.hys)
        if scatter:
            axs.extend(self.ax)

        for ax in axs:
            ax.legend(**kwargs)
