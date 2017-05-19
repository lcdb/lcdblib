import os
import sys
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
from lcdblib.plotting import results_table
from lcdblib.stats import fisher
import pybedtools
import pandas
import argh
from argh import arg


def plot(de_results, regions=None, peaks=None, selected=None, x='baseMean',
         y='log2FoldChange', disable_logx=False, logy=False, pval_col='padj',
         alpha=0.1, lfc_cutoff=0, plot_filename=None,
         disable_raster_points=False, genes_to_label=None, label_column=None,
         report=None, gene_lists=None
    ):
    """
    M-A plot showing up- and downregulated genes with optional labeling and
    Fishers exact tests.

    If --plot-filename is not specified, then the plot will be displayed and
    points can be clicked for interactive exploration.

    If --peaks and --regions are specified, then results from Fishers exact
    tests will be printed to stdout, or to --report if specified.

    Parameters
    ----------
    de_results : str or pandas.DataFrame
        If str, it's the filename of a TSV of differential expression results,
        with first column as gene ID. It will be parsed into a dataframe where
        the index is gene ID. When called as a library, an already-created
        pandas.DataFrame can optionally be provided instead.

    regions : str or pybedtools.BedTool
        Gene regions in which to look for intersections with peaks. BED file
        where the 4th column contains gene IDs that are also present in first
        column of `de_results`. Typically this would be a BED file of promoters
        or gene bodies. When called as a library, a pybedtools.BedTool object
        can optionally be provided instead.

    peaks : str or pybedtools.BedTool
        BED file to be intersected with `regions`. When called as a library,
        a pybedtools.BedTool object can optionally be provided instead.

    selected : str or list-like
        Replaces `regions` `peaks` arguments; useful for when you already know
        which genes you want to select (e.g., upregulated from a different
        experiment). If a string, assume it's a filename and use the first
        column which will be used as an index into the `de_results` dataframe.
        When called as a library, if `selected` is not a string it will be used
        as an index into the dataframe.

    x : str
        Column to use for x-axis. Default of "baseMean" expects DESeq2
        results

    y : str
        Column to use for y-axis. Default of "log2FoldChange" expects DESeq2
        results

    disable_logx : bool
        Disable default behavior of transforming x values using log10

    logy : bool
        Transform y values using log2

    pval-col : str
        Column to use for statistical significance. Default "padj" expectes
        DESeq2 results.

    alpha : float
        Threshold for calling significance. Applied to `pval_col`

    lfc_cutoff : float
        Log2fold change cutoff to be applied to y values. Threshold is applied
        post-transformation, if any specified (e.g., `logy` argument).

    plot_filename : str
        File to save plot. Format auto-detected by extension. Output directory
        will be created if needed.

    disable_raster_points : bool
        Disable the default behavior of rasterizing points in a PDF. Use
        sparingly, since drawing 30k+ individual points in a PDF may slow down
        your machine.

    genes_to_label: str or list-like
        Optional file containing genes to label with text. First column must be
        a subset of the first column of `de_results`. Lines starting with '#'
        and subsequent tab-separated columns will be ignored. When called as
        a library, a list-like object of gene IDs can be provided.

    label_column : str
        Optional column from which to take gene labels found in
        `genes_to_label` (e.g., "symbol"). If the value in this column is
        missing, fall back to the index. Use this if your gene IDs are long
        Ensembl IDs but you want the gene symbols to show up on the plot.

    report : str
        Where to write out Fisher's exact test results. Default is stdout

    gene_lists : str
        Prefix to gene lists. If specified, gene lists corresponding to the
        cells of the 2x2 Fishers exact test will be written to
        {prefix}.up.tsv and {prefix}.dn.tsv. These are subsets of `de_results`
        where genes are up and have a peak in region (or are selected), or
        downregulated and have a peak in region (or are selected),
        respectively.
    """
    rasterized = not disable_raster_points
    rt = results_table.DESeq2Results(de_results, import_kwargs=dict(index_col=0))
    up = rt.upregulated(alpha=alpha, lfc=lfc_cutoff)
    dn = rt.downregulated(alpha=alpha, lfc=-lfc_cutoff)
    un = ~(up | dn)
    sns.set_context('talk')
    sns.set_style('white')
    general_kwargs=dict(marker='.', alpha=0.4, color='0.5', picker=5,
                        label="_", linewidth=0, rasterized=rasterized)
    genes_to_highlight = [
        (
            up,
            dict(color='#990000', marker='o', s=40, alpha=0.5, label='up (%s)' % sum(up))
        ),
        (
            dn,
            dict(color='#005c99', marker='o', s=40, alpha=0.5, label='down (%s)' % sum(dn))
        ),
    ]

    if genes_to_label:
        _genes_to_label = []
        if isinstance(genes_to_label, str):
            for i in open(genes_to_label):
                if i.startswith('#'):
                    continue
                _genes_to_label.append(i.split('\t')[0].strip())
        else:
            _genes_to_label = genes_to_label

        ind = rt.data.index.isin(_genes_to_label)

        # Don't add labels if a coordinate is null
        ind = ind & ~(rt.data[y].isnull() | rt.data[x].isnull())

        if label_column:
            names = rt.data.loc[ind, label_column]

            # Fill in nulls with index to avoid labeling all genes with no
            # symbol as "nan"
            n = names.isnull()
            names[n] = rt.index[n]

        else:
            names = rt.index[ind]

        names = list(names)

        genes_to_highlight.append(
            (
                ind,
                dict(rasterized=rasterized, names=names, facecolor='None',
                     alpha=1.0, s=160, linewidth=1, zorder=100, label='_')
            )
        )

    if not disable_logx:
        xfunc=np.log10
    else:
        xfunc=None

    if report is None:
        output = sys.stdout
    else:
        output = open(report, 'w')

    if selected and (peaks or regions):
        raise ValueError(
            "`selected` is mutually exclusive with `peaks` and `regions`")

    do_fisher = False

    if selected:
        do_fisher = True
        if isinstance(selected, str):
            selected = list(pandas.read_table(selected, index_col=0).index)
        selected_genes = rt.index.isin(selected)
        row_names = ['selected', 'not selected']

    elif peaks is not None and regions is not None:
        do_fisher = True
        row_names = ['has peak', 'no peak']
        regions = pybedtools.BedTool(regions)
        peaks = pybedtools.BedTool(peaks)
        with_peak = list(set([i.name for i in regions.intersect(peaks, u=True)]))
        in_region = peaks.intersect(regions, u=True)

        selected_genes = rt.index.isin(with_peak)

        npeaks = len(peaks)
        nregions = len(regions)
        npeaks_in_region = len(in_region)

        output.write('Total peaks: {}\n'.format(npeaks))
        output.write('Peaks in regions: {0} ({1:.2f}%)\n\n'.format(npeaks_in_region, npeaks_in_region / npeaks * 100))

    if do_fisher:

        genes_to_highlight.append(
            (
                selected_genes,
                dict(color='#ff9900', alpha=0.8, s=30,
                    rasterized=rasterized, label='{0} ({1})'.format(row_names[0], sum(selected_genes)))
            )
        )

        output.write(
            fisher.fisher_tables(
                table=fisher.table_from_bool(selected_genes, up),
                row_names=row_names,
                col_names=['upregulated', 'not'],
                title='Upregulated (lfc>{0}; padj<{1})'.format(lfc_cutoff, alpha)))

        output.write('\n\n')
        output.write(
            fisher.fisher_tables(
                table=fisher.table_from_bool(selected_genes, dn),
                row_names=row_names,
                col_names=['downregulated', 'not'],
                title='Downregulated (lfc<-{0}; padj<{1})'.format(lfc_cutoff, alpha)))

        if gene_lists is not None:
            rt.data[selected_genes & up].to_csv(gene_lists + '.up.tsv', sep='\t')
            rt.data[selected_genes & dn].to_csv(gene_lists + '.dn.tsv', sep='\t')

        if report is not None:
            output.close()

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    ax = rt.scatter(
        ax=ax,
        x=x,
        y=y,
        xfunc=xfunc,
        genes_to_highlight=genes_to_highlight,
        general_kwargs=general_kwargs,
        offset_kwargs=dict(x=-0.1),
        label_kwargs=dict(
            horizontalalignment='right',
            verticalalignment='center',
            style='italic',
            size=10,
            zorder=500,
            bbox=dict(facecolor='w', alpha=0.5, edgecolor='none'))
    )

    ax.legend(loc='best', prop=dict(size=10))

    if plot_filename:
        dirname = os.path.dirname(plot_filename)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        fig.savefig(plot_filename)
    return ax

if __name__ == "__main__":
    argh.dispatch_command(plot)
    plt.show()
