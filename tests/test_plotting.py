import matplotlib
matplotlib.use('agg')
from lcdblib.plotting import results_table
import pytest
from textwrap import dedent
import pandas
from io import StringIO
import pybedtools
import gffutils


@pytest.fixture
def deseq_results():
    data = dedent("""
          baseMean       log2FoldChange  lfcSE      stat        pvalue         padj
    g1    6088.817976     3.253          0.081666   39.844574   0.000000e+00   0.000000e+00
    g2    108593.082921  -3.802          0.095365  -39.874870   0.000000e+00   0.000000e+00
    g3    3035.710964     3.788          0.086147   43.976691   0.000000e+00   0.000000e+00
    g4    3416.026005    -3.311          0.093745  -35.324175   2.499055e-273  1.047854e-269
    g5    67599.565846    2.517          0.073089   34.441404   6.056132e-260  2.031469e-256
    g6    4.009425       -0.096          0.476756  -0.203350    8.388611e-01   8.942725e-01
    g7    0.643246        0.157          0.304968   0.515483    6.062155e-01   NaN
    g8    0.157905       -0.022          0.212231  -0.107794    9.141590e-01   NaN
    g9    0.000000        NaN            NaN        NaN         NaN            NaN
    g10   17.039983      -1.009          0.391316  -2.578720    9.916715e-03   2.265983e-02
    g11   45.651073      -2.473          0.314267  -7.869736    3.553900e-15   3.592888e-14
    g12   0.000000        NaN            NaN        NaN         NaN            NaN
    g13   21.854535       0.017          0.364342   0.048802    9.610770e-01   9.754423e-01
    g14   64.838568      -0.838          0.243228  -3.446669    5.675434e-04   1.661518e-03
    g15   0.000000        NaN            NaN        NaN         NaN            NaN
    g16   0.000000        NaN            NaN        NaN         NaN            NaN
    g17   0.000000        NaN            NaN        NaN         NaN            NaN
    g18   3.691108        1.141          0.473275   2.411347    0.015894       0.034530
    g19   0.000000        NaN            NaN        NaN         NaN            NaN
    g20   0.000000        NaN            NaN        NaN         NaN            NaN
    g21   2.968828       -0.876          0.462677  -1.894393    0.058173       0.107702
    g22   0.000000        NaN            NaN        NaN         NaN            NaN
    g23   0.157905       -0.022          0.212231  -0.107794    0.914159       NaN
    g24   13.196089      -0.890          0.426685  -2.087471    0.036846       0.072439
    g25   47.856093      -0.354          0.272996  -1.296944    0.194650       0.295608
    """)
    return pandas.read_table(StringIO(data), delim_whitespace=True)


@pytest.fixture
def gff():
    return pybedtools.BedTool(
        """
        chr1  . gene 1 10    . + . ID="g1"
        chr1  . gene 20 30   . + . ID="g2"
        chr1  . gene 35 40   . + . ID="g3"
        chr1  . gene 45 50   . + . ID="g4"
        chr1  . gene 60 80   . + . ID="g5"
        chr1  . gene 90 100  . + . ID="g6"
        chr1  . gene 101 110 . + . ID="g7"
        chr1  . gene 105 120 . + . ID="g8"
        chr1  . gene 150 350 . + . ID="g9"
        chr1  . gene 400 500 . + . ID="g10"
        chr2  . gene 1 10    . + . ID="g11"
        chr2  . gene 20 30   . + . ID="g12"
        chr2  . gene 35 40   . + . ID="g13"
        chr2  . gene 45 50   . + . ID="g14"
        chr2  . gene 60 80   . + . ID="g15"
        chr2  . gene 90 100  . + . ID="g16"
        chr2  . gene 101 110 . + . ID="g17"
        chr2  . gene 105 120 . + . ID="g18"
        chr2  . gene 150 350 . + . ID="g19"
        chr2  . gene 400 500 . + . ID="g20"
        chr3  . gene 1 100   . + . ID="g21"
        chr3  . gene 10 50   . + . ID="g22"
        chr3  . gene 25 30   . + . ID="g23"
        chr3  . gene 45 50   . + . ID="g24"
        chr3  . gene 75 100  . + . ID="g25"
        """, from_string=True).fn


@pytest.fixture
def gffdb(gff):
    return gffutils.create_db(gff, ":memory:")


@pytest.fixture
def rt(deseq_results):
    "an easy-to-type fixture, because laziness"
    return results_table.ResultsTable(deseq_results)


def test_up_down(deseq_results):
    r = results_table.DESeq2Results(deseq_results)
    assert list(r.upregulated(idx=False).index) == [ 'g1', 'g3', 'g5', 'g18']
    assert list(r.upregulated(alpha=0.01, idx=False).index) == [ 'g1', 'g3', 'g5']
    assert list(r.upregulated(alpha=0.01, lfc=3.5, idx=False).index) == ['g3']
    assert list(r.downregulated(idx=False).index) == [ 'g2', 'g4', 'g10', 'g11', 'g14']
    assert list(r.downregulated(alpha=0.01, idx=False).index) == ['g2', 'g4', 'g11', 'g14']
    assert list(r.downregulated(alpha=0.01, lfc=-3.5, idx=False).index) == ['g2']


def test_reindex_to(rt):
    r2 = rt.reindex_to(
        pybedtools.BedTool(
            """
            chr1 0 100 g23
            chr1 0 100 g10
            """, from_string=True))
    assert list(r2.index) == ['g23', 'g10']

    r2 = rt.reindex_to(
        pybedtools.BedTool(
            """
            chr1 0 100 . . . g6
            chr1 0 100 . . . g17
            """, from_string=True), attribute=6)
    assert list(r2.index) == ['g6', 'g17']

    # check that ID works
    r2 = rt.reindex_to(
        pybedtools.BedTool(
            """
            chr1 ucb gene 631 899  . + . ID=g5;Name=asdf;
            chr1 ucb gene 631 913  . + . ID=g2;Parent=other;Name=zz;
            """, from_string=True))
    assert list(r2.index) == ['g5', 'g2']

    r2 = rt.reindex_to(
        pybedtools.BedTool(
            """
            chr1 ucb gene 631 899  . + . Name=g5;
            chr1 ucb gene 631 913  . + . Name=g2;Parent=other;
            """, from_string=True), 'name')
    assert list(r2.index) == ['g5', 'g2']

    r2 = rt.reindex_to(
        pybedtools.BedTool(
            """
            chr1 ucb gene 631 899  . + . XYZ=g5;
            chr1 ucb gene 631 913  . + . XYZ=g2;Parent=other;
            """, from_string=True), 'XYZ')
    assert list(r2.index) == ['g5', 'g2']


def test_dataframe_access(rt):
    """
    underlying dataframe should be easily accessible
    """
    assert rt.padj is rt.data.padj
    assert rt['padj'] is rt.data.padj


def test_subsetting(rt):
    assert all(rt[:10].data == rt.data[:10])
    assert all(rt.update(rt.data[:10]).data == rt.data[:10])


def test_copy(rt):
    e = rt.copy()
    e.data['padj'] = False
    assert sum(e.padj == 0)
    assert sum(rt.padj > 0)


def test_features(gffdb, rt):
    rt.attach_db(gffdb)
    assert len(pybedtools.BedTool(rt.features())) == 25


def test_scatter(rt):
    # just a smoke test
    rt.scatter(
        x='padj',
        y='baseMean'
    )
