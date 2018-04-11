import os
import subprocess
import pytest
from textwrap import dedent
import gzip

from lcdblib.utils import chrom_convert


class TestImport:
    """ Test conversion table import and construction of a mapping dictionary. """
    def test_import_conversion_FlyBase_to_UCSC(self):
        mapping = chrom_convert.import_conversion('FlyBase', 'UCSC')
        assert mapping['2L'] == 'chr2L'
        assert mapping['3L'] == 'chr3L'
        assert mapping['mitochondrion_genome'] == 'chrM'

    def test_import_conversion_FlyBase_to_GenBank(self):
        mapping = chrom_convert.import_conversion('FlyBase', 'GenBank')
        assert mapping['2L'] == 'AE014134.6'
        assert mapping['3L'] == 'AE014296.5'
        assert mapping['mitochondrion_genome'] == 'KJ947872.2'

    def test_import_conversion_FlyBase_to_RefSeq(self):
        mapping = chrom_convert.import_conversion('FlyBase', 'RefSeq')
        assert mapping['2L'] == 'NT_033779.5'
        assert mapping['3L'] == 'NT_037436.4'
        assert mapping['mitochondrion_genome'] == 'NC_024511.2'

    def test_import_conversion_UCSC_to_FlyBase(self):
        mapping = chrom_convert.import_conversion('UCSC', 'FlyBase')
        assert mapping['chr2L'] == '2L'
        assert mapping['chr3L'] == '3L'
        assert mapping['chrM'] == 'mitochondrion_genome'

    def test_import_conversion_UCSC_to_GenBank(self):
        mapping = chrom_convert.import_conversion('UCSC', 'GenBank')
        assert mapping['chr2L'] == 'AE014134.6'
        assert mapping['chr3L'] == 'AE014296.5'
        assert mapping['chrM'] == 'KJ947872.2'

    def test_import_conversion_UCSC_to_RefSeq(self):
        mapping = chrom_convert.import_conversion('UCSC', 'RefSeq')
        assert mapping['chr2L'] == 'NT_033779.5'
        assert mapping['chr3L'] == 'NT_037436.4'
        assert mapping['chrM'] == 'NC_024511.2'

    def test_import_conversion_GenBank_to_FlyBase(self):
        mapping = chrom_convert.import_conversion('GenBank', 'FlyBase')
        assert mapping['AE014134.6'] == '2L'
        assert mapping['AE014296.5'] == '3L'
        assert mapping['KJ947872.2'] == 'mitochondrion_genome'

    def test_import_conversion_GenBank_to_UCSC(self):
        mapping = chrom_convert.import_conversion('GenBank', 'UCSC')
        assert mapping['AE014134.6'] == 'chr2L'
        assert mapping['AE014296.5'] == 'chr3L'
        assert mapping['KJ947872.2'] == 'chrM'

    def test_import_conversion_GenBank_RefSeq(self):
        mapping = chrom_convert.import_conversion('GenBank', 'RefSeq')
        assert mapping['AE014134.6'] == 'NT_033779.5'
        assert mapping['AE014296.5'] == 'NT_037436.4'
        assert mapping['KJ947872.2'] == 'NC_024511.2'

    def test_import_conversion_RefSeq_to_FlyBase(self):
        mapping = chrom_convert.import_conversion('RefSeq', 'FlyBase')
        assert mapping['NT_033779.5'] == '2L'
        assert mapping['NT_037436.4'] == '3L'
        assert mapping['NC_024511.2'] == 'mitochondrion_genome'

    def test_import_conversion_RefSeq_to_UCSC(self):
        mapping = chrom_convert.import_conversion('RefSeq', 'UCSC')
        assert mapping['NT_033779.5'] == 'chr2L'
        assert mapping['NT_037436.4'] == 'chr3L'
        assert mapping['NC_024511.2'] == 'chrM'

    def test_import_conversion_RefSeq_GenBank(self):
        mapping = chrom_convert.import_conversion('RefSeq', 'GenBank')
        assert mapping['NT_033779.5'] == 'AE014134.6'
        assert mapping['NT_037436.4'] == 'AE014296.5'
        assert mapping['NC_024511.2'] == 'KJ947872.2'


@pytest.fixture(scope='session')
def mapper():
        return chrom_convert.import_conversion('UCSC', 'FlyBase')


@pytest.fixture(scope='session')
def inputs(tmpdir_factory):
    import shutil
    import pybedtools

    d = str(tmpdir_factory.mktemp('files'))

    shutil.copy(pybedtools.example_filename('x.bam'), os.path.join(d, 'x.bam'))
    shutil.copy(pybedtools.example_filename('x.bed'), os.path.join(d, 'x.bed'))
    shutil.copy(pybedtools.example_filename('gdc.gff'), os.path.join(d, 'x.gff'))

    # Make SAM file
    cmd = ['samtools', 'view', '-h', os.path.join(d, 'x.bam'), '-O', 'SAM', '-o', os.path.join(d, 'x.sam')]
    subprocess.run(cmd)

    # FASTA with basic header
    fa = dedent("""\
        >chr2L
        TAATTAAAACAGATCCTGAGAAAATTTCCACAATTATGAAGTATCCGATTCCACAAAACATTAGAGAGCTTCGAAGTTTT
        CTAGGCCTCACCGGCTACTACCGTAAATTTGTCCGAAATTATGCAAAAATTGCCAAACCTCTAACCAAATACTTAGGAGG
        AAATAATGGAAAAATTTCTAGGAGAATGTCTACAAAAATTAAAATACAGTTAGATGACCCAGCTGTTAAAGCTTTTAACG
        AACTTAAAGATAATTTAATAGCACAAGTGGAATTAGTTCAACCTGATTATAACAAAAAAATTCACTTTAACGACAGACGC
        """)

    fname = os.path.join(d, 'dm6_basic.fa')
    with open(fname, 'w') as fh:
        fh.write(fa)

    # Make GZIP FILE
    with gzip.open(os.path.join(d, 'dm6_basic.fa.gz'), 'wb') as fh:
        fh.write(fa.encode())

    # FASTA with full header
    fa = dedent("""\
        >chr2L type=golden_path_region; loc=chr2L:1..14983; ID=chr2L REFSEQ:NW_001845015; length=14983; release=r6.09; species=Dmel;
        TAATTAAAACAGATCCTGAGAAAATTTCCACAATTATGAAGTATCCGATTCCACAAAACATTAGAGAGCTTCGAAGTTTT
        CTAGGCCTCACCGGCTACTACCGTAAATTTGTCCGAAATTATGCAAAAATTGCCAAACCTCTAACCAAATACTTAGGAGG
        AAATAATGGAAAAATTTCTAGGAGAATGTCTACAAAAATTAAAATACAGTTAGATGACCCAGCTGTTAAAGCTTTTAACG
        AACTTAAAGATAATTTAATAGCACAAGTGGAATTAGTTCAACCTGATTATAACAAAAAAATTCACTTTAACGACAGACGC
        """)

    fname = os.path.join(d, 'dm6_full.fa')
    with open(fname, 'w') as fh:
        fh.write(fa)

    return d


def test_pysam_convert_BAM(inputs, mapper):
    bam = os.path.join(inputs, 'x.bam')
    oname = os.path.join(inputs, 'x_convert.bam')
    chrom_convert.pysam_convert(bam, oname, 'BAM', mapper)
    chrom = subprocess.run(('samtools view {} | head -n1'.format(os.path.join(inputs, 'x_convert.bam'))),
                           stdout=subprocess.PIPE, shell=True).stdout.decode().split('\t')[2]
    # FIXME:
    # assert chrom == '2L'



def test_pysam_convert_SAM(inputs, mapper):
    sam = os.path.join(inputs, 'x.sam')
    oname = os.path.join(inputs, 'x_convert.sam')
    chrom_convert.pysam_convert(sam, oname, 'SAM', mapper)
    chrom = subprocess.run(('tail -n1 {}'.format(os.path.join(inputs, 'x_convert.sam'))),
                           stdout=subprocess.PIPE, shell=True).stdout.decode().split('\t')[2]
    # FIXME:
    # assert chrom == '2L'

def test_pysam_convert_SAM_PIPE(inputs, mapper):
    sam = os.path.join(inputs, 'x.sam')
    cmd = 'cat {sam} | chrom_convert -i - --from UCSC --to FlyBase \
            --fileType SAM | grep -v -e "^@" | head -n1 | cut -f3'.format(sam=sam)

    out = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

    # FIXME:
    # assert out.stdout.decode('UTF-8').strip() == '2L'

def test_pybedtools_convert_BED(inputs, mapper):
    bed = os.path.join(inputs, 'x.bed')
    oname = os.path.join(inputs, 'x_convert.bed')
    chrom_convert.pybedtools_convert(bed, oname, mapper)
    chrom = subprocess.run(('head -n1 {}'.format(os.path.join(inputs, 'x_convert.bed'))),
                           stdout=subprocess.PIPE, shell=True).stdout.decode().split('\t')[0]
    assert chrom == '2L'


def test_pybedtools_convert_GFF(inputs, mapper):
    gff = os.path.join(inputs, 'x.gff')
    oname = os.path.join(inputs, 'x_convert.gff')
    chrom_convert.pybedtools_convert(gff, oname, mapper)
    chrom = subprocess.run(('head -n1 {}'.format(os.path.join(inputs, 'x_convert.gff'))),
                           stdout=subprocess.PIPE, shell=True).stdout.decode().split('\t')[0]
    assert chrom == '2L'


def test_pybedtools_convert_GFF_PIPE(inputs, mapper):
    gff = os.path.join(inputs, 'x.gff')
    cmd = 'cat {gff} | chrom_convert -i - --from UCSC --to FlyBase \
            --fileType GFF | head -n1 | cut -f1'.format(gff=gff)

    out = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

    assert out.stdout.decode('UTF-8').strip() == '2L'


def test_fasta_convert_basic_header(inputs, mapper):
    fname = os.path.join(inputs, 'dm6_basic.fa')
    oname = os.path.join(inputs, 'dm6_basic_convert.fa')
    chrom_convert.fasta_convert(fname, oname, mapper)
    with open(oname, 'r') as fh:
        header = '>2L'
        assert header == fh.readline().strip()


def test_fasta_convert_basic_header_PIPE(inputs, mapper):
    fasta = os.path.join(inputs, 'dm6_basic.fa')
    cmd = 'cat {fasta} | chrom_convert -i - --from UCSC --to FlyBase \
            --fileType FASTA | head -n1'.format(fasta=fasta)

    out = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

    assert out.stdout.decode('UTF-8').strip() == '>2L'


def test_fasta_convert_full_header(inputs, mapper):
    fname = os.path.join(inputs, 'dm6_full.fa')
    oname = os.path.join(inputs, 'dm6_full_convert.fa')
    chrom_convert.fasta_convert(fname, oname, mapper)
    with open(oname, 'r') as fh:
        header = '>2L type=golden_path_region; loc=2L:1..14983; ID=2L REFSEQ:NW_001845015; length=14983; release=r6.09; species=Dmel;'
        assert header == fh.readline().strip()


def test_fasta_convert_basic_header(inputs, mapper):
    fname = os.path.join(inputs, 'dm6_basic.fa.gz')
    oname = os.path.join(inputs, 'dm6_basic_convert.fa')
    chrom_convert.fasta_convert(fname, oname, mapper)

    with open(oname, 'r') as fh:
        header = '>2L'
        assert header == fh.readline().strip()
