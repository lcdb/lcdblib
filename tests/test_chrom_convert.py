import os
import pytest

from lcdblib.utils import chrom_convert

def test_import_conversion_FlyBase_to_UCSC():
    mapping = chrom_convert.import_conversion('FlyBase', 'UCSC')
    assert mapping['2L'] == 'chr2L'
    assert mapping['3L'] == 'chr3L'
    assert mapping['mitochondrion_genome'] == 'chrM'

def test_import_conversion_FlyBase_to_GenBank():
    mapping = chrom_convert.import_conversion('FlyBase', 'GenBank')
    assert mapping['2L'] == 'AE014134.6'
    assert mapping['3L'] == 'AE014296.5'
    assert mapping['mitochondrion_genome'] == 'KJ947872.2'

def test_import_conversion_FlyBase_to_RefSeq():
    mapping = chrom_convert.import_conversion('FlyBase', 'RefSeq')
    assert mapping['2L'] == 'NT_033779.5'
    assert mapping['3L'] == 'NT_037436.4'
    assert mapping['mitochondrion_genome'] == 'NC_024511.2'

def test_import_conversion_UCSC_to_FlyBase():
    mapping = chrom_convert.import_conversion('UCSC', 'FlyBase')
    assert mapping['chr2L'] == '2L'
    assert mapping['chr3L'] == '3L'
    assert mapping['chrM'] == 'mitochondrion_genome'

def test_import_conversion_UCSC_to_GenBank():
    mapping = chrom_convert.import_conversion('UCSC', 'GenBank')
    assert mapping['chr2L'] == 'AE014134.6'
    assert mapping['chr3L'] == 'AE014296.5'
    assert mapping['chrM'] == 'KJ947872.2'

def test_import_conversion_UCSC_to_RefSeq():
    mapping = chrom_convert.import_conversion('UCSC', 'RefSeq')
    assert mapping['chr2L'] == 'NT_033779.5'
    assert mapping['chr3L'] == 'NT_037436.4'
    assert mapping['chrM'] == 'NC_024511.2'

def test_import_conversion_GenBank_to_FlyBase():
    mapping = chrom_convert.import_conversion('GenBank', 'FlyBase')
    assert mapping['AE014134.6'] == '2L'
    assert mapping['AE014296.5'] == '3L'
    assert mapping['KJ947872.2'] == 'mitochondrion_genome'

def test_import_conversion_GenBank_to_UCSC():
    mapping = chrom_convert.import_conversion('GenBank', 'UCSC')
    assert mapping['AE014134.6'] == 'chr2L'
    assert mapping['AE014296.5'] == 'chr3L'
    assert mapping['KJ947872.2'] == 'chrM'

def test_import_conversion_GenBank_RefSeq():
    mapping = chrom_convert.import_conversion('GenBank', 'RefSeq')
    assert mapping['AE014134.6'] == 'NT_033779.5'
    assert mapping['AE014296.5'] == 'NT_037436.4'
    assert mapping['KJ947872.2'] == 'NC_024511.2'

def test_import_conversion_RefSeq_to_FlyBase():
    mapping = chrom_convert.import_conversion('RefSeq', 'FlyBase')
    assert mapping['NT_033779.5'] == '2L'
    assert mapping['NT_037436.4'] == '3L'
    assert mapping['NC_024511.2'] == 'mitochondrion_genome'

def test_import_conversion_RefSeq_to_UCSC():
    mapping = chrom_convert.import_conversion('RefSeq', 'UCSC')
    assert mapping['NT_033779.5'] == 'chr2L'
    assert mapping['NT_037436.4'] == 'chr3L'
    assert mapping['NC_024511.2'] == 'chrM'

def test_import_conversion_RefSeq_GenBank():
    mapping = chrom_convert.import_conversion('RefSeq', 'GenBank')
    assert mapping['NT_033779.5'] == 'AE014134.6'
    assert mapping['NT_037436.4'] == 'AE014296.5'
    assert mapping['NC_024511.2'] == 'KJ947872.2'


