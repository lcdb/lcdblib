import pytest
from lcdblib.snakemake import aligners

def test_hisat2_prefixes():
    files = [
        'a/b/c.1.ht2',
        'a/b/c.2.ht2',
        'a/b/c.3.ht2',
        'a/b/c.4.ht2',
        'a/b/c.5.ht2',
        'a/b/c.6.ht2',
        'a/b/c.7.ht2',
        'a/b/c.8.ht2']

    assert aligners.hisat2_index_from_prefix('a/b/c') == files
    assert aligners.prefix_from_hisat2_index(files) == 'a/b/c'
    assert aligners.prefix_from_hisat2_index(files[0]) == 'a/b/c'

    with pytest.raises(ValueError):
        aligners.prefix_from_hisat2_index(['a/b.1.ht2', 'z/b.2.ht2'])


def test_bowtie2_prefixes():
    files = [
        'a/b/c.1.bt2',
        'a/b/c.2.bt2',
        'a/b/c.3.bt2',
        'a/b/c.4.bt2',
        'a/b/c.rev.1.bt2',
        'a/b/c.rev.2.bt2',
    ]

    assert aligners.bowtie2_index_from_prefix('a/b/c') == files
    assert aligners.prefix_from_bowtie2_index(files) == 'a/b/c'
    assert aligners.prefix_from_bowtie2_index(files[0]) == 'a/b/c'

    with pytest.raises(ValueError):
        aligners.prefix_from_bowtie2_index(['a/b.1.bt2', 'z/b.2.bt2'])
