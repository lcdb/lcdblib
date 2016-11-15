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

    with pytest.raises(ValueError):
        aligners.prefix_from_hisat2_index(['a/b.1.ht2', 'z/b.2.ht2'])



