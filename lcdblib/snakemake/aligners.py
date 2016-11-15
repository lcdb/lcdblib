"""
Helper functions for working with aligners within Snakefiles
"""

def hisat2_index_from_prefix(prefix):
    """
    Given a prefix, return a list of the corresponding hisat2 index files.
    """
    return ['{prefix}.{n}.ht2'.format(prefix=prefix, n=n) for n in range(1, 9)]


def prefix_from_hisat2_index(index_files):
    """
    Given a list of index files for hisat2, return the corresponding prefix.
    """
    prefixes = list(set(map(lambda x: '.'.join(x.split('.')[:-2]), index_files)))
    if len(prefixes) != 1:
        raise ValueError("More than one prefix detected from '{0}'".format(prefixes))
    return prefixes[0]
