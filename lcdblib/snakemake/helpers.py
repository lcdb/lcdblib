import collections
from snakemake.shell import shell
from snakemake.io import expand

def fill_patterns(patterns, fill):
    """
    Fills in a dictionary of patterns with the dictionary `fill`.

    >>> patterns = dict(a='{sample}_R{N}.fastq')
    >>> fill = dict(sample=['one', 'two'], N=[1, 2])
    >>> sorted(fill_patterns(patterns, fill)['a'])
    ['one_R1.fastq', 'one_R2.fastq', 'two_R1.fastq', 'two_R2.fastq']
    """

    def update(d, u):
        for k, v in u.items():
            if isinstance(v, collections.Mapping):
                r = update(d.get(k, {}), v)
                d[k] = r
            else:
                d[k] = expand(u[k], **fill)
        return d
    d = {}
    return update(d, patterns)


def rscript(string, scriptname, log=None):
    """
    Saves the string as `scriptname` and then runs it

    Parameters
    ----------
    string : str
        Filled-in template to be written as R script

    scriptname : str
        File to save script to

    log : str
        File to redirect stdout and stderr to. If None, no redirection occurs.
    """
    with open(scriptname, 'w') as fout:
        fout.write(string)
    if log:
        _log = '> {0} 2>&1'.format(log)
    else:
        _log = ""
    shell('Rscript {scriptname} {_log}')
