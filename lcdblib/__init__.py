from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

# The following assumes versioneer was configured with the "pep440-pre" style
# in setup.cfg
toks = __version__.split('.dev')
if len(toks) == 2:
    __conda_version__, __conda_build__ = toks
else:
    __conda_build__ = '0'
    __conda_version__ = toks[0]
