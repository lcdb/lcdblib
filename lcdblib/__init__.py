from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

# The following assumes versioneer was configured with the "pep440-pre" style
# in setup.cfg
__conda_version__, __conda_build__ = __version__.split('.dev')
if __conda_build__ == "":
    __conda_build__ = '0'
