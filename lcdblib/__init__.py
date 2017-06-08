# -*- coding: utf-8 -*-
from ._version import get_version
version, build = get_version()
__version__ = '{0}.{1}'.format(version, build)
