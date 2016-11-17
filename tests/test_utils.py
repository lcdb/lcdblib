from lcdblib.utils import imports
import sys
import os
from textwrap import dedent

def test_resolve_name(tmpdir):
    t = str(tmpdir)
    os.makedirs(os.path.join(t, 'a', 'b'))
    with open(os.path.join(t, 'a', '__init__.py'), 'w'):
        pass
    with open(os.path.join(t, 'a', 'b', '__init__.py'), 'w'):
        pass
    with open(os.path.join(t, 'a', 'b', 'c.py'), 'w') as fout:
        fout.write(dedent(
            '''
            def x():
                return 'this is the result'
            '''))
    orig_path = sys.path[:]
    sys.path.insert(0, t)
    x = imports.resolve_name('a.b.c.x')
    assert x() == 'this is the result'
    sys.path = orig_path
