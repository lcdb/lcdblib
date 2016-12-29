from lcdblib.utils import imports, utils
import sys
import os
from textwrap import dedent


def test_temp_env():
    assert 'new_key_testing' not in os.environ
    with utils.temp_env(dict(new_key_testing='asdf')):
        assert os.environ['new_key_testing'] == 'asdf'
    assert 'new_key_testing' not in os.environ


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
