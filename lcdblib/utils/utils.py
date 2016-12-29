from collections.abc import Iterable
import contextlib

@contextlib.contextmanager
def temp_env(env):
    """
    Context manager to temporarily set os.environ.
    """
    env = dict(env)
    orig = os.environ.copy()
    _env = {k: str(v) for k, v in env.items()}
    os.environ.update(_env)
    try:
        yield
    finally:
        os.environ.clear()
        os.environ.update(orig)

def flatten(iter):
    """
    Flatten an arbitrarily nested iterable whose innermost items are strings
    into a flat list of strings.
    """
    if isinstance(iter, dict):
        iter = iter.values()
    def gen():
        for item in iter:
            if isinstance(item, dict):
                item = item.values()
            if isinstance(item, Iterable) and not isinstance(item, str):
                yield from flatten(item)
            else:
                yield item
    return list(gen())


def test_flatten():
    assert sorted(flatten({
        'a': {
            'b': {
                'c': ['a','b','c'],
            },
        },
        'x': ['e', 'f', 'g'],
        'y': {
            'z': 'd'
        },
    })) == ['a', 'b', 'c', 'd', 'e', 'f', 'g']
