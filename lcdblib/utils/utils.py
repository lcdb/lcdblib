from collections.abc import Iterable

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
